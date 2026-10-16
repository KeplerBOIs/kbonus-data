"""SQLite index of every downloaded light curve and the detector channel of each quarter.

The per-quarter extensions of a KBonus light curve record where the source fell
on the detector (MODULE, OUTPUT, CHANNEL) and aperture metrics (CROWDSAP,
FLFRCSAP, ...).  None of this is in the catalogs, so the index is built by
reading the headers of every file.  It also stores the byte offset of each
per-quarter extension, so that a quarter can later be read with a single read
instead of walking every header in the file (see :mod:`kbonus_data.channels`).

Tables::

    files         one row per source: target, path, kic, gaia_id, ra, dec, mags
    observations  one row per source and quarter: channel, module, output,
                  aperture metrics, number of cadences, extension offsets
    errors        files that could not be read in the last build

The database lives at ``<data dir>/index/kbonus_index.sqlite`` by default.  It
stores file paths relative to the data directory, never absolute ones, so the
data directory can be moved or copied with its index, or the index used on its
own with the data elsewhere: ``Index(root=<data dir>, path=<index file>)``.  Building it
reads ~20 small blocks from each of ~575,000 files; on the NFS-mounted data
directory that is I/O-bound at a few tens of files per second, so the build
is resumable and only scans files it has not seen.
"""

import os
import re
import sqlite3
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import MaskedColumn, Table

from .paths import _FILENAME_RE, _version_key, data_dir

INDEX_NAME = "index/kbonus_index.sqlite"
SCHEMA_VERSION = 1

BLOCK = 2880
READ_SIZE = 32768  # one read covers an aperture HDU and the next quarter's header

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS files (
    id INTEGER PRIMARY KEY,
    target TEXT UNIQUE NOT NULL,   -- kic-008311864 or gaia-dr3-<id>
    path TEXT NOT NULL,            -- relative to the data directory
    version TEXT,
    size INTEGER,
    mtime REAL,
    kic INTEGER,
    gaia_id INTEGER,
    ra REAL,
    dec REAL,
    kepmag REAL,
    gmag REAL,
    bpmag REAL,
    rpmag REAL,
    tpforg INTEGER,
    qdetect TEXT,                  -- detection flag per quarter, one character each
    n_quarters INTEGER
);
CREATE TABLE IF NOT EXISTS observations (
    file_id INTEGER NOT NULL REFERENCES files(id),
    quarter INTEGER NOT NULL,
    season INTEGER,
    module INTEGER,
    output INTEGER,
    channel INTEGER,
    n_cadences INTEGER,
    qdetect INTEGER,
    tpf_row INTEGER,
    tpf_col INTEGER,
    row REAL,
    col REAL,
    sap TEXT,
    flfrcsap REAL,
    crowdsap REAL,
    npixsap INTEGER,
    psffrac REAL,
    pertrati REAL,
    pertstd REAL,
    iterneg INTEGER,
    hdu_offset INTEGER,            -- byte offset of the LIGHTCURVE_Q<n> header
    hdu_size INTEGER,              -- bytes of header + padded data
    aperture_offset INTEGER,
    aperture_size INTEGER,
    PRIMARY KEY (file_id, quarter)
) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS obs_quarter_channel ON observations (quarter, channel);
CREATE INDEX IF NOT EXISTS obs_channel ON observations (channel);
CREATE INDEX IF NOT EXISTS files_kic ON files (kic);
CREATE TABLE IF NOT EXISTS errors (path TEXT PRIMARY KEY, error TEXT, time TEXT);
"""

# header keyword -> observations column
QUARTER_KEYWORDS = {
    "QUARTER": "quarter",
    "SEASON": "season",
    "MODULE": "module",
    "OUTPUT": "output",
    "CHANNEL": "channel",
    "NAXIS2": "n_cadences",
    "TPF_ROW": "tpf_row",
    "TPF_COL": "tpf_col",
    "ROW": "row",
    "COLUMN": "col",
    "SAP": "sap",
    "FLFRCSAP": "flfrcsap",
    "CROWDSAP": "crowdsap",
    "NPIXSAP": "npixsap",
    "PSFFRAC": "psffrac",
    "PERTRATI": "pertrati",
    "PERTSTD": "pertstd",
    "ITERNEG": "iterneg",
}
OBS_COLUMNS = (
    "file_id",
    *QUARTER_KEYWORDS.values(),
    "qdetect",
    "hdu_offset",
    "hdu_size",
    "aperture_offset",
    "aperture_size",
)
FILE_COLUMNS = (
    "target", "path", "version", "size", "mtime", "kic", "gaia_id", "ra", "dec",
    "kepmag", "gmag", "bpmag", "rpmag", "tpforg", "qdetect", "n_quarters",
)  # fmt: skip

_END_CARD = b"END" + b" " * 77
_GAIA_ID_RE = re.compile(r"(\d{10,})")


def index_path(root=None):
    return data_dir(root) / INDEX_NAME


# --- reading headers ---------------------------------------------------------


class _Reader:
    """``pread`` with a one-buffer cache, so nearby headers cost a single read."""

    def __init__(self, fd):
        self.fd, self.start, self.buf = fd, 0, b""

    def read(self, offset, n):
        if not (self.start <= offset and offset + n <= self.start + len(self.buf)):
            self.buf, self.start = os.pread(self.fd, max(n, READ_SIZE), offset), offset
        return self.buf[offset - self.start : offset - self.start + n]


def _padded(n):
    return -(-n // BLOCK) * BLOCK


def walk_hdus(path):
    """Yield ``(header, header offset, total HDU size)`` for each HDU of a FITS file.

    Reads only header blocks, skipping over data.  Much cheaper than
    ``astropy.io.fits.open`` on network file systems, where astropy's reads
    pull in most of each file.
    """
    fd = os.open(path, os.O_RDONLY)
    try:
        size = os.fstat(fd).st_size
        reader = _Reader(fd)
        offset = 0
        while offset + BLOCK <= size:
            start, cards = offset, []
            while True:
                block = reader.read(offset, BLOCK)
                if len(block) < BLOCK:
                    raise OSError(f"{path}: truncated header at byte {offset}")
                offset += BLOCK
                cards.append(block)
                if any(block[i : i + 80] == _END_CARD for i in range(0, BLOCK, 80)):
                    break
            header = fits.Header.fromstring(b"".join(cards).decode("ascii"))
            naxis = header.get("NAXIS", 0)
            n = int(np.prod([header[f"NAXIS{i}"] for i in range(1, naxis + 1)])) if naxis else 0
            nbytes = abs(header.get("BITPIX", 8)) // 8 * header.get("GCOUNT", 1) * (header.get("PCOUNT", 0) + n)
            offset += _padded(nbytes)
            yield header, start, offset - start
    finally:
        os.close(fd)


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _float(value):
    """Float or None; the files use 1e20 for a missing magnitude."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) and abs(value) < 1e19 else None


def scan_file(path):
    """Return ``(file row, [observation rows])`` read from the headers of one file.

    Rows are dicts keyed by column name; observation rows lack ``file_id``.
    """
    path = Path(path)
    hdus = list(walk_hdus(path))
    primary = hdus[0][0]
    qdetect = str(primary.get("QDETECT", "") or "")
    observations, apertures = {}, {}
    for header, offset, nbytes in hdus[1:]:
        name = header.get("EXTNAME", "")
        if name.startswith("LIGHTCURVE_Q"):
            obs = {column: header.get(keyword) for keyword, column in QUARTER_KEYWORDS.items()}
            quarter = obs["quarter"] = int(name.removeprefix("LIGHTCURVE_Q"))
            flag = qdetect[quarter] if quarter < len(qdetect) else None
            obs.update(qdetect=_int(flag), hdu_offset=offset, hdu_size=nbytes)
            observations[quarter] = obs
        elif name.startswith("APERTURE_Q"):
            apertures[int(name.removeprefix("APERTURE_Q"))] = (offset, nbytes)
    for quarter, obs in observations.items():
        obs["aperture_offset"], obs["aperture_size"] = apertures.get(quarter, (None, None))

    m = _FILENAME_RE.search(path.name)
    gaia = _GAIA_ID_RE.search(str(primary.get("GAIAID", "")))
    stat = path.stat()
    row = {
        "target": m.group(1),
        "version": m.group(2),
        "size": stat.st_size,
        "mtime": stat.st_mtime,
        "kic": _int(primary.get("KEPLERID")) or None,
        "gaia_id": int(gaia.group(1)) if gaia else None,
        "ra": _float(primary.get("RA_OBJ")),
        "dec": _float(primary.get("DEC_OBJ")),
        "kepmag": _float(primary.get("KEPMAG")),
        "gmag": _float(primary.get("GMAG")),
        "bpmag": _float(primary.get("BPMAG")),
        "rpmag": _float(primary.get("RPMAG")),
        "tpforg": _int(primary.get("TPFORG")),
        "qdetect": qdetect or None,
        "n_quarters": len(observations),
    }
    return row, [observations[q] for q in sorted(observations)]


# --- building ----------------------------------------------------------------


def connect(path, readonly=False):
    path = Path(path)
    if readonly:
        if not path.exists():
            raise FileNotFoundError(
                f"index {path.resolve()} not found. Point root= (or $KBONUS_DATA_DIR) at the data "
                "directory that holds index/, pass path= to the index file, or build it with: "
                "kbonus-data index build"
            )
        # wait out the commits of a build that is running at the same time
        return sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=60)
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path)
    con.executescript(SCHEMA)
    con.execute("INSERT OR IGNORE INTO meta VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))
    return con


def lc_files(root=None, prefixes=None):
    """Return the relative path of each target's highest-version light curve file."""
    root = data_dir(root)
    best = {}
    lcs = root / "lcs"
    for directory in sorted(prefixes or (e.name for e in os.scandir(lcs) if e.is_dir())):
        for entry in os.scandir(lcs / directory):
            m = _FILENAME_RE.search(entry.name)
            if not m:
                continue
            current = best.get(m.group(1))
            if current is None or _version_key(Path(entry.name)) > _version_key(Path(current)):
                best[m.group(1)] = f"lcs/{directory}/{entry.name}"
    return sorted(best.values())


def _insert(con, row, observations):
    old = con.execute("SELECT id FROM files WHERE target = ?", (row["target"],)).fetchone()
    if old:
        con.execute("DELETE FROM observations WHERE file_id = ?", old)
        con.execute("DELETE FROM files WHERE id = ?", old)
    cursor = con.execute(
        f"INSERT INTO files ({', '.join(FILE_COLUMNS)}) VALUES ({', '.join('?' * len(FILE_COLUMNS))})",
        [row[c] for c in FILE_COLUMNS],
    )
    file_id = cursor.lastrowid
    con.executemany(
        f"INSERT INTO observations ({', '.join(OBS_COLUMNS)}) VALUES ({', '.join('?' * len(OBS_COLUMNS))})",
        [[file_id] + [obs[c] for c in OBS_COLUMNS[1:]] for obs in observations],
    )
    con.execute("DELETE FROM errors WHERE path = ?", (row["path"],))


def build_index(root=None, db_path=None, workers=12, prefixes=None, limit=None, rebuild=False, batch=1000, report_every=60):
    """Scan light curve headers into the index.  Resumable.

    Files already in the index (by path) are skipped unless ``rebuild``.  A
    target whose path changed (e.g. a new version was downloaded) is rescanned
    and replaces its old entry.  Files that fail are recorded in the
    ``errors`` table and retried on the next run.  Returns the number of
    files scanned.
    """
    root = data_dir(root)
    db_path = Path(db_path or index_path(root))
    con = connect(db_path)
    print(f"index {db_path}", flush=True)
    paths = lc_files(root, prefixes)
    done = set() if rebuild else {p for (p,) in con.execute("SELECT path FROM files")}
    todo = [p for p in paths if p not in done]
    print(f"{len(paths)} light curve files, {len(paths) - len(todo)} indexed, {len(todo)} to scan", flush=True)
    if limit is not None:
        todo = todo[:limit]

    def scan(relative):
        try:
            row, observations = scan_file(root / relative)
            row["path"] = relative
            return relative, row, observations, None
        except Exception as exc:  # recorded and retried next run
            return relative, None, None, repr(exc)

    start = last = time.monotonic()
    scanned = failed = 0
    with ThreadPoolExecutor(workers) as pool:
        for i in range(0, len(todo), batch):
            for relative, row, observations, error in pool.map(scan, todo[i : i + batch]):
                if error is None:
                    _insert(con, row, observations)
                    scanned += 1
                else:
                    con.execute(
                        "INSERT OR REPLACE INTO errors VALUES (?, ?, ?)",
                        (relative, error, time.strftime("%Y-%m-%dT%H:%M:%S")),
                    )
                    failed += 1
            con.execute(
                "INSERT OR REPLACE INTO meta VALUES ('updated', ?)", (time.strftime("%Y-%m-%dT%H:%M:%S"),)
            )
            con.commit()
            now = time.monotonic()
            if now - last >= report_every or i + batch >= len(todo):
                last = now
                n = scanned + failed
                rate = n / (now - start)
                eta = (len(todo) - n) / rate / 3600 if rate else float("nan")
                print(
                    f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {n}/{len(todo)} files  "
                    f"{rate:.1f} files/s  {failed} failed  ETA {eta:.1f} h",
                    flush=True,
                )
    con.close()
    if failed:
        print(f"{failed} files failed; see the errors table. Re-run to retry them.")
    return scanned


# --- querying ----------------------------------------------------------------


def _to_table(cursor):
    names = [d[0] for d in cursor.description]
    rows = cursor.fetchall()
    table = Table()
    for name, values in zip(names, zip(*rows) if rows else [()] * len(names)):
        missing = [v is None for v in values]
        if any(missing):
            present = [v for v in values if v is not None]
            fill = "" if present and isinstance(present[0], str) else 0
            table[name] = MaskedColumn([fill if m else v for v, m in zip(values, missing)], mask=missing)
        else:
            table[name] = np.array(values) if rows else np.array([], dtype=float)
    return table


# keyword -> SQL column, for range and equality filters in Index.query
_FILTER_COLUMNS = {
    "quarter": "o.quarter",
    "channel": "o.channel",
    "module": "o.module",
    "output": "o.output",
    "season": "o.season",
    "qdetect": "o.qdetect",
    "n_cadences": "o.n_cadences",
    "flfrcsap": "o.flfrcsap",
    "crowdsap": "o.crowdsap",
    "npixsap": "o.npixsap",
    "psffrac": "o.psffrac",
    "pertstd": "o.pertstd",
    "kepmag": "f.kepmag",
    "gmag": "f.gmag",
    "kic": "f.kic",
    "ra": "f.ra",
    "dec": "f.dec",
}


class Index:
    """Read-only access to the light curve index.

    >>> idx = Index()
    >>> stars = idx.query(quarter=5, channel=41, kepmag=(None, 15), crowdsap=(0.9, None))

    Filters take a single value (equality), a list (membership) or a
    ``(low, high)`` tuple (inclusive range; ``None`` leaves a side open).
    """

    def __init__(self, root=None, path=None):
        self.root = data_dir(root)
        self.path = Path(path or index_path(self.root))
        self.con = connect(self.path, readonly=True)

    def close(self):
        self.con.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def sql(self, query, params=()):
        """Run any SQL query and return the result as an astropy Table."""
        return _to_table(self.con.execute(query, params))

    def query(self, where=None, params=(), order_by="COALESCE(f.kepmag, f.gmag)", limit=None, **filters):
        """Return one row per (source, quarter) matching the filters.

        Columns are those of ``observations`` joined with ``files``.  ``where``
        adds a raw SQL condition (tables are aliased ``o`` and ``f``) with
        ``params`` for its placeholders.  Rows are brightest first (Kepler
        magnitude, else Gaia G for sources without a KIC match); a ``kepmag``
        filter drops sources with no Kepler magnitude.
        """
        clauses, values = [], []
        for key, value in filters.items():
            try:
                column = _FILTER_COLUMNS[key]
            except KeyError:
                raise TypeError(f"unknown filter {key!r}; use one of {sorted(_FILTER_COLUMNS)}") from None
            if value is None:
                continue
            if isinstance(value, tuple):
                low, high = value
                if low is not None:
                    clauses.append(f"{column} >= ?")
                    values.append(low)
                if high is not None:
                    clauses.append(f"{column} <= ?")
                    values.append(high)
            elif isinstance(value, (list, set, np.ndarray)):
                value = list(value)
                clauses.append(f"{column} IN ({', '.join('?' * len(value))})")
                values.extend(value)
            else:
                clauses.append(f"{column} = ?")
                values.append(value)
        if where:
            clauses.append(f"({where})")
            values.extend(params)
        query = "SELECT f.target, f.kic, f.gaia_id, f.ra, f.dec, f.kepmag, f.gmag, f.path, o.* " \
                "FROM observations o JOIN files f ON f.id = o.file_id"  # fmt: skip
        if clauses:
            query += " WHERE " + " AND ".join(clauses)
        if order_by:
            query += f" ORDER BY {order_by}"
        if limit is not None:
            query += f" LIMIT {int(limit)}"
        return self.sql(query, values)

    def quarters(self, target):
        """Return the per-quarter detector placement of one target."""
        from .paths import normalize_target

        return self.query(where="f.target = ?", params=(normalize_target(target),), order_by="o.quarter")

    def channel_counts(self, quarter=None):
        """Number of sources per (quarter, channel)."""
        query = "SELECT quarter, channel, COUNT(*) AS n FROM observations"
        params = ()
        if quarter is not None:
            query += " WHERE quarter = ?"
            params = (quarter,)
        return self.sql(query + " GROUP BY quarter, channel ORDER BY quarter, channel", params)

    def summary(self):
        """Return a dict of counts and build metadata."""
        one = lambda q: self.con.execute(q).fetchone()[0]  # noqa: E731
        return {
            "path": str(self.path),
            "files": one("SELECT COUNT(*) FROM files"),
            "observations": one("SELECT COUNT(*) FROM observations"),
            "errors": one("SELECT COUNT(*) FROM errors"),
            "updated": (self.con.execute("SELECT value FROM meta WHERE key = 'updated'").fetchone() or [None])[0],
        }
