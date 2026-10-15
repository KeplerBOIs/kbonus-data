"""Read every light curve on one detector channel in one quarter.

The result is a :class:`ChannelQuarter`: 2-D arrays of shape
``(n_stars, n_cadences)`` aligned on cadence number, plus a table describing
each star.  This is the input for building systematics basis vectors from the
stars that share a channel.

Sources are found with the index (:mod:`kbonus_data.index`), which also gives
the byte offset of each per-quarter extension, so each star costs one read.
"""

import os
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import numpy as np
from astropy.io import fits
from astropy.table import Table

from .focalplane import module_output
from .index import Index
from .lightcurve import read_quarter
from .paths import data_dir

# columns of the LIGHTCURVE_Q<n> tables:
#   CADENCENO TIME FLUX FLUX_ERR SAP_FLUX SAP_FLUX_ERR CENTROID_COLUMN CENTROID_ROW
#   SAP_QUALITY SAP_BKG RED_CHI2 PSF_FLUX_NOVA PSF_FLUX_ERR_NOVA
DEFAULT_COLUMNS = (
    "FLUX",
    "FLUX_ERR",
    "SAP_FLUX",
    "SAP_FLUX_ERR",
    "SAP_BKG",
    "CENTROID_COLUMN",
    "CENTROID_ROW",
    "SAP_QUALITY",
)

# kind -> (flux column, error column) in the per-quarter tables
FLUX_KINDS = {
    "sap": ("SAP_FLUX", "SAP_FLUX_ERR"),
    "psf": ("FLUX", "FLUX_ERR"),
    "psf_nova": ("PSF_FLUX_NOVA", "PSF_FLUX_ERR_NOVA"),
}


@dataclass
class ChannelQuarter:
    """Light curves of the stars on one channel in one quarter.

    ``data[name]`` is an ``(n_stars, n_cadences)`` array for each per-quarter
    column read, also available as a lower-case attribute (``cq.sap_flux``).
    Row ``i`` is the star in ``targets[i]``; column ``j`` is cadence
    ``cadenceno[j]``.  ``cadenceno`` holds every cadence at which any star has
    data.  Where a star has no row for a cadence, ``present`` is False, float
    columns are NaN and integer columns (``SAP_QUALITY``) are 0.

    ``time`` is the median BKJD of each cadence over the stars; per-star
    barycentric times differ by up to ~2 minutes across the field (read the
    ``TIME`` column to get them).
    """

    quarter: int
    channel: int
    targets: Table
    cadenceno: np.ndarray
    time: np.ndarray
    present: np.ndarray
    data: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.targets)

    def __getattr__(self, name):
        data = self.__dict__.get("data", {})
        if name.upper() in data:
            return data[name.upper()]
        raise AttributeError(f"{type(self).__name__} has no attribute or column {name!r}")

    @property
    def shape(self):
        return self.present.shape

    @property
    def module(self):
        return module_output(self.channel)[0]

    @property
    def output(self):
        return module_output(self.channel)[1]

    def get_flux(self, kind="sap"):
        """Return ``(flux, flux_err)`` for ``kind`` in 'sap', 'psf', 'psf_nova'."""
        try:
            names = FLUX_KINDS[kind]
        except KeyError:
            raise ValueError(f"kind must be one of {sorted(FLUX_KINDS)}, not {kind!r}") from None
        missing = [n for n in names if n not in self.data]
        if missing:
            raise KeyError(f"columns {missing} were not read; pass them in columns=")
        return self.data[names[0]], self.data[names[1]]

    def save(self, path):
        """Save to a compressed ``.npz`` (reload with :meth:`load`)."""
        np.savez_compressed(
            path,
            quarter=self.quarter,
            channel=self.channel,
            targets=_plain(self.targets),
            cadenceno=self.cadenceno,
            time=self.time,
            present=self.present,
            columns=np.array(list(self.data)),
            **{f"data_{name}": value for name, value in self.data.items()},
        )

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as f:
            return cls(
                quarter=int(f["quarter"]),
                channel=int(f["channel"]),
                targets=Table(f["targets"]),
                cadenceno=f["cadenceno"],
                time=f["time"],
                present=f["present"],
                data={str(name): f[f"data_{name}"] for name in f["columns"]},
            )


def _plain(table):
    """Structured array with masked values filled (NaN, -1 or ""), for np.savez."""
    columns = {}
    for name in table.colnames:
        column = table[name]
        if hasattr(column, "mask") and np.any(column.mask):
            fill = np.nan if column.dtype.kind == "f" else ("" if column.dtype.kind in "US" else -1)
            if column.dtype.kind in "iu" and fill == -1:
                column = column.astype(np.int64)
            columns[name] = column.filled(fill)
        else:
            columns[name] = np.asarray(column)
    return Table(columns).as_array()


def read_table_at(path, offset, size, extname):
    """Read one binary table HDU at a known byte offset, as a native-order array."""
    fd = os.open(path, os.O_RDONLY)
    try:
        buf = os.pread(fd, size, offset)
    finally:
        os.close(fd)
    hdu = fits.BinTableHDU.fromstring(buf)
    if hdu.name != extname:
        raise ValueError(f"{path}: expected {extname} at byte {offset}, found {hdu.name}")
    return np.array(hdu.data, dtype=hdu.data.dtype.newbyteorder("="))


def read_channel(
    quarter,
    channel,
    root=None,
    columns=DEFAULT_COLUMNS,
    dtype=None,
    workers=16,
    index=None,
    progress=False,
    **filters,
):
    """Read the light curves of every star on ``channel`` in ``quarter``.

    ``filters`` select stars and are passed to :meth:`Index.query`, e.g.
    ``kepmag=(None, 15)``, ``crowdsap=(0.9, None)``, ``limit=500``.
    ``columns`` are per-quarter table columns to read; ``dtype`` (e.g.
    ``np.float32``) converts the float columns to save memory.  ``index`` is an
    open :class:`Index` to reuse.
    """
    root = data_dir(root)
    columns = tuple(c.upper() for c in columns)
    extname = f"LIGHTCURVE_Q{quarter}"
    if index is None:
        with Index(root) as idx:
            targets = idx.query(quarter=quarter, channel=channel, **filters)
    else:
        targets = index.query(quarter=quarter, channel=channel, **filters)

    keep = ("CADENCENO", "TIME", *columns)

    def read(row):
        path = root / row["path"]
        try:
            table = read_table_at(path, int(row["hdu_offset"]), int(row["hdu_size"]), extname)
        except (OSError, ValueError):  # file replaced since indexing: fall back to a full walk
            table = read_quarter(path, quarter).data
        return {name: table[name] for name in dict.fromkeys(keep)}

    start = time.monotonic()
    tables = []
    with ThreadPoolExecutor(workers) as pool:
        for i, table in enumerate(pool.map(read, targets), 1):
            tables.append(table)
            if progress and (i % 500 == 0 or i == len(targets)):
                rate = i / (time.monotonic() - start)
                print(f"{i}/{len(targets)} light curves  {rate:.1f}/s", flush=True)

    if tables:
        cadenceno = np.unique(np.concatenate([t["CADENCENO"] for t in tables]))
    else:
        cadenceno = np.array([], dtype=np.int32)
    shape = (len(tables), len(cadenceno))
    present = np.zeros(shape, dtype=bool)
    times = np.full(shape, np.nan)
    data = {}
    for name in columns:
        kind = tables[0][name].dtype if tables else np.dtype(float)
        if kind.kind == "f":
            data[name] = np.full(shape, np.nan, dtype=dtype or kind)
        else:
            data[name] = np.zeros(shape, dtype=kind)
    for i, table in enumerate(tables):
        j = np.searchsorted(cadenceno, table["CADENCENO"])
        present[i, j] = True
        times[i, j] = table["TIME"]
        for name in columns:
            data[name][i, j] = table[name]
    with np.errstate(all="ignore"):
        median_time = np.nanmedian(times, axis=0) if len(tables) else np.array([])
    return ChannelQuarter(
        quarter=int(quarter),
        channel=int(channel),
        targets=targets,
        cadenceno=cadenceno,
        time=median_time,
        present=present,
        data=data,
    )
