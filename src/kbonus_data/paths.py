"""Target naming, local storage layout and archive URLs for KBonus-BKG.

Local layout under the data directory::

    catalogs/<catalog file>
    lcs/<prefix>/<light curve file>
    logs/

``<prefix>`` is the first four characters of the zero-padded KIC number or of
the Gaia DR3 source id, the same prefix MAST uses.  MAST adds one more
directory per source (``lcs/<prefix>/<id>/<file>``); that level is dropped
locally because it holds exactly one file each.
"""

import os
import re
from pathlib import Path

ARCHIVE_BASE = "https://archive.stsci.edu/hlsps/kbonus-bkg"
DATA_DIR_ENV = "KBONUS_DATA_DIR"
DEFAULT_DATA_DIR = "."  # the current directory
HLSP_VERSION = "1.0"

_FILENAME_RE = re.compile(
    r"hlsp_kbonus-bkg_kepler_kepler_(kic-\d{9}|gaia-dr3-\d+)_kepler_v([\d.]+)_lc\.fits$"
)
_KIC_RE = re.compile(r"^(?:kic)?[\s_-]*(\d{1,9})$", re.IGNORECASE)
_GAIA_RE = re.compile(r"^(?:gaia)?[\s_-]*(?:dr3)?[\s_-]*(\d{10,})$", re.IGNORECASE)


def data_dir(root=None):
    """Return the data directory: ``root``, else $KBONUS_DATA_DIR, else the default."""
    return Path(root or os.environ.get(DATA_DIR_ENV) or DEFAULT_DATA_DIR)


def normalize_target(target):
    """Return the canonical target string, ``kic-008311864`` or ``gaia-dr3-<id>``.

    Accepts a KIC number (int or str), a Gaia DR3 source id, strings such as
    ``"KIC 8311864"`` or ``"Gaia DR3 2049134501869730432"``, or a light curve
    file name.  Bare numbers above 9 digits are taken to be Gaia source ids.
    """
    text = str(target).strip()
    m = _FILENAME_RE.search(text)
    if m:
        return m.group(1)
    m = _KIC_RE.match(text)
    if m:
        return f"kic-{int(m.group(1)):09d}"
    m = _GAIA_RE.match(text)
    if m:
        return f"gaia-dr3-{m.group(1)}"
    raise ValueError(f"Cannot interpret {target!r} as a KIC or Gaia DR3 target")


def target_id(target):
    """Return the id part of a target: 9-digit KIC string or Gaia source id."""
    return normalize_target(target).rsplit("-", 1)[1]


def prefix(target):
    """Return the 4-character directory prefix for a target."""
    return target_id(target)[:4]


def lc_filename(target, version=HLSP_VERSION):
    return f"hlsp_kbonus-bkg_kepler_kepler_{normalize_target(target)}_kepler_v{version}_lc.fits"


def lc_url(target, version=HLSP_VERSION, base=ARCHIVE_BASE):
    """Return the archive URL of a light curve.  ``target`` may be a file name."""
    m = _FILENAME_RE.search(str(target))
    if m:
        version = m.group(2)
    name = lc_filename(target, version)
    return f"{base}/lcs/{prefix(name)}/{target_id(name)}/{name}"


def local_lc_path(target, root=None, version=HLSP_VERSION):
    """Return where a light curve is stored locally (whether or not it exists)."""
    m = _FILENAME_RE.search(str(target))
    if m:
        version = m.group(2)
    name = lc_filename(target, version)
    return data_dir(root) / "lcs" / prefix(name) / name


def _version_key(path):
    return tuple(int(p) for p in _FILENAME_RE.search(path.name).group(2).split("."))


def find_lc(target, root=None):
    """Return the path of a downloaded light curve, the highest version present.

    ``target`` may also be a path to an existing file, which is returned as is.
    """
    if isinstance(target, (str, os.PathLike)) and os.path.isfile(target):
        return Path(target)
    directory = data_dir(root) / "lcs" / prefix(target)
    matches = [p for p in directory.glob(lc_filename(target, "*")) if _FILENAME_RE.search(p.name)]
    if not matches:
        raise FileNotFoundError(
            f"No light curve for {normalize_target(target)} in {directory}. "
            f"Fetch it with: kbonus-data download lcs --targets {normalize_target(target)}"
        )
    return max(matches, key=_version_key)
