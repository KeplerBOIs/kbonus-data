"""Read and download KBonus-BKG Kepler light curves."""

from .catalog import CATALOGS, catalog_path, lc_filenames, load_catalog
from .lightcurve import (
    BKJD_OFFSET,
    LightCurve,
    QuarterData,
    available_quarters,
    read_lightcurve,
    read_quarter,
)
from .paths import data_dir, find_lc, lc_filename, lc_url, local_lc_path, normalize_target

__version__ = "0.1.0"

__all__ = [
    "BKJD_OFFSET",
    "CATALOGS",
    "LightCurve",
    "QuarterData",
    "available_quarters",
    "catalog_path",
    "data_dir",
    "find_lc",
    "lc_filename",
    "lc_filenames",
    "lc_url",
    "load_catalog",
    "local_lc_path",
    "normalize_target",
    "read_lightcurve",
    "read_quarter",
]
