"""Read and download KBonus-BKG Kepler light curves."""

from .catalog import CATALOGS, catalog_path, lc_filenames, load_catalog
from .channels import ChannelQuarter, read_channel
from .focalplane import CHANNELS, MODULES, QUARTERS, channel, module_output, season
from .index import Index, build_index, index_path
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
    "CHANNELS",
    "ChannelQuarter",
    "Index",
    "LightCurve",
    "MODULES",
    "QUARTERS",
    "QuarterData",
    "available_quarters",
    "build_index",
    "catalog_path",
    "channel",
    "data_dir",
    "find_lc",
    "index_path",
    "lc_filename",
    "lc_filenames",
    "lc_url",
    "load_catalog",
    "local_lc_path",
    "module_output",
    "normalize_target",
    "read_channel",
    "read_lightcurve",
    "read_quarter",
    "season",
]
