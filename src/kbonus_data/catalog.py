"""KBonus-BKG catalogs: the full source catalog and the three sub-catalogs."""

import warnings

from astropy.io import fits
from astropy.table import Table
from astropy.units import UnitsWarning

from .paths import ARCHIVE_BASE, data_dir, lc_filename

# name -> archive file name.  "source" lists every light curve.
CATALOGS = {
    "source": "hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits",
    "kois-nns": "hlsp_kbonus-bkg_kepler_kepler_kois-nns-catalog_kepler_v1.0_cat.fits",
    "mstars": "hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits",
    "wd": "hlsp_kbonus-bkg_kepler_kepler_wd-catalog_kepler_v1.0_cat.fits",
}


def _check(name):
    if name not in CATALOGS:
        raise ValueError(f"catalog must be one of {sorted(CATALOGS)}, not {name!r}")


def catalog_url(name, base=ARCHIVE_BASE):
    _check(name)
    return f"{base}/catalogs/{CATALOGS[name]}"


def catalog_path(name, root=None):
    _check(name)
    return data_dir(root) / "catalogs" / CATALOGS[name]


def load_catalog(name="source", root=None):
    """Load a catalog as an astropy Table."""
    path = catalog_path(name, root)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Fetch it with: kbonus-data download catalogs")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UnitsWarning)  # 'e-/s' is not a FITS unit
        return Table.read(path, hdu=1)


def lc_filenames(name="source", root=None):
    """Return the light curve file name of each catalog row.

    Rows are not unique in the kois-nns catalog: a source that neighbours two
    KOIs is listed twice.
    """
    path = catalog_path(name, root)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Fetch it with: kbonus-data download catalogs")
    with fits.open(path) as hdul:
        # the fname column holds the target, e.g. "kic-003357387"
        return [lc_filename(str(value).strip()) for value in hdul[1].data["fname"]]
