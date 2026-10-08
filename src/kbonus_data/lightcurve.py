"""Read KBonus-BKG light curve FITS files."""

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.io import fits

from .paths import find_lc, normalize_target

STITCHED_EXT = "LIGHTCURVE_STITCHED"
BKJD_OFFSET = 2454833.0  # TIME column is BJD - 2454833
LONG_CADENCE_DAYS = 1765.4635 / 86400.0  # Kepler long cadence integration time

# kind -> (flux column, error column) in the stitched table
FLUX_KINDS = {
    "sap": ("sap_flux", "sap_flux_err"),
    "psf": ("flux", "flux_err"),
    "psf_flat": ("psf_flat_flux", "psf_flat_flux_err"),
}


def _native(column):
    """Copy a FITS column into a native-byte-order array (FITS is big-endian)."""
    return np.array(column, dtype=column.dtype.newbyteorder("="))


@dataclass
class LightCurve:
    """The stitched (all quarters) light curve of one source.

    ``time`` is BKJD (BJD - 2454833).  ``flux`` is the PSF photometry,
    ``sap_flux`` the aperture photometry and ``psf_flat_flux`` the flattened
    PSF photometry, all in e-/s.
    """

    path: Path
    header: fits.Header
    time: np.ndarray
    cadenceno: np.ndarray
    quarter: np.ndarray
    flux: np.ndarray
    flux_err: np.ndarray
    sap_flux: np.ndarray
    sap_flux_err: np.ndarray
    psf_flat_flux: np.ndarray
    psf_flat_flux_err: np.ndarray
    sap_quality: np.ndarray
    flatten_mask: np.ndarray

    def __len__(self):
        return len(self.time)

    @property
    def kic(self):
        """KIC number, or None for sources without a KIC match."""
        kic = self.header.get("KEPLERID")
        return int(kic) if kic else None

    @property
    def gaia_id(self):
        return self.header.get("GAIAID")

    @property
    def label(self):
        return self.header.get("LABEL")

    @property
    def ra(self):
        return self.header.get("RA_OBJ")

    @property
    def dec(self):
        return self.header.get("DEC_OBJ")

    @property
    def bjd(self):
        return self.time + BKJD_OFFSET

    def get_flux(self, kind="sap"):
        """Return ``(flux, flux_err)`` for ``kind`` in 'sap', 'psf', 'psf_flat'."""
        try:
            flux_name, err_name = FLUX_KINDS[kind]
        except KeyError:
            raise ValueError(f"kind must be one of {sorted(FLUX_KINDS)}, not {kind!r}") from None
        return getattr(self, flux_name), getattr(self, err_name)

    def to_phot(self, kind="sap", time_zero=2454900.0, itime=LONG_CADENCE_DAYS, phot=None):
        """Fill a bls_cuda ``transitPy5.phot_class``-style object.

        ``phot`` is the object to fill (e.g. ``tpy5.phot_class()``); if omitted
        a plain namespace with the same attributes is returned, so bls_cuda is
        not required.  Times are BJD - ``time_zero``.
        """
        flux, ferr = self.get_flux(kind)
        if phot is None:
            phot = SimpleNamespace()
        n = len(self)
        phot.time = self.bjd - time_zero
        phot.flux = flux.copy()
        phot.ferr = ferr.copy()
        phot.itime = np.full(n, itime)
        phot.tflag = np.zeros(n)  # 1 = in transit
        phot.flux_f = flux.copy()  # placeholder for detrended flux
        phot.icut = np.zeros(n)  # 0 = keep, 1 = toss
        return phot


@dataclass
class QuarterData:
    """One quarter of a light curve: table, per-quarter header and SAP aperture."""

    quarter: int
    data: np.ndarray  # structured array, columns as in the LIGHTCURVE_Q<n> table
    header: fits.Header
    aperture: np.ndarray | None  # pixel mask in the shape of the originating TPF


def read_lightcurve(target, root=None):
    """Read the stitched light curve of ``target``.

    ``target`` is a path to a FITS file, or a KIC number / Gaia DR3 id / target
    string that is looked up in the data directory (``root``, else
    $KBONUS_DATA_DIR).
    """
    path = find_lc(target, root)
    with fits.open(path) as hdul:
        hdu = hdul[STITCHED_EXT] if STITCHED_EXT in hdul else hdul[1]
        data = hdu.data
        columns = {
            name: _native(data[name.upper()])
            for name in (
                "time",
                "cadenceno",
                "quarter",
                "flux",
                "flux_err",
                "sap_flux",
                "sap_flux_err",
                "psf_flat_flux",
                "psf_flat_flux_err",
                "sap_quality",
                "flatten_mask",
            )
        }
        return LightCurve(path=path, header=hdul[0].header.copy(), **columns)


def available_quarters(target, root=None):
    """Return the quarters for which ``target`` has a per-quarter extension."""
    with fits.open(find_lc(target, root)) as hdul:
        return sorted(
            int(hdu.name.removeprefix("LIGHTCURVE_Q"))
            for hdu in hdul
            if hdu.name.startswith("LIGHTCURVE_Q")
        )


def read_quarter(target, quarter, root=None):
    """Read the per-quarter table, header and aperture mask for one quarter."""
    with fits.open(find_lc(target, root)) as hdul:
        try:
            hdu = hdul[f"LIGHTCURVE_Q{quarter}"]
        except KeyError:
            raise KeyError(f"{normalize_target(hdul.filename())} has no quarter {quarter}") from None
        dtype = hdu.data.dtype.newbyteorder("=")
        data = np.array(hdu.data, dtype=dtype)
        aperture_name = f"APERTURE_Q{quarter}"
        aperture = hdul[aperture_name].data.copy() if aperture_name in hdul else None
        return QuarterData(quarter=int(quarter), data=data, header=hdu.header.copy(), aperture=aperture)
