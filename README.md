# kbonus-data

Read and download [KBonus-BKG](https://archive.stsci.edu/hlsp/kbonus-bkg) Kepler
light curves (Martinez-Palomera, Hedges & Dotson 2023,
[doi:10.17909/7jbr-w430](https://dx.doi.org/10.17909/7jbr-w430)). This is the
shared data-access library for the [KeplerBOIs](https://github.com/KeplerBOIs)
repos.

## Install

```sh
pip install git+https://github.com/KeplerBOIs/kbonus-data
```

or, for development:

```sh
git clone --recurse-submodules https://github.com/KeplerBOIs/kbonus-data
pip install -e 'kbonus-data[test]'
```

In another repo's `pyproject.toml`:

```toml
dependencies = ["kbonus-data @ git+https://github.com/KeplerBOIs/kbonus-data"]
```

The library needs only `numpy`, `astropy` and `requests`. The
[`bls_cuda`](https://github.com/jasonfrowe/bls_cuda) submodule is used only by
the example notebook.

## Data directory

Data lives under `$KBONUS_DATA_DIR` (default `/opt/data2/rowe/kbonus`); every
function also takes a `root=` argument.

```
catalogs/                 the four catalog FITS files
lcs/<prefix>/<file>       light curves
logs/                     lists of failed downloads
```

`<prefix>` is the first four characters of the zero-padded KIC number
(`0083` for KIC 8311864) or of the Gaia DR3 source id (`2049`), the same
prefix MAST uses. That gives 164 directories of at most ~10,000 files. MAST's
extra per-source directory (`lcs/0083/008311864/<file>`) is dropped because
each one holds a single file.

## Reading light curves

```python
import kbonus_data as kb

lc = kb.read_lightcurve(8311864)             # KIC number, looked up in the data directory
lc = kb.read_lightcurve("gaia-dr3-2049134501869730432")
lc = kb.read_lightcurve("/path/to/hlsp_kbonus-bkg_..._lc.fits")

lc.time            # BKJD = BJD - 2454833
lc.sap_flux        # aperture photometry, e-/s (with lc.sap_flux_err)
lc.flux            # PSF photometry (with lc.flux_err)
lc.psf_flat_flux   # flattened PSF photometry (with lc.psf_flat_flux_err)
lc.quarter, lc.cadenceno, lc.sap_quality, lc.flatten_mask
lc.header          # primary header; also lc.kic, lc.gaia_id, lc.ra, lc.dec

flux, err = lc.get_flux("psf")               # "sap", "psf" or "psf_flat"
```

Arrays are native byte order, so they can go straight into numba or CUDA code.

Per-quarter tables (centroids, background, PSF chi-square) and aperture masks:

```python
kb.available_quarters(8311864)               # [0, 1, ..., 17]
q = kb.read_quarter(8311864, 5)
q.data["CENTROID_ROW"], q.header["CROWDSAP"], q.aperture
```

For `bls_cuda`:

```python
import transitPy5 as tpy5
phot = lc.to_phot(kind="sap", phot=tpy5.phot_class())   # time = BJD - 2454900
```

Catalogs:

```python
cat = kb.load_catalog("source")              # astropy Table, 606,900 rows
kb.load_catalog("kois-nns")                  # also "mstars", "wd"
```

Sources are named by KIC number when they have a KIC match and by Gaia DR3 id
otherwise; the catalog's `fname` column gives the name for each row.

## Downloading

The full set is 606,900 files, 4.8 TB. There is no bundle for the whole set and
no AWS mirror, so files are fetched one by one from MAST.

```sh
kbonus-data download catalogs                    # 165 MB, lists every light curve
kbonus-data download lcs --subset kois-nns       # KOIs and neighbours first (32,319 files)
kbonus-data download lcs --workers 16            # everything
kbonus-data status                               # progress per catalog
kbonus-data download lcs --targets 8311864 11904151
```

The download is resumable: a file is done when it exists under its final name.
Each file is written to `<name>.part`, checked against the server's
Content-Length and renamed, so an interrupted run leaves no truncated files.
Re-run the same command to continue or to retry failures (listed in `logs/`).

A pilot of 450 files ran at 23 MB/s with 8 workers and 33 MB/s with 16, which
puts the full download at roughly 40 to 60 hours. Run it detached, e.g.

```sh
tmux new -s kbonus 'kbonus-data download lcs --workers 16 2>&1 | tee -a kbonus_download.log'
```

## Tests

```sh
pytest
```
