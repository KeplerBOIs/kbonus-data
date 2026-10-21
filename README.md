# kbonus-data

Read and download [KBonus-BKG](https://archive.stsci.edu/hlsp/kbonus-bkg) Kepler
light curves (Martinez-Palomera, Hedges & Dotson 2023,
[doi:10.17909/7jbr-w430](https://dx.doi.org/10.17909/7jbr-w430)). This is the
shared data-access library for the [KeplerBOIs](https://github.com/KeplerBOIs)
repos.

## Install

For work on the KeplerBOIs project, set up the shared workspace described in
the [KeplerBOIs workspace guide](https://github.com/KeplerBOIs/.github#readme):
the repos are cloned side by side and installed, editable, into one virtual
environment, so changes here are seen at once by every other repo.

To install only this package:

```sh
pip install git+https://github.com/KeplerBOIs/kbonus-data
```

or, to work on it on its own:

```sh
git clone https://github.com/KeplerBOIs/kbonus-data
pip install -e 'kbonus-data[test]'
```

The library needs only `numpy`, `astropy` and `requests`. The example
notebook (`notebooks/`) also uses `matplotlib` and
[`pytfit5`](https://github.com/jasonfrowe/bls_cuda) for the BLS search; both
are in the workspace environment.

Other KeplerBOIs repos list it in their `pyproject.toml` by name only:

```toml
dependencies = ["kbonus-data"]
```

not as `kbonus-data @ git+https://...`. In the workspace, pip then uses the
local editable copy; a direct URL would make pip replace it with a frozen copy
from GitHub. kbonus-data is not on PyPI, so to install such a repo on its own,
install kbonus-data from GitHub first.

## Data directory

Everything lives in one data directory. Our copy is at
`/opt/data2/rowe/kbonus`, which the examples use; put yours wherever there is
room (4.8 TB for the full set) and substitute its path.

### Layout

```
/opt/data2/rowe/kbonus/              the data directory
├── catalogs/                        4 catalog FITS files, 165 MB
│   ├── hlsp_kbonus-bkg_kepler_kepler_source-catalog_kepler_v1.0_cat.fits     every light curve, 606,900 rows
│   ├── hlsp_kbonus-bkg_kepler_kepler_kois-nns-catalog_kepler_v1.0_cat.fits   KOIs and their neighbours
│   ├── hlsp_kbonus-bkg_kepler_kepler_mstars-catalog_kepler_v1.0_cat.fits     M dwarfs
│   └── hlsp_kbonus-bkg_kepler_kepler_wd-catalog_kepler_v1.0_cat.fits         white dwarfs
├── index/
│   └── kbonus_index.sqlite          quarter/channel index, 1.2 GB
├── lcs/                             light curves, 606,900 files, 4.8 TB
│   ├── 0007/ … 0129/                120 directories of KIC sources (541,399)
│   │   └── hlsp_kbonus-bkg_kepler_kepler_kic-008300021_kepler_v1.0_lc.fits
│   └── 2049/ … 2143/                44 directories of Gaia-only sources (65,501)
│       └── hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2049134501869730432_kepler_v1.0_lc.fits
└── logs/
    ├── download.log                 progress of the download (from the tee in the command below)
    ├── failed_<time>.txt            files a download run could not fetch
    └── index_build.log              progress of the index build
```

Each light curve sits in `lcs/<prefix>/`, where `<prefix>` is the first four
characters of the 9-digit zero-padded KIC number (`0083` for KIC 8311864,
`008311864`) or of the Gaia DR3 source id (`2049`). This is the prefix MAST
uses; it keeps directories under ~10,000 files. MAST's extra per-source level
(`lcs/0083/008311864/<file>`) is dropped because each holds a single file.

What needs which directory:

| directory | written by | needed for |
|---|---|---|
| `catalogs/` | `kbonus-data download catalogs` | `load_catalog`, and the list of files to download |
| `lcs/` | `kbonus-data download lcs` | reading any light curve |
| `index/` | `kbonus-data index build` | `read_channel` and `Index` (quarter/channel selection) |
| `logs/` | the download and index commands | nothing; for your information |

Other files and directories in the data directory are ignored (ours has a
`LCs_v111/` holding one test file), as are `.part` files left by an
interrupted download.

### Telling the package where it is

The data directory is, in order of precedence:

1. the `root=` argument of any function, or `--root` on the command line;
2. the `KBONUS_DATA_DIR` environment variable;
3. the current directory.

The simplest is to set `KBONUS_DATA_DIR` at the top of a script or notebook.
It is read each time a function needs it, so set it before the first call
(before or after `import kbonus_data`):

```python
import os
os.environ["KBONUS_DATA_DIR"] = "/opt/data2/rowe/kbonus"   # your data directory

import kbonus_data as kb
lc = kb.read_lightcurve(8311864)             # reads /opt/data2/rowe/kbonus/lcs/0083/...
```

or for the shell, e.g. in `~/.bashrc` or `~/.zshrc`:

```sh
export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus
```

or for one call: `kb.read_lightcurve(8311864, root="/opt/data2/rowe/kbonus")`,
`kbonus-data --root /opt/data2/rowe/kbonus status`.

Nothing in the package or the index refers to a particular machine. The
index stores paths relative to the data directory, so the directory can be
moved or copied whole and used from its new location.

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

For [`pytfit5`](https://github.com/jasonfrowe/bls_cuda) (bls_cuda):

```python
import pytfit5.transitPy5 as tpy5
phot = lc.to_phot(kind="sap", phot=tpy5.phot_class())   # time = BJD - 2454900
```

Catalogs:

```python
cat = kb.load_catalog("source")              # astropy Table, 606,900 rows
kb.load_catalog("kois-nns")                  # also "mstars", "wd"
```

Sources are named by KIC number when they have a KIC match and by Gaia DR3 id
otherwise; the catalog's `fname` column gives the name for each row.

## Light curves by quarter and channel

Where a star fell on the detector (module, output, channel) changes every
quarter and is recorded only in each file's per-quarter headers. An SQLite
index of those headers lets you select every star on a channel in a quarter
and read them as one array, e.g. to build systematics basis vectors.

```python
import os
os.environ["KBONUS_DATA_DIR"] = "/opt/data2/rowe/kbonus"   # your data directory

import kbonus_data as kb

cq = kb.read_channel(5, 41)                  # Q5, channel 41 (module 13, output 1)
cq = kb.read_channel(5, 41, kepmag=(None, 15), crowdsap=(0.9, None), limit=1000)

cq.flux, cq.flux_err                         # PSF flux, (n_stars, n_cadences)
cq.sap_flux, cq.sap_bkg, cq.centroid_row     # any per-quarter column, see columns=
cq.sap_quality                               # integer; 0 where a star has no data
cq.present                                   # False where a star has no row for a cadence (flux is NaN)
cq.cadenceno, cq.time                        # cadence grid; time is the median BKJD over stars
cq.targets                                   # one row per star: target, kic, kepmag, crowdsap, flfrcsap, ...

cq.save("q05_ch41.npz"); cq = kb.ChannelQuarter.load("q05_ch41.npz")
```

The cadence grid holds every cadence at which any star on the channel has
data. Rows are brightest first. Filters take a value, a list or a
`(low, high)` range; `kepmag` drops sources with no KIC match, which are
ordered by Gaia G. `dtype=np.float32` halves memory. A channel holds ~1,500
to ~16,500 stars in a quarter (median ~3,900) and reads at ~20 stars/s from the
NFS data directory, so about 3 minutes for a typical channel. Module 3
(channels 5–8) failed during Q4 and has no data from then on.

The index can also be queried directly:

```python
idx = kb.Index()                                         # index/kbonus_index.sqlite in the data directory
idx.query(quarter=5, channel=41, crowdsap=(0.95, None))   # astropy Table, one row per star-quarter
idx.query(quarter=5, module=13)                          # all four outputs of a module
idx.quarters(8311864)                                    # channel of one star in every quarter
idx.channel_counts(5)                                    # stars per channel
idx.sql("SELECT channel, AVG(crowdsap) FROM observations WHERE quarter = 5 GROUP BY channel")

kb.channel(13, 1), kb.module_output(41), kb.season(5)    # 41, (13, 1), 3
```

From the shell:

```sh
export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus           # or give --root on each command
kbonus-data index status --quarter 5                     # stars per channel
kbonus-data channel 5 41 --kepmag-max 15 --out q05_ch41.npz
```

### Where the index lives

The index stores each file's path relative to the data directory
(`lcs/0083/hlsp_..._lc.fits`), never an absolute path. You give the data
location when you open it, so the same index file works wherever the light
curves are:

```python
# default: the index inside the data directory, <root>/index/kbonus_index.sqlite
idx = kb.Index(root="/opt/data2/rowe/kbonus")

# the index kept somewhere else, e.g. a copy on fast local disk
idx = kb.Index(root="/opt/data2/rowe/kbonus", path="/scratch/kbonus_index.sqlite")
cq = kb.read_channel(5, 41, index=idx)       # light curves are read from idx.root
```

```sh
kbonus-data --root /opt/data2/rowe/kbonus channel 5 41 --db /scratch/kbonus_index.sqlite
```

To move the data, copy the whole directory, index included, and point
`$KBONUS_DATA_DIR` or `root=` at the new location; nothing needs rebuilding.
Someone who has the light curves but not the index can copy just
`kbonus_index.sqlite` (1.2 GB) instead of spending 11 hours rebuilding it. It
is only valid for the same set of files, which is all of KBonus-BKG v1.0.

Our index can be downloaded from the
[KeplerBOIs data page](https://kona.ubishops.ca/KeplerBOIs/) (534 MB
gzipped). Querying it with `Index` needs no light curves:

```sh
mkdir -p "$KBONUS_DATA_DIR/index"
curl -o "$KBONUS_DATA_DIR/index/kbonus_index.sqlite.gz" \
     https://kona.ubishops.ca/KeplerBOIs/index/kbonus_index.sqlite.gz
gunzip "$KBONUS_DATA_DIR/index/kbonus_index.sqlite.gz"
```

### What the index holds

`files`, one row per source (the highest version on disk):

| column | |
|---|---|
| `target` | `kic-008311864` or `gaia-dr3-<id>` |
| `path` | file path relative to the data directory |
| `kic`, `gaia_id` | identifiers; `kic` is NULL for sources without a KIC match |
| `ra`, `dec` | degrees |
| `kepmag`, `gmag`, `bpmag`, `rpmag` | magnitudes; NULL when missing (the files use 1e20) |
| `qdetect` | detection flag per quarter, one character for each of Q0–Q17 |
| `n_quarters`, `version`, `size`, `mtime`, `tpforg` | |

`observations`, one row per source and quarter, indexed on `(quarter, channel)`:

| column | |
|---|---|
| `quarter`, `season` | season is 0–3 (`kb.season(q)`) |
| `channel`, `module`, `output` | detector placement in that quarter |
| `n_cadences` | rows in the quarter's table |
| `qdetect` | that quarter's detection flag |
| `crowdsap`, `flfrcsap`, `npixsap`, `sap` | aperture metrics: fraction of flux in the aperture from the target, fraction of the target's flux captured, pixel count, aperture mode |
| `psffrac`, `pertrati`, `pertstd`, `iterneg` | PSF-model metrics from the per-quarter header |
| `row`, `col`, `tpf_row`, `tpf_col` | mean SAP centroid and origin of the pixel cutout |
| `hdu_offset`, `hdu_size`, `aperture_offset`, `aperture_size` | where the quarter's table and aperture sit in the file |

`Index.query` joins the two and returns an astropy Table; `Index.sql` runs
any query. The `errors` table lists files the last build could not read.

Each quarter's table (`columns=` in `read_channel`) has `CADENCENO`, `TIME`,
`FLUX`, `FLUX_ERR` (PSF photometry), `SAP_FLUX`, `SAP_FLUX_ERR`,
`CENTROID_COLUMN`, `CENTROID_ROW`, `SAP_QUALITY`, `SAP_BKG`, `RED_CHI2`,
`PSF_FLUX_NOVA` and `PSF_FLUX_ERR_NOVA`.

The [preconditioning](https://github.com/KeplerBOIs/preconditioning) repo
uses `read_channel` to build detrending basis vectors per quarter and
channel; its README shows the workflow.

### Building the index

```sh
export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus
mkdir -p $KBONUS_DATA_DIR/logs
tmux new -s kbonus-index "kbonus-data --root $KBONUS_DATA_DIR index build --workers 12 2>&1 | tee -a $KBONUS_DATA_DIR/logs/index_build.log"
```

The index is `index/kbonus_index.sqlite` in the data directory (1.2 GB;
606,900 sources, 6.5 million source-quarters). The
build reads only header blocks, ~20 small reads per file, plus the byte offset
of each per-quarter extension so that `read_channel` needs one read per star.
On our NFS-mounted copy it is limited by the storage server's random reads
to ~16 files/s, 11 hours for the full set; more workers do not help. It
is resumable and scans only files not yet indexed, so re-run it after
downloading more light curves. A newer version of a file replaces the old
entry. Files that could not be read are listed in the `errors` table and
retried on the next run.

## Downloading

The full set is 606,900 files, 4.8 TB. There is no bundle for the whole set and
no AWS mirror, so files are fetched one by one from MAST.

Files go into the data directory (see [Layout](#layout)). Set it before
downloading; our copy was downloaded to `/opt/data2/rowe/kbonus`:

```sh
export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus       # where the data will go
kbonus-data download catalogs                      # creates $KBONUS_DATA_DIR/catalogs/
```

`kbonus-data --root /opt/data2/rowe/kbonus download catalogs` does the same
without the variable. With neither, files go into the current directory
(`./catalogs/`, `./lcs/`, `./logs/`).

The commands (all take `--root`):

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
mkdir -p $KBONUS_DATA_DIR/logs
tmux new -s kbonus "kbonus-data --root $KBONUS_DATA_DIR download lcs --workers 16 2>&1 | tee -a $KBONUS_DATA_DIR/logs/download.log"
```

## Tests

```sh
pytest
```
