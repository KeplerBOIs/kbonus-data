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

All data (catalogs, light curves, logs and the index) lives in one data
directory. It is, in order of precedence:

1. the `root=` argument of any function, or `--root` on the command line;
2. `$KBONUS_DATA_DIR`;
3. the current directory.

Nothing in the package or the index refers to a particular machine, so the
directory can be put anywhere and moved later. On our server the copy is at
`/opt/data2/rowe/kbonus`, and the examples below use that location:

```sh
export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus      # in ~/.bashrc or ~/.zshrc
```

```python
import kbonus_data as kb
lc = kb.read_lightcurve(8311864, root="/opt/data2/rowe/kbonus")   # or rely on $KBONUS_DATA_DIR
```

The directory holds:

```
catalogs/                 the four catalog FITS files
lcs/<prefix>/<file>       light curves
logs/                     lists of failed downloads
index/kbonus_index.sqlite the quarter/channel index (see below)
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

## Light curves by quarter and channel

Where a star fell on the detector (module, output, channel) changes every
quarter and is recorded only in each file's per-quarter headers. An SQLite
index of those headers lets you select every star on a channel in a quarter
and read them as one array, e.g. to build systematics basis vectors.

```python
import kbonus_data as kb

cq = kb.read_channel(5, 41)                  # Q5, channel 41 (module 13, output 1), from $KBONUS_DATA_DIR
cq = kb.read_channel(5, 41, root="/opt/data2/rowe/kbonus")   # or name the data directory
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
idx = kb.Index("/opt/data2/rowe/kbonus")                 # or kb.Index() to use $KBONUS_DATA_DIR
idx.query(quarter=5, channel=41, crowdsap=(0.95, None))   # astropy Table, one row per star-quarter
idx.query(quarter=5, module=13)                          # all four outputs of a module
idx.quarters(8311864)                                    # channel of one star in every quarter
idx.channel_counts(5)                                    # stars per channel
idx.sql("SELECT channel, AVG(crowdsap) FROM observations WHERE quarter = 5 GROUP BY channel")

kb.channel(13, 1), kb.module_output(41), kb.season(5)    # 41, (13, 1), 3
```

From the shell:

```sh
kbonus-data --root /opt/data2/rowe/kbonus index status --quarter 5      # stars per channel
kbonus-data --root /opt/data2/rowe/kbonus channel 5 41 --kepmag-max 15 --out q05_ch41.npz
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
tmux new -s kbonus-index 'kbonus-data index build --workers 12 2>&1 | tee -a $KBONUS_DATA_DIR/logs/index_build.log'
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

Files go into the data directory: the current directory unless you give
`--root` or set `$KBONUS_DATA_DIR`. To download into the current directory:

```sh
mkdir kbonus && cd kbonus
kbonus-data download catalogs                    # writes ./catalogs, ./lcs, ./logs
```

Our copy was downloaded to `/opt/data2/rowe/kbonus`, giving the location
explicitly:

```sh
kbonus-data --root /opt/data2/rowe/kbonus download catalogs
# or: export KBONUS_DATA_DIR=/opt/data2/rowe/kbonus, then omit --root
```

The commands, all of which take `--root`:

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
tmux new -s kbonus 'kbonus-data --root /opt/data2/rowe/kbonus download lcs --workers 16 2>&1 | tee -a /opt/data2/rowe/kbonus/logs/download.log'
```

## Tests

```sh
pytest
```
