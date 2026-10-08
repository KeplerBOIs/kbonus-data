import numpy as np
import pytest
from astropy.io import fits

import kbonus_data as kb
from kbonus_data import download

KIC_NAME = "hlsp_kbonus-bkg_kepler_kepler_kic-008311864_kepler_v1.0_lc.fits"
GAIA_NAME = "hlsp_kbonus-bkg_kepler_kepler_gaia-dr3-2049134501869730432_kepler_v1.0_lc.fits"

STITCHED_COLUMNS = [
    ("TIME", "D"),
    ("CADENCENO", "J"),
    ("QUARTER", "K"),
    ("FLUX", "D"),
    ("FLUX_ERR", "D"),
    ("SAP_FLUX", "D"),
    ("SAP_FLUX_ERR", "D"),
    ("PSF_FLAT_FLUX", "D"),
    ("PSF_FLAT_FLUX_ERR", "D"),
    ("SAP_QUALITY", "J"),
    ("FLATTEN_MASK", "J"),
]


def write_lc(path, n=20):
    """Write a small file with the same HDU structure as a KBonus light curve."""
    primary = fits.PrimaryHDU()
    primary.header["KEPLERID"] = 8311864
    primary.header["GAIAID"] = "Gaia DR3 1"
    primary.header["RA_OBJ"] = 290.0
    values = {name: np.arange(n) + i for i, (name, _) in enumerate(STITCHED_COLUMNS)}
    values["TIME"] = 120.5 + np.arange(n) * 0.02
    values["QUARTER"] = np.repeat([0, 1], n // 2)
    stitched = fits.BinTableHDU.from_columns(
        [fits.Column(name=name, format=fmt, array=values[name]) for name, fmt in STITCHED_COLUMNS],
        name="LIGHTCURVE_STITCHED",
    )
    hdus = [primary, stitched]
    for q in (0, 1):
        table = fits.BinTableHDU.from_columns(
            [
                fits.Column(name="CADENCENO", format="J", array=np.arange(n // 2)),
                fits.Column(name="TIME", format="D", array=np.arange(n // 2) * 0.02),
            ],
            name=f"LIGHTCURVE_Q{q}",
        )
        table.header["QUARTER"] = q
        hdus += [table, fits.ImageHDU(np.ones((3, 4), dtype=np.uint8), name=f"APERTURE_Q{q}")]
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path)
    return values


@pytest.mark.parametrize(
    "target, expected",
    [
        (8311864, "kic-008311864"),
        ("8311864", "kic-008311864"),
        ("KIC 8311864", "kic-008311864"),
        ("kic-008311864", "kic-008311864"),
        (2049134501869730432, "gaia-dr3-2049134501869730432"),
        ("Gaia DR3 2049134501869730432", "gaia-dr3-2049134501869730432"),
        (KIC_NAME, "kic-008311864"),
        ("/some/dir/" + GAIA_NAME, "gaia-dr3-2049134501869730432"),
    ],
)
def test_normalize_target(target, expected):
    assert kb.normalize_target(target) == expected


def test_normalize_target_rejects_garbage():
    with pytest.raises(ValueError):
        kb.normalize_target("not a star")


def test_urls_and_local_paths(tmp_path):
    base = "https://archive.stsci.edu/hlsps/kbonus-bkg/lcs"
    assert kb.lc_url(8311864) == f"{base}/0083/008311864/{KIC_NAME}"
    assert kb.lc_url(GAIA_NAME) == f"{base}/2049/2049134501869730432/{GAIA_NAME}"
    assert kb.local_lc_path(8311864, tmp_path) == tmp_path / "lcs" / "0083" / KIC_NAME
    assert kb.local_lc_path(GAIA_NAME, tmp_path) == tmp_path / "lcs" / "2049" / GAIA_NAME


def test_data_dir_env(tmp_path, monkeypatch):
    monkeypatch.setenv("KBONUS_DATA_DIR", str(tmp_path))
    assert kb.data_dir() == tmp_path
    assert kb.data_dir("/elsewhere") != tmp_path


def test_find_lc_prefers_highest_version(tmp_path):
    for version in ("1.0", "1.1.1"):
        write_lc(kb.local_lc_path(8311864, tmp_path, version=version))
    assert kb.find_lc(8311864, tmp_path).name.endswith("_v1.1.1_lc.fits")
    with pytest.raises(FileNotFoundError):
        kb.find_lc(757039, tmp_path)


def test_read_lightcurve(tmp_path):
    path = kb.local_lc_path(8311864, tmp_path)
    values = write_lc(path)
    for lc in (kb.read_lightcurve(8311864, tmp_path), kb.read_lightcurve(path)):
        assert len(lc) == 20
        assert lc.kic == 8311864
        assert lc.ra == 290.0
        np.testing.assert_array_equal(lc.time, values["TIME"])
        np.testing.assert_array_equal(lc.sap_flux, values["SAP_FLUX"])
        np.testing.assert_array_equal(lc.quarter, values["QUARTER"])
        assert lc.time.dtype.isnative and lc.quarter.dtype.isnative
        np.testing.assert_array_equal(lc.bjd, values["TIME"] + 2454833.0)
        flux, err = lc.get_flux("psf_flat")
        np.testing.assert_array_equal(flux, values["PSF_FLAT_FLUX"])
        np.testing.assert_array_equal(err, values["PSF_FLAT_FLUX_ERR"])
    with pytest.raises(ValueError):
        lc.get_flux("bogus")


def test_to_phot(tmp_path):
    path = kb.local_lc_path(8311864, tmp_path)
    values = write_lc(path)
    lc = kb.read_lightcurve(path)
    phot = lc.to_phot(kind="psf")
    np.testing.assert_allclose(phot.time, values["TIME"] - 67.0)
    np.testing.assert_array_equal(phot.flux, values["FLUX"])
    np.testing.assert_array_equal(phot.ferr, values["FLUX_ERR"])
    assert phot.itime.shape == phot.tflag.shape == phot.icut.shape == phot.flux_f.shape == (20,)

    class Phot:
        pass

    mine = Phot()
    assert lc.to_phot(phot=mine) is mine
    np.testing.assert_array_equal(mine.flux, values["SAP_FLUX"])


def test_quarters(tmp_path):
    path = kb.local_lc_path(8311864, tmp_path)
    write_lc(path)
    assert kb.available_quarters(path) == [0, 1]
    q = kb.read_quarter(path, 1)
    assert q.quarter == 1 and q.header["QUARTER"] == 1
    assert q.data.dtype.names == ("CADENCENO", "TIME") and len(q.data) == 10
    assert q.aperture.shape == (3, 4)
    with pytest.raises(KeyError):
        kb.read_quarter(path, 5)


class FakeResponse:
    def __init__(self, body, status=200, length=None):
        self.body, self.status_code = body, status
        self.headers = {"Content-Length": str(len(body) if length is None else length)}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk):
        yield self.body


class FakeSession:
    def __init__(self, responses):
        self.responses, self.urls = list(responses), []

    def get(self, url, **kwargs):
        self.urls.append(url)
        return self.responses.pop(0)


def test_fetch_is_atomic_and_retries_short_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(download.time, "sleep", lambda s: None)
    dest = tmp_path / "a" / "file.fits"
    session = FakeSession([FakeResponse(b"abc", length=10), FakeResponse(b"abcdef")])
    assert download.fetch(session, "http://x/file.fits", dest) == 6
    assert dest.read_bytes() == b"abcdef"
    assert not list(dest.parent.glob("*.part"))

    truncated = FakeSession([FakeResponse(b"abc", length=10)])
    with pytest.raises(OSError):
        download.fetch(truncated, "http://x/other.fits", tmp_path / "other.fits", retries=1)
    assert not (tmp_path / "other.fits").exists()
    assert not list(tmp_path.glob("*.part"))

    with pytest.raises(FileNotFoundError):
        download.fetch(FakeSession([FakeResponse(b"", status=404)]), "http://x/no", tmp_path / "no")


def test_download_lightcurves_skips_present_and_logs_failures(tmp_path, monkeypatch):
    write_lc(kb.local_lc_path(KIC_NAME, tmp_path))
    fetched = []

    def fake_fetch(session, url, dest, retries=5):
        if "gaia" in url:
            raise OSError("boom")
        fetched.append(url)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b"x")
        return 1

    monkeypatch.setattr(download, "fetch", fake_fetch)
    other = kb.lc_filename(757039)
    failed = download.download_lightcurves([KIC_NAME, other, GAIA_NAME], tmp_path, workers=2)
    assert fetched == [kb.lc_url(757039)]
    assert [name for name, _ in failed] == [GAIA_NAME]
    assert download.downloaded_lc_names(tmp_path) == {KIC_NAME, other}
    assert GAIA_NAME in next((tmp_path / "logs").glob("failed_*.txt")).read_text()
