import numpy as np
import pytest
from astropy.io import fits

import kbonus_data as kb
from kbonus_data.index import scan_file, walk_hdus

QUARTER_COLUMNS = [
    ("CADENCENO", "J"),
    ("TIME", "D"),
    ("FLUX", "D"),
    ("FLUX_ERR", "D"),
    ("SAP_FLUX", "D"),
    ("SAP_FLUX_ERR", "D"),
    ("CENTROID_COLUMN", "D"),
    ("CENTROID_ROW", "D"),
    ("SAP_QUALITY", "J"),
    ("SAP_BKG", "D"),
    ("RED_CHI2", "D"),
    ("PSF_FLUX_NOVA", "D"),
    ("PSF_FLUX_ERR_NOVA", "D"),
]


def write_lc(root, target, placements, kepmag=14.0, version="1.0", cadences=None):
    """Write a KBonus-like file.  ``placements`` maps quarter -> channel.

    Flux in quarter q of the star is ``kepmag * 1000 + cadenceno``, so reads
    can be checked.  ``cadences`` maps quarter -> cadence numbers (default
    ``100 * q + arange(10)``).
    """
    primary = fits.PrimaryHDU()
    if target.startswith("kic"):
        primary.header["KEPLERID"] = int(target.split("-")[1])
    primary.header["GAIAID"] = "Gaia DR3 2049134501869730432"
    primary.header["RA_OBJ"] = 290.0
    primary.header["DEC_OBJ"] = 44.0
    primary.header["KEPMAG"] = kepmag
    primary.header["GMAG"] = kepmag + 0.1 if kepmag < 1e19 else 17.0
    primary.header["QDETECT"] = "".join("1" if q in placements else "0" for q in range(18))
    stitched = fits.BinTableHDU.from_columns(
        [fits.Column(name="TIME", format="D", array=np.arange(1000.0))], name="LIGHTCURVE_STITCHED"
    )
    hdus = [primary, stitched]
    for quarter, channel in sorted(placements.items()):
        cad = np.asarray(cadences[quarter] if cadences else 100 * quarter + np.arange(10))
        values = {name: np.zeros(len(cad)) for name, _ in QUARTER_COLUMNS}
        values["CADENCENO"] = cad
        values["TIME"] = 130.0 + cad * 0.0204
        values["FLUX"] = kepmag * 1000 + cad
        values["SAP_FLUX"] = kepmag * 1000 + cad + 0.5
        values["SAP_QUALITY"] = cad % 3
        table = fits.BinTableHDU.from_columns(
            [fits.Column(name=n, format=f, array=values[n]) for n, f in QUARTER_COLUMNS],
            name=f"LIGHTCURVE_Q{quarter}",
        )
        module, output = kb.module_output(channel)
        for key, value in [
            ("SEASON", kb.season(quarter)),
            ("MODULE", module),
            ("OUTPUT", output),
            ("CHANNEL", channel),
            ("QUARTER", quarter),
            ("CROWDSAP", 0.95),
            ("FLFRCSAP", 0.9),
            ("NPIXSAP", 4),
            ("SAP", "optimal"),
        ]:
            table.header[key] = value
        hdus += [table, fits.ImageHDU(np.ones((3, 4), dtype=np.uint8), name=f"APERTURE_Q{quarter}")]
    path = kb.local_lc_path(target, root, version=version)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


@pytest.fixture
def root(tmp_path):
    write_lc(tmp_path, "kic-008311864", {2: 41, 3: 42, 6: 41}, kepmag=12.0)
    write_lc(tmp_path, "kic-000757039", {2: 41, 3: 5}, kepmag=15.0)
    write_lc(tmp_path, "kic-000757040", {2: 41}, kepmag=16.5, cadences={2: [202, 203, 210, 211]})
    write_lc(tmp_path, "gaia-dr3-2049134501869730432", {2: 7}, kepmag=1e20)
    return tmp_path


def test_focalplane():
    assert kb.channel(2, 1) == 1 and kb.channel(24, 4) == 84
    assert kb.channel(13, 1) == 41 and kb.module_output(41) == (13, 1)
    assert all(kb.channel(*kb.module_output(c)) == c for c in kb.CHANNELS)
    assert [kb.season(q) for q in range(8)] == [3, 3, 0, 1, 2, 3, 0, 1]
    with pytest.raises(ValueError):
        kb.channel(5, 1)


def test_walk_hdus_matches_astropy(root):
    path = kb.find_lc(8311864, root)
    walked = list(walk_hdus(path))
    with fits.open(path) as hdul:
        assert [h.get("EXTNAME", "PRIMARY") for h, _, _ in walked] == [h.name for h in hdul]
        for (_, offset, size), hdu in zip(walked, hdul):
            info = hdu.fileinfo()
            assert offset == info["hdrLoc"]
            assert offset + size >= info["datLoc"] + info["datSpan"]


def test_scan_file(root):
    row, observations = scan_file(kb.find_lc(8311864, root))
    assert row["target"] == "kic-008311864" and row["kic"] == 8311864
    assert row["gaia_id"] == 2049134501869730432 and row["kepmag"] == 12.0
    assert [o["quarter"] for o in observations] == [2, 3, 6]
    assert [o["channel"] for o in observations] == [41, 42, 41]
    assert observations[0]["module"] == 13 and observations[0]["n_cadences"] == 10
    assert observations[0]["qdetect"] == 1 and observations[0]["aperture_offset"] > 0
    gaia_row, _ = scan_file(kb.find_lc("gaia-dr3-2049134501869730432", root))
    assert gaia_row["kic"] is None and gaia_row["kepmag"] is None  # 1e20 = missing


def test_build_and_query(root, capsys):
    assert kb.build_index(root, workers=2) == 4
    assert kb.build_index(root, workers=2) == 0  # resumable: nothing new
    with kb.Index(root) as idx:
        assert idx.summary()["files"] == 4 and idx.summary()["observations"] == 7
        stars = idx.query(quarter=2, channel=41)
        assert list(stars["target"]) == ["kic-008311864", "kic-000757039", "kic-000757040"]  # brightest first
        assert list(idx.query(quarter=2, channel=41, kepmag=(13, None))["kic"]) == [757039, 757040]
        assert len(idx.query(quarter=2, channel=[41, 7])) == 4
        assert list(idx.query(quarter=2, order_by="COALESCE(f.kepmag, f.gmag)")["target"])[-1].startswith("gaia")
        assert list(idx.quarters(8311864)["channel"]) == [41, 42, 41]
        counts = idx.channel_counts(2)
        assert dict(zip(counts["channel"], counts["n"])) == {7: 1, 41: 3}
        with pytest.raises(TypeError):
            idx.query(bogus=1)


def test_new_version_replaces_old(root):
    kb.build_index(root, workers=2)
    write_lc(root, "kic-000757039", {4: 60}, kepmag=15.0, version="1.1")
    assert kb.build_index(root, workers=2) == 1
    with kb.Index(root) as idx:
        rows = idx.quarters(757039)
        assert list(rows["quarter"]) == [4] and rows["path"][0].endswith("_v1.1_lc.fits")
        assert idx.summary()["files"] == 4


def test_read_channel(root, tmp_path):
    kb.build_index(root, workers=2)
    cq = kb.read_channel(2, 41, root, workers=2)
    assert len(cq) == 3 and (cq.module, cq.output) == (13, 1)
    np.testing.assert_array_equal(cq.cadenceno, np.arange(200, 212))
    assert cq.shape == (3, 12) and cq.flux.shape == (3, 12)
    np.testing.assert_array_equal(cq.flux[0, :10], 12000 + np.arange(200, 210))
    assert np.all(np.isnan(cq.flux[0, 10:])) and not cq.present[0, 10:].any()
    third = cq.present[2]
    np.testing.assert_array_equal(cq.cadenceno[third], [202, 203, 210, 211])
    np.testing.assert_array_equal(cq.sap_flux[2, third], 16500.5 + cq.cadenceno[third])
    assert cq.sap_quality.dtype.kind == "i" and cq.sap_quality[2, 0] == 0
    np.testing.assert_allclose(cq.time, 130.0 + cq.cadenceno * 0.0204)
    flux, err = cq.get_flux("psf")
    assert flux is cq.data["FLUX"]

    bright = kb.read_channel(2, 41, root, kepmag=(None, 13), columns=["SAP_FLUX"], dtype=np.float32)
    assert len(bright) == 1 and bright.sap_flux.dtype == np.float32
    with pytest.raises(KeyError):
        bright.get_flux("psf")

    out = tmp_path / "cq.npz"
    cq.save(out)
    back = kb.ChannelQuarter.load(out)
    assert (back.quarter, back.channel) == (2, 41)
    assert list(back.targets["target"]) == list(cq.targets["target"])
    np.testing.assert_array_equal(back.flux, cq.flux)
    np.testing.assert_array_equal(back.present, cq.present)

    assert len(kb.read_channel(2, 80, root)) == 0


def test_cli(root, tmp_path, capsys):
    from kbonus_data.cli import main

    assert main(["--root", str(root), "index", "build", "--workers", "2"]) == 0
    assert main(["--root", str(root), "index", "status", "--quarter", "2"]) == 0
    assert "ch 41       3" in capsys.readouterr().out
    out = tmp_path / "x.npz"
    assert main(["--root", str(root), "channel", "2", "41", "--out", str(out), "--kepmag-max", "15"]) == 0
    assert len(kb.ChannelQuarter.load(out)) == 2


def test_default_data_dir_is_current_directory(tmp_path, monkeypatch):
    monkeypatch.delenv("KBONUS_DATA_DIR", raising=False)
    monkeypatch.chdir(tmp_path)
    assert kb.data_dir().resolve() == tmp_path.resolve()
    assert kb.local_lc_path(8311864).resolve() == (tmp_path / "lcs" / "0083" / kb.lc_filename(8311864)).resolve()


def test_index_is_relocatable(root, tmp_path_factory, monkeypatch):
    import shutil

    tmp_path = tmp_path_factory.mktemp("other")  # outside the data directory

    monkeypatch.delenv("KBONUS_DATA_DIR", raising=False)
    kb.build_index(root, workers=2)
    with kb.Index(root) as idx:
        assert all(not p.startswith("/") for p in idx.sql("SELECT path FROM files")["path"])

    moved = tmp_path / "elsewhere" / "kbonus"
    shutil.copytree(root, moved)
    shutil.rmtree(root)  # nothing may still point at the old location
    with kb.Index(moved) as idx:
        assert len(kb.read_channel(2, 41, index=idx)) == 3  # data location comes from the index
    assert len(kb.read_channel(2, 41, moved)) == 3

    db = tmp_path / "my_index.sqlite"
    shutil.move(moved / "index" / "kbonus_index.sqlite", db)  # index kept apart from the data
    with kb.Index(root=moved, path=db) as idx:
        cq = kb.read_channel(2, 41, index=idx)
        np.testing.assert_array_equal(cq.flux[0, :10], 12000 + np.arange(200, 210))

    from kbonus_data.cli import main

    out = tmp_path / "x.npz"
    assert main(["--root", str(moved), "channel", "2", "41", "--db", str(db), "--out", str(out)]) == 0
    assert len(kb.ChannelQuarter.load(out)) == 3
