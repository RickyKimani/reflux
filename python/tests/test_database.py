"""Component data and fluids (M10a step 2): a Database on a data root (D9),
its compounds as the files hold them, and fluids with their constants, kᵢⱼ
and α, all checked against the data files themselves."""

import os
import shutil
import tomllib

import numpy as np
import pytest

import reflux
from reflux import _native
from support import ROOT, isolated

DATA = ROOT / "data"
COMPONENTS = DATA / "components"
VENDOR = DATA / "vendor" / "chemsep"


def component_file(cas):
    return tomllib.loads((COMPONENTS / f"{cas}.toml").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


def test_a_database_opens_on_the_data_root(db):
    assert reflux.data_root() == str(DATA) == db.root and db.extra == []
    index = tomllib.loads((COMPONENTS / "index.toml").read_text(encoding="utf-8"))
    compounds = db.compounds()
    assert [c["cas"] for c in compounds] == sorted(index["compounds"])
    assert compounds[0] == {"cas": "100-18-5", **index["compounds"]["100-18-5"]}
    with pytest.raises(reflux.DataError, match="is not a data root: it has no components"):
        reflux.Database(ROOT / "python")


def test_a_compound_is_its_file(db):
    # The file the library writes, read back, is the file in the repository.
    for q in ("methane", "METHANE", "74-82-8", "VNWKTOKETHGBQD-UHFFFAOYSA-N", " Methane "):
        assert db.compound(q) == component_file("74-82-8")
    assert db.compound("water") == component_file("7732-18-5")
    with pytest.raises(reflux.NotFound) as e:
        db.compound("watr")
    assert str(e.value).startswith("no compound `watr` in the database. Did you mean Water "
                                   "(7732-18-5),")
    assert isinstance(e.value, LookupError) and e.value.code == _native.NOT_FOUND
    with pytest.raises(reflux.NotFound, match="^no compound `` in the database$"):
        db.compound("")
    with pytest.raises(ValueError, match="query has a NUL character"):
        db.compound("meth\0ane")
    with pytest.raises(TypeError):
        db.compound(b"methane")


def test_your_own_directory_comes_after_the_roots(db, tmp_path):
    # Methane with another Tc replaces the root's: the same CAS number.
    mine = tmp_path / "mine"
    mine.mkdir()
    methane = (COMPONENTS / "74-82-8.toml").read_text(encoding="utf-8")
    assert methane.count("tc = 190.56\n") == 1
    (mine / "74-82-8.toml").write_text(methane.replace("tc = 190.56\n", "tc = 190.6\n"),
                                       encoding="utf-8")
    with reflux.Database(extra=[mine]) as d:
        assert d.extra == [str(mine)]
        assert d.compound("methane")["constants"]["tc"] == 190.6
        assert d.fluid("PR", ["methane"]).tc[0] == 190.6
        assert len(d.compounds()) == len(db.compounds())


def test_text_longer_than_its_first_buffer(db, monkeypatch):
    # The library says how much room text takes, and the package asks again.
    want = db.compound("methane")
    monkeypatch.setattr(_native, "TEXT_SIZE", 16)
    assert db.compound("methane") == want


def test_a_fluid_is_its_data(db):
    fl = db.fluid("pr", ["methane", "106-97-8", "Nitrogen"], alpha="mathias_copeman",
                  overrides={("nitrogen", "74-82-8"): 0.0311})
    assert fl.equation == "PR" and fl.nc == 3
    assert fl.names == ["Methane", "N-butane", "Nitrogen"]
    assert fl.cas == ["74-82-8", "106-97-8", "7727-37-9"]
    files = [component_file(c) for c in fl.cas]
    for name, key in (("tc", "tc"), ("pc", "pc"), ("omega", "omega"), ("molar_mass", "mw")):
        assert getattr(fl, name).tolist() == [f["constants"][key] for f in files]
    # kᵢⱼ: methane/n-butane's first listing in pr.ipd, line 133; the override;
    # n-butane/nitrogen from the bank.
    pairs = fl.setup["pairs"]
    assert [p["source"] for p in pairs] == ["file", "override", "file"]
    assert pairs[0]["line"] == 133 and fl.kij[0, 1] == fl.kij[1, 0] == 0.0244
    assert fl.kij[0, 2] == fl.kij[2, 0] == 0.0311
    line = VENDOR.joinpath("pr.ipd").read_text().splitlines()[pairs[2]["line"] - 1].split()
    assert {line[0], line[1]} == {"106-97-8", "7727-37-9"} and fl.kij[1, 2] == float(line[2])
    assert np.all(np.diag(fl.kij) == 0)
    assert fl.setup["components"][0]["alpha"].startswith("Mathias–Copeman, mc-PR.lib line ")
    # Read-only: the fluid in the library wouldn't change.
    with pytest.raises(ValueError):
        fl.kij[0, 1] = 0.5
    with pytest.raises(ValueError):
        fl.tc[0] = 200.0
    # No kᵢⱼ, and SRK's.
    assert not db.fluid("PR", ["methane", "n-butane"], kij=False).kij.any()
    assert db.fluid("SRK", ["methane", "n-butane"]).setup["kij_bank"] == "srk.ipd"


@pytest.mark.parametrize("args, kwargs, error, message", [
    (("PRX", ["methane"]), {}, ValueError, "no equation 'PRX': PR, SRK, RK or vdW"),
    (("PR", ["methane"]), {"alpha": "soave"}, ValueError, "no alpha 'soave'"),
    (("PR", "methane"), {}, TypeError, "names should be a list"),
    (("SRK", ["methane"]), {"alpha": "mathias_copeman"}, reflux.BadArgument,
     "fitted for PR, not SRK"),
    (("PR", []), {}, reflux.BadArgument, "a fluid needs at least one compound"),
    (("PR", ["methane", "74-82-8"]), {}, reflux.BadArgument, "which is in the list already"),
    (("PR", ["methane", "n-butane"]), {"overrides": {("methane", "methane"): 0.1}},
     reflux.BadArgument, "pairs `methane` with itself"),
    (("PR", ["methane", "n-butane"]), {"overrides": {("methane", "n-butane"): float("nan")}},
     reflux.BadArgument, r"pair_k\[0\] is nan"),
    (("PR", ["methane", "unobtanium"]), {}, reflux.NotFound, "no compound `unobtanium`"),
    (("PR", ["3967-55-3"]), {}, reflux.DataError,
     r"\(3967-55-3\) has no Tc, Pc, ω, so no equation of state can use it"),
])
def test_a_fluids_mistakes_are_refused(db, args, kwargs, error, message):
    with pytest.raises(error, match=message):
        db.fluid(*args, **kwargs)


def test_closed_is_closed():
    with reflux.Database() as d:
        fl = d.fluid("PR", ["methane", "n-butane"])
        handle = fl.handle
    # The fluid outlives its database, in the library too.
    assert "Database(" in repr(d) and repr(d).endswith(", closed)")
    with pytest.raises(ValueError, match="the database is closed"):
        d.compound("methane")
    setup = _native.text_out(_native.lib.reflux_fluid_setup, handle)
    assert tomllib.loads(setup) == fl.setup
    with fl:
        pass
    fl.close()
    with pytest.raises(ValueError, match="the fluid is closed"):
        _ = fl.handle
    with pytest.raises(reflux.BadHandle, match=f"handle {handle} was freed already"):
        _native.free(handle)
    with pytest.raises(TypeError, match="make a Fluid with Database.fluid"):
        reflux.Fluid()


def test_a_fluid_only_python_held_is_freed():
    import gc
    with reflux.Database() as d:
        fl = d.fluid("PR", ["methane"])
    handle = fl.handle
    del fl
    gc.collect()
    with pytest.raises(reflux.BadHandle, match="freed already"):
        _native.free(handle)


def moved_copy(root):
    """A data root away from the source tree, with methane and n-butane, whose
    pr.ipd lists 0.0313 for the pair where the source tree's lists 0.0244, and
    whose mc-PR.lib has no methane (D9)."""
    comps = root / "components"
    vendor = root / "vendor" / "chemsep"
    comps.mkdir(parents=True)
    vendor.mkdir(parents=True)
    for cas in ("74-82-8", "106-97-8"):
        shutil.copy(COMPONENTS / f"{cas}.toml", comps)
    ipd = (VENDOR / "pr.ipd").read_bytes()
    assert ipd.count(b"106-97-8   .244e-1") == 1
    (vendor / "pr.ipd").write_bytes(ipd.replace(b"106-97-8   .244e-1", b"106-97-8   .313e-1"))
    lib = (VENDOR / "mc-PR.lib").read_bytes().splitlines(keepends=True)
    (vendor / "mc-PR.lib").write_bytes(b"".join(x for x in lib if not x.startswith(b"74-82-8\t")))
    return root


def test_a_copy_of_the_data_is_used_whole(tmp_path):
    root = moved_copy(tmp_path / "data")
    with reflux.Database(root) as d:
        assert len(d.compounds()) == 2
        fl = d.fluid("PR", ["methane", "n-butane"], alpha="mathias_copeman")
    assert fl.kij[0, 1] == 0.0313
    assert fl.setup["components"][0]["alpha"] == "the equation's own: mc-PR.lib has no Methane"


def test_reflux_data_moves_the_default(tmp_path):
    # Set before the process starts, and set by Python after the import: the
    # library reads the same environment.
    root = moved_copy(tmp_path / "data")
    code = """
import os, reflux
print(reflux.data_root())
print(reflux.Database().fluid("PR", ["methane", "n-butane"]).kij[0, 1])
os.environ["REFLUX_DATA"] = {other!r}
print(reflux.data_root())
"""
    out = isolated(code.format(other=str(DATA)), env={"REFLUX_DATA": str(root)}).split("\n")
    assert out[:3] == [str(root), "0.0313", str(DATA)]
