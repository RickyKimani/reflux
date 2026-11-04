"""The C API against an independent reference (M1 spike's scenarios, now
through the package): testdata/m1_thermo.json, PR and SRK states computed by
Caleb Bell's `thermo` with Smith, Van Ness & Abbott's constants. Those
constants go in a directory of the caller's own compounds, so the fluid is
built through the database as any other; then Z, V, ln φ and the residual H
and S must agree to M1's 1e-9. And 10,000 points come back in one call."""

import json
import re
import shutil

import numpy as np
import pytest

import reflux
from support import ROOT

FIXTURE = json.loads((ROOT / "testdata" / "m1_thermo.json").read_text(encoding="utf-8"))
CAS = {"methane": "74-82-8", "n-butane": "106-97-8", "carbon dioxide": "124-38-9"}
R = 8.314462618


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    """The database, with the fixture's constants for its three compounds."""
    mine = tmp_path_factory.mktemp("m1")
    for c in FIXTURE["components"]:
        cas = CAS[c["name"]]
        text = (ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8")
        for key in ("tc", "pc", "omega", "mw"):
            text, n = re.subn(rf"^{key} = .*$", f"{key} = {c[key]!r}", text, flags=re.M)
            assert n == 1, (cas, key)
        (mine / f"{cas}.toml").write_text(text, encoding="utf-8")
    with reflux.Database(extra=[mine]) as d:
        yield d
    shutil.rmtree(mine, ignore_errors=True)


def agree(got, want, tol, floor):
    return abs(got - want) <= tol * max(floor, abs(want))


@pytest.mark.parametrize("case", FIXTURE["cases"],
                         ids=lambda c: f"{c['eos']}-{len(c['components'])}-{c['t']}-{c['p']}")
def test_the_same_models_as_thermo(db, case):
    names = [FIXTURE["components"][i]["name"] for i in case["components"]]
    nc = len(names)
    k = case["kij"]
    overrides = {(names[i], names[j]): k[i * nc + j] for i in range(nc) for j in range(i + 1, nc)}
    with db.fluid(case["eos"], [CAS[n] for n in names], kij=False, overrides=overrides) as fl:
        assert fl.tc.tolist() == [FIXTURE["components"][i]["tc"] for i in case["components"]]
        t, p, n = case["t"], case["p"], case["n"]
        for ph in case["phases"]:
            at = f"{case['eos']} {names} {ph['phase']} at {t} K, {p} Pa"
            r = fl.props(t, p, n, root=ph["phase"])
            assert r.status == reflux.Status.OK, at
            assert agree(r.z, ph["z"], 1e-9, 0.0), f"{at}: Z {r.z} vs {ph['z']}"
            assert agree(r.v, ph["v"], 1e-9, 0.0), f"{at}: V"
            for i, lp in enumerate(ph["ln_phi"]):
                assert agree(r.ln_phi[i], lp, 1e-9, 1.0), f"{at}: ln φ[{i}]"
            # Residual H and S, as the difference from the same T and z at
            # 1 mPa, where the residual part is below 1e-11 of these: the
            # ideal gas's H doesn't depend on P, and its S goes as −R ln P.
            ig = fl.props(t, 1e-3, n, root="vapor")
            h_dep = r.h - ig.h
            s_dep = r.s - ig.s + R * np.log(p / 1e-3)
            assert agree(h_dep, ph["h_dep"], 1e-9, R * t), f"{at}: Hʳ {h_dep} vs {ph['h_dep']}"
            assert agree(s_dep, ph["s_dep"], 1e-9, R), f"{at}: Sʳ {s_dep} vs {ph['s_dep']}"


def test_ten_thousand_points_in_one_call(db):
    # The spike's vectorised scenario, through the package: n-butane vapour at
    # 5 bar from 300 to 500 K.
    with db.fluid("PR", ["n-butane"]) as fl:
        T = np.linspace(300.0, 500.0, 10_000)
        r = fl.props(T, 5e5, root="vapor")
        assert np.all(r.ok) and np.all((r.z > 0.5) & (r.z < 1.0))
        assert np.all(np.diff(r.z) > 0)
