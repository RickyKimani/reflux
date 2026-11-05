"""Stability and the PT flash at many points (M10a step 4). Against `thermo`'s
splits (testdata/m3_thermo.json) through the package; a grid is its points
one by one; a split's balances close and its fugacities agree, checked from
the arrays; the whole is its phases; each phase is in its label's place; and
stability agrees with the flash (D6, D7)."""

import numpy as np
import pytest

import reflux
from reflux import Phase, Status
from support import M3, m3_compounds, m3_fluid
FLASH = ("beta", "t", "p", "h", "s", "v")


def same(a, b):
    """Equal bit for bit, NaN included."""
    return np.array_equal(a, b, equal_nan=True)


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


@pytest.fixture(scope="module")
def fl(db):
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        yield f


@pytest.fixture(scope="module")
def m3(tmp_path_factory):
    """The database, with M3's constants for its eleven compounds."""
    with reflux.Database(extra=[m3_compounds(tmp_path_factory.mktemp("m3"))]) as d:
        yield d


def test_the_splits_are_thermos(m3):
    # β and the compositions to 1e-8 and the phases' volumes to 1e-9, as
    # flash/thermo_test.v holds the V API to: all 83 of M3's PT states.
    worst_split, worst_v, splits = 0.0, 0.0, 0
    for c in M3["pt"]:
        with m3_fluid(m3, c) as f:
            r = f.flash_pt(c["t"], c["p"], c["z"])
        at = f"{c['eos']} {c['components']} z {c['z']} at {c['t']} K, {c['p']} Pa"
        assert r.status == Status.OK, f"{at}: {r.reason(())}"
        assert r.phases == len(c["phases"]), at
        if r.phases == 2:
            splits += 1
        # The fixture's phases are in order of mass density. A split's less
        # dense phase is the vapour, whatever thermo labels it, and a single
        # phase is in its label's place, as thermo labels it too.
        if r.phases == 2:
            pairs = [(r.vapor, c["phases"][0]), (r.liquid, c["phases"][1])]
        else:
            want = c["phases"][0]
            ph, other = (r.vapor, r.liquid) if want["label"] == "vapor" else (r.liquid, r.vapor)
            assert ph.present and not other.present, at
            pairs = [(ph, want)]
        for ph, want in pairs:
            worst_split = max(worst_split, abs(ph.fraction - want["fraction"]),
                              *np.abs(ph.x - want["z"]))
            worst_v = max(worst_v, abs(ph.v - want["v"]) / want["v"])
    assert splits == 39
    assert worst_split <= 1e-8, worst_split
    assert worst_v <= 1e-9, worst_v


def test_a_grid_is_its_points_one_by_one(fl):
    T = np.array([250.0, 300.0, 350.0])[:, None]
    P = np.array([2e5, 20e5, 60e5])
    r = fl.flash_pt(T, P, [0.5, 0.5])
    assert r.beta.shape == (3, 3) and r.x.shape == (3, 3, 2, 2)
    for i in range(3):
        for j in range(3):
            one = fl.flash_pt(T[i, 0], P[j], [0.5, 0.5])
            assert one.beta.shape == () and one.x.shape == (2, 2)
            for name in FLASH:
                assert same(getattr(r, name)[i, j], getattr(one, name)), name
            assert same(r.x[i, j], one.x) and same(r.fraction[i, j], one.fraction)
            assert r.phases[i, j] == one.phases


def test_a_splits_balances_close(db):
    # From the arrays alone: the fractions add to 1, β is the vapour's, each
    # compound's moles balance, and its fugacity is the same in both phases,
    # computed again through props.
    with db.fluid("PR", ["methane", "ethane", "n-butane"]) as f:
        z = np.array([0.3, 0.3, 0.4])
        T = np.linspace(240.0, 320.0, 9)
        r = f.flash_pt(T, 20e5, z)
        two = r.phases == 2
        assert two.sum() >= 5
        v, l_ = r.vapor, r.liquid
        assert np.all(np.abs(v.fraction[two] + l_.fraction[two] - 1) < 1e-15)
        assert same(r.beta[two], v.fraction[two])
        balance = v.fraction[two, None] * v.x[two] + l_.fraction[two, None] * l_.x[two]
        assert np.all(np.abs(balance - z) < 1e-14)
        pv = f.props(T[two], 20e5, v.x[two])
        pl = f.props(T[two], 20e5, l_.x[two])
        assert np.all(np.abs(pv.ln_phi + np.log(v.x[two]) - pl.ln_phi - np.log(l_.x[two])) < 1e-11)
        assert np.all(v.v[two] > l_.v[two])


def test_the_whole_is_its_phases(fl):
    r = fl.flash_pt(300.0, np.linspace(2e5, 60e5, 12), [0.5, 0.5])
    v, l_ = r.vapor, r.liquid
    for name in ("h", "s", "v"):
        whole = getattr(r, name)
        two = r.phases == 2
        # Summed vapour first, as the flash gives its phases.
        assert same(whole[two], v.fraction[two] * getattr(v, name)[two]
                    + l_.fraction[two] * getattr(l_, name)[two])
        one = r.phases == 1
        alone = np.where(v.present[one], getattr(v, name)[one], getattr(l_, name)[one])
        assert same(whole[one], alone)
    # A single phase's properties are props' at the stable root.
    single = np.flatnonzero(r.phases == 1)
    assert single.size > 0
    p = np.linspace(2e5, 60e5, 12)[single]
    pr = fl.props(300.0, p, [0.5, 0.5])
    assert same(r.h[single], pr.h) and same(r.v[single], pr.v)


def test_a_single_phase_takes_its_labels_place(fl):
    r = fl.flash_pt([400.0, 250.0], [5e5, 50e5], [[0.5, 0.5], [0.2, 0.8]])
    assert r.phases.tolist() == [1, 1] and r.beta.tolist() == [1.0, 0.0]
    assert r.vapor.present.tolist() == [True, False]
    assert r.liquid.present.tolist() == [False, True]
    assert np.all(np.isnan(r.x[0, Phase.LIQUID])) and np.all(np.isnan(r.x[1, Phase.VAPOR]))
    assert r.vapor.x[0].tolist() == [0.5, 0.5] and r.liquid.x[1].tolist() == [0.2, 0.8]
    # The feed is normalised, as amounts.
    twice = fl.flash_pt(400.0, 5e5, [1.0, 1.0])
    assert twice.vapor.x.tolist() == [0.5, 0.5]


def test_stability_agrees_with_the_flash(fl):
    T = np.linspace(220.0, 380.0, 9)[:, None]
    P = np.linspace(2e5, 80e5, 9)
    st = fl.stability(T, P, [0.5, 0.5])
    r = fl.flash_pt(T, P, [0.5, 0.5])
    assert np.all(st.ok) and np.all(r.ok)
    assert np.array_equal(~st.stable, r.phases == 2)
    assert np.all(st.tm[~st.stable] < -1e-10) and np.all(st.tm[st.stable] > -1e-10)
    assert st.w.shape == (9, 9, 2) and np.allclose(st.w.sum(axis=-1), 1, atol=1e-12)
    assert st.stable.dtype == bool


def test_a_flash_that_fails_says_why(fl, db):
    r = fl.flash_pt([300.0, np.nan, 1e10], [20e5, 20e5, 1e5], [0.5, 0.5])
    assert r.status.tolist() == [Status.OK, Status.BAD_INPUT, Status.NOT_FINITE]
    assert r.phases.tolist() == [2, 0, 0] and np.all(np.isnan(r.x[1:]))
    for ph in (r.liquid, r.vapor):
        for values in (ph.fraction, ph.v, ph.h, ph.s):
            assert np.all(np.isnan(values[1:]))
    assert not r.liquid.present[1:].any() and not r.vapor.present[1:].any()
    assert r.reason(1) == "T = nan K: a temperature must be finite and above 0"
    assert r.reason(2).startswith("h is nan at T = 1e+10 K")
    # Water, n-decane and methane with a large water–hydrocarbon kᵢⱼ have a
    # third phase at 350 K and 5 bar.
    with db.fluid("PR", ["water", "n-decane", "methane"], kij=False,
                  overrides={("water", "n-decane"): 0.5, ("water", "methane"): 0.5,
                             ("n-decane", "methane"): 0.04}) as three:
        r3 = three.flash_pt(350.0, 5e5, [0.5, 0.3, 0.2])
        assert r3.status == Status.MORE_PHASES
        assert "unstable" in r3.reason(()) and "three or more phases" in r3.reason(())
    st = fl.stability([300.0, -1.0], 1e5, [0.5, 0.5])
    assert st.status.tolist() == [Status.OK, Status.BAD_INPUT]
    assert not st.stable[1] and np.isnan(st.tm[1])
    assert st.reason(1) == "T = -1.0 K: a temperature must be finite and above 0"
    with pytest.raises(ValueError, match="z's last axis"):
        fl.flash_pt(300.0, 1e5, [1.0])
