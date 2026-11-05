"""The other flashes (M10a step 5): PH, PS, TV, UV and the vapour-fraction
flashes, with bubble and dew points, as specs of the one flash call. Each
lands back on the PT state it was given; the boundaries and saturation states
are thermo's (testdata/m3_thermo.json); a vapour-fraction flash gives its
incipient phase; and each kind of failure has its status."""

import numpy as np
import pytest

import reflux
from reflux import Status
from support import M3, m3_compounds, m3_fluid


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


@pytest.fixture(scope="module")
def m3(tmp_path_factory):
    with reflux.Database(extra=[m3_compounds(tmp_path_factory.mktemp("m3"))]) as d:
        yield d


def close(got, want, rel):
    return np.all(np.abs(got - want) <= rel * np.abs(want))


def test_every_flash_lands_on_the_pt_state(db):
    # A grid of PT states, splits and single phases, then each flash given
    # what the PT flash found: T, P and β come back to 1e-8.
    with db.fluid("PR", ["methane", "ethane", "n-butane"]) as f:
        z = [0.3, 0.3, 0.4]
        T = np.linspace(230.0, 400.0, 6)[:, None]
        P = np.array([5e5, 20e5, 50e5])
        pt = f.flash_pt(T, P, z)
        assert np.all(pt.ok) and np.any(pt.phases == 2) and np.any(pt.phases == 1)
        u = pt.h - pt.p * pt.v
        for r in (f.flash_ph(P, pt.h, z), f.flash_ps(P, pt.s, z), f.flash_tv(T, pt.v, z),
                  f.flash_uv(u, pt.v, z)):
            assert np.all(r.ok), r.reason(np.unravel_index(np.argmin(r.ok), r.ok.shape))
            assert close(r.t, pt.t, 1e-8) and close(r.p, pt.p, 1e-8)
            assert np.all(np.abs(r.beta - pt.beta) <= 1e-8)
            assert np.array_equal(r.phases, pt.phases)
        two = pt.phases == 2
        Tb = np.broadcast_to(T, pt.t.shape)[two]
        Pb = np.broadcast_to(P, pt.p.shape)[two]
        tvf = f.flash_tvf(Tb, pt.beta[two], z)
        pvf = f.flash_pvf(Pb, pt.beta[two], z)
        assert np.all(tvf.ok) and np.all(pvf.ok)
        # Near the critical region two pressures (or temperatures) can have the
        # same β, and Newton finds the one its start leads to: here one of ten
        # splits gives 80.6 bar back for 50. So each answer must be a state
        # with that β, which a PT flash there confirms, and the rest land on
        # the state they came from.
        for r in (tvf, pvf):
            back = f.flash_pt(r.t, r.p, z)
            assert np.all(np.abs(back.beta - pt.beta[two]) <= 1e-8)
        assert np.sum(np.abs(tvf.p - Pb) <= 1e-8 * Pb) >= 9
        assert np.sum(np.abs(pvf.t - Tb) <= 1e-8 * Tb) >= 9


def test_the_boundaries_are_thermos(m3):
    # thermo's 40 bubble and dew points, through the package: T or P to 1e-11
    # relative and the incipient phase to 1e-11, as flash/vf_test.v holds the
    # V API to.
    worst_tp, worst_w, runs = 0.0, 0.0, 0
    for c in M3["boundaries"]:
        if c.get("error") or c.get("rejected"):
            continue
        with m3_fluid(m3, c) as f:
            if c["spec"] == "p":
                r = f.flash_pvf(c["p"], c["vf"], c["z"])
                worst_tp = max(worst_tp, abs(r.t - c["t"]) / c["t"])
            else:
                r = f.flash_tvf(c["t"], c["vf"], c["z"])
                worst_tp = max(worst_tp, abs(r.p - c["p"]) / c["p"])
        assert r.status == Status.OK and r.phases == 2
        incipient = r.vapor if c["vf"] == 0.0 else r.liquid
        assert incipient.fraction == 0.0
        worst_w = max(worst_w, *np.abs(incipient.x - c["incipient"]))
        runs += 1
    assert runs == 40
    assert worst_tp <= 1e-11 and worst_w <= 1e-11, (worst_tp, worst_w)


def test_saturation_is_thermos(m3):
    # A pure compound's answer is its saturation point whatever β: thermo's 36
    # saturation states to 1e-11, with two distinct phases.
    worst = 0.0
    for c in M3["saturation"]:
        with m3_fluid(m3, c) as f:
            if c["spec"] == "p":
                r = f.flash_pvf(c["p"], c["vf"])
                worst = max(worst, abs(r.t - c["t"]) / c["t"])
            else:
                r = f.flash_tvf(c["t"], c["vf"])
                worst = max(worst, abs(r.p - c["p"]) / c["p"])
        assert r.status == Status.OK and r.beta == c["vf"] and r.phases == 2
        assert r.vapor.v > r.liquid.v
    assert worst <= 1e-11, worst


def test_bubble_and_dew_points_give_the_incipient_phase(db):
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        z = np.array([0.2, 0.8])
        T = np.linspace(250.0, 350.0, 5)
        bub, dew = f.bubble_p(T, z), f.dew_p(T, z)
        assert np.all(bub.ok) and np.all(dew.ok)
        assert np.all(bub.phases == 2) and np.all(bub.beta == 0.0)
        assert np.all(bub.vapor.fraction == 0.0) and np.all(bub.liquid.fraction == 1.0)
        assert np.all(bub.liquid.x == z) and np.allclose(bub.vapor.x.sum(axis=-1), 1, atol=1e-15)
        assert np.all(dew.beta == 1.0) and np.all(dew.vapor.x == z)
        assert np.all(bub.p > dew.p)
        # The named ones are the vapour-fraction flashes at 0 and 1.
        assert np.array_equal(bub.p, f.flash_tvf(T, 0.0, z).p)
        bt, dt = f.bubble_t(10e5, z), f.dew_t(10e5, z)
        assert bt.ok and dt.ok and bt.t < dt.t
        assert bt.t == f.flash_pvf(10e5, 0.0, z).t


def test_each_kind_of_failure_has_its_status(m3, db):
    # M3's methane/n-butane, 40/60, whose envelope's top is near 93 bar: there's
    # no dew point at 90 bar and the iteration can't converge (85 to 93 bar
    # alike), and neither point exists at 100 bar (94 to 106 bar alike).
    case = {"components": [0, 3], "kij": [0.0, 0.02, 0.02, 0.0]}
    with m3_fluid(m3, case) as f:
        z = [0.4, 0.6]
        r = f.dew_t([90e5, 100e5], z)
        assert r.status.tolist() == [Status.NO_CONVERGENCE, Status.OUTSIDE_ENVELOPE]
        assert r.reason(0).startswith("flash: the dew point at 9e+06 Pa")
        assert f.bubble_t(100e5, z).status == Status.OUTSIDE_ENVELOPE
        assert np.all(r.phases == 0) and np.all(np.isnan(r.x))
        # Further above it, a vapour-fraction flash fails as its iteration
        # goes: mostly outside the envelope, and now and then on PR's
        # liquid–liquid branch near 43 K, with its phases swapped or a third
        # phase. Those are isolated points, which rounding moves from one
        # platform to another (M10a's exit review: Ubuntu's 110 bar wasn't
        # Windows'), so here each failure need only be one of those kinds,
        # with its own reason. capi/flash_test.v checks each kind's number.
        why = {Status.OUTSIDE_ENVELOPE: "There is very likely no such point",
               Status.PHASES_SWAPPED: "where the phase of fraction β is the denser",
               Status.MORE_PHASES: "is unstable (tm =",
               Status.NO_CONVERGENCE: "flash: "}
        seen = set()
        for z1 in (0.4, 0.8, 0.95):
            g = f.flash_pvf(np.arange(105.0, 200.0, 5.0)[:, None] * 1e5,
                            [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0], [z1, 1 - z1])
            assert set(g.status.ravel().tolist()) <= {Status.OK, *why}, z1
            for s in set(g.status.ravel().tolist()) - {Status.OK}:
                at = tuple(int(i) for i in np.argwhere(g.status == s)[0])
                assert why[s] in g.reason(at), (z1, Status(s), g.reason(at))
                seen.add(s)
        # Outside the envelope above all; the rest as rounding has it.
        assert Status.OUTSIDE_ENVELOPE in seen
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        # An enthalpy no temperature reaches is the flash's own bad input.
        far = f.flash_ph(1e5, 1e10, [0.5, 0.5])
        assert far.status == Status.BAD_INPUT and far.reason(()).startswith("flash: ")
        bad = f.flash_tvf([300.0, 300.0], [1.5, -0.1], [0.5, 0.5])
        assert bad.status.tolist() == [Status.BAD_INPUT] * 2
        assert bad.reason(0) == "β = 1.5: a vapour fraction must be between 0 and 1"
        assert f.flash_tv(300.0, -1.0, [0.5, 0.5]).reason(()) == (
            "V = -1.0 m3/mol: a molar volume must be finite and above 0")
        assert f.flash_uv(np.nan, 1e-3, [0.5, 0.5]).reason(()) == "U = nan J/mol: it must be finite"
