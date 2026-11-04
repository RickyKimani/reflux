"""Properties and vapour pressure at many points (M10a step 3). Broadcasting
gives each point what a call for it alone gives, bit for bit; the energies
keep their definitions; the stable root has the lower Gibbs energy; at a
vapour pressure the liquid's fugacity is the vapour's; and a point without an
answer says why, while the others are unaffected (D6, D7)."""

import numpy as np
import pytest

import reflux
from reflux import Phase, Status, points
from support import ROOT


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


@pytest.fixture(scope="module")
def fl(db):
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        yield f


NAMES = ("z", "v", "h", "s", "g", "u", "a", "cp", "cv", "w", "mu_jt")


def same(a, b):
    """Equal bit for bit, NaN included."""
    return np.array_equal(a, b, equal_nan=True)


def test_a_grid_is_its_points_one_by_one(fl):
    T = np.array([250.0, 300.0, 400.0])[:, None]
    P = np.array([1e5, 20e5, 60e5])
    z = np.array([0.4, 0.6])
    grid = fl.props(T, P, z)
    assert grid.status.shape == (3, 3) and grid.ln_phi.shape == (3, 3, 2)
    for i in range(3):
        for j in range(3):
            one = fl.props(T[i, 0], P[j], z)
            assert one.status.shape == ()
            for name in NAMES:
                assert getattr(grid, name)[i, j] == getattr(one, name), name
            assert same(grid.ln_phi[i, j], one.ln_phi)
            assert grid.phase[i, j] == one.phase and grid.status[i, j] == Status.OK
    # A composition per point, the same everywhere, is the one composition.
    rows = fl.props(T, P, np.broadcast_to(z, (3, 3, 2)).copy())
    assert all(same(getattr(rows, name), getattr(grid, name)) for name in NAMES)
    # Amounts, not fractions: twice the moles is the same per mole, to rounding.
    twice = fl.props(T, P, 2 * z)
    for name in NAMES:
        assert np.allclose(getattr(twice, name), getattr(grid, name), rtol=1e-13, atol=0), name


def test_inputs_already_in_shape_arent_copied():
    a = np.linspace(200.0, 400.0, 6).reshape(2, 3)
    t, n = points.flat(a, (2, 3), 6)
    assert n == 6 and np.shares_memory(t, a)
    # A scalar, or one value, is one value for every point.
    assert points.flat(np.float64(300.0), (2, 3), 6)[1] == 1
    # What has to be broadcast or converted is copied.
    assert not np.shares_memory(points.flat(a[:, :1], (2, 3), 6)[0], a)
    assert not np.shares_memory(points.flat(a.T, (3, 2), 6)[0], a)
    z = np.array([[0.4, 0.6], [0.5, 0.5]])
    rows, k = points.flat_rows(z, (2,), 2, 2)
    assert k == 2 and np.shares_memory(rows, z)
    assert points.flat_rows(z[:1], (2,), 2, 2)[1] == 1


def test_the_energies_keep_their_definitions(fl):
    T = np.array([200.0, 300.0, 450.0])
    P = np.array([2e5, 30e5, 80e5])
    r = fl.props(T[:, None], P, [0.3, 0.7])
    T, P = T[:, None], P[None, :]
    # g = h − Ts, a = u − Ts and u = h − Pv, to the rounding of the larger term.
    for got, want, scale in ((r.g, r.h - T * r.s, np.abs(T * r.s)), (r.a, r.u - T * r.s,
                             np.abs(T * r.s)), (r.u, r.h - P * r.v, np.abs(r.h))):
        assert np.all(np.abs(got - want) <= 1e-12 * (scale + np.abs(want)))
    assert np.all(r.cp > r.cv) and np.all(r.w > 0)
    # Z = Pv/RT.
    assert np.allclose(r.z, P * r.v / (8.314462618 * T), rtol=1e-12, atol=0)
    # At 1 Pa, an ideal gas: Z = 1 and ln φ = 0, to a few parts in 10⁶.
    low = fl.props(300.0, 1.0, [0.3, 0.7])
    assert abs(low.z - 1) < 1e-6 and np.all(np.abs(low.ln_phi) < 1e-6)
    assert low.phase == Phase.VAPOR


def test_the_stable_root_has_the_lower_gibbs_energy(fl):
    # Pressures where the cubic has both roots at 300 K, and the stable one
    # changes between them.
    T, P, z = 300.0, np.linspace(2e5, 10e5, 8), [0.1, 0.9]
    liq = fl.props(T, P, z, root="liquid")
    vap = fl.props(T, P, z, root="vapor")
    stable = fl.props(T, P, z, root="stable")
    assert np.all(liq.v < vap.v) and np.all(liq.phase == Phase.LIQUID)
    assert np.all(vap.phase == Phase.VAPOR)
    lower = np.where(liq.g < vap.g, liq.v, vap.v)
    assert same(stable.v, lower)
    assert np.any(stable.v == liq.v) and np.any(stable.v == vap.v)
    # At 20 bar there's only the liquid root, and every root asked for is it.
    one = [fl.props(T, 20e5, z, root=r) for r in ("liquid", "vapor", "stable")]
    assert one[0].v == one[1].v == one[2].v and one[1].phase == Phase.LIQUID


def test_at_its_vapour_pressure_a_liquid_and_its_vapour_agree(db):
    with db.fluid("PR", ["n-butane"], alpha="mathias_copeman") as pure:
        T = np.linspace(200.0, 420.0, 12)
        ps = pure.psat(T)
        assert ps.p.shape == (12, 1) and np.all(ps.ok)
        liq = pure.props(T, ps.p[:, 0], root="liquid")
        vap = pure.props(T, ps.p[:, 0], root="vapor")
        assert np.all(liq.v < vap.v)
        # Equal fugacities, to the iteration's tolerance (1e-10 in ln f).
        assert np.all(np.abs(liq.ln_phi - vap.ln_phi) < 1e-9)
        # Psat rises with T.
        assert np.all(np.diff(ps.p[:, 0]) > 0)


def test_vapour_pressure_by_compound(fl):
    T = np.array([[150.0, 300.0], [400.0, 430.0]])
    both = fl.psat(T)
    assert both.p.shape == (2, 2, 2)
    # Methane's Tc is 190.56 K and n-butane's 425.12 K.
    assert both.status.tolist() == [[[0, 0], [7, 0]], [[7, 0], [7, 7]]]
    assert np.isnan(both.p[0, 1, 0]) and both.p[0, 1, 1] > 0
    for which in (1, -1, "n-butane", "N-BUTANE", "106-97-8"):
        one = fl.psat(T, which)
        assert one.p.shape == (2, 2) and same(one.p, both.p[..., 1])
    assert both.reason((0, 1, 0)) == ("Methane: T = 300.0 K is at or above its Tc, 190.56 K, "
                                      "so it has no vapour pressure")
    assert both.reason((0, 0, 0)) == ""
    # The last axis says whose reason: n-butane's at 430 K, not methane's.
    assert both.reason((1, 1, 1)).startswith("N-butane: T = 430.0 K is at or above its Tc")
    assert both.reason((1, 1, -2)).startswith("Methane: T = 430.0 K")
    assert fl.psat(np.nan, 0).status == Status.BAD_INPUT
    with pytest.raises(KeyError, match="'water' isn't in the fluid: Methane, N-butane"):
        fl.psat(300.0, "water")
    with pytest.raises(IndexError):
        fl.psat(300.0, 2)


def test_a_pure_fluids_vapour_pressure_keeps_its_compound_axis(db):
    with db.fluid("PR", ["n-butane"]) as pure:
        s = pure.psat([300.0, 430.0])
        assert s.p.shape == (2, 1) and s.status[:, 0].tolist() == [Status.OK, Status.SUPERCRITICAL]
        assert s.reason((1, 0)).startswith("N-butane: T = 430.0 K is at or above its Tc")
        assert s.reason((0, -1)) == ""
        one = pure.psat([300.0, 430.0], 0)
        assert one.p.shape == (2,) and one.reason(-1) == s.reason((1, 0))


def test_a_point_without_an_answer_leaves_the_others(fl):
    T = np.array([300.0, np.nan, 350.0, 300.0, 300.0])
    P = np.array([10e5, 10e5, 0.0, 10e5, 10e5])
    z = np.array([[0.4, 0.6], [0.4, 0.6], [0.4, 0.6], [-0.1, 1.1], [0.0, 0.0]])
    r = fl.props(T, P, z)
    assert r.status.tolist() == [0, 1, 1, 1, 1] and r.ok.tolist() == [True] + [False] * 4
    assert np.all(np.isnan(r.h[1:])) and np.all(np.isnan(r.ln_phi[1:]))
    assert np.all(r.phase[1:] == Phase.NONE)
    alone = fl.props(300.0, 10e5, [0.4, 0.6])
    assert r.h[0] == alone.h
    assert r.reason(0) == ""
    assert r.reason(1) == "T = nan K: a temperature must be finite and above 0"
    assert r.reason(2) == "P = 0.0 Pa: a pressure must be finite and above 0"
    assert r.reason(-2) == "z[0] = -0.1: amounts must be finite and not negative"
    assert r.reason(4) == "z adds up to 0.0: it must be finite and above 0"
    with pytest.raises(IndexError):
        r.reason(5)
    with pytest.raises(IndexError):
        r.reason((0, 0))


def test_beyond_the_model_is_no_answer(fl):
    # Above about 1e10 K the ideal-gas Cp correlations overflow, so H and S
    # are NaN; at 1e-300 Pa the vapour's volume is past what a float holds.
    r = fl.props([300.0, 1e10, 300.0], [1e5, 1e5, 1e-300], [0.5, 0.5])
    assert r.status.tolist() == [Status.OK, Status.NOT_FINITE, Status.NO_VOLUME_ROOT]
    assert r.reason(1).startswith("h is nan at T = 1e+10 K")
    assert r.reason(2).startswith("no volume root at T = 300.0 K")
    assert np.all(np.isnan(r.z[1:])) and np.all(r.phase[1:] == Phase.NONE)


def test_mistakes_in_the_call(fl, db):
    with pytest.raises(ValueError, match="no root 'gas'"):
        fl.props(300.0, 1e5, [0.5, 0.5], root="gas")
    with pytest.raises(TypeError, match="z is needed: the fluid has 2 compounds"):
        fl.props(300.0, 1e5)
    with pytest.raises(ValueError, match="z's last axis should have the fluid's 2 compounds"):
        fl.props(300.0, 1e5, [0.2, 0.3, 0.5])
    with pytest.raises(ValueError):
        fl.props([300.0, 310.0, 320.0], [1e5, 2e5], [0.5, 0.5])
    empty = fl.props(np.empty(0), 1e5, [0.5, 0.5])
    assert empty.h.shape == (0,) and empty.ln_phi.shape == (0, 2)
    with db.fluid("PR", ["methane"]) as closed:
        pass
    with pytest.raises(ValueError, match="the fluid is closed"):
        closed.props(300.0, 1e5)


def test_absolute_says_whether_energies_have_a_basis(fl, tmp_path):
    assert fl.absolute
    assert all("dhf" in c and "s0" in c for c in fl.setup["components"])
    # Methane of your own without ΔH°f: its energies sit on an arbitrary zero,
    # though differences at one composition still hold.
    text = (ROOT / "data" / "components" / "74-82-8.toml").read_text(encoding="utf-8")
    assert text.count("dhf = -74520.0\n") == 1
    (tmp_path / "74-82-8.toml").write_text(text.replace("dhf = -74520.0\n", ""), encoding="utf-8")
    with reflux.Database(extra=[tmp_path]) as d, d.fluid("PR", ["methane", "n-butane"]) as mine:
        assert not mine.absolute and "dhf" not in mine.setup["components"][0]
        assert "s0" in mine.setup["components"][0]
        r = mine.props(300.0, [1e5, 2e5], [0.5, 0.5])
        assert np.all(r.ok)
        ref = fl.props(300.0, [1e5, 2e5], [0.5, 0.5])
        assert abs((r.h[1] - r.h[0]) - (ref.h[1] - ref.h[0])) < 1e-9 * abs(ref.h[0])
