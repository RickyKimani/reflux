"""A phase's transport from Python (M5b's D17): Fluid.transport's shapes,
statuses and reasons. Its values are the V API's bit for bit in
test_reference.py; here, IAPWS-95 water's against R12-08's, R15-11's and
R1-76's own tables."""

import numpy as np
import pytest

import reflux
from reflux import Phase, Status


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


def test_values_broadcast_with_a_status_each(db):
    with db.fluid("PR", ["n-butane", "n-pentane"]) as f:
        tr = f.transport([[280.0], [300.0]], [5e5, 8e5, 10e5], [0.4, 0.6])
        assert tr.viscosity.shape == tr.thermal_conductivity.shape == (2, 3)
        assert tr.surface_tension.shape == (2, 3)
        assert tr.status.shape == (2, 3, 3) and tr.ok.all()
        # A liquid's viscosity falls with T and rises with P.
        assert (tr.viscosity[0] > tr.viscosity[1]).all()
        assert (np.diff(tr.viscosity, axis=1) > 0).all()
        assert tr.reason((1, 2)) == ""
        assert repr(tr) == "Transport((2, 3), liquid, 6 of 6 with every value)"


def test_a_vapour_has_no_surface_tension(db):
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        for phase in ("vapor", Phase.VAPOR):
            tr = f.transport(400.0, 1e5, [0.5, 0.5], phase=phase)
            assert tr.status.tolist() == [Status.OK, Status.OK, Status.NO_VALUE]
            assert np.isnan(tr.surface_tension) and tr.reason(()) == ""
            # A dilute gas: about 11 µPa·s and 0.035 W/(m·K).
            assert 1e-5 < tr.viscosity < 1.3e-5 and 0.03 < tr.thermal_conductivity < 0.04


def test_a_point_without_values_says_why(db):
    with db.fluid("PR", ["methane", "n-butane"]) as f:
        tr = f.transport([-1.0, 250.0, 250.0], 20e5, [[0.3, 0.7], [0.3, 0.7], [0.0, 1.0]])
        assert tr.status[0].tolist() == [Status.BAD_INPUT] * 3
        assert tr.status[1].tolist() == [Status.NO_VALUE] * 3
        assert tr.ok[2].all()
        assert np.isnan(tr.viscosity[:2]).all()
        assert tr.reason(0) == "T = -1.0 K: a temperature must be finite and above 0"
        assert tr.reason(1).endswith("Methane, above its Tc, is 0.3 of the liquid's moles, more "
                                     "than 0.1 (D9)")
        with pytest.raises(IndexError):
            tr.reason(3)


def test_a_phase_is_liquid_or_vapor(db):
    with db.fluid("PR", ["methane"]) as f:
        for bad in ("gas", "stable", Phase.NONE, 1):
            with pytest.raises(ValueError, match="no phase"):
                f.transport(300.0, 1e5, phase=bad)


def test_iapws95_water_is_iapws(db):
    # M5a's functions, which water/transport_test.v holds to IAPWS's
    # tables, give at 300 K and 0.1 MPa μ 853.74 µPa·s, λ 0.6095 W/(m·K)
    # and σ 71.69 mN/m.
    with db.fluid("IAPWS-95", ["water"]) as w:
        tr = w.transport(300.0, 1e5)
        assert tr.ok.all()
        assert abs(tr.viscosity - 853.74e-6) < 0.01e-6
        assert abs(tr.thermal_conductivity - 0.6095) < 1e-4
        assert abs(tr.surface_tension - 71.69e-3) < 0.01e-3
        # Steam at 400 K and 1 bar.
        st = w.transport(400.0, 1e5, phase="vapor")
        assert st.status.tolist() == [Status.OK, Status.OK, Status.NO_VALUE]
        assert 1.3e-5 < st.viscosity < 1.4e-5
