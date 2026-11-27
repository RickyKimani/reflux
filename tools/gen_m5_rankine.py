# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "numpy==2.2.6", "scipy==1.18.1"]
# ///
"""Write testdata/m5_rankine.json: M5a's Rankine cycle
(testdata/flowsheets/rankine.toml) worked by hand from the iapws package's
IAPWS-95 states, for criterion 8.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_rankine.py

Per kg, on IAPWS's reference, with reflux's unit models (M4's D4): the
pump's work v1 (P2 - P1)/eta, from the condensate's v; the boiler to T3 at
P2; the turbine's eta times its isentropic drop, h(P4, s3); the condenser
back to the condensate. A one-phase state is the package's IAPWS95 at the
pair given. A wet state is the lever rule on the package's saturated liquid
and vapour at its P (by `saturated`): the package's own (P, s) inside the dome can land off
the saturation line (at 10 kPa, s = 6.7266 kJ/(kg K) gave T 0.15 K above
Tsat and h 0.94 kJ/kg high; s = 4 and 7.5 gave 322 K). It reproduces its
output byte for byte on the same day.
"""
import datetime
import json
import pathlib
import warnings

import iapws
from iapws import IAPWS95

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5_rankine.json"

T1, P1 = 318.15, 0.01  # K, MPa: the condensate
P2 = 8.0  # MPa: the pump's outlet
T3 = 773.15  # K: the boiler's outlet
ETA_PUMP, ETA_TURBINE = 0.75, 0.85
FLOW = 50.0 * 0.018015268  # kg/s: 50 mol/s


def saturated(p):
    """The saturated liquid and vapour at p (MPa), at the package's own
    Psat(T) = p by Newton in T with Clapeyron's slope. Its IAPWS95(P=p, x=0)
    stops early: at 10 kPa its T is 1.5e-8 K high, P(T) 7.6e-10 over p;
    CoolProp's Tsat is within 1.4e-10 K of the Newton's."""
    liquid = IAPWS95(P=p, x=0.0)
    t = liquid.T
    for _ in range(3):
        liquid, vapour = IAPWS95(T=t, x=0.0), IAPWS95(T=t, x=1.0)
        slope = (vapour.h - liquid.h) / (t * (vapour.v - liquid.v)) * 1e-3  # MPa/K
        t -= (liquid.P - p) / slope
    return IAPWS95(T=t, x=0.0), IAPWS95(T=t, x=1.0)


def main():
    s1 = IAPWS95(T=T1, P=P1)
    w_pump = s1.v * (P2 - P1) * 1e3 / ETA_PUMP  # kJ/kg: m3/kg x MPa = MJ/kg
    s2 = IAPWS95(P=P2, h=s1.h + w_pump)
    s3 = IAPWS95(T=T3, P=P2)
    liquid, vapour = saturated(P1)
    x4s = (s3.s - liquid.s) / (vapour.s - liquid.s)
    h4s = liquid.h + x4s * (vapour.h - liquid.h)
    h4 = s3.h - ETA_TURBINE * (s3.h - h4s)
    x4 = (h4 - liquid.h) / (vapour.h - liquid.h)
    states = {}
    for name, st in (("condensate", s1), ("feedwater", s2), ("steam", s3)):
        states[name] = {"T": float(st.T), "P": float(st.P), "h": float(st.h), "s": float(st.s), "x": float(st.x)}
    states["exhaust"] = {
        "T": float(liquid.T),
        "P": P1,
        "h": float(h4),
        "s": float(liquid.s + x4 * (vapour.s - liquid.s)),
        "x": float(x4),
    }
    states["isentropic exhaust"] = {"T": float(liquid.T), "P": P1, "h": float(h4s), "s": float(s3.s), "x": float(x4s)}
    q_boiler = s3.h - s2.h
    w_turbine = s3.h - h4
    q_condenser = s1.h - h4
    out = {
        "provenance": {
            "script": "tools/gen_m5_rankine.py",
            "iapws": iapws.__version__,
            "date": datetime.date.today().isoformat(),
            "units": "K, MPa, kJ/kg, kJ/(kg K); powers kW at the flow",
        },
        "flow_kg_s": FLOW,
        "states": states,
        "per_kg": {
            "pump": float(w_pump),
            "boiler": float(q_boiler),
            "turbine": float(w_turbine),
            "condenser": float(q_condenser),
        },
        "efficiency": float((w_turbine - w_pump) / q_boiler),
    }
    OUT.write_bytes((json.dumps(out, indent=1, sort_keys=True) + "\n").encode("utf-8"))
    print(f"wrote {OUT}: efficiency {out['efficiency']:.6f}, x4 {x4:.6f}")


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
