# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "mpmath==1.3.0", "CoolProp==8.0.0"]
# ///
"""Write testdata/m5_iapws95.json: the references for M5a's IAPWS-95
(docs/m5-plan.md, criteria 1 to 3).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_iapws95.py [part ...]

Parts (all by default; a part not named keeps what the file has):
  coefficients  the iapws package's IAPWS-95 coefficients (its _constants
                and Fi0), to hold water/coefficients.v's transcription from
                the release against;
  derivatives   phi^r and its five derivatives at (tau, delta) points, by
                mpmath.diff at 50 digits of an implementation written here
                from the release's Eq. (6), with each derivative's scale, the
                sum of its terms' magnitudes, which sets what rounding allows;
  states        (T, P) states from the triple point to 2273.15 K and 1 Pa to
                1000 MPa: the iapws package's and CoolProp's densities, h
                and s (each on its own reference), and the iapws package's
                phase;
  saturation    the saturation line by both: p, rho', rho'', h', h'', s', s''.

It reproduces its output byte for byte on the same day.
"""
import datetime
import json
import math
import pathlib
import sys
import warnings

import CoolProp
import iapws
import mpmath
from CoolProp.CoolProp import PropsSI
from iapws import IAPWS95

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5_iapws95.json"
C = IAPWS95._constants
TC = 647.096
RHOC = 322.0


# --- phi^r in mpmath, from the release's Eq. (6) -------------------------------

def term_functions():
    """The 56 terms of phi^r, each a function of (tau, delta) in mpmath."""
    mp = mpmath
    out = []
    for n, d, t in zip(C["nr1"], C["d1"], C["t1"]):
        out.append(lambda tau, delta, n=n, d=d, t=t: mp.mpf(n) * delta ** d * tau ** mp.mpf(t))
    for n, c, d, t in zip(C["nr2"], C["c2"], C["d2"], C["t2"]):
        out.append(lambda tau, delta, n=n, c=c, d=d, t=t:
                   mp.mpf(n) * delta ** d * tau ** t * mp.exp(-delta ** c))
    for n, d, t, a, b, g, e in zip(C["nr3"], C["d3"], C["t3"], C["alfa3"], C["beta3"],
                                   C["gamma3"], C["epsilon3"]):
        out.append(lambda tau, delta, n=n, d=d, t=t, a=a, b=b, g=g, e=e:
                   mp.mpf(n) * delta ** d * tau ** t
                   * mp.exp(-a * (delta - e) ** 2 - b * (tau - mp.mpf(g)) ** 2))
    for n, a, b, A, B, Cc, D, bt in zip(C["nr4"], C["a4"], C["b4"], C["A"], C["B"], C["C"],
                                         C["D"], C["beta4"]):
        def nonanalytic(tau, delta, n=n, a=a, b=b, A=A, B=B, Cc=Cc, D=D, bt=bt):
            d2 = (delta - 1) ** 2
            theta = (1 - tau) + mp.mpf(A) * d2 ** (1 / (2 * mp.mpf(bt)))
            big = theta ** 2 + mp.mpf(B) * d2 ** mp.mpf(a)
            return mp.mpf(n) * big ** mp.mpf(b) * delta * mp.exp(-Cc * d2 - D * (tau - 1) ** 2)
        out.append(nonanalytic)
    return out


TERMS = term_functions()
ORDERS = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (1, 1)]  # (tau, delta): f, d, dd, t, tt, dt


def derivatives(tau, delta):
    """phi^r's six derivatives and their scales at (tau, delta): the
    derivatives at 50 digits, the scales (a magnitude) at 20."""
    tm, dm = mpmath.mpf(tau), mpmath.mpf(delta)
    values, scales = [], []
    for ot, od in ORDERS:
        mpmath.mp.dps = 50
        v = mpmath.diff(lambda a, b: mpmath.fsum(f(a, b) for f in TERMS), (tm, dm), (ot, od))
        values.append(float(v))
        mpmath.mp.dps = 20
        scales.append(float(mpmath.fsum(abs(mpmath.diff(f, (tm, dm), (ot, od))) for f in TERMS)))
    mpmath.mp.dps = 50
    return values, scales


def points():
    """(T, delta) over the range and the release's extrapolations, and rays
    into the critical point."""
    out = []
    for t in [130.0, 200.0, 250.0]:  # the gas below the triple point, to 130 K
        for d in [1e-8, 1e-6, 1e-4]:
            out.append((TC / t, d))
    for t in [251.165, 273.16, 300.0, 373.15, 450.0, 550.0, 623.15, 640.0, 700.0, 800.0, 1000.0,
              1273.0, 2000.0, 5000.0]:
        for d in [1e-8, 1e-5, 1e-3, 0.01, 0.1, 0.3, 0.6, 1.0, 1.3, 1.8, 2.3, 2.8, 3.1, 3.4, 3.7, 3.95]:
            out.append((TC / t, d))
    for dt, dd in [(0.0, 1.0), (1.0, 0.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 1.0), (-1.0, 0.5)]:
        for k in range(2, 11):
            s = 10.0 ** -k
            out.append((1.0 + dt * s, 1.0 + dd * s))
    # delta = 1 exactly, Table 5's 0/0, away from tau = 1.
    for tau in [0.5, 0.9, 1.1, 2.0]:
        out.append((tau, 1.0))
    return out


def part_coefficients():
    keys = ["nr1", "d1", "t1", "nr2", "c2", "d2", "t2", "nr3", "d3", "t3", "alfa3", "beta3",
            "gamma3", "epsilon3", "nr4", "a4", "b4", "A", "B", "C", "D", "beta4"]
    return {
        "residual": {k: [float(x) for x in C[k]] for k in keys},
        "ideal": {k: [float(x) for x in IAPWS95.Fi0[k]] for k in ["ao_log", "ao_pow", "ao_exp",
                                                                  "titao"]},
        "R_molar": C["R"],
        "M": IAPWS95.M,
    }


def part_derivatives():
    rows = []
    for tau, delta in points():
        v, s = derivatives(tau, delta)
        rows.append({"tau": tau, "delta": delta, "values": v, "scales": s})
    return rows


def coolprop(t, p):
    """CoolProp's state at (T, P), on its own reference, or why there's none.
    Every row has every key, so that a test can decode it as one struct."""
    try:
        return {"ok": True, "rho": PropsSI("Dmolar", "T", t, "P", p, "Water"),
                "h": PropsSI("Hmolar", "T", t, "P", p, "Water"),
                "s": PropsSI("Smolar", "T", t, "P", p, "Water"), "error": ""}
    except ValueError as e:
        return {"ok": False, "rho": 0.0, "h": 0.0, "s": 0.0, "error": str(e).splitlines()[0][:160]}


def part_states():
    rows = []
    temps = [273.16, 280.0, 300.0, 350.0, 373.15, 400.0, 450.0, 500.0, 550.0, 600.0, 623.15, 640.0,
             645.0, 647.0, 647.096, 650.0, 660.0, 700.0, 800.0, 900.0, 1073.15, 1273.0, 1500.0,
             2000.0, 2273.15]
    pressures = [10 ** (k / 4) for k in range(0, 37)]  # 1 Pa to 1 GPa
    for t in temps:
        for p in pressures:
            w = IAPWS95(T=t, P=p / 1e6)
            ok = getattr(w, "status", 1) == 1 and w.rho is not None
            row = {"T": t, "P": p}
            m = IAPWS95.M / 1000
            if ok:
                row["iapws"] = {"ok": True, "rho": w.rho / m, "h": w.h * 1e3 * m,
                                "s": w.s * 1e3 * m, "phase": w.phase}
            else:
                row["iapws"] = {"ok": False, "rho": 0.0, "h": 0.0, "s": 0.0, "phase": ""}
            row["coolprop"] = coolprop(t, p)
            rows.append(row)
    return rows


def part_saturation():
    rows = []
    temps = [273.16, 275.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 625.0, 640.0, 645.0,
             646.0, 646.5, 647.0, 647.05, 647.09, 647.095]
    m = IAPWS95.M / 1000
    for t in temps:
        liq, vap = IAPWS95(T=t, x=0), IAPWS95(T=t, x=1)
        row = {"T": t,
               "iapws": {"ok": True, "p": liq.P * 1e6, "rho_l": liq.rho / m, "rho_v": vap.rho / m,
                         "h_l": liq.h * 1e3 * m, "h_v": vap.h * 1e3 * m,
                         "s_l": liq.s * 1e3 * m, "s_v": vap.s * 1e3 * m, "error": ""}}
        try:
            row["coolprop"] = {
                "ok": True,
                "p": PropsSI("P", "T", t, "Q", 0, "Water"),
                "rho_l": PropsSI("Dmolar", "T", t, "Q", 0, "Water"),
                "rho_v": PropsSI("Dmolar", "T", t, "Q", 1, "Water"),
                "h_l": PropsSI("Hmolar", "T", t, "Q", 0, "Water"),
                "h_v": PropsSI("Hmolar", "T", t, "Q", 1, "Water"),
                "s_l": PropsSI("Smolar", "T", t, "Q", 0, "Water"),
                "s_v": PropsSI("Smolar", "T", t, "Q", 1, "Water"), "error": ""}
        except ValueError as e:
            row["coolprop"] = {"ok": False, "p": 0.0, "rho_l": 0.0, "rho_v": 0.0, "h_l": 0.0,
                               "h_v": 0.0, "s_l": 0.0, "s_v": 0.0,
                               "error": str(e).splitlines()[0][:160]}
        rows.append(row)
    return rows


PARTS = {"coefficients": part_coefficients, "derivatives": part_derivatives,
         "states": part_states, "saturation": part_saturation}


def main():
    warnings.simplefilter("ignore")
    names = sys.argv[1:] or list(PARTS)
    for n in names:
        if n not in PARTS:
            sys.exit(f"no part {n}: {', '.join(PARTS)}")
    old = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    out = {
        "source": {"generator": "tools/gen_m5_iapws95.py", "iapws": iapws.__version__,
                   "coolprop": CoolProp.__version__, "mpmath": mpmath.__version__,
                   "generated": datetime.date.today().isoformat(),
                   "release": "IAPWS R6-95(2018)"},
        "note": ("Per mole, SI. derivatives: phi^r's f, d, dd, t, tt, dt (tau, delta) by "
                 "mpmath.diff at 50 digits; scales are the sums of the terms' |derivative|. "
                 "states and saturation: each source on its own reference for h and s; the "
                 "iapws package's are IAPWS's (u = s = 0 for the triple point's liquid)."),
    }
    for n in PARTS:
        out[n] = PARTS[n]() if n in names else old.get(n)
    OUT.write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {', '.join(names)}")


if __name__ == "__main__":
    main()
