# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath==1.3.0"]
# ///
"""M8 step 8's reference (criterion 6): an isothermal first-order PFR's
conversion and its derivatives, at 50 digits.

n-butane to isobutane in a PFR at its inlet's T, rate k(T)·p of n-butane per
volume, k = A·exp(−Eₐ/RT): the moles don't change, so dF_A/dV = −kP·F_A/F and

    X = 1 − exp(−a),  a = k(T)·P·V/F,
    dX/dV = (kP/F)·exp(−a),
    dX/dT = exp(−a)·(PV/F)·k·Eₐ/(RT²).

Each closed-form derivative is checked here against mpmath's own numerical
derivative of X, to 30 digits, before it's written. R is the 2019 SI's,
N_A·k_B, as reflux's built-in constant.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_sensitivities.py

writes testdata/m8_sensitivities.json.
"""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
R = mp.mpf("1.380649e-23") * mp.mpf("6.02214076e23")
A = mp.mpf(70)  # mol/(m³·s·Pa)
EA = mp.mpf(50000)  # J/mol
P = mp.mpf(100000)  # Pa
F = mp.mpf(1)  # mol/s


def conversion(t, v):
    k = A * mp.exp(-EA / (R * t))
    return 1 - mp.exp(-k * P * v / F)


def case(t, v):
    t, v = mp.mpf(t), mp.mpf(v)
    k = A * mp.exp(-EA / (R * t))
    a = k * P * v / F
    x = 1 - mp.exp(-a)
    dxdv = (k * P / F) * mp.exp(-a)
    dxdt = mp.exp(-a) * (P * v / F) * k * EA / (R * t**2)
    assert abs(dxdv - mp.diff(lambda w: conversion(t, w), v)) < mp.mpf(10) ** -30 * abs(dxdv)
    assert abs(dxdt - mp.diff(lambda s: conversion(s, v), t)) < mp.mpf(10) ** -30 * abs(dxdt)
    return {
        "T": float(t),
        "V": float(v),
        "X": float(x),
        "dX/dV": float(dxdv),
        "dX/dT": float(dxdt),
    }


def main():
    cases = [case(t, v) for t in ["380", "400", "420"] for v in ["0.5", "1", "2"]]
    out = {
        "provenance": "tools/gen_m8_sensitivities.py: mpmath 1.3.0 at 50 digits; an isothermal "
        "first-order PFR, n-butane to isobutane, k = A exp(-Ea/RT) per volume on n-butane's "
        "partial pressure, A = 70 mol/(m^3 s Pa), Ea = 50 kJ/mol, P = 1 bar, F = 1 mol/s of "
        "n-butane; each derivative checked against mpmath's numerical one to 1e-30",
        "cases": cases,
    }
    path = Path(__file__).resolve().parent.parent / "testdata" / "m8_sensitivities.json"
    path.write_bytes((json.dumps(out, indent=2) + "\n").encode())
    print(f"wrote {len(cases)} cases to {path}")


if __name__ == "__main__":
    main()
