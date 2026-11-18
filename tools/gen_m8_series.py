# /// script
# requires-python = ">=3.11"
# dependencies = ["mpmath==1.3.0"]
# ///
"""M8 step 10's reference (criterion 8): series reactions in an isothermal
PFR, A -> B -> C, and three optima of B's outflow, at 50 digits.

n-pentane (A) to isopentane (B) to neopentane (C), each rate first order
in its reactant's partial pressure, k_i(T) p per volume, k_i = A_i exp(-E_i/RT).
The moles don't change, so with K_i = k_i P/F (per m^3) and V the volume,

    F_A = F exp(-K1 V),
    F_B = F K1/(K2 - K1) (exp(-K1 V) - exp(-K2 V)),
    F_C = F - F_A - F_B,

and B peaks at V* = ln(K2/K1)/(K2 - K1), where F_B/F = (K1/K2)^(K2/(K2 - K1))
(Levenspiel, Chemical Reaction Engineering, ch. 8). E1 > E2, so a higher T
favours the first reaction, and B's peak rises with T.

The three optima:
- interior: B over V in [0.1, 2] m^3 at 400 K, at V*;
- on a bound: B over V and T, T in [380, 420] K: T on its upper bound, V at
  V*(420 K), and the bound's multiplier, by the envelope theorem dB/dT at
  (V*, 420 K) with V held;
- on a constraint: B over V at 400 K with F_C <= 0.05 mol/s: V where
  F_C = 0.05, and the constraint's multiplier, dB/dV over dC/dV there, the
  rate B rises as the limit on C is relaxed.

Each is checked here: the interior optimum's derivative is zero, B's peak
rises with T over the whole range, and the constraint binds below V*.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_series.py

writes testdata/m8_series.json.
"""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50
R = mp.mpf("1.380649e-23") * mp.mpf("6.02214076e23")
A1, E1 = mp.mpf(70), mp.mpf(50000)  # mol/(m^3 s Pa), J/mol
A2, E2 = mp.mpf("0.09"), mp.mpf(30000)
P = mp.mpf(100000)  # Pa
F = mp.mpf(1)  # mol/s
C_MAX = mp.mpf("0.05")  # mol/s


def ks(t):
    return A1 * mp.exp(-E1 / (R * t)) * P / F, A2 * mp.exp(-E2 / (R * t)) * P / F


def flows(t, v):
    k1, k2 = ks(t)
    a = F * mp.exp(-k1 * v)
    b = F * k1 / (k2 - k1) * (mp.exp(-k1 * v) - mp.exp(-k2 * v))
    return a, b, F - a - b


def peak(t):
    k1, k2 = ks(t)
    return mp.log(k2 / k1) / (k2 - k1)


def main():
    t = mp.mpf(400)
    v_star = peak(t)
    b_star = flows(t, v_star)[1]
    k1, k2 = ks(t)
    assert abs(b_star - F * (k1 / k2) ** (k2 / (k2 - k1))) < mp.mpf(10) ** -40
    assert abs(mp.diff(lambda v: flows(t, v)[1], v_star)) < mp.mpf(10) ** -40
    assert mp.mpf("0.1") < v_star < 2
    interior = {"T": float(t), "V": float(v_star), "B": float(b_star)}

    # On T's upper bound: B's peak rises with T across [380, 420] K.
    t_hi = mp.mpf(420)
    best = lambda s: flows(s, peak(s))[1]
    for s in [380 + 2 * i for i in range(21)]:
        assert mp.diff(best, s) > 0
    v_hi = peak(t_hi)
    assert mp.mpf("0.1") < v_hi < 2
    dbdt = mp.diff(lambda s: flows(s, v_hi)[1], t_hi)
    assert abs(dbdt - mp.diff(best, t_hi)) < mp.mpf(10) ** -30  # the envelope
    bound = {"T": float(t_hi), "V": float(v_hi), "B": float(flows(t_hi, v_hi)[1]),
             "T_multiplier": float(dbdt)}

    # On the constraint: C reaches its limit before B peaks.
    assert flows(t, v_star)[2] > C_MAX
    v_c = mp.findroot(lambda v: flows(t, v)[2] - C_MAX, mp.mpf("0.2"))
    assert mp.mpf("0.1") < v_c < v_star
    db = mp.diff(lambda v: flows(t, v)[1], v_c)
    dc = mp.diff(lambda v: flows(t, v)[2], v_c)
    assert db > 0 and dc > 0
    constraint = {"T": float(t), "V": float(v_c), "B": float(flows(t, v_c)[1]),
                  "C_max": float(C_MAX), "multiplier": float(db / dc)}

    out = {
        "provenance": "tools/gen_m8_series.py: mpmath 1.3.0 at 50 digits; n-pentane to "
        "isopentane to neopentane in an isothermal PFR, each rate k_i p per volume, "
        "A1 = 70 mol/(m^3 s Pa), E1 = 50 kJ/mol, A2 = 0.09 mol/(m^3 s Pa), E2 = 30 kJ/mol, "
        "P = 1 bar, F = 1 mol/s of n-pentane; B's optimum inside V's bounds, on T's upper "
        "bound, and on F_C <= 0.05 mol/s",
        "interior": interior,
        "bound": bound,
        "constraint": constraint,
    }
    path = Path(__file__).resolve().parent.parent / "testdata" / "m8_series.json"
    path.write_bytes((json.dumps(out, indent=2) + "\n").encode())
    print(f"wrote {path}: {json.dumps(out['interior'])}, {json.dumps(out['bound'])}, "
          f"{json.dumps(out['constraint'])}")


if __name__ == "__main__":
    main()
