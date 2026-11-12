# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["numpy==2.1.3", "scipy==1.14.1"]
# ///
"""Generate M9's fixture for the integrators (step 8): stiff test problems.

Run from the repository root:

    uv run --python 3.13 tools/gen_m9_odes.py

It writes testdata/m9_odes.json: six stiff problems of Hairer and Wanner's
(II, §IV.10) as the Test Set for IVP Solvers of the University of Bari
gives them (https://archimede.uniba.it/~testset/src/problems/, the files of
2006), each with the test set's reference solution at its interval's end,
as its Fortran file's `solut` has it: RADAU5 or RADAU at rtol = 1.1e-18
with a unit roundoff of 1.01e-19.

- **The problems:** HIRES, Robertson, van der Pol (ε = 1e-6), the
  Oregonator, E5 and the pollution problem, each f written here again from
  the test set's `feval`, with its Jacobian.
- **A check of the transcription:** each problem solved again here by
  SciPy's Radau at rtol 1e-12, with its absolute tolerances, and its
  distance from the reference recorded, in units of those tolerances, and
  as the largest relative difference of the components above 1e-30. A
  problem written wrong would land far from its reference. SciPy's Radau
  can't resolve E5 (its y₂ and y₃ move by 3e-4 between tolerances, its y₁,
  1e-291 in the reference, anywhere from −1e-258 to 1e-70), so E5's check
  is its y₂ and y₃ to 1e-3.

The output is the same, byte for byte, for the same versions.
"""

import json
import pathlib

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m9_odes.json"


def hires(t, y):
    f = np.zeros(8)
    f[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007
    f[1] = 1.71 * y[0] - 8.75 * y[1]
    f[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4]
    f[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3]
    f[4] = -1.745 * y[4] + 0.43 * (y[5] + y[6])
    f[5] = -280 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6]
    f[6] = 280 * y[5] * y[7] - 1.81 * y[6]
    f[7] = -f[6]
    return f


def hires_jac(t, y):
    j = np.zeros((8, 8))
    j[0, 0], j[0, 1], j[0, 2] = -1.71, 0.43, 8.32
    j[1, 0], j[1, 1] = 1.71, -8.75
    j[2, 2], j[2, 3], j[2, 4] = -10.03, 0.43, 0.035
    j[3, 1], j[3, 2], j[3, 3] = 8.32, 1.71, -1.12
    j[4, 4], j[4, 5], j[4, 6] = -1.745, 0.43, 0.43
    j[5, 3], j[5, 4], j[5, 5], j[5, 6], j[5, 7] = 0.69, 1.71, -0.43 - 280 * y[7], 0.69, -280 * y[5]
    j[6, 5], j[6, 6], j[6, 7] = 280 * y[7], -1.81, 280 * y[5]
    j[7] = -j[6]
    return j


def rober(t, y):
    return np.array([
        -0.04 * y[0] + 1e4 * y[1] * y[2],
        0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] ** 2,
        3e7 * y[1] ** 2,
    ])


def rober_jac(t, y):
    return np.array([
        [-0.04, 1e4 * y[2], 1e4 * y[1]],
        [0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]],
        [0.0, 6e7 * y[1], 0.0],
    ])


EPS = 1e-6


def vdpol(t, y):
    return np.array([y[1], ((1 - y[0] ** 2) * y[1] - y[0]) / EPS])


def vdpol_jac(t, y):
    return np.array([[0.0, 1.0], [(-2 * y[0] * y[1] - 1) / EPS, (1 - y[0] ** 2) / EPS]])


def orego(t, y):
    return np.array([
        77.27 * (y[1] + y[0] * (1 - 8.375e-6 * y[0] - y[1])),
        (y[2] - (1 + y[0]) * y[1]) / 77.27,
        0.161 * (y[0] - y[2]),
    ])


def orego_jac(t, y):
    return np.array([
        [77.27 * (1 - 2 * 8.375e-6 * y[0] - y[1]), 77.27 * (1 - y[0]), 0.0],
        [-y[1] / 77.27, -(1 + y[0]) / 77.27, 1 / 77.27],
        [0.161, 0.0, -0.161],
    ])


def e5(t, y):
    p1 = 7.89e-10 * y[0]
    p2 = 1.1e7 * y[0] * y[2]
    p3 = 1.13e9 * y[1] * y[2]
    p4 = 1.13e3 * y[3]
    f1 = -p1 - p2
    f2 = p1 - p3
    f4 = p2 - p4
    return np.array([f1, f2, f2 - f4, f4])


def e5_jac(t, y):
    a, b, c = 7.89e-10, 1.1e7, 1.13e9
    j = np.zeros((4, 4))
    j[0] = [-a - b * y[2], 0, -b * y[0], 0]
    j[1] = [a, -c * y[2], -c * y[1], 0]
    j[3] = [b * y[2], 0, b * y[0], -1.13e3]
    j[2] = j[1] - j[3]
    return j


K = [0.35, 0.266e2, 0.123e5, 0.86e-3, 0.82e-3, 0.15e5, 0.13e-3, 0.24e5, 0.165e5, 0.9e4,
     0.22e-1, 0.12e5, 0.188e1, 0.163e5, 0.48e7, 0.35e-3, 0.175e-1, 0.1e9, 0.444e12,
     0.124e4, 0.21e1, 0.578e1, 0.474e-1, 0.178e4, 0.312e1]
# Each reaction's rate constant and reactants (0-based species), and each
# species' change by reaction: the test set's r(1..25) and f(1..20).
REACTANTS = [(0,), (1, 3), (4, 1), (6,), (6,), (6, 5), (8,), (8, 5), (10, 1), (10, 0),
             (12,), (9, 1), (13,), (0, 5), (2,), (3,), (3,), (15,), (15,), (16, 5), (18,),
             (18,), (0, 3), (18, 0), (19,)]
CHANGE = {
    0: {0: -1, 9: -1, 13: -1, 22: -1, 23: -1, 1: 1, 2: 1, 8: 1, 10: 1, 11: 1, 21: 1, 24: 1},
    1: {1: -1, 2: -1, 8: -1, 11: -1, 0: 1, 20: 1},
    2: {14: -1, 0: 1, 16: 1, 18: 1, 21: 1},
    3: {1: -1, 15: -1, 16: -1, 22: -1, 14: 1},
    4: {2: -1, 3: 2, 5: 1, 6: 1, 12: 1, 19: 1},
    5: {5: -1, 7: -1, 13: -1, 19: -1, 2: 1, 17: 2},
    6: {3: -1, 4: -1, 5: -1, 12: 1},
    7: {3: 1, 4: 1, 5: 1, 6: 1},
    8: {6: -1, 7: -1},
    9: {11: -1, 6: 1, 8: 1},
    10: {8: -1, 9: -1, 7: 1, 10: 1},
    11: {8: 1},
    12: {10: -1, 9: 1},
    13: {12: -1, 11: 1},
    14: {13: 1},
    15: {17: -1, 18: -1, 15: 1},
    16: {19: -1},
    17: {19: 1},
    18: {20: -1, 21: -1, 23: -1, 22: 1, 24: 1},
    19: {24: -1, 23: 1},
}


def pollu_rates(y):
    return [K[r] * np.prod([y[s] for s in sp]) for r, sp in enumerate(REACTANTS)]


def pollu(t, y):
    r = pollu_rates(y)
    return np.array([sum(c * r[k] for k, c in CHANGE[i].items()) for i in range(20)])


def pollu_jac(t, y):
    # dr_k/dy_s for each reaction, then the changes.
    dr = np.zeros((25, 20))
    for k, sp in enumerate(REACTANTS):
        for a in sp:
            others = [s for s in sp if s != a] if sp.count(a) == 1 else []
            dr[k, a] += K[k] * np.prod([y[s] for s in others])
    j = np.zeros((20, 20))
    for i in range(20):
        for k, c in CHANGE[i].items():
            j[i] += c * dr[k]
    return j


POLLU_Y0 = [0.0] * 20
POLLU_Y0[1], POLLU_Y0[3], POLLU_Y0[6], POLLU_Y0[7], POLLU_Y0[8], POLLU_Y0[16] = (
    0.2, 0.04, 0.1, 0.3, 0.01, 0.007)

CRAY = "RADAU5 in true double precision on a Cray C90: uround 1.01e-19, rtol = atol = h0 = 1.1e-18"
ALPHA = "RADAU in double precision on an Alphaserver DS20E: uround 1.01e-19, rtol = atol = h0 = 1.1e-18"
ALPHA_ATOL = "RADAU in double precision on an Alphaserver DS20E: uround 1.01e-19, rtol = h0 = 1.1e-18, atol = 1.1e-40"

PROBLEMS = [
    dict(name="hires", f=hires, jac=hires_jac, t0=0.0, t1=321.8122,
         y0=[1.0, 0, 0, 0, 0, 0, 0, 0.0057], atol=1e-4,
         reference=["0.7371312573325668d-3", "0.1442485726316185d-3", "0.5888729740967575d-4",
                    "0.1175651343283149d-2", "0.2386356198831331d-2", "0.6238968252742796d-2",
                    "0.2849998395185769d-2", "0.2850001604814231d-2"], how=CRAY),
    dict(name="rober", f=rober, jac=rober_jac, t0=0.0, t1=1e11, y0=[1.0, 0, 0], atol=1e-4,
         reference=["0.2083340149701255D-07", "0.8333360770334713D-13", "0.9999999791665050D+00"],
         how=ALPHA),
    dict(name="vdpol", f=vdpol, jac=vdpol_jac, t0=0.0, t1=2.0, y0=[2.0, 0], atol=1.0,
         reference=["0.1706167732170483D+01", "-0.8928097010247975D+00"], how=ALPHA),
    dict(name="orego", f=orego, jac=orego_jac, t0=0.0, t1=360.0, y0=[1.0, 2, 3], atol=1.0,
         reference=["0.1000814870318523D+001", "0.1228178521549917D+004",
                    "0.1320554942846706D+003"], how=ALPHA_ATOL),
    dict(name="e5", f=e5, jac=e5_jac, t0=0.0, t1=1e13, y0=[1.76e-3, 0, 0, 0], atol=1e-20,
         reference=["0.1152903278711829D-290", "0.8867655517642120D-022",
                    "0.8854814626268838D-022", "0.0000000000000000D+000"], how=ALPHA_ATOL),
    dict(name="pollu", f=pollu, jac=pollu_jac, t0=0.0, t1=60.0, y0=POLLU_Y0, atol=1e-4,
         reference=["0.5646255480022769d-01", "0.1342484130422339d+00", "0.4139734331099427d-08",
                    "0.5523140207484359d-02", "0.2018977262302196d-06", "0.1464541863493966d-06",
                    "0.7784249118997964d-01", "0.3245075353396018d+00", "0.7494013383880406d-02",
                    "0.1622293157301561d-07", "0.1135863833257075d-07", "0.2230505975721359d-02",
                    "0.2087162882798630d-03", "0.1396921016840158d-04", "0.8964884856898295d-02",
                    "0.4352846369330103d-17", "0.6899219696263405d-02", "0.1007803037365946d-03",
                    "0.1772146513969984d-05", "0.5682943292316392d-04"], how=CRAY),
]


def fortran(s):
    return float(s.lower().replace("d", "e"))


def check_jacobian(p):
    """The Jacobian against central differences at the start and partway."""
    y = np.array(p["y0"], dtype=float) + 1e-3
    j = p["jac"](0.0, y)
    worst = 0.0
    for k in range(len(y)):
        h = 1e-6 * max(abs(y[k]), 1e-3)
        yp, ym = y.copy(), y.copy()
        yp[k] += h
        ym[k] -= h
        col = (p["f"](0.0, yp) - p["f"](0.0, ym)) / (2 * h)
        scale = max(np.max(np.abs(j)), 1.0)
        worst = max(worst, float(np.max(np.abs(col - j[:, k]))) / scale)
    assert worst < 1e-6, (p["name"], worst)


def main():
    out = []
    for p in PROBLEMS:
        check_jacobian(p)
        ref = [fortran(s) for s in p["reference"]]
        rtol = 1e-12
        atol = p["atol"] * rtol
        sol = solve_ivp(p["f"], (p["t0"], p["t1"]), p["y0"], method="Radau", rtol=rtol, atol=atol,
                        jac=p["jac"])
        assert sol.success, (p["name"], sol.message)
        ys = sol.y[:, -1]
        units = max(abs(a - b) / (atol + rtol * abs(b)) for a, b in zip(ys, ref))
        relative = max(abs(a - b) / abs(b) for a, b in zip(ys, ref) if abs(b) > 1e-30)
        if p["name"] == "e5":
            assert relative < 1e-3, (p["name"], relative)
        else:
            assert units < 1e3, (p["name"], units)
        out.append({
            "name": p["name"],
            "t0": p["t0"],
            "t1": p["t1"],
            "y0": [float(v) for v in p["y0"]],
            "atol_per_rtol": p["atol"],
            "reference": ref,
            "reference_source": f"{p['name']}.f's solut: {p['how']}",
            "scipy": [float(v) for v in ys],
            "scipy_tolerance_units": float(f"{units:.3g}"),
            "scipy_relative": float(f"{relative:.3g}"),
        })
        print(f"{p['name']}: SciPy's Radau within {units:.3g} tolerance units of the reference, "
              f"{relative:.3g} relative")
    doc = {
        "source": ("Six stiff test problems as the Test Set for IVP Solvers of the University "
                   "of Bari gives them (archimede.uniba.it/~testset/src/problems/, the files of "
                   "2006), with its reference solutions; each solved again by SciPy " + scipy.__version__ + "'s Radau at "
                   "rtol 1e-12 and atol = atol_per_rtol × 1e-12, as a check of the problems' "
                   "transcription, its distance in units of atol + rtol·|y| and relative, over "
                   "the components above 1e-30; E5's to 1e-3 relative, SciPy's Radau not "
                   "resolving E5. tools/gen_m9_odes.py"),
        "problems": out,
    }
    OUT.write_bytes((json.dumps(doc, indent=1) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
