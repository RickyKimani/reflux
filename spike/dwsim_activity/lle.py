# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""Do ChemSep's NRTL and UNIQUAC parameters for water-1-butanol predict two
liquids? M6 step 0, beside run.py's DWSIM flashes. For each model, the
mixture's Gibbs energy of mixing over x, and the split by the common
tangent where it isn't convex.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/dwsim_activity/lle.py
"""
import math

from thermo.nrtl import NRTL
from thermo.uniquac import UNIQUAC

R_CAL = 1.98721  # DECHEMA's R, cal/(mol K), which DWSIM's NRTL uses (run.py)
# ChemSep's banks, water first: nrtl.ipd and uniquac.ipd, cal/mol.
NRTL_A = (2633.6951, 504.0381, 0.4447)
UNIQUAC_A = (548.2453, 89.0444)
R_Q = ([0.92, 3.4543], [1.4, 3.052])  # water, 1-butanol: the component files' r and q


def nrtl(x, t):
    a12, a21, alpha = NRTL_A
    return NRTL(xs=x, T=t, tau_bs=[[0.0, a12 / R_CAL], [a21 / R_CAL, 0.0]],
                alpha_cs=[[0.0, alpha], [alpha, 0.0]]).gammas()


def uniquac(x, t):
    a12, a21 = UNIQUAC_A
    return UNIQUAC(xs=x, T=t, rs=R_Q[0], qs=R_Q[1],
                   tau_bs=[[0.0, -a12 / R_CAL], [-a21 / R_CAL, 0.0]]).gammas()


def g_mix(model, x1, t):
    x = [x1, 1.0 - x1]
    g = model(x, t)
    return sum(xi * (math.log(xi) + math.log(gi)) for xi, gi in zip(x, g))


def split(model, t, n=20000):
    """The two liquids as the longest edge of G of mixing's lower convex hull,
    the common tangent, on a grid of n; None where G is convex."""
    xs = [i / n for i in range(1, n)]
    pts = [(x, g_mix(model, x, t)) for x in xs]
    hull = []
    for p in pts:  # Andrew's monotone chain, the lower hull
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    gaps = [(hull[i + 1][0] - hull[i][0], hull[i][0], hull[i + 1][0]) for i in range(len(hull) - 1)]
    width, xb, xa = max(gaps)
    if width < 10.0 / n:
        return None
    # Equal activities at the two ends, as a check on the grid's answer.
    ga, gb = model([xa, 1 - xa], t), model([xb, 1 - xb], t)
    e = max(abs(math.log(xa * ga[0] / (xb * gb[0]))), abs(math.log((1 - xa) * ga[1] / ((1 - xb) * gb[1]))))
    return e, xa, xb


for label, model in [("NRTL", nrtl), ("UNIQUAC", uniquac)]:
    for t in [298.15, 333.15, 360.0]:
        s = split(model, t)
        if s is None:
            print(f"{label} at {t} K: one liquid, G of mixing convex")
        else:
            e, xa, xb = s
            print(f"{label} at {t} K: two liquids, water {xa:.4f} and {xb:.4f} (activity gap {e:.1e} on a 1/2000 grid)")
print("Experiment at 298 K (Sorensen and Arlt's DECHEMA LLE volume): about 0.98 water in the aqueous phase,"
      " 0.51 in the organic.")
