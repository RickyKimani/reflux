# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2", "mpmath==1.3.0", "scipy==1.14.1"]
# ///
"""Write M6's fit references (step 11, criterion 9): testdata/m6_fits.json and
the data in testdata/m6_fit/.

- methanol_water.csv, a published set: methanol-water's bubble points at
  101.32 kPa (Yang, Ma and Yin, J. Chem. Eng. Data 56 (2011) 3747, from the
  ThermoML archive, as testdata/m6_published.json has them; its 13 points
  inside (0, 1), the numbers as printed). It determines NRTL's two energies,
  but not alpha besides, nor UNIQUAC's two, each nearly a combination of
  the other (correlations of -1.0);
- ethanol_water_nrtl.csv and ethanol_water_uniquac.csv, drawn from thermo's
  model at ChemSep's pair, ethanol-water at 323.15, 343.15 and 363.15 K and
  x1 = 0.05 to 0.95, P with 0.5 % and y with 0.003 of seeded noise, as M8's
  kij data were: more non-ideal, over three T, where alpha and UNIQUAC's
  energies are determined.

SciPy fits them as reflux's fit does
(M6's D16, M8's D15): each bubble point's P and y less the model's, over
its weights (1 % of P and 0.005 in y for the published set, the noise's
own for the drawn), by least_squares (lm, tolerances 1e-15) with the
Jacobian by central differences of reflux's steps (an energy +- 1e-4 of
itself or of 2500 J/mol, alpha +- 1e-4), then s^2 (J^T J)^-1 for the
standard errors and Student's t by mpmath for the 95 % intervals. A step
to parameters with no bubble point gets residuals of 1e10, which lm
refuses.

The model is thermo's: its NRTL (tau = A/RT, A in J/mol, R the SI's exact
value) or UNIQUAC (tau = exp(-A/RT)) in GibbsExcessLiquid on the 'Poynting'
basis, each compound's ChemSep vapour pressure and liquid density, the
vapour the ideal gas; each bubble point from its fugacities, P = P sum x_i
phi_i(T, P, x) by substitution to 1e-14, y_i = x_i phi_i. Three fits: NRTL's
two energies with alpha at the bank's 0.2999, then with alpha too, and
UNIQUAC's two from zero energies. NRTL's with alpha starts from the bank's
pair, as a file does: from zero energies and alpha 0.3 lm finds a minimum at
alpha = -0.59, s^2 = 4.2, where the noise gives about 1.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_fits.py

It reproduces its output byte for byte on the same day.
"""
import datetime
import json
import math
import pathlib
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m6_activity as ga  # noqa: E402
import gen_m6_gamma_phi as gp  # noqa: E402
import gen_m8_fits as m8  # noqa: E402

from thermo import GibbsExcessLiquid  # noqa: E402
from thermo.nrtl import NRTL  # noqa: E402
from thermo.uniquac import UNIQUAC  # noqa: E402

ROOT = ga.ROOT
R = ga.R_THERMO  # the SI's exact R, as reflux's expressions have it
METHANOL_WATER = ["67-56-1", "7732-18-5"]
ETHANOL_WATER = ["64-17-5", "7732-18-5"]
DOI = "10.1021/je200341c"
# (name, model, the data's file, compounds, parameters, P's weight relative,
# y's weight, the start: the bank's pair, or this)
FITS = [
    ("nrtl", "NRTL", "methanol_water.csv", METHANOL_WATER, ["A[methanol, water]", "A[water, methanol]"],
     0.01, 0.005, None),
    ("nrtl_alpha", "NRTL", "ethanol_water_nrtl.csv", ETHANOL_WATER,
     ["A[ethanol, water]", "A[water, ethanol]", "alpha[ethanol, water]"], 0.005, 0.003, None),
    ("uniquac", "UNIQUAC", "ethanol_water_uniquac.csv", ETHANOL_WATER, ["A[ethanol, water]", "A[water, ethanol]"],
     0.005, 0.003, [0.0, 0.0]),
]


def data():
    pub = json.loads((ROOT / "testdata" / "m6_published.json").read_text(encoding="utf-8"))
    s = [v for v in pub["vle"] if v["doi"] == DOI][0]
    assert s["point_test"]["consistent"]
    pts = [pt for pt in s["points"] if 0.0 < float(pt["x1"]) < 1.0]
    return s["P_kPa"], pts


def liquid(model, theta, cas):
    recs = [ga.record(c) for c in cas]
    a12, a21 = theta[0], theta[1]
    if model == "NRTL":
        alpha = theta[2] if len(theta) > 2 else 0.2999
        ge = NRTL(xs=[0.5, 0.5], T=330.0, tau_bs=[[0.0, a12 / R], [a21 / R, 0.0]],
                  alpha_cs=[[0.0, alpha], [alpha, 0.0]])
    else:
        ge = UNIQUAC(xs=[0.5, 0.5], T=330.0, rs=[r["parameters"]["uniquac_r"] for r in recs],
                     qs=[r["parameters"]["uniquac_q"] for r in recs], tau_bs=[[0.0, -a12 / R], [-a21 / R, 0.0]])
    return GibbsExcessLiquid(VaporPressures=[gp.vapour_pressure(r) for r in recs],
                             VolumeLiquids=[gp.liquid_volume(r) for r in recs],
                             HeatCapacityGases=[gp.heat_capacity(r) for r in recs], GibbsExcessModel=ge,
                             equilibrium_basis="Poynting", caloric_basis="Poynting", T=330.0, P=1e5, zs=[0.5, 0.5])


def bubble(liq, t, x):
    """The bubble point at T and x on the ideal-gas vapour: P = P sum x_i
    phi_i, which only the Poynting term moves with P, by substitution."""
    p = 1e5
    for _ in range(200):
        ph = [math.exp(v) for v in liq.to(T=t, P=p, zs=x).lnphis()]
        s = sum(xi * f for xi, f in zip(x, ph))
        p_new = p * s
        # thermo's fugacities are good to about 1e-15, so s settles within
        # 1e-14 of 1, and P with it.
        if abs(s - 1.0) < 1e-14:
            return p_new, [xi * f for xi, f in zip(x, ph)]
        p = p_new
    raise RuntimeError(f"no bubble point at T = {t}, x = {x}")


def step(name, v):
    return 1e-4 * max(abs(v), 2500.0) if name.startswith("A[") else 1e-4


def bank_pair(model, cas):
    """The bank's two energies, A in cal/mol over DECHEMA's R as J/mol, and
    NRTL's alpha."""
    b = ga.bank("nrtl.ipd" if model == "NRTL" else "uniquac.ipd", 3 if model == "NRTL" else 2)
    v = ga.oriented(b, cas[0], cas[1])
    return [v[0] / ga.R_DECHEMA * R, v[1] / ga.R_DECHEMA * R] + ([v[2]] if model == "NRTL" else [])


def drawn(model, cas, seed):
    """Bubble points of thermo's model at the bank's pair, with noise."""
    rng = np.random.default_rng(seed)
    liq = liquid(model, bank_pair(model, cas), cas)
    rows = []
    for t in (323.15, 343.15, 363.15):
        for i in range(10):
            x1 = 0.05 + 0.1 * i
            p, y = bubble(liq, t, [x1, 1.0 - x1])
            rows.append([t, p * (1 + 0.005 * rng.standard_normal()), x1, float(y[0] + 0.003 * rng.standard_normal())])
    return rows


def fit(model, cas, names, rows, p_weight, y_weight, start):
    def residuals(theta):
        liq = liquid(model, theta, cas)
        out = []
        for t, p_meas, x1, y_meas in rows:
            try:
                p, y = bubble(liq, t, [x1, 1.0 - x1])
            except (RuntimeError, OverflowError, ValueError):
                return np.full(2 * len(rows), 1e10)
            out += [(p - p_meas) / (p_weight * p_meas), (y[0] - y_meas) / y_weight]
        return np.array(out)

    def jacobian(theta):
        cols = []
        for q, name in enumerate(names):
            h = step(name, theta[q])
            up, down = np.array(theta, dtype=float), np.array(theta, dtype=float)
            up[q] += h
            down[q] -= h
            cols.append((residuals(up) - residuals(down)) / (up[q] - down[q]))
        return np.array(cols).T

    if start is None:
        start = bank_pair(model, cas)[:len(names)]
    res = optimize.least_squares(residuals, start, jac=jacobian, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    theta = [float(t) for t in res.x]
    r = residuals(theta)
    j = jacobian(theta)
    n, k = len(r), len(theta)
    ssr = float(np.sum(r ** 2))
    cov = ssr / (n - k) * np.linalg.inv(j.T @ j)
    se = [float(math.sqrt(cov[q, q])) for q in range(k)]
    t = m8.t_quantile(0.975, n - k)
    return {"model": model, "parameters": names, "start": start, "value": theta, "standard_error": se,
            "lower": [theta[q] - t * se[q] for q in range(k)], "upper": [theta[q] + t * se[q] for q in range(k)],
            "ssr": ssr, "residuals": n, "dof": n - k, "t": t, "scipy_nfev": int(res.nfev)}


def write_csv(path, alcohol, rows):
    lines = [f"T [K],P [Pa],x[{alcohol}],y[{alcohol}]"] + [",".join(repr(float(v)) for v in r) for r in rows]
    path.write_bytes(("\n".join(lines) + "\n").encode())


def main():
    p_kpa, pts = data()
    out_dir = ROOT / "testdata" / "m6_fit"
    out_dir.mkdir(exist_ok=True)
    lines = ["T [K],P [kPa],x[methanol],y[methanol]"] + [f"{pt['T']},{p_kpa},{pt['x1']},{pt['y1']}" for pt in pts]
    (out_dir / "methanol_water.csv").write_bytes(("\n".join(lines) + "\n").encode())
    rows = {"methanol_water.csv": [[float(pt["T"]), float(p_kpa) * 1e3, float(pt["x1"]), float(pt["y1"])] for pt in pts],
            "ethanol_water_nrtl.csv": drawn("NRTL", ETHANOL_WATER, 20261124),
            "ethanol_water_uniquac.csv": drawn("UNIQUAC", ETHANOL_WATER, 20261125)}
    write_csv(out_dir / "ethanol_water_nrtl.csv", "ethanol", rows["ethanol_water_nrtl.csv"])
    write_csv(out_dir / "ethanol_water_uniquac.csv", "ethanol", rows["ethanol_water_uniquac.csv"])
    out = {
        "source": {"generator": "tools/gen_m6_fits.py", "thermo": "0.6.1", "scipy": "1.14.1",
                   "generated": datetime.date.today().isoformat()},
        "note": ("SciPy's fits of testdata/m6_fit's data (methanol_water.csv, je200341c's 13 points at 101.32 kPa; "
                 "the ethanol-water sets drawn from thermo's models at the bank's pair, given as 'true') with "
                 "thermo's NRTL and UNIQUAC on the Poynting basis and the ideal gas, as reflux's fit makes them: "
                 "central differences, s^2 (J^T J)^-1, Student's t by mpmath. Energies in J/mol."),
        "data": {"doi": DOI, "P_kPa": p_kpa, "points": len(pts)},
        "fits": {},
    }
    for name, model, csv, cas, names, pw, yw, start in FITS:
        out["fits"][name] = fit(model, cas, names, rows[csv], pw, yw, start)
        out["fits"][name]["data"] = csv
        if csv != "methanol_water.csv":
            out["fits"][name]["true"] = bank_pair(model, cas)[:len(names)]
        f = out["fits"][name]
        print(name, [round(v, 6) for v in f["value"]], "+-", [round(v, 6) for v in f["standard_error"]],
              "s^2", round(f["ssr"] / f["dof"], 4), "nfev", f["scipy_nfev"])
    (ROOT / "testdata" / "m6_fits.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
