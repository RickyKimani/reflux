# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = ["numpy==2.1.3", "scipy==1.14.1", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""M8 steps 11 and 12's references (criterion 9): data to fit, and SciPy's
fits of them, which reflux's must match.

1. Rate data. Ethylene's hydrogenation on a nickel site, LHHW,

       r1 = k1 p_E p_H / (1 + K_E p_E + K_H p_H)^2,

   and ethane's hydrogenolysis, a power law, r2 = k2 p_A p_H^0.5, on the ideal
   gas (p_i = y_i P, bar), k = A exp(-Ea/RT), K_E = A exp(-dH/RT), K_H fixed:
   the constants in TRUE, at 40 states (five T from 500 to 600 K, eight
   compositions and pressures), each rate with 5% Gaussian noise, seeded.
   testdata/m8_fit/ethylene.csv is the data.

   SciPy fits the seven constants by least_squares (lm, tolerances 1e-15),
   each rate weighted by 5% of its measurement, the Arrhenius and van 't
   Hoff pairs centred on the data's mean T as reflux fits them. The
   covariance is s^2 (J^T J)^-1 with J taken again in the file's own form,
   A and Ea, by mpmath's derivatives at 30 digits: not by reflux's
   transformation, which it checks. Each interval is the value plus or
   minus t(0.975, n - p) standard errors.

2. VLE data. Ethane and propane at 250 K, nine liquid compositions, each
   bubble point by thermo's Peng-Robinson (PRMIX) on reflux's own
   constants (data/components) with k_ij = 0.012, solved here to 1e-14 by
   successive substitution on thermo's ln phi, then P with 0.5% noise and
   y with 0.003, seeded. testdata/m8_fit/ethane_propane.csv is the data.
   SciPy fits k_ij on the same bubble points, P weighted by 0.5% of its
   measurement and y by 0.003.

3. Student's t quantiles, for the intervals, by mpmath at 40 digits.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_fits.py

writes the two CSV files and testdata/m8_fits.json.
"""
import json
import tomllib
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import optimize, stats
from thermo.eos_mix import PRMIX

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m8_fits.json"
DATA = ROOT / "testdata" / "m8_fit"
R = 1.380649e-23 * 6.02214076e23

# The rates' true constants, in the file's units: A in mol/(kg s bar^2) and
# mol/(kg s bar^1.5), Ea and dH in J/mol, K_E's A and K_H in 1/bar.
TRUE = {"k1.A": 1e5, "k1.Ea": 60000.0, "K_E.A": 1e-4, "K_E.dH": -50000.0, "K_H": 0.5,
        "k2.A": 2e4, "k2.Ea": 80000.0}
NAMES = list(TRUE)
# Where reflux's file starts the fit, as its [reaction] and [site] give them.
START = {"k1.A": 3e5, "k1.Ea": 55000.0, "K_E.A": 3e-4, "K_E.dH": -45000.0, "K_H": 0.3,
         "k2.A": 1e4, "k2.Ea": 75000.0}
SPECIES = ["ethylene", "hydrogen", "ethane", "methane"]


def rates(c, t, p, y):
    """The two rates at T (K), P (bar) and y, from constants in the file's form."""
    pe, ph, pa = y[0] * p, y[1] * p, y[2] * p
    k1 = c["k1.A"] * np.exp(-c["k1.Ea"] / (R * t))
    ke = c["K_E.A"] * np.exp(-c["K_E.dH"] / (R * t))
    k2 = c["k2.A"] * np.exp(-c["k2.Ea"] / (R * t))
    r1 = k1 * pe * ph / (1 + ke * pe + c["K_H"] * ph) ** 2
    r2 = k2 * pa * ph ** 0.5
    return r1, r2


def rate_data():
    rng = np.random.default_rng(20261118)
    rows = []
    for t in [500.0, 525.0, 550.0, 575.0, 600.0]:
        for k in range(8):
            p = 1.0 + 0.5 * k  # bar
            raw = rng.uniform(0.1, 1.0, 4)
            y = raw / raw.sum()
            r1, r2 = rates(TRUE, t, p, y)
            rows.append([t, p, *y, r1 * (1 + 0.05 * rng.standard_normal()),
                         r2 * (1 + 0.05 * rng.standard_normal())])
    return rows


def centred(theta, t_ref):
    """The file's form from the centred: A = k_ref exp(Ea/(R T_ref))."""
    k1r, ea1, ker, dhe, kh, k2r, ea2 = theta
    return {"k1.A": k1r * np.exp(ea1 / (R * t_ref)), "k1.Ea": ea1,
            "K_E.A": ker * np.exp(dhe / (R * t_ref)), "K_E.dH": dhe, "K_H": kh,
            "k2.A": k2r * np.exp(ea2 / (R * t_ref)), "k2.Ea": ea2}


def fit_rates(rows):
    t_ref = sum(r[0] for r in rows) / len(rows)
    meas = np.array([[r[6], r[7]] for r in rows])

    def residuals(theta):
        c = centred(theta, t_ref)
        out = []
        for row, m in zip(rows, meas):
            r1, r2 = rates(c, row[0], row[1], row[2:6])
            out += [(r1 - m[0]) / (0.05 * m[0]), (r2 - m[1]) / (0.05 * m[1])]
        return np.array(out)

    s = START
    theta0 = [s["k1.A"] * np.exp(-s["k1.Ea"] / (R * t_ref)), s["k1.Ea"],
              s["K_E.A"] * np.exp(-s["K_E.dH"] / (R * t_ref)), s["K_E.dH"], s["K_H"],
              s["k2.A"] * np.exp(-s["k2.Ea"] / (R * t_ref)), s["k2.Ea"]]
    res = optimize.least_squares(residuals, theta0, method="lm", xtol=1e-15, ftol=1e-15,
                                 gtol=1e-15, x_scale="jac", max_nfev=100000)
    assert res.success, res.message
    best = centred(res.x, t_ref)
    n, p = 2 * len(rows), len(NAMES)
    ssr = float(np.sum(res.fun ** 2))
    s2 = ssr / (n - p)
    # J in the file's form, by mpmath at 30 digits.
    mp.mp.dps = 30
    jac = np.zeros((n, p))
    for k, name in enumerate(NAMES):
        def column(v, name=name):
            c = {q: mp.mpf(best[q]) for q in NAMES}
            c[name] = v
            out = []
            for row, m in zip(rows, meas):
                pe, ph, pa = (mp.mpf(row[2]) * row[1], mp.mpf(row[3]) * row[1],
                              mp.mpf(row[4]) * row[1])
                t = mp.mpf(row[0])
                k1 = c["k1.A"] * mp.exp(-c["k1.Ea"] / (R * t))
                ke = c["K_E.A"] * mp.exp(-c["K_E.dH"] / (R * t))
                k2 = c["k2.A"] * mp.exp(-c["k2.Ea"] / (R * t))
                r1 = k1 * pe * ph / (1 + ke * pe + c["K_H"] * ph) ** 2
                r2 = k2 * pa * mp.sqrt(ph)
                out += [r1 / (0.05 * m[0]), r2 / (0.05 * m[1])]
            return out
        h = mp.mpf(best[name]) * mp.mpf(10) ** -12
        up = column(mp.mpf(best[name]) + h)
        down = column(mp.mpf(best[name]) - h)
        jac[:, k] = [float((u - d) / (2 * h)) for u, d in zip(up, down)]
    cov = s2 * np.linalg.inv(jac.T @ jac)
    se = np.sqrt(np.diag(cov))
    q = float(stats.t.ppf(0.975, n - p))
    corr = cov / np.outer(se, se)
    return {
        "t_ref": t_ref,
        "names": NAMES,
        "true": [TRUE[k] for k in NAMES],
        "start": [START[k] for k in NAMES],
        "value": [float(best[k]) for k in NAMES],
        "standard_error": [float(v) for v in se],
        "lower": [float(best[k] - q * e) for k, e in zip(NAMES, se)],
        "upper": [float(best[k] + q * e) for k, e in zip(NAMES, se)],
        "correlation": [[float(v) for v in row] for row in corr],
        "ssr": ssr,
        "residuals": n,
        "dof": n - p,
        "t": q,
        "scipy_nfev": int(res.nfev),
    }


def constants(cas):
    data = tomllib.loads((ROOT / "data" / "components" / f"{cas}.toml").read_text())["constants"]
    return data["tc"], data["pc"], data["omega"]


ETHANE, PROPANE = constants("74-84-0"), constants("74-98-6")
TCS, PCS, WS = [ETHANE[0], PROPANE[0]], [ETHANE[1], PROPANE[1]], [ETHANE[2], PROPANE[2]]
K_TRUE = 0.012
T_VLE = 250.0


def ln_phis(t, p, z, kij, liquid):
    eos = PRMIX(T=t, P=p, Tcs=TCS, Pcs=PCS, omegas=WS, zs=list(z), kijs=[[0, kij], [kij, 0]])
    if liquid:
        zz = eos.Z_l if hasattr(eos, "Z_l") else eos.Z_g
    else:
        zz = eos.Z_g if hasattr(eos, "Z_g") else eos.Z_l
    return np.array(eos.fugacity_coefficients(zz))


def bubble(t, x, kij, p0=None, y0=None):
    """thermo's PR bubble point at T and x, by successive substitution on K."""
    p = p0 if p0 is not None else 10e5
    y = np.array(y0) if y0 is not None else np.array(x)
    for _ in range(10000):
        k = np.exp(ln_phis(t, p, x, kij, True) - ln_phis(t, p, y, kij, False))
        s = float(np.sum(k * x))
        y_new = k * x / s
        p_new = p * s
        if abs(s - 1) < 1e-14 and np.max(np.abs(y_new - y)) < 1e-14:
            return p_new, y_new
        p, y = p_new, y_new
    raise RuntimeError(f"no bubble point at x = {x}")


def vle_data():
    rng = np.random.default_rng(20261119)
    rows = []
    p, y = None, None
    for i in range(9):
        x1 = 0.1 * (i + 1)
        x = np.array([x1, 1 - x1])
        p, y = bubble(T_VLE, x, K_TRUE, p, y)
        rows.append([T_VLE, p * (1 + 0.005 * rng.standard_normal()), x1,
                     float(y[0] + 0.003 * rng.standard_normal())])
    return rows


def fit_vle(rows):
    def residuals(theta):
        out = []
        for t, p_meas, x1, y_meas in rows:
            p, y = bubble(t, np.array([x1, 1 - x1]), theta[0])
            out += [(p - p_meas) / (0.005 * p_meas), (y[0] - y_meas) / 0.003]
        return np.array(out)

    res = optimize.least_squares(residuals, [0.0], method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                                 diff_step=1e-6)
    assert res.success, res.message
    kij = float(res.x[0])
    n, p = 2 * len(rows), 1
    ssr = float(np.sum(res.fun ** 2))
    h = 1e-6
    jac = (residuals([kij + h]) - residuals([kij - h])) / (2 * h)
    se = float(np.sqrt(ssr / (n - p) / np.dot(jac, jac)))
    q = float(stats.t.ppf(0.975, n - p))
    return {"true": K_TRUE, "value": kij, "standard_error": se, "lower": kij - q * se,
            "upper": kij + q * se, "ssr": ssr, "residuals": n, "dof": n - p, "t": q,
            "constants": {"tc": TCS, "pc": PCS, "omega": WS}}


def t_quantile(p, nu):
    """Student's t quantile by mpmath at 40 digits: the root of the CDF,
    1 - I_x(nu/2, 1/2)/2 at x = nu/(nu + t^2), from SciPy's value. SciPy's own
    is off by 1.2e-12 at nu = 1, p = 0.6, where it's tan(pi/10) exactly."""
    mp.mp.dps = 40
    p, nu = mp.mpf(str(p)), mp.mpf(str(nu))

    def cdf(t):
        return 1 - mp.betainc(nu / 2, mp.mpf(1) / 2, 0, nu / (nu + t * t), regularized=True) / 2

    t = mp.findroot(lambda t: cdf(t) - p, mp.mpf(float(stats.t.ppf(float(p), float(nu)))))
    if nu == 1:
        assert abs(t - mp.tan(mp.pi * (p - mp.mpf(1) / 2))) < mp.mpf(10) ** -30
    return float(t)


def write_csv(path, head, rows):
    lines = [",".join(head)] + [",".join(repr(float(v)) for v in row) for row in rows]
    path.write_bytes(("\n".join(lines) + "\n").encode())


def main():
    DATA.mkdir(exist_ok=True)
    rows = rate_data()
    head = ["T [K]", "P [bar]"] + [f"y[{s}]" for s in SPECIES] + [
        "r[hydrogenation] [mol/(kg*s)]", "r[hydrogenolysis] [mol/(kg*s)]"]
    write_csv(DATA / "ethylene.csv", head, rows)
    vle = vle_data()
    write_csv(DATA / "ethane_propane.csv", ["T [K]", "P [Pa]", "x[ethane]", "y[ethane]"], vle)
    quantiles = []
    for p in [0.6, 0.9, 0.975, 0.995]:
        for nu in [1, 2, 3, 4, 5, 7, 10, 17, 33, 73, 100, 1000, 1e6]:
            quantiles.append({"p": p, "nu": nu, "t": t_quantile(p, nu)})
    out = {
        "provenance": "tools/gen_m8_fits.py: SciPy 1.14.1 (least_squares, lm, tolerances 1e-15; "
        "stats.t for the intervals), NumPy 2.1.3, mpmath 1.3.0 (the t quantiles at 40 digits), thermo 0.6.1 (PRMIX on reflux's constants); "
        "synthetic data seeded 20261118 (rates, 5% noise) and 20261119 (VLE, 0.5% in P, 0.003 "
        "in y)",
        "rates": fit_rates(rows),
        "vle": fit_vle(vle),
        "t_quantiles": quantiles,
    }
    OUT.write_bytes((json.dumps(out, indent=1) + "\n").encode())
    r, v = out["rates"], out["vle"]
    print(f"rates: {dict(zip(r['names'], r['value']))}, s^2 = {r['ssr'] / r['dof']:.3f}")
    print(f"vle: kij = {v['value']} +- {v['standard_error']}, s^2 = {v['ssr'] / v['dof']:.3f}")


if __name__ == "__main__":
    main()
