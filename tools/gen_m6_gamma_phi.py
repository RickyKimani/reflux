# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2", "mpmath==1.3.0", "scipy==1.14.1"]
# ///
"""Write M6's gamma-phi references (steps 5 and 6, criteria 3, 5 and 7):
testdata/m6_gamma_phi.json.

- "liquid": thermo 0.6.1's GibbsExcessLiquid at seeded states, on ChemSep's
  own coefficients (each compound's vapour pressure and liquid density from
  its component file, the activity model from the banks and tables): ln phi,
  its T-, P- and n-derivatives, and H, S and Cp less the ideal gas's. Bases
  'Psat' and 'Poynting' with the ideal gas, and 'Poynting&PhiSat' with
  Peng-Robinson's pure phi at Psat.
- "henry": gases in water by Henry's law, ChemSep's henrys_law.ipd, basis
  'Psat' with the ideal solution (a solute's own gamma is 1 in reflux's
  convention, and water alone is ideal).
- "henry_mixed": CO2 and hydrogen in water and ethylene glycol, a pair no
  bank lists, so an ideal solution in reflux too ('ideal'). ChemSep's bank
  has water alone as a solvent, so each gas's constant in the glycol is
  synthetic, given here: water's with A less 0.33 and ln 8 (CO2) or ln 3
  (hydrogen), and D more 1e-3/K.
  thermo mixes them by its 'solvents_with_parameters' mode; it has no
  composition derivatives, so dn is Richardson's extrapolation of central
  differences of its ln phi.

The activity models are built as tools/gen_m6_activity.py builds them, the
banks' A over DECHEMA's R. The vapour pressure is thermo's DIPPR101 on
ChemSep's coefficients, and the liquid volume its DIPPR105 or 106 reciprocal
on them, A scaled to mol/m^3.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_gamma_phi.py

It reproduces its output byte for byte.
"""
import datetime
import json
import math
import pathlib
import random
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m6_activity as ga  # noqa: E402
import mpmath as mp  # noqa: E402

from thermo import (VaporPressure, VolumeLiquid, HeatCapacityGas, GibbsExcessLiquid,  # noqa: E402
                    ChemicalConstantsPackage, PropertyCorrelationsPackage, IdealGas, FlashVL)
from scipy.optimize import brentq  # noqa: E402
from thermo.activity import IdealSolution  # noqa: E402
from thermo.nrtl import NRTL  # noqa: E402
from thermo.uniquac import UNIQUAC  # noqa: E402
from thermo.unifac import UNIFAC  # noqa: E402
from thermo.eos import PR  # noqa: E402

ROOT = ga.ROOT
SEED = 7
R_DECHEMA = ga.R_DECHEMA

LIQUID_SYSTEMS = [
    ("NRTL", ["67-56-1", "7732-18-5"]),
    ("UNIQUAC", ["67-56-1", "7732-18-5"]),
    ("NRTL", ["64-17-5", "7732-18-5"]),
    ("NRTL", ["67-64-1", "67-56-1", "7732-18-5"]),
    ("UNIFAC", ["67-56-1", "7732-18-5"]),
]
HENRY_SYSTEMS = [["7732-18-5", "1333-74-0"], ["7732-18-5", "7727-37-9", "7782-44-7"]]


def correlation(rec, key):
    c = rec["correlations"][key]
    return c["eq"], c["coeffs"], c["t_range"][0], c["t_range"][1]


def vapour_pressure(rec):
    eq, k, tmin, tmax = correlation(rec, "vapor_pressure")
    assert eq == 101, (rec["name"], eq)
    vp = VaporPressure(Tc=rec["constants"]["tc"], Pc=rec["constants"]["pc"], omega=rec["constants"]["omega"])
    vp.add_correlation(name="CHEMSEP", model="DIPPR101", Tmin=tmin, Tmax=tmax, A=k[0], B=k[1], C=k[2], D=k[3], E=k[4])
    vp.method = "CHEMSEP"
    return vp


def liquid_volume(rec):
    eq, k, tmin, tmax = correlation(rec, "liquid_density")
    c = rec["constants"]
    vl = VolumeLiquid(Tc=c["tc"], Pc=c["pc"], omega=c["omega"])
    if eq == 105:
        vl.add_correlation(name="CHEMSEP", model="DIPPR105_reciprocal", Tmin=tmin, Tmax=tmax,
                           A=k[0] * 1e3, B=k[1], C=k[2], D=k[3])
    else:
        assert eq == 106, (rec["name"], eq)
        vl.add_correlation(name="CHEMSEP", model="DIPPR106_inv", Tmin=tmin, Tmax=tmax, Tc=c["tc"],
                           A=k[0] * 1e3, B=k[1], C=k[2], D=k[3], E=k[4])
    vl.method = "CHEMSEP"
    return vl


def mp_volume_d2(rec, t):
    """d^2V/dT^2 of the liquid volume, by mpmath at 40 digits from the
    correlation's closed form."""
    mp.mp.dps = 40
    eq, k, _, _ = correlation(rec, "liquid_density")
    a = [mp.mpf(v) for v in k]
    if eq == 105:
        def v(temp):
            return 1 / (a[0] * 1000 / a[1] ** (1 + (1 - temp / a[2]) ** a[3]))
    else:
        tc = mp.mpf(rec["constants"]["tc"])

        def v(temp):
            tr = temp / tc
            return 1 / (a[0] * 1000 * (1 - tr) ** (a[1] + tr * (a[2] + tr * (a[3] + tr * a[4]))))
    return float(mp.diff(v, mp.mpf(t), 2))


def cp_corrected(liq, recs, vls, vps, t, p, x):
    """thermo's Cp with its numerical d^2V/dT^2 (0.05% to 0.3% off, by mpmath)
    replaced by mpmath's: in Poynting's ln V(P - Psat)/RT the second
    T-derivative carries V''(P - Psat)/RT, so Cp moves by
    -T*sum x_i (V''_exact - V''_thermo)(P - Psat_i)."""
    cp = liq.Cp() - liq.Cp_ideal_gas()
    for xi, rec, vl, vp in zip(x, recs, vls, vps):
        d2 = mp_volume_d2(rec, t) - vl.T_dependent_property_derivative(t, order=2)
        cp -= t * xi * d2 * (p - vp.T_dependent_property(t))
    return cp


def placeholder_pressure():
    """A solute's vapour pressure, which thermo computes and Henry's law then replaces."""
    vp = VaporPressure()
    vp.add_correlation(name="PLACEHOLDER", model="constant", Tmin=10.0, Tmax=5000.0, value=1e5)
    vp.method = "PLACEHOLDER"
    return vp


def heat_capacity(rec):
    """Any ideal-gas Cp: the references are H, S and Cp less the ideal gas's."""
    cp = HeatCapacityGas()
    cp.add_correlation(name="FLAT", model="constant", Tmin=10.0, Tmax=5000.0, value=40.0)
    cp.method = "FLAT"
    return cp


def activity(name, cas, t, x):
    recs = [ga.record(c) for c in cas]
    nc = len(cas)
    if name == "UNIFAC":
        mains, subs_files, version, key = ga.TABLES["UNIFAC"]
        subs, params = ga.unifac_table(mains, subs_files)
        interaction = {m: {n: v[0] for n, v in row.items()} for m, row in params.items()}
        groups = [{int(k): v for k, v in r["groups"][key].items()} for r in recs]
        return UNIFAC.from_subgroups(T=t, xs=x, chemgroups=groups, subgroups=subs,
                                     interaction_data=interaction, version=version)
    b = ga.bank("nrtl.ipd" if name == "NRTL" else "uniquac.ipd", 3 if name == "NRTL" else 2)
    p = {(i, j): ga.oriented(b, cas[i], cas[j]) for i in range(nc) for j in range(nc) if i != j}
    if name == "NRTL":
        bs = [[0.0 if i == j else p[i, j][0] / R_DECHEMA for j in range(nc)] for i in range(nc)]
        al = [[0.0 if i == j else p[i, j][2] for j in range(nc)] for i in range(nc)]
        return NRTL(xs=x, T=t, tau_bs=bs, alpha_cs=al)
    bs = [[0.0 if i == j else -p[i, j][0] / R_DECHEMA for j in range(nc)] for i in range(nc)]
    return UNIQUAC(xs=x, T=t, rs=[r["parameters"]["uniquac_r"] for r in recs],
                   qs=[r["parameters"]["uniquac_q"] for r in recs], tau_bs=bs)


def pure_pr(rec, t, p):
    c = rec["constants"]
    return PR(Tc=c["tc"], Pc=c["pc"], omega=c["omega"], T=t, P=p)


def phi_sat_terms(liq, recs, vps, t):
    """For 'Poynting&PhiSat': thermo's ln phi_sat with its two T-derivatives,
    and Peng-Robinson's ln phi of the vapour at the correlation's Psat, with
    its total T-derivatives along Psat(T) by Richardson's extrapolation of
    central differences. thermo takes phi_sat at the equation's own vapour
    pressure, and reflux at the correlation's, which is what makes a pure
    liquid's fugacity its vapour's at Psat; the test moves thermo's answer
    from the one to the other by these terms."""
    phi, dphi, d2phi = liq.phis_sat(), liq.dphis_sat_dT(), liq.d2phis_sat_dT2()
    own = []
    for i in range(len(recs)):
        d = dphi[i] / phi[i]
        own.append([math.log(phi[i]), d, d2phi[i] / phi[i] - d * d])

    def g(rec, vp, temp):
        return pure_pr(rec, temp, vp.T_dependent_property(temp)).lnphi_g

    at = []
    for rec, vp in zip(recs, vps):
        h = 0.05
        f = {k: g(rec, vp, t + k * h) for k in (-2, -1, 0, 1, 2)}
        d1 = (f[1] - f[-1]) / (2 * h)
        d2 = (f[2] - f[-2]) / (4 * h)
        dd1 = (f[1] - 2 * f[0] + f[-1]) / (h * h)
        dd2 = (f[2] - 2 * f[0] + f[-2]) / (4 * h * h)
        at.append([f[0], (4 * d1 - d2) / 3, (4 * dd1 - dd2) / 3])
    return {"thermo": own, "at_psat": at}


def answer(liq, x, t, volume=True):
    n = len(x)
    # GibbsExcessLiquid has no composition derivatives of its own, and ln phi
    # depends on composition through ln gamma alone (a solute's H through
    # the solvents', and water is the one solvent here), so n*d(ln phi_i)/dn_j
    # is the activity model's, as gen_m6_activity.py takes it.
    hess = liq.GibbsExcessModel.d2nGE_dninjs()
    dn = [[hess[i][j] / (ga.R_THERMO * t) for j in range(n)] for i in range(n)]
    return {
        "ln_phi": list(liq.lnphis()),
        "dt": list(liq.dlnphis_dT()),
        "dp": list(liq.dlnphis_dP()),
        "dn": [dn[i][j] for i in range(n) for j in range(n)],
        "h": liq.H() - liq.H_ideal_gas(),
        "s": liq.S() - liq.S_ideal_gas(),
        "cp": liq.Cp() - liq.Cp_ideal_gas(),
        "v": liq.V() if volume else None,
    }


def liquid_cases(rng):
    cases = []
    for name, cas in LIQUID_SYSTEMS:
        recs = [ga.record(c) for c in cas]
        vps = [vapour_pressure(r) for r in recs]
        vls = [liquid_volume(r) for r in recs]
        cps = [heat_capacity(r) for r in recs]
        bases = ["Psat", "Poynting"] + (["Poynting&PhiSat"] if name == "NRTL" and len(cas) == 2 else [])
        for basis in bases:
            for _ in range(4):
                t = round(rng.uniform(300.0, 360.0), 2)
                p = round(10 ** rng.uniform(4.0, 6.5), 0)
                raw = [rng.random() for _ in cas]
                x = [v / sum(raw) for v in raw]
                x[-1] = 1.0 - sum(x[:-1])
                kw = {}
                if basis == "Poynting&PhiSat":
                    kw["eos_pure_instances"] = [pure_pr(r, t, p) for r in recs]
                liq = GibbsExcessLiquid(VaporPressures=vps, VolumeLiquids=vls, HeatCapacityGases=cps,
                                        GibbsExcessModel=activity(name, cas, t, x), equilibrium_basis=basis,
                                        caloric_basis=basis, T=t, P=p, zs=x, **kw)
                case = {"model": name, "basis": basis, "cas": cas, "T": t, "P": p, "x": x,
                        "thermo": answer(liq, x, t)}
                if basis != "Psat":
                    case["thermo"]["cp_thermo"] = case["thermo"]["cp"]
                    case["thermo"]["cp"] = cp_corrected(liq, recs, vls, vps, t, p, x)
                if basis == "Poynting&PhiSat":
                    case["phi_sat"] = phi_sat_terms(liq, recs, vps, t)
                cases.append(case)
    return cases


def henry_lines():
    """henrys_law.ipd's constants in water, by gas: A, B, C, D (E is 1 throughout), tmin, tmax."""
    out = {}
    for f in ga.ipd_lines("henrys_law.ipd"):
        if len(f) < 11 or f[1] != "7732-18-5":
            continue
        a, b, c, d, e = (float(v) for v in f[6:11])
        assert int(f[3]) == 121 and e == 1.0, f
        out[f[0]] = (a, b, c, d, float(f[4]), float(f[5]))
    return out


def henry_cases(rng):
    lines = henry_lines()
    cases = []
    for cas in HENRY_SYSTEMS:
        recs = [ga.record(c) for c in cas]
        nc = len(cas)
        henry = [False] + [True] * (nc - 1)
        # thermo takes each pair's (A, B, C, D, E, F) of ln H = A + B/T + C ln T
        # + D T + E/T^2 + F T^2, solute i in solvent j at [i][j].
        abcdef = [[(0.0,) * 6 for _ in range(nc)] for _ in range(nc)]
        for i in range(1, nc):
            a, b, c, d, _, _ = lines[cas[i]]
            abcdef[i][0] = (a, b, c, d, 0.0, 0.0)
        for _ in range(6):
            t = round(rng.uniform(290.0, 360.0), 2)
            p = round(10 ** rng.uniform(5.0, 6.5), 0)
            raw = [0.0] + [rng.random() * 1e-3 for _ in cas[1:]]
            raw[0] = 1.0 - sum(raw[1:])
            x = raw
            liq = GibbsExcessLiquid(VaporPressures=[vapour_pressure(recs[0])] + [placeholder_pressure() for _ in recs[1:]],
                                    VolumeLiquids=[liquid_volume(recs[0])] + [None] * (nc - 1),
                                    HeatCapacityGases=[heat_capacity(r) for r in recs],
                                    GibbsExcessModel=IdealSolution(T=t, xs=x), equilibrium_basis="Psat",
                                    caloric_basis="Psat", henry_components=henry,
                                    henry_abcdef=abcdef,
                                    T=t, P=p, zs=x)
            # A solute adds no volume in reflux, and has none in thermo here.
            cases.append({"cas": cas, "T": t, "P": p, "x": x, "thermo": answer(liq, x, t, volume=False)})
    return cases


FLASH_SYSTEMS = [
    ("NRTL", ["67-56-1", "7732-18-5"]),
    ("NRTL", ["64-17-5", "7732-18-5"]),
    ("NRTL", ["67-64-1", "67-56-1"]),
    ("NRTL", ["67-64-1", "67-56-1", "7732-18-5"]),
    ("UNIFAC", ["67-56-1", "7732-18-5"]),
]


def parts(name, cas):
    """thermo's FlashVL on the gamma-phi liquid, basis 'Poynting' (D4's
    default), and the ideal gas; and the liquid phase alone."""
    recs = [ga.record(c) for c in cas]
    consts = ChemicalConstantsPackage(
        Tcs=[r["constants"]["tc"] for r in recs], Pcs=[r["constants"]["pc"] for r in recs],
        omegas=[r["constants"]["omega"] for r in recs], MWs=[r["constants"]["mw"] * 1e3 for r in recs],
        CASs=cas)
    vps = [vapour_pressure(r) for r in recs]
    vls = [liquid_volume(r) for r in recs]
    cps = [heat_capacity(r) for r in recs]
    corr = PropertyCorrelationsPackage(constants=consts, VaporPressures=vps, VolumeLiquids=vls,
                                       HeatCapacityGases=cps, skip_missing=True)
    x0 = [1.0 / len(cas)] * len(cas)
    liq = GibbsExcessLiquid(VaporPressures=vps, VolumeLiquids=vls, HeatCapacityGases=cps,
                            GibbsExcessModel=activity(name, cas, 330.0, x0), equilibrium_basis="Poynting",
                            caloric_basis="Poynting", T=330.0, P=1e5, zs=x0)
    gas = IdealGas(HeatCapacityGases=cps, T=330.0, P=1e5, zs=x0)
    return FlashVL(consts, corr, liquid=liq, gas=gas), liq


def phis(liq, t, p, x):
    return [math.exp(v) for v in liq.to(T=t, P=p, zs=x).lnphis()]


def boundary(liq, given, value, z, vf):
    """A bubble (vf 0) or dew (vf 1) point by thermo's liquid fugacities and
    the ideal gas, y_i = x_i phi_i^L: thermo's own flash at a vapour
    fraction fails on these mixtures. At a bubble point x = z and
    sum z_i phi_i(x) = 1; at a dew point y = z, x_i = z_i/phi_i(x)
    normalised, found by substitution, and sum z_i/phi_i(x) = 1. Brent's
    method in T, or in ln P, to 1e-13."""
    nc = len(z)

    def x_of(t, p):
        if vf == 0.0:
            return z
        x = list(z)
        for _ in range(500):
            ph = phis(liq, t, p, x)
            new = [z[i] / ph[i] for i in range(nc)]
            tot = sum(new)
            new = [v / tot for v in new]
            if max(abs(new[i] - x[i]) for i in range(nc)) < 1e-15:
                return new
            x = new
        return x

    def gap(u):
        t, p = (value, math.exp(u)) if given == "T" else (u, value)
        x = x_of(t, p)
        ph = phis(liq, t, p, x)
        return sum(x[i] * ph[i] for i in range(nc)) - 1.0 if vf == 0.0 else             sum(z[i] / ph[i] for i in range(nc)) - 1.0

    if given == "T":
        lo, hi = math.log(1e2), math.log(1e8)
        grid = [lo + (hi - lo) * k / 400 for k in range(401)]
    else:
        grid = [250.0 + k * 0.5 for k in range(601)]
    vals = [gap(u) for u in grid]
    for u0, u1, g0, g1 in zip(grid, grid[1:], vals, vals[1:]):
        if (g0 > 0) != (g1 > 0):
            u = brentq(gap, u0, u1, xtol=1e-13, rtol=1e-15)
            t, p = (value, math.exp(u)) if given == "T" else (u, value)
            x = x_of(t, p)
            ph = phis(liq, t, p, x)
            y = [x[i] * ph[i] for i in range(nc)]
            if vf == 0.0:
                return {"T": t, "P": p, "VF": 0.0, "x": list(z), "y": y}
            return {"T": t, "P": p, "VF": 1.0, "x": x, "y": list(z)}
    raise ValueError("no boundary")


def split(liq, t, p, z):
    """A PT split from thermo's liquid fugacities and the ideal gas, by
    successive substitution on K_i = phi_i^L(x) with Rachford-Rice, to 1e-15
    in x: what FlashVL's PT flash should give, which it doesn't always near
    an azeotrope (it called a state between its own liquid's bubble and dew
    points all vapour)."""
    nc = len(z)
    x = list(z)
    for _ in range(2000):
        k = phis(liq, t, p, x)

        def rr(b):
            return sum(z[i] * (k[i] - 1.0) / (1.0 + b * (k[i] - 1.0)) for i in range(nc))
        beta = brentq(rr, 1e-14, 1.0 - 1e-14, xtol=1e-16, rtol=1e-15)
        new = [z[i] / (1.0 + beta * (k[i] - 1.0)) for i in range(nc)]
        if max(abs(new[i] - x[i]) for i in range(nc)) < 1e-15:
            x = new
            break
        x = new
    k = phis(liq, t, p, x)
    return {"T": t, "P": p, "VF": beta, "x": x, "y": [x[i] * k[i] for i in range(nc)]}


def flash_answer(res):
    out = {"T": res.T, "P": res.P, "VF": res.VF}
    if res.gas is not None:
        out["y"] = list(res.gas.zs)
    if res.liquids:
        out["x"] = list(res.liquids[0].zs)
    return out


def flash_cases(rng):
    cases = []
    for name, cas in FLASH_SYSTEMS:
        f, liq = parts(name, cas)
        for _ in range(6):
            raw = [rng.random() + 0.05 for _ in cas]
            z = [v / sum(raw) for v in raw]
            z[-1] = 1.0 - sum(z[:-1])
            p = round(10 ** rng.uniform(4.3, 5.7), 0)
            bub = boundary(liq, "P", p, z, 0.0)
            dew = boundary(liq, "P", p, z, 1.0)
            # a T in the two-phase range, one below it and one above
            for t in [round(0.5 * (bub["T"] + dew["T"]), 3), round(bub["T"] - 15.0, 3), round(dew["T"] + 15.0, 3)]:
                vl = flash_answer(f.flash(T=t, P=p, zs=z))
                inside = bub["T"] < t < dew["T"]
                case = {"model": name, "cas": cas, "kind": "pt", "z": z, "spec": [t, p],
                        "thermo": split(liq, t, p, z) if inside else vl}
                if inside and not 0.0 < vl["VF"] < 1.0:
                    case["flashvl"] = vl
                cases.append(case)
            cases.append({"model": name, "cas": cas, "kind": "bubble_t", "z": z, "spec": [p], "thermo": bub})
            cases.append({"model": name, "cas": cas, "kind": "dew_t", "z": z, "spec": [p], "thermo": dew})
            t = round(bub["T"], 1)
            cases.append({"model": name, "cas": cas, "kind": "bubble_p", "z": z, "spec": [t],
                          "thermo": boundary(liq, "T", t, z, 0.0)})
            cases.append({"model": name, "cas": cas, "kind": "dew_p", "z": z, "spec": [t],
                          "thermo": boundary(liq, "T", t, z, 1.0)})
    return cases


def azeotropes():
    """Each binary's azeotropes at 1 atm, where its bubble point's y1 - x1
    changes sign, by Brent's method to 1e-14 in x."""
    out = []
    for name, cas in FLASH_SYSTEMS:
        if len(cas) != 2:
            continue
        _, liq = parts(name, cas)

        def gap(x1):
            r = boundary(liq, "P", 101325.0, [x1, 1.0 - x1], 0.0)
            return r["y"][0] - x1

        xs = [k / 50.0 for k in range(1, 50)]
        gaps = [gap(x) for x in xs]
        found = []
        for a, b, g0, g1 in zip(xs, xs[1:], gaps, gaps[1:]):
            if (g0 > 0) != (g1 > 0):
                x1 = brentq(gap, a, b, xtol=1e-14, rtol=1e-15)
                found.append({"x1": x1, "T": boundary(liq, "P", 101325.0, [x1, 1.0 - x1], 0.0)["T"]})
        out.append({"model": name, "cas": cas, "P": 101325.0, "azeotropes": found})
    return out


HENRY_MIXED = ["7732-18-5", "107-21-1", "124-38-9", "1333-74-0"]


def henry_mixed(rng):
    """Two gases in two solvents: water's constants from the bank, the
    glycol's synthetic, mixed over the solvents present by thermo."""
    lines = henry_lines()
    cas = HENRY_MIXED
    recs = [ga.record(c) for c in cas]
    nc = len(cas)
    henry = [False, False, True, True]
    abcdef = [[(0.0,) * 6 for _ in range(nc)] for _ in range(nc)]
    constants = []
    for i in (2, 3):
        a, b, c, d, tmin, tmax = lines[cas[i]]
        abcdef[i][0] = (a, b, c, d, 0.0, 0.0)
        constants.append({"solute": i, "solvent": 0, "abcde": [a, b, c, d, 1.0], "tmin": tmin, "tmax": tmax})
        a2, d2 = a - math.log(8.0 if i == 2 else 3.0) - 0.33, d + 1e-3
        abcdef[i][1] = (a2, b, c, d2, 0.0, 0.0)
        constants.append({"solute": i, "solvent": 1, "abcde": [a2, b, c, d2, 1.0], "tmin": tmin, "tmax": tmax})
    states = []
    for _ in range(8):
        t = round(rng.uniform(290.0, 360.0), 2)
        p = round(10 ** rng.uniform(5.0, 6.5), 0)
        w = rng.uniform(0.1, 0.9)
        gas = [rng.random() * 1e-3 for _ in range(2)]
        rest = 1.0 - sum(gas)
        x = [w * rest, 0.0, gas[0], gas[1]]
        x[1] = 1.0 - x[0] - x[2] - x[3]
        liq = GibbsExcessLiquid(VaporPressures=[vapour_pressure(r) for r in recs[:2]] + [placeholder_pressure() for _ in recs[2:]],
                                VolumeLiquids=[liquid_volume(r) for r in recs[:2]] + [None, None],
                                HeatCapacityGases=[heat_capacity(r) for r in recs],
                                GibbsExcessModel=IdealSolution(T=t, xs=x), equilibrium_basis="Psat",
                                caloric_basis="Psat", henry_components=henry, henry_abcdef=abcdef,
                                henry_mode="solvents_with_parameters", T=t, P=p, zs=x)
        ans = answer(liq, x, t, volume=False)
        ans["dn"] = differences(liq, t, p, x)
        states.append({"T": t, "P": p, "x": x, "thermo": ans})
    return {"cas": cas, "constants": constants, "states": states}


def differences(liq, t, p, x):
    """n*d(ln phi_i)/dn_j at one mole in all, by Richardson's extrapolation of
    central differences of thermo's ln phi in n_j: a gas's step is half its
    amount at most, which ln phi doesn't depend on (each depends on the
    solvents' proportions alone), so the step costs nothing."""
    n = len(x)

    def lnphis(moles):
        total = sum(moles)
        return liq.to(T=t, P=p, zs=[m / total for m in moles]).lnphis()

    out = [[0.0] * n for _ in range(n)]
    for j in range(n):
        h = min(1e-3, 0.5 * x[j])
        d = []
        for step in (h, h / 2.0):
            up = list(x)
            dn = list(x)
            up[j] += step
            dn[j] -= step
            fu, fd = lnphis(up), lnphis(dn)
            d.append([(fu[i] - fd[i]) / (2.0 * step) for i in range(n)])
        for i in range(n):
            out[i][j] = (4.0 * d[1][i] - d[0][i]) / 3.0
    return [out[i][j] for i in range(n) for j in range(n)]


def main():
    rng = random.Random(SEED)
    out = {
        "source": {"generator": "tools/gen_m6_gamma_phi.py", "thermo": "0.6.1", "seed": SEED,
                   "generated": datetime.date.today().isoformat(), "R_thermo": ga.R_THERMO},
        "note": ("thermo's GibbsExcessLiquid on ChemSep's coefficients: ln phi, its derivatives (dn is "
                 "n*d(ln phi_i)/dn_j), and H, S, Cp less the ideal gas's (J/mol, J/(mol K)), V in m^3/mol. "
                 "With Poynting's factor, cp has thermo's numerical d2V/dT2 replaced by mpmath's, and "
                 "cp_thermo is thermo's own."),
        "liquid": liquid_cases(rng),
        "henry": henry_cases(rng),
        "flash": flash_cases(rng),
        "azeotropes": azeotropes(),
        "henry_mixed": henry_mixed(random.Random(SEED + 1)),
    }
    (ROOT / "testdata" / "m6_gamma_phi.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8",
                                                          newline="\n")
    print("wrote testdata/m6_gamma_phi.json:", len(out["liquid"]), "liquid states,", len(out["henry"]), "Henry's,",
          len(out["flash"]), "flashes,", sum(len(a["azeotropes"]) for a in out["azeotropes"]), "azeotropes,",
          len(out["henry_mixed"]["states"]), "in a mixed solvent")


if __name__ == "__main__":
    main()
