# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""M6 step 0's spike: DWSIM 9's activity-coefficient packages, by script.

Run from the repository root, on Windows with DWSIM 9 installed:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/dwsim_activity/run.py

It finds, for docs/spike-dwsim-activity.md:
- each package's settings (its Parameters), and NRTL's and UNIQUAC's
  parameters for methanol-water, against ChemSep's banks;
- which form and which R DWSIM's NRTL uses, from its gamma against thermo's
  NRTL with tau = A/RT for each candidate R;
- how its liquid enthalpy is taken: from the activity model and Psat's slope,
  or from a heat of vaporization;
- a PT flash of methanol-water by each package;
- whether its three-phase flash runs on water-1-butanol, and what it gives.
"""
import math
import os
import sys

import pythonnet

pythonnet.load("netfx")
import clr  # noqa: E402

HOME = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
sys.path.append(HOME)
HERE = os.getcwd()
os.chdir(HOME)  # DWSIM finds its databases relative to the working directory
for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings",
             "DWSIM.SharedClasses", "DWSIM.Thermodynamics", "DWSIM.UnitOperations",
             "DWSIM.Inspector", "System.Buffers"]:
    clr.AddReference(name)

from System import Array, Double  # noqa: E402
from DWSIM.Automation import Automation3  # noqa: E402
from DWSIM.Interfaces.Enums import FlashCalculationType  # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402
from DWSIM.Thermodynamics.PropertyPackages import State  # noqa: E402

P_ATM = 101325.0
vec = lambda xs: Array[Double](list(xs))  # noqa: E731
auto = Automation3()
print(auto.GetVersion())


def sheet(compounds, package, z):
    sim = auto.CreateFlowsheet()
    for c in compounds:
        comp = sim.AvailableCompounds[c]
        sim.SelectedCompounds.Add(c, comp)
    pp = sim.CreateAndAddPropertyPackage(package).__implementation__
    s = sim.AddObject(ObjectType.MaterialStream, 50, 50, "feed").GetAsObject()
    s.SetTemperature(330.0)
    s.SetPressure(P_ATM)
    s.SetMolarFlow(1.0)
    s.SetOverallComposition(vec(z))
    pp.CurrentMaterialStream = s
    return sim, pp


def members(obj, words):
    names = [n for n in dir(obj) if any(w.lower() in n.lower() for w in words)]
    return sorted(set(names))


def show_parameters(pp):
    try:
        keys = sorted(str(k) for k in pp.Parameters.Keys)
        for k in keys:
            print(f"    {k} = {pp.Parameters[k]}")
    except Exception as e:  # noqa: BLE001
        print("    no Parameters:", e)


def pair_parameters(pp, a, b):
    """The binary parameters a package holds for a pair, whatever it calls them."""
    for attr in ["m_uni", "m_pr", "m_nrtl", "m_uniquac", "m_wilson"]:
        holder = getattr(pp, attr, None)
        if holder is None:
            continue
        ip = getattr(holder, "InteractionParameters", None)
        if ip is None:
            print(f"    {attr}: {members(holder, ['inter', 'param'])}")
            continue
        for x, y in [(a, b), (b, a)]:
            if ip.ContainsKey(x) and ip[x].ContainsKey(y):
                d = ip[x][y]
                fields = {n: getattr(d, n) for n in dir(d) if n[0] in "ABCa" and not n.startswith("__")
                          and isinstance(getattr(d, n), (int, float))}
                return attr, x, y, fields
    return None


def psat_dwsim(pp, i, t):
    return pp.AUX_PVAPi(i, t)


def liquid_gamma(pp, x, t, p):
    """gamma from DWSIM's liquid phi, phi = gamma Psat / P for an ideal vapour with no Poynting."""
    phi = list(pp.DW_CalcFugCoeff(vec(x), t, p, State.Liquid))
    return [phi[i] * p / psat_dwsim(pp, i, t) for i in range(len(x))]


def nrtl_gamma(x, t, a12, a21, alpha, r_cal):
    tau12, tau21 = a12 / (r_cal * t), a21 / (r_cal * t)
    g12, g21 = math.exp(-alpha * tau12), math.exp(-alpha * tau21)
    x1, x2 = x
    l1 = x2 ** 2 * (tau21 * (g21 / (x1 + x2 * g21)) ** 2 + tau12 * g12 / (x2 + x1 * g12) ** 2)
    l2 = x1 ** 2 * (tau12 * (g12 / (x2 + x1 * g12)) ** 2 + tau21 * g21 / (x1 + x2 * g21) ** 2)
    return [math.exp(l1), math.exp(l2)]


def flash_pt(pp, t, p, z):
    r = pp.CalculateEquilibrium(FlashCalculationType.PressureTemperature, p, t, vec(z), None, 0.0)
    return {"beta": r.GetVaporPhaseMoleFraction(), "y": list(r.GetVaporPhaseMoleFractions()),
            "x1": list(r.GetLiquidPhase1MoleFractions()), "l2": r.GetLiquidPhase2MoleFraction(),
            "x2": list(r.GetLiquidPhase2MoleFractions()) if r.GetLiquidPhase2MoleFraction() > 0 else [],
            "iterations": r.IterationsTaken}


MW = ["Methanol", "Water"]
CHEMSEP_NRTL = {"A12": -189.0469, "A21": 792.8020, "alpha": 0.2999}  # nrtl.ipd, methanol/water, cal/mol
CHEMSEP_UNIQUAC = {"A12": -337.1298, "A21": 549.2958}

print("\n== methanol-water, each package")
x = [0.3, 0.7]
for package in ["NRTL", "UNIQUAC", "Wilson", "UNIFAC", "Modified UNIFAC (Dortmund)", "UNIFAC-LL"]:
    print(f"\n-- {package}")
    try:
        sim, pp = sheet(MW, package, [0.5, 0.5])
    except Exception as e:  # noqa: BLE001
        print("   can't make it:", e)
        continue
    print("   class", pp.GetType().FullName)
    show_parameters(pp)
    got = pair_parameters(pp, "Methanol", "Water")
    print("   pair parameters:", got)
    for c in MW:
        comp = sim.SelectedCompounds[c]
        print(f"   {c}: db {comp.OriginalDB}, Psat eq {comp.VaporPressureEquation}, "
              f"Hvap eq {getattr(comp, 'HVaporizationEquation', '?')}")
    t = 330.0
    gam = liquid_gamma(pp, x, t, P_ATM)
    print(f"   gamma from phi_L at {t} K, x = {x}: {gam}")
    if package == "NRTL":
        for label, r_cal in [("1.98721", 1.98721), ("R/4.184", 8.314462618 / 4.184),
                             ("R/4.1868", 8.314462618 / 4.1868)]:
            ref = nrtl_gamma(x, t, CHEMSEP_NRTL["A12"], CHEMSEP_NRTL["A21"], CHEMSEP_NRTL["alpha"], r_cal)
            print(f"   NRTL tau = A/({label} T): {ref}, ln-gamma gaps {[math.log(a / b) for a, b in zip(gam, ref)]}")
    # The liquid's enthalpy: pure water, H_L - H_ig, against -Hvap and Clapeyron.
    for i, c in enumerate(MW):
        pure = [1.0 if j == i else 0.0 for j in range(2)]
        comp = sim.SelectedCompounds[c]
        mw = comp.Molar_Weight
        try:
            dep = pp.DW_CalcEnthalpyDeparture(vec(pure), t, P_ATM, State.Liquid) * mw  # J/mol
        except Exception as e:  # noqa: BLE001
            dep = float("nan")
            print("   no enthalpy departure:", e)
        hvap = pp.RET_HVAPM(vec(pure), t) * mw if hasattr(pp, "RET_HVAPM") else float("nan")
        h = 1e-3 * t
        dlnp = (math.log(psat_dwsim(pp, i, t + h)) - math.log(psat_dwsim(pp, i, t - h))) / (2 * h)
        clapeyron = 8.314462618 * t * t * dlnp  # ideal vapour, liquid volume neglected
        print(f"   {c}: H_L - H_ig = {dep:.2f} J/mol; Hvap correlation {hvap:.2f}; Clapeyron (ideal vapour) {clapeyron:.2f}")
    try:
        print("   PT flash 350 K, 1 atm, z = 0.5/0.5:", flash_pt(pp, 350.0, P_ATM, [0.5, 0.5]))
    except Exception as e:  # noqa: BLE001
        print("   PT flash failed:", e)

print("\n== the settings and the gamma functions, NRTL and UNIQUAC")
from thermo.uniquac import UNIQUAC  # noqa: E402

R_CAL = 1.98721
for package in ["NRTL", "UNIQUAC"]:
    sim, pp = sheet(MW, package, [0.5, 0.5])
    print(f"-- {package}")
    for n in members(pp, ["enth", "poynting", "vapor", "ideal", "mode", "liqdens", "calc", "fugac"]):
        try:
            v = getattr(pp, n)
            if not callable(v):
                print(f"   {n} = {v}")
        except Exception:  # noqa: BLE001
            pass
    print("   m_uni members:", members(pp.m_uni, ["gam", "calc", "ln"]))
    t = 330.0
    # DWSIM's own liquid volumes, for its Poynting term.
    vols = []
    for i in range(2):
        pure = [1.0 if j == i else 0.0 for j in range(2)]
        try:
            vols.append(pp.AUX_LIQDENSi(sim.SelectedCompounds[MW[i]], t))
        except Exception as e:  # noqa: BLE001
            vols.append(None)
            print("   no AUX_LIQDENSi:", e)
    print("   liquid densities (kg/m3):", vols)
    gam = liquid_gamma(pp, x, t, P_ATM)
    if all(v is not None for v in vols):
        mws = [sim.SelectedCompounds[c].Molar_Weight / 1000 for c in MW]
        poy = [mws[i] / vols[i] * (P_ATM - psat_dwsim(pp, i, t)) / (8.314462618 * t) for i in range(2)]
        print("   Poynting VL(P - Psat)/RT:", poy)
        gam_np = [g / math.exp(pv) for g, pv in zip(gam, poy)]
        print("   gamma without Poynting:", gam_np)
        if package == "NRTL":
            for label, r_cal in [("1.98721", 1.98721), ("R/4.184", 8.314462618 / 4.184),
                                 ("R/4.1868", 8.314462618 / 4.1868)]:
                ref = nrtl_gamma(x, t, CHEMSEP_NRTL["A12"], CHEMSEP_NRTL["A21"], CHEMSEP_NRTL["alpha"], r_cal)
                print(f"   NRTL A/({label} T) gaps in ln gamma: {[math.log(a / b) for a, b in zip(gam_np, ref)]}")
        else:
            r = [sim.SelectedCompounds[c].UNIQUAC_R for c in MW]
            q = [sim.SelectedCompounds[c].UNIQUAC_Q for c in MW]
            print(f"   DWSIM's r {r}, q {q}")
            for label, tau in [
                ("exp(-A/RT)", lambda a: math.exp(-a / (R_CAL * t))),
                ("exp(-A/T)", lambda a: math.exp(-a / t)),
            ]:
                taus = [[1.0, tau(CHEMSEP_UNIQUAC["A12"])], [tau(CHEMSEP_UNIQUAC["A21"]), 1.0]]
                ref = UNIQUAC(xs=x, T=t, rs=r, qs=q, tau_coeffs=None, ABCDEF=None,
                              tau_as=[[math.log(v) for v in row] for row in taus]).gammas()
                print(f"   UNIQUAC tau = {label}: {ref}, gaps {[math.log(a / b) for a, b in zip(gam_np, ref)]}")

print("\n== DWSIM's butanols:", [k for k in auto.AvailableCompounds.Keys if "utanol" in k.lower()][:12])

print("\n== water-1-butanol, three phases")
BUTANOL = next((k for k in auto.AvailableCompounds.Keys if k.lower() in ("1-butanol", "n-butanol")), "1-Butanol")
for package in ["NRTL", "UNIQUAC"]:
    try:
        sim, pp = sheet(["Water", BUTANOL], package, [0.5, 0.5])
    except Exception as e:  # noqa: BLE001
        print(package, "can't make it:", e)
        continue
    print(f"-- {package}: pair {pair_parameters(pp, 'Water', BUTANOL)}")
    print("   flash members:", members(pp, ["flash"])[:30])
    algos = None
    try:
        import DWSIM.Thermodynamics.PropertyPackages.Auxiliary.FlashAlgorithms as FA  # noqa: N812

        algos = [n for n in dir(FA) if not n.startswith("_")]
        print("   flash algorithms:", algos)
    except Exception as e:  # noqa: BLE001
        print("   no flash algorithms module:", e)
    for label, make in [("default", None), ("NestedLoops3PV3", "NestedLoops3PV3"),
                        ("GibbsMinimization3P", "GibbsMinimization3P"), ("NestedLoopsSLE", None)]:
        if make is None and label != "default":
            continue
        if make is not None:
            if algos is None or make not in algos:
                print(f"   {label}: not available")
                continue
            try:
                pp.FlashBase = getattr(FA, make)()
            except Exception as e:  # noqa: BLE001
                print(f"   {label}: can't set it as FlashBase ({e}); members {members(pp, ['algorithm', 'flashbase'])}")
                continue
        for t in [298.15, 360.0, 366.0]:
            try:
                print(f"   {label}, {t} K, 1 atm: {flash_pt(pp, t, P_ATM, [0.5, 0.5])}")
            except Exception as e:  # noqa: BLE001
                print(f"   {label}, {t} K: failed: {e}")
sys.stdout.flush()
os.chdir(HERE)
os._exit(0)
