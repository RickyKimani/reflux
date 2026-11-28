# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""Write M5b's DWSIM references (step 10, criterion 9): testdata/m5b_dwsim.json.

DWSIM 9.0.5's transport of mixtures by state, on its "Peng-Robinson (PR)"
and "NRTL" packages at their defaults: a material stream at each case's T,
P and composition, calculated alone, and each phase's viscosity, thermal
conductivity and surface tension as it reports them (Phases 2, 3 and 4:
vapour, liquid1 and liquid2), with the phase's fraction, composition,
density and molar mass.

Each number is worked out again from DWSIM's own formulas (its source,
DanWBR/dwsim's PropertyPackage.vb, AUX_LIQVISCm, AUX_VAPVISCm, AUX_CONDTL,
AUX_CONDTG, AUX_SURFTM, as step 0's spike/transport/dwsim.py found them):
- a liquid's viscosity: the mole-fraction mean of each compound's
  correlation times Lucas's (1981) factor at |P - Psat|, a compound above
  its Tc taking its Henry's constant (AUX_KHenry, whatever the solvent) in
  place of Psat;
- a vapour's: the mole-fraction mean, then Jossi, Stiel and Thodos'
  residual on the mixture's mole-fraction-mean Tc, Pc, Vc and M, at the
  vapour's density;
- a liquid's conductivity: Li's rule on each compound's critical volume;
- a vapour's: the mole-fraction mean;
- the surface tension: the mole-fraction mean over the liquid's compounds
  below their Tc, each term divided by one less the moles of the compounds
  above their Tc met so far in the compounds' order (so the answer depends
  on that order: the "hexane with methane" case lists the "methane in
  hexane" case's compounds the other way round);
- a second liquid's viscosity and conductivity: the first's (the functions
  called without the phase).
Twice: on DWSIM's own compound values ("reproduced", the relative
difference from what it reports), and on reflux's component files' (the
"files" values, ChemSep LITE's correlations, as tools/gen_m5b_published.py
reads them), so that a difference from reflux's is put down to the data or
to the method. Each compound's values by both are recorded too.

And its water: the "Steam Tables (IAPWS-IF97)" package's viscosity and
conductivity, which are IAPS's 1985 formulations (its source), and its
"CoolProp" package's, which are CoolProp's IAPWS (not ChemSep's
correlations, as step 0 read the source: those are 0.7% to 20% off them):
recorded, at each state's T and density, for transport/dwsim_test.v to
hold against M5a's IAPWS functions.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5b_dwsim.py

It needs DWSIM 9.0.5 at %LOCALAPPDATA%\\DWSIM (Windows), and takes about a
minute. It reproduces its output byte for byte on the same day but the
water-1-butanol case's phases: DWSIM's liquid-liquid split moves in its
last digits from run to run (M6 step 12 found it), about 1e-13 of each
number, which no test's bound sees.
"""
import datetime
import json
import math
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m5b_published as files  # noqa: E402  (the component files' correlations)

import pythonnet  # noqa: E402

pythonnet.load("netfx")
import clr  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5b_dwsim.json"
HOME = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
sys.path.append(HOME)
os.chdir(HOME)
for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings", "DWSIM.SharedClasses",
             "DWSIM.Thermodynamics", "DWSIM.UnitOperations", "DWSIM.Inspector", "System.Buffers"]:
    clr.AddReference(name)

from System import Array, Double  # noqa: E402
from DWSIM.Automation import Automation3  # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402

vec = lambda xs: Array[Double](list(xs))  # noqa: E731
auto = Automation3()
BUTANOL = next(k for k in auto.AvailableCompounds.Keys if k.lower() in ("1-butanol", "n-butanol"))

CASES = [
    # (label, package, compounds by DWSIM's names, z, T K, P Pa)
    ("benzene-toluene liquid, 1 atm", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 300.0,
     101325.0),
    ("benzene-toluene liquid, 100 bar", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 300.0, 1e7),
    ("benzene-toluene, two phases", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 370.0, 101325.0),
    ("methanol-water, two phases", "NRTL", ["Methanol", "Water"], [0.5, 0.5], 350.0, 101325.0),
    ("methanol-water liquid", "NRTL", ["Methanol", "Water"], [0.3, 0.7], 300.0, 101325.0),
    ("reforming gas, 30 bar", "Peng-Robinson (PR)",
     ["Methane", "Hydrogen", "Water", "Carbon dioxide", "Carbon monoxide"], [0.15, 0.5, 0.2, 0.05, 0.1], 1000.0,
     3e6),
    ("methane-propane, 100 bar", "Peng-Robinson (PR)", ["Methane", "Propane"], [0.8, 0.2], 320.0, 1e7),
    ("methane in hexane, 20 bar", "Peng-Robinson (PR)", ["Methane", "N-hexane"], [0.05, 0.95], 300.0, 2e6),
    ("hexane with methane, 20 bar", "Peng-Robinson (PR)", ["N-hexane", "Methane"], [0.95, 0.05], 300.0, 2e6),
    ("water-1-butanol, two liquids", "NRTL", ["Water", BUTANOL], [0.7, 0.3], 300.0, 101325.0),
]
PHASES = {2: "vapour", 3: "liquid1", 4: "liquid2"}
WATER_STATES = [(300.0, 101325.0), (373.15, 1e5), (500.0, 1e6), (500.0, 1e7), (700.0, 1e7)]


def sheet(package, compounds):
    sim = auto.CreateFlowsheet()
    for c in compounds:
        sim.SelectedCompounds.Add(c, sim.AvailableCompounds[c])
    pp = sim.CreateAndAddPropertyPackage(package).__implementation__
    ms = sim.AddObject(ObjectType.MaterialStream, 50, 50, "s").GetAsObject()
    ms.SetPropertyPackage(pp)
    return sim, pp, ms


def num(v):
    return None if v is None else float(v)


def lucas_abs(tr, p, pc, pvap, w):
    """Lucas's compressed-liquid factor at |P - Psat|, DWSIM's."""
    dpr = abs((p - pvap) / pc)
    tr = min(tr, 1.0)
    a = 0.9991 - (0.0004674 / (1.0523 * tr ** -0.03877 - 1.0513))
    d = (0.3257 / (1.0039 - tr ** 2.573) ** 0.2906) - 0.2086
    cc = (-0.07921 + 2.1616 * tr - 13.404 * tr ** 2 + 44.1706 * tr ** 3 - 84.8291 * tr ** 4 + 96.1209 * tr ** 5
          - 59.8127 * tr ** 6 + 15.6719 * tr ** 7)
    return (1 + d * (dpr / 2.118) ** a) / (1 + cc * w * dpr)


def jst(eta0, t, v, tc, pc, vc, mm):
    """Jossi, Stiel and Thodos at the vapour's volume v (m3/mol), DWSIM's
    arithmetic: Pc in atm, M in g/mol, viscosity in micropoise inside."""
    pc = pc / 101325
    rr = vc / v
    xit = (tc / mm ** 3 / pc ** 4) ** (1 / 6)
    eta0 = eta0 * 1e6 / 100 * 1000
    tmp = ((1.023 + 0.23364 * rr + 0.58533 * rr ** 2 - 0.40758 * rr ** 3 + 0.093324 * rr ** 4) ** 4 - 1) / xit + eta0
    return tmp / 10000 / 1000


def dwsim_method(pure, t, p, x, v_vap):
    """DWSIM's formulas over a phase of mole fractions x, on each compound's
    values in `pure` (tc, pc K Pa, omega, vc m3/mol, mw g/mol, and at T: mu_l,
    k_l, sigma (None at or above Tc), psat, mu_v, k_v)."""
    out = {}
    n = len(x)
    out["mu_l"] = sum(xi * c["mu_l"] * lucas_abs(t / c["tc"], p, c["pc"], c["henry"] if t >= c["tc"] else c["psat"],
                                                 c["omega"]) for xi, c in zip(x, pure) if xi > 0)
    s = sum(xi * c["vc"] for xi, c in zip(x, pure) if xi > 0)
    phi = [xi * c["vc"] / s for xi, c in zip(x, pure)]
    out["k_l"] = sum(phi[i] * phi[j] * 2 / (1 / pure[i]["k_l"] + 1 / pure[j]["k_l"]) for i in range(n)
                     for j in range(n) if x[i] > 0 and x[j] > 0)
    sigma, dropped = 0.0, 0.0
    for xi, c in zip(x, pure):
        if c["sigma"] is None:
            dropped += xi
        else:
            sigma += xi / (1 - dropped) * c["sigma"]
    out["sigma"] = sigma
    eta0 = sum(xi * c["mu_v"] for xi, c in zip(x, pure))
    if v_vap is not None:
        tcm = sum(xi * c["tc"] for xi, c in zip(x, pure))
        pcm = sum(xi * c["pc"] for xi, c in zip(x, pure))
        vcm = sum(xi * c["vc"] for xi, c in zip(x, pure))
        mmm = sum(xi * c["mw"] for xi, c in zip(x, pure))
        out["mu_v"] = jst(eta0, t, v_vap, tcm, pcm, vcm, mmm)
    out["k_v"] = sum(xi * c["k_v"] for xi, c in zip(x, pure))
    return out


def dwsim_pure(pp, sim, compounds, t, p):
    """Each compound's values by DWSIM's own constants and functions."""
    out = []
    for k, name in enumerate(compounds):
        c = sim.SelectedCompounds[name]
        tc = float(c.Critical_Temperature)
        mu_l = pp.CalcCSTDepProp(str(c.LiquidViscosityEquation), c.Liquid_Viscosity_Const_A,
                                 c.Liquid_Viscosity_Const_B, c.Liquid_Viscosity_Const_C, c.Liquid_Viscosity_Const_D,
                                 c.Liquid_Viscosity_Const_E, t, tc)
        k_l = pp.CalcCSTDepProp(str(c.LiquidThermalConductivityEquation), c.Liquid_Thermal_Conductivity_Const_A,
                                c.Liquid_Thermal_Conductivity_Const_B, c.Liquid_Thermal_Conductivity_Const_C,
                                c.Liquid_Thermal_Conductivity_Const_D, c.Liquid_Thermal_Conductivity_Const_E, t, tc)
        out.append({
            "name": name, "cas": str(c.CAS_Number), "tc": tc, "pc": float(c.Critical_Pressure),
            "omega": float(c.Acentric_Factor), "vc": float(c.Critical_Volume) / 1000.0, "mw": float(c.Molar_Weight),
            "mu_l": float(mu_l), "k_l": float(k_l), "psat": float(pp.AUX_PVAPi(k, t)),
            "sigma": float(pp.AUX_SURFTi(c, t)) if t < tc else None,
            "henry": float(pp.AUX_KHenry(name, t)) if t >= tc else None,
            "mu_v": float(pp.AUX_VAPVISCi(c, t)), "k_v": float(pp.AUX_VAPTHERMCONDi(c, t, p)),
        })
    return out


def files_pure(cas, t, dw):
    """Each compound's values by reflux's component files, at T: ChemSep
    LITE's correlations, read as tools/gen_m5b_published.py reads them. A
    vapour's conductivity is its correlation's, at no pressure. A compound
    above its Tc takes DWSIM's Henry's constant, which is its method's, not
    data reflux has."""
    out = []
    for c, d in zip(cas, dw):
        f = files.Pure(c)
        out.append({
            "cas": c, "tc": f.tc, "pc": f.pc, "omega": f.w, "vc": f.vc, "mw": f.m,
            "mu_l": f.get("liquid_viscosity", t), "k_l": f.get("liquid_thermal_conductivity", t),
            "psat": f.get("vapor_pressure", t),
            "sigma": f.get("surface_tension", t) if t < f.tc else None,
            "henry": d["henry"] if t >= f.tc else None,
            "mu_v": f.get("vapor_viscosity", t), "k_v": f.get("vapor_thermal_conductivity", t),
        })
    return out


def rel(a, b):
    return None if a is None or b is None or b == 0 else (a - b) / b


def case(label, package, compounds, z, t, p):
    sim, pp, ms = sheet(package, compounds)
    ms.SetTemperature(t)
    ms.SetPressure(p)
    ms.SetMolarFlow(1.0)
    ms.SetOverallComposition(vec(z))
    pp.CurrentMaterialStream = ms
    ms.Calculate(True, True)
    pp.CurrentMaterialStream = ms  # its AUX_ functions read the stream's compounds
    dw = dwsim_pure(pp, sim, compounds, t, p)
    cas = [c["cas"] for c in dw]
    fl = files_pure(cas, t, dw)
    phases = []
    for i, name in PHASES.items():
        ph = ms.Phases[i]
        q = ph.Properties
        frac = num(q.molarfraction) or 0.0
        if frac <= 0:
            continue
        x = [float(ph.Compounds[c].MoleFraction or 0.0) for c in compounds]
        rho, mw = num(q.density), num(q.molecularWeight)
        v = mw / rho / 1000.0 if name == "vapour" else None
        got = {"mu": num(q.viscosity), "k": num(q.thermalConductivity), "sigma": num(q.surfaceTension)}
        keys = {"mu": "mu_v" if name == "vapour" else "mu_l", "k": "k_v" if name == "vapour" else "k_l",
                "sigma": "sigma"}
        mine = dwsim_method(dw, t, p, x, v)
        on_files = dwsim_method(fl, t, p, x, v)
        phases.append({
            "phase": name, "fraction": frac, "x": x, "density": rho, "molar_mass": mw,
            "dwsim": got,
            "reproduced": {k: rel(got[k], mine[keys[k]]) for k in got if got[k] is not None and keys[k] in mine},
            "on_files": {k: on_files[keys[k]] for k in got if keys[k] in on_files},
        })
    # A second liquid's viscosity and conductivity are the first's: the
    # functions called without the phase. Its "reproduced" is against the
    # first liquid's composition.
    by = {ph["phase"]: ph for ph in phases}
    if "liquid2" in by and "liquid1" in by:
        first = dwsim_method(dw, t, p, by["liquid1"]["x"], None)
        l2 = by["liquid2"]
        l2["reproduced"]["mu"] = rel(l2["dwsim"]["mu"], first["mu_l"])
        l2["reproduced"]["k"] = rel(l2["dwsim"]["k"], first["k_l"])
    return {"label": label, "package": package, "compounds": compounds, "cas": cas, "z": z, "T": t, "P": p,
            "pure_dwsim": dw, "pure_files": fl, "phases": phases}


def water():
    out = {}
    for key, package in [("steam_tables", "Steam Tables (IAPWS-IF97)"), ("coolprop", "CoolProp")]:
        rows = []
        for t, p in WATER_STATES:
            sim, pp, ms = sheet(package, ["Water"])
            ms.SetTemperature(t)
            ms.SetPressure(p)
            ms.SetMolarFlow(1.0)
            ms.SetOverallComposition(vec([1.0]))
            pp.CurrentMaterialStream = ms
            ms.Calculate(True, True)
            for i, name in PHASES.items():
                q = ms.Phases[i].Properties
                if (num(q.molarfraction) or 0.0) > 0:
                    rows.append({"T": t, "P": p, "phase": name, "density": num(q.density),
                                 "mu": num(q.viscosity), "k": num(q.thermalConductivity)})
        out[key] = rows
    return out


def main():
    cases = [case(*c) for c in CASES]
    for c in cases:
        for ph in c["phases"]:
            worst = max((abs(v) for v in ph["reproduced"].values() if v is not None), default=0.0)
            print(f"{c['label']:34} {ph['phase']:8} reproduced to {worst:.1e}")
    out = {
        "dwsim": str(auto.GetVersion()),
        "date": datetime.date.today().isoformat(),
        "generator": "tools/gen_m5b_dwsim.py",
        "units": "SI: Pa s, W/(m K), N/m, K, Pa, kg/m3; molar masses in g/mol, critical volumes in m3/mol",
        "cases": cases,
        "water": water(),
    }
    OUT.write_bytes((json.dumps(out, indent=1) + "\n").encode("utf-8"))
    print(f"wrote {OUT}: {len(cases)} cases")


if __name__ == "__main__":
    main()
