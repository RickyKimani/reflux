# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5"]
# ///
"""M5b step 0: DWSIM 9.0.5's transport properties of mixtures, by script.

For each case, a material stream on one of DWSIM's packages at its
defaults, calculated alone: each phase's viscosity, thermal conductivity
and surface tension as DWSIM reports them, and the same numbers worked out
again here from DWSIM's own compound constants by the formulas its source
uses (DanWBR/dwsim, DWSIM.Thermodynamics/PropertyPackages/PropertyPackage.vb,
AUX_LIQVISCm, AUX_VAPVISCm, AUX_CONDTL, AUX_CONDTG, AUX_SURFTM, and
Models/FluidProperties.vb's helpers):

- the liquid's viscosity: the mole-fraction mean of each compound's
  correlation, each times Lucas's (1981) pressure factor at |P - Psat|
  (its defaults: LiquidViscosity_MixingRule MoleAverage,
  LiquidViscosity_CorrectExpDataForPressure True);
- the vapour's: the mole-fraction mean of each compound's correlation, then
  Jossi, Stiel and Thodos' dense-gas correction on the mixture's
  mole-fraction-mean Tc, Pc, Vc and M (AUX_TCM, AUX_PCM, AUX_VCM, AUX_MMM);
- the liquid's conductivity: Li's rule with each compound's critical volume
  where Li has its liquid molar volume;
- the vapour's: the mole-fraction mean of each compound's correlation;
- the surface tension: the mole-fraction mean of each compound's, over the
  liquid's composition.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/transport/dwsim.py

Not part of the simulator: it answers the plan's questions (can a script
read DWSIM's mixture transport, and are its methods what its source says).
"""
import math
import os
import sys
import time

import pythonnet

pythonnet.load("netfx")
import clr  # noqa: E402

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

CASES = [
    # (label, package, compounds, z, T K, P Pa)
    ("benzene-toluene liquid, 1 atm", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 300.0, 101325.0),
    ("benzene-toluene liquid, 100 bar", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 300.0, 1e7),
    ("benzene-toluene at its bubble region", "Peng-Robinson (PR)", ["Benzene", "Toluene"], [0.4, 0.6], 370.0,
     101325.0),
    ("methanol-water, two phases", "NRTL", ["Methanol", "Water"], [0.5, 0.5], 350.0, 101325.0),
    ("methanol-water liquid", "NRTL", ["Methanol", "Water"], [0.3, 0.7], 300.0, 101325.0),
    ("reforming gas, 30 bar", "Peng-Robinson (PR)",
     ["Methane", "Hydrogen", "Water", "Carbon dioxide", "Carbon monoxide"],
     [0.15, 0.5, 0.2, 0.05, 0.1], 1000.0, 3e6),
    ("methane-propane, 100 bar (dense gas)", "Peng-Robinson (PR)", ["Methane", "Propane"], [0.8, 0.2], 320.0, 1e7),
    ("water-1-butanol, two liquids", "NRTL", ["Water", "BUTANOL"], [0.7, 0.3], 300.0, 101325.0),
]


BUTANOL = next(k for k in auto.AvailableCompounds.Keys if k.lower() in ("1-butanol", "n-butanol"))


def sheet(package, compounds):
    compounds[:] = [BUTANOL if c == "BUTANOL" else c for c in compounds]
    sim = auto.CreateFlowsheet()
    for c in compounds:
        comp = sim.AvailableCompounds[c]
        sim.SelectedCompounds.Add(c, comp)
    pp = sim.CreateAndAddPropertyPackage(package).__implementation__
    ms = sim.AddObject(ObjectType.MaterialStream, 50, 50, "s").GetAsObject()
    ms.SetPropertyPackage(pp)
    return sim, pp, ms


def props(ph):
    p = ph.Properties
    g = lambda v: None if v is None else float(v)  # noqa: E731
    return {"frac": g(p.molarfraction), "mu": g(p.viscosity), "k": g(p.thermalConductivity),
            "sigma": g(p.surfaceTension), "rho": g(p.density), "mw": g(p.molecularWeight)}


def fractions(ph, compounds):
    return [float(ph.Compounds[c].MoleFraction or 0.0) for c in compounds]


def eval_eq(pp, eq, a, b, c, d, e, t, tc):
    return pp.CalcCSTDepProp(str(eq), a, b, c, d, e, t, tc)


def lucas(tr, p, pc, pvap, w):
    dpr = abs((p - pvap) / pc)
    tr = min(tr, 1.0)
    a = 0.9991 - (0.0004674 / (1.0523 * tr ** -0.03877 - 1.0513))
    d = (0.3257 / (1.0039 - tr ** 2.573) ** 0.2906) - 0.2086
    cc = (-0.07921 + 2.1616 * tr - 13.404 * tr ** 2 + 44.1706 * tr ** 3 - 84.8291 * tr ** 4 + 96.1209 * tr ** 5
          - 59.8127 * tr ** 6 + 15.6719 * tr ** 7)
    return (1 + d * (dpr / 2.118) ** a) / (1 + cc * w * dpr)


def jst(eta0, t, v, tc, pc, vc, mm):
    pc = pc / 101325
    rr = vc / v
    xit = (tc / mm ** 3 / pc ** 4) ** (1 / 6)
    eta0 = eta0 * 1e6 / 100 * 1000
    tmp = ((1.023 + 0.23364 * rr + 0.58533 * rr ** 2 - 0.40758 * rr ** 3 + 0.093324 * rr ** 4) ** 4 - 1) / xit + eta0
    return tmp / 10000 / 1000


def ours(pp, sim, compounds, t, p, x, y, rho_v, mw_v):
    """The same numbers from DWSIM's own compound constants, by its formulas."""
    cps = [sim.SelectedCompounds[c] for c in compounds]
    out = {}
    if x is not None:
        mu = 0.0
        for k, (xi, c) in enumerate(zip(x, cps)):
            if xi <= 0:
                continue
            m = eval_eq(pp, c.LiquidViscosityEquation, c.Liquid_Viscosity_Const_A, c.Liquid_Viscosity_Const_B,
                        c.Liquid_Viscosity_Const_C, c.Liquid_Viscosity_Const_D, c.Liquid_Viscosity_Const_E, t,
                        c.Critical_Temperature)
            pv = pp.AUX_PVAPi(k, t)
            mu += xi * m * lucas(t / c.Critical_Temperature, p, c.Critical_Pressure, pv, c.Acentric_Factor)
        out["mu_l"] = mu
        ks = [eval_eq(pp, c.LiquidThermalConductivityEquation, c.Liquid_Thermal_Conductivity_Const_A,
                      c.Liquid_Thermal_Conductivity_Const_B, c.Liquid_Thermal_Conductivity_Const_C,
                      c.Liquid_Thermal_Conductivity_Const_D, c.Liquid_Thermal_Conductivity_Const_E, t,
                      c.Critical_Temperature) for c in cps]
        vcs = [c.Critical_Volume for c in cps]
        s = sum(xi * v for xi, v in zip(x, vcs) if xi > 0)
        phi = [xi * v / s for xi, v in zip(x, vcs)]
        out["k_l"] = sum(phi[i] * phi[j] * 2 / (1 / ks[i] + 1 / ks[j]) for i in range(len(x)) for j in range(len(x))
                         if x[i] > 0 and x[j] > 0)
        out["sigma"] = sum(xi * pp.AUX_SURFTi(c, t) for xi, c in zip(x, cps) if t < c.Critical_Temperature)
    if y is not None:
        eta0 = sum(yi * pp.AUX_VAPVISCi(c, t) for yi, c in zip(y, cps))
        tcm = sum(yi * c.Critical_Temperature for yi, c in zip(y, cps))
        pcm = sum(yi * c.Critical_Pressure for yi, c in zip(y, cps))
        mmm = sum(yi * c.Molar_Weight for yi, c in zip(y, cps))
        vcm = sum(yi * c.Critical_Volume for yi, c in zip(y, cps)) / 1000  # AUX_VCM: a mole-fraction mean, m3/mol
        out["mu_v0"] = eta0
        out["mu_v"] = jst(eta0, t, mw_v / rho_v / 1000, tcm, pcm, vcm, mmm)
        out["k_v"] = sum(yi * pp.AUX_VAPTHERMCONDi(c, t, p) for yi, c in zip(y, cps))
    return out


def rel(a, b):
    if a is None or b is None or b == 0:
        return float("nan")
    return (a - b) / b


def main():
    print("DWSIM", auto.GetVersion())
    for label, package, compounds, z, t, p in CASES:
        sim, pp, ms = sheet(package, compounds)
        ms.SetTemperature(t)
        ms.SetPressure(p)
        ms.SetMolarFlow(1.0)
        ms.SetOverallComposition(vec(z))
        pp.CurrentMaterialStream = ms
        t0 = time.perf_counter()
        n = 20
        for _ in range(n):
            ms.Calculate(True, True)
        dt = (time.perf_counter() - t0) / n
        ph = {i: props(ms.Phases[i]) for i in (0, 1, 2, 3, 4)}
        print(f"\n== {label} ({package}), {t} K, {p:g} Pa: a stream's calculation {dt * 1e3:.2f} ms")
        names = {0: "mixture", 1: "liquid (overall)", 2: "vapour", 3: "liquid1", 4: "liquid2"}
        for i, d in ph.items():
            if (d["frac"] or 0) > 0 or i == 0:
                print(f"   {names[i]:17} frac {d['frac']}: mu {d['mu']}, k {d['k']}, sigma {d['sigma']}, "
                      f"rho {d['rho']}")
        x = fractions(ms.Phases[3], compounds) if (ph[3]["frac"] or 0) > 0 else None
        y = fractions(ms.Phases[2], compounds) if (ph[2]["frac"] or 0) > 0 else None
        pp.CurrentMaterialStream = ms  # its AUX_ functions read the stream's compounds
        o = ours(pp, sim, compounds, t, p, x, y, ph[2]["rho"], ph[2]["mw"])
        if x is not None:
            print(f"   liquid1 x {x}: from its source mu {o['mu_l']:.10g} ({rel(ph[3]['mu'], o['mu_l']):.2e}), "
                  f"k {o['k_l']:.10g} ({rel(ph[3]['k'], o['k_l']):.2e}), sigma {o['sigma']:.10g} "
                  f"(liquid1's {rel(ph[3]['sigma'], o['sigma']):.2e})")
        if y is not None:
            print(f"   vapour y {y}: from its source mu {o['mu_v']:.10g} ({rel(ph[2]['mu'], o['mu_v']):.2e}; "
                  f"before JST {o['mu_v0']:.6g}), k {o['k_v']:.10g} ({rel(ph[2]['k'], o['k_v']):.2e})")


if __name__ == "__main__":
    main()
