# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "thermo==0.6.1"]
# ///
"""M3 step 0 spike: drive DWSIM 9 from Python, with no GUI.

    uv run --python 3.13 spike/dwsim/run.py

It needs Windows and DWSIM 9 installed for the current user, in
%LOCALAPPDATA%\\DWSIM. DWSIM's classic build runs on .NET Framework, which
pythonnet reaches as "netfx". pythonnet 3.0.5 stops at Python 3.13.

The script:
  1. builds a flowsheet with a methane/n-butane stream on Peng-Robinson;
  2. reads the constants, Cp coefficients, kij and flash settings DWSIM uses;
  3. runs a PT flash through the flowsheet solver, and again through the
     property package directly;
  4. runs a PH flash at the enthalpy the PT flash gave, which should land back
     on the same temperature;
  5. tightens DWSIM's flash tolerances, and measures how well its answer
     satisfies its own equilibrium condition;
  6. times DWSIM's flashes;
  7. explains what is left: the PT flash in thermo's PR, and in a small PR
     here with DWSIM's rounded constants; DWSIM's enthalpy split into its
     ideal-gas and departure parts.

DWSIM reports enthalpy in kJ/kg. Multiplying by the molar mass in g/mol
gives J/mol.
"""

import math
import os
import sys
import time

T, P, Z = 250.0, 30e5, [0.4, 0.6]  # K, Pa, mole fractions
NAMES = ["Methane", "N-butane"]  # DWSIM's ChemSep entries ("n-Butane" is CoolProp's)
T_REF = 298.15  # K, DWSIM's ideal-gas enthalpy reference
SQRT2 = math.sqrt(2)

# Peng-Robinson's constants: exact, and as DWSIM's PengRobinson class has them.
EXACT = dict(oa=0.45723552892138218, ob=0.077796073903888456, r=8.314462618)
DWSIM_PR = dict(oa=0.45724, ob=0.0778, r=8.314)


def load_dwsim():
    import pythonnet

    pythonnet.load("netfx")
    import clr

    home = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
    sys.path.append(home)
    os.chdir(home)  # DWSIM finds its databases relative to the working directory
    for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings",
                 "DWSIM.SharedClasses", "DWSIM.Thermodynamics", "DWSIM.UnitOperations",
                 "DWSIM.Inspector", "System.Buffers"]:
        clr.AddReference(name)


def timed(f, *args):
    t0 = time.perf_counter()
    out = f(*args)
    return out, time.perf_counter() - t0


def main():
    _, t_load = timed(load_dwsim)
    from System import Array, Double
    from DWSIM.Automation import Automation3
    from DWSIM.Interfaces.Enums import FlashCalculationType, FlashSetting
    from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType
    from DWSIM.Thermodynamics.PropertyPackages import State

    auto, t_auto = timed(Automation3)
    sim, t_sheet = timed(auto.CreateFlowsheet)
    print(auto.GetVersion())
    print(f"start-up: {t_load:.2f} s loading assemblies, {t_auto:.2f} s Automation3, "
          f"{t_sheet:.2f} s CreateFlowsheet")

    for name in NAMES:
        sim.SelectedCompounds.Add(name, sim.AvailableCompounds[name])
    pp = sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)").__implementation__
    s = sim.AddObject(ObjectType.MaterialStream, 50, 50, "feed").GetAsObject()
    s.SetTemperature(T)
    s.SetPressure(P)
    s.SetMolarFlow(1.0)
    s.SetOverallComposition(Array[Double](Z))
    pp.CurrentMaterialStream = s  # the package reads its compounds from this stream
    vec = lambda xs: Array[Double](list(xs))

    # What DWSIM uses.
    print("\nconstants (DWSIM's database):")
    comps = []
    for name in NAMES:
        c = sim.SelectedCompounds[name]
        cp = [c.Ideal_Gas_Heat_Capacity_Const_A, c.Ideal_Gas_Heat_Capacity_Const_B,
              c.Ideal_Gas_Heat_Capacity_Const_C, c.Ideal_Gas_Heat_Capacity_Const_D,
              c.Ideal_Gas_Heat_Capacity_Const_E]
        comps.append(dict(name=name, cas=c.CAS_Number, tc=c.Critical_Temperature,
                          pc=c.Critical_Pressure, omega=c.Acentric_Factor, mw=c.Molar_Weight, cp=cp))
        print(f"  {name:9} {c.CAS_Number:9} {c.OriginalDB}: Tc {c.Critical_Temperature} K, "
              f"Pc {c.Critical_Pressure} Pa, omega {c.Acentric_Factor}, M {c.Molar_Weight} g/mol, "
              f"Cp eq {c.IdealgasCpEquation} {cp} J/(kmol K)")
    K = pp.RET_VKij()
    kij = [[K.GetValue(i, j) for j in range(2)] for i in range(2)]
    print(f"  kij {kij}")
    settings = pp.FlashSettings
    print(f"  flash: {pp.FlashBase.GetType().Name}, approach {pp.FlashCalculationApproach}, "
          f"PT tolerances {settings[FlashSetting.PTFlash_External_Loop_Tolerance]} / "
          f"{settings[FlashSetting.PTFlash_Internal_Loop_Tolerance]}, "
          f"PH tolerances {settings[FlashSetting.PHFlash_External_Loop_Tolerance]} / "
          f"{settings[FlashSetting.PHFlash_Internal_Loop_Tolerance]}")
    mw = sum(z * c["mw"] for z, c in zip(Z, comps))

    def flash_pt(t=T):
        r = pp.CalculateEquilibrium(FlashCalculationType.PressureTemperature, P, t, vec(Z), None, 0.0)
        out = dict(beta=r.GetVaporPhaseMoleFraction(), y=list(r.GetVaporPhaseMoleFractions()),
                   x=list(r.GetLiquidPhase1MoleFractions()), h=r.CalculatedEnthalpy * mw,
                   iters=r.IterationsTaken)
        fv = pp.DW_CalcFugCoeff(vec(out["y"]), t, P, State.Vapor)
        fl = pp.DW_CalcFugCoeff(vec(out["x"]), t, P, State.Liquid)
        out["residual"] = max(abs(math.log(out["y"][i] * fv[i]) - math.log(out["x"][i] * fl[i]))
                              for i in range(len(Z)))
        return out

    def flash_ph(h):
        r = pp.CalculateEquilibrium(FlashCalculationType.PressureEnthalpy, P, h / mw, vec(Z), None, 0.0)
        return dict(t=r.CalculatedTemperature, beta=r.GetVaporPhaseMoleFraction(), iters=r.IterationsTaken)

    # PT through the flowsheet solver.
    errs, t_fs = timed(auto.CalculateFlowsheet2, sim)
    assert not errs, [str(e) for e in errs]
    vap, liq = s.Phases[2], s.Phases[3]
    sheet = dict(beta=vap.Properties.molarfraction,
                 y=[vap.Compounds[n].MoleFraction for n in NAMES],
                 x=[liq.Compounds[n].MoleFraction for n in NAMES],
                 h=s.Phases[0].Properties.molar_enthalpy)  # kJ/kmol = J/mol
    print(f"\nPT, flowsheet solver ({t_fs:.3f} s, first call):")
    print(f"  beta {sheet['beta']!r}, y {sheet['y']}, x {sheet['x']}, H {sheet['h']!r} J/mol")

    # PT through the property package. The flowsheet solver leaves the
    # package without a current stream, so set it again.
    pp.CurrentMaterialStream = s
    pt, t_pt = timed(flash_pt)
    same = all(pt[k] == sheet[k] for k in ("beta", "y", "x"))
    print(f"\nPT, property package ({t_pt * 1e3:.1f} ms, {pt['iters']} iterations):")
    print(f"  beta {pt['beta']!r}, y {pt['y']}, x {pt['x']}, H {pt['h']!r} J/mol")
    print(f"  same as the flowsheet: {same}; H differs by {pt['h'] - sheet['h']:.3g} J/mol; "
          f"|ln fV - ln fL| {pt['residual']:.3g}")

    # PH back to the same state.
    ph, t_ph = timed(flash_ph, pt["h"])
    print(f"\nPH at that H ({t_ph * 1e3:.1f} ms, {ph['iters']} iterations): T {ph['t']!r} K "
          f"(off by {ph['t'] - T:.3g} K), beta off by {ph['beta'] - pt['beta']:.3g}")

    # Tighter tolerances. Past 1e-8 nothing changes: the residual stays put.
    print("\ntighter tolerances, PT:")
    for tol in ("1E-6", "1E-8", "1E-12"):
        for key in (FlashSetting.PTFlash_External_Loop_Tolerance, FlashSetting.PTFlash_Internal_Loop_Tolerance,
                    FlashSetting.PHFlash_External_Loop_Tolerance, FlashSetting.PHFlash_Internal_Loop_Tolerance):
            settings[key] = tol
        tight = flash_pt()
        print(f"  {tol:5}: beta {tight['beta']!r} ({tight['iters']} iterations, moved "
              f"{tight['beta'] - pt['beta']:.3g}), |ln fV - ln fL| {tight['residual']:.3g}")
    ph = flash_ph(tight["h"])
    print(f"  PH at 1E-12: T off by {ph['t'] - T:.3g} K ({ph['iters']} iterations)")

    # Timings, warm, at the tight tolerances.
    n = 200
    t0 = time.perf_counter()
    for i in range(n):
        pp.CalculateEquilibrium(FlashCalculationType.PressureTemperature, P, T + 0.05 * (i % 40), vec(Z), None, 0.0)
    t_pt_warm = (time.perf_counter() - t0) / n
    t0 = time.perf_counter()
    for i in range(n):
        pp.CalculateEquilibrium(FlashCalculationType.PressureEnthalpy, P, tight["h"] / mw * (1 + 1e-4 * (i % 40)),
                                vec(Z), None, 0.0)
    t_ph_warm = (time.perf_counter() - t0) / n
    m = 10
    betas = set()
    t0 = time.perf_counter()
    for i in range(m):
        s.SetTemperature(T + 0.1 * (i + 1))
        assert not auto.CalculateFlowsheet2(sim)
        betas.add(s.Phases[2].Properties.molarfraction)
    t_fs_warm = (time.perf_counter() - t0) / m
    assert len(betas) == m, "a flowsheet solve was skipped"
    print(f"\ntimings, warm: PT {t_pt_warm * 1e3:.2f} ms, PH {t_ph_warm * 1e3:.2f} ms "
          f"(property package, {n} each); flowsheet solve {t_fs_warm * 1e3:.1f} ms ({m}, each one recalculated)")

    pp.CurrentMaterialStream = s
    explain(comps, kij, tight, pp, State, vec)


def pr_phase(comps, kij, x, t, p, vapor, k):
    """Z, ln phi and the enthalpy departure (J/mol) of one PR phase."""
    import numpy as np

    tc = np.array([c["tc"] for c in comps])
    pc = np.array([c["pc"] for c in comps])
    w = np.array([c["omega"] for c in comps])
    x, kij, r = np.array(x), np.array(kij), k["r"]
    m = 0.37464 + 1.54226 * w - 0.26992 * w**2
    s = 1 + m * (1 - np.sqrt(t / tc))
    ac = k["oa"] * r**2 * tc**2 / pc
    ai, dai = ac * s**2, -m * ac * s / np.sqrt(t * tc)
    bi = k["ob"] * r * tc / pc
    root = np.sqrt(np.outer(ai, ai))
    aij = root * (1 - kij)
    a, b = x @ aij @ x, x @ bi
    da = x @ (0.5 * root * np.add.outer(dai / ai, dai / ai) * (1 - kij)) @ x
    A, B = a * p / (r * t) ** 2, b * p / (r * t)
    roots = sorted(z.real for z in np.roots([1, -(1 - B), A - 3 * B**2 - 2 * B, -(A * B - B**2 - B**3)])
                   if abs(z.imag) < 1e-12 and z.real > B)
    z = roots[-1] if vapor else roots[0]
    log_term = math.log((z + (1 + SQRT2) * B) / (z + (1 - SQRT2) * B))
    ln_phi = bi / b * (z - 1) - math.log(z - B) - A / (2 * SQRT2 * B) * (2 * (aij @ x) / a - bi / b) * log_term
    h_dep = r * t * (z - 1) + (t * da - a) / (2 * SQRT2 * b) * log_term
    return z, ln_phi, h_dep


def pr_flash(comps, kij, k):
    """A PT flash converged to round-off: successive substitution with an
    exact Rachford-Rice by bisection."""
    import numpy as np

    zs = np.array(Z)
    tc = np.array([c["tc"] for c in comps])
    pc = np.array([c["pc"] for c in comps])
    w = np.array([c["omega"] for c in comps])
    kv = pc / P * np.exp(5.373 * (1 + w) * (1 - tc / T))  # Wilson
    for _ in range(200):
        lo, hi = 0.0, 1.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if np.sum(zs * (kv - 1) / (1 + mid * (kv - 1))) > 0 else (lo, mid)
        beta = 0.5 * (lo + hi)
        x = zs / (1 + beta * (kv - 1))
        y = kv * x
        x, y = x / x.sum(), y / y.sum()
        kv = np.exp(pr_phase(comps, kij, x, T, P, False, k)[1] - pr_phase(comps, kij, y, T, P, True, k)[1])
    return beta, list(y), list(x)


def cp_integral(cp, t1, t2):
    """ChemSep eq 16, Cp = A + exp(B/T + C + D T + E T^2) in J/(kmol K),
    integrated from t1 to t2, in J/mol."""
    from scipy.integrate import quad

    a, b, c, d, e = cp
    rest, _ = quad(lambda t: math.exp(b / t + c + d * t + e * t * t), t1, t2, epsabs=0, epsrel=1e-13)
    return (a * (t2 - t1) + rest) / 1000


def explain(comps, kij, dw, pp, State, vec):
    from thermo import CEOSGas, CEOSLiquid, ChemicalConstantsPackage, FlashVL, PRMIX
    from thermo import HeatCapacityGas, PropertyCorrelationsPackage

    tcs = [c["tc"] for c in comps]
    pcs = [c["pc"] for c in comps]
    ws = [c["omega"] for c in comps]
    constants = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=ws, MWs=[c["mw"] for c in comps],
                                         CASs=[c["cas"] for c in comps])
    # A PT flash doesn't use Cp; FlashVL only needs the objects to exist.
    cps = [HeatCapacityGas(CASRN=c["cas"]) for c in comps]
    correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=cps, skip_missing=True)
    eos = dict(Tcs=tcs, Pcs=pcs, omegas=ws, kijs=kij)
    flasher = FlashVL(constants, correlations, liquid=CEOSLiquid(PRMIX, eos, HeatCapacityGases=cps),
                      gas=CEOSGas(PRMIX, eos, HeatCapacityGases=cps))

    def residual(y, x):
        g, l = PRMIX(T=T, P=P, zs=y, **eos), PRMIX(T=T, P=P, zs=x, **eos)
        return max(abs(math.log(a) - math.log(b)) for a, b in zip(g.fugacities_g, l.fugacities_l))

    print("\nthe same PT flash elsewhere (DWSIM at 1E-12 minus each):")
    for tol in (None, 1e-30):
        if tol:
            flasher.PT_SS_TOL = tol  # its error is a sum of squares, so 1e-13 is about 3e-7 in K
        res = flasher.flash(T=T, P=P, zs=Z)
        print(f"  thermo {__import__('thermo').__version__}, PT_SS_TOL {flasher.PT_SS_TOL:g}: beta {res.VF!r}, "
              f"|ln fV - ln fL| {residual(res.gas.zs, res.liquid0.zs):.3g}; DWSIM minus this: "
              f"beta {dw['beta'] - res.VF:.3g}, y {max(abs(a - b) for a, b in zip(dw['y'], res.gas.zs)):.3g}, "
              f"x {max(abs(a - b) for a, b in zip(dw['x'], res.liquid0.zs)):.3g}")
    for label, k in (("exact constants", EXACT), ("DWSIM's constants", DWSIM_PR)):
        beta, y, x = pr_flash(comps, kij, k)
        print(f"  PR here, {label:17}: beta {beta!r}; DWSIM minus this: beta {dw['beta'] - beta:.3g}, "
              f"y {max(abs(a - b) for a, b in zip(dw['y'], y)):.3g}")

    print("\nDWSIM's enthalpy, part by part:")
    for i, c in enumerate(comps):
        pure = [0.0] * len(comps)
        pure[i] = 1.0
        diffs = []
        for t in (250.0, 400.0, 600.0):
            mine = cp_integral(c["cp"], T_REF, t)
            diffs.append(f"{t:.0f} K {pp.RET_Hid(T_REF, t, vec(pure)) * c['mw'] - mine:+.3g}")
        print(f"  {c['name']:9} ideal gas from {T_REF} K, DWSIM minus exact quadrature (J/mol): {', '.join(diffs)}")
    for label, comp, vapor, st in (("vapour", dw["y"], True, State.Vapor), ("liquid", dw["x"], False, State.Liquid)):
        m = sum(a * c["mw"] for a, c in zip(comp, comps))
        theirs = pp.DW_CalcEnthalpyDeparture(vec(comp), T, P, st) * m
        exact = pr_phase(comps, kij, comp, T, P, vapor, EXACT)[2]
        rounded = pr_phase(comps, kij, comp, T, P, vapor, DWSIM_PR)[2]
        print(f"  {label} departure {theirs!r} J/mol; DWSIM minus PR with exact constants {theirs - exact:.3g} "
              f"({(theirs - exact) / abs(exact):.2g}), with DWSIM's {theirs - rounded:.3g} "
              f"({(theirs - rounded) / abs(rounded):.2g})")


if __name__ == "__main__":
    main()
