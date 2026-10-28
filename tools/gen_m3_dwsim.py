# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Generate testdata/m3_dwsim.json: DWSIM 9's flashes, for M3's exit criterion 6.

Run from the repository root:

    uv run --python 3.13 tools/gen_m3_dwsim.py

It needs Windows and DWSIM 9 installed for the current user, in
%LOCALAPPDATA%\\DWSIM, and pythonnet 3.0.5 stops at Python 3.13. CI never
runs it: like the CoolProp fixtures, it's run by hand and its output is
committed. It grew from spike/dwsim/run.py and does what "What step 5 must
do" in docs/spike-dwsim.md asks.

Each case has up to three answers:

dwsim
    DWSIM's "Peng-Robinson (PR)" package, called directly
    (CalculateEquilibrium), with its PT and PH tolerances at 1e-8, past
    which nothing changes. residual is DWSIM's own fugacity residual at its
    answer, computed with its own fugacity coefficients, so that its
    convergence floor shows here and not in reflux's error.
exact
    The same Peng-Robinson with DWSIM's constants and kij, converged to
    round-off: thermo 0.6.1's FlashVL at PT_SS_TOL 1e-30, and for H and S
    thermo's departures plus an exact quadrature (mpmath) of DWSIM's
    ideal-gas Cp. This is what reflux must reproduce to 1e-8.
rounded
    The same again with the constants DWSIM's PR actually uses: Omega_a
    0.45724, Omega_b 0.0778 and R 8.314. With DWSIM's own ideal-gas H and S
    (RET_Hid, RET_Sid) it reproduces DWSIM's H and S, which explains them;
    the script stops if it doesn't.

DWSIM's ideal-gas entropy has a mixing term of -R sum(z_i ln z_i) M/M_i,
where -R sum(z_i ln z_i) belongs (M is the mixture's molar mass): the mass
basis applied to mole fractions. The script measures it and records it, as
s_mixing_error. Within one program it cancels from a PS flash's round trip,
which is one reason PH and PS are compared by state: from the same (T, P)
and z, each program computes its own H or S and flashes back, and dwsim
records where DWSIM lands.
"""

import datetime
import functools
import json
import math
import os
import pathlib
import platform
import sys

import mpmath as mp
import thermo
from fluids.constants import R as R_EXACT
from thermo import PRMIX

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m3_fixtures as g3  # noqa: E402  the thermo helpers: flashers, residuals, polishing

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m3_dwsim.json"
T_REF = 298.15  # K, DWSIM's ideal-gas reference for H and S
P_REF = 101325.0  # Pa, DWSIM's ideal-gas reference pressure for S
R_DWSIM = 8.314
TOL = "1E-8"
mp.mp.dps = 30

SYSTEMS = [
    {"name": "methane/n-butane", "compounds": ["Methane", "N-butane"], "z": [0.4, 0.6]},
    {"name": "natural gas", "compounds": ["Methane", "Ethane", "Propane", "N-butane", "N-pentane"],
     "z": [0.70, 0.10, 0.08, 0.07, 0.05]},
    {"name": "CO2/n-butane", "compounds": ["Carbon dioxide", "N-butane"], "z": [0.5, 0.5]},
    {"name": "nitrogen/methane", "compounds": ["Nitrogen", "Methane"], "z": [0.3, 0.7]},
    {"name": "n-butane", "compounds": ["N-butane"], "z": [1.0]},
]
C1C4, GAS, CO2C4, N2C1, C4 = range(5)

CASES = [
    # PT: two-phase, vapour, liquid, near the critical locus, and pure.
    {"system": C1C4, "kind": "pt", "t": 250.0, "p": 30e5},
    {"system": C1C4, "kind": "pt", "t": 350.0, "p": 10e5},
    {"system": C1C4, "kind": "pt", "t": 220.0, "p": 60e5},
    {"system": C1C4, "kind": "pt", "t": 340.0, "p": 110e5},
    {"system": GAS, "kind": "pt", "t": 250.0, "p": 40e5},
    {"system": GAS, "kind": "pt", "t": 300.0, "p": 50e5},
    {"system": CO2C4, "kind": "pt", "t": 320.0, "p": 40e5},
    {"system": N2C1, "kind": "pt", "t": 150.0, "p": 30e5},
    {"system": C4, "kind": "pt", "t": 300.0, "p": 1e5},
    {"system": C4, "kind": "pt", "t": 300.0, "p": 10e5},
    # PH and PS, by state: (T, P), or (P, vapour fraction) for the pure
    # component's saturated state.
    {"system": C1C4, "kind": "ph", "t": 250.0, "p": 30e5},
    {"system": C1C4, "kind": "ph", "t": 350.0, "p": 10e5},
    {"system": GAS, "kind": "ph", "t": 250.0, "p": 40e5},
    {"system": CO2C4, "kind": "ph", "t": 320.0, "p": 40e5},
    {"system": C4, "kind": "ph", "p": 5e5, "vf": 0.5},
    {"system": C1C4, "kind": "ps", "t": 250.0, "p": 30e5},
    {"system": C1C4, "kind": "ps", "t": 220.0, "p": 60e5},
    {"system": GAS, "kind": "ps", "t": 250.0, "p": 40e5},
    {"system": C4, "kind": "ps", "p": 5e5, "vf": 0.5},
    # Bubble and dew points, at a given P or T.
    {"system": C1C4, "kind": "bubble", "p": 30e5},
    {"system": C1C4, "kind": "bubble", "t": 250.0},
    {"system": GAS, "kind": "bubble", "p": 40e5},
    {"system": C1C4, "kind": "dew", "p": 30e5},
    {"system": C1C4, "kind": "dew", "t": 300.0},
    {"system": GAS, "kind": "dew", "p": 40e5},
]


class PRMIXRounded(PRMIX):
    """thermo's PR with DWSIM's Omega_a and Omega_b. R cancels from A and B,
    so the phase equilibrium needs no change of R; departures scale with R,
    and rounded_dep converts them."""

    c1, c2 = 0.45724, 0.0778
    c1R2 = c1 * R_EXACT * R_EXACT
    c2R = c2 * R_EXACT
    c1R2_c2R = c1R2 / c2R


def rounded_dep(ph):
    """A PRMIXRounded phase's H and S departures with DWSIM's R."""
    return ph.H_dep() * R_DWSIM / R_EXACT, ph.S_dep() * R_DWSIM / R_EXACT


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


@functools.lru_cache(maxsize=None)
def ig_integrals(cp, t):
    """ChemSep eq 16, Cp = A + exp(B/T + C + D T + E T^2) in J/(kmol K),
    from T_REF to t: (int Cp dT in J/mol, int Cp/T dT in J/(mol K)), by
    mpmath at 30 digits."""
    a, b, c, d, e = (mp.mpf(v) for v in cp)

    def f(x):
        return a + mp.exp(b / x + c + d * x + e * x * x)

    return float(mp.quad(f, [T_REF, t])) / 1000, float(mp.quad(lambda x: f(x) / x, [T_REF, t])) / 1000


def exact_hs(comps, pairs, t, p):
    """Molar H and S of thermo phases (fraction, phase) on DWSIM's
    reference state, done right: the exact Cp integral, R exact, and ideal
    mixing -R sum(x ln x)."""
    h = s = 0.0
    for beta, ph in pairs:
        hi = si = 0.0
        for x, c in zip(ph.zs, comps):
            ih, isx = ig_integrals(tuple(c["cp"]), t)
            hi += x * ih
            si += x * isx - (R_EXACT * x * math.log(x) if x > 0 else 0.0)
        si -= R_EXACT * math.log(p / P_REF)
        h += beta * (hi + ph.H_dep())
        s += beta * (si + ph.S_dep())
    return h, s


def split_of(pairs, n):
    """beta, y and x from phases sorted by density: a single phase is the
    vapour (beta 1) or the liquid (beta 0) by thermo's label, with zeros for
    the absent phase, as DWSIM reports it."""
    if len(pairs) == 2:
        return pairs[0][0], list(pairs[0][1].zs), list(pairs[1][1].zs)
    ph = pairs[0][1]
    zero = [0.0] * n
    return (1.0, list(ph.zs), zero) if ph.__class__.__name__ == "CEOSGas" else (0.0, zero, list(ph.zs))


class Dwsim:
    """One DWSIM flowsheet and property package per system, with its PT and PH
    tolerances at tol, or at DWSIM's defaults for tol None (the benchmarks
    time both)."""

    def __init__(self, auto, system, tol=TOL):
        from System import Array, Double
        from DWSIM.Interfaces.Enums import FlashCalculationType, FlashSetting
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType
        from DWSIM.Thermodynamics.PropertyPackages import State

        self.FCT, self.State = FlashCalculationType, State
        self.vec = lambda xs: Array[Double](list(xs))
        sim = auto.CreateFlowsheet()
        for name in system["compounds"]:
            c = sim.AvailableCompounds[name]
            assert c.OriginalDB == "ChemSep", (name, c.OriginalDB)  # "n-Butane" is CoolProp's, "N-butane" ChemSep's
            sim.SelectedCompounds.Add(name, c)
        self.pp = sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)").__implementation__
        stream = sim.AddObject(ObjectType.MaterialStream, 50, 50, "feed").GetAsObject()
        stream.SetTemperature(300.0)
        stream.SetPressure(1e5)
        stream.SetMolarFlow(1.0)
        stream.SetOverallComposition(self.vec(system["z"]))
        self.pp.CurrentMaterialStream = stream  # the package reads its compounds from this
        for key in (FlashSetting.PTFlash_External_Loop_Tolerance, FlashSetting.PTFlash_Internal_Loop_Tolerance,
                    FlashSetting.PHFlash_External_Loop_Tolerance, FlashSetting.PHFlash_Internal_Loop_Tolerance):
            if tol is not None:
                self.pp.FlashSettings[key] = tol
        self.settings = {str(k): str(self.pp.FlashSettings[k]) for k in self.pp.FlashSettings.Keys}
        n = len(system["compounds"])
        self.comps = []
        for name in system["compounds"]:
            c = sim.SelectedCompounds[name]
            assert c.IdealgasCpEquation == "16", (name, c.IdealgasCpEquation)
            self.comps.append({
                "name": name, "cas": c.CAS_Number, "db": c.OriginalDB, "tc": c.Critical_Temperature,
                "pc": c.Critical_Pressure, "omega": c.Acentric_Factor, "mw": c.Molar_Weight / 1000,
                "cp_eq": 16, "cp": [c.Ideal_Gas_Heat_Capacity_Const_A, c.Ideal_Gas_Heat_Capacity_Const_B,
                                    c.Ideal_Gas_Heat_Capacity_Const_C, c.Ideal_Gas_Heat_Capacity_Const_D,
                                    c.Ideal_Gas_Heat_Capacity_Const_E],
                # The vapour-pressure correlation, which DWSIM's PR package
                # uses for a pure component's saturation instead of the EOS.
                "psat_eq": int(c.VaporPressureEquation),
                "psat": [c.Vapor_Pressure_Constant_A, c.Vapor_Pressure_Constant_B, c.Vapor_Pressure_Constant_C,
                         c.Vapor_Pressure_Constant_D, c.Vapor_Pressure_Constant_E],
            })
        k = self.pp.RET_VKij()
        self.kij = [[k.GetValue(i, j) for j in range(n)] for i in range(n)]
        self.mws = [sim.SelectedCompounds[name].Molar_Weight for name in system["compounds"]]

    def mw(self, x):
        return sum(a * b for a, b in zip(x, self.mws))  # g/mol, so kJ/kg times this is J/mol

    def flash(self, kind, a, b, z, estimate=0.0):
        r = self.pp.CalculateEquilibrium(getattr(self.FCT, kind), a, b, self.vec(z), None, estimate)
        out = {"t": r.CalculatedTemperature, "p": r.CalculatedPressure, "beta": r.GetVaporPhaseMoleFraction(),
               "y": list(r.GetVaporPhaseMoleFractions()), "x": list(r.GetLiquidPhase1MoleFractions()),
               "iterations": r.IterationsTaken}
        assert r.GetLiquidPhase2MoleFraction() == 0.0 and r.GetSolidPhaseMoleFraction() == 0.0
        out["k"] = list(r.Kvalues) if r.Kvalues is not None else []
        out["h_mass"], out["s_mass"] = r.CalculatedEnthalpy, r.CalculatedEntropy
        return out

    def ln_phi(self, x, t, p, state):
        return [math.log(v) for v in self.pp.DW_CalcFugCoeff(self.vec(x), t, p, state)]

    def residual(self, y, x, t, p):
        """DWSIM's own max |ln fV - ln fL| at a split."""
        lv, ll = self.ln_phi(y, t, p, self.State.Vapor), self.ln_phi(x, t, p, self.State.Liquid)
        return max(abs(math.log(a) + la - math.log(b) - lb) for a, la, b, lb in zip(y, lv, x, ll) if a > 0 and b > 0)

    def hs_phase(self, x, t, p, vapor):
        """DWSIM's molar H and S of one phase: its ideal gas plus its departure."""
        st = self.State.Vapor if vapor else self.State.Liquid
        m = self.mw(x)
        h = (self.pp.RET_Hid(T_REF, t, self.vec(x)) + self.pp.DW_CalcEnthalpyDeparture(self.vec(x), t, p, st)) * m
        s = (self.pp.RET_Sid(T_REF, t, p, self.vec(x)) + self.pp.DW_CalcEntropyDeparture(self.vec(x), t, p, st)) * m
        return h, s

    def mixing_error(self, x):
        """DWSIM's ideal mixing entropy minus -R sum(x ln x), J/(mol K), with R
        8.314: it uses -R sum(x_i ln x_i) M/M_i."""
        m = self.mw(x)
        dwsim = (self.pp.RET_Sid(T_REF, T_REF, P_REF, self.vec(x)) * m)
        return dwsim + R_DWSIM * sum(a * math.log(a) for a in x if a > 0)


def pt_case(dw, case, z):
    t, p = case["t"], case["p"]
    r = dw.flash("PressureTemperature", p, t, z)
    m = dw.mw(z)
    out = {"t": t, "p": p, "beta": r["beta"], "y": r["y"], "x": r["x"], "h": r["h_mass"] * m, "s": r["s_mass"] * m,
           "iterations": r["iterations"]}
    if 0.0 < r["beta"] < 1.0:
        out["residual"] = dw.residual(r["y"], r["x"], t, p)
    return out


def exact_pt(dw, t, p, z, eos_cls=PRMIX):
    f = g3.make_flasher(eos_cls, dw.comps, dw.kij)
    res, tol = g3.thermo_pt(f, t, p, z)
    return g3.by_density(res), tol


def run_pt(dw, case, z):
    t, p = case["t"], case["p"]
    d = pt_case(dw, case, z)
    pairs, tol = exact_pt(dw, t, p, z)
    beta, y, x = split_of(pairs, len(z))
    h, s = exact_hs(dw.comps, pairs, t, p)
    exact = {"beta": beta, "y": y, "x": x, "h": h, "s": s, "tol": tol}
    if len(pairs) == 2:
        exact["residual"] = g3.fugacity_residual(pairs[0][1], pairs[1][1])
    rpairs, _ = exact_pt(dw, t, p, z, PRMIXRounded)
    rb, ry, rx = split_of(rpairs, len(z))
    # Each of DWSIM's phases, rebuilt at its own composition: DWSIM's ideal
    # gas plus the rounded PR's departures. (Its mixture H can't be rebuilt
    # the same way: its y, x and beta balance the feed only to its
    # convergence, 7e-8 at the first case, and it weights the phases by
    # mass.) rounded's h and s are the mole-weighted sums.
    rh = rs = 0.0
    phases = [(d["beta"], d["y"], True), (1 - d["beta"], d["x"], False)]
    for frac, comp, vapor in phases:
        if frac == 0.0:
            continue
        m = dw.mw(comp)
        cls = g3.CEOSGas if vapor else g3.CEOSLiquid
        ph = cls(PRMIXRounded, dict(Tcs=[c["tc"] for c in dw.comps], Pcs=[c["pc"] for c in dw.comps],
                                    omegas=[c["omega"] for c in dw.comps], kijs=dw.kij),
                 HeatCapacityGases=None).to(T=t, P=p, zs=comp)
        hdep, sdep = rounded_dep(ph)
        ph_h = dw.pp.RET_Hid(T_REF, t, dw.vec(comp)) * m + hdep
        ph_s = dw.pp.RET_Sid(T_REF, t, p, dw.vec(comp)) * m + sdep
        theirs = dw.hs_phase(comp, t, p, vapor)
        for got, want, scale in ((ph_h, theirs[0], R_EXACT * t), (ph_s, theirs[1], R_EXACT)):
            err = abs(got - want) / (abs(want) + scale)
            assert err < 1e-12, f"the rounded PR doesn't reproduce DWSIM's phase: {err:.3g} at {case}"
        rh += frac * ph_h
        rs += frac * ph_s
    rounded = {"beta": rb, "y": ry, "x": rx, "h": rh, "s": rs}
    return {"dwsim": d, "exact": exact, "rounded": rounded,
            "s_mixing_error": sum(frac * dw.mixing_error(comp) for frac, comp, _ in phases if frac > 0)}


def run_by_state(dw, case, z):
    """PH or PS by state: each program computes its own H or S at the state
    and flashes back."""
    spec_key = "h" if case["kind"] == "ph" else "s"
    kind = "PressureEnthalpy" if spec_key == "h" else "PressureEntropy"
    p = case["p"]
    m = dw.mw(z)
    if "vf" in case:
        # The pure component's saturated state, at (P, vapour fraction).
        sat = dw.flash("PressureVaporFraction", p, case["vf"], z)
        t_dw = sat["t"]
        hl, sl = dw.hs_phase(z, t_dw, p, False)
        hv, sv = dw.hs_phase(z, t_dw, p, True)
        vf = case["vf"]
        spec = (1 - vf) * (hl if spec_key == "h" else sl) + vf * (hv if spec_key == "h" else sv)
        f = g3.make_flasher(PRMIX, dw.comps, dw.kij)
        res = f.flash(P=p, VF=vf, zs=z)
        pairs = g3.by_density(res)
        if len(pairs) == 1:
            raise RuntimeError("expected both saturated phases")
        t_state = res.T
        state = {"t": t_state, "beta": vf, "y": list(z), "x": list(z)}
        h, s = exact_hs(dw.comps, pairs, t_state, p)
        dwsim_state = {"t": t_dw, "iterations": sat["iterations"]}
    else:
        t_state = case["t"]
        d = pt_case(dw, case, z)
        spec = d[spec_key]
        pairs, _ = exact_pt(dw, t_state, p, z)
        beta, y, x = split_of(pairs, len(z))
        state = {"t": t_state, "beta": beta, "y": y, "x": x}
        h, s = exact_hs(dw.comps, pairs, t_state, p)
        dwsim_state = {"t": t_state, "beta": d["beta"], "y": d["y"], "x": d["x"]}
    back = dw.flash(kind, p, spec / m, z)
    landing = {"t": back["t"], "beta": back["beta"], "y": back["y"], "x": back["x"], "iterations": back["iterations"]}
    if 0.0 < back["beta"] < 1.0:
        landing["residual"] = dw.residual(back["y"], back["x"], back["t"], p)
    return {"dwsim": {"state": dwsim_state, "spec": spec, "landing": landing},
            "exact": {"state": state, "spec": h if spec_key == "h" else s}}


def run_boundary(dw, case, z):
    vf = 0.0 if case["kind"] == "bubble" else 1.0
    spec = "p" if "p" in case else "t"
    value = case[spec]
    f = g3.make_flasher(PRMIX, dw.comps, dw.kij)
    res = f.flash(zs=z, VF=vf, **{spec.upper(): value})
    feed, new = (res.liquid0, res.gas) if vf == 0.0 else (res.gas, res.liquid0)
    raw = g3.fugacity_residual(feed, new)
    t, p, w = g3.polish_boundary(f, z, vf, spec, res.T, res.P, list(new.zs))
    feed_phase, new_phase = (f.liquid, f.gas) if vf == 0.0 else (f.gas, f.liquid)
    residual = g3.fugacity_residual(feed_phase.to(T=t, P=p, zs=z), new_phase.to(T=t, P=p, zs=w))
    shift = (t - res.T) / res.T if spec == "p" else (p - res.P) / res.P
    problem = g3.boundary_problem(f, z, spec, t, p, w, residual, shift)
    assert problem is None, (case, problem)
    exact = {"t": t, "p": p, "incipient": w, "residual": residual, "thermo_residual": raw}
    rf = g3.make_flasher(PRMIXRounded, dw.comps, dw.kij)
    rres = rf.flash(zs=z, VF=vf, **{spec.upper(): value})
    rt, rp, _ = g3.polish_boundary(rf, z, vf, spec, rres.T, rres.P, list((rres.gas if vf == 0.0 else rres.liquid0).zs))
    kind = "PressureVaporFraction" if spec == "p" else "TemperatureVaporFraction"
    out = {"exact": exact, "rounded": {"t": rt, "p": rp}}
    try:
        r = dw.flash(kind, value, vf, z)
        estimate = None
    except Exception as e:
        # DWSIM's bubble point at a given T fails from its own start
        # ("Unable to calculate compressility factor"), so give it one.
        first = f"{type(e).__name__}: {str(e).splitlines()[0]}"[:200]
        estimate = p * 1.05 if spec == "t" else t * 1.02
        try:
            r = dw.flash(kind, value, vf, z, estimate)
        except Exception as e2:
            out["dwsim"] = {"error": first, "with_estimate": f"{type(e2).__name__}: {str(e2).splitlines()[0]}"[:200]}
            return out
        out["dwsim_without_estimate"] = first
    # DWSIM reports the incipient phase only through its K-values.
    ks = r["k"]
    w = [k * zi for k, zi in zip(ks, z)] if vf == 0.0 else [zi / k for k, zi in zip(ks, z)]
    total = math.fsum(w)
    wn = [v / total for v in w]
    feed_state, new_state = (dw.State.Liquid, dw.State.Vapor) if vf == 0.0 else (dw.State.Vapor, dw.State.Liquid)
    lf = dw.ln_phi(z, r["t"], r["p"], feed_state)
    ln = dw.ln_phi(wn, r["t"], r["p"], new_state)
    res_dw = max(abs(math.log(a) + la - math.log(b) - lb) for a, la, b, lb in zip(z, lf, wn, ln))
    out["dwsim"] = {"t": r["t"], "p": r["p"], "incipient": wn, "sum": total, "residual": res_dw,
                    "iterations": r["iterations"]}
    if estimate is not None:
        out["dwsim"]["estimate"] = estimate
    return out


def main():
    load_dwsim()
    from DWSIM.Automation import Automation3

    auto = Automation3()
    dws = [Dwsim(auto, s) for s in SYSTEMS]
    systems = [{"name": s["name"], "components": dw.comps, "kij": [x for row in dw.kij for x in row], "z": s["z"]}
               for s, dw in zip(SYSTEMS, dws)]
    cases = []
    for case in CASES:
        dw, z = dws[case["system"]], SYSTEMS[case["system"]]["z"]
        out = dict(case)
        if case["kind"] == "pt":
            out.update(run_pt(dw, case, z))
        elif case["kind"] in ("ph", "ps"):
            out.update(run_by_state(dw, case, z))
        else:
            out.update(run_boundary(dw, case, z))
        cases.append(out)
        print(f"  {case}: done", flush=True)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_dwsim.py",
            "python": platform.python_version(),
            "dwsim": str(auto.GetVersion()),
            "thermo": thermo.__version__,
            "mpmath": mp.__version__,
        },
        "note": "DWSIM 9's Peng-Robinson (PR) package against the same PR with DWSIM's constants and kij. dwsim is "
        "DWSIM's answer at PT and PH tolerances 1e-8, residual its own max |ln fV - ln fL|; exact is the PR "
        "converged to round-off (thermo, and mpmath for the ideal-gas Cp integrals); rounded is the PR with "
        "DWSIM's Omega_a, Omega_b and R. h in J/mol relative to the ideal gas at t_ref, s in J/(mol K) relative to "
        "the ideal gas at t_ref and p_ref_s; cp is ChemSep eq 16 in J/(kmol K); mw in kg/mol. Single phases have "
        "beta 1 (vapour) or 0 (liquid), and zeros for the absent phase. PH and PS are by state: spec is each "
        "program's own H or S at the state, and landing is where DWSIM's flash from it ends.",
        "t_ref": T_REF,
        "p_ref_s": P_REF,
        "dwsim_pr": {"omega_a": PRMIXRounded.c1, "omega_b": PRMIXRounded.c2, "r": R_DWSIM},
        "flash_settings": dws[0].settings,
        "systems": systems,
        "cases": cases,
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    report(data)


def report(data):
    print(f"wrote {OUT.relative_to(ROOT)}: {len(data['cases'])} cases, DWSIM {data['source']['dwsim']}")
    for c in data["cases"]:
        d, e = c.get("dwsim", {}), c.get("exact", {})
        if c["kind"] == "pt":
            dev = max([abs(d["beta"] - e["beta"])] + [abs(a - b) for a, b in zip(d["y"] + d["x"], e["y"] + e["x"])])
            print(f"  pt {c['system']} {c['t']} K {c['p'] / 1e5:g} bar: split {dev:.2g} (rounded "
                  f"{abs(d['beta'] - c['rounded']['beta']):.2g}), H {(d['h'] - e['h']) / abs(e['h']):.2g}, "
                  f"S {d['s'] - e['s']:.3g} J/(mol K) (mixing {c['s_mixing_error']:.3g}), residual "
                  f"{d.get('residual', 0):.2g}")
        elif c["kind"] in ("ph", "ps"):
            land, st = d["landing"], e["state"]
            print(f"  {c['kind']} {c['system']}: DWSIM lands {land['t'] - st['t']:.3g} K from the exact state "
                  f"({land['t'] - d['state']['t']:.3g} K from its own), beta {land['beta'] - st['beta']:.2g}; spec "
                  f"DWSIM minus exact {d['spec'] - e['spec']:.3g}")
        elif "error" in d:
            print(f"  {c['kind']} {c['system']}: DWSIM failed: {d['error']}; {d['with_estimate']}")
        else:
            r = c["rounded"]
            print(f"  {c['kind']} {c['system']}: DWSIM minus exact T {d['t'] - e['t']:.3g} K, P "
                  f"{(d['p'] - e['p']) / e['p']:.2g}; minus rounded T {d['t'] - r['t']:.2g} K, P "
                  f"{(d['p'] - r['p']) / r['p']:.2g}; residual {d['residual']:.2g}, sum {d['sum'] - 1:.2g}"
                  + (f", estimate {d['estimate']:.6g}" if "estimate" in d else ""))


if __name__ == "__main__":
    main()
