# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""Generate testdata/m9_dwsim.json: DWSIM 9's reactors, for M9 step 10
(D13, criteria 7 and 8).

Run from the repository root, on Windows with DWSIM 9 installed:

    uv run --python 3.13 tools/gen_m9_dwsim.py [part ...]

with parts among conversion, equilibrium, gibbs, cstr, pfr and ideal_gas; with none,
all of them. A part not named keeps what the file has. pfr takes about a
quarter of an hour, two of its solves at dV = 0.001.

Each part builds its case here from numbers, as M9 step 0's spike does
(spike/dwsim_pfr/run.py, whose PFR, rate model and viscosity this imports),
and checks every reactor's Calculated and ErrorMessage: DWSIM's failures are
silent. REFLUX_TUTORIAL names the owner's reforming file and
REFLUX_TUTORIAL2 the methanol one; if set, each is loaded, solved and
checked to be the case built here, and the check is printed, not recorded.
Neither file is in the repository.

conversion
    The owner's Tutorial 2: CO + 2 H2 = CH3OH in an isothermal conversion
    reactor on DWSIM's PR78, CO's conversion 85 - 0.18*(T - 273.15) % as
    DWSIM takes it, 1500 kg/h of CO and H2 at 0.2 and 0.8. The base case at
    320 C and 70 bar, the temperature sweep (200 to 320 C at 70 bar) and the
    pressure sweep (30 and 50 bar at 320 C). Each point's duty taken apart:
    DWSIM's own departures and ideal-gas enthalpies (RET_Hid) and its heat of
    reaction, which must add up to its duty; the same parts as reflux has
    them, from thermo's PR with exact constants and an exact quadrature of
    reflux's Cp coefficients; and thermo's PR78 with DWSIM's rounded
    constants, which explains most of DWSIM's departures. And the adiabatic
    reactor, whose outlet shows the T DWSIM evaluates its expression at.
equilibrium
    The same reaction in DWSIM's equilibrium reactor, K from the Gibbs
    energies, at 320 C and 30, 50 and 70 bar on PR, at DWSIM's default
    tolerances and at 1e-10. Beside it, the equilibrium as reflux defines
    it (fugacities in bar, K at 1 bar from the component files' data, PR
    with exact constants), and with fugacities in atm against the same K on
    DWSIM's PR, which is what DWSIM's answers turn out to be.
gibbs
    The owner's reforming feed in DWSIM's Gibbs reactor on PR, at 1000 K
    and 1 atm and at 900 K and 10 bar, beside the same two equilibria
    computed here both ways.
cstr
    DWSIM 9.0.5's CSTR reacts nothing in the vapour, silently (recorded), so
    the reference is the exact answer of DWSIM's own equations: the owner's
    three rates as typed, in a CSTR of 1 m3 with 0.386 kg of catalyst per
    m3, isothermal at 1000 K and 1 atm, solved by mpmath at 30 digits.
pfr
    The owner's reforming PFR (criterion 8): DWSIM at the file's dV = 0.05
    and at 0.001, isobaric and with Ergun; and the exact answer of DWSIM's
    equations without its sections, by SciPy's Radau at 1e-12: isobaric,
    and with Ergun on DWSIM's viscosity method and PR's density with exact
    constants and reflux's kij, which is reflux's density. DWSIM's own
    density uses R = 8.314, which moves its velocity by 5.6e-5.
ideal_gas
    Each reforming compound's ideal-gas enthalpy at the PFR's 1000 K,
    DWSIM's and reflux's, with DWSIM's Cp coefficients: its water is eq 16
    with other constants, where reflux's file has ChemSep's eq 107, 72 J/mol
    apart at 1000 K, which is most of the PFR's duty's difference.

Units throughout: K, Pa, mol/s, W, J/mol. The output is the same, byte for
byte, for the same versions, but for "generated" and the times.
"""

import datetime
import json
import math
import os
import pathlib
import platform
import sys
import time
import tomllib

import mpmath as mp
import numpy as np
from fluids.constants import R as R_EXACT
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from thermo import PR78MIX, PRMIX

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "spike" / "dwsim_pfr"))
import run as spike  # noqa: E402  DWSIM's PFR, its rate model and viscosity (M9 step 0)

OUT = ROOT / "testdata" / "m9_dwsim.json"
PARTS = set(sys.argv[1:]) or {"conversion", "equilibrium", "gibbs", "cstr", "pfr", "ideal_gas"}
TUTORIAL2 = os.environ.get("REFLUX_TUTORIAL2")
R_DWSIM = 8.314
mp.mp.dps = 30

CO, H2, MEOH = "Carbon monoxide", "Hydrogen", "Methanol"
METHANOL = {"names": [CO, H2, MEOH], "cas": ["630-08-0", "1333-74-0", "67-56-1"]}
T2 = {"t": 593.15, "p": 70e5, "mass": 1500.0 / 3600.0, "z": [0.2, 0.8, 0.0],
      "expression": "85-0.18*(T-273.15)"}
SYNTHESIS = {CO: -1.0, H2: -2.0, MEOH: 1.0}
REFORMING = {"names": spike.COMPS, "cas": spike.CAS}


# Component data, as reflux's component files have them.

def component(cas):
    return tomllib.loads((ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8"))


def reflux_kij(cas):
    """reflux's PR kij: the first entry pr.ipd lists for each pair."""
    k = [[0.0] * len(cas) for _ in cas]
    seen = set()
    for line in (ROOT / "data" / "vendor" / "chemsep" / "pr.ipd").read_text(encoding="latin-1").splitlines():
        parts = line.split()
        if len(parts) < 3 or parts[0] not in cas or parts[1] not in cas:
            continue
        i, j = cas.index(parts[0]), cas.index(parts[1])
        if (min(i, j), max(i, j)) in seen:
            continue
        seen.add((min(i, j), max(i, j)))
        k[i][j] = k[j][i] = float(parts[2])
    return k


def cp_integrals(cor, t, t0=298.15):
    """A component file's ideal-gas Cp, ChemSep's eq 16 or eq 107 (Aly and
    Lee's), from t0 to t, exactly: (∫Cp dT in J/mol, ∫Cp/T dT in J/(mol K))."""
    a, b, c, d, e = (mp.mpf(v) for v in cor["coeffs"])
    if cor["eq"] == 16:
        def cp(x):
            return (a + mp.exp(b / x + c + d * x + e * x * x)) / 1000
    else:
        assert cor["eq"] == 107, cor["eq"]

        def cp(x):
            return (a + b * (c / x / mp.sinh(c / x)) ** 2 + d * (e / x / mp.cosh(e / x)) ** 2) / 1000

    return mp.quad(cp, [t0, t]), mp.quad(lambda x: cp(x) / x, [t0, t])


def ln_k(cas, nu, t):
    """ln K at 1 bar from the component files' dhf, s0 and Cp, as reflux has
    it (D7): -ΔG°/RT, with H and S along Cp from 298.15 K."""
    dh = ds = mp.mpf(0)
    for c, v in zip(cas, nu):
        if v == 0:
            continue
        doc = component(c)
        h, s = cp_integrals(doc["correlations"]["ideal_gas_cp"], t)
        dh += v * (mp.mpf(doc["constants"]["dhf"]) + h)
        ds += v * (mp.mpf(doc["constants"]["s0"]) + s)
    return float(-(dh - mp.mpf(t) * ds) / (mp.mpf(R_EXACT) * t))


# Thermo's Peng-Robinson, as reflux has it and as DWSIM has it.

class Rounded:
    """DWSIM's PR constants: Omega_a 0.45724, Omega_b 0.0778 (and R 8.314,
    which only scales departures, applied by the caller)."""
    c1, c2 = 0.45724, 0.0778
    c1R2 = c1 * R_EXACT * R_EXACT
    c2R = c2 * R_EXACT
    c1R2_c2R = c1R2 / c2R


class PRRounded(Rounded, PRMIX):
    pass


class PR78Rounded(Rounded, PR78MIX):
    pass


class Fluid:
    """A PR mixture: its constants, kij and equation."""

    def __init__(self, cas, kij, cls=PRMIX, scale=1.0):
        self.docs = [component(c) for c in cas]
        self.tcs = [d["constants"]["tc"] for d in self.docs]
        self.pcs = [d["constants"]["pc"] for d in self.docs]
        self.omegas = [d["constants"]["omega"] for d in self.docs]
        self.mws = [d["constants"]["mw"] for d in self.docs]
        self.kij, self.cls, self.scale = kij, cls, scale

    def eos(self, t, p, n):
        tot = sum(n)
        return self.cls(Tcs=self.tcs, Pcs=self.pcs, omegas=self.omegas, kijs=self.kij, zs=[x / tot for x in n],
                        T=t, P=p)

    @staticmethod
    def gas(eos):
        """The one root, which thermo may label liquid (M9 step 9)."""
        return "g" if hasattr(eos, "Z_g") else "l"

    def departure(self, t, p, n):
        """H - H_ideal of the flow n, in W."""
        e = self.eos(t, p, n)
        return getattr(e, f"H_dep_{self.gas(e)}") * sum(n) * self.scale

    def ln_phi(self, t, p, n):
        e = self.eos(t, p, n)
        return list(e.fugacity_coefficients(getattr(e, f"Z_{self.gas(e)}")))

    def z(self, t, p, n):
        e = self.eos(t, p, n)
        return getattr(e, f"Z_{self.gas(e)}")


def equilibrium(fluid, nu_rows, ln_ks, feed, t, p, p_ref, guess):
    """The extents at which each reaction's ln Q, from fugacities over p_ref,
    is its ln K: solved by SciPy from a guess near the answer (DWSIM's), then
    polished by Newton to 1e-13. A step that would leave a species at or
    below zero is refused by a large residual."""
    def residual(x):
        n = list(feed)
        for row, xi in zip(nu_rows, x):
            n = [a + v * xi for a, v in zip(n, row)]
        if min(n) <= 0:
            return [1e6] * len(x)
        tot = sum(n)
        lp = fluid.ln_phi(t, p, n)
        return [sum(v * (math.log(n[i] / tot) + lp[i] + math.log(p / p_ref)) for i, v in enumerate(row) if v)
                - k for row, k in zip(nu_rows, ln_ks)]

    x = [float(v) for v in fsolve(residual, guess, xtol=1e-14, full_output=False)]
    for _ in range(30):
        r = residual(x)
        size = max(abs(v) for v in r)
        if size < 1e-13:
            break
        jac = np.zeros((len(x), len(x)))
        for j in range(len(x)):
            xp = list(x)
            xp[j] += 1e-7 * max(1e-3, abs(x[j]))
            jac[:, j] = [(a - b) / (xp[j] - x[j]) for a, b in zip(residual(xp), r)]
        dx = np.linalg.solve(jac, [-v for v in r])
        step = 1.0
        while step > 1e-6:
            trial = [a + step * float(d) for a, d in zip(x, dx)]
            if max(abs(v) for v in residual(trial)) < size:
                x = trial
                break
            step /= 2
        else:
            break
    assert max(abs(v) for v in residual(x)) < 1e-12, residual(x)
    n = list(feed)
    for row, xi in zip(nu_rows, x):
        n = [a + v * xi for a, v in zip(n, row)]
    return n


# DWSIM.

def dictionary(pairs):
    from System import Double, String
    from System.Collections.Generic import Dictionary

    d = Dictionary[String, Double]()
    for k, v in pairs.items():
        d[k] = v
    return d


class Reactor:
    """A feed, one reactor, its outlets and an energy stream, in a new
    flowsheet: the vapour outlet on port 0 and, where there is one, the
    liquid on port 1."""

    def __init__(self, auto, names, package, kind):
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType as OT
        from System import Array, Double

        self.auto, self.names = auto, names
        self.vec = lambda xs: Array[Double](list(xs))
        self.sim = auto.CreateFlowsheet()
        for c in names:
            compound = self.sim.AvailableCompounds[c]
            assert compound.OriginalDB == "ChemSep", (c, compound.OriginalDB)
            self.sim.SelectedCompounds.Add(c, compound)
        self.pp = self.sim.CreateAndAddPropertyPackage(package).__implementation__
        self.x = 0
        add = self.add
        self.feed = add(OT.MaterialStream, "feed")
        self.r = add(getattr(OT, kind), "reactor")
        self.top = add(OT.MaterialStream, "vapour")
        self.q = add(OT.EnergyStream, "q")
        self.sim.ConnectObjects(self.feed.GraphicObject, self.r.GraphicObject, -1, 0)
        self.sim.ConnectObjects(self.r.GraphicObject, self.top.GraphicObject, 0, -1)
        self.bottom = None
        if kind in ("RCT_Conversion", "RCT_Equilibrium", "RCT_Gibbs"):
            self.bottom = add(OT.MaterialStream, "liquid")
            self.sim.ConnectObjects(self.r.GraphicObject, self.bottom.GraphicObject, 1, -1)
        try:
            self.sim.ConnectObjects(self.r.GraphicObject, self.q.GraphicObject, 2, -1)
        except Exception:
            self.sim.ConnectObjects(self.q.GraphicObject, self.r.GraphicObject, -1, 1)
        self.pp.CurrentMaterialStream = self.feed.GetAsObject()
        k = self.pp.RET_VKij()
        self.kij = [[k.GetValue(i, j) for j in range(len(names))] for i in range(len(names))]

    def add(self, kind, name):
        self.x += 100
        return self.sim.AddObject(kind, self.x, 100, name)

    @property
    def o(self):
        return self.r.GetAsObject()

    def reaction(self, rx):
        self.sim.AddReaction(rx)
        self.sim.AddReactionToSet(rx.ID, "DefaultSet", True, 0)
        self.o.ReactionSetID = "DefaultSet"
        return rx

    def mode(self, name):
        spike.enum(self.o, "ReactorOperationMode", name)

    def set_feed(self, t, p, flows=None, mass=None, z=None):
        f = self.feed.GetAsObject()
        f.SetTemperature(t)
        f.SetPressure(p)
        if mass is not None:
            f.SetMassFlow(mass)
            f.SetOverallComposition(self.vec(z))
        else:
            f.SetMolarFlow(sum(flows))
            f.SetOverallComposition(self.vec(x / sum(flows) for x in flows))

    def solve(self):
        t0 = time.perf_counter()
        errs = self.auto.CalculateFlowsheet2(self.sim)
        dt = time.perf_counter() - t0
        bad = [str(e) for e in errs] if errs else []
        o = self.o
        if not o.Calculated or (o.ErrorMessage or ""):
            bad.append(f"reactor: calculated {o.Calculated}, {o.ErrorMessage!r}")
        if bad:
            raise RuntimeError(f"DWSIM failed: {bad}")
        return dt

    def flows(self, ms):
        """A stream's flows, and whether DWSIM left them NaN, as it leaves an
        empty liquid outlet at some states."""
        m = ms.GetAsObject().Phases[0]
        f = [m.Compounds[c].MolarFlow for c in self.names]
        nan = any(x is None or math.isnan(x) for x in f)
        return [0.0 if x is None or math.isnan(x) else x for x in f], nan

    def outlet(self):
        top, _ = self.flows(self.top)
        out = {"t": self.top.GetAsObject().Phases[0].Properties.temperature,
               "p": self.top.GetAsObject().Phases[0].Properties.pressure, "flows": top}
        if self.bottom is not None:
            bottom, nan = self.flows(self.bottom)
            assert max(bottom) == 0.0, f"a liquid: {bottom}"
            out["liquid_nan"] = nan
        out["duty"] = (self.o.DeltaQ or 0.0) * 1000.0
        return out

    def conversions(self):
        c = self.o.ComponentConversions
        return {k: (c[k] if math.isfinite(c[k]) else str(c[k]).replace("inf", "Infinity"))
                for k in self.names if c.ContainsKey(k)}

    def parts(self, t, p, n):
        """DWSIM's departure and ideal-gas enthalpy of a flow, in W."""
        from DWSIM.Thermodynamics.PropertyPackages import State

        tot = sum(n)
        x = [v / tot for v in n]
        m = sum(a * self.sim.SelectedCompounds[c].Molar_Weight for a, c in zip(x, self.names))
        self.pp.CurrentMaterialStream = self.feed.GetAsObject()  # a solve leaves it unset
        dep = self.pp.DW_CalcEnthalpyDeparture(self.vec(x), t, p, State.Vapor) * m * tot
        ig = self.pp.RET_Hid(298.15, t, self.vec(x)) * m * tot
        return dep, ig


def load():
    return spike.load()


# The parts.

def conversion(auto):
    def sheet(mode):
        sh = Reactor(auto, METHANOL["names"], "Peng-Robinson 1978 (PR78)", "RCT_Conversion")
        rx = sh.reaction(sh.sim.CreateConversionReaction("Methanol Synthesis", "", dictionary(SYNTHESIS), CO,
                                                         "Vapor", T2["expression"]))
        sh.mode(mode)
        return sh, rx

    sh, rx = sheet("Isothermic")
    cas = METHANOL["cas"]
    exact = Fluid(cas, reflux_kij(cas))
    pr78 = Fluid(cas, sh.kij, PR78Rounded, R_DWSIM / R_EXACT)
    docs = [component(c) for c in cas]
    nu = [SYNTHESIS[c] for c in METHANOL["names"]]
    points = []
    for t, p in [(T2["t"], T2["p"]), (473.15, 70e5), (513.15, 70e5), (553.15, 70e5), (593.15, 30e5), (593.15, 50e5)]:
        sh.set_feed(t, p, mass=T2["mass"], z=T2["z"])
        dt = sh.solve()
        feed, _ = sh.flows(sh.feed)
        out = sh.outlet()
        extent = feed[0] - out["flows"][0]
        dep_in, ig_in = sh.parts(t, p, feed)
        dep_out, ig_out = sh.parts(t, p, out["flows"])
        closure = dep_out - dep_in + ig_out - ig_in + rx.ReactionHeat * extent - out["duty"]
        assert abs(closure) <= 1e-6 * abs(out["duty"]), f"DWSIM's duty isn't its parts: {closure} W"
        ig_reflux = sum(v * extent * float(cp_integrals(d["correlations"]["ideal_gas_cp"], t)[0])
                        for v, d in zip(nu, docs))
        points.append({
            "t": t, "p": p, "feed": feed, "outlet": out, "conversions": sh.conversions(), "extent": extent,
            "dwsim": {"departure": dep_out - dep_in, "ideal_gas": ig_out - ig_in,
                      "reaction_heat": rx.ReactionHeat},
            "thermo": {"departure_pr": exact.departure(t, p, out["flows"]) - exact.departure(t, p, feed),
                       "departure_pr78_dwsim_constants": pr78.departure(t, p, out["flows"])
                       - pr78.departure(t, p, feed),
                       "ideal_gas_reflux_cp": ig_reflux},
            "seconds": round(dt, 3),
        })
        print(f"  conversion at {t} K, {p / 1e5:g} bar: duty {out['duty']:.6f} W", flush=True)
    ad, _ = sheet("Adiabatic")
    ad.set_feed(T2["t"], T2["p"], mass=T2["mass"], z=T2["z"])
    ad.solve()
    adiabatic = {"feed": ad.flows(ad.feed)[0], "outlet": ad.outlet(), "conversions": ad.conversions()}
    print(f"  adiabatic: outlet {adiabatic['outlet']['t']:.4f} K at CO conversion {adiabatic['conversions'][CO]}")
    if TUTORIAL2:
        check_tutorial2(auto, points[0])
    return {"package": "Peng-Robinson 1978 (PR78)", "expression": T2["expression"], "mass_flow": T2["mass"],
            "z": T2["z"], "names": METHANOL["names"], "cas": cas, "kij": sh.kij, "points": points,
            "adiabatic": adiabatic}


def check_tutorial2(auto, base):
    fs = auto.LoadFlowsheet(TUTORIAL2)
    objs = {o.GraphicObject.Tag: o.__implementation__ for o in fs.SimulationObjects.Values}
    rx = next(r for r in fs.Reactions.Values)
    reactor = next(o for o in objs.values() if o.GetType().Name == "Reactor_Conversion")
    same = (rx.Expression == T2["expression"] and rx.BaseReactant == CO
            and {c: rx.Components[c].StoichCoeff for c in rx.Components.Keys} == SYNTHESIS
            and list(fs.PropertyPackages.Values)[0].__implementation__.ComponentName == "Peng-Robinson 1978 (PR78)")
    errs = auto.CalculateFlowsheet2(fs)
    top = objs[reactor.GraphicObject.OutputConnectors[0].AttachedConnector.AttachedTo.Tag]
    flows = [top.Phases[0].Compounds[c].MolarFlow for c in METHANOL["names"]]
    worst = max(abs(a - b) / max(abs(b), 1e-300) for a, b in zip(flows, base["outlet"]["flows"]))
    print(f"  the owner's Tutorial 2: the same reaction and package {same}, errors {list(errs) if errs else []}, "
          f"outlet within {worst:.1e} of this script's, duty {reactor.DeltaQ * 1000:.6f} W against "
          f"{base['outlet']['duty']:.6f}")


def equilibrium_part(auto):
    sh = Reactor(auto, METHANOL["names"], "Peng-Robinson (PR)", "RCT_Equilibrium")
    rx = sh.reaction(sh.sim.CreateEquilibriumReaction("Methanol Synthesis", "", dictionary(SYNTHESIS), CO, "Vapor",
                                                      "Fugacity", "", 0.0, ""))
    assert str(rx.KExprType) == "Gibbs", rx.KExprType
    sh.mode("Isothermic")
    cas = METHANOL["cas"]
    exact = Fluid(cas, reflux_kij(cas))
    dwsim = Fluid(cas, sh.kij, PRRounded)
    nu = [SYNTHESIS[c] for c in METHANOL["names"]]
    default_tol = (sh.o.ExternalLoopTolerance, sh.o.InternalLoopTolerance)
    points = []
    for p in [30e5, 50e5, 70e5]:
        row = {"t": T2["t"], "p": p}
        for label, tol in [("default", None), ("tight", 1e-10)]:
            if tol:
                sh.o.ExternalLoopTolerance = sh.o.InternalLoopTolerance = tol
            else:
                sh.o.ExternalLoopTolerance, sh.o.InternalLoopTolerance = default_tol
            sh.set_feed(T2["t"], p, mass=T2["mass"], z=T2["z"])
            sh.solve()
            row[label] = {"outlet": sh.outlet(), "conversions": sh.conversions()}
        feed, _ = sh.flows(sh.feed)
        row["feed"] = feed
        lk = ln_k(cas, nu, T2["t"])
        guess = [row["tight"]["outlet"]["flows"][2]]
        row["ln_k"] = lk
        row["bar"] = equilibrium(exact, [nu], [lk], feed, T2["t"], p, 1e5, guess)
        row["atm"] = equilibrium(dwsim, [nu], [lk], feed, T2["t"], p, 101325.0, guess)
        points.append(row)
        x = lambda n: (feed[0] - n[0]) / feed[0]  # noqa: E731
        print(f"  equilibrium at {p / 1e5:g} bar: DWSIM {row['tight']['conversions'][CO]:.6f}, in atm "
              f"{x(row['atm']):.6f}, in bar {x(row['bar']):.6f}", flush=True)
    return {"package": "Peng-Robinson (PR)", "names": METHANOL["names"], "cas": cas, "kij": sh.kij,
            "default_tolerances": list(default_tol), "points": points}


def imbalance(names, feed, out):
    """The largest of C, H and O's imbalances, over what the feed has."""
    atoms = {"Methane": (1, 4, 0), "Water": (0, 2, 1), "Hydrogen": (0, 2, 0), "Carbon monoxide": (1, 0, 1),
             "Carbon dioxide": (1, 0, 2)}
    worst = 0.0
    for e in range(3):
        a = sum(f * atoms[c][e] for f, c in zip(feed, names))
        b = sum(f * atoms[c][e] for f, c in zip(out, names))
        worst = max(worst, abs(b - a) / a)
    return worst


def affinities(fluid, rows, lks, n, t, p, p_ref):
    """Each reaction's ln(Q/K) at n: zero at equilibrium."""
    tot = sum(n)
    lp = fluid.ln_phi(t, p, n)
    return [sum(v * (math.log(n[i] / tot) + lp[i] + math.log(p / p_ref)) for i, v in enumerate(row) if v) - k
            for row, k in zip(rows, lks)]


def gibbs(auto):
    sh = Reactor(auto, REFORMING["names"], "Peng-Robinson (PR)", "RCT_Gibbs")
    sh.mode("Isothermic")
    for c in REFORMING["names"]:
        if not sh.o.ComponentIDs.Contains(c):
            sh.o.ComponentIDs.Add(c)
    cas = REFORMING["cas"]
    exact = Fluid(cas, reflux_kij(cas))
    dwsim = Fluid(cas, sh.kij, PRRounded)
    # Two independent reactions span the equilibrium of these five species.
    names = REFORMING["names"]
    rows = [[{"Methane": -1, "Water": -1, "Carbon monoxide": 1, "Hydrogen": 3}.get(c, 0) for c in names],
            [{"Carbon monoxide": -1, "Water": -1, "Carbon dioxide": 1, "Hydrogen": 1}.get(c, 0) for c in names]]
    feed = list(spike.FEED)
    points = []
    for t, p in [(1000.0, 101325.0), (900.0, 1e6)]:
        sh.set_feed(t, p, flows=feed)
        sh.feed.GetAsObject().Calculate(True, True)
        sh.o.CreateElementMatrix()
        lks = [ln_k(cas, r, t) for r in rows]
        answers = {}
        for solver, ipopt in [("ipopt", True), ("other", False)]:
            sh.o.UseIPOPTSolver = ipopt
            # Its first solve from a feed without CO or CO2 refuses with
            # "invalid initial estimates", and a second, from what the first
            # left, runs.
            retried = False
            try:
                dt = sh.solve()
            except RuntimeError as e:
                assert "invalid initial estimates" in str(e), e
                retried = True
                dt = sh.solve()
            out = sh.outlet()
            answers[solver] = {"outlet": out, "retried": retried, "seconds": round(dt, 3),
                               "element_imbalance": imbalance(names, feed, out["flows"]),
                               "ln_q_over_k_atm": affinities(dwsim, rows, lks, out["flows"], t, p, 101325.0)}
        got = answers["other"]["outlet"]["flows"]
        co, co2 = got[names.index("Carbon monoxide")], got[names.index("Carbon dioxide")]
        guess = [co + co2, co2]  # the methane reformed, and the CO shifted
        points.append({
            "t": t, "p": p, "feed": feed, "dwsim": answers, "elements": list(sh.o.Elements), "ln_k": lks,
            "bar": equilibrium(exact, rows, lks, feed, t, p, 1e5, guess),
            "atm": equilibrium(dwsim, rows, lks, feed, t, p, 101325.0, guess),
        })
        print(f"  Gibbs at {t} K, {p / 1e5:g} bar: DWSIM's IPOPT {answers['ipopt']['outlet']['flows']} "
              f"(elements off by {answers['ipopt']['element_imbalance']:.1e}), its other method "
              f"{answers['other']['outlet']['flows']} ({answers['other']['element_imbalance']:.1e}); atm "
              f"{points[-1]['atm']}; bar {points[-1]['bar']}", flush=True)
    return {"package": "Peng-Robinson (PR)", "names": names, "cas": cas, "kij": sh.kij, "reactions": rows,
            "points": points}


def cstr(auto):
    # DWSIM's CSTR, which reacts nothing in the vapour: recorded.
    sh = Reactor(auto, REFORMING["names"], "Peng-Robinson (PR)", "RCT_CSTR")
    for r in spike.REACTIONS:
        d = dictionary(dict(r["comps"]))
        rx = sh.sim.CreateHetCatReaction(r["name"], "", d, r["base"], "Vapor", "Partial Pressure", "atm",
                                         "kmol/[kg.h]", r["num"], r["den"])
        rx.Tmin, rx.Tmax = 0.0, 2000.0
        sh.reaction(rx)
    sh.mode("Isothermic")
    o = sh.o
    o.Volume, o.Headspace, o.CatalystAmount = spike.VOLUME, spike.VOLUME, spike.LOADING * spike.VOLUME
    o.IsothermalTemperature = spike.T_FEED
    o.Tolerance = 1e-12
    sh.set_feed(spike.T_FEED, spike.P_FEED, flows=spike.FEED)
    sh.solve()
    dw = sh.outlet()
    reacted = max(abs(a - b) for a, b in zip(dw["flows"], spike.FEED))
    # The exact answer of DWSIM's equations: F - F0 = W·V·Σ ν r(y·P), r as
    # typed in kmol/(kg·h) on partial pressures in atm.
    model = spike.Model()
    w = spike.LOADING * spike.VOLUME
    feed = [mp.mpf(x) for x in spike.FEED]
    p_atm = mp.mpf(spike.P_FEED) / 101325
    t = mp.mpf(spike.T_FEED)

    def rates(n):
        tot = sum(n)
        out = []
        for r in model.rx:
            env = {"exp": mp.exp, "T": t}
            env.update({k: n[i] / tot * p_atm for k, i in r["names"].items()})
            out.append(eval(r["num"], {"__builtins__": {}}, env) / eval(r["den"], {"__builtins__": {}}, env))
        return out

    # Two extents span what the three reactions can do (C, H and O are
    # conserved in five species): x1 the methane reformed, x2 the CO2 made.
    names = REFORMING["names"]
    i = {c: names.index(c) for c in names}

    def flows(x1, x2):
        n = list(feed)
        n[i["Methane"]] -= x1
        n[i["Water"]] -= x1 + x2
        n[i["Hydrogen"]] += 3 * x1 + x2
        n[i["Carbon monoxide"]] += x1 - x2
        n[i["Carbon dioxide"]] += x2
        return n

    def residual(x1, x2, scale=1):
        n = flows(x1, x2)
        made = [mp.mpf(0)] * len(n)
        for r, rate in zip(model.rx, rates(n)):
            for k, v in enumerate(r["nu"]):
                if v:
                    made[k] += scale * w * rate * 1000 / 3600 * v / abs(r["nub"])  # per mole of the base reactant
        return [-made[i["Methane"]] - x1, made[i["Carbon dioxide"]] - x2]

    # The rates are steep near equilibrium, so the catalyst is brought in
    # by steps from a millionth, each solve starting from the last.
    x = [mp.mpf(0), mp.mpf(0)]
    for k in range(241):
        scale = mp.mpf(10) ** (mp.mpf(k) / 40 - 6)
        x = list(mp.findroot(lambda a, b: residual(a, b, scale), x, tol=mp.mpf("1e-40") if k == 240 else
                             mp.mpf("1e-24"), maxsteps=200))
    sol = flows(*x)
    assert max(abs(v) for v in residual(*x)) < mp.mpf("1e-25"), residual(*x)
    exact = [float(x) for x in sol]
    print(f"  CSTR: DWSIM moves the feed by {reacted:.1e} mol/s; the exact answer {exact}", flush=True)
    return {"dwsim": {"outlet": dw, "reacted": reacted}, "names": REFORMING["names"], "cas": REFORMING["cas"],
            "feed": spike.FEED,
            "t": spike.T_FEED, "p": spike.P_FEED, "volume": spike.VOLUME, "loading": spike.LOADING, "exact": exact}


def ideal_gas(auto):
    """Each reforming compound's ideal-gas H(T) − H(298.15 K) at the PFR's
    T, DWSIM's (RET_Hid) and an exact quadrature of reflux's Cp, with
    DWSIM's Cp equation and coefficients: DWSIM's copy of ChemSep gives
    water eq 16 where reflux's file has eq 107."""
    from System import Array, Double

    sh = spike.Sheet(auto, reactor=False)
    names, t = REFORMING["names"], spike.T_FEED
    out = {"names": names, "cas": REFORMING["cas"], "t": t, "dwsim": [], "reflux": [], "dwsim_cp": []}
    for i, (name, cas) in enumerate(zip(names, REFORMING["cas"])):
        c = sh.sim.SelectedCompounds[name]
        x = Array[Double]([1.0 if j == i else 0.0 for j in range(len(names))])
        sh.pp.CurrentMaterialStream = sh.feed.GetAsObject()
        out["dwsim"].append(sh.pp.RET_Hid(298.15, t, x) * c.Molar_Weight)
        out["reflux"].append(float(cp_integrals(component(cas)["correlations"]["ideal_gas_cp"], t)[0]))
        out["dwsim_cp"].append([c.Ideal_Gas_Heat_Capacity_Const_A, c.Ideal_Gas_Heat_Capacity_Const_B,
                                c.Ideal_Gas_Heat_Capacity_Const_C, c.Ideal_Gas_Heat_Capacity_Const_D,
                                c.Ideal_Gas_Heat_Capacity_Const_E])
    print(f"  ideal gas at {t} K, DWSIM less reflux: "
          f"{[round(a - b, 3) for a, b in zip(out['dwsim'], out['reflux'])]} J/mol", flush=True)
    return out


def viscosity(fluid, t, p, n):
    """DWSIM's gas viscosity (M9 step 0): the mole-fraction mean of ChemSep's
    eq 102, corrected by Jossi, Stiel and Thodos with the mixture's Tc, Pc,
    Vc and M as mole-fraction means and its molar volume from PR."""
    tot = sum(n)
    y = [v / tot for v in n]
    mean = sum(a * spike.dippr102(d["correlations"]["vapor_viscosity"]["coeffs"], t) for a, d in zip(y, fluid.docs))
    crit = [sum(a * d["constants"][k] for a, d in zip(y, fluid.docs)) for k in ("tc", "pc", "vc")]
    mm = sum(a * d["constants"]["mw"] for a, d in zip(y, fluid.docs)) * 1e3
    v = fluid.z(t, p, n) * R_EXACT * t / p
    return spike.jossi_stiel_thodos(mean, v, crit[0], crit[1], crit[2], mm)


def pfr_exact_ergun(fluid, rtol=1e-12):
    """DWSIM's equations without its sections, with Ergun's drop applied
    continuously: dP/dV = -(Ergun's drop per length)/A, A = V/L, on DWSIM's
    viscosity method and PR's density, the rates at p_i = y_i·P."""
    model = spike.Model()
    area = spike.VOLUME / spike.LENGTH
    eps, dp = spike.VOID, spike.DP
    mws = fluid.mws

    def f(_v, y):
        n, p = list(y[:-1]), y[-1]
        tot = sum(n)
        dn = model.rhs(n, spike.T_FEED, p / 101325, tot)
        q = fluid.z(spike.T_FEED, p, n) * tot * R_EXACT * spike.T_FEED / p
        rho = sum(a * m for a, m in zip(n, mws)) / q
        u = q / area
        mu = viscosity(fluid, spike.T_FEED, p, n)
        drop = 150 * mu * (1 - eps) ** 2 * u / (eps ** 3 * dp ** 2) + 1.75 * rho * (1 - eps) * u ** 2 / (eps ** 3 * dp)
        return dn + [-drop / area]

    n0 = [x if x > 0.0 else 1e-10 for x in spike.FEED]
    sol = solve_ivp(f, (0.0, spike.VOLUME), n0 + [spike.P_FEED], method="Radau", rtol=rtol, atol=1e-16)
    assert sol.success, sol.message
    return list(sol.y[:-1, -1]), float(sol.y[-1, -1])


def pfr(auto):
    sh = spike.Sheet(auto)
    o = sh.o
    out = {"feed": spike.FEED, "t": spike.T_FEED, "p": spike.P_FEED, "names": spike.COMPS, "cas": spike.CAS,
           "volume": spike.VOLUME, "length": spike.LENGTH, "loading": spike.LOADING, "void": spike.VOID,
           "particle": spike.DP, "dwsim": []}
    for dv in [0.05, 0.001]:
        for ergun in [False, True]:
            o.dV = dv
            o.UseUserDefinedPressureDrop = not ergun
            o.UserDefinedPressureDrop = 0.0
            dt, _ = sh.solve()
            r = sh.outlet()
            out["dwsim"].append({"dv": dv, "ergun": ergun, "flows": r["flows"], "x": r["x"], "drop": r["dp"],
                                 "duty": r["q"] * 1000.0, "seconds": round(dt, 1)})
            print(f"  PFR, DWSIM at dV {dv}, Ergun {ergun}: X {r['x']:.8f}, drop {r['dp']:.4f} Pa, {dt:.0f} s",
                  flush=True)
    flows, nfev = spike.continuous(spike.Model())
    out["exact_isobaric"] = {"flows": flows, "x": (spike.FEED[0] - flows[0]) / spike.FEED[0]}
    fluid = Fluid(spike.CAS, reflux_kij(spike.CAS))
    flows, p_out = pfr_exact_ergun(fluid)
    out["exact_ergun"] = {"flows": flows, "x": (spike.FEED[0] - flows[0]) / spike.FEED[0], "drop": spike.P_FEED - p_out}
    out["kij"] = fluid.kij
    print(f"  PFR, exact: isobaric X {out['exact_isobaric']['x']:.10f}; with Ergun X {out['exact_ergun']['x']:.10f}, "
          f"drop {out['exact_ergun']['drop']:.6f} Pa", flush=True)
    return out


def main():
    previous = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    auto = load()
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m9_dwsim.py",
            "python": platform.python_version(),
            "dwsim": str(auto.GetVersion()),
        },
        "note": "DWSIM 9's reactors on cases built from numbers (see the generator's docstring). Units: K, Pa, "
                "mol/s, W into the material, J/mol; flows in the order of names. dwsim parts are DWSIM's own; "
                "thermo parts are thermo 0.6.1's PR (exact constants, reflux's kij from pr.ipd) or PR78 with "
                "DWSIM's rounded constants and R; bar and atm are equilibria with fugacities over 1 bar and over "
                "1 atm, K at 1 bar from the component files' dhf, s0 and Cp; exact answers are DWSIM's equations "
                "solved without its approximations.",
    }
    for part, fn in [("conversion", conversion), ("equilibrium", equilibrium_part), ("gibbs", gibbs),
                     ("cstr", cstr), ("pfr", pfr), ("ideal_gas", ideal_gas)]:
        data[part] = fn(auto) if part in PARTS else previous.get(part, {})
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}, DWSIM {data['source']['dwsim']}")


if __name__ == "__main__":
    main()
