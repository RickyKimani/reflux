# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""M9 step 0: DWSIM 9's plug-flow reactor, headless.

Run from the repository root, on Windows with DWSIM 9 installed:

    uv run --python 3.13 spike/dwsim_pfr/run.py [sections]

with sections among 1 to 6: all but 6 by default, since 6, the owner's 50-point
study, takes about a quarter of an hour. Section 3 takes about half an hour.

It answers step 0's questions (docs/m9-plan.md) and prints every number
docs/spike-dwsim-pfr.md quotes. The case is the owner's methane steam
reforming PFR, built here from its numbers: Xu and Froment's three rates as
DWSIM's text boxes have them, on an isothermal PFR of 1 m³ with Ergun's
pressure drop. If REFLUX_TUTORIAL names the owner's own DWSIM file, section 1
loads it too and checks that this script builds the same case. The file isn't
in the repository.

1. The case at the file's dV = 0.05: its outlet, and the file's.
2. How DWSIM reads its inputs: the rate's units and the catalyst loading,
   its partial pressures, Ergun's ΔP section by section, the void fraction,
   and its gas viscosity. Each is checked by reproducing DWSIM's numbers.
3. dV from 0.05 down, with each of DWSIM's four ODE solvers, against the
   exact answer of DWSIM's own model.
4. Failures: are they silent?
5. Where the time goes.
6. The owner's study: the feed from 700 to 1000 °C in 50 points.
"""

import math
import os
import pathlib
import sys
import time
import tomllib

from scipy.integrate import solve_ivp

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader

SECTIONS = set(sys.argv[1:]) or {"1", "2", "3", "4", "5"}
TUTORIAL = os.environ.get("REFLUX_TUTORIAL")

# The owner's case, as its DWSIM file has it. Compounds in the file's order.
COMPS = ["Methane", "Hydrogen", "Water", "Carbon dioxide", "Carbon monoxide"]
SHORT = ["CH4", "H2", "H2O", "CO2", "CO"]
CAS = ["74-82-8", "1333-74-0", "7732-18-5", "124-38-9", "630-08-0"]
FEED = [10.0 / 3.6, 0.1 / 3.6, 10.0 / 3.6, 0.0, 0.0]  # mol/s: 10, 0.1 and 10 kmol/h
T_FEED, P_FEED = 1000.0, 101325.0
VOLUME, LENGTH, DV = 1.0, 1.0, 0.05  # m³, m, and the section as a fraction of the volume
LOADING, VOID, DP = 0.386, 0.4, 0.002  # kg of catalyst per m³ of reactor, the bed's void fraction, m
# Each reaction's compounds in the file's order, which is what names them:
# DWSIM calls a reaction's reactants R1, R2, ... and its products P1, P2, ...
# in the order of its compound list.
REACTIONS = [
    {"name": "Overall Reaction", "base": "Methane",
     "comps": [("Methane", -1.0), ("Hydrogen", 4.0), ("Water", -2.0), ("Carbon dioxide", 1.0)],
     "num": "1.02E+15*exp(-243900/8.314/T)/P1^3.5*(R1*R2^2-P1^4*P2/(exp(34.218-31266/T)))",
     "den": "(1+6.65E-4*exp(38280/8.314/T)*R1+1.77E+5*exp(-88680/8.314/T)*R2/P1"
            "+6.12E-9*exp(82900/8.314/T)*P1)^2"},
    {"name": "Steam Reforming", "base": "Methane",
     "comps": [("Methane", -1.0), ("Hydrogen", 3.0), ("Water", -1.0), ("Carbon monoxide", 1.0)],
     "num": "4.22E+15*exp(-240100/8.314/T)/P1^2.5*(R1*R2-P1^3*P2/(exp(30.42-27106/T)))",
     "den": "(1+8.23E-5*exp(70650/8.314/T)*P2+6.65E-4*exp(38280/8.314/T)*R1"
            "+1.77E+5*exp(-88680/8.314/T)*R2/P1+6.12E-9*exp(82900/8.314/T)*P1)^2"},
    {"name": "Water Gas Shift", "base": "Carbon monoxide",
     "comps": [("Hydrogen", 1.0), ("Water", -1.0), ("Carbon dioxide", 1.0), ("Carbon monoxide", -1.0)],
     "num": "1.96E+6*exp(-67130/8.314/T)/P1*(R1*R2-P1*P2/(exp(-3.798+4160/T)))",
     "den": "(1+1.77E+5*exp(-88680/8.314/T)*R1/P1+6.12E-9*exp(82900/8.314/T)*P1"
            "+8.23E-5*exp(70650/8.314/T)*R2)^2"},
]
SOLVERS = ["implicit Runge-Kutta 5 (Radau)", "explicit Runge-Kutta 4(5)", "Adams-Moulton", "Gear's BDF"]


def load():
    gd.load_dwsim()
    import clr

    clr.AddReference("DWSIM.MathOps.DotNumerics")
    clr.AddReference("Ciloci.Flee")
    from DWSIM.Automation import Automation3

    return Automation3()


def enum(obj, prop, name):
    setattr(obj, prop, getattr(getattr(obj, prop).__class__, name))


class Sheet:
    """The owner's case in a new DWSIM flowsheet: feed, PFR, product."""

    def __init__(self, auto, reactor=True):
        from System import Array, Double
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType

        self.auto, self.OT = auto, ObjectType
        self.vec = lambda xs: Array[Double](list(xs))
        self.sim = auto.CreateFlowsheet()
        for c in COMPS:
            compound = self.sim.AvailableCompounds[c]
            assert compound.OriginalDB == "ChemSep", (c, compound.OriginalDB)
            self.sim.SelectedCompounds.Add(c, compound)
        # DWSIM's default flash tolerances, as the owner's file has them: the
        # gas is one phase, and only the stability test runs.
        self.pp = self.sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)").__implementation__
        self.x = 0
        self.feed = self.add("MaterialStream", "feed")
        self.set_feed(T_FEED)
        self.pp.CurrentMaterialStream = self.feed.GetAsObject()
        if not reactor:
            return
        self.reactions = [self.reaction(r) for r in REACTIONS]
        self.pfr = self.add("RCT_PFR", "PFR")
        self.product = self.add("MaterialStream", "product")
        self.q = self.add("EnergyStream", "Q")
        self.sim.ConnectObjects(self.feed.GraphicObject, self.pfr.GraphicObject, -1, 0)
        self.sim.ConnectObjects(self.pfr.GraphicObject, self.product.GraphicObject, 0, -1)
        self.sim.ConnectObjects(self.q.GraphicObject, self.pfr.GraphicObject, -1, 1)
        o = self.pfr.GetAsObject()
        enum(o, "ReactorOperationMode", "Isothermic")
        enum(o, "ReactorSizingType", "Length")
        o.ReactionSetID = "DefaultSet"
        o.Volume, o.Length, o.dV = VOLUME, LENGTH, DV
        o.CatalystLoading, o.CatalystVoidFraction, o.CatalystParticleDiameter = LOADING, VOID, DP
        o.InternalSolver = 0

    def add(self, kind, name):
        self.x += 100
        return self.sim.AddObject(getattr(self.OT, kind), self.x, 100, name)

    def set_feed(self, t, flows=FEED, p=P_FEED):
        o = self.feed.GetAsObject()
        o.SetTemperature(t)
        o.SetPressure(p)
        o.SetMolarFlow(sum(flows))
        o.SetOverallComposition(self.vec(f / sum(flows) for f in flows))

    def reaction(self, r):
        from System import Double, String
        from System.Collections.Generic import Dictionary

        d = Dictionary[String, Double]()
        for c, nu in r["comps"]:
            d[c] = nu
        rx = self.sim.CreateHetCatReaction(r["name"], "", d, r["base"], "Vapor", "Partial Pressure", "atm",
                                           "kmol/[kg.h]", r["num"], r["den"])
        self.default_t_range = (rx.Tmin, rx.Tmax)
        rx.Tmin, rx.Tmax = 0.0, 2000.0  # the file's
        self.sim.AddReaction(rx)
        if not self.sim.ReactionSets.ContainsKey("DefaultSet"):
            raise RuntimeError("no default reaction set")
        self.sim.AddReactionToSet(rx.ID, "DefaultSet", True, 0)
        return rx

    @property
    def o(self):
        return self.pfr.GetAsObject()

    def solve(self, check=True):
        """Solve; with check, refuse any failure, which DWSIM's solver doesn't report."""
        t0 = time.perf_counter()
        errs = self.auto.CalculateFlowsheet2(self.sim)
        dt = time.perf_counter() - t0
        errs = [str(e) for e in errs] if errs else []
        if check:
            bad = errs + self.unit_errors()
            if bad:
                raise RuntimeError(f"DWSIM failed: {bad}")
        return dt, errs

    def unit_errors(self):
        o = self.o
        return [] if o.Calculated and not (o.ErrorMessage or "") else [
            f"PFR: calculated {o.Calculated}, {o.ErrorMessage!r}"]

    def outlet(self):
        m = self.product.GetAsObject()
        ph = m.Phases[0].Properties
        flows = [m.Phases[0].Compounds[c].MolarFlow for c in COMPS]
        if ph.temperature is None or ph.pressure is None or any(f is None for f in flows):
            return None
        return {"t": ph.temperature, "p": ph.pressure, "flows": flows,
                "x": (FEED[0] - flows[0]) / FEED[0], "dp": self.o.DeltaP, "q": self.o.DeltaQ}

    def profile(self):
        out = []
        for row in self.o.Profile:
            by = {it.Compound: (it.MolarFlow, it.MolarConcentration) for it in row.Item4}
            out.append({"x": row.Item1, "t": row.Item2, "p": row.Item3, "flows": [by[c][0] for c in COMPS],
                        "conc": [by[c][1] for c in COMPS]})
        return out


class Probe:
    """DWSIM's own properties at a state: a stream it calculates alone, on a
    flowsheet of its own with the same package."""

    def __init__(self, auto):
        self.sh = Sheet(auto, reactor=False)
        self.ms = self.sh.feed.GetAsObject()
        self.ms.SetPropertyPackage(self.sh.pp)

    def at(self, t, p, flows):
        self.sh.set_feed(t, flows, p)
        self.sh.pp.CurrentMaterialStream = self.ms
        self.ms.Calculate(True, True)
        v, m = self.ms.Phases[2].Properties, self.ms.Phases[0].Properties
        return {"z": v.compressibilityFactor, "qv": v.volumetric_flow, "rho": m.density, "mu_v": v.viscosity,
                "mu": m.viscosity, "xv": v.massfraction, "q": m.volumetric_flow, "mw": m.molecularWeight}


# DWSIM's model in Python, from the same strings.

def variables(r):
    """DWSIM's names for a reaction's compounds: R1, R2, ... and P1, P2, ..."""
    names, ir, ip = {}, 0, 0
    for c, nu in r["comps"]:
        if nu < 0:
            ir += 1
            names[f"R{ir}"] = COMPS.index(c)
        elif nu > 0:
            ip += 1
            names[f"P{ip}"] = COMPS.index(c)
    return names


class Model:
    """DWSIM's PFR equations: dN/dV = -W Σ_j r_j ν_ij/ν_bj, the rate in
    kmol/(kg·h) converted as DWSIM converts it, at partial pressures in atm."""

    def __init__(self, loading=LOADING):
        self.loading = loading
        self.rx = []
        for r in REACTIONS:
            nu = [0.0] * len(COMPS)
            for c, v in r["comps"]:
                nu[COMPS.index(c)] = v
            # ^ is Flee's power; Python's ** binds the same way in these strings.
            self.rx.append({"names": variables(r), "nu": nu, "nub": dict(r["comps"])[r["base"]],
                            "num": compile(r["num"].replace("^", "**"), r["name"], "eval"),
                            "den": compile(r["den"].replace("^", "**"), r["name"], "eval")})
        self.calls = 0

    def rates(self, p_atm, t):
        out = []
        for r in self.rx:
            env = {"exp": math.exp, "T": t}
            env.update({k: p_atm[i] for k, i in r["names"].items()})
            out.append(eval(r["num"], {"__builtins__": {}}, env) / eval(r["den"], {"__builtins__": {}}, env))
        return out

    def rhs(self, n, t, conv, qf):
        """As ODEFunc: C = N/Qf, p = C·conv with conv = Z·8.314·T in atm."""
        self.calls += 1
        p = [0.0 if ni < 0.0 else ni / qf * conv for ni in n]
        ri = [0.0] * len(COMPS)
        for r, rate in zip(self.rx, self.rates(p, t)):
            rsi = rate / 60 / 60 * 1000  # kmol/(kg·h) to mol/(kg·s), as DWSIM's converter does
            for i, nu in enumerate(r["nu"]):
                if nu != 0.0:
                    ri[i] += rsi * nu / r["nub"] * self.loading
        return [-x for x in ri]


def sections(dv):
    """DWSIM's sections, as Calculate_Internal steps through them."""
    out, cur, d0 = [], 0.0, dv
    while True:
        out.append(d0)
        cur += d0 * VOLUME
        if cur + d0 * VOLUME > VOLUME:
            d0 = (VOLUME - cur) / VOLUME
        if cur >= VOLUME:
            return out


def continuous(model, t=T_FEED, p=P_FEED, factor=1.0, rtol=1e-12):
    """DWSIM's equations without its sections, isobaric: p_i = y_i·P·factor,
    integrated to rtol by SciPy's Radau."""
    def f(_v, n):
        tot = sum(n)
        return model.rhs(list(n), t, p / 101325 * factor, tot)

    n0 = [x if x > 0.0 else 1e-10 for x in FEED]  # DWSIM starts a zero flow at 1e-10 mol/s
    sol = solve_ivp(f, (0.0, VOLUME), n0, method="Radau", rtol=rtol, atol=1e-16)
    assert sol.success, sol.message
    return list(sol.y[:, -1]), sol.nfev


def reproduce(model, starts, dv, integrator="exact"):
    """DWSIM's sections, each integrated with its start's Qf and conv frozen
    as ODEFunc has them. starts gives (conv, qf) for each section, from
    DWSIM's own state there. integrator "exact" is SciPy's Radau at 1e-12;
    "dwsim" is DWSIM's own (DotNumerics) at its default tolerances."""
    n = [x if x > 0.0 else 1e-10 for x in FEED]
    out = [n]
    for k, d in enumerate(sections(dv)):
        conv, qf = starts[k]
        h = d * VOLUME
        n0 = [x if x > 0.0 else 1e-10 for x in n]
        if integrator == "exact":
            sol = solve_ivp(lambda _v, y: model.rhs(list(y), T_FEED, conv, qf), (0.0, h), n0, method="Radau",
                            rtol=1e-12, atol=1e-16)
            assert sol.success, sol.message
            n = list(sol.y[:, -1])
        else:
            n = dotnumerics(model, n0, h, conv, qf)
        out.append(n)
    return out


def dotnumerics(model, n0, h, conv, qf):
    from System import Array, Double
    from DotNumerics.ODE import OdeFunction, OdeImplicitRungeKutta5, OdeSolution

    last = {}

    def f(_x, y):
        return Array[Double](model.rhs(list(y), T_FEED, conv, qf))

    def keep(_x, y):
        last["y"] = list(y)

    s = OdeImplicitRungeKutta5()
    s.InitializeODEs(OdeFunction(f), len(n0), 0.0, Array[Double](n0))
    s.Solve(Array[Double](n0), 0.0, 0.01 * h, h, OdeSolution(keep))
    return last["y"]


def section_starts(probe, prof, dv):
    """Each section's Qf and conv as ODEFunc sees them: from the state DWSIM
    flashed before it, at the pressure the profile records there."""
    out = []
    assert len(prof) == len(sections(dv)) + 1, (len(prof), len(sections(dv)))
    for k in range(len(sections(dv))):
        row = prof[k]
        st = probe.at(row["t"], row["p"], row["flows"] if k else FEED)
        out.append((st["z"] * 8.314 * row["t"] / 101325, st["qv"]))
    return out


def rel(a, b):
    return max(abs(x - y) / max(abs(y), 1e-300) for x, y in zip(a, b))


def show(label, o):
    f = "  ".join(f"{s} {v:.10f}" for s, v in zip(SHORT, o["flows"]))
    print(f"  {label}: X(CH4) {o['x']:.8f}, ΔP {o['dp']:.4f} Pa, Q {o['q']:.6f} kW; mol/s {f}")


def section_1(auto):
    print("\n1. The owner's case, at dV = 0.05")
    sh = Sheet(auto)
    print(f"  CreateHetCatReaction's default T range: {sh.default_t_range} K (set to the file's 0 and 2000)")
    from DWSIM.GlobalSettings import Settings
    print(f"  DWSIM's inspector enabled: {Settings.InspectorEnabled}")
    dt, _ = sh.solve()
    out = sh.outlet()
    show(f"built here, solved in {dt:.2f} s", out)
    print(f"  outlet T {out['t']:.6f} K, P {out['p']:.4f} Pa; {len(sh.profile())} profile points")
    if TUTORIAL:
        fs = auto.LoadFlowsheet(TUTORIAL)
        pfr = next(o for o in fs.SimulationObjects.Values if o.GetType().Name == "Reactor_PFR").__implementation__
        for r in REACTIONS:
            rx = next(x for x in fs.Reactions.Values if x.Name == r["name"])
            same = (rx.RateEquationNumerator == r["num"] and rx.RateEquationDenominator == r["den"]
                    and [(c, rx.Components[c].StoichCoeff) for c in rx.Components.Keys] == r["comps"]
                    and rx.BaseReactant == r["base"] and (rx.Tmin, rx.Tmax) == (0.0, 2000.0))
            print(f"  the file's {r['name']!r}: the same as this script's: {same}")
        print(f"  the file's PFR: V {pfr.Volume}, L {pfr.Length}, dV {pfr.dV}, loading {pfr.CatalystLoading}, "
              f"void {pfr.CatalystVoidFraction}, dp {pfr.CatalystParticleDiameter}, solver {pfr.InternalSolver}, "
              f"mode {pfr.ReactorOperationMode}")
        t0 = time.perf_counter()
        errs = auto.CalculateFlowsheet2(fs)
        dt = time.perf_counter() - t0
        prod = fs.SimulationObjects[pfr.GraphicObject.OutputConnectors[0].AttachedConnector.AttachedTo.Name]
        prod = prod.__implementation__
        flows = [prod.Phases[0].Compounds[c].MolarFlow for c in COMPS]
        print(f"  the owner's file solved in {dt:.2f} s, errors {list(errs) if errs else []}, PFR calculated "
              f"{pfr.Calculated}; its outlet's flows differ from this script's by {rel(flows, out['flows']):.2e} "
              f"(relative), its ΔP by {abs(pfr.DeltaP - out['dp']):.2e} Pa")
    return sh


def section_2(auto, sh):
    print("\n2. How DWSIM reads its inputs, from the dV = 0.05 solve")
    probe = Probe(auto)
    prof = sh.profile()
    out = sh.outlet()
    # The partial pressures: Z·8.314·T·C_i, C_i = N_i/Qf, so y_i·P times 8.314
    # over the R in DWSIM's density.
    st = probe.at(T_FEED, P_FEED, FEED)
    ratio = st["z"] * 8.314 * T_FEED * sum(FEED) / (st["qv"] * P_FEED)
    print(f"  the feed: Z {st['z']:.10f}, Qf {st['qv']:.10f} m³/s, so Σp_i/P = Z·8.314·T·N/(Qf·P) = {ratio:.12f} "
          f"(8.314/R = {8.314 / 8.31446261815324:.12f})")
    # The rates and the loading: DWSIM's sections reproduced with its own
    # Qf and Z at each section's start.
    model = Model()
    starts = section_starts(probe, prof, DV)
    for integ in ["dwsim", "exact"]:
        model.calls = 0
        rep = reproduce(model, starts, DV, integ)
        worst = max(rel(rep[k], prof[k]["flows"]) for k in range(1, len(rep)))
        print(f"  sections reproduced with {'DWSIM' if integ == 'dwsim' else 'an exact'} integrator: every "
              f"profile point's flows within {worst:.2e} (relative), the outlet's within "
              f"{rel(rep[-1], out['flows']):.2e}; {model.calls} rate evaluations")
    for w in [LOADING * 2]:
        rep = reproduce(Model(w), starts, DV, "dwsim")
        print(f"  with the loading doubled in the reproduction, the outlet moves by {rel(rep[-1], out['flows']):.2e}")
    # Ergun's ΔP, section by section, from DWSIM's state after each section.
    d = math.sqrt(4 * VOLUME / LENGTH / math.pi)
    area = math.pi * d * d / 4
    drops = [prof[k + 1]["p"] - prof[k + 2]["p"] for k in range(len(prof) - 2)] + [prof[-1]["p"] - out["p"]]
    worst, total = 0.0, 0.0
    terms = None
    for k, (dl, drop) in enumerate(zip(sections(DV), drops)):
        row = prof[k + 1]
        s = probe.at(row["t"], row["p"], row["flows"])
        u = s["qv"] / area
        mu = s["mu_v"] * s["xv"]
        lk = dl * LENGTH
        visc = 150 * mu * lk / DP ** 2 * (1 - VOID) ** 2 / VOID ** 3 * u
        inert = 1.75 * lk * s["rho"] / DP * (1 - VOID) / VOID ** 3 * u ** 2
        worst = max(worst, abs(visc + inert - drop) / drop)
        total += visc + inert
        if k == 0:
            terms = (visc, inert, u, s["rho"], mu)
    print(f"  Ergun per section, at the state after it with L = dV·length: every section's drop within "
          f"{worst:.2e} of DWSIM's; their sum {total:.6f} Pa, DWSIM's ΔP {out['dp']:.6f} Pa")
    print(f"  the first section: viscous term {terms[0]:.4f} Pa, inertial {terms[1]:.4f} Pa; u {terms[2]:.6f} m/s, "
          f"ρ {terms[3]:.6f} kg/m³, μ {terms[4]:.6e} Pa·s")
    print(f"  the pressure each section's rates use: the profile's P at its start, P_0 for the first two "
          f"({prof[0]['p']:.4f}, {prof[1]['p']:.4f}, {prof[2]['p']:.4f} Pa, ...)")
    # The void fraction, with no pressure drop: it's in Ergun alone.
    o = sh.o
    o.UseUserDefinedPressureDrop, o.UserDefinedPressureDrop = True, 0.0
    res = []
    for e in [VOID, 0.3, 0.6]:
        o.CatalystVoidFraction = e
        sh.solve()
        res.append(sh.outlet())
    print(f"  no pressure drop, void fraction 0.4, 0.3, 0.6: outlets differ by {rel(res[1]['flows'], res[0]['flows']):.1e} "
          f"and {rel(res[2]['flows'], res[0]['flows']):.1e}")
    show("isobaric", res[0])
    show("with Ergun", out)
    print(f"  Ergun's 4.9 kPa moves the outlet's flows by {rel(out['flows'], res[0]['flows']):.2e} (relative), "
          f"X(CH4) by {out['x'] - res[0]['x']:+.3e}")
    o.CatalystVoidFraction, o.UseUserDefinedPressureDrop = VOID, False
    viscosity(probe, prof)


def dippr102(c, t):
    a, b, cc, d = c
    return a * t ** b / (1 + cc / t + d / (t * t))


def jossi_stiel_thodos(eta0, v, tc, pc, vc, mm):
    """Jossi, Stiel and Thodos' dense-gas correction as DWSIM documents it:
    [(η − η0)ξ + 1]^¼ = 1.023 + 0.23364ρr + 0.58533ρr² − 0.40758ρr³ +
    0.093324ρr⁴, η in µP, ξ = (Tc/(M³Pc⁴))^⅙ with Pc in atm and M in g/mol,
    ρr = Vc/V. At ρr = 0 it isn't zero: 1.023⁴ − 1 = 0.0952."""
    xi = (tc / (mm ** 3 * (pc / 101325) ** 4)) ** (1 / 6)
    rr = vc / v
    rhs = 1.023 + 0.23364 * rr + 0.58533 * rr ** 2 - 0.40758 * rr ** 3 + 0.093324 * rr ** 4
    return eta0 + (rhs ** 4 - 1) / xi * 1e-7


def viscosity(probe, prof):
    """DWSIM's gas viscosity against ChemSep's correlations, as reflux's
    component files have them, three mixing rules, and DWSIM's own method:
    the mole-fraction mean corrected by Jossi, Stiel and Thodos with the
    mixture's constants as mole-fraction means."""
    pure, mw, crit = [], [], []
    for cas in CAS:
        doc = tomllib.loads((ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8"))
        cor = doc["correlations"]["vapor_viscosity"]
        assert cor["eq"] == 102, cas
        pure.append((cor["coeffs"], cor["t_range"]))
        mw.append(doc["constants"]["mw"])
        crit.append((doc["constants"]["tc"], doc["constants"]["pc"], doc["constants"]["vc"],
                     doc["constants"]["mw"] * 1e3))
    print("  pure gases at 1000 K, 1 atm: DWSIM's viscosity against ChemSep's correlation (eq. 102)")
    for i, s in enumerate(SHORT):
        flows = [1e-12] * len(COMPS)
        flows[i] = 1.0
        st = probe.at(T_FEED, P_FEED, flows)
        mu0 = dippr102(pure[i][0], T_FEED)
        jst = jossi_stiel_thodos(mu0, st["mw"] / st["rho"] / 1000, *crit[i])
        print(f"    {s:4} DWSIM {st['mu_v']:.6e}, ChemSep {mu0:.6e} Pa·s (its range {pure[i][1][0]}-"
              f"{pure[i][1][1]} K); corrected as DWSIM does, {(jst - st['mu_v']) / st['mu_v']:+.1e} from DWSIM's")
    for label, row in [("the feed", {"t": T_FEED, "p": P_FEED, "flows": FEED}), ("mid-reactor", prof[len(prof) // 2]),
                       ("the outlet", prof[-1])]:
        st = probe.at(row["t"], row["p"], row["flows"])
        y = [f / sum(row["flows"]) for f in row["flows"]]
        mu = [dippr102(c, row["t"]) for c, _ in pure]
        linear = sum(a * b for a, b in zip(y, mu))
        hz = sum(a * b * math.sqrt(m) for a, b, m in zip(y, mu, mw)) / sum(a * math.sqrt(m) for a, m in zip(y, mw))
        wilke = 0.0
        for i in range(len(y)):
            den = sum(y[j] * (1 + math.sqrt(mu[i] / mu[j]) * (mw[j] / mw[i]) ** 0.25) ** 2
                      / math.sqrt(8 * (1 + mw[i] / mw[j])) for j in range(len(y)))
            wilke += y[i] * mu[i] / den
        mean = [sum(a * c[k] for a, c in zip(y, crit)) for k in range(4)]
        jst = jossi_stiel_thodos(linear, st["mw"] / st["rho"] / 1000, *mean)
        print(f"  {label}: DWSIM's gas viscosity {st['mu_v']:.6e} Pa·s; from ChemSep's, mole-fraction mean "
              f"{linear:.6e}, Herning-Zipperer {hz:.6e}, Wilke {wilke:.6e}; the mean corrected as DWSIM does "
              f"{jst:.6e}, {(jst - st['mu_v']) / st['mu_v']:+.1e} from DWSIM's")


def section_3(auto, sh):
    print("\n3. dV and DWSIM's solvers")
    probe = Probe(auto)
    model = Model()
    o = sh.o
    # The exact answer of DWSIM's model with no sections, isobaric, for p_i
    # = y_i·P·(8.314/R_density) as section 2 found.
    st = probe.at(T_FEED, P_FEED, FEED)
    factor = st["z"] * 8.314 * T_FEED * sum(FEED) / (st["qv"] * P_FEED)
    exact, nfev = continuous(model, factor=factor)
    loose, _ = continuous(model, factor=factor, rtol=1e-10)
    x_exact = (FEED[0] - exact[0]) / FEED[0]
    print(f"  the exact isobaric answer: X(CH4) {x_exact:.10f}, {nfev} rate evaluations; at rtol 1e-10 it moves "
          f"by {rel(loose, exact):.1e}")
    print("  " + "  ".join(f"{s} {v:.10f}" for s, v in zip(SHORT, exact)))
    grid = [0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    for solver in range(4):
        o.InternalSolver = solver
        print(f"  solver {solver}, {SOLVERS[solver]}:")
        for dv in grid if solver == 0 else [0.05, 0.01]:
            o.dV = dv
            o.UseUserDefinedPressureDrop = False
            dt, _ = sh.solve()
            ergun = sh.outlet()
            o.UseUserDefinedPressureDrop, o.UserDefinedPressureDrop = True, 0.0
            dt2, _ = sh.solve()
            iso = sh.outlet()
            line = (f"    dV {dv:<6} Ergun: X {ergun['x']:.8f}, ΔP {ergun['dp']:9.4f} Pa, {dt:6.2f} s | isobaric: "
                    f"X {iso['x']:.8f}, {dt2:6.2f} s, off the exact answer by {rel(iso['flows'], exact):.2e} "
                    f"(X by {iso['x'] - x_exact:+.2e})")
            if solver == 0 and dv >= 0.005:
                starts = section_starts(probe, sh.profile(), dv)
                rep = reproduce(model, starts, dv, "exact")
                line += f"; its sections, exactly: {rel(rep[-1], iso['flows']):.2e} from it"
            print(line)
    o.InternalSolver, o.dV, o.UseUserDefinedPressureDrop = 0, DV, False


def section_4(auto):
    print("\n4. Failures")

    def run(label, change):
        sh = Sheet(auto)
        change(sh)
        dt, errs = sh.solve(check=False)
        out = sh.outlet()
        what = "no outlet" if out is None else f"outlet X(CH4) {out['x']:.6f}"
        print(f"  {label}: {dt:.2f} s; CalculateFlowsheet2's errors {errs}; PFR calculated {sh.o.Calculated}, "
              f"its message {sh.o.ErrorMessage!r}; product calculated {sh.product.GetAsObject().Calculated}, "
              f"{what}")

    def bad_name(sh):
        sh.reactions[1].RateEquationNumerator = REACTIONS[1]["num"].replace("R2", "R3")

    def bad_syntax(sh):
        sh.reactions[1].RateEquationDenominator = REACTIONS[1]["den"] + ")"

    def below_t_range(sh):
        for rx in sh.reactions:
            rx.Tmax = 900.0

    def negative_p(sh):
        sh.o.CatalystParticleDiameter = 2e-5

    def no_point(sh):
        sh.reactions[2].RateEquationNumerator = "1E+3*" + REACTIONS[2]["num"]

    def huge_rate(sh):
        sh.reactions[2].RateEquationNumerator = "1.0E+300*" + REACTIONS[2]["num"]

    def zero_denominator(sh):
        sh.reactions[2].RateEquationDenominator = "0*" + REACTIONS[2]["den"]

    run("an unknown name in a rate (R3)", bad_name)
    run("a rate that doesn't parse", bad_syntax)
    run("a number written 1E+3, with no decimal point", no_point)
    run("every reaction's Tmax at 900 K, under the feed's 1000 K", below_t_range)
    run("a pressure drop larger than the pressure (particles of 0.02 mm)", negative_p)
    run("a rate scaled by 1.0E+300", huge_rate)
    run("a denominator of zero", zero_denominator)


def section_5(auto, sh):
    print("\n5. Where the time goes")
    times = [sh.solve()[0] for _ in range(3)]
    print(f"  the case re-solved: {', '.join(f'{t:.2f}' for t in times)} s")
    # Flee compiles each rate's numerator and denominator on every call of
    # ODEFunc: what that costs alone.
    import clr
    from System import Double, Math
    from Ciloci.Flee import ExpressionContext

    from System.Globalization import CultureInfo

    def odefunc_flee():
        """What one ODEFunc call does with Flee, reaction by reaction: a new
        context, System.Math imported, the variables, both expressions
        compiled and evaluated."""
        for r in REACTIONS:
            c = ExpressionContext()
            c.Imports.AddType(clr.GetClrType(Math))
            c.Options.ParseCulture = CultureInfo.InvariantCulture
            c.Variables.Clear()
            c.Variables.Add("T", 1000.0)
            for k in ["R1", "R2", "P1", "P2"]:
                c.Variables.Add(k, 0.3)
            c.CompileGeneric[Double](r["num"]).Evaluate()
            c.CompileGeneric[Double](r["den"]).Evaluate()

    n = 300
    odefunc_flee()
    t0 = time.perf_counter()
    for _ in range(n):
        odefunc_flee()
    per = (time.perf_counter() - t0) / n
    ctx = ExpressionContext()
    ctx.Imports.AddType(clr.GetClrType(Math))
    ctx.Variables.Add("T", 1000.0)
    for k in ["R1", "R2", "P1", "P2"]:
        ctx.Variables.Add(k, 0.3)
    t0 = time.perf_counter()
    for _ in range(n):
        for r in REACTIONS:
            ctx.CompileGeneric[Double](r["num"]).Evaluate()
            ctx.CompileGeneric[Double](r["den"]).Evaluate()
    reuse = (time.perf_counter() - t0) / n
    e = ctx.CompileGeneric[Double](REACTIONS[1]["num"])
    t0 = time.perf_counter()
    for _ in range(n * 10):
        e.Evaluate()
    ev = (time.perf_counter() - t0) / (n * 10)
    print(f"  Flee: one ODEFunc call's contexts, compiles and evaluations of the six expressions: {per * 1e3:.3f} ms "
          f"({reuse * 1e3:.3f} ms with one context kept); evaluating one compiled expression from Python: "
          f"{ev * 1e6:.2f} µs")
    probe = Probe(auto)
    prof = sh.profile()
    model = Model()
    model.calls = 0
    reproduce(model, section_starts(probe, prof, DV), DV, "dwsim")
    print(f"  DotNumerics' Radau at DWSIM's tolerances takes {model.calls} rate evaluations over the 20 sections, "
          f"and ODEFunc runs once more after each, so Flee's work alone is about "
          f"{(model.calls + len(sections(DV))) * per:.1f} s of a solve")
    t0 = time.perf_counter()
    for row in prof:
        probe.at(row["t"], row["p"], row["flows"])
    print(f"  a stream's flash and properties: {(time.perf_counter() - t0) / len(prof) * 1e3:.1f} ms")


def section_6(auto, sh):
    print("\n6. The owner's study: the feed from 700 to 1000 °C in 50 points, at dV = 0.05")
    ts = [973.15 + 300.0 * i / 49 for i in range(50)]
    t0 = time.perf_counter()
    rows = []
    for t in ts:
        sh.set_feed(t)
        dt, _ = sh.solve()
        rows.append((t, sh.outlet()["x"], dt))
    total = time.perf_counter() - t0
    for t, x, dt in rows[::7] + [rows[-1]]:
        print(f"  {t - 273.15:7.2f} °C: X(CH4) {x:.6f}, {dt:.2f} s")
    print(f"  50 solves in {total:.1f} s ({total / 60:.1f} min), {min(r[2] for r in rows):.2f} to "
          f"{max(r[2] for r in rows):.2f} s each")
    sh.set_feed(T_FEED)


def main():
    t0 = time.perf_counter()
    auto = load()
    print(f"{auto.GetVersion()}; loaded in {time.perf_counter() - t0:.1f} s")
    sh = section_1(auto)
    if "2" in SECTIONS:
        section_2(auto, sh)
    if "3" in SECTIONS:
        section_3(auto, sh)
    if "4" in SECTIONS:
        section_4(auto)
    if "5" in SECTIONS:
        section_5(auto, sh)
    if "6" in SECTIONS:
        section_6(auto, sh)


if __name__ == "__main__":
    main()
