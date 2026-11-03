# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""M4 step 0: DWSIM 9's flowsheets, headless.

Run from the repository root, on Windows with DWSIM 9 installed:

    uv run --python 3.13 spike/dwsim_flowsheet/run.py [sections]

with sections among 1 to 4 (all by default; section 2 takes a few minutes,
mostly Broyden failing).

It answers step 0's questions (docs/m4-plan.md) and prints every number
docs/spike-dwsim-flowsheet.md quotes:
- the reference flowsheet: feed, mixer, cooler, flash drum, splitter, a pump
  back to the mixer through a recycle block, and a valve to a second drum;
- how tightly DWSIM's recycle converges, with each acceleration method, and
  what the overall balances close to then;
- DWSIM's unit models, each against DWSIM's own flashes: pump, compressor,
  expander, mixer, valve, heater, and the drum's modes;
- the open-loop fallback (D12): the loop cut, fed DWSIM's own converged
  recycle stream, must hand it back.
"""

import math
import sys
import time

sys.path.insert(0, "tools")
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader

COMPS = ["Methane", "Ethane", "Propane", "N-butane", "N-pentane"]
Z = [0.70, 0.10, 0.08, 0.07, 0.05]
FEED = dict(t=300.0, p=30e5, f=100.0)  # K, Pa, mol/s
COOL_T = 230.0
PUMP_P = 31e5
VALVE_P = 5e5
RATIOS = [0.5, 0.5, 0.0]  # DWSIM's splitter wants one ratio per port, three in all
TOL_FIELDS = ["Temperatura", "Pressao", "VazaoMassica", "FracaoVapor", "Entalpia", "Entropia", "Composicao"]
TIGHT = 1e-8  # the scale on DWSIM's recycle tolerances past which nothing improves
SECTIONS = set(sys.argv[1:]) or {"1", "2", "3", "4"}


def load():
    gd.load_dwsim()
    from DWSIM.Automation import Automation3

    return Automation3()


class Sheet:
    """One DWSIM flowsheet with PR, its flash tolerances at 1e-8 as in M3."""

    def __init__(self, auto):
        from System import Array, Double
        from DWSIM.Interfaces.Enums import FlashCalculationType, FlashSetting
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType

        self.auto, self.OT, self.FCT = auto, ObjectType, FlashCalculationType
        self.vec = lambda xs: Array[Double](list(xs))
        self.sim = auto.CreateFlowsheet()
        for c in COMPS:
            self.sim.SelectedCompounds.Add(c, self.sim.AvailableCompounds[c])
        self.pp = self.sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)").__implementation__
        for key in (FlashSetting.PTFlash_External_Loop_Tolerance, FlashSetting.PTFlash_Internal_Loop_Tolerance,
                    FlashSetting.PHFlash_External_Loop_Tolerance, FlashSetting.PHFlash_Internal_Loop_Tolerance):
            self.pp.FlashSettings[key] = "1E-8"
        self.x = 0
        self.mw = [self.sim.SelectedCompounds[c].Molar_Weight for c in COMPS]  # g/mol

    def add(self, kind, name):
        self.x += 100
        return self.sim.AddObject(getattr(self.OT, kind), self.x, 100, name)

    def link(self, a, b, ia=-1, ib=-1):
        self.sim.ConnectObjects(a.GraphicObject, b.GraphicObject, ia, ib)

    def feed(self, name, t, p, f, z):
        s = self.add("MaterialStream", name)
        o = s.GetAsObject()
        o.SetTemperature(t)
        o.SetPressure(p)
        o.SetMolarFlow(f)
        o.SetOverallComposition(self.vec(z))
        if self.pp.CurrentMaterialStream is None:
            self.pp.CurrentMaterialStream = o  # the package reads its compounds from this
        return s

    def solve(self):
        t0 = time.perf_counter()
        errs = self.auto.CalculateFlowsheet2(self.sim)
        return time.perf_counter() - t0, [str(e) for e in errs] if errs else []

    def unit_errors(self, units):
        """DWSIM's flowsheet solver doesn't report a unit's failure: look in each."""
        out = []
        for u in units:
            o = u.GetAsObject()
            if not o.Calculated or (o.ErrorMessage or ""):
                out.append(f"{u.GraphicObject.Tag}: calculated {o.Calculated}, {o.ErrorMessage!r}")
        return out

    def state(self, stream):
        m = stream.GetAsObject()
        ph = m.Phases[0].Properties
        if ph.molarflow is None:
            return None
        vf = m.Phases[2].Properties.molarfraction
        return dict(t=ph.temperature, p=ph.pressure, f=ph.molarflow, h=ph.molar_enthalpy, s=ph.molar_entropy,
                    vf=vf if vf is not None else 0.0, v=(ph.molecularWeight / ph.density) * 1e-3,
                    flows=[m.Phases[0].Compounds[c].MolarFlow for c in COMPS],
                    z=[m.Phases[0].Compounds[c].MoleFraction for c in COMPS])

    def flash(self, kind, a, b, z):
        """DWSIM's own flash, through the package, on a mass basis like its API."""
        r = self.pp.CalculateEquilibrium(getattr(self.FCT, kind), a, b, self.vec(z), None, 0.0)
        mw = sum(zi * m for zi, m in zip(z, self.mw))
        return dict(t=r.CalculatedTemperature, p=r.CalculatedPressure, h=r.CalculatedEnthalpy * mw,
                    s=r.CalculatedEntropy * mw, vf=r.GetVaporPhaseMoleFraction())


def reference(auto, tol_scale=None, accel=None, max_iter=50, open_loop=None, legacy=True):
    """The reference flowsheet. tol_scale multiplies DWSIM's recycle tolerances
    (None leaves its defaults). open_loop, a stream state, cuts the loop and
    feeds that state in where the recycle's outlet was."""
    sh = Sheet(auto)
    feed = sh.feed("feed", FEED["t"], FEED["p"], FEED["f"], Z)
    o = {"feed": feed}
    for kind, name in [("NodeIn", "MIX"), ("MaterialStream", "s1"), ("Cooler", "COOL"), ("EnergyStream", "Qcool"),
                       ("MaterialStream", "s2"), ("Vessel", "DRUM"), ("MaterialStream", "vapour"),
                       ("MaterialStream", "liquid"), ("NodeOut", "SPLIT"), ("MaterialStream", "rec_a"),
                       ("MaterialStream", "liq_product"), ("Pump", "PUMP"), ("EnergyStream", "Qpump"),
                       ("MaterialStream", "rec_b"), ("Valve", "VALVE"), ("MaterialStream", "s3"),
                       ("Vessel", "DRUM2"), ("MaterialStream", "vapour2"), ("MaterialStream", "liquid2")]:
        o[name] = sh.add(kind, name)
    if open_loop is None:
        o["RECYCLE"] = sh.add("OT_Recycle", "RECYCLE")
        o["rec_c"] = sh.add("MaterialStream", "rec_c")
        rc = o["rec_c"].GetAsObject()  # DWSIM needs a first estimate for the torn stream
        rc.SetTemperature(COOL_T)
        rc.SetPressure(PUMP_P)
        rc.SetMolarFlow(20.0)
        rc.SetOverallComposition(sh.vec([0.2, 0.15, 0.25, 0.22, 0.18]))
    else:
        o["rec_c"] = sh.feed("rec_c", open_loop["t"], open_loop["p"], open_loop["f"], open_loop["z"])
    L = sh.link
    L(feed, o["MIX"], -1, 0)
    L(o["rec_c"], o["MIX"], -1, 1)
    L(o["MIX"], o["s1"])
    L(o["s1"], o["COOL"])
    L(o["COOL"], o["s2"])
    L(o["COOL"], o["Qcool"])
    L(o["s2"], o["DRUM"])
    L(o["DRUM"], o["vapour"], 0, -1)
    L(o["DRUM"], o["liquid"], 1, -1)
    L(o["liquid"], o["SPLIT"])
    L(o["SPLIT"], o["rec_a"], 0, -1)
    L(o["SPLIT"], o["liq_product"], 1, -1)
    L(o["rec_a"], o["PUMP"])
    L(o["Qpump"], o["PUMP"])
    L(o["PUMP"], o["rec_b"])
    if open_loop is None:
        L(o["rec_b"], o["RECYCLE"])
        L(o["RECYCLE"], o["rec_c"])
    L(o["liq_product"], o["VALVE"])
    L(o["VALVE"], o["s3"])
    L(o["s3"], o["DRUM2"])
    L(o["DRUM2"], o["vapour2"], 0, -1)
    L(o["DRUM2"], o["liquid2"], 1, -1)
    c = o["COOL"].GetAsObject()
    c.CalcMode = c.CalcMode.__class__.OutletTemperature
    c.OutletTemperature = COOL_T
    c.DeltaP = 0.0
    s = o["SPLIT"].GetAsObject()
    s.Ratios.Clear()
    for r in RATIOS:
        s.Ratios.Add(r)
    pu = o["PUMP"].GetAsObject()
    pu.CalcMode = pu.CalcMode.__class__.OutletPressure
    pu.Pout = PUMP_P
    va = o["VALVE"].GetAsObject()
    va.CalcMode = va.CalcMode.__class__.OutletPressure
    va.OutletPressure = VALVE_P
    if open_loop is None:
        r = o["RECYCLE"].GetAsObject()
        r.MaximumIterations = max_iter
        r.LegacyMode = legacy
        if accel is not None:
            r.AccelerationMethod = getattr(r.AccelerationMethod.__class__, accel)
        if tol_scale is not None:
            cp = r.ConvergenceParameters
            for f in TOL_FIELDS:
                setattr(cp, f, getattr(cp, f) * tol_scale)
    return sh, o


UNITS = ["MIX", "COOL", "DRUM", "SPLIT", "PUMP", "VALVE", "DRUM2"]
PRODUCTS = ["vapour", "vapour2", "liquid2"]


def balances(sh, o):
    """Overall component and energy balances, relative to the feed."""
    fd = sh.state(o["feed"])
    prods = [sh.state(o[n]) for n in PRODUCTS]
    comp = max(abs(sum(p["flows"][i] for p in prods) - fd["flows"][i]) for i in range(len(COMPS))) / fd["f"]
    q_cool = o["COOL"].GetAsObject().DeltaQ * 1e3  # kW to W, heat removed
    w_pump = o["PUMP"].GetAsObject().DeltaQ * 1e3
    h_in = fd["f"] * fd["h"] - q_cool + w_pump
    h_out = sum(p["f"] * p["h"] for p in prods)
    total = abs(sum(p["f"] for p in prods) - fd["f"]) / fd["f"]
    return comp, abs(h_out - h_in) / abs(fd["f"] * fd["h"]), q_cool, w_pump, total


def show(sh, o, names):
    for n in names:
        st = sh.state(o[n])
        if st is None:
            print(f"    {n:12} not calculated")
            continue
        print(f"    {n:12} T {st['t']:.6f} K  P {st['p'] / 1e5:.6f} bar  F {st['f']:.8f} mol/s  vf {st['vf']:.8f}")


def main():
    t_load = time.perf_counter()
    auto = load()
    print(f"{auto.GetVersion()}; loaded in {time.perf_counter() - t_load:.1f} s")

    # 1. The reference flowsheet at DWSIM's defaults.
    print("\n1. The reference flowsheet, at DWSIM's default recycle tolerances")
    sh, o = reference(auto)
    t_cold, errs = sh.solve()
    r = o["RECYCLE"].GetAsObject()
    print(f"  solved in {t_cold:.2f} s, {r.IterationsTaken} recycle iterations, converged {r.Converged}, "
          f"reported errors {errs}, unit errors {sh.unit_errors([o[u] for u in UNITS + ['RECYCLE']])}")
    show(sh, o, ["feed", "s2", "vapour", "liquid", "rec_c", "liq_product", "vapour2", "liquid2"])
    comp, energy, q, w, total = balances(sh, o)
    print(f"  overall balances: total flow {total:.3g}, worst component {comp:.3g}, energy {energy:.3g} "
          f"(of the feed's); cooler {q / 1e3:.4f} kW, pump {w:.4f} W")
    t0 = time.perf_counter()
    for i in range(5):
        o["feed"].GetAsObject().SetTemperature(FEED["t"] + 0.01 * (i + 1))
        sh.solve()
    print(f"  warm re-solves: {(time.perf_counter() - t0) / 5 * 1e3:.0f} ms each")

    # 2. How tight the recycle converges, by method.
    if "2" in SECTIONS:
        print("\n2. The recycle's tolerances scaled down: iterations, converged, overall balances, time")
        runs = [(a, True, sc) for a in ["None", "Wegstein", "Dominant_Eigenvalue"]
                for sc in [1.0, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10]]
        runs += [(a, False, sc) for a in ["None", "Wegstein", "Dominant_Eigenvalue"] for sc in [1.0, 1e-4, 1e-8]]
        runs += [("GlobalBroyden", lg, 1e-4) for lg in [True, False]]
        for accel, legacy, scale in runs:
            sh, o = reference(auto, tol_scale=scale, accel=accel, max_iter=500, legacy=legacy)
            dt, errs = sh.solve()
            r = o["RECYCLE"].GetAsObject()
            bad = sh.unit_errors([o[u] for u in UNITS])
            label = f"{accel}{'' if legacy else ', not legacy'}"
            if bad or sh.state(o["vapour"]) is None:
                print(f"  {label:32} x{scale:.0e}: failed: {errs or bad}")
                continue
            comp, energy, _, _, total = balances(sh, o)
            print(f"  {label:32} x{scale:.0e}: {r.IterationsTaken:4} iterations, converged {r.Converged}, "
                  f"total {total:.2e}, components {comp:.2e}, energy {energy:.2e}, {dt:.2f} s")

    # The reference for what follows: DWSIM's recycle tolerances scaled by
    # 1e-8, where its balances have stopped improving.
    sh, o = reference(auto, tol_scale=TIGHT, max_iter=500)
    dt, _ = sh.solve()
    r = o["RECYCLE"].GetAsObject()
    comp, energy, _, _, total = balances(sh, o)
    print(f"\nThe tight reference (tolerances x{TIGHT:.0e}): {r.IterationsTaken} iterations, {dt:.2f} s, "
          f"total {total:.2e}, components {comp:.2e}, energy {energy:.2e}; "
          f"unit errors {sh.unit_errors([o[u] for u in UNITS])}")
    show(sh, o, ["s2", "vapour", "rec_b", "vapour2", "liquid2"])
    t0 = time.perf_counter()
    for i in range(5):
        o["feed"].GetAsObject().SetTemperature(FEED["t"] + 0.01 * (i + 1))
        sh.solve()
    print(f"  warm re-solves: {(time.perf_counter() - t0) / 5 * 1e3:.0f} ms each")
    o["feed"].GetAsObject().SetTemperature(FEED["t"])
    sh.solve()

    # 3. DWSIM's unit models, against its own flashes.
    if "3" in SECTIONS:
        units(auto, sh, o)
    if "4" in SECTIONS:
        open_loop(auto, sh, o)


def units(auto, sh, o):
    print("\n3. Unit models, from the tight reference")
    a, b = sh.state(o["rec_a"]), sh.state(o["rec_b"])
    pu = o["PUMP"].GetAsObject()
    w = pu.DeltaQ * 1e3 / a["f"]  # J/mol
    eta = pu.Efficiency / 100.0
    incompressible = a["v"] * (b["p"] - a["p"]) / eta
    hs = sh.flash("PressureEntropy", b["p"], a["s"] / sum(zi * m for zi, m in zip(a["z"], sh.mw)), a["z"])["h"]
    print(f"  pump: work {w:.6f} J/mol; v·ΔP/η {incompressible:.6f} (v {a['v'] * 1e6:.4f} cm³/mol); "
          f"isentropic ΔH_s/η {(hs - a['h']) / eta:.6f}; outlet H − inlet H {b['h'] - a['h']:.6f}")
    s1, s2 = sh.state(o["liq_product"]), sh.state(o["s3"])
    print(f"  valve: H in {s1['h']:.6f}, out {s2['h']:.6f} J/mol (isenthalpic if equal)")
    # DWSIM's own PS and PH flashes, from a state's own S or H at its P: they
    # should come back to the state.
    for label, st in [("liquid (rec_a)", a), ("two-phase (feed)", sh.state(o["feed"]))]:
        mwi = sum(zi * m for zi, m in zip(st["z"], sh.mw))
        ps = sh.flash("PressureEntropy", st["p"], st["s"] / mwi, st["z"])
        ph = sh.flash("PressureEnthalpy", st["p"], st["h"] / mwi, st["z"])
        print(f"  round trips, {label}: PS lands {ps['t'] - st['t']:+.3e} K and {ps['h'] - st['h']:+.4f} J/mol "
              f"away; PH {ph['t'] - st['t']:+.3e} K")
    # Three liquid volumes for the pump's inlet: from DWSIM's reported density
    # (which its pump uses), from the phase's own Z, and the one DWSIM's H and
    # S imply. At constant T, dH − T dS = v dP exactly, so between two
    # compressed-liquid states 1 bar apart, ΔH − T·ΔS = v̄·ΔP. 31 and 32 bar
    # keep both states clear of the bubble point at 30 bar.
    liq = o["rec_a"].GetAsObject().Phases[3].Properties
    v_density = a["v"]
    v_z = liq.compressibilityFactor * 8.314462618 * a["t"] / a["p"]
    p1, p2 = a["p"] + 1e5, a["p"] + 2e5
    f1 = sh.flash("PressureTemperature", p1, a["t"], a["z"])
    f2 = sh.flash("PressureTemperature", p2, a["t"], a["z"])
    v_hs = ((f2["h"] - f1["h"]) - a["t"] * (f2["s"] - f1["s"])) / (p2 - p1)
    print(f"  pump inlet liquid volume: {v_density * 1e6:.3f} cm³/mol from DWSIM's density "
          f"({liq.density:.3f} kg/m³), {v_z * 1e6:.3f} from its Z ({liq.compressibilityFactor:.6f}), "
          f"{v_hs * 1e6:.3f} implied by its H and S (vapour fractions {f1['vf']}, {f2['vf']})")
    c1, c2 = sh.state(o["s1"]), sh.state(o["s2"])
    q = o["COOL"].GetAsObject().DeltaQ * 1e3 / c1["f"]
    print(f"  cooler: duty {q:.6f} J/mol; H in − H out {c1['h'] - c2['h']:.6f}")
    m_in = [sh.state(o["feed"]), sh.state(o["rec_c"])]
    print(f"  mixer: outlet P {c1['p']:.1f} Pa, inlets {[s['p'] for s in m_in]} (DWSIM's default is the minimum); "
          f"ΣF·H in {sum(s['f'] * s['h'] for s in m_in):.6f}, out {c1['f'] * c1['h']:.6f} W")
    # Compressor and expander, on their own flowsheet.
    cs = Sheet(auto)
    gas = cs.feed("gas", 350.0, 20e5, 10.0, Z)
    comp_ = cs.add("Compressor", "COMP")
    qc = cs.add("EnergyStream", "Qc")
    hot = cs.add("MaterialStream", "hot")
    exp_ = cs.add("Expander", "EXP")
    qe = cs.add("EnergyStream", "Qe")
    cold = cs.add("MaterialStream", "cold")
    cs.link(gas, comp_)
    cs.link(qc, comp_)
    cs.link(comp_, hot)
    cs.link(hot, exp_)
    cs.link(exp_, cold)
    cs.link(exp_, qe)
    co, eo = comp_.GetAsObject(), exp_.GetAsObject()
    co.CalcMode = co.CalcMode.__class__.OutletPressure
    co.POut = 60e5
    eo.CalcMode = eo.CalcMode.__class__.OutletPressure
    eo.POut = 15e5
    cs.solve()
    print(f"  compressor and expander: unit errors {cs.unit_errors([comp_, exp_])}")
    g, h, k = cs.state(gas), cs.state(hot), cs.state(cold)
    for label, unit, i, out in [("compressor", co, g, h), ("expander", eo, h, k)]:
        mwi = sum(zi * m for zi, m in zip(i["z"], cs.mw))
        hs = cs.flash("PressureEntropy", out["p"], i["s"] / mwi, i["z"])["h"]
        e = unit.AdiabaticEfficiency / 100.0
        model = (hs - i["h"]) / e if label == "compressor" else e * (hs - i["h"])
        print(f"  {label}: outlet H − inlet H {out['h'] - i['h']:.6f} J/mol; isentropic ΔH_s {hs - i['h']:.6f}, "
              f"so the model gives {model:.6f}; power {unit.DeltaQ * 1e3 / i['f']:.6f} J/mol; "
              f"outlet T {out['t']:.6f} K; efficiency {unit.AdiabaticEfficiency}%")

    # The drum's modes, with two inlets at different T and P.
    for mode in ["Legacy", "Adiabatic"]:
        ds = Sheet(auto)
        f1 = ds.feed("f1", 300.0, 30e5, 10.0, Z)
        f2 = ds.feed("f2", 220.0, 25e5, 10.0, Z)
        dr = ds.add("Vessel", "D")
        v_, l_ = ds.add("MaterialStream", "v"), ds.add("MaterialStream", "l")
        ds.link(f1, dr, -1, 0)
        ds.link(f2, dr, -1, 1)
        ds.link(dr, v_, 0, -1)
        ds.link(dr, l_, 1, -1)
        do = dr.GetAsObject()
        do.CalculationMode = getattr(do.CalculationMode.__class__, mode)
        ds.solve()
        a1, a2, vv, ll = ds.state(f1), ds.state(f2), ds.state(v_), ds.state(l_)
        h_in = a1["f"] * a1["h"] + a2["f"] * a2["h"]
        h_out = vv["f"] * vv["h"] + ll["f"] * ll["h"]
        ph = ds.flash("PressureEnthalpy", min(a1["p"], a2["p"]),
                      h_in / (a1["f"] + a2["f"]) / sum(zi * m for zi, m in zip(Z, ds.mw)), Z)
        print(f"  drum, {mode}: unit errors {ds.unit_errors([dr])}; outlet T {vv['t']:.6f} K, P {vv['p']:.1f} Pa; "
              f"an adiabatic mix at the lower P gives {ph['t']:.6f} K; H out − H in {h_out - h_in:.6f} W; "
              f"flows in {a1['f'] + a2['f']:.6f}, out {vv['f'] + ll['f']:.6f} mol/s (v {vv['f']:.6f}, l {ll['f']:.6f})")


def open_loop(auto, sh, o):
    """The loop cut, fed the converged recycle stream: it must come back."""
    rb = sh.state(o["rec_b"])
    print(f"\n4. The open loop: the loop cut, fed the converged recycle stream ({rb['f']:.10f} mol/s)")
    op, oo = reference(auto, open_loop=rb)
    dt, errs = op.solve()
    back = op.state(oo["rec_b"])
    d_f = max(abs(x - y) for x, y in zip(back["flows"], rb["flows"])) / rb["f"]
    print(f"  solved in {dt:.2f} s, unit errors {op.unit_errors([oo[u] for u in UNITS])}")
    print(f"  handed back: flows within {d_f:.2e} of the stream fed, T within {abs(back['t'] - rb['t']):.2e} K, "
          f"P within {abs(back['p'] - rb['p']):.2e} Pa")


if __name__ == "__main__":
    main()
