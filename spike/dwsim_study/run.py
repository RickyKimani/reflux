# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""M8 step 0: DWSIM 9's sensitivity analysis, Adjust, optimiser and regression,
headless.

Run from the repository root, on Windows with DWSIM 9 installed:

    uv run --python 3.13 spike/dwsim_study/run.py [sections]

with sections among 1 to 6 and 50: 1 to 6 by default. Section 50 is the
owner's 50-point study, twice (DWSIM's own form, then a loop of solves), and
takes about twenty minutes.

It answers step 0's questions (docs/m8-plan.md) and prints every number
docs/spike-dwsim-study.md quotes.

None of the four has an API. Each is a window of DWSIM's classic app, and its
work is that window's button handler. So this script makes the app's main
window and its flowsheet window without showing them (DWSIM's first
automation class does that), fills each form's own case, and calls the
handler its button would. Nothing here re-implements DWSIM's study, optimiser
or regression: what runs is DWSIM's code. A single Adjust is also replayed,
with the two MathNet calls its panel makes, to count its flowsheet solves at
each tolerance (section 3).

1. The sensitivity analysis on the owner's reforming PFR (M9's Tutorial 1,
   built by spike/dwsim_pfr/run.py): a few points by the form, and the same
   points by a loop of solves.
2. A failed point in a study: what the form records.
3. Adjust on the reference flowsheet (M4's, with its recycle): a single
   Adjust as its control panel solves it, and the simultaneous adjust solver.
4. The optimiser on a two-stage compression: every method, from the form.
5. The regression: a PR kij from P-x-y data.
6. The solver's modes: does "background parallel" calculate in parallel?
50. The owner's study, 50 points.
"""

import importlib.util
import math
import os
import pathlib
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
SECTIONS = set(sys.argv[1:]) or {"1", "2", "3", "4", "5", "6"}
sys.argv = sys.argv[:1]  # the spikes imported below read their own sections from it


def module(name, path):
    """Another spike's script, under a name of its own: each is a run.py."""
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


pfr = module("dwsim_pfr", "spike/dwsim_pfr/run.py")  # the owner's case
ref = module("dwsim_flowsheet", "spike/dwsim_flowsheet/run.py")  # the reference flowsheet
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader

T_LO, T_HI = 973.15, 1273.15  # the owner's study: the feed from 700 to 1000 °C
VAPOUR_TARGET = 60.0  # mol/s: the spec of section 3, on the reference flowsheet's drum vapour
SPEC_LO, SPEC_HI = 200.0, 260.0  # K: the cooler's outlet T may vary between these

# Section 4's flowsheet: a lean gas compressed from 1 to 16 bar in two stages,
# cooled back to 300 K between them. The work is least near √(1·16) = 4 bar.
GAS = ["Methane", "Ethane", "Propane"]
GAS_Z = [0.90, 0.07, 0.03]
GAS_T, GAS_P, GAS_F, P_OUT = 300.0, 1e5, 10.0, 16e5
MID_LO, MID_HI, MID_START = 1.5e5, 12e5, 8e5


class Dwsim:
    """Both of DWSIM's automation classes. Automation3 builds flowsheets; the
    first, Automation, makes the classic app's main window, unshown, which the
    forms need."""

    def __init__(self):
        gd.load_dwsim()
        import clr

        for name in ["DWSIM", "System.Windows.Forms", "DWSIM.MathOps.DotNumerics", "Ciloci.Flee",
                     "MathNet.Numerics", "DWSIM.FlowsheetSolver"]:
            clr.AddReference(name)
        from DWSIM.Automation import Automation, Automation3
        from System.Reflection import BindingFlags

        self.clr = clr
        self.all = (BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static)
        self.a3 = Automation3()
        t0 = time.perf_counter()
        self.classic = Automation()
        self.classic_time = time.perf_counter() - t0

    def call(self, obj, name, *args):
        """A method by name, private ones too: a button's handler."""
        from System import Array, Object

        m = obj.GetType().GetMethod(name, self.all)
        assert m is not None, name
        return m.Invoke(obj, Array[Object](list(args)))

    def field(self, obj, name):
        return obj.GetType().GetField(name, self.all).GetValue(obj)

    def set_field(self, obj, name, value):
        obj.GetType().GetField(name, self.all).SetValue(obj, value)

    def to_classic(self, sim):
        """A flowsheet built by Automation3, saved and loaded into the classic
        app as a FormFlowsheet, and made its active simulation, which the
        forms read. A FormFlowsheet made bare can't have objects added."""
        path = os.path.join(tempfile.mkdtemp(prefix="reflux_dwsim_"), "case.dwxmz")
        self.a3.SaveFlowsheet(sim, path, True)
        fs = self.classic.LoadFlowsheet(path).__implementation__
        my = self.clr.GetClrType(type(fs)).Assembly.GetType("DWSIM.My.MyProject")
        app = my.GetProperty("Application", self.all).GetValue(None)
        app.GetType().GetProperty("ActiveSimulation", self.all).SetValue(app, fs)
        fs.Options.NumberFormat = "G17"  # the forms' tables, at full precision
        return fs

    def pump(self, done, timeout):
        """Run the window messages the forms post, until done() or timeout."""
        from System.Windows.Forms import Application

        t0 = time.perf_counter()
        while time.perf_counter() - t0 < timeout:
            Application.DoEvents()
            if done():
                Application.DoEvents()
                return True
            time.sleep(0.02)
        return False


def by_tag(fs):
    return {o.GraphicObject.Tag: o for o in fs.SimulationObjects.Values}


def solve_classic(dw, fs):
    """Solve as CalculateFlowsheet2 does, and list what failed: the solver's
    errors and every unit that isn't calculated, since a unit fails silently."""
    t0 = time.perf_counter()
    errs = dw.classic.CalculateFlowsheet2(fs)
    dt = time.perf_counter() - t0
    bad = [str(e).splitlines()[0][:160] for e in errs] if errs else []
    for o in fs.SimulationObjects.Values:
        if not o.GraphicObject.Calculated or (o.ErrorMessage or ""):
            bad.append(f"{o.GraphicObject.Tag}: calculated {o.GraphicObject.Calculated}, {o.ErrorMessage!r}"[:200])
    return dt, bad


def nullable(x):
    from System import Double, Nullable

    return Nullable[Double](x)


def run_study(dw, fs, iv, lo, hi, points, deps, timeout):
    """DWSIM's sensitivity analysis, by its own form. iv is (tag, property
    id, unit) and deps a list of the same. Returns the rows of the form's
    results table, the seconds from its Run button to "Done!", and its log."""
    from DWSIM.SharedClasses.Flowsheet.Optimization import SAVariable, SensitivityAnalysisCase
    from System import EventArgs

    objs = by_tag(fs)
    case = SensitivityAnalysisCase()
    case.name = "study"
    case.iv1.objectID, case.iv1.objectTAG = objs[iv[0]].Name, iv[0]
    case.iv1.propID, case.iv1.unit = iv[1], iv[2]
    case.iv1.lowerlimit, case.iv1.upperlimit, case.iv1.points = nullable(lo), nullable(hi), points
    for k, (tag, pid, unit) in enumerate(deps):
        v = SAVariable()
        v.objectID, v.objectTAG, v.propID, v.unit, v.id, v.name = objs[tag].Name, tag, pid, unit, f"dv{k}", f"y{k}"
        case.depvariables.Add(v.id, v)
    fs.Collections.OPT_SensAnalysisCollection.Clear()
    fs.Collections.OPT_SensAnalysisCollection.Add(case)
    form = dw.field(fs, "FormSensAnalysis0")
    dw.call(form, "FormSensAnalysis_Load", form, EventArgs.Empty)
    form.lbCases.SelectedIndex = 0
    assert form.btnRun.Enabled and form.dgDepVariables.Rows.Count == len(deps), "the form didn't take the case"
    t0 = time.perf_counter()
    dw.call(form, "btnRun_Click", form, EventArgs.Empty)
    finished = dw.pump(lambda: "Done!" in form.tbStats.Text, timeout)
    dt = time.perf_counter() - t0
    g = form.dgvResults
    rows = [[float(str(g.Rows[r].Cells[c].Value)) for c in range(g.Columns.Count) if c != 1]
            for r in range(g.Rows.Count) if g.Rows[r].Cells[0].Value is not None]
    return rows, dt, form.tbStats.Text.replace("\r", ""), finished


def pfr_classic(dw, dv=pfr.DV):
    """The owner's case in the classic app."""
    sh = pfr.Sheet(dw.a3)
    sh.o.dV = dv
    return dw.to_classic(sh.sim)


CH4 = ("product", "PROP_MS_104/Methane", "mol/s")  # the product's methane, as a molar flow
FEED_T = ("feed", "PROP_MS_0", "K")
FEED_P = ("feed", "PROP_MS_1", "Pa")


def loop(dw, fs, iv, values, dep):
    """The study as M9 step 11 ran it: set the value, solve, read."""
    objs = by_tag(fs)
    out = []
    t0 = time.perf_counter()
    for v in values:
        objs[iv[0]].SetPropertyValue(iv[1], v)
        _, bad = solve_classic(dw, fs)
        out.append((v, float(objs[dep[0]].GetPropertyValue(dep[1])), bad))
    return out, time.perf_counter() - t0


def section_1(dw, points=4):
    print(f"\n1. The sensitivity analysis: the owner's study, {points} points, by DWSIM's own form")
    fs = pfr_classic(dw)
    dt0, bad = solve_classic(dw, fs)
    print(f"  the case loaded into the classic app and solved in {dt0:.2f} s, failures {bad}")
    rows, dt, log, finished = run_study(dw, fs, FEED_T, T_LO, T_HI, points, [CH4], 120 * points)
    print(f"  the form's Run handler: finished {finished}, {len(rows)} rows in {dt:.1f} s ({dt / points:.2f} s a point)")
    objs = by_tag(fs)
    print(f"  after it, the feed is back at {float(objs['feed'].GetPropertyValue('PROP_MS_0')):.2f} K: the form "
          f"restores a snapshot of the flowsheet")
    _, bad = solve_classic(dw, fs)
    print(f"  a solve of that restored flowsheet: failures {bad}")
    values = [T_LO + i * (T_HI - T_LO) / (points - 1) for i in range(points)]
    fresh = pfr_classic(dw)
    solve_classic(dw, fresh)
    lp, dt_loop = loop(dw, fresh, FEED_T, values, CH4)
    worst = 0.0
    for row, (t, n, bad) in zip(rows, lp):
        worst = max(worst, abs(row[1] - n) / n)
        x = (pfr.FEED[0] - row[1]) / pfr.FEED[0]
        if points <= 8 or row is rows[0] or row is rows[-1]:
            print(f"    {row[0] - 273.15:7.2f} °C: CH4 out {row[1]:.10f} mol/s (X {x:.6f}); the loop's {n:.10f}{' ' + str(bad) if bad else ''}")
    print(f"  the same points by a loop of set, solve and read: {dt_loop:.1f} s ({dt_loop / points:.2f} s a point); "
          f"its answers differ from the form's by at most {worst:.1e} (relative)")
    return dt, dt_loop


def section_2(dw):
    print("\n2. A failed point: the feed's pressure swept across where Ergun's drop exceeds it")
    for lo, hi in [(101325.0, 11325.0), (11325.0, 101325.0)]:
        fs = pfr_classic(dw, dv=0.25)
        solve_classic(dw, fs)
        rows, dt, log, finished = run_study(dw, fs, FEED_P, lo, hi, 4, [CH4, ("product", "PROP_MS_1", "Pa")], 600)
        print(f"  from {lo:.0f} to {hi:.0f} Pa: finished {finished}, {len(rows)} rows, {dt:.1f} s; the form's log: "
              f"{' | '.join(x for x in log.splitlines() if x)[:400]}")
        alone = {}
        for row in rows:
            one = pfr_classic(dw, dv=0.25)
            by_tag(one)["feed"].SetPropertyValue("PROP_MS_1", row[0])
            _, bad = solve_classic(dw, one)
            prod = by_tag(one)["product"]
            alone[row[0]] = (float(prod.GetPropertyValue(CH4[1])), bad)
            print(f"    {row[0]:9.1f} Pa: the study records CH4 {row[1]:.8f} mol/s, outlet P {row[2]:.1f} Pa; solved "
                  f"alone from a fresh flowsheet: CH4 {alone[row[0]][0]:.8f}, failures {bad}")


def spec_sheet(dw, tol_scale):
    """The reference flowsheet, its recycle's tolerances scaled."""
    sh, o = ref.reference(dw.a3, tol_scale=tol_scale, max_iter=500)
    dt, errs = sh.solve()
    assert not errs and not sh.unit_errors([o[u] for u in ref.UNITS]), (errs, sh.unit_errors([o[u] for u in ref.UNITS]))
    return sh, o


def section_3(dw):
    from MathNet.Numerics.RootFinding import Brent, Secant
    from scipy.optimize import brentq
    from System import Double, Func

    print(f"\n3. Adjust: the reference flowsheet's cooler outlet T, until the drum's vapour is {VAPOUR_TARGET} mol/s")
    exact = None
    for scale, label in [(ref.TIGHT, f"tight (x{ref.TIGHT:.0e})"), (None, "DWSIM's defaults")]:
        sh, o = spec_sheet(dw, scale)
        cool, vap = o["COOL"].GetAsObject(), o["vapour"].GetAsObject()
        calls = [0]

        def f(t):
            calls[0] += 1
            cool.SetPropertyValue("PROP_CL_2", float(t))
            _, errs = sh.solve()
            bad = sh.unit_errors([o[u] for u in ref.UNITS])
            if errs or bad:
                raise RuntimeError(f"at {t} K: {errs or bad}")
            return float(vap.GetPropertyValue("PROP_MS_3")) - VAPOUR_TARGET

        def reset():
            calls[0] = 0
            return f(ref.COOL_T)

        f0 = reset()
        if exact is None:
            from DWSIM.UnitOperations.SpecialOps import Adjust

            adj = Adjust()
            print(f"  an Adjust's defaults: tolerance {adj.Tolerance}, {adj.MaximumIterations} iterations, step "
                  f"{adj.StepSize}, bounds {adj.MinVal} and {adj.MaxVal} (its panel then takes 0.2 and 2 times the "
                  f"value: {0.2 * ref.COOL_T:.0f} to {2 * ref.COOL_T:.0f} K here)")
            print(f"  at {ref.COOL_T} K the vapour is {f0 + VAPOUR_TARGET:.8f} mol/s")
            t0 = time.perf_counter()
            exact = brentq(f, SPEC_LO, SPEC_HI, xtol=1e-12, rtol=8.9e-16)
            print(f"  the answer, by SciPy's brentq to 1e-12 K on the tight recycle: {exact:.10f} K, {calls[0]} solves, "
                  f"{time.perf_counter() - t0:.1f} s; the residual there {f(exact):+.2e} mol/s")
        print(f"  recycle tolerances: {label}")
        runs = [("secant", 1e-4, 10), ("Brent", 1e-4, 10), ("secant", 1e-10, 100), ("Brent", 1e-10, 100)]
        for method, tol, maxit in runs:
            reset()
            calls[0] = 0
            t0 = time.perf_counter()
            try:
                # As EditingForm_Adjust_ControlPanel.btnIniciar_Click calls them.
                if method == "secant":
                    x = Secant.FindRoot(Func[Double, Double](f), ref.COOL_T, ref.COOL_T * 1.01, SPEC_LO, SPEC_HI, tol, maxit)
                else:
                    x = Brent.FindRoot(Func[Double, Double](f), SPEC_LO, SPEC_HI, tol, maxit)
                n = calls[0]
                res = f(x)
                print(f"    {method:6} tol {tol:.0e}, at most {maxit}: {x:.10f} K ({x - exact:+.2e} from the answer), "
                      f"residual {res:+.2e} mol/s, {n} flowsheet solves, {time.perf_counter() - t0:.1f} s")
            except Exception as e:  # MathNet's NonConvergenceException, or a failed solve
                print(f"    {method:6} tol {tol:.0e}, at most {maxit}: {type(e).__name__}: {str(e).splitlines()[0][:150]} "
                      f"after {calls[0]} solves")
        if scale is None:
            continue
        # The panel's default bounds, 0.2 and 2 times the value.
        reset()
        calls[0] = 0
        try:
            x = Brent.FindRoot(Func[Double, Double](f), 0.2 * ref.COOL_T, 2 * ref.COOL_T, 1e-4, 10)
            print(f"    Brent between the panel's default bounds: {x:.8f} K, {calls[0]} solves")
        except Exception as e:
            print(f"    Brent between the panel's default bounds, {0.2 * ref.COOL_T:.0f} and {2 * ref.COOL_T:.0f} K: "
                  f"{type(e).__name__}: {str(e).splitlines()[0][:200]} after {calls[0]} solves")
        simultaneous(dw, exact)
    panel(dw)


def simultaneous(dw, exact):
    """The solver's own simultaneous adjust: an Adjust object on the flowsheet."""
    from System import Action, String
    from DWSIM.Interfaces import IFlowsheet

    for tol in [1e-4, 1e-10]:
        sh, o = spec_sheet(dw, ref.TIGHT)
        adj = sh.add("OT_Adjust", "ADJ").GetAsObject()
        cool, vap = o["COOL"].GetAsObject(), o["vapour"].GetAsObject()
        m, c = adj.ManipulatedObjectData, adj.ControlledObjectData
        m.ID, m.Name, m.PropertyName = cool.Name, "COOL", "PROP_CL_2"
        c.ID, c.Name, c.PropertyName = vap.Name, "vapour", "PROP_MS_3"
        adj.ManipulatedObject, adj.ControlledObject = cool, vap
        adj.AdjustValue, adj.Tolerance, adj.SimultaneousAdjust = VAPOUR_TARGET, tol, True
        msgs = []
        sh.sim.SetMessageListener(Action[String, IFlowsheet.MessageType](lambda s, t: msgs.append(str(s))))
        t0 = time.perf_counter()
        _, errs = sh.solve()
        dt = time.perf_counter() - t0
        its = [x for x in msgs if "Iteration #" in x or "NSSE" in x]
        other = [x[:160] for x in msgs if x not in its and ("djust" in x or "rror" in x)]
        t = float(cool.GetPropertyValue("PROP_CL_2"))
        res = float(vap.GetPropertyValue("PROP_MS_3")) - VAPOUR_TARGET
        print(f"    simultaneous adjust, tol {tol:.0e}: enabled {sh.sim.FlowsheetOptions.SimultaneousAdjustSolverEnabled}; "
              f"{len(its)} iterations, so {3 * len(its) - 2 if its else 0} flowsheet solves by its source (one, and two "
              f"for the derivative, each); {t:.10f} K ({t - exact:+.2e}), residual {res:+.2e} mol/s, {dt:.1f} s; "
              f"errors {errs}; other messages {other[:3]}")
        if its:
            print(f"      its last message: {its[-1][:120]}")


def panel(dw):
    """Can a single Adjust's control panel be run unshown?"""
    from System import Activator, EventArgs

    sh, o = spec_sheet(dw, ref.TIGHT)
    adj = sh.add("OT_Adjust", "ADJ").GetAsObject()
    cool, vap = o["COOL"].GetAsObject(), o["vapour"].GetAsObject()
    m, c = adj.ManipulatedObjectData, adj.ControlledObjectData
    m.ID, m.Name, m.PropertyName = cool.Name, "COOL", "PROP_CL_2"
    c.ID, c.Name, c.PropertyName = vap.Name, "vapour", "PROP_MS_3"
    adj.ManipulatedObject, adj.ControlledObject = cool, vap
    adj.AdjustValue = VAPOUR_TARGET
    _, errs = sh.solve()
    t = float(cool.GetPropertyValue("PROP_CL_2"))
    print(f"  an Adjust that isn't simultaneous, after a flowsheet solve: the cooler is still at {t} K (errors {errs}): "
          f"the solver doesn't run it; only its control panel's Start button does")
    try:
        fs = dw.to_classic(sh.sim)
        a = by_tag(fs)["ADJ"]
        t_panel = dw.clr.GetClrType(type(a.__implementation__)).Assembly.GetType(
            "DWSIM.UnitOperations.EditingForm_Adjust_ControlPanel")
        p = Activator.CreateInstance(t_panel)
        dw.set_field(p, "myADJ", a.__implementation__)
        dw.call(p, "UI_AdjControlPanelForm_Load", p, EventArgs.Empty)
        dw.call(p, "btnIniciar_Click", p, EventArgs.Empty)
        dw.pump(lambda: p.btnIniciar.Enabled, 120)
        tc = float(by_tag(fs)["COOL"].GetPropertyValue("PROP_CL_2"))
        res = float(by_tag(fs)["vapour"].GetPropertyValue("PROP_MS_3")) - VAPOUR_TARGET
        print(f"  its panel, made unshown and its Start handler called (the secant, its default): the cooler ends at "
              f"{tc:.8f} K, residual {res:+.2e} mol/s, {p.Grid1.Rows.Count} rows in its table, status "
              f"{p.lblStatus.Text!r}")
    except Exception as e:
        inner = getattr(e, "InnerException", None)
        print(f"  its panel, made unshown: {type(e).__name__}: {str(inner.Message if inner else e).splitlines()[0][:200]}")


class Compression:
    """Section 4's flowsheet, on Automation3."""

    def __init__(self, dw):
        sh = ref.Sheet(dw.a3)
        sh.sim.SelectedCompounds.Clear()
        for cname in GAS:
            sh.sim.SelectedCompounds.Add(cname, sh.sim.AvailableCompounds[cname])
        gas = sh.feed("gas", GAS_T, GAS_P, GAS_F, GAS_Z)
        names = [("Compressor", "C1"), ("EnergyStream", "W1"), ("MaterialStream", "hot"), ("Cooler", "COOL"),
                 ("EnergyStream", "Q"), ("MaterialStream", "cooled"), ("Compressor", "C2"), ("EnergyStream", "W2"),
                 ("MaterialStream", "out")]
        o = {"gas": gas}
        for kind, name in names:
            o[name] = sh.add(kind, name)
        sh.link(gas, o["C1"])
        sh.link(o["W1"], o["C1"])
        sh.link(o["C1"], o["hot"])
        sh.link(o["hot"], o["COOL"])
        sh.link(o["COOL"], o["cooled"])
        sh.link(o["COOL"], o["Q"])
        sh.link(o["cooled"], o["C2"])
        sh.link(o["W2"], o["C2"])
        sh.link(o["C2"], o["out"])
        for name, p in [("C1", MID_START), ("C2", P_OUT)]:
            c = o[name].GetAsObject()
            c.CalcMode = c.CalcMode.__class__.OutletPressure
            c.POut = p
        c = o["COOL"].GetAsObject()
        c.CalcMode = c.CalcMode.__class__.OutletTemperature
        c.OutletTemperature, c.DeltaP = GAS_T, 0.0
        self.sh, self.o = sh, o

    def work(self, p_mid):
        self.o["C1"].GetAsObject().SetPropertyValue("PROP_CO_4", float(p_mid))
        _, errs = self.sh.solve()
        bad = self.sh.unit_errors([self.o[u] for u in ["C1", "COOL", "C2"]])
        assert not errs and not bad, (errs, bad)
        return sum(float(self.o[w].GetAsObject().GetPropertyValue("PROP_ES_0")) for w in ["W1", "W2"])


class Drum:
    """Section 4's second flowsheet: the reference flowsheet's feed, cooled
    and flashed. The liquid's propane fraction peaks at a cooler temperature
    between the bounds: one property, so an objective with no expression."""

    def __init__(self, dw):
        sh = ref.Sheet(dw.a3)
        feed = sh.feed("feed", ref.FEED["t"], ref.FEED["p"], ref.FEED["f"], ref.Z)
        o = {"feed": feed}
        for kind, name in [("Cooler", "COOL"), ("EnergyStream", "Q"), ("MaterialStream", "cold"), ("Vessel", "DRUM"),
                           ("MaterialStream", "vapour"), ("MaterialStream", "liquid")]:
            o[name] = sh.add(kind, name)
        sh.link(feed, o["COOL"])
        sh.link(o["COOL"], o["cold"])
        sh.link(o["COOL"], o["Q"])
        sh.link(o["cold"], o["DRUM"])
        sh.link(o["DRUM"], o["vapour"], 0, -1)
        sh.link(o["DRUM"], o["liquid"], 1, -1)
        c = o["COOL"].GetAsObject()
        c.CalcMode = c.CalcMode.__class__.OutletTemperature
        c.OutletTemperature, c.DeltaP = DRUM_START, 0.0
        self.sh, self.o = sh, o

    def value(self, t):
        self.o["COOL"].GetAsObject().SetPropertyValue("PROP_CL_2", float(t))
        _, errs = self.sh.solve()
        bad = self.sh.unit_errors([self.o[u] for u in ["COOL", "DRUM"]])
        assert not errs and not bad, (errs, bad)
        return float(self.o["liquid"].GetAsObject().GetPropertyValue(PROPANE))


DRUM_LO, DRUM_HI, DRUM_START = 200.0, 290.0, 280.0  # K
PROPANE = "PROP_MS_102/Propane"  # a stream's propane, as a mole fraction


def optimise(dw, sim, method, variables, expression, maximise):
    """DWSIM's optimiser, by its own form, with one of its methods. variables
    are (name, tag, property id, unit, kind, lower, upper, start). Returns
    the form, the flowsheet's objects, the evaluations it recorded, the
    seconds, and why its objective function fails if it made none."""
    from DWSIM.SharedClasses.Flowsheet.Optimization import (OPTObjectiveFunctionType, OptimizationCase, OPTType,
                                                           OPTVariable, OPTVariableType)
    from System import Array, Double, EventArgs, Object

    fs = dw.to_classic(sim)
    objs = by_tag(fs)
    case = OptimizationCase()
    case.name = "opt"
    case.solvm = getattr(case.solvm.__class__, method)
    case.type = OPTType.Maximization if maximise else OPTType.Minimization
    if expression:
        case.objfunctype, case.expression = OPTObjectiveFunctionType.Expression, expression
    else:
        case.objfunctype = OPTObjectiveFunctionType.Variable
    start = None
    for name, tag, pid, unit, kind, lo, hi, x0 in variables:
        v = OPTVariable()
        v.id, v.name, v.objectID, v.objectTAG, v.propID, v.unit = name, name, objs[tag].Name, tag, pid, unit
        v.type = getattr(OPTVariableType, kind)
        v.lowerlimit, v.upperlimit = nullable(lo), nullable(hi)
        if kind == "Independent":
            v.initialvalue = v.currentvalue = start = x0
        case.variables.Add(v.id, v)
    fs.Collections.OPT_OptimizationCollection.Clear()
    fs.Collections.OPT_OptimizationCollection.Add(case)
    form = dw.field(fs, "FormOptimization0")
    t0 = time.perf_counter()
    dw.call(form, "FormOptimization_Load", form, EventArgs.Empty)
    dw.set_field(form, "selectedoptcase", case)
    dw.call(form, "PopulateForm", case)
    dw.call(form, "btnRun_Click", form, EventArgs.Empty)
    dw.pump(lambda: form.btnRun.Enabled, 5)
    dt = time.perf_counter() - t0
    used = dw.field(form, "selectedoptcase")
    why = ""
    if used.results.Count == 0:
        # The form catches what its objective throws and writes it to a log
        # window. Call the objective once more, to read it.
        try:
            m = [x for x in form.GetType().GetMethods(dw.all) if x.Name == "FunctionValue"][0]
            m.Invoke(form, Array[Object]([Array[Double]([start])]))
            why = "its objective runs when called alone"
        except Exception as e:
            while getattr(e, "InnerException", None) is not None:
                e = e.InnerException
            where = str(getattr(e, "StackTrace", "")).strip().splitlines()
            line = where[-1].split("\\")[-1] if where else ""
            why = f"{type(e).__name__}: {str(getattr(e, 'Message', e)).splitlines()[0][:90]} ({line})"
    return objs, used.results.Count, dt, why


def section_4(dw):
    from scipy.optimize import minimize_scalar
    from DWSIM.SharedClasses.Flowsheet.Optimization import OptimizationCase

    case0 = OptimizationCase()
    methods = list(case0.solvm.GetType().GetEnumNames())
    print("\n4. The optimiser, by its own form")
    print(f"  a case's defaults: tolerance {case0.tolerance}, {case0.maxits} iterations, the derivative's relative "
          f"step {case0.epsilon}, {case0.numdevscheme} points, barrier multiplier {case0.barriermultiplier}")
    print(f"  its methods: {', '.join(methods)}")

    print("  a. The pressure between two compression stages, for the least work: an expression, w1 + w2")
    cp = Compression(dw)
    t0 = time.perf_counter()
    w0 = cp.work(MID_START)
    r = minimize_scalar(cp.work, bounds=(MID_LO, MID_HI), method="bounded", options={"xatol": 1e-4, "maxiter": 500})
    p_opt, w_opt = r.x, r.fun
    print(f"    at the start, {MID_START / 1e5} bar, the work is {w0:.8f} kW; the optimum by SciPy on DWSIM's solves: "
          f"{p_opt / 1e5:.8f} bar, {w_opt:.10f} kW, {r.nfev} solves ({(time.perf_counter() - t0) / (r.nfev + 1) * 1e3:.0f} "
          f"ms each); the ideal gas's √(P1·P2) is {math.sqrt(GAS_P * P_OUT) / 1e5:.4f} bar")
    h = 50.0
    curv = (cp.work(p_opt + h) - 2 * w_opt + cp.work(p_opt - h)) / h ** 2
    print(f"    the work's curvature there: {curv:.3e} kW/Pa², so 1e-6 kW above the least is "
          f"{math.sqrt(2e-6 / curv):.0f} Pa away")
    for bounds, label in [((MID_LO, MID_HI), "the optimum inside the bounds"), ((6e5, MID_HI), "the optimum on a bound")]:
        print(f"    {label}: {bounds[0] / 1e5} to {bounds[1] / 1e5} bar, from {MID_START / 1e5} bar")
        target = min(max(p_opt, bounds[0]), bounds[1])
        w_target = cp.work(target)
        for method in methods:
            cp.work(MID_START)
            vs = [("x", "C1", "PROP_CO_4", "Pa", "Independent", bounds[0], bounds[1], MID_START),
                  ("w1", "W1", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0),
                  ("w2", "W2", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0)]
            objs, n, dt, why = optimise(dw, cp.sh.sim, method, vs, "w1+w2", False)
            p = float(objs["C1"].GetPropertyValue("PROP_CO_4"))
            w = sum(float(objs[t].GetPropertyValue("PROP_ES_0")) for t in ["W1", "W2"])
            print(f"      {method:24} ends at {p / 1e5:.6f} bar ({p - target:+.3e} Pa from the optimum), work "
                  f"{w - w_target:+.2e} kW above the least; {n} evaluations, {dt:.1f} s{'; ' + why if why else ''}")

    print("  b. A cooler's outlet T, for the most propane in the drum's liquid: one property, no expression")
    dr = Drum(dw)
    scan = [(t, dr.value(t)) for t in [200.0, 215.0, 230.0, 245.0, 260.0, 275.0, 290.0]]
    print("    the liquid's propane fraction: " + ", ".join(f"{t:.0f} K {x:.5f}" for t, x in scan))
    r = minimize_scalar(lambda t: -dr.value(t), bounds=(DRUM_LO, DRUM_HI), method="bounded",
                        options={"xatol": 1e-9, "maxiter": 500})
    t_opt, x_opt = r.x, -r.fun
    h = 0.05
    curv = (dr.value(t_opt + h) - 2 * x_opt + dr.value(t_opt - h)) / h ** 2
    print(f"    the optimum by SciPy on DWSIM's solves: {t_opt:.8f} K, {x_opt:.10f}, {r.nfev} solves; the curvature "
          f"there {curv:.3e} /K², so 1e-8 below the most is {math.sqrt(-2e-8 / curv):.4f} K away")
    for maximise in [True, False]:
        # Downhill from the start the fraction falls to the upper bound, a
        # local minimum; the least over the whole range is at the lower one.
        target = t_opt if maximise else DRUM_HI
        x_target = dr.value(target)
        print(f"    asked to {'maximise' if maximise else 'minimise'} it, from {DRUM_START} K between {DRUM_LO} and "
              f"{DRUM_HI} K: the answer is {target:.6f} K, {x_target:.8f}"
              + ("" if maximise else f" (downhill from the start; at {DRUM_LO} K it's {scan[0][1]:.8f})"))
        for method in methods:
            dr.value(DRUM_START)
            vs = [("x", "COOL", "PROP_CL_2", "K", "Independent", DRUM_LO, DRUM_HI, DRUM_START),
                  ("y", "liquid", PROPANE, "", "Dependent", 0.0, 0.0, 0.0)]
            try:
                objs, n, dt, why = optimise(dw, dr.sh.sim, method, vs, "", maximise)
            except Exception as e:
                inner = getattr(e, "InnerException", None)
                print(f"      {method:24} {type(e).__name__}: {str(inner.Message if inner else e).splitlines()[0][:160]}")
                continue
            t = float(objs["COOL"].GetPropertyValue("PROP_CL_2"))
            x = float(objs["liquid"].GetPropertyValue(PROPANE))
            print(f"      {method:24} ends at {t:.6f} K ({t - target:+.3e} K from the answer), the fraction "
                  f"{x - x_target:+.2e} from its; {n} evaluations, {dt:.1f} s{'; ' + why if why else ''}")


def section_5(dw):
    """A PR kij from P-x-y data, made by thermo's PR with a kij of 0.13."""
    from thermo import ChemicalConstantsPackage, CEOSGas, CEOSLiquid, FlashVL, PRMIX
    from System import Activator, EventArgs
    from System.Collections import ArrayList

    print("\n5. The regression: a PR kij for CO2 and ethane, from P-x-y data")
    kij, t = 0.13, 250.0
    consts, props = ChemicalConstantsPackage.from_IDs(["carbon dioxide", "ethane"])
    kw = dict(Tcs=consts.Tcs, Pcs=consts.Pcs, omegas=consts.omegas, kijs=[[0, kij], [kij, 0]])
    gas = CEOSGas(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases)
    liq = CEOSLiquid(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases)
    flasher = FlashVL(consts, props, liquid=liq, gas=gas)
    data = []
    for x1 in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        r = flasher.flash(T=t, VF=0.0, zs=[x1, 1 - x1])
        data.append((x1, r.gas.zs[0], r.P))
    print(f"  the data: thermo 0.6.1's PR bubble points at {t} K with kij = {kij}: "
          + ", ".join(f"x {x:.1f} y {y:.4f} P {p / 1e5:.4f} bar" for x, y, p in data[::4]))
    asm = dw.clr.AddReference("DWSIM")
    try:
        form = Activator.CreateInstance(asm.GetType("DWSIM.FormDataRegression"))
        dw.call(form, "FormDataRegression_Load", form, EventArgs.Empty)
        combos = {n: [str(i) for i in dw.field(form, "_" + n).Items] for n in
                  ["cbModel", "cbDataType", "cbRegMethod", "cbObjFunc", "cbTunit", "cbPunit"]}
        print(f"  the form's choices: models {combos['cbModel']}; data {combos['cbDataType']}; methods "
              f"{combos['cbRegMethod']}; objectives {combos['cbObjFunc']}")
        RegressionCase = asm.GetType("DWSIM.DWSIM.Optimization.DatRegression.RegressionCase")
        settings = asm.GetType("DWSIM.My.MySettings").GetProperty("Default", dw.all).GetValue(None)
        print(f"  the app's EnableParallelProcessing setting, which picks the objective's flash: "
              f"{settings.EnableParallelProcessing}")
        punits = [u for u in combos["cbPunit"] if u in ("Pa", "bar")]
        print(f"  the form's pressure units: {combos['cbPunit']}")
        runs = [(m, settings.EnableParallelProcessing, u) for u in punits for m in combos["cbRegMethod"]]
        runs += [(combos["cbRegMethod"][0], not settings.EnableParallelProcessing, u) for u in punits]
        for method, parallel, punit in runs:
            settings.EnableParallelProcessing = parallel
            case = Activator.CreateInstance(RegressionCase)
            case.comp1, case.comp2, case.model = "Carbon dioxide", "Ethane", "Peng-Robinson"
            case.datatype = case.datatype.__class__.Pxy
            case.method, case.objfunction = method, combos["cbObjFunc"][0]
            case.tunit, case.punit, case.cunit = "K", punit, ""
            case.iepar1, case.llim1, case.ulim1 = 0.0, -0.5, 0.5
            for name in ["tp", "x1p", "x2p", "yp", "pp", "checkp", "ts", "tl"]:
                setattr(case, name, ArrayList())
            for x1, y1, p in data:
                case.tp.Add(t)
                case.x1p.Add(x1)
                case.x2p.Add(0.0)
                case.yp.Add(y1)
                case.pp.Add(p if punit == "Pa" else p / 1e5)
                case.ts.Add(0.0)
                case.tl.Add(0.0)
                case.checkp.Add(True)
            t0 = time.perf_counter()
            try:
                dw.call(form, "LoadCase", case, True)
                dw.call(form, "btnDoReg_Click", form, EventArgs.Empty)
                dw.pump(lambda: dw.field(form, "_btnDoReg").Enabled, 300)
            except Exception as e:
                inner = getattr(e, "InnerException", None)
                print(f"    {method:30} {type(e).__name__}: {str(inner.Message if inner else e).splitlines()[0][:200]}")
                continue
            dt = time.perf_counter() - t0
            pars = dw.field(form, "regressedparameters")
            found = {str(k): float(pars[k]) for k in pars.Keys} if pars is not None else {}
            text = dw.field(form, "_tbRegResults").Text.replace("\r", "")
            lines = [x for x in text.splitlines() if x.strip()]
            tail = lines[-2:]
            good = [x for x in lines if "Function Value" in x]
            print(f"    {method:30} P in {punit}, parallel {parallel}: {found}, {dw.field(form, 'itn')} iterations, "
                  f"{len(good)} of them with a value, the last: {good[-1][:90] if good else None}; {dt:.1f} s; the "
                  f"form's log ends: {' | '.join(x[:90] for x in tail)}")
            if (method, parallel, punit) == runs[0]:
                cc = dw.field(form, "currcase")
                head = [x for x in text.splitlines() if x.strip()][:3]
                print(f"      the case as the form stored it: {cc.comp1}, {cc.comp2}, {cc.model}, {cc.datatype}, units "
                      f"{cc.tunit!r} {cc.punit!r}, {cc.tp.Count} points, T {list(cc.tp)[:2]}, P {list(cc.pp)[:2]}, x "
                      f"{list(cc.x1p)[:2]}, y {list(cc.yp)[:2]}, used {list(cc.checkp)[:2]}, kij from {cc.iepar1} in "
                      f"{cc.llim1} to {cc.ulim1}, tolerance {cc.tolerance}, {cc.maxits} iterations; its log begins: "
                      f"{' | '.join(x[:110] for x in head)}")
    except Exception as e:
        inner = getattr(e, "InnerException", None)
        print(f"  the form: {type(e).__name__}: {str(inner.Message if inner else e).splitlines()[0][:300]}")


def section_6(dw):
    from DWSIM.GlobalSettings import Settings

    print("\n6. The solver's modes: three of the owner's PFRs side by side, each its own feed")
    print(f"  DWSIM's settings here: SolverMode {Settings.SolverMode}, EnableParallelProcessing "
          f"{Settings.EnableParallelProcessing}, MaxDegreeOfParallelism {Settings.MaxDegreeOfParallelism}, "
          f"{os.cpu_count()} logical processors")
    sh = pfr.Sheet(dw.a3)
    sh.o.dV = 0.25
    units = [sh.pfr]
    for k in range(2):
        feed = sh.add("MaterialStream", f"feed{k}")
        o = feed.GetAsObject()
        o.SetTemperature(pfr.T_FEED)
        o.SetPressure(pfr.P_FEED)
        o.SetMolarFlow(sum(pfr.FEED))
        o.SetOverallComposition(sh.vec(f / sum(pfr.FEED) for f in pfr.FEED))
        unit, prod, q = sh.add("RCT_PFR", f"PFR{k}"), sh.add("MaterialStream", f"product{k}"), sh.add("EnergyStream", f"Q{k}")
        sh.sim.ConnectObjects(feed.GraphicObject, unit.GraphicObject, -1, 0)
        sh.sim.ConnectObjects(unit.GraphicObject, prod.GraphicObject, 0, -1)
        sh.sim.ConnectObjects(q.GraphicObject, unit.GraphicObject, -1, 1)
        u, src = unit.GetAsObject(), sh.o
        pfr.enum(u, "ReactorOperationMode", "Isothermic")
        pfr.enum(u, "ReactorSizingType", "Length")
        u.ReactionSetID = "DefaultSet"
        u.Volume, u.Length, u.dV = src.Volume, src.Length, src.dV
        u.CatalystLoading, u.CatalystVoidFraction, u.CatalystParticleDiameter = (
            src.CatalystLoading, src.CatalystVoidFraction, src.CatalystParticleDiameter)
        u.InternalSolver = 0
        units.append(unit)
    one = pfr.Sheet(dw.a3)
    one.o.dV = 0.25
    one.solve()
    t1 = min(one.solve()[0] for _ in range(2))
    print(f"  one PFR alone: {t1:.2f} s a solve")
    keep = Settings.SolverMode, Settings.EnableParallelProcessing
    for mode, par in [(1, False), (1, True), (2, False), (2, True), (0, True)]:
        Settings.SolverMode, Settings.EnableParallelProcessing = mode, par
        sh.set_feed(pfr.T_FEED + 1.0)  # so that every unit is calculated again
        for k in range(2):
            by = {str(o.GraphicObject.Tag): o for o in sh.sim.SimulationObjects.Values}
            by[f"feed{k}"].SetPropertyValue("PROP_MS_0", pfr.T_FEED + 1.0)
        dt, errs = sh.solve(check=False)
        sh.set_feed(pfr.T_FEED)
        for k in range(2):
            by[f"feed{k}"].SetPropertyValue("PROP_MS_0", pfr.T_FEED)
        dt, errs = sh.solve(check=False)
        ok = all(u.GetAsObject().Calculated for u in units)
        print(f"  SolverMode {mode}, EnableParallelProcessing {par}: three PFRs in {dt:.2f} s ({dt / t1:.2f} times one), "
              f"all calculated {ok}, errors {errs}")
    Settings.SolverMode, Settings.EnableParallelProcessing = keep


def section_50(dw):
    print("\n50. The owner's study, 50 points")
    dt, dt_loop = section_1(dw, points=50)
    print(f"  the form: {dt / 60:.2f} min; the loop: {dt_loop / 60:.2f} min")


def main():
    t0 = time.perf_counter()
    dw = Dwsim()
    print(f"{dw.a3.GetVersion()}; loaded in {time.perf_counter() - t0:.1f} s, the classic app's main window, "
          f"unshown, in {dw.classic_time:.1f} s of it")
    for name, f in [("1", section_1), ("2", section_2), ("3", section_3), ("4", section_4), ("5", section_5),
                    ("6", section_6), ("50", section_50)]:
        if name in SECTIONS:
            f(dw)


if __name__ == "__main__":
    main()
