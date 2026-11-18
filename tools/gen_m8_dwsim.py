# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""Write testdata/m8_dwsim.json: M8's references from DWSIM 9 (D18, as
amended in step 0; criteria 8 and 11), each answer of DWSIM's own model.

Run from the repository root, on Windows with DWSIM 9 installed:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_dwsim.py [part ...]

The parts are spec, compression and series, all three by default; a part
not named keeps what the file has. About twenty minutes in all, most of it
the series reactions'.

- spec: the reference flowsheet (M4's, with its recycle, its tolerances
  scaled by 1e-8 as M4 needed), the cooler's outlet T between 200 and 260 K
  until the first drum's vapour is 60 mol/s. Its answer by SciPy's brentq on
  DWSIM's solves to 1e-12 K, and an Adjust's two methods at 1e-10 as its
  control panel calls them (MathNet's secant and Brent), with their solves;
  the vapour at 230 K and its slope at the answer.
- compression: a lean gas from 1 to 16 bar in two stages, each 75%
  efficient, cooled to 300 K between them, the middle pressure for the least
  work: its optimum by SciPy on DWSIM's solves, and by DWSIM's optimiser's
  own form with bounded Brent, the one method of its 11 that finds one
  (step 0), inside the bounds (1.5 to 12 bar) and on a bound (6 to 12); the
  work at five pressures, and its slope at 6 bar.
- series: criterion 8's series reactions, n-pentane to isopentane to
  neopentane, in DWSIM's PFR, each typed as a heterogeneous catalytic rate
  per kg of catalyst on partial pressures in atm and kmol/(kg h), with 1 kg
  of catalyst per m3, so per volume as reflux's; isothermal at 400 K,
  isobaric, 1 mol/s of n-pentane at 1 bar. Isopentane's outflow at five
  volumes, and its most over the volume by SciPy on DWSIM's solves and by
  the optimiser's form with bounded Brent, 0.1 to 2 m3. Each solve takes
  DWSIM about 16 s, so this part about 15 minutes.

Steps 0's spikes build the flowsheets: spike/dwsim_study/run.py (the forms,
the compression), spike/dwsim_flowsheet/run.py (the reference flowsheet)
and spike/dwsim_pfr/run.py (DWSIM's PFR). DWSIM's answers are the same from
run to run, so the file reproduces byte for byte.
"""

import importlib.util
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m8_dwsim.json"
PARTS = set(sys.argv[1:]) or {"spec", "compression", "series"}
sys.argv = sys.argv[:1]
sys.path.insert(0, str(ROOT / "tools"))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


study = module("dwsim_study", "spike/dwsim_study/run.py")
ref, pfr = study.ref, study.pfr

R = 8.31446261815324
SERIES = {"T": 400.0, "P": 1e5, "F": 1.0, "A1": 70.0, "E1": 50000.0, "A2": 0.09, "E2": 30000.0,
          "lo": 0.1, "hi": 2.0, "start": 1.0}


def spec(dw):
    from MathNet.Numerics.RootFinding import Brent, Secant
    from scipy.optimize import brentq
    from System import Double, Func

    sh, o = study.spec_sheet(dw, ref.TIGHT)
    cool, vap = o["COOL"].GetAsObject(), o["vapour"].GetAsObject()
    calls = [0]

    def vapour(t):
        calls[0] += 1
        cool.SetPropertyValue("PROP_CL_2", float(t))
        _, errs = sh.solve()
        bad = sh.unit_errors([o[u] for u in ref.UNITS])
        if errs or bad:
            raise RuntimeError(f"at {t} K: {errs or bad}")
        return float(vap.GetPropertyValue("PROP_MS_3"))

    def f(t):
        return vapour(t) - study.VAPOUR_TARGET

    at_230 = vapour(ref.COOL_T)
    exact = brentq(f, study.SPEC_LO, study.SPEC_HI, xtol=1e-12, rtol=8.9e-16)
    h = 0.01
    slope = (vapour(exact + h) - vapour(exact - h)) / (2 * h)
    methods = {}
    for name in ["secant", "brent"]:
        f(ref.COOL_T)
        calls[0] = 0
        if name == "secant":
            x = Secant.FindRoot(Func[Double, Double](f), ref.COOL_T, ref.COOL_T * 1.01, study.SPEC_LO,
                                study.SPEC_HI, 1e-10, 100)
        else:
            x = Brent.FindRoot(Func[Double, Double](f), study.SPEC_LO, study.SPEC_HI, 1e-10, 100)
        n = calls[0]
        methods[name] = {"T": x, "solves": n, "residual": f(x)}
        print(f"  spec, {name}: {x:.10f} K in {n} solves")
    print(f"  spec: {exact:.10f} K; vapour at 230 K {at_230:.8f} mol/s; slope {slope:.6f} mol/(s K)")
    return {"T": exact, "vapour_at_230": at_230, "slope": slope, "target": study.VAPOUR_TARGET,
            "bounds": [study.SPEC_LO, study.SPEC_HI], "methods": methods}


def compression(dw):
    from scipy.optimize import minimize_scalar

    cp = study.Compression(dw)
    for c in ["C1", "C2"]:
        cp.o[c].GetAsObject().AdiabaticEfficiency = 75.0
    work = lambda p: 1e3 * cp.work(p)  # W, from DWSIM's kW
    r = minimize_scalar(work, bounds=(study.MID_LO, study.MID_HI), method="bounded",
                        options={"xatol": 1e-3, "maxiter": 500})
    p_opt, w_opt = float(r.x), float(r.fun)
    grid = [2e5, 4e5, 6e5, 8e5, 12e5]
    at = {f"{p:.0f}": work(p) for p in grid}
    # Each stage on its own: its work and its outlet's T.
    stages = {}
    for p in grid:
        cp.work(p)
        stages[f"{p:.0f}"] = {name: {"work": 1e3 * float(cp.o[w].GetAsObject().GetPropertyValue("PROP_ES_0")),
                                     "T": float(cp.o[out].GetAsObject().GetPropertyValue("PROP_MS_0"))}
                              for name, w, out in [("c1", "W1", "hot"), ("c2", "W2", "out")]}
    flash = {str(k): str(cp.sh.pp.FlashSettings[k]) for k in cp.sh.pp.FlashSettings.Keys}
    h = 50.0
    slope_6 = (work(6e5 + h) - work(6e5 - h)) / (2 * h)
    curvature = (work(p_opt + h) - 2 * w_opt + work(p_opt - h)) / h ** 2
    forms = {}
    for label, lo, hi in [("inside", study.MID_LO, study.MID_HI), ("bound", 6e5, study.MID_HI)]:
        cp.work(study.MID_START)
        vs = [("x", "C1", "PROP_CO_4", "Pa", "Independent", lo, hi, study.MID_START),
              ("w1", "W1", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0),
              ("w2", "W2", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0)]
        objs, n, _, why = study.optimise(dw, cp.sh.sim, "AL_BRENT_B", vs, "w1+w2", False)
        p = float(objs["C1"].GetPropertyValue("PROP_CO_4"))
        w = 1e3 * sum(float(objs[t].GetPropertyValue("PROP_ES_0")) for t in ["W1", "W2"])
        forms[label] = {"bounds": [lo, hi], "P": p, "work": w, "evaluations": n, "why": why}
        print(f"  compression, bounded Brent {label}: {p / 1e5:.8f} bar, {w:.6f} W, {n} evaluations")
    print(f"  compression: {p_opt / 1e5:.8f} bar, {w_opt:.6f} W, curvature {curvature:.3e} W/Pa2, "
          f"slope at 6 bar {slope_6:.6e} W/Pa")
    return {"P": p_opt, "work": w_opt, "curvature": curvature, "work_at": at, "stages_at": stages,
            "slope_at_6_bar": slope_6, "efficiency": 0.75, "forms": forms, "flash_settings": flash}


class Series:
    """Criterion 8's series reactions in DWSIM's PFR."""

    def __init__(self, dw):
        from System import Array, Double, String
        from System.Collections.Generic import Dictionary
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType

        s = SERIES
        self.dw, self.sim = dw, dw.a3.CreateFlowsheet()
        comps = ["N-pentane", "Isopentane", "Neopentane"]
        for c in comps:
            compound = self.sim.AvailableCompounds[c]
            assert compound.OriginalDB == "ChemSep", (c, compound.OriginalDB)
            self.sim.SelectedCompounds.Add(c, compound)
        self.sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)")
        objs = {}
        x = 0
        for kind, name in [("MaterialStream", "feed"), ("RCT_PFR", "PFR"), ("MaterialStream", "product"),
                           ("EnergyStream", "Q")]:
            x += 100
            objs[name] = self.sim.AddObject(getattr(ObjectType, kind), x, 100, name)
        feed = objs["feed"].GetAsObject()
        feed.SetTemperature(s["T"])
        feed.SetPressure(s["P"])
        feed.SetMolarFlow(s["F"])
        feed.SetOverallComposition(Array[Double]([1.0, 0.0, 0.0]))
        # k·p per m³ in mol/(m³ s) with p in Pa is, per kg on 1 kg/m³, k·101325·3.6
        # kmol/(kg h) with p in atm.
        for name, a, e, comp in [("first", s["A1"], s["E1"], [("N-pentane", -1.0), ("Isopentane", 1.0)]),
                                 ("second", s["A2"], s["E2"], [("Isopentane", -1.0), ("Neopentane", 1.0)])]:
            d = Dictionary[String, Double]()
            for c, nu in comp:
                d[c] = nu
            rx = self.sim.CreateHetCatReaction(name, "", d, comp[0][0], "Vapor", "Partial Pressure", "atm",
                                               "kmol/[kg.h]", f"{a}*exp(-{e}/{R}/T)*101325*3.6*R1", "1")
            rx.Tmin, rx.Tmax = 0.0, 2000.0
            self.sim.AddReaction(rx)
            self.sim.AddReactionToSet(rx.ID, "DefaultSet", True, 0)
        self.sim.ConnectObjects(objs["feed"].GraphicObject, objs["PFR"].GraphicObject, -1, 0)
        self.sim.ConnectObjects(objs["PFR"].GraphicObject, objs["product"].GraphicObject, 0, -1)
        self.sim.ConnectObjects(objs["Q"].GraphicObject, objs["PFR"].GraphicObject, -1, 1)
        o = objs["PFR"].GetAsObject()
        pfr.enum(o, "ReactorOperationMode", "Isothermic")
        pfr.enum(o, "ReactorSizingType", "Length")
        o.ReactionSetID = "DefaultSet"
        o.Volume, o.Length = s["start"], 1.0
        o.CatalystLoading = 1.0
        o.UseUserDefinedPressureDrop, o.UserDefinedPressureDrop = True, 0.0
        o.InternalSolver = 0
        self.objs, self.o = objs, o

    def outflow(self, v):
        self.o.Volume = float(v)
        errs = self.dw.a3.CalculateFlowsheet2(self.sim)
        o = self.o
        if errs or not o.Calculated or (o.ErrorMessage or ""):
            raise RuntimeError(f"at {v} m3: {[str(e) for e in errs] if errs else o.ErrorMessage}")
        product = self.objs["product"].GetAsObject()
        return float(product.Phases[0].Compounds["Isopentane"].MolarFlow)


def series(dw):
    from scipy.optimize import minimize_scalar

    s = SERIES
    sh = Series(dw)
    at = {}
    for v in [0.1, 0.5, 0.655, 1.0, 2.0]:
        at[f"{v}"] = sh.outflow(v)
        print(f"  series: {v} m3, {at[f'{v}']!r} mol/s", flush=True)
    r = minimize_scalar(lambda v: -sh.outflow(v), bounds=(s["lo"], s["hi"]), method="bounded",
                        options={"xatol": 1e-7, "maxiter": 500})
    v_opt, b_opt = float(r.x), float(-r.fun)
    # Near the peak, where it's flattest: does DWSIM's outflow move smoothly?
    near = {}
    for v in [0.64, 0.645, 0.65, 0.655, 0.66, v_opt, 0.665, 0.67]:
        near[f"{v!r}"] = sh.outflow(v)
        print(f"  series, near the peak: {v!r} m3, {near[f'{v!r}']!r} mol/s", flush=True)
    sh.outflow(s["start"])
    vs = [("x", "PFR", "PROP_PF_2", "m3", "Independent", s["lo"], s["hi"], s["start"]),
          ("y", "product", "PROP_MS_104/Isopentane", "mol/s", "Dependent", 0.0, 0.0, 0.0)]
    objs, n, _, why = study.optimise(dw, sh.sim, "AL_BRENT_B", vs, "", True)
    v = float(objs["PFR"].GetPropertyValue("PROP_PF_2"))
    b = float(objs["product"].GetPropertyValue("PROP_MS_104/Isopentane"))
    print(f"  series: {v_opt:.8f} m3, {b_opt:.10f} mol/s by SciPy; bounded Brent {v:.8f} m3, {b:.10f} mol/s, "
          f"{n} evaluations; dV {sh.o.dV}")
    return {"setup": s, "dV": float(sh.o.dV), "outflow_at": at, "near_peak": near, "V": v_opt, "B": b_opt,
            "form": {"V": v, "B": b, "evaluations": n, "why": why}}


def main():
    old = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    dw = study.Dwsim()
    out = {
        "provenance": "tools/gen_m8_dwsim.py: DWSIM 9.0.5 through pythonnet 3.0.5, its forms driven unshown "
        "(docs/spike-dwsim-study.md); SciPy 1.14.1 for the answers of DWSIM's own model. spec: the reference "
        "flowsheet, its recycle's tolerances scaled by 1e-8; compression: DWSIM's own compounds and kij, as "
        "testdata/flowsheets/compression_dwsim.toml has them; series: DWSIM's PFR on typed rates per kg on "
        "1 kg/m3",
        "dwsim": str(dw.a3.GetVersion()),
    }
    for name, f in [("spec", spec), ("compression", compression), ("series", series)]:
        out[name] = f(dw) if name in PARTS else old.get(name)
    OUT.write_bytes((json.dumps(out, indent=1) + "\n").encode())
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
