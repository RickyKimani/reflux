# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Generate DWSIM 9's references for M4 (D12, exit criteria 2 and 5).

Run from the repository root, on Windows with DWSIM 9 installed, and V on the
PATH for the open loop:

    uv run --python 3.13 tools/gen_m4_dwsim.py [units] [reference]

with the parts to run (both by default). It writes:

testdata/m4_dwsim/components/*.toml
    DWSIM's constants for the five compounds, as reflux component files: its
    Tc, Pc, omega, molar mass and ideal-gas Cp (ChemSep's eq. 16), with no
    formation data, so that H is from the ideal gas at 298.15 K as DWSIM's is.
testdata/flowsheets/reference_dwsim.toml
    The reference flowsheet on those files and DWSIM's kij.
testdata/m4_dwsim.json
    Each unit operation alone, from a stated inlet, and the reference
    flowsheet, with every stream and duty as DWSIM 9 gives them; then the
    open loop: the loop cut and fed reflux's own converged recycle, which
    DWSIM must hand back.

It does what "What step 8 must do" in docs/spike-dwsim-flowsheet.md asks:
every unit's Calculated and ErrorMessage and every stream are checked before a
number is used, the splitter has a ratio per port, the torn stream an
estimate, and the recycle's tolerances are scaled by 1e-8 with acceleration
off. CI never runs it: like M3's, its output is committed.
"""

import datetime
import json
import math
import pathlib
import platform
import subprocess
import sys
import tempfile
import time
import tomllib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m4_dwsim.json"
COMPONENTS_DIR = ROOT / "testdata" / "m4_dwsim" / "components"
FLOWSHEET = ROOT / "testdata" / "flowsheets" / "reference_dwsim.toml"
REFERENCE = ROOT / "testdata" / "flowsheets" / "reference.toml"

# DWSIM's names, and the file's: lookup ignores case, but the file's are the
# reference flowsheet's.
COMPS = ["Methane", "Ethane", "Propane", "N-butane", "N-pentane"]
NAMES = ["methane", "ethane", "propane", "n-butane", "n-pentane"]
Z = [0.70, 0.10, 0.08, 0.07, 0.05]
# A liquid like the reference drum's, subcooled at 220 K and 30 bar.
ZL = [0.287, 0.174, 0.201, 0.195, 0.143]
TOL = "1E-8"
TIGHT = 1e-8  # the scale on the recycle's tolerances, past which nothing improves
TOL_FIELDS = ["Temperatura", "Pressao", "VazaoMassica", "FracaoVapor", "Entalpia", "Entropia", "Composicao"]
PARTS = set(sys.argv[1:]) or {"units", "reference"}


class Sheet:
    """One DWSIM flowsheet with PR, its flash tolerances at 1e-8 as in M3."""

    def __init__(self, auto):
        from System import Array, Double
        from DWSIM.Interfaces.Enums import FlashSetting
        from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType

        self.auto, self.OT = auto, ObjectType
        self.vec = lambda xs: Array[Double](list(xs))
        self.sim = auto.CreateFlowsheet()
        for c in COMPS:
            compound = self.sim.AvailableCompounds[c]
            assert compound.OriginalDB == "ChemSep", (c, compound.OriginalDB)
            self.sim.SelectedCompounds.Add(c, compound)
        self.pp = self.sim.CreateAndAddPropertyPackage("Peng-Robinson (PR)").__implementation__
        for key in (FlashSetting.PTFlash_External_Loop_Tolerance, FlashSetting.PTFlash_Internal_Loop_Tolerance,
                    FlashSetting.PHFlash_External_Loop_Tolerance, FlashSetting.PHFlash_Internal_Loop_Tolerance):
            self.pp.FlashSettings[key] = TOL
        self.x = 0
        self.units = []

    def add(self, kind, name):
        self.x += 100
        obj = self.sim.AddObject(getattr(self.OT, kind), self.x, 100, name)
        if kind not in ("MaterialStream", "EnergyStream"):
            self.units.append(obj)
        return obj

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
        """Solve, and refuse any failure: DWSIM's solver reports none of a unit's."""
        t0 = time.perf_counter()
        errs = self.auto.CalculateFlowsheet2(self.sim)
        dt = time.perf_counter() - t0
        bad = [str(e) for e in errs] if errs else []
        for u in self.units:
            o = u.GetAsObject()
            if not o.Calculated or (o.ErrorMessage or ""):
                bad.append(f"{u.GraphicObject.Tag}: calculated {o.Calculated}, {o.ErrorMessage!r}")
        if bad:
            raise RuntimeError(f"DWSIM failed: {bad}")
        return dt

    def state(self, stream):
        """A stream as DWSIM has it, refused if any of it wasn't calculated."""
        m = stream.GetAsObject()
        ph = m.Phases[0].Properties
        values = [ph.temperature, ph.pressure, ph.molarflow, ph.molar_enthalpy, ph.molar_entropy, ph.density,
                  ph.molecularWeight]
        if any(v is None for v in values):
            raise RuntimeError(f"{stream.GraphicObject.Tag} wasn't calculated")
        vf = m.Phases[2].Properties.molarfraction or 0.0
        lf = m.Phases[3].Properties.molarfraction or 0.0
        flows = [m.Phases[0].Compounds[c].MolarFlow for c in COMPS]
        out = {"t": ph.temperature, "p": ph.pressure, "f": ph.molarflow, "flows": flows,
               "z": [m.Phases[0].Compounds[c].MoleFraction for c in COMPS], "beta": vf, "h": ph.molar_enthalpy,
               "s": ph.molar_entropy, "v": ph.molecularWeight / ph.density * 1e-3}
        out["y"] = [m.Phases[2].Compounds[c].MoleFraction for c in COMPS] if vf > 0 else []
        out["x"] = [m.Phases[3].Compounds[c].MoleFraction for c in COMPS] if lf > 0 else []
        if any(v is None for v in flows + out["z"] + out["y"] + out["x"]):
            raise RuntimeError(f"{stream.GraphicObject.Tag} has a property that wasn't calculated")
        return out

    def constants(self):
        """DWSIM's constants and kij. The package needs a stream to read them from."""
        if self.pp.CurrentMaterialStream is None:
            self.feed("constants", 300.0, 1e5, 1.0, Z)
        out = []
        for name, c in zip(NAMES, COMPS):
            d = self.sim.SelectedCompounds[c]
            assert d.IdealgasCpEquation == "16", (c, d.IdealgasCpEquation)
            # g/mol to kg/mol by moving the decimal point, as reflux's importer does:
            # dividing by 1000 would round twice.
            mw = float(repr(d.Molar_Weight) + "e-3")
            out.append({"name": name, "dwsim_name": c, "cas": d.CAS_Number, "tc": d.Critical_Temperature,
                        "pc": d.Critical_Pressure, "omega": d.Acentric_Factor, "mw": mw,
                        "cp": [d.Ideal_Gas_Heat_Capacity_Const_A, d.Ideal_Gas_Heat_Capacity_Const_B,
                               d.Ideal_Gas_Heat_Capacity_Const_C, d.Ideal_Gas_Heat_Capacity_Const_D,
                               d.Ideal_Gas_Heat_Capacity_Const_E]})
        k = self.pp.RET_VKij()
        kij = [[k.GetValue(i, j) for j in range(len(COMPS))] for i in range(len(COMPS))]
        settings = {str(key): str(self.pp.FlashSettings[key]) for key in self.pp.FlashSettings.Keys}
        return out, kij, settings


def enum(obj, prop, name):
    setattr(obj, prop, getattr(getattr(obj, prop).__class__, name))


# Each unit alone. spec is reflux's, as its file format has it; a stream is
# (T K, P Pa, F mol/s, z). The sign is how DWSIM's DeltaQ, in kW, goes into
# the material: a cooler's is heat removed, an expander's power produced.
UNIT_CASES = [
    {"name": "cooler to a T", "dwsim": "Cooler", "sign": -1, "inlets": [(300.0, 30e5, 100.0, Z)],
     "spec": {"type": "heater", "mode": "outlet_t", "value": 230.0, "pressure_drop": 0.0}},
    {"name": "heater by its duty", "dwsim": "Heater", "sign": 1, "inlets": [(230.0, 30e5, 100.0, Z)],
     "spec": {"type": "heater", "mode": "duty", "value": 500e3, "pressure_drop": 0.0}},
    {"name": "heater to a vapour fraction", "dwsim": "Heater", "sign": 1, "inlets": [(230.0, 30e5, 100.0, Z)],
     "spec": {"type": "heater", "mode": "outlet_vapour_fraction", "value": 0.9, "pressure_drop": 0.0},
     "exact": "vapour_fraction"},
    {"name": "heater by a rise, with a pressure drop", "dwsim": "Heater", "sign": 1,
     "inlets": [(230.0, 30e5, 100.0, Z)],
     "spec": {"type": "heater", "mode": "delta_t", "value": 25.0, "pressure_drop": 0.5e5}},
    {"name": "valve", "dwsim": "Valve", "sign": 1, "inlets": [(220.0, 30e5, 34.8, ZL)],
     "spec": {"type": "valve", "mode": "outlet_p", "value": 5e5}},
    {"name": "pump", "dwsim": "Pump", "sign": 1, "inlets": [(220.0, 30e5, 34.8, ZL)],
     "spec": {"type": "pump", "mode": "outlet_p", "value": 40e5, "efficiency": 0.75}, "exact": "pump"},
    {"name": "compressor", "dwsim": "Compressor", "sign": 1, "inlets": [(350.0, 20e5, 10.0, Z)],
     "spec": {"type": "compressor", "mode": "outlet_p", "value": 60e5, "efficiency": 0.75}},
    {"name": "expander", "dwsim": "Expander", "sign": -1, "inlets": [(350.0, 60e5, 10.0, Z)],
     "spec": {"type": "expander", "mode": "outlet_p", "value": 20e5, "efficiency": 0.8}, "exact": "expander"},
    {"name": "mixer", "dwsim": "NodeIn", "sign": 1, "inlets": [(300.0, 30e5, 100.0, Z), (220.0, 25e5, 34.8, ZL)],
     "spec": {"type": "mixer"}},
    {"name": "splitter", "dwsim": "NodeOut", "sign": 1, "inlets": [(230.0, 30e5, 100.0, Z)],
     "spec": {"type": "splitter", "fractions": [0.3, 0.7]}},
    {"name": "adiabatic drum", "dwsim": "Vessel", "sign": 1, "inlets": [(250.0, 30e5, 100.0, Z)],
     "spec": {"type": "flash_drum", "mode": "duty", "value": 0.0}},
    {"name": "drum at a T", "dwsim": "Vessel", "sign": 1, "inlets": [(250.0, 30e5, 100.0, Z)],
     "spec": {"type": "flash_drum", "mode": "outlet_t", "value": 220.0}},
    # The separator's EnergyImb is what goes in less what comes out.
    {"name": "component separator", "dwsim": "ComponentSeparator", "sign": -1, "inlets": [(300.0, 30e5, 100.0, Z)],
     "spec": {"type": "component_separator", "fractions": [0.9, 0.5, 0.1, 0.0, 0.0]}},
]


def outlets_of(case):
    kind = case["spec"]["type"]
    if kind in ("splitter", "flash_drum", "component_separator"):
        return 2
    return 1


def configure(case, o):
    """A unit's specification, in DWSIM's terms."""
    spec = case["spec"]
    kind = spec["type"]
    if kind == "heater":
        o.DeltaP = spec["pressure_drop"]
        mode = spec["mode"]
        if mode == "outlet_t":
            enum(o, "CalcMode", "OutletTemperature")
            o.OutletTemperature = spec["value"]
        elif mode == "duty":
            enum(o, "CalcMode", "HeatAdded")
            o.DeltaQ = spec["value"] / 1e3
        elif mode == "outlet_vapour_fraction":
            enum(o, "CalcMode", "OutletVaporFraction")
            o.OutletVaporFraction = spec["value"]
        else:
            enum(o, "CalcMode", "TemperatureChange")
            o.DeltaT = spec["value"]
    elif kind == "valve":
        enum(o, "CalcMode", "OutletPressure")
        o.OutletPressure = spec["value"]
    elif kind == "pump":
        enum(o, "CalcMode", "OutletPressure")
        o.Pout = spec["value"]
        o.Efficiency = spec["efficiency"] * 100
    elif kind in ("compressor", "expander"):
        enum(o, "CalcMode", "OutletPressure")
        enum(o, "ProcessPath", "Adiabatic")
        o.POut = spec["value"]
        o.AdiabaticEfficiency = spec["efficiency"] * 100
    elif kind == "mixer":
        enum(o, "PressureCalculation", "Minimum")
    elif kind == "splitter":
        o.Ratios.Clear()
        for r in spec["fractions"] + [0.0]:  # a ratio for each of its three ports
            o.Ratios.Add(r)
    elif kind == "flash_drum":
        if spec["mode"] == "duty":
            enum(o, "CalculationMode", "Adiabatic")
        else:
            enum(o, "CalculationMode", "HeatingCoolingIsothermic")
            o.OverrideT = True
            o.FlashTemperature = spec["value"]
    elif kind == "component_separator":
        configure_separator(o, spec["fractions"])


def configure_separator(o, fractions):
    from DWSIM.UnitOperations.UnitOperations.Auxiliary import ComponentSeparationSpec, SeparationSpec

    o.SpecifiedStreamIndex = 0  # the fractions are of what goes to the first outlet
    o.ComponentSepSpecs.Clear()
    for c, frac in zip(COMPS, fractions):
        o.ComponentSepSpecs.Add(c, ComponentSeparationSpec(c, SeparationSpec.PercentInletMolarFlow, frac * 100, "%"))


class Models:
    """The same Peng-Robinson two ways, from M3's tool (gen_m3_dwsim.py): exact,
    with DWSIM's constants and kij and an exact quadrature of its ideal-gas Cp;
    and rounded, with DWSIM's Omega_a, Omega_b and R and its own ideal gas,
    which reproduces each of DWSIM's phases to 1e-12."""

    def __init__(self, auto):
        system = gd.SYSTEMS[1]  # M3's natural gas: these five compounds, in this order
        assert system["compounds"] == COMPS, system
        self.dw = gd.Dwsim(auto, system)
        self.flasher = gd.g3.make_flasher(gd.PRMIX, self.dw.comps, self.dw.kij)

    def exact(self, t, p, z):
        """Molar H and S, and the phases, at (T, P, z)."""
        pairs, _ = gd.exact_pt(self.dw, t, p, z)
        h, s = gd.exact_hs(self.dw.comps, pairs, t, p)
        return h, s, pairs

    def rounded_phase(self, comp, t, p, vapor):
        dw = self.dw
        cls = gd.g3.CEOSGas if vapor else gd.g3.CEOSLiquid
        args = dict(Tcs=[c["tc"] for c in dw.comps], Pcs=[c["pc"] for c in dw.comps],
                    omegas=[c["omega"] for c in dw.comps], kijs=dw.kij)
        ph = cls(gd.PRMIXRounded, args, HeatCapacityGases=None).to(T=t, P=p, zs=comp)
        hdep, _ = gd.rounded_dep(ph)
        return dw.pp.RET_Hid(gd.T_REF, t, dw.vec(comp)) * dw.mw(comp) + hdep

    def enthalpies(self, st):
        """A DWSIM state's molar H by the exact PR, at its T, P and z, and by the
        rounded PR at DWSIM's own phases."""
        st["h_exact"], _, _ = self.exact(st["t"], st["p"], st["z"])
        phases = [(st["beta"], st["y"], True), (1.0 - st["beta"], st["x"], False)]
        st["h_rounded"] = sum(f * self.rounded_phase(c, st["t"], st["p"], v) for f, c, v in phases if c and f > 0)
        return st

    @staticmethod
    def beta(pairs):
        """The vapour fraction of phases sorted by density, as reflux has it."""
        if len(pairs) == 2:
            return pairs[0][0]
        return 1.0 if pairs[0][1].__class__.__name__ == "CEOSGas" else 0.0

    def vapour_fraction(self, p, vf, z, lo, hi):
        """The exact PR's state at (P, vapour fraction), by bisection on its PT
        flash's: thermo 0.6.1's own flash at a fraction between 0 and 1 fails on
        a mixture, with `dew` used before it's set."""
        t = self.solve_t(lambda t: self.beta(self.exact(t, p, z)[2]) - vf, lo, hi)
        h, _, pairs = self.exact(t, p, z)
        return {"t": t, "h": h, "beta": self.beta(pairs)}

    def solve_t(self, f, lo, hi):
        """T where f, increasing in T, is 0, by bisection to 1e-10 K."""
        assert f(lo) < 0 < f(hi), (lo, hi)
        while hi - lo > 1e-10:
            mid = 0.5 * (lo + hi)
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    def expander(self, inlet, p_out, efficiency):
        """The exact PR's expander: the isentropic outlet by S, then the actual one by H."""
        z = inlet["z"]
        h_in, s_in, _ = self.exact(inlet["t"], inlet["p"], z)
        t_s = self.solve_t(lambda t: self.exact(t, p_out, z)[1] - s_in, 200.0, inlet["t"])
        h_s = self.exact(t_s, p_out, z)[0]
        h_out = h_in + efficiency * (h_s - h_in)
        t = self.solve_t(lambda t: self.exact(t, p_out, z)[0] - h_out, 200.0, inlet["t"])
        return {"t_s": t_s, "t": t, "beta": self.beta(self.exact(t, p_out, z)[2]), "h_in": h_in, "h_out": h_out}


def run_unit(auto, models, case):
    sh = Sheet(auto)
    feeds = [sh.feed(f"in{i}", *s) for i, s in enumerate(case["inlets"])]
    unit = sh.add(case["dwsim"], "UNIT")
    outs = [sh.add("MaterialStream", f"out{i}") for i in range(outlets_of(case))]
    for i, f in enumerate(feeds):
        sh.link(f, unit, -1, i)
    for i, s in enumerate(outs):
        sh.link(unit, s, i, -1)
    if case["dwsim"] in ("Heater", "Pump", "Compressor", "Vessel"):
        sh.link(sh.add("EnergyStream", "Q"), unit)
    elif case["dwsim"] in ("Cooler", "Expander", "ComponentSeparator"):
        sh.link(unit, sh.add("EnergyStream", "Q"))
    o = unit.GetAsObject()
    configure(case, o)
    sh.solve()
    inlets = [models.enthalpies(sh.state(f)) for f in feeds]
    outlets = [models.enthalpies(sh.state(s)) for s in outs]
    if case["dwsim"] == "ComponentSeparator":
        dq = o.EnergyImb
    elif hasattr(o, "DeltaQ") and o.DeltaQ is not None:
        dq = o.DeltaQ
    else:
        dq = 0.0
    rate = case["sign"] * dq * 1e3
    flow_in = sum(s["f"] for s in inlets)
    h_in = sum(s["f"] * s["h"] for s in inlets)
    h_out = sum(s["f"] * s["h"] for s in outlets)
    comp = max(abs(sum(s["flows"][i] for s in inlets) - sum(s["flows"][i] for s in outlets)) for i in range(len(COMPS)))
    out = {k: case[k] for k in ("name", "dwsim", "spec")}
    # energy is DWSIM's own number for the unit's duty or work, and from_streams
    # what its streams' enthalpies say it is: they differ for the separator.
    out.update({"inlets": inlets, "outlets": outlets, "energy": rate, "from_streams": h_out - h_in,
                "dwsim_delta_q": dq,
                "balance": {"component": comp / flow_in, "energy": (h_out - h_in - rate) / (flow_in * 8.314462618 *
                                                                                          max(s["t"] for s in inlets))}})
    # The exact PR's own answer, for the units that are further from DWSIM than
    # criterion 2's targets, and what explains DWSIM's.
    kind, i = case.get("exact"), inlets[0]
    if kind == "vapour_fraction":
        out["exact"] = models.vapour_fraction(i["p"] - case["spec"]["pressure_drop"], case["spec"]["value"], i["z"],
                                              i["t"], outlets[0]["t"] + 20.0)
    elif kind == "pump":
        _, _, pairs = models.exact(i["t"], i["p"], i["z"])
        assert len(pairs) == 1, "the pump's inlet isn't one phase"
        out["exact"] = {"v": pairs[0][1].V()}
    elif kind == "expander":
        out["exact"] = models.expander(i, case["spec"]["value"], case["spec"]["efficiency"])
        o = outlets[0]
        # DWSIM's ideal mixing entropy less -R sum(x ln x), as M3 measured it.
        out["mixing_error"] = {"inlet": models.dw.mixing_error(i["z"]), "vapour": models.dw.mixing_error(o["y"]),
                               "liquid": models.dw.mixing_error(o["x"])}
    return out


def reference(auto, tol_scale=TIGHT, max_iter=500, open_loop=None):
    """The reference flowsheet, as step 0's spike built it. open_loop, a stream
    (T, P, F, z), cuts the loop and feeds that where the recycle came in."""
    sh = Sheet(auto)
    o = {"feed": sh.feed("feed", 300.0, 30e5, 100.0, Z)}
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
        rc.SetTemperature(230.0)
        rc.SetPressure(31e5)
        rc.SetMolarFlow(20.0)
        rc.SetOverallComposition(sh.vec([0.2, 0.15, 0.25, 0.22, 0.18]))
    else:
        o["rec_c"] = sh.feed("rec_c", *open_loop)
    L = sh.link
    L(o["feed"], o["MIX"], -1, 0)
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
    enum(c, "CalcMode", "OutletTemperature")
    c.OutletTemperature = 230.0
    c.DeltaP = 0.0
    s = o["SPLIT"].GetAsObject()
    s.Ratios.Clear()
    for r in [0.5, 0.5, 0.0]:
        s.Ratios.Add(r)
    pu = o["PUMP"].GetAsObject()
    enum(pu, "CalcMode", "OutletPressure")
    pu.Pout = 31e5
    pu.Efficiency = 75.0
    va = o["VALVE"].GetAsObject()
    enum(va, "CalcMode", "OutletPressure")
    va.OutletPressure = 5e5
    for d in ("DRUM", "DRUM2"):
        enum(o[d].GetAsObject(), "CalculationMode", "Adiabatic")
    if open_loop is None:
        r = o["RECYCLE"].GetAsObject()
        r.MaximumIterations = max_iter
        enum(r, "AccelerationMethod", "None")
        cp = r.ConvergenceParameters
        for f in TOL_FIELDS:
            setattr(cp, f, getattr(cp, f) * tol_scale)
    return sh, o


# DWSIM's streams by reflux's names. rec_b leaves the pump; in the closed loop
# the recycle block hands it on as rec_c.
STREAMS = {"feed": "feed", "mixed": "s1", "cold": "s2", "gas": "vapour", "liquid": "liquid", "back": "rec_a",
           "on": "liq_product", "recycle": "rec_b", "let down": "s3", "flash gas": "vapour2", "product": "liquid2"}
PRODUCTS = ["gas", "flash gas", "product"]


def run_reference(auto, models):
    """The reference flowsheet, closed. Its time isn't kept: the fixture is
    reproduced byte for byte, and step 9 times it."""
    sh, o = reference(auto)
    sh.solve()
    r = o["RECYCLE"].GetAsObject()
    if not r.Converged:
        raise RuntimeError(f"the recycle didn't converge in {r.IterationsTaken} iterations")
    streams = {name: models.enthalpies(sh.state(o[dw])) for name, dw in STREAMS.items()}
    energy = {"cooler": -o["COOL"].GetAsObject().DeltaQ * 1e3, "pump": o["PUMP"].GetAsObject().DeltaQ * 1e3,
              "drum": (o["DRUM"].GetAsObject().DeltaQ or 0.0) * 1e3,
              "second": (o["DRUM2"].GetAsObject().DeltaQ or 0.0) * 1e3}
    feed = streams["feed"]
    prods = [streams[n] for n in PRODUCTS]
    comp = max(abs(sum(p["flows"][i] for p in prods) - feed["flows"][i]) for i in range(len(COMPS))) / feed["f"]
    h_in = feed["f"] * feed["h"] + sum(energy.values())
    h_out = sum(p["f"] * p["h"] for p in prods)
    rt = feed["f"] * 8.314462618 * feed["t"]
    return {"iterations": r.IterationsTaken, "tolerance_scale": TIGHT, "max_iterations": 500,
            "streams": streams, "energy": energy,
            "balance": {"component": comp, "energy": abs(h_out - h_in) / rt,
                        "total": abs(sum(p["f"] for p in prods) - feed["f"]) / feed["f"]}}


def reflux_recycle():
    """reflux's converged recycle, from reflux solving the generated file."""
    with tempfile.TemporaryDirectory() as tmp:
        out = pathlib.Path(tmp) / "results.json"
        subprocess.run(["v", "run", ".", "solve", str(FLOWSHEET.relative_to(ROOT)), "--json", str(out)], cwd=ROOT,
                       check=True, capture_output=True)
        report = json.loads(out.read_text(encoding="utf-8"))
    s = next(s for s in report["streams"] if s["name"] == "recycle")
    return s


def run_open_loop(auto, fed):
    sh, o = reference(auto, open_loop=(fed["T"], fed["P"], fed["molar_flow"], fed["mole_fractions"]))
    sh.solve()
    back = sh.state(o["rec_b"])
    return {"fed": {"t": fed["T"], "p": fed["P"], "flows": fed["molar_flows"]}, "handed_back": back,
            "flows": max(abs(a - b) for a, b in zip(back["flows"], fed["molar_flows"])) / fed["molar_flow"],
            "t": back["t"] - fed["T"], "p": back["p"] - fed["P"]}


def write_components(constants):
    """DWSIM's constants as reflux component files, in data-format.md's format,
    with the ideal-gas Cp's range from the built-in file."""
    COMPONENTS_DIR.mkdir(parents=True, exist_ok=True)
    for c in constants:
        builtin = tomllib.loads((ROOT / "data" / "components" / f"{c['cas']}.toml").read_text(encoding="utf-8"))
        t_range = builtin["correlations"]["ideal_gas_cp"]["t_range"]
        text = "\n".join([
            f"# {c['dwsim_name']} ({c['cas']}), with DWSIM 9's constants for its Peng-Robinson.",
            "# Generated by tools/gen_m4_dwsim.py: don't edit.",
            "format = 1",
            f'name = "{c["dwsim_name"]}"',
            f'source = "DWSIM 9, {c["dwsim_name"]} from its ChemSep database (tools/gen_m4_dwsim.py)"',
            "",
            "[ids]",
            f'cas = "{c["cas"]}"',
            "",
            "[constants]",
            f"tc = {c['tc']!r}",
            f"pc = {c['pc']!r}",
            f"omega = {c['omega']!r}",
            f"mw = {c['mw']!r}",
            "",
            "[correlations.ideal_gas_cp]",
            "eq = 16",
            f"coeffs = [{', '.join(repr(x) for x in c['cp'])}]",
            f"t_range = [{', '.join(repr(float(x)) for x in t_range)}]",
            'units = "J/kmol/K"',
            "",
        ])
        (COMPONENTS_DIR / f"{c['cas']}.toml").write_text(text, encoding="utf-8", newline="\n")
    (COMPONENTS_DIR / "README.md").write_text(
        "# DWSIM 9's constants\n\n"
        "Written by [tools/gen_m4_dwsim.py](../../../tools/gen_m4_dwsim.py) from DWSIM 9's copy of the ChemSep\n"
        "database, for [reference_dwsim.toml](../../flowsheets/reference_dwsim.toml). They're ChemSep's data,\n"
        "under its Artistic License 2.0, as [data/components](../../../data/components/README.md) are.\n",
        encoding="utf-8", newline="\n")


def write_flowsheet(kij):
    """The reference flowsheet on DWSIM's constants and kij."""
    ref = REFERENCE.read_text(encoding="utf-8")
    pairs = []
    for i in range(len(NAMES)):
        for j in range(i + 1, len(NAMES)):
            assert kij[i][j] == kij[j][i], (i, j)
            if kij[i][j] != 0.0:
                pairs.append(f'{{ a = "{NAMES[i]}", b = "{NAMES[j]}", k = {kij[i][j]!r} }}')
    text = "\n".join([
        "# M4's reference flowsheet on DWSIM 9's constants and kij, to compare with",
        "# DWSIM (testdata/m4_dwsim.json). Generated by tools/gen_m4_dwsim.py: don't edit.",
        "format = 1",
        "",
        "[components]",
        f"names = [{', '.join(json.dumps(n) for n in NAMES)}]",
        'data = ["../m4_dwsim/components"]',
        "",
        "[model]",
        'equation = "PR"',
        "kij = false",
        f"kij_overrides = [{', '.join(pairs)}]",
        "",
        "",
    ]) + ref[ref.index("[solver]"):]
    FLOWSHEET.write_text(text, encoding="utf-8", newline="\n")


def main():
    gd.load_dwsim()
    from DWSIM.Automation import Automation3

    auto = Automation3()
    models = Models(auto)
    constants, kij, settings = Sheet(auto).constants()
    assert kij == models.dw.kij, "M3's DWSIM and this one have different kij"
    write_components(constants)
    write_flowsheet(kij)
    previous = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m4_dwsim.py",
            "python": platform.python_version(),
            "dwsim": str(auto.GetVersion()),
        },
        "note": "DWSIM 9's Peng-Robinson (PR) package at flash tolerances 1e-8. Streams: t K, p Pa, f and flows "
        "mol/s, h J/mol from the ideal gas at 298.15 K, s J/(mol K), v m3/mol from DWSIM's density, beta the "
        "vapour fraction, y and x the phases' mole fractions (empty when absent). energy is W into the material. "
        "spec is reflux's, as its file format has it. The reference flowsheet's recycle tolerances are DWSIM's "
        "scaled by tolerance_scale, with no acceleration; open_loop feeds the loop reflux's converged recycle.",
        "components": constants,
        "kij": kij,
        "flash_settings": settings,
        "units": previous.get("units", []),
        "reference": previous.get("reference", {}),
        "open_loop": previous.get("open_loop", {}),
    }
    if "units" in PARTS:
        data["units"] = []
        for case in UNIT_CASES:
            data["units"].append(run_unit(auto, models, case))
            print(f"  {case['name']}: done", flush=True)
    if "reference" in PARTS:
        data["reference"] = run_reference(auto, models)
        print(f"  the reference flowsheet: {data['reference']['iterations']} iterations", flush=True)
        data["open_loop"] = run_open_loop(auto, reflux_recycle())
        print("  the open loop: done", flush=True)
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    report(data)


def report(data):
    print(f"wrote {OUT.relative_to(ROOT)}, DWSIM {data['source']['dwsim']}")
    for u in data["units"]:
        b = u["balance"]
        print(f"  {u['name']}: energy {u['energy']:.6g} W, balances {b['component']:.2g} and {b['energy']:.2g}; "
              + "; ".join(f"out T {s['t']:.6f} K P {s['p'] / 1e5:.6f} bar beta {s['beta']:.6f}" for s in u["outlets"]))
    ref = data["reference"]
    if ref:
        b = ref["balance"]
        print(f"  reference: {ref['iterations']} iterations, balances {b['component']:.2g}, {b['energy']:.2g}, "
              f"total {b['total']:.2g}; energy {ref['energy']}")
    ol = data["open_loop"]
    if ol:
        print(f"  open loop: handed back within {ol['flows']:.2g} in flows, {ol['t']:.2g} K, {ol['p']:.2g} Pa")


if __name__ == "__main__":
    main()
