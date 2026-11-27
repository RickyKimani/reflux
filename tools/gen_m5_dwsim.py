# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5"]
# ///
"""Write M5a's DWSIM references (step 9, D14, criterion 8):
testdata/m5_dwsim.json.

DWSIM 9.0.5's two packages for pure water, by script, each at its defaults:
"CoolProp", which is IAPWS-95 by the CoolProp library, and "Steam Tables
(IAPWS-IF97)". For each:

- a PT flash at (T, P) states over IF97's five regions, near the critical
  point and in the compressed liquid: h and s per kg (DWSIM's kJ/kg and
  kJ/(kg K), on IAPWS's reference), the vapour fraction, or the error it
  raised;
- a PH flash back from each state's h, and at qualities inside the dome
  (h by the package's own saturated phases at 1 bar and 10 bar);
- its bubble temperature at five pressures;
- the Rankine cycle of testdata/flowsheets/rankine.toml as a DWSIM flowsheet
  (its Pump, Heater, Expander and Cooler), every unit checked calculated and
  without an error message (a failed unit is silent in DWSIM): each stream's
  T, P, mass flow, h, s and vapour fraction, and each unit's DeltaQ, kW.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_dwsim.py

On Windows with DWSIM 9 installed, in about a minute. Its numbers
reproduce byte for byte on the same day.
"""
import datetime
import json
import math
import os
import pathlib
import sys

import pythonnet

pythonnet.load("netfx")
import clr  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5_dwsim.json"
HOME = os.path.join(os.environ["LOCALAPPDATA"], "DWSIM")
sys.path.append(HOME)
os.chdir(HOME)  # DWSIM finds its databases relative to the working directory
for name in ["CapeOpen", "DWSIM.Automation", "DWSIM.Interfaces", "DWSIM.GlobalSettings", "DWSIM.SharedClasses",
             "DWSIM.Thermodynamics", "DWSIM.UnitOperations", "DWSIM.Inspector", "System.Buffers"]:
    clr.AddReference(name)

from System import Array, Double  # noqa: E402
from System.Diagnostics import FileVersionInfo  # noqa: E402
from DWSIM.Automation import Automation3  # noqa: E402
from DWSIM.Interfaces.Enums import FlashCalculationType  # noqa: E402
from DWSIM.Interfaces.Enums.GraphicObjects import ObjectType  # noqa: E402

vec = lambda xs: Array[Double](list(xs))  # noqa: E731
auto = Automation3()
PACKAGES = {"coolprop": "CoolProp", "steam_tables": "Steam Tables (IAPWS-IF97)"}

# (T K, P Pa): a grid over IF97's regions 1, 2 and 3, the compressed liquid
# to 100 MPa, and region 5.
STATES = [(t, p) for t in (280.0, 300.0, 350.0, 400.0, 500.0, 600.0, 620.0, 640.0, 650.0, 660.0, 700.0, 800.0,
                           1000.0)
          for p in (1e4, 1e5, 1e6, 1e7, 2.5e7, 5e7, 1e8)]
STATES += [(647.0, 2.2e7), (648.0, 2.21e7), (650.0, 2.55837018e7), (750.0, 7.83095639e7),
           (1500.0, 5e5), (1500.0, 3e7), (2000.0, 3e7)]
BUBBLE_P = [1e4, 1e5, 1e6, 1e7, 2e7]
DOME_P = [1e5, 1e6]
DOME_X = [0.1, 0.5, 0.9]


def short(e):
    return str(e).splitlines()[0][:200]


def package(sim, name):
    pp = sim.CreateAndAddPropertyPackage(name).__implementation__
    return pp


def flash(pp, kind, a, b):
    """A flash: (T, P, h, s, beta), or the error DWSIM raised."""
    spec = getattr(FlashCalculationType, kind)
    r = pp.CalculateEquilibrium(spec, a, b, vec([1.0]), None, 0.0)
    return {"T": r.CalculatedTemperature, "P": r.CalculatedPressure, "h": r.CalculatedEnthalpy,
            "s": r.CalculatedEntropy, "beta": r.GetVaporPhaseMoleFraction()}


def guarded(f):
    try:
        out = f()
    except Exception as e:  # noqa: BLE001
        return {"error": short(e)}
    if any(v is None or (isinstance(v, float) and not math.isfinite(v)) for v in out.values()):
        return {"error": f"not finite: {out}"}
    return out


def states(pp):
    rows = []
    for t, p in STATES:
        r = guarded(lambda: flash(pp, "PressureTemperature", p, t))
        row = {"T": t, "P": p, **{k: r[k] for k in ("h", "s", "beta") if k in r}}
        if "error" in r:
            row["error"] = r["error"]
        else:
            back = guarded(lambda: flash(pp, "PressureEnthalpy", p, r["h"]))
            row["ph_T"] = back.get("T")
            if "error" in back:
                row["ph_error"] = back["error"]
        rows.append(row)
    return rows


def dome(pp):
    rows = []
    for p in DOME_P:
        liq = guarded(lambda: flash(pp, "PressureVaporFraction", p, 0.0))
        vap = guarded(lambda: flash(pp, "PressureVaporFraction", p, 1.0))
        if "error" in liq or "error" in vap:
            rows.append({"P": p, "error": liq.get("error") or vap.get("error")})
            continue
        for x in DOME_X:
            h = liq["h"] + x * (vap["h"] - liq["h"])
            r = guarded(lambda: flash(pp, "PressureEnthalpy", p, h))
            rows.append({"P": p, "x": x, "h": h, "h_liquid": liq["h"], "h_vapour": vap["h"], **r})
    return rows


def bubble(pp):
    rows = []
    for p in BUBBLE_P:
        r = guarded(lambda: flash(pp, "PressureVaporFraction", p, 0.0))
        rows.append({"P": p, **({"error": r["error"]} if "error" in r else {"T": r["T"]})})
    return rows


def rankine(pp_name):
    """testdata/flowsheets/rankine.toml in DWSIM: 50 mol/s of condensate at
    318.15 K and 10 kPa, pumped to 8 MPa (75%), heated to 773.15 K,
    expanded to 10 kPa (85%, adiabatic) and cooled to 318.15 K."""
    from DWSIM.UnitOperations.UnitOperations import Pump, Heater, Expander, Cooler  # noqa: F401

    sim = auto.CreateFlowsheet()
    sim.SelectedCompounds.Add("Water", sim.AvailableCompounds["Water"])
    pp = package(sim, pp_name)
    x = [0]
    units = []

    def add(kind, name):
        x[0] += 100
        obj = sim.AddObject(getattr(ObjectType, kind), x[0], 100, name)
        if kind not in ("MaterialStream", "EnergyStream"):
            units.append(obj)
        return obj

    def link(a, b):
        sim.ConnectObjects(a.GraphicObject, b.GraphicObject, -1, -1)

    def enum(o, prop, value):
        setattr(o, prop, getattr(getattr(o, prop).__class__, value))

    feed = add("MaterialStream", "condensate")
    f = feed.GetAsObject()
    f.SetTemperature(318.15)
    f.SetPressure(1e4)
    f.SetMolarFlow(50.0)
    f.SetOverallComposition(vec([1.0]))
    pp.CurrentMaterialStream = f
    names = ["condensate", "feedwater", "steam", "exhaust", "return"]
    streams = [feed] + [add("MaterialStream", n) for n in names[1:]]
    kinds = [("Pump", "pump"), ("Heater", "boiler"), ("Expander", "turbine"), ("Cooler", "condenser")]
    made = {}
    for i, (kind, name) in enumerate(kinds):
        u = add(kind, name)
        made[name] = u
        link(streams[i], u)
        link(u, streams[i + 1])
        q = add("EnergyStream", "Q" + name)
        if kind in ("Pump", "Heater"):
            link(q, u)
        else:
            link(u, q)
    pump = made["pump"].GetAsObject()
    enum(pump, "CalcMode", "OutletPressure")
    pump.Pout = 8e6
    pump.Efficiency = 75.0
    boiler = made["boiler"].GetAsObject()
    enum(boiler, "CalcMode", "OutletTemperature")
    boiler.OutletTemperature = 773.15
    boiler.DeltaP = 0.0
    turbine = made["turbine"].GetAsObject()
    enum(turbine, "CalcMode", "OutletPressure")
    enum(turbine, "ProcessPath", "Adiabatic")
    turbine.POut = 1e4
    turbine.AdiabaticEfficiency = 85.0
    condenser = made["condenser"].GetAsObject()
    enum(condenser, "CalcMode", "OutletTemperature")
    condenser.OutletTemperature = 318.15
    condenser.DeltaP = 0.0
    errs = auto.CalculateFlowsheet2(sim)
    bad = [short(e) for e in errs] if errs else []
    for u in units:
        o = u.GetAsObject()
        if not o.Calculated or (o.ErrorMessage or ""):
            bad.append(f"{u.GraphicObject.Tag}: calculated {o.Calculated}, {o.ErrorMessage!r}")
    if bad:
        return {"error": bad}
    out = {"streams": {}, "units": {}}
    for name, s in zip(names, streams):
        m = s.GetAsObject()
        ph = m.Phases[0].Properties
        out["streams"][name] = {"T": ph.temperature, "P": ph.pressure, "mass_flow": ph.massflow,
                                "molar_flow": ph.molarflow, "h": ph.enthalpy, "s": ph.entropy,
                                "beta": m.Phases[2].Properties.molarfraction or 0.0}
    for name in ("pump", "boiler", "turbine", "condenser"):
        out["units"][name] = made[name].GetAsObject().DeltaQ
    return out


def main():
    sim = auto.CreateFlowsheet()
    sim.SelectedCompounds.Add("Water", sim.AvailableCompounds["Water"])
    stream = sim.AddObject(ObjectType.MaterialStream, 50, 50, "w").GetAsObject()
    stream.SetTemperature(300.0)
    stream.SetPressure(1e5)
    stream.SetMolarFlow(1.0)
    stream.SetOverallComposition(vec([1.0]))
    out = {
        "provenance": {
            "script": "tools/gen_m5_dwsim.py",
            "dwsim": str(auto.GetVersion()),
            "coolprop_dll": FileVersionInfo.GetVersionInfo(os.path.join(HOME, "CoolProp.dll")).FileVersion,
            "date": datetime.date.today().isoformat(),
            "water_molar_mass": sim.SelectedCompounds["Water"].Molar_Weight,
            "units": "K, Pa, kJ/kg, kJ/(kg K), kg/s, mol/s; a unit's DeltaQ in kW",
        },
    }
    for key, name in PACKAGES.items():
        pp = package(sim, name)
        pp.CurrentMaterialStream = stream
        out[key] = {"package": name, "states": states(pp), "dome": dome(pp), "bubble": bubble(pp),
                    "rankine": rankine(name)}
    OUT.write_bytes((json.dumps(out, indent=1) + "\n").encode("utf-8"))
    for key in PACKAGES:
        bad = [r for r in out[key]["states"] if "error" in r]
        print(f"{key}: {len(out[key]['states'])} states, {len(bad)} refused; rankine "
              f"{'ok' if 'error' not in out[key]['rankine'] else out[key]['rankine']['error']}")


if __name__ == "__main__":
    main()
