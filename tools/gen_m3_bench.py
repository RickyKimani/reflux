# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Write testdata/m3_bench.json: M3's benchmarks (decision D10, exit
criterion 8), reflux, thermo and DWSIM timed on the same flashes.
tools/benchmarks.vsh writes docs/benchmarks.md from it.

Run from the repository root, on Windows with DWSIM 9 installed:

    uv run --python 3.13 tools/gen_m3_bench.py

The flashes: PT with two phases and with one, PH, and the bubble point at a
given P, for a 2-, a 5- and a 10-component mixture, all Peng-Robinson with
DWSIM's constants and kij (read back from DWSIM, as gen_m3_dwsim.py does).
Each case is 40 nearby states, so that nothing is merely repeated: PT and PH
at T0 + 0.05 K * k, the bubble point at P0 * (1 + 1e-3 * k). PH is by state:
each program computes its own H at the PT state, untimed, and the PH flash
from it is timed.

reflux
    tools/bench_case.v, built here with gcc -prod and run once per case,
    at reflux's default tolerances.
thermo
    0.6.1's FlashVL at its default tolerances, with thermo's own ideal-gas
    Cp for each compound (which PT and bubble flashes don't use).
dwsim_default, dwsim_tight
    DWSIM 9's "Peng-Robinson (PR)" package, called directly
    (CalculateEquilibrium) through pythonnet, warm, at DWSIM's default PT
    and PH tolerances and at 1e-8.

Each flash is timed on its own (perf_counter_ns here, the monotonic clock
in reflux), after one untimed pass over the states; the median and the
minimum are in microseconds. The programs time each case one after the
other, within seconds: on a laptop the same flash's median moves by up to
3x from one session to another, so timed minutes apart the ratios moved by
2x. answer is
each program's result at the first state: beta for PT (the less dense
phase's fraction), T for PH and the bubble point, and for PH also the H
(J/mol) the program computed at the state. DWSIM's H is relative to the
ideal gas at 298.15 K, with the same Cp reflux is given; thermo's has its
own Cp and reference. Timings vary from run to
run, so unlike the fixture generators this one can't reproduce its output
byte for byte; its answers it can.
"""

import datetime
import json
import os
import pathlib
import platform
import statistics
import subprocess
import sys
import tempfile
import time

import thermo
from thermo import PRMIX, ChemicalConstantsPackage, CEOSGas, CEOSLiquid, FlashVL, HeatCapacityGas
from thermo import PropertyCorrelationsPackage

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader and its flowsheet per system
import gen_m3_fixtures as g3  # noqa: E402  by_density

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m3_bench.json"
STATES = 40
CALLS = {"thermo": 200, "dwsim_default": 200, "dwsim_tight": 200}
TIGHT = "1E-8"

SYSTEMS = [
    {"name": "methane/n-butane", "compounds": ["Methane", "N-butane"], "z": [0.4, 0.6],
     "pt2": (250.0, 30e5), "pt1": (350.0, 10e5), "bubble_p": 30e5},
    {"name": "natural gas", "compounds": ["Methane", "Ethane", "Propane", "N-butane", "N-pentane"],
     "z": [0.70, 0.10, 0.08, 0.07, 0.05], "pt2": (250.0, 40e5), "pt1": (350.0, 30e5), "bubble_p": 40e5},
    {"name": "gas condensate",
     "compounds": ["Nitrogen", "Carbon dioxide", "Methane", "Ethane", "Propane", "N-butane", "N-pentane",
                   "N-hexane", "N-heptane", "N-decane"],
     "z": [0.01, 0.02, 0.60, 0.10, 0.07, 0.05, 0.04, 0.04, 0.04, 0.03],
     "pt2": (280.0, 50e5), "pt1": (500.0, 20e5), "bubble_p": 30e5},
]
KINDS = ["pt2", "pt1", "ph", "bubble"]


def states(system, kind):
    """The case's states as (T, P) pairs; T is None for the bubble point."""
    if kind == "bubble":
        p0 = system["bubble_p"]
        return [(None, p0 * (1 + 1e-3 * k)) for k in range(STATES)]
    t0, p0 = system["pt1" if kind == "pt1" else "pt2"]
    return [(t0 + 0.05 * k, p0) for k in range(STATES)]


def timed(call, args, calls):
    """Median and minimum of calls timed calls, cycling through args after
    one untimed pass over them."""
    for a in args:
        call(a)
    ns = []
    for i in range(calls):
        a = args[i % len(args)]
        t0 = time.perf_counter_ns()
        call(a)
        ns.append(time.perf_counter_ns() - t0)
    return {"median_us": statistics.median(ns) / 1e3, "min_us": min(ns) / 1e3, "calls": calls}


def thermo_flasher(comps, kij):
    """FlashVL at its default tolerances, with DWSIM's constants and kij and
    thermo's own ideal-gas Cp."""
    tcs, pcs, ws = [c["tc"] for c in comps], [c["pc"] for c in comps], [c["omega"] for c in comps]
    constants = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=ws, MWs=[c["mw"] * 1000 for c in comps],
                                         CASs=[c["cas"] for c in comps])
    cps = [HeatCapacityGas(CASRN=c["cas"]) for c in comps]
    correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=cps, skip_missing=True)
    args = dict(Tcs=tcs, Pcs=pcs, omegas=ws, kijs=kij)
    return FlashVL(constants, correlations, liquid=CEOSLiquid(PRMIX, args, HeatCapacityGases=cps),
                   gas=CEOSGas(PRMIX, args, HeatCapacityGases=cps))


def bench_thermo(f, kind, sts, z):
    if kind in ("pt2", "pt1"):
        res = f.flash(T=sts[0][0], P=sts[0][1], zs=z)
        pairs = g3.by_density(res)
        beta = pairs[0][0] if len(pairs) == 2 else (1.0 if pairs[0][1].__class__.__name__ == "CEOSGas" else 0.0)
        return timed(lambda s: f.flash(T=s[0], P=s[1], zs=z), sts, CALLS["thermo"]), {"beta": beta}
    if kind == "ph":
        specs = [(p, f.flash(T=t, P=p, zs=z).H()) for t, p in sts]
        answer = {"t": f.flash(P=specs[0][0], H=specs[0][1], zs=z).T, "h": specs[0][1]}
        return timed(lambda s: f.flash(P=s[0], H=s[1], zs=z), specs, CALLS["thermo"]), answer
    answer = {"t": f.flash(P=sts[0][1], VF=0.0, zs=z).T}
    return timed(lambda s: f.flash(P=s[1], VF=0.0, zs=z), sts, CALLS["thermo"]), answer


def bench_dwsim(dw, kind, sts, z, calls):
    """Times CalculateEquilibrium itself; the answer comes from a separate
    call that also reads the result back."""
    zv = dw.vec(z)
    calc = dw.pp.CalculateEquilibrium
    if kind in ("pt2", "pt1"):
        fct = dw.FCT.PressureTemperature
        answer = {"beta": dw.flash("PressureTemperature", sts[0][1], sts[0][0], z)["beta"]}
        return timed(lambda s: calc(fct, s[1], s[0], zv, None, 0.0), sts, calls), answer
    if kind == "ph":
        fct = dw.FCT.PressureEnthalpy
        specs = [(p, dw.flash("PressureTemperature", p, t, z)["h_mass"]) for t, p in sts]
        # h_mass is in kJ/kg and mw in g/mol, so their product is J/mol.
        answer = {"t": dw.flash("PressureEnthalpy", specs[0][0], specs[0][1], z)["t"],
                  "h": specs[0][1] * dw.mw(z)}
        return timed(lambda s: calc(fct, s[0], s[1], zv, None, 0.0), specs, calls), answer
    fct = dw.FCT.PressureVaporFraction
    answer = {"t": dw.flash("PressureVaporFraction", sts[0][1], 0.0, z)["t"]}
    return timed(lambda s: calc(fct, s[1], 0.0, zv, None, 0.0), sts, calls), answer


def run(name, bench):
    """One program's result for one case, or the error it raised."""
    try:
        timing, answer = bench()
    except Exception as e:  # a failed flash is a result to record, not a reason to stop
        return {"error": f"{type(e).__name__}: {str(e).splitlines()[0] if str(e) else ''}"[:200]}
    print(f"    {name}: {timing['median_us']:.1f} us median, answer {answer}", flush=True)
    return {**timing, "answer": answer}


def build_reflux(tmp):
    """bench_case.exe, built with gcc -prod, and how it was built."""
    exe = pathlib.Path(tmp) / "bench_case.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_case.v"], cwd=ROOT,
                   check=True, env=env)
    build = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True, check=True).stdout)
    return exe, build


def bench_reflux(exe, tmp, system, case):
    """reflux's timing of one case, from bench_case.exe."""
    path = pathlib.Path(tmp) / "case.json"
    path.write_text(json.dumps({"system": system, "flash": case}), encoding="utf-8")
    out = subprocess.run([str(exe), str(path)], capture_output=True, text=True, check=True).stdout
    timing = json.loads(out)
    if timing.get("error"):
        return {"error": timing["error"]}
    print(f"    reflux: {timing['median_us']:.1f} us median, answer {timing['answer']}", flush=True)
    return timing


def cpu_name():
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                             capture_output=True, text=True, timeout=60).stdout
        return out.strip() or platform.processor()
    except (OSError, subprocess.SubprocessError):
        return platform.processor()


def main():
    gd.load_dwsim()
    from DWSIM.Automation import Automation3

    auto = Automation3()
    tmp = tempfile.mkdtemp()
    exe, build = build_reflux(tmp)
    tight = [gd.Dwsim(auto, s, TIGHT) for s in SYSTEMS]
    default = [gd.Dwsim(auto, s, None) for s in SYSTEMS]
    systems, cases = [], []
    for i, s in enumerate(SYSTEMS):
        dw = tight[i]
        n = len(s["compounds"])
        systems.append({"name": s["name"], "components": dw.comps, "kij": [x for row in dw.kij for x in row],
                        "z": s["z"]})
        f = thermo_flasher(dw.comps, dw.kij)
        for kind in KINDS:
            sts = states(s, kind)
            print(f"  {s['name']} ({n}), {kind}:", flush=True)
            case = {"system": i, "kind": kind, "t": [t for t, _ in sts] if kind != "bubble" else [],
                    "p": [p for _, p in sts]}
            case["reflux"] = bench_reflux(exe, tmp, systems[-1], case)
            case["thermo"] = run("thermo", lambda: bench_thermo(f, kind, sts, s["z"]))
            case["dwsim_default"] = run("dwsim_default",
                                        lambda: bench_dwsim(default[i], kind, sts, s["z"], CALLS["dwsim_default"]))
            case["dwsim_tight"] = run("dwsim_tight",
                                      lambda: bench_dwsim(tight[i], kind, sts, s["z"], CALLS["dwsim_tight"]))
            cases.append(case)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_bench.py",
            "python": platform.python_version(),
            "thermo": thermo.__version__,
            "dwsim": str(auto.GetVersion()),
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": build,
        "note": "M3's benchmarks: reflux's, thermo's and DWSIM's timings and answers, each case timed by the "
        "three programs one after the other. "
        "Peng-Robinson with DWSIM's constants and kij for all three. Each case is 40 states (t, p; no t for the "
        "bubble point), each flash timed on its own after a warm pass, median_us and min_us in microseconds. "
        "answer is at the first state: beta for PT (the less dense phase's fraction, 1 or 0 for one phase), T "
        "for PH (by state: from the program's own H at the PT state, h, in J/mol) and for the bubble point at P. "
        "dwsim_default is at DWSIM's default PT and PH tolerances, dwsim_tight at 1e-8.",
        "states": STATES,
        "dwsim_settings": {"default": default[0].settings, "tight": tight[0].settings},
        "systems": systems,
        "cases": cases,
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
