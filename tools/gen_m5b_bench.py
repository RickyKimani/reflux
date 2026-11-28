# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""Write testdata/m5b_bench.json: M5b's speed (step 11, criterion 10), a
phase's transport in reflux, DWSIM 9 and thermo, and IF97 behind reflux's
flash beside IAPWS-95's, water.if97's functions and DWSIM's Steam Tables
package, in one session. tools/benchmarks.vsh adds it to docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on
the PATH:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5b_bench.py

A few minutes.

- reflux: tools/bench_m5b.v built with gcc -prod, before the others' runs
  and again after: a phase's transport (transport.Mixture.phase, all its
  properties at once, 2,000 rounds) at three states, benzene-toluene's
  liquid at 300 K and 1 atm, the reforming gas at 1000 K and 30 bar and
  methane-propane at 320 K and 100 bar; IF97's and IAPWS-95's PT and PH
  flashes and water.if97's state_pt and state_ph at M5a's three states
  (200 rounds); and a solve's report with each phase's transport and
  without (200 rounds).
- DWSIM: its "Peng-Robinson (PR)" package's transport functions on a
  stream calculated at each state, the liquid's AUX_LIQVISCm, AUX_CONDTL
  and AUX_SURFTM and the vapour's AUX_VAPVISCm and AUX_CONDTG, 200 rounds;
  and its "Steam Tables (IAPWS-IF97)" package's CalculateEquilibrium by PT
  and by PH, 20 rounds.
- thermo: its mixture properties at their default methods, from
  ChemicalConstantsPackage.from_IDs (its own pure data, not reflux's), each
  property's mixture_property, 200 rounds.

A round is a case's calls; each program runs one untimed round first.
Timings vary from run to run, so this can't reproduce its output byte for
byte.
"""
import datetime
import json
import os
import pathlib
import platform
import subprocess
import sys
import tempfile
import time
import warnings

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m5b_bench.json"
sys.path.insert(0, str(ROOT / "tools"))
cwd = os.getcwd()
import gen_m5_dwsim as gd  # noqa: E402  loads DWSIM, and changes to its directory

os.chdir(cwd)
import thermo  # noqa: E402
from thermo import ChemicalConstantsPackage  # noqa: E402

REFLUX_ROUNDS, DWSIM_ROUNDS, THERMO_ROUNDS, FLASH_ROUNDS = 200, 200, 200, 20
# Each case: its name, compounds by DWSIM's names and by thermo's, mole
# fractions, T K, P Pa, and the phase.
CASES = [
    ("benzene-toluene liquid", ["Benzene", "Toluene"], ["benzene", "toluene"], [0.4, 0.6], 300.0, 101325.0,
     "liquid"),
    ("reforming gas", ["Methane", "Hydrogen", "Water", "Carbon dioxide", "Carbon monoxide"],
     ["methane", "hydrogen", "water", "carbon dioxide", "carbon monoxide"], [0.15, 0.5, 0.2, 0.05, 0.1], 1000.0,
     3e6, "vapour"),
    ("methane-propane dense gas", ["Methane", "Propane"], ["methane", "propane"], [0.8, 0.2], 320.0, 1e7,
     "vapour"),
]
WATER_STATES = [(300.0, 3e6), (500.0, 1e6), (700.0, 3e7)]  # T K, P Pa


def cpu_name():
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                             capture_output=True, text=True, timeout=60).stdout
        return out.strip() or platform.processor()
    except (OSError, subprocess.SubprocessError):
        return platform.processor()


def build(tmp):
    exe = pathlib.Path(tmp) / "bench_m5b.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_m5b.v"], cwd=ROOT,
                   check=True, env=env)
    info = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True, check=True).stdout)
    return exe, info


def reflux_run(exe):
    out = subprocess.run([str(exe), str(REFLUX_ROUNDS)], cwd=ROOT, capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


def timed(name, kind, call, rounds):
    """A case's answers, then each round's microseconds; or the error its
    first call raised."""
    try:
        answers = [float(v) for v in call()]
    except Exception as e:  # noqa: BLE001  a failed call is a result to record
        return {"name": name, "kind": kind, "error": f"{type(e).__name__}: {str(e).splitlines()[0][:200]}"}
    us = []
    for _ in range(rounds):
        t0 = time.perf_counter_ns()
        call()
        us.append((time.perf_counter_ns() - t0) / 1e3)
    return {"name": name, "kind": kind, "us": us, "answers": answers}


def dwsim_transport():
    out = []
    for name, compounds, _, z, t, p, phase in CASES:
        sim = gd.auto.CreateFlowsheet()
        for c in compounds:
            sim.SelectedCompounds.Add(c, sim.AvailableCompounds[c])
        pp = gd.package(sim, "Peng-Robinson (PR)")
        ms = sim.AddObject(gd.ObjectType.MaterialStream, 50, 50, "s").GetAsObject()
        ms.SetPropertyPackage(pp)
        ms.SetTemperature(t)
        ms.SetPressure(p)
        ms.SetMolarFlow(1.0)
        ms.SetOverallComposition(gd.vec(z))
        pp.CurrentMaterialStream = ms
        ms.Calculate(True, True)
        pp.CurrentMaterialStream = ms
        if phase == "liquid":
            def call(pp=pp, t=t, p=p):
                return [pp.AUX_LIQVISCm(t, p, 3), pp.AUX_CONDTL(t, 3), pp.AUX_SURFTM(t)]
        else:
            q = ms.Phases[2].Properties
            rho, mw = float(q.density), float(q.molecularWeight)

            def call(pp=pp, t=t, p=p, rho=rho, mw=mw):
                return [pp.AUX_VAPVISCm(t, rho, mw), pp.AUX_CONDTG(t, p)]
        out.append(timed(name, "phase", call, DWSIM_ROUNDS))
        print("DWSIM", name, out[-1].get("answers", out[-1].get("error")), flush=True)
    return out


def dwsim_water():
    sim = gd.auto.CreateFlowsheet()
    sim.SelectedCompounds.Add("Water", sim.AvailableCompounds["Water"])
    stream = sim.AddObject(gd.ObjectType.MaterialStream, 50, 50, "w").GetAsObject()
    stream.SetTemperature(300.0)
    stream.SetPressure(1e5)
    stream.SetMolarFlow(1.0)
    stream.SetOverallComposition(gd.vec([1.0]))
    pp = gd.package(sim, "Steam Tables (IAPWS-IF97)")
    pp.CurrentMaterialStream = stream
    pt = timed("Steam Tables", "pt",
               lambda: [gd.flash(pp, "PressureTemperature", p, t)["h"] for t, p in WATER_STATES], FLASH_ROUNDS)
    hs = pt.get("answers", [])
    ph = timed("Steam Tables", "ph",
               lambda: [gd.flash(pp, "PressureEnthalpy", p, h)["T"] for (_, p), h in zip(WATER_STATES, hs)],
               FLASH_ROUNDS)
    for c in (pt, ph):  # per call, as reflux's are
        if "us" in c:
            c["us"] = [u / len(WATER_STATES) for u in c["us"]]
    print("DWSIM Steam Tables", pt.get("answers"), ph.get("answers"), flush=True)
    return [pt, ph]


def thermo_transport():
    out = []
    for name, _, ids, z, t, p, phase in CASES:
        constants, correlations = ChemicalConstantsPackage.from_IDs(ids)
        mws = constants.MWs
        mass = [zi * m for zi, m in zip(z, mws)]
        ws = [m / sum(mass) for m in mass]
        if phase == "liquid":
            props = [correlations.ViscosityLiquidMixture, correlations.ThermalConductivityLiquidMixture,
                     correlations.SurfaceTensionMixture]
        else:
            props = [correlations.ViscosityGasMixture, correlations.ThermalConductivityGasMixture]

        def call(props=props, t=t, p=p, z=z, ws=ws):
            return [m.mixture_property(t, p, z, ws) for m in props]

        out.append(timed(name, "phase", call, THERMO_ROUNDS))
        print("thermo", name, out[-1].get("answers", out[-1].get("error")), flush=True)
    return out


def main():
    warnings.simplefilter("ignore")
    t_start = time.perf_counter()
    tmp = tempfile.mkdtemp()
    exe, info = build(tmp)
    print("reflux, before", flush=True)
    before = reflux_run(exe)
    dwsim = dwsim_transport()
    steam = dwsim_water()
    th = thermo_transport()
    print("reflux, after", flush=True)
    after = reflux_run(exe)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m5b_bench.py",
            "python": platform.python_version(),
            "dwsim": str(gd.auto.GetVersion()),
            "thermo": thermo.__version__,
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": info,
        "note": ("M5b's speed: a phase's transport at three states in reflux (gcc -prod), DWSIM 9's PR package's "
                 "functions and thermo's mixture properties; IF97 and IAPWS-95 behind reflux's flash, IF97's "
                 "functions and DWSIM's Steam Tables; a solve's report with and without transport. Each round's "
                 "microseconds, per call for the water cases, and each case's answers."),
        "reflux_before": before,
        "reflux_after": after,
        "dwsim": dwsim,
        "dwsim_steam_tables": steam,
        "thermo": th,
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}", flush=True)
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
