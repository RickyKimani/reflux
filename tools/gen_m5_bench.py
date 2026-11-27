# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "CoolProp==8.0.0", "iapws==1.5.4", "numpy==2.2.6", "scipy==1.18.1"]
# ///
"""Write testdata/m5_bench.json: M5a's speed (step 10, criterion 10), water
by IAPWS-95 and IF97 in reflux, DWSIM 9, CoolProp and the iapws package in
one session. tools/benchmarks.vsh adds it to docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on
the PATH:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_bench.py

A few minutes. The states are 300 K and 3 MPa (a liquid), 500 K and 1 MPa
(a vapour) and 700 K and 30 MPa (dense, above Tc and pc). A PT call gives
h; a PH call gives T back from the h each program's own PT call gave. A
round is a case's three calls; each program runs one untimed round first.

- reflux: tools/bench_m5.v built with gcc -prod, 200 rounds, before the
  others' runs and again after: IAPWS-95's PT and PH flashes (flash.pt and
  flash.ph) and IF97's state_pt and state_ph.
- DWSIM: its "CoolProp" package (IAPWS-95) and its "Steam Tables
  (IAPWS-IF97)", CalculateEquilibrium by PT and by PH, 20 rounds.
- CoolProp: an AbstractState, "HEOS" (IAPWS-95) and "IF97", update by
  PT_INPUTS and HmassP_INPUTS, 200 rounds.
- iapws: IAPWS95(T=, P=) and IAPWS95(P=, h=), 20 rounds; IAPWS97 the same.

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
OUT = ROOT / "testdata" / "m5_bench.json"
sys.path.insert(0, str(ROOT / "tools"))
cwd = os.getcwd()
import gen_m5_dwsim as gd  # noqa: E402  loads DWSIM, and changes to its directory

os.chdir(cwd)
import CoolProp  # noqa: E402
import CoolProp.CoolProp as CP  # noqa: E402
import iapws  # noqa: E402
from iapws import IAPWS95, IAPWS97  # noqa: E402

STATES = [(300.0, 3e6), (500.0, 1e6), (700.0, 3e7)]  # T K, P Pa
REFLUX_ROUNDS, FAST_ROUNDS, OTHER_ROUNDS = 200, 200, 20


def cpu_name():
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                             capture_output=True, text=True, timeout=60).stdout
        return out.strip() or platform.processor()
    except (OSError, subprocess.SubprocessError):
        return platform.processor()


def build(tmp):
    exe = pathlib.Path(tmp) / "bench_m5.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_m5.v"], cwd=ROOT,
                   check=True, env=env)
    info = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True, check=True).stdout)
    return exe, info


def reflux_round(exe):
    out = subprocess.run([str(exe), str(REFLUX_ROUNDS)], cwd=ROOT, capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


def timed(model, kind, calls, rounds):
    """A case's answers, then each round's microseconds per call; or the
    error its first call raised."""
    try:
        answers = [float(c()) for c in calls]
    except Exception as e:  # noqa: BLE001  a failed call is a result to record
        return {"model": model, "kind": kind, "error": f"{type(e).__name__}: {str(e).splitlines()[0][:200]}"}
    us = []
    for _ in range(rounds):
        t0 = time.perf_counter_ns()
        for c in calls:
            c()
        us.append((time.perf_counter_ns() - t0) / 1e3 / len(calls))
    return {"model": model, "kind": kind, "us": us, "answers": answers}


def pair(model, pt_call, ph_call, rounds):
    """A model's PT case (h, kJ/kg) and its PH case back from its own h."""
    pt = timed(model, "pt", [lambda t=t, p=p: pt_call(t, p) for t, p in STATES], rounds)
    if "error" in pt:
        return [pt, {"model": model, "kind": "ph", "error": "no h: " + pt["error"]}]
    ph = timed(model, "ph", [lambda p=p, h=h: ph_call(p, h) for (_, p), h in zip(STATES, pt["answers"])], rounds)
    return [pt, ph]


def dwsim_cases():
    sim = gd.auto.CreateFlowsheet()
    sim.SelectedCompounds.Add("Water", sim.AvailableCompounds["Water"])
    stream = sim.AddObject(gd.ObjectType.MaterialStream, 50, 50, "w").GetAsObject()
    stream.SetTemperature(300.0)
    stream.SetPressure(1e5)
    stream.SetMolarFlow(1.0)
    stream.SetOverallComposition(gd.vec([1.0]))
    out = []
    for model, name in [("IAPWS-95", "CoolProp"), ("IF97", "Steam Tables (IAPWS-IF97)")]:
        pp = gd.package(sim, name)
        pp.CurrentMaterialStream = stream
        out += pair(model, lambda t, p, pp=pp: gd.flash(pp, "PressureTemperature", p, t)["h"],
                    lambda p, h, pp=pp: gd.flash(pp, "PressureEnthalpy", p, h)["T"], OTHER_ROUNDS)
        print("DWSIM", model, [(c["kind"], c.get("answers", c.get("error"))) for c in out[-2:]], flush=True)
    return out


def coolprop_cases():
    out = []
    for model, backend in [("IAPWS-95", "HEOS"), ("IF97", "IF97")]:
        st = CP.AbstractState(backend, "Water")

        def pt(t, p, st=st):
            st.update(CP.PT_INPUTS, p, t)
            return st.hmass() * 1e-3

        def ph(p, h, st=st):
            st.update(CP.HmassP_INPUTS, h * 1e3, p)
            return st.T()

        out += pair(model, pt, ph, FAST_ROUNDS)
        print("CoolProp", model, [(c["kind"], c.get("answers", c.get("error"))) for c in out[-2:]], flush=True)
    return out


def iapws_cases():
    out = []
    for model, cls in [("IAPWS-95", IAPWS95), ("IF97", IAPWS97)]:
        out += pair(model, lambda t, p, cls=cls: cls(T=t, P=p / 1e6).h,
                    lambda p, h, cls=cls: cls(P=p / 1e6, h=h).T, OTHER_ROUNDS)
        print("iapws", model, [(c["kind"], c.get("answers", c.get("error"))) for c in out[-2:]], flush=True)
    return out


def main():
    warnings.simplefilter("ignore")
    t_start = time.perf_counter()
    tmp = tempfile.mkdtemp()
    exe, info = build(tmp)
    print("reflux, before", flush=True)
    before = reflux_round(exe)
    dwsim = dwsim_cases()
    coolprop = coolprop_cases()
    iapws_ = iapws_cases()
    print("reflux, after", flush=True)
    after = reflux_round(exe)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m5_bench.py",
            "python": platform.python_version(),
            "dwsim": str(gd.auto.GetVersion()),
            "coolprop": CoolProp.__version__,
            "iapws": iapws.__version__,
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": info,
        "note": ("M5a's speed: water at 300 K and 3 MPa, 500 K and 1 MPa, and 700 K and 30 MPa, by IAPWS-95 and "
                 "IF97, a PT call (h, kJ/kg on IAPWS's reference, or CoolProp's own) and a PH call back (T, K), in "
                 "reflux (gcc -prod), DWSIM 9's two packages, CoolProp's AbstractState and the iapws package; each "
                 "round's microseconds per call, and each call's answer."),
        "reflux_before": before,
        "reflux_after": after,
        "dwsim": dwsim,
        "coolprop": coolprop,
        "iapws": iapws_,
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}", flush=True)
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
