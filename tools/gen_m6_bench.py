# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "thermo==0.6.1", "chemicals==1.5.2", "mpmath==1.3.0", "scipy==1.14.1"]
# ///
"""Write testdata/m6_bench.json: M6's speed (step 13, D18, criterion 12),
gamma-phi flashes by reflux, DWSIM 9 and thermo in one session, and the
cubic's beside them. tools/benchmarks.vsh adds it to docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on
the PATH:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_bench.py

A few minutes. Every case is methanol and water at 1 atm: PT flashes at
350 K (z = 0.5), 357 K (0.3) and 346 K (0.7), two-phase by NRTL, and the
bubble temperatures of liquids of 0.1, 0.5 and 0.9 methanol; by NRTL and
UNIFAC as gamma-phi fluids (the vapour the ideal gas, Poynting on), and by
Peng-Robinson. A round is a case's three calls; each program runs one
untimed round first.

- reflux: tools/bench_m6.v built with gcc -prod, 200 rounds, before
  DWSIM's and thermo's runs and again after.
- DWSIM: its NRTL, UNIFAC and Peng-Robinson packages at their defaults,
  CalculateEquilibrium by PT and by P and a vapour fraction of 0, 20 rounds.
- thermo: FlashVL on GibbsExcessLiquid ('Poynting', the ideal gas) with
  NRTL and UNIFAC on ChemSep's parameters (tools/gen_m6_gamma_phi.py's), and
  on PRMIX with pr.ipd's kij, -0.0778, as reflux's and DWSIM's take, 20
  rounds; its flash at a vapour fraction of 0, or its error where it has
  none.

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

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m6_bench.json"
sys.path.insert(0, str(ROOT / "tools"))
import gen_m3_bench as gb  # noqa: E402  the machine's name
import gen_m6_activity as ga  # noqa: E402
import gen_m6_gamma_phi as gp  # noqa: E402

os_cwd = os.getcwd()
import gen_m6_dwsim as gd  # noqa: E402  loads DWSIM, and changes to its directory

os.chdir(os_cwd)
from thermo import (ChemicalConstantsPackage, PropertyCorrelationsPackage, FlashVL, CEOSGas,  # noqa: E402
                    CEOSLiquid, PRMIX, HeatCapacityGas)

P = 101325.0
PT = [(350.0, 0.5), (357.0, 0.3), (346.0, 0.7)]
BUBBLE = [0.1, 0.5, 0.9]
CAS = ["67-56-1", "7732-18-5"]
REFLUX_ROUNDS, OTHER_ROUNDS = 200, 20


def build(tmp):
    exe = pathlib.Path(tmp) / "bench_m6.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_m6.v"], cwd=ROOT,
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
        answers = [c() for c in calls]
    except Exception as e:  # noqa: BLE001  a failed flash is a result to record
        return {"model": model, "kind": kind, "error": f"{type(e).__name__}: {str(e).splitlines()[0][:200]}"}
    us = []
    for _ in range(rounds):
        t0 = time.perf_counter_ns()
        for c in calls:
            c()
        us.append((time.perf_counter_ns() - t0) / 1e3 / len(calls))
    return {"model": model, "kind": kind, "us": us, "answers": answers}


def dwsim_cases():
    out = []
    names = ["Methanol", "Water"]
    for model, package in [("NRTL", "NRTL"), ("UNIFAC", "UNIFAC"), ("PR", "Peng-Robinson (PR)")]:
        pk = gd.Package(names, package, [0.5, 0.5], tol=None)
        calc, fct, vec = pk.pp.CalculateEquilibrium, gd.FlashCalculationType, gd.vec
        pt = [lambda t=t, z=z: calc(fct.PressureTemperature, P, t, vec([z, 1 - z]), None, 0.0)
              .GetVaporPhaseMoleFraction() for t, z in PT]
        bub = [lambda x=x: calc(fct.PressureVaporFraction, P, 0.0, vec([x, 1 - x]), None, 0.0).CalculatedTemperature
               for x in BUBBLE]
        out.append(timed(model, "pt", pt, OTHER_ROUNDS))
        out.append(timed(model, "bubble", bub, OTHER_ROUNDS))
        print("DWSIM", model, [(c["kind"], c.get("answers", c.get("error"))) for c in out[-2:]], flush=True)
    return out


def thermo_pr():
    recs = [ga.record(c) for c in CAS]
    tcs = [r["constants"]["tc"] for r in recs]
    pcs = [r["constants"]["pc"] for r in recs]
    ws = [r["constants"]["omega"] for r in recs]
    consts = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=ws, MWs=[r["constants"]["mw"] * 1e3 for r in recs],
                                      CASs=CAS)
    cps = [gp.heat_capacity(r) for r in recs]
    corr = PropertyCorrelationsPackage(constants=consts, HeatCapacityGases=cps, skip_missing=True)
    k = -0.0778  # pr.ipd: Methanol/Water p730
    kw = dict(Tcs=tcs, Pcs=pcs, omegas=ws, kijs=[[0.0, k], [k, 0.0]])
    return FlashVL(consts, corr, liquid=CEOSLiquid(PRMIX, kw, HeatCapacityGases=cps),
                   gas=CEOSGas(PRMIX, kw, HeatCapacityGases=cps))


def thermo_cases():
    out = []
    for model in ["NRTL", "UNIFAC", "PR"]:
        f = thermo_pr() if model == "PR" else gp.parts(model, CAS)[0]
        pt = [lambda t=t, z=z: f.flash(T=t, P=P, zs=[z, 1 - z]).VF for t, z in PT]
        bub = [lambda x=x: f.flash(P=P, VF=0.0, zs=[x, 1 - x]).T for x in BUBBLE]
        out.append(timed(model, "pt", pt, OTHER_ROUNDS))
        out.append(timed(model, "bubble", bub, OTHER_ROUNDS))
        print("thermo", model, [(c["kind"], c.get("answers", c.get("error"))) for c in out[-2:]], flush=True)
    return out


def main():
    t_start = time.perf_counter()
    tmp = tempfile.mkdtemp()
    exe, info = build(tmp)
    print("reflux, before", flush=True)
    before = reflux_round(exe)
    dwsim = dwsim_cases()
    thermo = thermo_cases()
    print("reflux, after", flush=True)
    after = reflux_round(exe)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m6_bench.py",
            "python": platform.python_version(),
            "thermo": __import__("importlib.metadata").metadata.version("thermo"),
            "dwsim": gd.VERSION,
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": gb.cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": info,
        "note": ("M6's speed: methanol-water at 1 atm, PT flashes at 350, 357 and 346 K and bubble temperatures of "
                 "0.1, 0.5 and 0.9 methanol, by NRTL and UNIFAC as gamma-phi fluids and by PR, in reflux (gcc "
                 "-prod), DWSIM 9 at its defaults and thermo's FlashVL; each round's microseconds per call, and "
                 "each call's answer, the vapour fraction or the bubble T."),
        "reflux_before": before,
        "reflux_after": after,
        "dwsim": dwsim,
        "thermo": thermo,
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}", flush=True)
    sys.stdout.flush()
    os._exit(0)


if __name__ == "__main__":
    main()
