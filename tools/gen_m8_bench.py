# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""Write testdata/m8_bench.json: M8's speed (step 14, criterion 11), reflux
against DWSIM 9 in one session. tools/benchmarks.vsh adds it to
docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on
the PATH:

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_bench.py

About twenty minutes, nearly all of it DWSIM's study.

- The owner's study: the reforming PFR's feed from 700 to 1000 C in 50
  points. DWSIM by its own sensitivity analysis's form (step 0's spike runs
  it unshown), the owner's file's setup (dV 0.05, Radau); reflux on
  testdata/flowsheets/reformer_dwsim.toml, the same case on DWSIM's rates as
  typed, cold, warm (from each point's neighbour, as a study runs) and warm
  on 2 and 4 workers.
- One solve of the published reforming model, testdata/flowsheets/reformer.toml,
  whose K(T) step 1 keeps until T moves: reflux alone.
- An optimisation: the two-stage compression's least work, by DWSIM's
  optimiser's form with bounded Brent, the one method of its 11 that finds
  it (step 0), and by reflux, testdata/flowsheets/compression_dwsim.toml.
- A spec: the reference flowsheet's cooler outlet T until the first drum's
  vapour is 60 mol/s, by an Adjust's secant at its default tolerance, 1e-4,
  as its control panel calls it, and by reflux's spec on reference_dwsim.toml.

reflux is tools/bench_m8.v built with gcc -prod, each time the best of
three runs after an untimed one, all of it before DWSIM's runs and again
after: a laptop's clock moves DWSIM's solve between 9.4 and 29 s (M9 step 0),
so each comparison is with reflux's runs on either side. Timings vary from
run to run, so this can't reproduce its output byte for byte.
"""

import datetime
import importlib.util
import json
import os
import pathlib
import platform
import subprocess
import sys
import tempfile
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m8_bench.json"
sys.argv = sys.argv[:1]
sys.path.insert(0, str(ROOT / "tools"))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


study = module("dwsim_study", "spike/dwsim_study/run.py")
ref = study.ref
import gen_m3_bench as gb  # noqa: E402  the machine's name

REPEATS = 3
OWNER = "testdata/flowsheets/reformer_dwsim.toml"
STUDIES = [("cold", 1, 0), ("warm", 1, 1), ("warm, 2 workers", 2, 1), ("warm, 4 workers", 4, 1)]


def build(tmp):
    exe = pathlib.Path(tmp) / "bench_m8.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_m8.v"], cwd=ROOT,
                   check=True, env=env)
    info = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True, check=True).stdout)
    return exe, info


def bench(exe, *args):
    r = subprocess.run([str(exe), *map(str, args)], cwd=ROOT, capture_output=True, text=True)
    out = json.loads(r.stdout)
    if out["error"]:
        raise RuntimeError(f"bench_m8 {' '.join(map(str, args))}: {out['error']}")
    return out


def reflux_round(exe):
    """Everything reflux runs, once: each time the best of REPEATS."""
    r = {"study": {}}
    for label, workers, warm in STUDIES:
        out = bench(exe, "study", OWNER, workers, warm, REPEATS)
        r["study"][label] = {"ms": out["ms"], "methane": out["values"]}
        print(f"  reflux, the owner's study {label}: {min(out['ms']):.0f} ms", flush=True)
    out = bench(exe, "solve", "testdata/flowsheets/reformer.toml", 10)
    r["published"] = {"ms": out["ms"]}
    out = bench(exe, "optimize", "testdata/flowsheets/compression_dwsim.toml", "least_work", REPEATS)
    r["optimize"] = {"ms": out["ms"], "P": out["values"][0], "work": out["objective"],
                     "iterations": out["iterations"], "evaluations": out["evaluations"], "solves": out["solves"]}
    out = bench(exe, "spec", "testdata/flowsheets/reference_dwsim.toml", REPEATS)
    r["spec"] = {"ms": out["ms"], "T": out["values"][0], "passes": out["passes"]}
    print(f"  reflux: published model {min(r['published']['ms']):.1f} ms; optimisation "
          f"{min(r['optimize']['ms']):.0f} ms; spec {min(r['spec']['ms']):.0f} ms", flush=True)
    return r


def dwsim_spec(dw):
    """An Adjust's secant at its default tolerance, on the tight recycle."""
    from MathNet.Numerics.RootFinding import Secant
    from System import Double, Func

    sh, o = study.spec_sheet(dw, ref.TIGHT)
    cool, vap = o["COOL"].GetAsObject(), o["vapour"].GetAsObject()
    calls = [0]

    def f(t):
        calls[0] += 1
        cool.SetPropertyValue("PROP_CL_2", float(t))
        _, errs = sh.solve()
        if errs:
            raise RuntimeError(errs)
        return float(vap.GetPropertyValue("PROP_MS_3")) - study.VAPOUR_TARGET

    f(ref.COOL_T)
    calls[0] = 0
    t0 = time.perf_counter()
    x = Secant.FindRoot(Func[Double, Double](f), ref.COOL_T, ref.COOL_T * 1.01, study.SPEC_LO, study.SPEC_HI,
                        1e-4, 10)
    return {"s": time.perf_counter() - t0, "T": x, "solves": calls[0]}


def dwsim_optimize(dw):
    cp = study.Compression(dw)
    for c in ["C1", "C2"]:
        cp.o[c].GetAsObject().AdiabaticEfficiency = 75.0
    cp.work(study.MID_START)
    vs = [("x", "C1", "PROP_CO_4", "Pa", "Independent", study.MID_LO, study.MID_HI, study.MID_START),
          ("w1", "W1", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0),
          ("w2", "W2", "PROP_ES_0", "kW", "Auxiliary", 0.0, 0.0, 0.0)]
    objs, n, dt, why = study.optimise(dw, cp.sh.sim, "AL_BRENT_B", vs, "w1+w2", False)
    return {"s": dt, "P": float(objs["C1"].GetPropertyValue("PROP_CO_4")), "evaluations": n, "why": why}


def dwsim_study(dw):
    fs = study.pfr_classic(dw)
    study.solve_classic(dw, fs)
    rows, dt, log, finished = study.run_study(dw, fs, study.FEED_T, study.T_LO, study.T_HI, 50, [study.CH4],
                                              120 * 50)
    assert finished and len(rows) == 50, log[-500:]
    return {"s": dt, "rows": [{"T": r[0], "methane": r[1]} for r in rows]}


def main():
    t_start = time.perf_counter()
    tmp = tempfile.mkdtemp()
    exe, info = build(tmp)
    dw = study.Dwsim()
    print("reflux, before DWSIM", flush=True)
    before = reflux_round(exe)
    print("DWSIM: the spec, the optimisation, then the owner's study", flush=True)
    spec = dwsim_spec(dw)
    print(f"  DWSIM's spec: {spec['s']:.1f} s, {spec['solves']} solves", flush=True)
    opt = dwsim_optimize(dw)
    print(f"  DWSIM's optimisation: {opt['s']:.1f} s, {opt['evaluations']} evaluations", flush=True)
    owner = dwsim_study(dw)
    print(f"  DWSIM's study: {owner['s']:.0f} s", flush=True)
    print("reflux, after DWSIM", flush=True)
    after = reflux_round(exe)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m8_bench.py",
            "python": platform.python_version(),
            "dwsim": str(dw.a3.GetVersion()),
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": gb.cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": info,
        "note": "M8's speed: the owner's 50-point study of the reforming PFR by DWSIM 9's own sensitivity analysis "
                "and by reflux (gcc -prod) cold, warm and on workers; one solve of the published model; the "
                "two-stage compression's least work by DWSIM's optimiser (bounded Brent) and by reflux; the "
                "reference flowsheet's spec by an Adjust's secant at 1e-4 and by reflux. DWSIM's times are "
                "seconds, reflux's milliseconds, each of REPEATS runs after an untimed one, all of reflux's "
                "before DWSIM's runs and again after.",
        "reflux_before": before,
        "reflux_after": after,
        "dwsim": {"study": owner, "optimize": opt, "spec": spec},
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
