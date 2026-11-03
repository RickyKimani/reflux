# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Write testdata/m4_bench.json: M4's flowsheet benchmark (decision D14, exit
criterion 9), the reference flowsheet on DWSIM's constants solved by reflux and
by DWSIM 9. tools/benchmarks.vsh adds it to docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on the
PATH:

    uv run --python 3.13 tools/gen_m4_bench.py

DWSIM builds the flowsheet as tools/gen_m4_dwsim.py does, at its default recycle
tolerances and at them scaled by 1e-8. After one untimed solve, so that .NET has
compiled it, it times first solves, each of a flowsheet just built, then solves
again, each with the feed 0.01 K warmer, which DWSIM starts from its last
answer. reflux, through tools/bench_case.v built with gcc -prod, times the same
from testdata/flowsheets/reference_dwsim.toml, by Wegstein and by direct
substitution, with no warm start (M8). The four are timed one after the other,
within a minute or two: on a laptop timings move with its power state. Timings
vary from run to run, so this can't reproduce its output byte for byte.
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

import thermo

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m3_bench as gb  # noqa: E402  reflux's build, and the machine
import gen_m3_dwsim as gd  # noqa: E402  DWSIM's loader
import gen_m4_dwsim as g4  # noqa: E402  the reference flowsheet in DWSIM

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m4_bench.json"
FLOWSHEET = "testdata/flowsheets/reference_dwsim.toml"
COLD = 5
WARM = 20


def dwsim_products(sh, o):
    return {name: sh.state(o[g4.STREAMS[name]])["f"] for name in g4.PRODUCTS}


def time_dwsim(auto, label, tol_scale, max_iter):
    """DWSIM's first solves, each of a flowsheet just built, and its solves
    again, from its last answer."""
    sh, o = g4.reference(auto, tol_scale=tol_scale, max_iter=max_iter)
    sh.solve()  # untimed: .NET compiles the flowsheet's code on its first solve
    cold, passes, products = [], None, None
    for _ in range(COLD):
        sh, o = g4.reference(auto, tol_scale=tol_scale, max_iter=max_iter)
        cold.append(sh.solve() * 1e3)
        if passes is None:
            passes = o["RECYCLE"].GetAsObject().IterationsTaken
            products = dwsim_products(sh, o)
    warm = []
    for k in range(1, WARM + 1):
        o["feed"].GetAsObject().SetTemperature(300.0 + 0.01 * k)
        warm.append(sh.solve() * 1e3)
    out = {"label": label, "load_ms": [], "cold_ms": cold, "warm_ms": warm, "passes": passes,
           "warm_passes": o["RECYCLE"].GetAsObject().IterationsTaken, "products": products, "error": ""}
    print(f"  {label}: first {sorted(cold)[len(cold) // 2]:.0f} ms, again {sorted(warm)[len(warm) // 2]:.0f} ms, "
          f"{passes} passes", flush=True)
    return out


def time_reflux(exe, method):
    out = subprocess.run([str(exe), "--flowsheet", FLOWSHEET, method, str(COLD), str(WARM)], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    timing = json.loads(out)
    print(f"  {timing['label']}: first {sorted(timing['cold_ms'])[COLD // 2]:.0f} ms, "
          f"{timing['passes']} passes {timing['error']}", flush=True)
    return timing


def main():
    gd.load_dwsim()
    from DWSIM.Automation import Automation3

    auto = Automation3()
    tmp = tempfile.mkdtemp()
    exe, build = gb.build_reflux(tmp)
    t0 = time.perf_counter()
    timings = [
        time_dwsim(auto, "DWSIM", 1.0, 50),
        time_dwsim(auto, "DWSIM at 1e-8", g4.TIGHT, 500),
        time_reflux(exe, "wegstein"),
        time_reflux(exe, "direct"),
    ]
    print(f"  timed in {time.perf_counter() - t0:.0f} s")
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "python": platform.python_version(),
            "thermo": thermo.__version__,
            "dwsim": str(auto.GetVersion()),
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": gb.cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": build,
        "note": "M4's flowsheet benchmark: the reference flowsheet on DWSIM's constants, solved by DWSIM 9 (its "
        "default recycle tolerances, and them scaled by 1e-8) and reflux (Wegstein, direct substitution). Times in "
        "milliseconds: cold_ms first solves of a flowsheet just built or loaded, warm_ms solves again with the feed "
        "0.01 K warmer each time, load_ms reflux's reading of the file, data and fluid. passes are the loop's in the "
        "first solve, warm_passes in the last solve again; products the products' flows after the first, mol/s.",
        "flowsheet": FLOWSHEET,
        "timings": timings,
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
