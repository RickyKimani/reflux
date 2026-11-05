# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy==2.5.4"]
# ///
"""Write testdata/m10a_bench.json: M10a's speed (exit criterion 10), what a
call from Python through the C API costs over V's own loop doing the same
work, for properties and a PT flash at 1, 100 and 10,000 points.
tools/benchmarks.vsh adds it to docs/benchmarks.md.

Run from the repository root, with V and gcc on the PATH:

    uv run --python <3.11 or later> tools/gen_m10a_bench.py

It builds the library with gcc -prod (tools/python.vsh build) and
tools/bench_case.v with the same compiler, then for each piece of work and
number of points times V's loop (bench_case --capi) and the package's call,
alternately: on a laptop timings move with its power state.
Each pair is timed for about four seconds, in eight rounds of V then Python, each
after a call that warms the caches. Timings vary from run to run, so this
can't reproduce its output byte for byte.
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

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "testdata" / "m10a_bench.json"
WORK = ("props", "flash_pt")
POINTS = (1, 100, 10_000)
SECONDS = 4.0
ROUNDS = 8
ENV = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")


def cpu_name():
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"],
                             capture_output=True, text=True, timeout=60).stdout
        return out.strip() or platform.processor()
    except (OSError, subprocess.SubprocessError):
        return platform.processor()


def states(n):
    """benchmarks.python_states: one point at 300 K, or n from 220 to 400 K."""
    return np.array([300.0]) if n == 1 else 220.0 + 180.0 * np.arange(n) / (n - 1)


def time_v(exe, what, n, reps):
    out = subprocess.run([str(exe), "--capi", what, str(n), str(reps)], capture_output=True, text=True,
                         check=True).stdout
    return json.loads(out)["us"]


def time_python(fl, what, n, reps):
    T = states(n)
    z = np.array([0.3, 0.3, 0.4])
    f = fl.props if what == "props" else fl.flash_pt
    r = f(T, 20e5, z)  # warms the caches, as V's loop's first call does
    assert np.all(r.ok), f"{what} at {n} points: {r.reason(int(np.argmin(r.ok)))}"
    us = []
    for _ in range(reps):
        t0 = time.perf_counter()
        f(T, 20e5, z)
        us.append((time.perf_counter() - t0) * 1e6)
    return us


def main():
    subprocess.run(["v", "run", "tools/python.vsh", "build"], cwd=ROOT, check=True, env=ENV)
    sys.path.insert(0, str(ROOT / "python"))
    import reflux

    with tempfile.TemporaryDirectory() as tmp:
        exe = pathlib.Path(tmp) / "bench_case.exe"
        subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_case.v"], cwd=ROOT,
                       check=True, env=ENV)
        build = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True,
                                          check=True).stdout)
        fl = reflux.Database().fluid("PR", ["methane", "ethane", "n-butane"])
        timings = []
        for what in WORK:
            for n in POINTS:
                # Enough calls for about a second, from one call's time, in
                # ROUNDS rounds of V then Python: the laptop's power state
                # drifts, and alternating lets it move both alike.
                one = time_v(exe, what, n, 1)[0]
                per = max(1, min(400, int(SECONDS * 1e6 / ROUNDS / max(one, 1.0))))
                v_us, py_us = [], []
                for _ in range(ROUNDS):
                    v_us += time_v(exe, what, n, per)
                    py_us += time_python(fl, what, n, per)
                reps = len(v_us)
                print(f"{what} at {n}: V {statistics.median(v_us):.1f} us, "
                      f"Python {statistics.median(py_us):.1f} us, {reps} calls", flush=True)
                timings.append({"what": what, "points": n, "v_us": [round(x, 3) for x in v_us],
                                "python_us": [round(x, 3) for x in py_us]})
    fixture = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m10a_bench.py",
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "machine": {"cpu": cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": build,
        "timings": timings,
    }
    OUT.write_text(json.dumps(fixture, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
