# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["pythonnet==3.0.5", "chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1", "scipy==1.14.1"]
# ///
"""Write testdata/m9_bench.json: M9's speed (criterion 9), the owner's
reforming PFR solved by DWSIM 9 and by reflux, one solve and the owner's
50-point study, the feed from 700 to 1000 C. tools/benchmarks.vsh adds it
to docs/benchmarks.md.

Run from the repository root, on Windows with DWSIM 9 installed and V on
the PATH:

    uv run --python 3.13 tools/gen_m9_bench.py

About a quarter of an hour, nearly all of it DWSIM's.

DWSIM solves the case as step 0's spike builds it (spike/dwsim_pfr/run.py),
the owner's file's setup: dV = 0.05, Ergun, its default integrator (Radau);
and, for one solve, its explicit Runge-Kutta 4(5), which step 0 found seven
times faster. reflux, through tools/bench_reactor.v built with gcc -prod,
solves testdata/flowsheets/reformer_dwsim.toml, the same case on DWSIM's
rates as typed with reflux's own viscosity, at its default rtol of 1e-8 and
at 1e-4, and testdata/flowsheets/reformer.toml, the published model.

The two programs take turns: a laptop's clock moves DWSIM's solve between
9.4 and 29 s (step 0). After an untimed solve each, which .NET needs to
compile and V's process to warm, the single solve is timed in five turns
(DWSIM, then each reflux file ten times), and each point of the study is
solved by DWSIM and then by each reflux file three times before the next.
A reflux time is the feed's flash and the flowsheet's solve, as a study's
point costs, and the file's loading is timed once apart. Timings vary from
run to run, so this can't reproduce its output byte for byte.
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
sys.path.insert(0, str(ROOT / "spike" / "dwsim_pfr"))
sys.path.insert(0, str(ROOT / "tools"))
import gen_m3_bench as gb  # noqa: E402  the machine's name
import run as spike  # noqa: E402  DWSIM's PFR on the owner's case (M9 step 0)

OUT = ROOT / "testdata" / "m9_bench.json"
TURNS = 5
SINGLE = 10
PER_POINT = 3
STUDY = [973.15 + 300.0 * i / 49 for i in range(50)]  # the owner's: 700 to 1000 C, as step 0 ran it
FILES = [
    {"label": "reflux", "file": "testdata/flowsheets/reformer_dwsim.toml", "rtol": None},
    {"label": "reflux at rtol 1e-4", "file": "testdata/flowsheets/reformer_dwsim.toml", "rtol": "1e-4"},
    {"label": "reflux, the published model", "file": "testdata/flowsheets/reformer.toml", "rtol": None},
]


def build(tmp):
    exe = pathlib.Path(tmp) / "bench_reactor.exe"
    env = dict(os.environ, V_C_ERROR_BUG_REPORT_DISABLED="1")
    subprocess.run(["v", "-prod", "-cc", "gcc", "-nocache", "-o", str(exe), "tools/bench_reactor.v"], cwd=ROOT,
                   check=True, env=env)
    info = json.loads(subprocess.run([str(exe), "--build"], capture_output=True, text=True, check=True).stdout)
    return exe, info


class Reflux:
    """bench_reactor.exe on one file, answering a point at a time."""

    def __init__(self, exe, tmp, spec):
        self.label = spec["label"]
        path = ROOT / spec["file"]
        if spec["rtol"]:
            text = path.read_text(encoding="utf-8").replace('pressure_drop = "ergun"',
                                                            f'pressure_drop = "ergun"\nrtol = {spec["rtol"]}')
            path = pathlib.Path(tmp) / f"rtol_{spec['rtol']}.toml"
            path.write_text(text, encoding="utf-8", newline="\n")
        self.p = subprocess.Popen([str(exe), str(path)], cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                  text=True, bufsize=1)
        self.load_ms = json.loads(self.p.stdout.readline())["load_ms"]

    def at(self, t, repeats):
        self.p.stdin.write(f"{t!r} {repeats}\n")
        self.p.stdin.flush()
        r = json.loads(self.p.stdout.readline())
        if r["error"]:
            raise RuntimeError(f"{self.label} at {t} K: {r['error']}")
        return r

    def close(self):
        self.p.stdin.close()
        self.p.wait(timeout=60)


def x_of(feed, product):
    return (feed[0] - product[0]) / feed[0]  # methane is first in both files and in DWSIM's order


def median(xs):
    s = sorted(xs)
    return s[len(s) // 2]


def main():
    auto = spike.load()
    tmp = tempfile.mkdtemp()
    exe, info = build(tmp)
    sh = spike.Sheet(auto)
    rf = [Reflux(exe, tmp, spec) for spec in FILES]
    t_start = time.perf_counter()
    # Untimed: .NET compiles the case's code on its first solve.
    sh.set_feed(spike.T_FEED)
    sh.solve()
    for r in rf:
        r.at(spike.T_FEED, 1)
    single = {"dwsim_s": [], "dwsim_rk_s": [], "reflux_ms": {r.label: [] for r in rf}}
    for _ in range(TURNS):
        sh.o.InternalSolver = 0
        dt, _ = sh.solve()
        single["dwsim_s"].append(dt)
        dwsim_x = sh.outlet()["x"]
        sh.o.InternalSolver = 1
        dt, _ = sh.solve()
        single["dwsim_rk_s"].append(dt)
        dwsim_rk_x = sh.outlet()["x"]
        sh.o.InternalSolver = 0
        for r in rf:
            a = r.at(spike.T_FEED, SINGLE)
            single["reflux_ms"][r.label].extend(a["ms"])
            single.setdefault("reflux_x", {})[r.label] = x_of(a["feed"], a["product"])
            single.setdefault("reflux_steps", {})[r.label] = a["steps"]
        print(f"  a single solve: DWSIM {single['dwsim_s'][-1]:.1f} s (explicit RK {single['dwsim_rk_s'][-1]:.1f} s); "
              + ", ".join(f"{r.label} {median(single['reflux_ms'][r.label][-SINGLE:]):.1f} ms" for r in rf),
              flush=True)
    single["dwsim_x"], single["dwsim_rk_x"] = dwsim_x, dwsim_rk_x
    points = []
    for t in STUDY:
        sh.set_feed(t)
        dt, _ = sh.solve()
        row = {"t": t, "dwsim_s": dt, "dwsim_x": sh.outlet()["x"], "reflux": {}}
        for r in rf:
            a = r.at(t, PER_POINT)
            row["reflux"][r.label] = {"ms": a["ms"], "x": x_of(a["feed"], a["product"]), "steps": a["steps"]}
        points.append(row)
        print(f"  {t - 273.15:7.2f} C: DWSIM {dt:.1f} s, X {row['dwsim_x']:.6f}; reflux "
              f"{median(row['reflux']['reflux']['ms']):.1f} ms, X {row['reflux']['reflux']['x']:.6f}", flush=True)
    for r in rf:
        r.close()
    sh.set_feed(spike.T_FEED)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m9_bench.py",
            "python": platform.python_version(),
            "dwsim": str(auto.GetVersion()),
            "pythonnet": __import__("importlib.metadata").metadata.version("pythonnet"),
        },
        "machine": {"cpu": gb.cpu_name(), "logical_cpus": os.cpu_count(), "os": platform.platform()},
        "reflux_build": info,
        "note": "M9's speed: the owner's reforming PFR (1 m3, 1000 K and 1 atm, Ergun), solved by DWSIM 9 as the "
                "owner's file sets it up (dV 0.05, its default integrator Radau, and for one solve its explicit "
                "Runge-Kutta) and by reflux (gcc -prod) on the same rates as typed, at rtol 1e-8 and 1e-4, and on "
                "the published model. dwsim_s are seconds a solve; reflux ms milliseconds a solve, the feed's flash "
                "included; load_ms the file's reading, data and fluid, once. x is methane's conversion; steps the "
                "PFR integrator's accepted steps. The study is the owner's: 50 feed temperatures, 700 to 1000 C, "
                "each solved by DWSIM and then by each reflux file, in turns.",
        "load_ms": {r.label: r.load_ms for r in rf},
        "single": single,
        "study": points,
        "seconds": round(time.perf_counter() - t_start, 1),
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    dw = sum(p["dwsim_s"] for p in points)
    rx = sum(median(p["reflux"]["reflux"]["ms"]) for p in points) / 1e3
    print(f"wrote {OUT.relative_to(ROOT)}: the study, DWSIM {dw:.0f} s, reflux {rx:.2f} s")


if __name__ == "__main__":
    main()
