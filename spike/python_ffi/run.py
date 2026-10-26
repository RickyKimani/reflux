# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy>=2,<3"]
# ///
"""M1 step 9 spike: call a reflux shared library from Python with NumPy.

    uv run spike/python_ffi/run.py <path to reflux_ffi.dll or .so>

Normally started by spike.vsh, which builds the library first. Each
scenario runs in its own Python process, so a crash in one is reported
as a crash instead of hiding the rest.
"""

import ctypes
import json
import math
import pathlib
import subprocess
import sys
import threading
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
BUTANE = (425.1, 37.96e5, 0.200)  # Tc (K), Pc (Pa), omega: as in the M1 tests
LIQUID, VAPOR = 0, 1

c_double_p = ctypes.POINTER(ctypes.c_double)
c_i32_p = ctypes.POINTER(ctypes.c_int32)


def load(path):
    lib = ctypes.CDLL(str(path))
    lib.reflux_scale.argtypes = [c_double_p, ctypes.c_int32, ctypes.c_double]
    lib.reflux_scale.restype = None
    lib.reflux_pr_z_many.argtypes = [
        ctypes.c_int32, c_double_p, c_double_p, ctypes.c_int32,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, c_double_p, c_i32_p,
    ]
    lib.reflux_pr_z_many.restype = ctypes.c_int32
    for name in ("reflux_model_new", "reflux_model_new_pinned"):
        f = getattr(lib, name)
        f.argtypes = [ctypes.c_double] * 3
        f.restype = ctypes.c_void_p
    lib.reflux_model_free.argtypes = [ctypes.c_void_p]
    lib.reflux_model_free.restype = None
    lib.reflux_model_z.argtypes = [ctypes.c_void_p, ctypes.c_double, ctypes.c_double, ctypes.c_int32]
    lib.reflux_model_z.restype = ctypes.c_double
    lib.reflux_gc_collect.restype = ctypes.c_uint64
    lib.reflux_churn.argtypes = [ctypes.c_int32]
    lib.reflux_churn.restype = ctypes.c_double
    return lib


def ptr(a, ctype=c_double_p):
    # NumPy hands over its own buffer: no copy. The array must be contiguous
    # and of the right dtype, or the C side would read garbage.
    assert a.flags["C_CONTIGUOUS"]
    return a.ctypes.data_as(ctype)


def z_many(lib, t, p, root):
    t = np.ascontiguousarray(t, dtype=np.float64)
    p = np.ascontiguousarray(p, dtype=np.float64)
    z = np.empty_like(t)
    status = np.empty(t.shape, dtype=np.int32)
    failures = lib.reflux_pr_z_many(t.size, ptr(t), ptr(p), root, *BUTANE, ptr(z), ptr(status, c_i32_p))
    return z, status, failures


def fixture_points():
    """PR n-butane points from testdata/m1_thermo.json (computed by thermo)."""
    fx = json.loads((ROOT / "testdata" / "m1_thermo.json").read_text())
    pts = []
    for c in fx["cases"]:
        if c["eos"] == "pr" and c["components"] == [1]:
            for ph in c["phases"]:
                pts.append((c["t"], c["p"], LIQUID if ph["phase"] == "liquid" else VAPOR, ph["z"]))
    return pts


# ---- scenarios -------------------------------------------------------------


def s_zero_copy(lib):
    a = np.arange(1.0, 6.0)
    lib.reflux_scale(ptr(a), a.size, 2.0)
    assert np.array_equal(a, [2.0, 4.0, 6.0, 8.0, 10.0]), a
    return "NumPy buffer modified in place"


def s_matches_thermo(lib):
    pts = fixture_points()
    worst = 0.0
    for t, p, root, z_ref in pts:
        z, status, _ = z_many(lib, [t], [p], root)
        assert status[0] == 0
        worst = max(worst, abs(z[0] - z_ref) / z_ref)
    assert worst < 1e-9, worst
    return f"{len(pts)} PR points agree with thermo to {worst:.1e}"


def s_vectorised(lib):
    n = 10_000
    t = np.linspace(300.0, 500.0, n)
    p = np.full(n, 5e5)
    t0 = time.perf_counter()
    z, status, failures = z_many(lib, t, p, VAPOR)
    dt = time.perf_counter() - t0
    assert failures == 0 and np.all(status == 0) and np.all((z > 0.5) & (z < 1.0))
    return f"{n} points in one call: {dt / n * 1e6:.1f} us/point"


def s_per_call_overhead(lib):
    n = 10_000
    h = lib.reflux_model_new_pinned(*BUTANE)
    t = np.linspace(300.0, 500.0, n)
    t0 = time.perf_counter()
    for ti in t:
        lib.reflux_model_z(h, float(ti), 5e5, VAPOR)
    dt_loop = time.perf_counter() - t0
    t0 = time.perf_counter()
    z_many(lib, t, np.full(n, 5e5), VAPOR)
    dt_vec = time.perf_counter() - t0
    lib.reflux_model_free(h)
    return f"Python loop {dt_loop / n * 1e6:.1f} us/point vs one vectorised call {dt_vec / n * 1e6:.1f} us/point"


def s_gc_stress(lib):
    t = np.linspace(300.0, 500.0, 1000)
    p = np.full(1000, 5e5)
    heaps = [lib.reflux_gc_collect()]
    ref = None
    for _ in range(5):
        for _ in range(100):
            z, _, failures = z_many(lib, t, p, VAPOR)
            assert failures == 0
            ref = z if ref is None else ref
            assert np.array_equal(z, ref)
        heaps.append(lib.reflux_gc_collect())
    # The collector grows its heap in steps until it fits the workload, then
    # it should stay put. A leak would keep growing round after round.
    assert heaps[-1] <= 1.1 * heaps[-2] + (64 << 10), heaps
    return "500 calls x 1000 points; heap after each 100 calls (KiB): " + ", ".join(str(h // 1024) for h in heaps)


def _handle_survives(lib, new):
    h = new(*BUTANE)
    z0 = lib.reflux_model_z(h, 400.0, 20e5, VAPOR)
    for _ in range(20):
        lib.reflux_gc_collect()
        lib.reflux_churn(200_000)
    z1 = lib.reflux_model_z(h, 400.0, 20e5, VAPOR)
    return z0, z1


def s_handle_unpinned(lib):
    # This scenario reports rather than asserts: whether an unpinned handle
    # survives is the question, and either answer is a finding.
    z0, z1 = _handle_survives(lib, lib.reflux_model_new)
    if z0 == z1:
        return f"survived this run (Z = {z1:.6f}), but nothing guarantees it"
    return f"UNSAFE: the collector freed it (Z {z0:.6f} -> {z1}); handles must be pinned"


def s_handle_pinned(lib):
    z0, z1 = _handle_survives(lib, lib.reflux_model_new_pinned)
    assert z0 == z1 and not math.isnan(z1), (z0, z1)
    return f"survived (Z = {z1:.6f})"


def s_python_thread(lib):
    result = {}

    def work():
        z, _, failures = z_many(lib, np.linspace(300.0, 500.0, 1000), np.full(1000, 5e5), VAPOR)
        lib.reflux_gc_collect()  # a collection started from that thread
        result["ok"] = failures == 0

    th = threading.Thread(target=work)
    th.start()
    th.join()
    assert result.get("ok"), result
    return "called and collected from a non-main Python thread"


SCENARIOS = {
    "zero_copy": s_zero_copy,
    "matches_thermo": s_matches_thermo,
    "vectorised": s_vectorised,
    "per_call_overhead": s_per_call_overhead,
    "gc_stress": s_gc_stress,
    "handle_unpinned": s_handle_unpinned,
    "handle_pinned": s_handle_pinned,
    "python_thread": s_python_thread,
}


def main():
    if len(sys.argv) == 3:  # child: run one scenario
        print(SCENARIOS[sys.argv[2]](load(sys.argv[1])), flush=True)
        return
    lib_path = pathlib.Path(sys.argv[1]).resolve()
    print(f"library: {lib_path.name}   python {sys.version.split()[0]}   numpy {np.__version__}")
    failed = 0
    for name in SCENARIOS:
        r = subprocess.run([sys.executable, __file__, str(lib_path), name], capture_output=True, text=True, timeout=600)
        if r.returncode == 0:
            print(f"  ok     {name:18} {r.stdout.strip()}")
        else:
            failed += 1
            last = (r.stderr.strip().splitlines() or ["(no output)"])[-1]
            kind = "FAIL" if "Error" in last else f"CRASH (exit {r.returncode})"
            print(f"  {kind:6} {name:18} {last}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
