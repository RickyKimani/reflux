"""Threads (D5). Python threads call in at once, each registered with the
collector for the length of its call: without that, a collection from one
aborts with "Collecting from unknown thread", on Windows as on Linux. And the
library's first call, which starts V's runtime and mustn't race (reflux.h), is
the package's own, made by the import. Each test runs in a process of its
own, so that a crash is reported, not hidden."""

from support import isolated

CALLS = """
import threading
import reflux
from reflux import _native

wrong = []
db = reflux.Database()
water = db.compound("water")


def work():
    for i in range(200):
        if reflux.collect() <= 0 or reflux.version() != reflux.__version__:
            wrong.append("collect or version")
        try:
            _native.free(10**15)
        except reflux.BadHandle as e:
            if str(e) != "1000000000000000 is not a handle":
                wrong.append(str(e))
        else:
            wrong.append("no BadHandle")
        # One database, shared: fluids made and freed, and files read.
        if i % 10 == 0:
            with db.fluid("PR", ["methane", "n-butane"]) as fl:
                if fl.kij[0, 1] != 0.0244:
                    wrong.append(f"kij {fl.kij[0, 1]}")
            if db.compound("water") != water:
                wrong.append("water")


threads = [threading.Thread(target=work) for _ in range(8)]
for t in threads:
    t.start()
for t in threads:
    t.join()
print(len(wrong), wrong[:3])
"""


def test_threads_call_and_collect_at_once():
    assert isolated(CALLS).strip() == "0 []"


def test_an_import_from_a_thread_that_ends():
    # The import makes the first call, so V's runtime starts in that thread,
    # and here the thread ends before anything else calls in.
    code = ("import threading\n"
            "t = threading.Thread(target=lambda: __import__('reflux'))\n"
            "t.start()\n"
            "t.join()\n") + CALLS
    assert isolated(code).strip().endswith("0 []")


EVERY = """
import threading
import numpy as np
import reflux

fl = reflux.Database().fluid("PR", ["methane", "ethane", "n-butane"])
z = np.array([0.3, 0.3, 0.4])
T = np.linspace(220.0, 400.0, 8)[:, None]
P = np.linspace(2e5, 60e5, 4)
pt0 = fl.flash_pt(T, P, z)
U = pt0.h - pt0.p * pt0.v


def arrays(r):
    out = [v for k, v in sorted(vars(r).items())
           if not k.startswith("_") and isinstance(v, np.ndarray)]
    for view in (getattr(r, "liquid", None), getattr(r, "vapor", None)):
        if view is not None:
            out += [view.fraction, view.x, view.v, view.h, view.s]
    return out


def run():
    # Every vectorised function, every flash among them.
    return [arrays(r) for r in (
        fl.props(T, P, z), fl.props(T, P, z, root="liquid"), fl.psat(T[:, 0]),
        fl.stability(T, P, z), fl.flash_pt(T, P, z), fl.flash_ph(P, pt0.h, z),
        fl.flash_ps(P, pt0.s, z), fl.flash_tv(T, pt0.v, z), fl.flash_uv(U, pt0.v, z),
        fl.flash_tvf(T[:, 0], 0.5, z), fl.flash_pvf(P, 0.5, z), fl.bubble_p(T[:, 0], z),
        fl.dew_t(P, z))]


def same(a, b):
    return all(len(x) == len(y) and all(np.array_equal(u, v, equal_nan=True) for u, v in zip(x, y))
               for x, y in zip(a, b))


serial = run()
wrong = []


def work():
    for _ in range(3):
        if not same(run(), serial):
            wrong.append("different")


threads = [threading.Thread(target=work) for _ in range(8)]
for t in threads:
    t.start()
for t in threads:
    t.join()
print(len(wrong), sum(len(x) for x in serial))
"""


def test_every_vectorised_function_under_threads():
    # Exit criterion 7: 8 threads calling every vectorised function at once,
    # three times each, get the serial answers bit for bit, NaN included.
    wrong, arrays = isolated(EVERY).split()
    assert wrong == "0" and int(arrays) > 100
