"""Memory (M10a step 6, exit criterion 6): over 100,000 calls of every kind,
fluids made and freed among them, the collector's heap levels off. A leak
would keep it growing round after round. It runs in a process of its own,
so that the heap is this test's alone (D11), and it replaces M1's spike's
collector stress, which made 500 calls."""

from support import isolated

CODE = """
import numpy as np
import reflux

db = reflux.Database()
fl = db.fluid("PR", ["methane", "n-butane"])
z = np.array([0.5, 0.5])
heaps = [reflux.collect()]
calls = 0
for _ in range(10):
    for i in range(10_000):
        k = i % 100
        if k == 0:
            # A fluid made and freed: its handle mustn't outlive it.
            with db.fluid("SRK", ["ethane", "propane"], overrides={("ethane", "propane"): 0.01}):
                pass
        elif k == 1:
            fl.flash_pt(300.0, 20e5, z)
        elif k == 2:
            fl.stability(300.0, 20e5, z)
        elif k == 3:
            db.compound("water")
        elif k == 4:
            fl.flash_ph(20e5, -1e4, z)
        else:
            fl.props(300.0 + 0.01 * k, 1e5, z)
        calls += 1
    heaps.append(reflux.collect())
print(calls, *heaps)
"""


def test_the_heap_levels_off_over_100000_calls():
    calls, *heaps = map(int, isolated(CODE, timeout=600).split())
    assert calls == 100_000
    # The collector grows its heap in steps until it fits the work, then
    # stays: the second half's largest is within 10% (and 1 MiB) of the
    # first half's, after the first round's warm-up.
    first, second = max(heaps[1:6]), max(heaps[6:])
    assert second <= 1.1 * first + (1 << 20), [h >> 10 for h in heaps]
