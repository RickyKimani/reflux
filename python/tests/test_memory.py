"""Memory (M10a step 6, exit criterion 6): over 100,000 calls of every kind,
fluids made and freed among them, what the library holds levels off. A leak
would keep it growing round after round. It runs in a process of its own,
so that the heap is this test's alone (D11), and it replaces M1's spike's
collector stress, which made 500 calls.

What it holds is the bytes live after a full collection, which
reflux.collect() returns (M8 step 0). Until then this test watched the
collector's heap size, free blocks included, and the collector grows that in
one large step when it chooses to and never gives it back: on Ubuntu it
doubled once, 5,580 KiB to about 11 MB, and failed a run with no leak in it."""

from support import isolated

CODE = """
import numpy as np
import reflux

db = reflux.Database()
fl = db.fluid("PR", ["methane", "n-butane"])
z = np.array([0.5, 0.5])
live = [reflux.collect()]
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
    live.append(reflux.collect())
print(calls, *live)
"""


def test_what_the_library_holds_levels_off_over_100000_calls():
    calls, *live = map(int, isolated(CODE, timeout=600).split())
    assert calls == 100_000
    # After the first round's warm-up, the second half's most is within
    # SLACK of the first half's.
    first, second = max(live[1:6]), max(live[6:])
    assert second <= first + SLACK, [n >> 10 for n in live]


# How far what's live may wander between rounds. The collector counts a
# block as live while anything in it is, so the count moves by blocks:
# measured over three runs, it stayed within 68 KiB from the first round
# on (2,836 KiB every round in one, 948 to 1,016 KiB in another), and the
# second half's most was never above the first's. A fluid's handle kept
# after each of a round's 100 fluids is freed would add more than this in
# one round.
SLACK = 256 << 10
