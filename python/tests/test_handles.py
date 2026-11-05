"""Handles (D4), M1's spike's question answered for a real object: a fluid
that only Python holds survives the library's collections, and stays the
same fluid while databases and fluids are made and freed around it. The
spike lost such a handle to the collector. It runs in a process of its own,
so that a crash is reported (D11)."""

from support import isolated

CODE = """
import reflux
from reflux import _native

db = reflux.Database()
fl = db.fluid("PR", ["methane", "n-butane"], overrides={("methane", "n-butane"): 0.0123})
first = _native.text_out(_native.lib.reflux_fluid_setup, fl.handle)
db.close()
del db
for _ in range(40):
    with reflux.Database() as d:
        d.fluid("SRK", ["ethane", "propane", "nitrogen"]).close()
    reflux.collect()
print(_native.text_out(_native.lib.reflux_fluid_setup, fl.handle) == first, fl.kij[0, 1])
"""


def test_a_fluid_survives_collections():
    assert isolated(CODE).split() == ["True", "0.0123"]
