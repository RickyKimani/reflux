"""The C API against the V API, bit for bit (M10a, D12 and exit criterion 3).
tools/capi_reference.v, built by tools/python.vsh with the library's own
compiler and flags, writes the V API's answers for a reference set of states
covering every function; here the package computes them again, through the
library, and every bit must agree: values, NaN where there's no answer,
statuses, phases and the text the library writes. Bits can differ between
compilers and between maths libraries, so the two must come from one build,
and the reference is never committed. Run through `v run tools/python.vsh
test`, which sets REFLUX_REFERENCE; without it, these tests are skipped."""

import json
import os
import struct

import numpy as np
import pytest

import reflux
from reflux import _native, points

PATH = os.environ.get("REFLUX_REFERENCE")
pytestmark = pytest.mark.skipif(not PATH, reason="tools/python.vsh test writes the V reference")


def value(h):
    """A double from the hex of its bits."""
    return struct.unpack("<d", int(h, 16).to_bytes(8, "little"))[0]


def values(hs):
    return np.array([value(h) for h in hs], dtype=np.float64)


def same(got, hs):
    """got's bits are exactly hs's."""
    got = np.ascontiguousarray(np.asarray(got, dtype=np.float64)).ravel().view(np.uint64)
    want = np.array([int(h, 16) for h in hs], dtype=np.uint64)
    return np.array_equal(got, want)


@pytest.fixture(scope="module")
def ref():
    with open(PATH, encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


def fluid_of(db, c):
    overrides = {(a, b): value(k) for (a, b), k in zip(c["overrides"], c["ks"])}
    return db.fluid(c["equation"], c["names"], alpha=c["alpha"], kij=c["kij"], overrides=overrides)


def test_the_text_is_the_v_apis(ref, db):
    for c in ref["fluids"]:
        with fluid_of(db, c) as f:
            assert _native.text_out(_native.lib.reflux_fluid_setup, f.handle) == c["setup"]
    for name, text in ref["compounds"].items():
        assert _native.text_out(_native.lib.reflux_compound, db.handle, name.encode()) == text


def test_properties_are_the_v_apis(ref, db):
    checked = 0
    for c in ref["fluids"]:
        with fluid_of(db, c) as f:
            for root in ("liquid", "vapor", "stable"):
                pts = [p for p in c["props"] if p["root"] == root]
                r = f.props(values([p["t"] for p in pts]), values([p["p"] for p in pts]),
                            np.array([values(p["z"]) for p in pts]), root=root)
                for i, p in enumerate(pts):
                    at = f"{c['names']} {root} point {i}"
                    assert r.status[i] == p["status"], at
                    assert same([getattr(r, n)[i] for n in _native.PROPS], p["values"]), at
                    assert same(r.ln_phi[i], p["ln_phi"]) and r.phase[i] == p["phase"], at
                    checked += 1
    assert checked == 4 * 3 * 20


def test_vapour_pressures_are_the_v_apis(ref, db):
    checked = 0
    for c in ref["fluids"]:
        with fluid_of(db, c) as f:
            for j in range(len(c["names"])):
                pts = [p for p in c["psat"] if p["component"] == j]
                r = f.psat(values([p["t"] for p in pts]), j)
                for i, p in enumerate(pts):
                    assert r.status[i] == p["status"] and same([r.p[i]], [p["p"]])
                    checked += 1
    assert checked == 4 * 11


def test_stability_is_the_v_apis(ref, db):
    for c in ref["fluids"]:
        with fluid_of(db, c) as f:
            pts = c["stability"]
            r = f.stability(values([p["t"] for p in pts]), values([p["p"] for p in pts]),
                            np.array([values(p["z"]) for p in pts]))
            for i, p in enumerate(pts):
                assert r.status[i] == p["status"] and r.stable[i] == (p["stable"] == 1)
                assert same([r.tm[i]], [p["tm"]]) and same(r.w[i], p["w"])


def test_every_flash_is_the_v_apis(ref, db):
    seen = {}
    for c in ref["fluids"]:
        with fluid_of(db, c) as f:
            for spec in _native.FLASHES:
                pts = [p for p in c["flash"] if p["spec"] == spec]
                if not pts:
                    continue
                r = points.flash(f, spec, values([p["a"] for p in pts]), values([p["b"] for p in pts]),
                                 np.array([values(p["z"]) for p in pts]))
                for i, p in enumerate(pts):
                    at = f"{c['names']} {spec} point {i}"
                    assert r.status[i] == p["status"] and r.phases[i] == p["n_phases"], at
                    assert same([getattr(r, n)[i] for n in _native.FLASH_VALUES], p["out"]), at
                    for view, name in ((r.liquid, "liquid"), (r.vapor, "vapor")):
                        assert same([view.fraction[i], view.v[i], view.h[i], view.s[i]], p[name]), at
                        assert same(view.x[i], p["x_" + name]), at
                    seen[spec] = seen.get(spec, 0) + (p["status"] == 0)
    # Every flash has answers among them, not only failures.
    assert set(seen) == set(_native.FLASHES) and all(n > 0 for n in seen.values()), seen
