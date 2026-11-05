"""M10a's exit script, python/examples/phase_diagrams.py (exit criterion 8): it
runs and draws its figure; its diagrams are phase diagrams, checked from the
claims alone (saturation at the pure ends, tie lines that close, a split
exactly where the feed is unstable); and its numbers are the V API's bit for
bit, against tools/capi_reference.v's "examples", on the same grids."""

import importlib.util
import json
import os
import struct
import subprocess
import sys

import numpy as np
import pytest

import reflux
from reflux import _native
from support import ROOT

SCRIPT = ROOT / "python" / "examples" / "phase_diagrams.py"
REFERENCE = os.environ.get("REFLUX_REFERENCE")

_spec = importlib.util.spec_from_file_location("phase_diagrams", SCRIPT)
diagrams = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(diagrams)


@pytest.fixture(scope="module")
def db():
    with reflux.Database() as d:
        yield d


@pytest.fixture(scope="module")
def txy(db):
    return diagrams.txy(db)


@pytest.fixture(scope="module")
def gas(db):
    return diagrams.phase_map(db)


def test_the_script_draws_its_figure(tmp_path):
    # tools/python.vsh, which sets REFLUX_REFERENCE, installs matplotlib, so
    # in CI it's never skipped; pytest alone may lack it.
    if not REFERENCE:
        pytest.importorskip("matplotlib")
    out = tmp_path / "figure.png"
    env = {**os.environ, "PYTHONPATH": str(ROOT / "python"), "PYTHONIOENCODING": "utf-8"}
    r = subprocess.run([sys.executable, str(SCRIPT), str(out)], env=env, capture_output=True,
                       text=True, encoding="utf-8", timeout=600)
    assert r.returncode == 0, f"exit {r.returncode}\n{r.stdout}\n{r.stderr}"
    lines = r.stdout.splitlines()
    assert lines[0] == "bubble points: 101 of 101 have answers"
    assert lines[1] == "dew points: 101 of 101 have answers"
    assert lines[2] == "map: 6480 of 6480 have answers"
    assert lines[-1] == f"wrote {out}"
    png = out.read_bytes()
    # A PNG's IHDR gives its size: 11 by 4.6 inches at 110 dots an inch.
    assert png[:8] == b"\x89PNG\r\n\x1a\n" and png[12:16] == b"IHDR"
    assert struct.unpack(">II", png[16:24]) == (1210, 506)


def test_its_t_x_y_is_a_phase_diagram(db, txy):
    bubble, dew = txy
    p = diagrams.TXY_P
    assert bubble.ok.all() and dew.ok.all()
    assert (bubble.phases == 2).all() and (dew.phases == 2).all()
    assert np.array_equal(bubble.p, np.full(101, p)) and np.array_equal(dew.p, np.full(101, p))
    # Methane is the lighter: the bubble point falls as there's more of it,
    # and between the pure ends the dew point is above it, by 5.1 K at least
    # (next to the ends).
    assert (np.diff(bubble.t) < 0).all() and (np.diff(dew.t) < 0).all()
    assert (dew.t[1:-1] - bubble.t[1:-1] > 1.0).all()
    # The pure ends are each compound's saturation at p, bubble and dew alike
    # (to the bit, measured): its vapour pressure at that T is p. Both solve
    # their fugacities to 1e-12, and ln Psat changes about ten times as fast
    # as ln T, so p is met to 5.2e-11 (n-butane) and 3.1e-11 (methane).
    with db.fluid("PR", diagrams.TXY_NAMES) as f:
        for end, compound in ((0, "n-butane"), (-1, "methane")):
            assert abs(dew.t[end] - bubble.t[end]) <= 1e-12 * bubble.t[end]
            psat = f.psat(bubble.t[end], compound)
            assert psat.ok and abs(psat.p - p) <= 2e-10 * p, (compound, psat.p)
        # A tie line closes: the first bubble of the liquid x is a vapour y
        # whose dew point is the same T, with x its first drop. Measured, to
        # 2.0e-13 in T and 4.2e-13 in x, the two flashes' own tolerance.
        y = bubble.vapor.x[1:-1]
        back = f.dew_t(p, y)
        assert back.ok.all()
        assert (np.abs(back.t - bubble.t[1:-1]) <= 1e-11 * bubble.t[1:-1]).all()
        assert np.abs(back.liquid.x - bubble.liquid.x[1:-1]).max() <= 1e-11


def test_its_map_splits_exactly_where_the_feed_is_unstable(db, gas):
    assert gas.status.shape == (80, 81)
    assert np.array_equal(gas.t[0], diagrams.MAP_T) and np.array_equal(gas.p[:, 0], diagrams.MAP_P)
    assert gas.ok.all(), [gas.reason(tuple(i)) for i in np.argwhere(~gas.ok)[:3]]
    with db.fluid("PR", diagrams.GAS_NAMES) as f:
        st = f.stability(diagrams.MAP_T, diagrams.MAP_P[:, None], diagrams.GAS_Z)
    assert st.ok.all()
    split = gas.phases == 2
    assert np.array_equal(split, ~st.stable)
    assert ((gas.beta > 0) & (gas.beta < 1))[split].all()
    assert np.isin(gas.beta[~split], [0.0, 1.0]).all()
    # Some of each: the map has an envelope in it, not one region.
    assert 1000 < split.sum() < gas.status.size - 1000


def hex_bits(points, key):
    """The points' key, each a list of hex bits, as a 2-d uint64 array."""
    return np.array([[int(h, 16) for h in p[key]] for p in points], dtype=np.uint64)


def assert_the_v_apis(r, points, a, b, z, at):
    """The FlashResult r is the V reference's points bit for bit, a point per
    row in C order, at the inputs a, b and z, broadcast to its points."""
    n, nc = len(points), r.x.shape[-1]
    assert r.status.size == n, at

    def u64(x, width):
        return np.ascontiguousarray(np.broadcast_to(x, r.status.shape + (width,)) if width > 1 else
                                    np.broadcast_to(x, r.status.shape), dtype=np.float64
                                    ).view(np.uint64).reshape(n, width)

    inputs = {"a": u64(a, 1), "b": u64(b, 1), "z": u64(z, nc)}
    for key, want in inputs.items():
        got = hex_bits(points, key) if key == "z" else hex_bits([{key: [p[key]]} for p in points], key)
        assert np.array_equal(want, got), f"{at}: the script's {key} isn't the reference's"
    outputs = {
        "out": np.stack([getattr(r, k) for k in _native.FLASH_VALUES], axis=-1),
        "liquid": np.stack([r.liquid.fraction, r.liquid.v, r.liquid.h, r.liquid.s], axis=-1),
        "vapor": np.stack([r.vapor.fraction, r.vapor.v, r.vapor.h, r.vapor.s], axis=-1),
        "x_liquid": r.liquid.x,
        "x_vapor": r.vapor.x,
    }
    for key, got in outputs.items():
        got = np.ascontiguousarray(got, dtype=np.float64).reshape(n, -1).view(np.uint64)
        bad = np.flatnonzero((got != hex_bits(points, key)).any(axis=1))
        assert bad.size == 0, f"{at}: {key} differs at points {bad[:5]} of {n}"
    assert np.array_equal(r.status.reshape(n), [p["status"] for p in points]), at
    assert np.array_equal(r.phases.reshape(n), [p["n_phases"] for p in points]), at


@pytest.mark.skipif(not REFERENCE, reason="tools/python.vsh test writes the V reference")
def test_its_numbers_are_the_v_apis(txy, gas):
    with open(REFERENCE, encoding="utf-8") as f:
        ref = json.load(f)["examples"]
    bubble, dew = txy
    z = np.stack([diagrams.X, 1 - diagrams.X], axis=-1)
    assert_the_v_apis(bubble, ref["bubble"], diagrams.TXY_P, 0.0, z, "bubble points")
    assert_the_v_apis(dew, ref["dew"], diagrams.TXY_P, 1.0, z, "dew points")
    assert_the_v_apis(gas, ref["map"], diagrams.MAP_T, diagrams.MAP_P[:, None], diagrams.GAS_Z,
                      "map")
