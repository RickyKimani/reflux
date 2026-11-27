# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "numpy==2.2.6", "scipy==1.18.1"]
# ///
"""Write testdata/m5_if97.json: the references for M5a's IF97
(docs/m5-plan.md, criteria 4 and 5).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_if97.py [part ...]

Parts (all by default; a part not named keeps what the file has):
  coefficients  the iapws package's IF97 coefficients, to hold
                water/if97/coefficients*.v's transcription from the releases
                against (D10): its constants module's arrays, the literal
                lists, tuples and dicts assigned in its functions, and, for
                the few functions that write theirs inline (B2bc, h3ab,
                region 3's n1, T3ef), every number in the function's body;
  states        the iapws package's IF97 at grids of states in every region:
                the basic equations, the backward equations as published,
                region 4, the boundaries, and its steam-table class from
                (p, T), (p, h), (p, s), (T, x) and (p, x).

Everything is in the package's units (MPa, K, kJ/kg, kJ/(kg K), m3/kg).
It reproduces its output byte for byte on the same day.
"""
import ast
import datetime
import inspect
import json
import math
import pathlib
import sys
import warnings

import iapws
import numpy as np
from iapws import _iapws97Constants as C
from iapws import iapws97 as m

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5_if97.json"

# Functions whose coefficients are written inline: every number in them.
INLINE = ["_P_2bc", "_hbc_P", "_h_3ab", "_Region3", "_tef_P"]


def plain(v):
    if isinstance(v, np.ndarray):
        return [plain(x) for x in v.tolist()]
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    if isinstance(v, dict):
        return {str(k): plain(x) for k, x in v.items()}
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v


def coefficients():
    # the (h, s) backward equations (SR2-01, SR4-04) aren't in M5a
    hs = ("_hs", "h1_s", "h3a_s", "h2ab_s", "h2c3b_s", "h13_s")
    const = {k: plain(v) for k, v in sorted(vars(C).items())
             if isinstance(v, np.ndarray) and "product" not in k and "less" not in k
             and not any(t in k for t in hs)}
    tree = ast.parse(inspect.getsource(m))
    local, inline = {}, {}
    for fn in tree.body:
        if not isinstance(fn, ast.FunctionDef):
            continue
        got = {}
        for node in ast.walk(fn):
            if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                try:
                    val = ast.literal_eval(node.value)
                except (ValueError, SyntaxError, TypeError):
                    continue
                if isinstance(val, (list, tuple, dict)) and val:
                    got[node.targets[0].id] = plain(val)
        if got:
            local[fn.name] = got
        if fn.name in INLINE:
            body = fn.body[1:] if isinstance(fn.body[0], ast.Expr) else fn.body  # no docstring
            nums = []
            for stmt in body:
                for node in ast.walk(stmt):
                    if isinstance(node, ast.Constant) and isinstance(node.value, float):
                        nums.append(node.value)
            inline[fn.name] = nums
    return {"const": const, "local": local, "inline": inline}


def dump(o, ind=0):
    """JSON with sorted keys, a key per line and a list of numbers, or each
    row of a list of rows, on one line."""
    pad = " " * (ind + 1)
    if isinstance(o, dict):
        items = [f"{pad}{json.dumps(k)}: {dump(o[k], ind + 1)}" for k in sorted(o)]
        return "{\n" + ",\n".join(items) + "\n" + " " * ind + "}"
    if isinstance(o, list) and o and all(isinstance(x, list) for x in o):
        rows = [pad + json.dumps(x, separators=(",", ":")) for x in o]
        return "[\n" + ",\n".join(rows) + "\n" + " " * ind + "]"
    return json.dumps(o, separators=(",", ":"))


def main():
    parts = sys.argv[1:] or ["coefficients", "states"]
    data = json.loads(OUT.read_text(encoding="utf-8")) if OUT.exists() else {}
    data["provenance"] = {
        "script": "tools/gen_m5_if97.py",
        "iapws": iapws.__version__,
        "numpy": np.__version__,
        "date": datetime.date.today().isoformat(),
        "units": "MPa, K, kJ/kg, kJ/(kg K), m3/kg",
    }
    if "coefficients" in parts:
        data["coefficients"] = coefficients()
    if "states" in parts:
        data["states"] = states()
    OUT.write_bytes((dump(data) + "\n").encode("utf-8"))
    print(f"wrote {OUT}")


PROPS = ["region", "T", "P", "v", "h", "s", "cp", "cv", "w", "x"]


def row(st):
    """A state of the package's steam-table class, NaN as null."""
    out = []
    for k in PROPS:
        v = getattr(st, k, None)
        v = None if v is None or (isinstance(v, float) and not math.isfinite(v)) else float(v)
        out.append(v)
    return out


def pt_grid():
    """(P, T) over the range: 273.15 K, then every 25 K from 275.15 K, off
    the regions' boundaries (623.15, 863.15 and 1073.15 K), then every 100 K
    to 2273.15 K, and pressures from 1 Pa to 100 MPa."""
    ts = [273.15] + [275.15 + 25.0 * k for k in range(32)] + [1173.15 + 100.0 * k for k in range(12)]
    ps = [1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 3.0, 5.0, 10.0, 15.0, 17.0, 20.0, 22.0, 23.0, 25.0,
          30.0, 40.0, 50.0, 70.0, 100.0]
    return [(p, t) for t in ts for p in ps if not (t > 1073.15 and p > 50.0)]


def state(**kw):
    """The package's state, or None where it refuses one (it raises for a
    state out of its bounds)."""
    try:
        st = iapws.IAPWS97(**kw)
    except (NotImplementedError, ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None
    return st if st.status == 1 else None


def states():
    out = {"columns": PROPS}
    # each row of "pt" is the (P, T) asked, then the package's state: its P
    # in region 3 is the one its density gives back, which differs by rounding
    asked = [(p, t, st) for p, t in pt_grid() if (st := state(P=p, T=t)) is not None]
    good = [st for _, _, st in asked]
    out["pt"] = [[p, t] + row(st) for p, t, st in asked]
    # (P, h) and (P, s) from those states' h and s at the P asked: each row
    # the state's index in "pt", then the package's answer
    out["ph"] = [[k] + row(b) for k, (p, _, st) in enumerate(asked) if (b := state(P=p, h=st.h)) is not None]
    out["ps"] = [[k] + row(b) for k, (p, _, st) in enumerate(asked) if (b := state(P=p, s=st.s)) is not None]
    # saturated states: (T, x) and (P, x), x of 0, 0.35 and 1
    tx = []
    for k in range(75):
        t = 273.15 + 5.0 * k
        for x in (0.0, 0.35, 1.0):
            st = state(T=t, x=x)
            if st is not None:
                tx.append(row(st))
    out["tx"] = tx
    px = []
    for p in [0.000611213, 0.001, 0.01, 0.1, 1.0, 5.0, 10.0, 15.0, 16.5, 17.0, 18.0, 19.0, 20.0, 21.0,
              21.5, 22.0]:
        for x in (0.0, 0.35, 1.0):
            st = state(P=p, x=x)
            if st is not None:
                px.append(row(st))
    out["px"] = px
    # the backward equations as published, at each grid state of their region
    back = {"T1_ph": [], "T1_ps": [], "T2_ph": [], "T2_ps": [], "T3_ph": [], "v3_ph": [], "T3_ps": [],
            "v3_ps": [], "v3_pt": []}
    for st in good:
        p, h, s, t = st.P, st.h, st.s, st.T
        if st.region == 1:
            back["T1_ph"].append([p, h, m._Backward1_T_Ph(p, h)])
            back["T1_ps"].append([p, s, m._Backward1_T_Ps(p, s)])
        elif st.region == 2:
            back["T2_ph"].append([p, h, m._Backward2_T_Ph(p, h)])
            back["T2_ps"].append([p, s, m._Backward2_T_Ps(p, s)])
        elif st.region == 3:
            back["T3_ph"].append([p, h, m._Backward3_T_Ph(p, h)])
            back["v3_ph"].append([p, h, m._Backward3_v_Ph(p, h)])
            back["T3_ps"].append([p, s, m._Backward3_T_Ps(p, s)])
            back["v3_ps"].append([p, s, m._Backward3_v_Ps(p, s)])
            back["v3_pt"].append([p, t, m._Backward3_v_PT(p, t)])
    # region 3's v(p, T) near the critical point, where the auxiliary
    # subregions are: a fine grid of (p, T) in region 3
    for k in range(41):
        p = 20.6 + 0.05 * k
        for j in range(41):
            t = 640.0 + 0.4 * j
            if p > m._P23_T(t) and t <= m._t_P(p) if t > 623.15 else False:
                back["v3_pt"].append([p, t, m._Backward3_v_PT(p, t)])
    out["backward"] = back
    # region 4 and the boundaries
    out["psat"] = [[273.15 + 2.5 * k, m._PSat_T(273.15 + 2.5 * k)] for k in range(150)]
    out["tsat"] = [[p, m._TSat_P(p)] for p in [0.000611213 * 1.1 ** k for k in range(110)] if p <= 22.064]
    out["p23"] = [[623.15 + 10.0 * k, m._P23_T(623.15 + 10.0 * k)] for k in range(25)]
    out["p3sat_h"] = [[h, m._PSat_h(h)] for h in [1670.86 + 20.0 * k for k in range(45)]]
    out["p3sat_s"] = [[s, m._PSat_s(s)] for s in [3.7783 + 0.025 * k for k in range(58)]]
    out["h3ab"] = [[p, m._h_3ab(p)] for p in [22.064 + 2.0 * k for k in range(40)]]
    return out


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
