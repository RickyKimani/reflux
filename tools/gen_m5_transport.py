# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "numpy==2.2.6", "scipy==1.18.1"]
# ///
"""Write testdata/m5_transport.json: the iapws package's water transport
properties, for M5a's step 6 (docs/m5-plan.md, criterion 7 and D11).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5_transport.py

- "iapws95": at (T, P) states from 275 K to 1273 K and 1 kPa to 100 MPa,
  and densely near the critical point, the package's IAPWS-95 density, and
  its viscosity and thermal conductivity there with their critical
  enhancements (R12-08, R15-11);
- "if97": the same states by its IF97 class: the density, and the
  industrial viscosity and thermal conductivity (R15-11's Eq. 25);
- "if97_critical": IF97's region 3 close to the critical point, 647 K to
  650 K and 21.8 MPa to 22.6 MPa, where ρ̄ crosses R15-11's Table 6 edges,
  and on the critical isochore 1 mK to 0.3 K above Tc (P from the package's
  region 3 at 322 kg/m³), where ζ = (p*/ρ*)(∂ρ/∂p)_T passes 1e3: the
  package's density, and its industrial viscosity and thermal conductivity
  there;
- "tension": the surface tension (R1-76) from 248.15 K to Tc.

Each row: T (K), P (MPa), rho (kg/m3), mu (Pa s), k (W/(m K)).
It reproduces its output byte for byte on the same day.
"""
import datetime
import json
import pathlib
import warnings

import iapws
import numpy as np
from iapws._iapws import _Tension
from iapws.iapws97 import _Region3

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "testdata" / "m5_transport.json"


def grid():
    ts = [275.0 + 25.0 * k for k in range(40)]
    ps = [0.001, 0.01, 0.1, 1.0, 5.0, 10.0, 15.0, 20.0, 22.0, 25.0, 30.0, 50.0, 70.0, 100.0]
    pts = [(t, p) for t in ts for p in ps]
    # near the critical point
    pts += [(640.0 + 1.0 * j, 20.0 + 0.5 * k) for j in range(21) for k in range(11)]
    return pts


def critical_grid():
    pts = [(647.0 + 0.25 * j, 21.8 + 0.02 * k) for j in range(13) for k in range(41)]
    for dt in (0.001, 0.003, 0.01, 0.03, 0.1, 0.3):
        t = 647.096 + dt
        pts.append((t, float(_Region3(322.0, t)["P"])))
    return pts


def region3(t, p):
    try:
        st = iapws.IAPWS97(T=t, P=p)
    except (NotImplementedError, ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None
    if st.region != 3 or st.mu is None or st.k is None or getattr(st, "x", 0) not in (0, 1):
        return None
    return [t, p, float(st.rho), float(st.mu), float(st.k)]


def state(cls, t, p):
    try:
        st = cls(T=t, P=p)
    except (NotImplementedError, ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None
    if getattr(st, "status", 1) != 1 or st.mu is None or st.k is None:
        return None
    if getattr(st, "x", 0) not in (0, 1):
        return None  # a two-phase mixture has no transport properties
    return [t, p, float(st.rho), float(st.mu), float(st.k)]


def dump(o):
    items = []
    for k in sorted(o):
        v = o[k]
        if isinstance(v, list) and v and isinstance(v[0], list):
            body = ",\n".join("  " + json.dumps(r, separators=(",", ":")) for r in v)
            items.append(f" {json.dumps(k)}: [\n{body}\n ]")
        else:
            items.append(f" {json.dumps(k)}: {json.dumps(v, separators=(',', ':'))}")
    return "{\n" + ",\n".join(items) + "\n}\n"


def main():
    out = {
        "provenance": {
            "script": "tools/gen_m5_transport.py",
            "iapws": iapws.__version__,
            "numpy": np.__version__,
            "date": datetime.date.today().isoformat(),
            "columns": ["T K", "P MPa", "rho kg/m3", "mu Pa s", "k W/(m K)"],
        },
        "iapws95": [r for t, p in grid() if (r := state(iapws.IAPWS95, t, p)) is not None],
        "if97": [r for t, p in grid() if (r := state(iapws.IAPWS97, t, p)) is not None],
        "if97_critical": [r for t, p in critical_grid() if (r := region3(t, p)) is not None],
        "tension": [[t, _Tension(t)] for t in [248.15 + 5.0 * k for k in range(80)] if t <= 647.096],
    }
    OUT.write_bytes(dump(out).encode("utf-8"))
    print(f"wrote {OUT}: {len(out['iapws95'])} IAPWS-95 states, {len(out['if97'])} IF97 states, "
          f"{len(out['if97_critical'])} near the critical point")


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
