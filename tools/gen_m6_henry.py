# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4"]
# ///
"""Write testdata/m6_henry_iapws.json: IAPWS's Henry's constants for gases
in water (G7-04, 2004; Fernandez-Prini, Alvarez and Harvey, 2003), by the
iapws package's implementation, for M6's Henry's law (D6 as amended after
step 6).

For each of the guideline's 14 gases in H2O: the coefficients A, B, C and
the range as the package has them, and ln(kH/Pa) at 25 temperatures evenly
across the range, with water's ln(p*/Pa) by Wagner and Pruss's equation
there (the package's, which the guideline takes p* from).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_henry.py

It reproduces its output byte for byte on the same day.
"""
import datetime
import inspect
import json
import math
import pathlib
import re
import warnings

import iapws
from iapws import _iapws

ROOT = pathlib.Path(__file__).resolve().parent.parent

# The guideline's names, with each gas's CAS number.
GASES = [("He", "7440-59-7"), ("Ne", "7440-01-9"), ("Ar", "7440-37-1"), ("Kr", "7439-90-9"),
         ("Xe", "7440-63-3"), ("H2", "1333-74-0"), ("N2", "7727-37-9"), ("O2", "7782-44-7"),
         ("CO", "630-08-0"), ("CO2", "124-38-9"), ("H2S", "7783-06-4"), ("CH4", "74-82-8"),
         ("C2H6", "74-84-0"), ("SF6", "2551-62-4")]


def tables():
    """The package's Table 2 and ranges, read from its source, where they're
    local to _Henry."""
    src = inspect.getsource(_iapws._Henry)
    limit, par = {}, {}
    for name, a, b in re.findall(r'"(\w+)": \(([-\d.]+), ([-\d.]+)\)', src):
        limit[name] = (float(a), float(b))
    for name, a, b, c in re.findall(r'"(\w+)": \(([-\d.]+), ([-\d.]+), ([-\d.]+)\)', src):
        par[name] = (float(a), float(b), float(c))
    return limit, par


def ln_psat(t):
    """Water's ln(p*/Pa) as _Henry computes it (eq. 4 of the guideline)."""
    tr = t / _iapws.Tc
    tau = 1 - tr
    ai = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
    bi = [1, 1.5, 3, 3.5, 4, 7.5]
    return math.log(_iapws.Pc * 1e6) + 1 / tr * sum(a * tau ** b for a, b in zip(ai, bi))


def main():
    warnings.simplefilter("ignore")
    limit, par = tables()
    out = {
        "source": {"generator": "tools/gen_m6_henry.py", "iapws": iapws.__version__,
                   "generated": datetime.date.today().isoformat(),
                   "guideline": "IAPWS G7-04, Table 2, H2O"},
        "note": "ln_kh is ln(kH/Pa), kH = lim f/x of the gas in water; ln_psat is water's ln(p*/Pa).",
        "gases": [],
    }
    for name, cas in GASES:
        lo, hi = limit[name]
        ts = [lo + (hi - lo) * k / 24 for k in range(25)]
        out["gases"].append({
            "name": name, "cas": cas, "abc": list(par[name]), "tmin": lo, "tmax": hi,
            "T": ts,
            "ln_kh": [math.log(_iapws._Henry(t, name) * 1e6) for t in ts],
            "ln_psat": [ln_psat(t) for t in ts],
        })
    (ROOT / "testdata" / "m6_henry_iapws.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8",
                                                           newline="\n")
    print("wrote testdata/m6_henry_iapws.json:", len(out["gases"]), "gases")


if __name__ == "__main__":
    main()
