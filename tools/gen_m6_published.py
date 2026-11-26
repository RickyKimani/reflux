# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["scipy==1.14.1", "numpy==2.1.3"]
# ///
"""Write M6's published data (step 10, criterion 8): testdata/m6_published.json.

Every data set is from the NIST TRC ThermoML Archive (doi:10.18434/mds2-2422),
which holds the data five journals published from 2003 to 2019, each file
fetched by its paper's DOI from https://trc.nist.gov/ThermoML/ and its
numbers kept as printed:

- methanol-water and ethanol-water VLE at 101.3 kPa: the binary sets that
  papers on other systems measured to check their stills, T-x and y-x joined
  by x;
- water-1-butanol: its mutual solubilities at 303.15 and 313.15 K, and its
  heteroazeotrope at 101.3 kPa, a vapour over two liquids (Gomis et al.'s
  papers);
- water + 1-butanol + p-xylene: LLE tie lines at 313.15 K (Gomis et al.,
  J. Chem. Thermodyn. 79 (2014) 242).

The archive's terms: NIST's open licence, NIST acknowledged as the source;
the files shared with the journals' permission; each value as the authors
gave it, with its uncertainty where the file has one (none of these do).

Each VLE set is given the point test for Gibbs-Duhem consistency (Fredenslund,
Gmehling and Rasmussen, 1977, after Van Ness, Byer and Gibbs): G^E/RT as a
four-term Legendre series in x1 - x2, fitted to the bubble pressures by
Barker's method with the vapour the ideal gas and each compound's vapour
pressure by its ChemSep correlation (data/components, equation 101), then
the mean |y_model - y| against 0.01. Barker's fit can have two minima where
a set covers part of the range, so it starts from zero and from the G^E the
measured vapour implies (modified Raoult's law), keeping the lower sum of
squares. Here by SciPy, independently of
fit/consistency.v, which components/published_test.v holds to it.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_published.py

The files are cached in $REFLUX_THERMOML_CACHE, or the system's temporary
directory. It reproduces its output byte for byte on the same day.
"""
import datetime
import json
import math
import os
import pathlib
import tempfile
import tomllib
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
from scipy.optimize import least_squares

ROOT = pathlib.Path(__file__).resolve().parent.parent
KEYS = {
    "OKKJLVBELUTLKV-UHFFFAOYSA-N": ("methanol", "67-56-1"),
    "LFQSCWFLJHTTHZ-UHFFFAOYSA-N": ("ethanol", "64-17-5"),
    "XLYOFNOQVPJJNP-UHFFFAOYSA-N": ("water", "7732-18-5"),
    "LRHPLDYGYMQRHN-UHFFFAOYSA-N": ("1-butanol", "71-36-3"),
    "URLKBWYHVLBVBO-UHFFFAOYSA-N": ("p-xylene", "106-42-3"),
}
# (DOI, the T-x set, the y-x set), each a binary of the alcohol and water at
# 101.3 kPa, chosen by scanning the archive's files of both compounds.
VLE = {
    "methanol": [
        ("10.1021/je200341c", 2, 3),
        ("10.1021/acs.jced.5b00849", 9, 8),
        ("10.1016/j.jct.2019.01.004", 5, 6),
        ("10.1016/j.jct.2019.05.005", 6, 7),
        ("10.1021/acs.jced.9b00306", 6, 7),
    ],
    "ethanol": [
        ("10.1021/je2005209", 1, 2),
        ("10.1021/je3007138", 1, 2),
        ("10.1016/j.jct.2013.08.020", 8, 9),
        ("10.1021/je2008704", 3, 4),
        ("10.1016/j.fluid.2017.03.010", 6, 7),
        ("10.1016/j.fluid.2014.07.022", 0, 1),
    ],
}
GOMIS_2014 = "10.1016/j.jct.2014.07.021"
# Water-1-butanol's mutual solubilities: (DOI, the organic phase's set, the
# aqueous phase's set), each a phase's water or 1-butanol fraction against T.
LLE = [
    ("10.1016/j.jct.2018.03.024", 4, 3),
    ("10.1021/je400501z", 3, 4),
    (GOMIS_2014, 4, 3),
]


def local(tag):
    return tag.split("}", 1)[-1]


def child(e, name):
    return [c for c in e if local(c.tag) == name]


def text(e, path):
    for name in path:
        found = child(e, name)
        if not found:
            return None
        e = found[0]
    return e.text


def fetch(doi):
    cache = pathlib.Path(os.environ.get("REFLUX_THERMOML_CACHE", pathlib.Path(tempfile.gettempdir()) / "reflux_thermoml"))
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / (doi.replace("/", "_") + ".xml")
    if not f.exists():
        req = urllib.request.Request(f"https://trc.nist.gov/ThermoML/{doi}.xml", headers={"User-Agent": "curl/8.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            f.write_bytes(r.read())
    return ET.parse(f).getroot()


def citation(root):
    c = child(root, "Citation")[0]
    return {
        "authors": [a.text for a in child(c, "sAuthor")],
        "title": text(c, ["sTitle"]),
        "journal": text(c, ["sPubName"]),
        "year": text(c, ["yrPubYr"]),
        "volume": text(c, ["sVol"]),
        "pages": text(c, ["sPage"]),
        "doi": text(c, ["sDOI"]),
    }


def compounds(root):
    return {text(c, ["RegNum", "nOrgNum"]): KEYS.get(text(c, ["sStandardInChIKey"]), ("?", "?"))[0]
            for c in child(root, "Compound")}


def data_set(root, k):
    """A set's constraints, variables and properties by (name, phase,
    compound), each row's values as printed, and its uncertainties."""
    names = compounds(root)
    d = child(root, "PureOrMixtureData")[k]
    cons = {}
    for c in child(d, "Constraint"):
        ct = list(child(child(c, "ConstraintID")[0], "ConstraintType")[0])[0]
        cons[ct.text] = text(c, ["nConstraintValue"])
    cols = {}
    for v in child(d, "Variable"):
        vt = child(v, "VariableID")[0]
        t = list(child(vt, "VariableType")[0])[0].text
        cols[("v", text(v, ["nVarNumber"]))] = (t, text(v, ["VarPhaseID", "eVarPhase"]), names.get(text(vt, ["RegNum", "nOrgNum"])))
    unc = {}
    for p in child(d, "Property"):
        pm = child(p, "Property-MethodID")[0]
        grp = [x for x in pm if local(x.tag) == "PropertyGroup"][0]
        pg = list(grp)[0]
        key = (text(pg, ["ePropName"]), text(p, ["PropPhaseID", "ePropPhase"]), names.get(text(pm, ["RegNum", "nOrgNum"])))
        cols[("p", text(p, ["nPropNumber"]))] = key
        found = [x.text for x in p.iter() if local(x.tag) in ("nExpandUncertValue", "nCombExpandUncertValue")]
        if found:
            unc["/".join(str(s) for s in key)] = found[0]
    rows = []
    for row in child(d, "NumValues"):
        r = {}
        for x in child(row, "VariableValue"):
            r[cols[("v", text(x, ["nVarNumber"]))]] = text(x, ["nVarValue"])
        for x in child(row, "PropertyValue"):
            r[cols[("p", text(x, ["nPropNumber"]))]] = text(x, ["nPropValue"])
        rows.append(r)
    comps = sorted(names.get(text(c, ["RegNum", "nOrgNum"])) for c in child(d, "Component"))
    return comps, cons, rows, unc


def only(row, name, phase=None, compound=None):
    hits = [v for (n, ph, c), v in row.items() if n == name and (phase is None or ph == phase)
            and (compound is None or c == compound)]
    assert len(hits) == 1, (name, phase, compound, row)
    return hits[0]


def vle_set(alcohol, doi, t_set, y_set):
    root = fetch(doi)
    comps, cons, trows, unc = data_set(root, t_set)
    comps2, cons2, yrows, unc2 = data_set(root, y_set)
    assert comps == comps2 == sorted([alcohol, "water"]), (doi, comps, comps2)
    p = cons.get("Pressure, kPa")
    assert p is not None and abs(float(p) - 101.3) < 0.1, (doi, cons)
    assert cons2.get("Pressure, kPa") in (None, p), (doi, cons2)
    ys = {}
    for r in yrows:
        x = only(r, "Mole fraction", "Liquid", alcohol)
        ys[x] = only(r, "Mole fraction", "Gas", alcohol)
    points = []
    for r in trows:
        x = only(r, "Mole fraction", "Liquid", alcohol)
        points.append({"x1": x, "T": only(r, "Boiling temperature at pressure P, K"), "y1": ys.get(x)})
    assert all(pt["y1"] is not None for pt in points), (doi, "a T-x point without its y")
    unc.update(unc2)
    return {"system": [alcohol, "water"], "cas": [KEYS_BY_NAME[alcohol], KEYS_BY_NAME["water"]], "doi": doi,
            "citation": citation(root), "sets": [t_set, y_set], "P_kPa": p, "uncertainty": unc, "points": points}


KEYS_BY_NAME = {name: cas for name, cas in KEYS.values()}


def psat_correlation(cas):
    rec = tomllib.loads((ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8"))
    vp = rec["correlations"]["vapor_pressure"]
    assert vp["eq"] == 101 and vp["units"] == "Pa", (cas, vp)
    return vp


def psat(vp, t):
    a, b, c, d, e = vp["coeffs"]
    return math.exp(a + b / t + c * math.log(t) + d * t ** e)


def legendre(z, n):
    l, d = [0.0] * n, [0.0] * n
    if n > 0:
        l[0] = 1.0
    if n > 1:
        l[1], d[1] = z, 1.0
    for k in range(2, n):
        l[k] = ((2 * k - 1) * z * l[k - 1] - (k - 1) * l[k - 2]) / k
        d[k] = d[k - 2] + (2 * k - 1) * l[k - 1]
    return l, d


def point_test(s, terms=4):
    """Barker's method on a four-term Legendre series, then the mean |dy|."""
    vp1, vp2 = psat_correlation(s["cas"][0]), psat_correlation(s["cas"][1])
    pts = [pt for pt in s["points"] if 0.0 < float(pt["x1"]) < 1.0]
    p = float(s["P_kPa"]) * 1e3
    x = np.array([float(pt["x1"]) for pt in pts])
    t = np.array([float(pt["T"]) for pt in pts])
    y = np.array([float(pt["y1"]) for pt in pts])
    p1 = np.array([psat(vp1, ti) for ti in t])
    p2 = np.array([psat(vp2, ti) for ti in t])
    c1 = np.zeros((len(x), terms))
    c2 = np.zeros((len(x), terms))
    for i, xi in enumerate(x):
        l, d = legendre(2 * xi - 1, terms)
        for k in range(terms):
            g = xi * (1 - xi) * l[k]
            dg = (1 - 2 * xi) * l[k] + 2 * xi * (1 - xi) * d[k]
            c1[i, k] = g + (1 - xi) * dg
            c2[i, k] = g - xi * dg

    def parts(a):
        return x * np.exp(c1 @ a) * p1, (1 - x) * np.exp(c2 @ a) * p2

    def res(a):
        u, v = parts(a)
        return (u + v - p) / p

    def jac(a):
        u, v = parts(a)
        return (u[:, None] * c1 + v[:, None] * c2) / p

    fit = least_squares(res, np.zeros(terms), jac=jac, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if np.all((y > 0) & (y < 1)):
        g = x * np.log(y * p / (x * p1)) + (1 - x) * np.log((1 - y) * p / ((1 - x) * p2))
        leg = np.array([legendre(2 * xi - 1, terms)[0] for xi in x])
        start, *_ = np.linalg.lstsq(leg, g / (x * (1 - x)), rcond=None)
        other = least_squares(res, start, jac=jac, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if other.cost < fit.cost:
            fit = other
    u, v = parts(fit.x)
    dy = u / (u + v) - y
    mean = float(np.mean(np.abs(dy)))
    return {"terms": terms, "points": len(pts), "coefficients": [float(a) for a in fit.x],
            "mean_dy": mean, "max_dy": float(np.max(np.abs(dy))), "consistent": mean <= 0.01}


def lle_points():
    out = []
    for doi, organic, aqueous in LLE:
        root = fetch(doi)
        sets = []
        for k in (organic, aqueous):
            comps, cons, rows, unc = data_set(root, k)
            assert comps == sorted(["1-butanol", "water"]), (doi, comps)
            sets.append((cons, rows, unc))
        (oc, orows, ounc), (ac, arows, aunc) = sets
        for ro, ra in zip(orows, arows):
            t = only(ro, "Temperature, K")
            assert t == only(ra, "Temperature, K"), (doi, ro, ra)
            if float(t) > 360.0:
                continue  # the heteroazeotrope's liquids, in vlle below
            p = oc.get("Pressure, kPa") or only(ro, "Pressure, kPa")

            def water(r, phase):
                hits = [(c, v) for (n, ph, c), v in r.items() if n == "Mole fraction" and ph == phase]
                assert len(hits) == 1, (doi, r)
                c, v = hits[0]
                return {"water": v} if c == "water" else {"1-butanol": v}

            out.append({"system": ["water", "1-butanol"], "doi": doi, "citation": citation(root), "sets": [organic, aqueous],
                        "T": t, "P_kPa": p, "organic": water(ro, "Liquid mixture 1"), "aqueous": water(ra, "Liquid mixture 2")})
    return out


def vlle():
    root = fetch(GOMIS_2014)
    _, cons, trows, _ = data_set(root, 5)
    _, cons2, yrows, _ = data_set(root, 6)
    liquids = [only(r, "Mole fraction", "Liquid", "1-butanol") for r in trows]
    ts = {only(r, "Boiling temperature at pressure P, K") for r in trows}
    ys = {only(r, "Mole fraction", "Gas", "1-butanol") for r in yrows}
    assert len(ts) == 1 and len(ys) == 1, (ts, ys)
    return [{"system": ["water", "1-butanol"], "doi": GOMIS_2014, "citation": citation(root), "sets": [5, 6],
             "P_kPa": cons["Pressure, kPa"], "T": ts.pop(), "organic": {"1-butanol": max(liquids, key=float)},
             "aqueous": {"1-butanol": min(liquids, key=float)}, "vapour": {"1-butanol": ys.pop()}}]


def ternary():
    root = fetch(GOMIS_2014)
    comps, cons, rows, unc = data_set(root, 0)
    assert comps == sorted(["1-butanol", "p-xylene", "water"])
    lines = []
    for r in rows:
        lines.append({"organic": {"water": only(r, "Mole fraction", "Liquid mixture 1", "water"),
                                  "1-butanol": only(r, "Mole fraction", "Liquid mixture 1", "1-butanol")},
                      "aqueous": {"water": only(r, "Mole fraction", "Liquid mixture 2", "water"),
                                  "1-butanol": only(r, "Mole fraction", "Liquid mixture 2", "1-butanol")}})
    return [{"system": ["water", "1-butanol", "p-xylene"], "doi": GOMIS_2014, "citation": citation(root), "sets": [0],
             "T": cons["Temperature, K"], "P_kPa": cons["Pressure, kPa"], "tie_lines": lines}]


def main():
    vle = []
    for alcohol, specs in VLE.items():
        for doi, ts, ys in specs:
            s = vle_set(alcohol, doi, ts, ys)
            s["point_test"] = point_test(s)
            vle.append(s)
    out = {
        "source": {"generator": "tools/gen_m6_published.py", "archive": "NIST TRC ThermoML Archive, doi:10.18434/mds2-2422",
                   "terms": "NIST's open licence, NIST acknowledged as the source; shared with the journals' permission",
                   "generated": datetime.date.today().isoformat()},
        "note": ("Each value as printed in the archive's file, a string; P in kPa and T in K. VLE: x1 and y1 are the "
                 "alcohol's mole fractions, T the bubble point at P. LLE: each liquid's water or 1-butanol fraction, "
                 "organic and aqueous; the ternary's p-xylene is the rest. point_test: Fredenslund et al.'s, by SciPy, "
                 "on the ChemSep vapour pressures in vapour_pressure."),
        "vapour_pressure": {cas: psat_correlation(cas) for cas in sorted(KEYS_BY_NAME.values())},
        "vle": vle,
        "lle": lle_points(),
        "vlle": vlle(),
        "ternary_lle": ternary(),
    }
    (ROOT / "testdata" / "m6_published.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n",
                                                        encoding="utf-8", newline="\n")
    for s in vle:
        pt = s["point_test"]
        print(s["system"][0], s["doi"], len(s["points"]), "points; mean |dy|", round(pt["mean_dy"], 5), "max",
              round(pt["max_dy"], 5), "consistent" if pt["consistent"] else "NOT consistent")
    print(len(out["lle"]), "LLE points,", len(out["vlle"]), "VLLE,", len(out["ternary_lle"][0]["tie_lines"]), "tie lines")


if __name__ == "__main__":
    main()
