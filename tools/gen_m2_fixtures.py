# /// script
# requires-python = ">=3.11"
# dependencies = ["chemicals==1.5.2", "mpmath==1.3.0", "CoolProp==8.0.0"]
# ///
"""Generate the M2 reference fixtures for reflux.

Run from the repository root:

    uv run tools/gen_m2_fixtures.py

It writes two files to testdata/:

m2_tdep.json
    Every equation form tdep supports, read by tdep/fixtures_test.v. Real
    coefficient sets come from the vendored ChemSep XML, plus made-up sets
    for the forms ChemSep never uses and for a few hard cases. Each is
    evaluated at 40 significant digits with mpmath: the value, dY/dT, ∫Y dT
    and ∫Y/T dT. The formulas below are written independently of reflux's V
    code. For the DIPPR forms they're also checked against Caleb Bell's
    `chemicals`: if a formula here disagrees, the script stops instead of
    writing a bad fixture.

m2_coolprop.json
    ChemSep's data against CoolProp's reference equations, for every
    CoolProp fluid that ChemSep also has (matched by CAS number). It holds
    Tc, Pc, ω and M, the ideal-gas Cp and the vapour pressure, and the
    deviation of ChemSep's values from CoolProp's at the time of generation.
    ChemSep's values are read from data/components/ and evaluated here, with
    the formulas below. components/coolprop_test.v recomputes the deviations
    through reflux and requires that they match and haven't grown. That's a
    regression guard on the data and its conversion, not a claim that
    ChemSep is as accurate as CoolProp.
"""

import datetime
import json
import math
import pathlib
import platform
import re
import tomllib
import xml.etree.ElementTree as ET

import chemicals
import chemicals.dippr as dippr
import CoolProp
import CoolProp.CoolProp as CP
import mpmath as mp

mp.mp.dps = 40
ROOT = pathlib.Path(__file__).resolve().parent.parent
VENDOR = ROOT / "data" / "vendor" / "chemsep"
COMPONENTS = ROOT / "data" / "components"
OUT = ROOT / "testdata" / "m2_tdep.json"
OUT_COOLPROP = ROOT / "testdata" / "m2_coolprop.json"

TAU_FORMS = {106, 114, 116}


def Y(eq, k, t, tc):
    """Y(T) for equation `eq` with coefficients k = [A..E], in mpmath."""
    A, B, C, D, E = (mp.mpf(x) for x in k)
    t = mp.mpf(t)
    if eq in (1, 2, 3, 4, 5):
        return sum(c * t**i for i, c in enumerate([A, B, C, D, E][:eq]))
    if eq == 10:
        return mp.exp(A - B / (t + C))
    if eq in (11, 12, 13, 14, 15):
        return mp.exp(sum(c * t**i for i, c in enumerate([A, B, C, D, E][: eq - 10])))
    if eq == 16:
        return A + mp.exp(B / t + C + D * t + E * t**2)
    if eq == 17:
        return A + mp.exp(B + C * t + D * t**2 + E * t**3)
    if eq == 45:
        return A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 + E * t**5 / 5
    if eq == 100:
        return A + B * t + C * t**2 + D * t**3 + E * t**4
    if eq == 101:
        return mp.exp(A + B / t + C * mp.log(t) + D * t**E)
    if eq == 102:
        return A * t**B / (1 + C / t + D / t**2)
    if eq == 104:
        return A + B / t + C / t**3 + D / t**8 + E / t**9
    if eq == 105:
        return A / B ** (1 + (1 - t / C) ** D)
    if eq == 106:
        tr = t / mp.mpf(tc)
        return A * (1 - tr) ** (B + C * tr + D * tr**2 + E * tr**3)
    if eq == 107:
        x, y = C / t, E / t
        return A + B * (x / mp.sinh(x)) ** 2 + D * (y / mp.cosh(y)) ** 2
    if eq == 114:
        tau = 1 - t / mp.mpf(tc)
        return A**2 / tau + B - 2 * A * C * tau - A * D * tau**2 - C**2 * tau**3 / 3 - C * D * tau**4 / 2 - D**2 * tau**5 / 5
    if eq == 115:
        return mp.exp(A + B / t + C * mp.log(t) + D * t**2 + E / t**2)
    if eq == 116:
        tau = 1 - t / mp.mpf(tc)
        return A + B * tau ** mp.mpf("0.35") + C * tau ** (mp.mpf(2) / 3) + D * tau + E * tau ** (mp.mpf(4) / 3)
    raise ValueError(eq)


def chemicals_value(eq, k, t, tc, order):
    """The same quantity from chemicals.dippr, for the DIPPR forms. None for
    the other forms, and where chemicals overflows (Aly–Lee far outside its
    range, the case reflux's evaluator is written to survive)."""
    try:
        return _chemicals_value(eq, k, t, tc, order)
    except OverflowError:
        return None


def _chemicals_value(eq, k, t, tc, order):
    A, B, C, D, E = k
    if eq == 100:
        return dippr.EQ100(t, A, B, C, D, E, order=order)
    if eq == 101:
        return dippr.EQ101(t, A, B, C, D, E, order=order)
    if eq == 102:
        return dippr.EQ102(t, A, B, C, D, order=order)
    if eq == 104:
        return dippr.EQ104(t, A, B, C, D, E, order=order)
    if eq == 105:
        return dippr.EQ105(t, A, B, C, D, order=order)
    if eq == 106:
        return dippr.EQ106(t, tc, A, B, C, D, E, order=order)
    if eq == 107:
        return dippr.EQ107(t, A, B, C, D, E, order=order)
    if eq == 114:
        return dippr.EQ114(t, tc, A, B, C, D, order=order)
    if eq == 115:
        return dippr.EQ115(t, A, B, C, D, E, order=order)
    if eq == 116:
        return dippr.EQ116(t, tc, A, B, C, D, E, order=order)
    return None


USED = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 10: 3, 11: 1, 12: 2, 13: 3, 14: 4, 15: 5, 102: 4, 105: 4, 114: 4}


def chemsep_sets():
    """Real coefficient sets, up to three per form. Records that use more of
    their coefficients come first, then records of different properties
    are preferred. The tau forms carry the compound's Tc."""
    bad = re.compile(r'value="([^"]*[<&][^"]*)"')
    candidates = {}
    forms = set()
    for fname in ("chemsep1.xml", "chemsep2.xml"):
        raw = (VENDOR / fname).read_text(encoding="latin-1")
        raw = bad.sub(lambda m: 'value="' + m.group(1).replace("&", "&amp;").replace("<", "&lt;") + '"', raw)
        for comp in ET.fromstring(raw).findall("compound"):
            name = comp.find("CompoundID").get("value")
            tc_el = comp.find("CriticalTemperature")
            tc = float(tc_el.get("value")) if tc_el is not None else 0.0
            for el in comp:
                eq_el = el.find("eqno")
                if eq_el is None:
                    continue
                eq = int(eq_el.get("value"))
                forms.add(eq)
                k = [float(el.find(c).get("value")) if el.find(c) is not None else 0.0 for c in "ABCDE"]
                tmin = float(el.find("Tmin").get("value"))
                tmax = float(el.find("Tmax").get("value"))
                # 114 and 116 aren't defined above Tc, and reflux refuses such
                # ranges. ChemSep's 106 records often end a little past Tc;
                # the fixture points stop well short of it (see evaluate).
                if not (tmin < tmax and (eq not in (114, 116) or tmax <= tc)):
                    continue
                candidates.setdefault(eq, []).append({"eq": eq, "coeffs": k, "tmin": tmin, "tmax": tmax,
                                                      "tc": tc if eq in TAU_FORMS else 0.0,
                                                      "origin": f"{name}/{el.tag}"})
    found = {}
    for eq, cs in candidates.items():
        cs.sort(key=lambda c: -sum(1 for x in c["coeffs"][: USED.get(eq, 5)] if x != 0))
        chosen = []
        for c in cs:  # one per property first
            if len(chosen) < 3 and all(s["origin"].split("/")[1] != c["origin"].split("/")[1] for s in chosen):
                chosen.append(c)
        for c in cs:  # then fill up with other compounds
            if len(chosen) < 3 and c not in chosen:
                chosen.append(c)
        # The records reflux leans on most, whatever the selection above did.
        chosen += [c for c in cs if c["origin"] in PINNED and c not in chosen]
        found[eq] = chosen
    return found, forms


PINNED = {
    "N-butane/VaporPressure",
    "N-butane/IdealGasHeatCapacityCp",
    "N-butane/LiquidDensity",
    "N-butane/HeatOfVaporization",
    "N-butane/LiquidHeatCapacityCp",
    "Water/VaporPressure",
    "Ethanol/IdealGasHeatCapacityCp",
}


# Made-up but plausible sets: for forms ChemSep 8.50 never uses, for forms
# whose only real records leave coefficients at zero, and one Aly–Lee set
# run far outside its range, where sinh and cosh overflow a float.
SYNTHETIC = [
    {"eq": 11, "coeffs": [2.5, 0, 0, 0, 0], "tmin": 200.0, "tmax": 600.0, "tc": 0.0},
    {"eq": 14, "coeffs": [-3.0, 0.02, -4e-5, 2e-8, 0], "tmin": 250.0, "tmax": 550.0, "tc": 0.0},
    {"eq": 15, "coeffs": [-3.0, 0.02, -4e-5, 2e-8, -1e-11], "tmin": 250.0, "tmax": 550.0, "tc": 0.0},
    {"eq": 17, "coeffs": [3.0e4, 9.0, 1e-3, -2e-6, 1e-9], "tmin": 300.0, "tmax": 900.0, "tc": 0.0},
    {"eq": 45, "coeffs": [3.0e4, 50.0, 0.1, -1e-4, 3e-8], "tmin": 250.0, "tmax": 800.0, "tc": 0.0},
    {"eq": 115, "coeffs": [-10.0, 1000.0, 0.5, -1e-6, 1e4], "tmin": 250.0, "tmax": 600.0, "tc": 0.0},
    # In ChemSep's second virial coefficients the D/T⁸ and E/T⁹ terms are
    # about 1e-14 of the total, too small for a test to see. Here every term
    # counts.
    {"eq": 104, "coeffs": [0.5, -1.0, -0.3, 0.05, -0.01], "tmin": 1.0, "tmax": 3.0, "tc": 0.0},
    {"eq": 107, "coeffs": [33298.3, 26618.3, 2572.35, 8630.72, 1138.84], "tmin": 2.0, "tmax": 1e5, "tc": 0.0,
     "ts": [2.0, 10.0, 100.0, 1000.0, 1e4, 1e5]},
]

SUPPORTED = [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 45, 100, 101, 102, 104, 105, 106, 107, 114, 115, 116]
CHEMICALS_INTEGRAL = {100, 104, 107, 114, 116}


def check(what, case, t, ours, theirs, rtol, scale=None):
    """Stop if this script and chemicals disagree, rather than write a bad fixture."""
    scale = abs(ours) if scale is None else scale
    if abs(ours - theirs) > rtol * max(scale, 1e-300):
        raise SystemExit(f"eq {case['eq']} ({case['origin']}) {what} at T={t}: mpmath {ours!r} vs chemicals {theirs!r}")


def evaluate(case):
    eq, k, tc = case["eq"], case["coeffs"], case["tc"]
    lo, hi = case["tmin"], case["tmax"]
    ts = case.pop("ts", None) or [lo + f * (hi - lo) for f in (0.05, 0.35, 0.65, 0.95)]
    if eq in TAU_FORMS and max(ts) >= tc:
        raise SystemExit(f"eq {eq} ({case['origin']}): fixture point at or above Tc = {tc}")
    f = lambda t: Y(eq, k, t, tc)
    g = lambda t: Y(eq, k, t, tc) / t
    t0 = ts[0]
    points = []
    integ = integ_t = mp.mpf(0)
    for i, t in enumerate(ts):
        if i > 0:
            # Piecewise between the points, so each piece is short.
            a = mp.mpf(ts[i - 1])
            d, e1 = mp.quad(f, [a, t], error=True)
            dt, e2 = mp.quad(g, [a, t], error=True)
            if e1 > 1e-25 * abs(d) or e2 > 1e-25 * abs(dt):
                raise SystemExit(f"eq {eq} ({case['origin']}): mpmath quadrature unsure on [{a}, {t}]")
            integ += d
            integ_t += dt
        y = f(t)
        dy = mp.diff(f, mp.mpf(t))
        p = {"t": t, "y": float(y), "dy": float(dy), "integral": float(integ), "integral_over_t": float(integ_t)}
        cy = chemicals_value(eq, k, t, tc, 0)
        cdy = chemicals_value(eq, k, t, tc, 1)
        if cy is not None and cdy is not None:
            check("value", case, t, float(y), cy, 1e-12)
            check("derivative", case, t, float(dy), cdy, 1e-9)
            if eq in CHEMICALS_INTEGRAL:
                hi_, lo_ = chemicals_value(eq, k, t, tc, -1), chemicals_value(eq, k, t0, tc, -1)
                if hi_ is not None and lo_ is not None:
                    check("integral", case, t, float(integ), hi_ - lo_, 1e-10,
                          max(abs(hi_), abs(lo_), abs(float(integ))))
            p["chemicals_y"] = cy
            p["chemicals_dy"] = cdy
        points.append(p)
    return {**case, "t0": t0, "points": points}


def usable(c):
    """Whether reflux can use a correlation record: the rules of
    components/build.v, restated."""
    if c is None or c["eq"] not in SUPPORTED or not c["coeffs"]:
        return False
    lo, hi = c["t_range"]
    return not (lo > hi or lo < 0 or hi <= 0)


def chemsep_value(c, t, tc, scale):
    """A ChemSep correlation record's value at t, in SI per mol."""
    k = (list(c["coeffs"]) + [0.0] * 5)[:5]
    return float(Y(c["eq"], k, t, c.get("tc", tc))) * scale


def points(lo, hi, value, reference):
    """Four temperatures across [lo, hi], with CoolProp's value and ChemSep's
    relative deviation from it. Points CoolProp can't evaluate are skipped."""
    out = []
    for f in (0.05, 0.35, 0.65, 0.95):
        t = lo + f * (hi - lo)
        try:
            ref = reference(t)
        except ValueError:
            continue
        out.append({"t": t, "coolprop": ref, "deviation": (value(t) - ref) / ref})
    return out


def coolprop_fluid(name, cas, rec):
    c = rec.get("constants", {})
    cor = rec.get("correlations", {})
    fluid = {"coolprop": name, "cas": cas, "name": rec["name"], "constants": {}}
    for key, param in (("tc", "Tcrit"), ("pc", "pcrit"), ("omega", "acentric"), ("mw", "molar_mass")):
        if key not in c:
            continue
        try:
            ref = CP.PropsSI(param, name)
        except ValueError:
            continue
        # ω can be near zero (argon's is −0.002), so its deviation is absolute.
        dev = c[key] - ref if key == "omega" else (c[key] - ref) / ref
        fluid["constants"][key] = {"coolprop": ref, "deviation": dev}
    t_low = max(CP.PropsSI("Ttriple", name), CP.PropsSI("Tmin", name))
    t_high = CP.PropsSI("Tmax", name)
    tc_ref = CP.PropsSI("Tcrit", name)
    cp = cor.get("ideal_gas_cp")
    if usable(cp):
        assert cp["units"] == "J/kmol/K", f"{cas} ideal_gas_cp in {cp['units']}"
        lo, hi = max(cp["t_range"][0], t_low), min(cp["t_range"][1], t_high)
        if lo < hi:
            fluid["cp0"] = {"key": "ideal_gas_cp", "points": points(lo, hi,
                lambda t: chemsep_value(cp, t, c.get("tc", 0.0), 1e-3),
                lambda t: CP.PropsSI("Cp0molar", "T", t, "Dmolar", 1e-3, name))}
    key = "vapor_pressure" if usable(cor.get("vapor_pressure")) else "antoine_vapor_pressure"
    vp = cor.get(key)
    if usable(vp):
        assert vp["units"] == "Pa", f"{cas} {key} in {vp['units']}"
        lo, hi = max(vp["t_range"][0], t_low), min(vp["t_range"][1], 0.98 * tc_ref)
        if lo < hi:
            fluid["psat"] = {"key": key, "points": points(lo, hi,
                lambda t: chemsep_value(vp, t, c.get("tc", 0.0), 1.0),
                lambda t: CP.PropsSI("P", "T", t, "Q", 0, name))}
    return fluid


def n_butane():
    """n-butane from CoolProp's reference equation, for the end-to-end test
    (M2 exit criterion 9): Psat at four temperatures, and the molar enthalpy
    at states whose differences are heating, vaporization, compression and a
    supercritical step. CoolProp's enthalpy zero is its own, so only the
    differences mean anything."""
    name = "n-Butane"
    psat = [{"t": t, "p": CP.PropsSI("P", "T", t, "Q", 0, name)} for t in (250.0, 300.0, 350.0, 400.0)]
    states = []
    for label, phase, t, p in (
        ("vapor 300 K 1 bar", "vapor", 300.0, 1e5),
        ("vapor 450 K 1 bar", "vapor", 450.0, 1e5),
        ("liquid 250 K 20 bar", "liquid", 250.0, 20e5),
        ("liquid 350 K 20 bar", "liquid", 350.0, 20e5),
        ("liquid 300 K 5 bar", "liquid", 300.0, 5e5),
        ("liquid 300 K 100 bar", "liquid", 300.0, 100e5),
        ("supercritical 450 K 50 bar", "vapor", 450.0, 50e5),
        ("supercritical 500 K 50 bar", "vapor", 500.0, 50e5),
    ):
        states.append({"label": label, "phase": phase, "t": t, "p": p,
                       "h": CP.PropsSI("Hmolar", "T", t, "P", p, name)})
    p300 = CP.PropsSI("P", "T", 300.0, "Q", 0, name)
    for q, phase in ((0, "liquid"), (1, "vapor")):
        states.append({"label": f"saturated {phase} 300 K", "phase": phase, "t": 300.0, "p": p300,
                       "h": CP.PropsSI("Hmolar", "T", 300.0, "Q", q, name)})
    return {"psat": psat, "states": states}


def coolprop_fixture():
    fluids = []
    for name in sorted(CP.get_global_param_string("fluids_list").split(","), key=str.lower):
        cas = CP.get_fluid_param_string(name, "CAS")
        path = COMPONENTS / f"{cas}.toml"
        if path.exists():
            fluids.append(coolprop_fluid(name, cas, tomllib.loads(path.read_text(encoding="utf-8"))))
    return {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m2_fixtures.py",
            "python": platform.python_version(),
            "coolprop": CoolProp.__version__,
            "chemsep": "data/components (ChemSep 8.50, converted)",
        },
        "note": "coolprop is CoolProp's value (HEOS backend); deviation is ChemSep minus CoolProp, "
        "relative, except for omega, where it's absolute. cp0 in J/(mol K) at low density, psat in Pa "
        "on the saturated liquid. key names the ChemSep correlation used. n_butane holds CoolProp's "
        "Psat (Pa) and molar enthalpy h (J/mol, CoolProp's own zero) for the end-to-end test.",
        "fluids": fluids,
        "n_butane": n_butane(),
    }


def summarize(data):
    """Print the worst deviation of each kind, for the record."""
    worst = {}
    for f in data["fluids"]:
        for k, v in f["constants"].items():
            worst.setdefault(k, []).append((abs(v["deviation"]), f["coolprop"]))
        for k in ("cp0", "psat"):
            for p in f.get(k, {}).get("points", []):
                worst.setdefault(k, []).append((abs(p["deviation"]), f"{f['coolprop']} at {p['t']:.1f} K"))
    for k, v in worst.items():
        v.sort(reverse=True)
        median = sorted(x for x, _ in v)[len(v) // 2]
        print(f"  {k:6} median {median:.2e}, worst: " + "; ".join(f"{w} {d:.2e}" for d, w in v[:3]))


def main():
    found, forms = chemsep_sets()
    cases = []
    for eq in SUPPORTED:
        for s in found.get(eq, []):
            cases.append(evaluate(s))
    for s in SYNTHETIC:
        cases.append(evaluate({**s, "origin": "synthetic"}))
    unsupported = sorted(e for e in forms if e not in SUPPORTED)
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m2_fixtures.py",
            "python": platform.python_version(),
            "chemicals": chemicals.__version__,
            "mpmath": mp.__version__,
            "chemsep": "data/vendor/chemsep (ChemSep 8.50)",
        },
        "note": "y, dy, integral and integral_over_t are mpmath at 40 digits; the integrals run from t0 "
        "to t. chemicals_y and chemicals_dy are chemicals.dippr, for the DIPPR forms.",
        "unsupported_forms_in_data": unsupported,
        "cases": cases,
    }
    OUT.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    forms = sorted({c["eq"] for c in cases})
    print(f"wrote {OUT.relative_to(ROOT)}: {len(cases)} cases covering forms {forms}")
    print(f"forms in the data that tdep doesn't support: {unsupported}")

    cp_data = coolprop_fixture()
    OUT_COOLPROP.write_text(json.dumps(cp_data, indent=1) + "\n", encoding="utf-8", newline="\n")
    n = len(cp_data["fluids"])
    print(f"wrote {OUT_COOLPROP.relative_to(ROOT)}: {n} CoolProp fluids also in ChemSep")
    summarize(cp_data)


if __name__ == "__main__":
    main()
