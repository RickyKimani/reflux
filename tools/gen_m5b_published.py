# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""Write M5b's measured mixtures (step 9, criterion 5): testdata/m5b_published.json.

Every data set is from the NIST TRC ThermoML Archive (doi:10.18434/mds2-2422),
which holds the data five journals published from 2003 to 2019, each file
fetched by its paper's DOI from https://trc.nist.gov/ThermoML/ and its points
kept as printed: T, P (101.325 kPa where a set gives none), the composition
(a mass fraction turned to mole fractions by the component files' molar
masses), the value in SI, and its expanded uncertainty as the file gives it.

- gas viscosity: H2 + CO2 at 300 kPa from 500 to 1100 K (Maltsev et al.,
  2004); CO2 + CH4 to 32 MPa (Locke et al., 2015); CH4 + C3H8 to 31 MPa
  (Stanwix et al., 2015); CH4 + C4H10 to 30 MPa (Locke et al., 2015);
- gas conductivity: CO2, and CH4 + CO2, to 15 MPa (Patek et al., 2005);
- liquid viscosity: heptane, octane and their mixtures to 10 MPa and 473 K
  (Abdulagatov et al., 2006) and from 3 to 5 MPa to 513 K (Yang et al.,
  2015); ethanol, and ethanol + heptane, at 1 atm and to 100 MPa
  (Zéberg-Mikkelsen et al., 2006); 1-butanol + 2-methylheptane
  (Torín-Ollarves et al., 2012); methanol
  + water (Semenov et al., 2019); diethylene glycol, and
  its mixtures with water (Begum et al., 2011); diethylethanolamine + water
  (Lampreia et al., 2007);
- liquid conductivity: 2-butoxyethanol, and its mixtures with water, to 150
  MPa (Mensah-Brown and Wakeham, 2010);
- surface tension: CH4 + C2H6 (Baidakov et al., 2013); O2, N2 and N2 + O2
  (Baidakov et al., 2008); diethylethanolamine + water (Lampreia et al.,
  2007).

The archive's terms: NIST's open licence, NIST acknowledged as the source;
the files shared with the journals' permission; each value as the authors
gave it.

For each set, the alternatives to M5b's methods (D4 to D8), each on the
same pure values, ChemSep's correlations from reflux's component files,
with `chemicals`' functions where it has them; a gas's density by
`thermo`'s Peng-Robinson with every kij 0, its lower-Gibbs root. Each
alternative's average absolute deviation, mean deviation and largest
deviation, as fractions of the measured values, and the points it gives
no value at. M5b's own are worked out in V, by transport/published_test.v,
on reflux's own Peng-Robinson, ChemSep's kij and the phase's label's root.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5b_published.py

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

from chemicals.interface import Diguilio_Teja
from chemicals.thermal_conductivity import DIPPR9H, DIPPR9I, Filippov, Lindsay_Bromley, Missenard, Stiel_Thodos_dense
from chemicals.viscosity import Herning_Zipperer, Lorentz_Bray_Clarke, Wilke
from thermo.eos_mix import PRMIX

ROOT = pathlib.Path(__file__).resolve().parent.parent
R = 8.31446261815324
# ThermoML's common names, to the component files' CAS numbers.
NAMES = {"methane": "74-82-8", "carbon dioxide": "124-38-9", "propane": "74-98-6", "butane": "106-97-8",
         "hydrogen": "1333-74-0", "methanol": "67-56-1", "water": "7732-18-5", "diethylene glycol": "111-46-6",
         "nitrogen": "7727-37-9", "oxygen": "7782-44-7", "ethane": "74-84-0", "2-butoxyethan-1-ol": "111-76-2",
         "(diethylamino)ethanol": "100-37-8", "heptane": "142-82-5", "octane": "111-65-9", "ethanol": "64-17-5",
         "butan-1-ol": "71-36-3", "2-methylheptane": "592-27-8"}

# Each set: its key in the tests, its paper's DOI, its place in the file,
# the property and the phase measured.
SETS = [
    ("H2 + CO2 gas viscosity, 300 kPa", "10.1021/je0342419", 2, "viscosity", "vapour"),
    ("CO2 + CH4 gas viscosity, to 32 MPa", "10.1016/j.jct.2015.03.007", 0, "viscosity", "vapour"),
    ("CH4 + C3H8 gas viscosity, to 31 MPa", "10.1021/je500857w", 0, "viscosity", "vapour"),
    ("CH4 + C4H10 gas viscosity, to 30 MPa", "10.1021/acs.jced.5b00635", 0, "viscosity", "vapour"),
    ("CO2 gas conductivity, to 15 MPa", "10.1007/s10765-005-5566-6", 0, "thermal_conductivity", "vapour"),
    ("CH4 + CO2 gas conductivity, to 12 MPa", "10.1007/s10765-005-5566-6", 1, "thermal_conductivity",
     "vapour"),
    ("heptane liquid viscosity, to 10 MPa", "10.1016/j.jct.2006.01.012", 2, "viscosity", "liquid"),
    ("octane liquid viscosity, to 10 MPa", "10.1016/j.jct.2006.01.012", 0, "viscosity", "liquid"),
    ("heptane + octane liquid viscosity, to 10 MPa", "10.1016/j.jct.2006.01.012", 4, "viscosity", "liquid"),
    ("heptane liquid viscosity, to 513 K", "10.1016/j.tca.2015.08.005", 0, "viscosity", "liquid"),
    ("octane liquid viscosity, to 513 K", "10.1016/j.tca.2015.08.005", 1, "viscosity", "liquid"),
    ("heptane + octane liquid viscosity, to 513 K", "10.1016/j.tca.2015.08.005", 2, "viscosity", "liquid"),
    ("ethanol liquid viscosity, to pressure", "10.1016/j.fluid.2006.01.030", 0, "viscosity", "liquid"),
    ("ethanol + heptane liquid viscosity, 1 atm", "10.1016/j.fluid.2006.01.030", 4, "viscosity", "liquid"),
    ("ethanol + heptane liquid viscosity, to pressure", "10.1016/j.fluid.2006.01.030", 5, "viscosity",
     "liquid"),
    ("1-butanol + 2-methylheptane liquid viscosity", "10.1016/j.jct.2011.08.012", 4, "viscosity", "liquid"),
    ("methanol + water liquid viscosity", "10.1016/j.jct.2019.05.013", 2, "viscosity", "liquid"),
    ("diethylene glycol liquid viscosity", "10.1021/je1009976", 0, "viscosity", "liquid"),
    ("diethylene glycol + water liquid viscosity", "10.1021/je1009976", 2, "viscosity", "liquid"),
    ("diethylethanolamine + water liquid viscosity", "10.1021/je700350b", 3, "viscosity", "liquid"),
    ("2-butoxyethanol liquid conductivity, to 150 MPa", "10.1021/je1007197", 0, "thermal_conductivity",
     "liquid"),
    ("2-butoxyethanol + water liquid conductivity, to 150 MPa", "10.1021/je1007197", 1,
     "thermal_conductivity", "liquid"),
    ("CH4 + C2H6 surface tension", "10.1016/j.fluid.2013.07.008", 0, "surface_tension", "liquid"),
    ("O2 surface tension", "10.1016/j.fluid.2008.06.016", 0, "surface_tension", "liquid"),
    ("N2 surface tension", "10.1016/j.fluid.2008.06.016", 1, "surface_tension", "liquid"),
    ("N2 + O2 surface tension", "10.1016/j.fluid.2008.06.016", 2, "surface_tension", "liquid"),
    ("diethylethanolamine + water surface tension", "10.1021/je700350b", 1, "surface_tension", "liquid"),
]


# --- reflux's component files, ChemSep's correlations ---

def record(cas):
    return tomllib.loads((ROOT / "data" / "components" / f"{cas}.toml").read_text(encoding="utf-8"))


def correlation(rec, key, t):
    d = rec["correlations"][key]
    a, b, c, dd, e = (d["coeffs"] + [0.0] * 5)[:5]
    eq = d["eq"]
    tc = d.get("tc", rec["constants"]["tc"])
    if eq in (1, 2, 3, 4, 5, 100):
        return a + t * (b + t * (c + t * (dd + t * e)))
    if eq == 10:
        return math.exp(a - b / (t + c))
    if eq == 16:
        return a + math.exp(b / t + c + dd * t + e * t * t)
    if eq == 101:
        return math.exp(a + b / t + c * math.log(t) + dd * t ** e)
    if eq == 102:
        return a * t ** b / (1 + c / t + dd / (t * t))
    if eq == 105:
        return a / b ** (1 + (1 - t / c) ** dd)
    if eq == 106:
        tr = t / tc
        return a * (1 - tr) ** (b + tr * (c + tr * (dd + tr * e)))
    if eq == 115:
        return math.exp(a + b / t + c * math.log(t) + dd * t * t + e / (t * t))
    raise ValueError(f"equation {eq}")


def scale(rec, key):
    u = rec["correlations"][key]["units"]
    return {"kmol/m3": 1e3, "Pa.s": 1.0, "W/m/K": 1.0, "N/m": 1.0, "Pa": 1.0}[u]


class Pure:
    def __init__(self, cas):
        r = record(cas)
        k = r["constants"]
        self.r, self.name = r, r["name"]
        self.tc, self.pc, self.vc, self.zc, self.w = k["tc"], k["pc"], k["vc"], k["zc"], k["omega"]
        self.mw = k["mw"]  # kg/mol
        self.m = k["mw"] * 1e3  # g/mol
        self.tb = k.get("tb")

    def get(self, key, t):
        return scale(self.r, key) * correlation(self.r, key, t)

    def vl(self, t):  # m3/mol, the liquid-density correlation's
        return 1.0 / self.get("liquid_density", t)


# --- the ThermoML files ---

def local(t):
    return t.split("}", 1)[-1]


def child(e, n):
    return [c for c in e if local(c.tag) == n]


def text(e, path):
    for n in path:
        f = child(e, n)
        if not f:
            return None
        e = f[0]
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
    authors = [e.text for e in c.iter() if local(e.tag) == "sAuthor"]
    journal = next((e.text for e in c.iter() if local(e.tag) == "sPubName"), "")
    year = next((e.text for e in c.iter() if local(e.tag) == "yrPubYr"), "")
    vol = next((e.text for e in c.iter() if local(e.tag) == "sVol"), "")
    pages = next((e.text for e in c.iter() if local(e.tag) == "sPage"), "")
    title = text(c, ["sTitle"]) or ""
    who = authors[0].split("[")[0].strip() + (" et al." if len(authors) > 1 else "")
    return f"{who}, {journal} {vol} ({year}) {pages}: {title.strip()}"


def data_set(doi, k):
    """Set k's points: T K, P Pa, x (mole fractions in the set's component
    order), value SI, its expanded uncertainty; and its components' names
    and CAS numbers, and whether a pressure was given."""
    root = fetch(doi)
    names = {text(c, ["RegNum", "nOrgNum"]): text(c, ["sCommonName"]) for c in child(root, "Compound")}
    d = child(root, "PureOrMixtureData")[k]
    comps = [names[text(c, ["RegNum", "nOrgNum"])] for c in child(d, "Component")]
    cas = [NAMES[c] for c in comps]
    fixed = {}
    for c in child(d, "Constraint"):
        cid = child(c, "ConstraintID")[0]
        kind = list(child(cid, "ConstraintType")[0])[0].text
        fixed[(kind, names.get(text(cid, ["RegNum", "nOrgNum"])))] = float(text(c, ["nConstraintValue"]))
    cols = {}
    for v in child(d, "Variable"):
        vid = child(v, "VariableID")[0]
        kind = list(child(vid, "VariableType")[0])[0].text
        cols[text(v, ["nVarNumber"])] = (kind, names.get(text(vid, ["RegNum", "nOrgNum"])))
    points, given = [], True
    for nv in child(d, "NumValues"):
        p = dict(fixed)
        for x in child(nv, "VariableValue"):
            p[cols[text(x, ["nVarNumber"])]] = float(text(x, ["nVarValue"]))
        pv = child(nv, "PropertyValue")[0]
        val = float(text(pv, ["nPropValue"]))
        u = text(pv, ["CombinedUncertainty", "nCombExpandUncertValue"])
        t = p[("Temperature, K", None)]
        if ("Pressure, kPa", None) in p:
            pr = 1e3 * p[("Pressure, kPa", None)]
        else:
            pr, given = 101325.0, False
        points.append({"T": t, "P": pr, "x": composition(p, comps, cas), "value": val,
                       "u": float(u) if u is not None else None})
    return comps, cas, points, given


def composition(p, comps, cas):
    n = len(comps)
    if n == 1:
        return [1.0]
    for (kind, c), v in p.items():
        if kind == "Mole fraction":
            i = comps.index(c)
            x = [0.0] * n
            x[i], x[1 - i] = v, 1 - v
            return x
        if kind == "Mass fraction":
            i = comps.index(c)
            mw = [Pure(s).mw for s in cas]
            w = [0.0] * n
            w[i], w[1 - i] = v, 1 - v
            moles = [w[j] / mw[j] for j in range(n)]
            return [m / sum(moles) for m in moles]
    raise ValueError(p)


# --- the alternatives ---

def pr_volume(ps, t, p, y):
    eos = PRMIX(Tcs=[c.tc for c in ps], Pcs=[c.pc for c in ps], omegas=[c.w for c in ps], zs=y, T=t, P=p,
                kijs=[[0.0] * len(ps) for _ in ps])
    if hasattr(eos, "V_l") and hasattr(eos, "V_g"):
        return eos.V_l if eos.G_dep_l < eos.G_dep_g else eos.V_g
    return eos.V_l if hasattr(eos, "V_l") else eos.V_g


def lucas_published(ps, t, p, y):
    """Lucas's high-pressure ratio as published (Poling et al., 5th ed., sec.
    9-6), not taken against its own value at zero pressure (D5's choice)."""
    tcm = sum(yi * c.tc for yi, c in zip(y, ps))
    pcm = R * tcm * sum(yi * c.zc for yi, c in zip(y, ps)) / sum(yi * c.vc for yi, c in zip(y, ps))
    tr, pr = t / tcm, p / pcm
    fp0 = fq0 = 0.0
    for yi, c in zip(y, ps):
        dip = c.r["constants"].get("dipole", 0.0) / 3.33564e-30
        mur = 52.46 * dip ** 2 * (c.pc / 1e5) / c.tc ** 2
        tri = t / c.tc
        if mur < 0.022:
            f = 1.0
        elif mur < 0.075:
            f = 1 + 30.55 * (0.292 - c.zc) ** 1.72
        else:
            f = 1 + 30.55 * (0.292 - c.zc) ** 1.72 * abs(0.96 + 0.1 * (tri - 0.7))
        fp0 += yi * f
        if c.name == "Hydrogen":
            s = 1.0 if tri > 12 else -1.0
            fq0 += yi * 1.22 * 0.76 ** 0.15 * (1 + 0.00385 * ((tri - 12) ** 2) ** (1 / c.m) * s)
        else:
            fq0 += yi
    z1 = (0.807 * tr ** 0.618 - 0.357 * math.exp(-0.449 * tr) + 0.340 * math.exp(-4.058 * tr) + 0.018) * fp0 * fq0
    if tr <= 1.0:
        if pr > 1.0:
            return None
        alpha = 3.262 + 14.98 * pr ** 5.508
        beta = 1.390 + 5.746 * pr
        z2 = 0.600 + 0.760 * pr ** alpha + (6.990 * pr ** beta - 0.6) * (1 - tr)
    else:
        a = 1.245e-3 / tr * math.exp(5.1726 * tr ** -0.3286)
        b = a * (1.6553 * tr - 1.2723)
        c = 0.4489 / tr * math.exp(3.0578 * tr ** -37.7332)
        d = 1.7368 / tr * math.exp(2.2310 * tr ** -7.6351)
        f = 0.9425 * math.exp(-0.1853 * tr ** 0.4489)
        z2 = z1 * (1 + a * pr ** 1.3088 / (b * pr ** f + 1 / (1 + c * pr ** d)))
    yy = z2 / z1
    fp = (1 + (fp0 - 1) * yy ** -3) / fp0
    fq = (1 + (fq0 - 1) * (1 / yy - 0.007 * math.log(yy) ** 4)) / fq0
    return yy * fp * fq


def jst(eta0, ps, t, v, y):
    """Jossi, Stiel and Thodos' residual on eta0, at every density, on the
    mixture's mole-fraction-mean Tc, Pc, Vc and M (DWSIM's)."""
    tc = sum(yi * c.tc for yi, c in zip(y, ps))
    pc = sum(yi * c.pc for yi, c in zip(y, ps)) / 101325
    vc = sum(yi * c.vc for yi, c in zip(y, ps))
    m = sum(yi * c.m for yi, c in zip(y, ps))
    rr = vc / v
    xi = tc ** (1 / 6) * m ** -0.5 * pc ** (-2 / 3)
    poly = 1.0230 + 0.23364 * rr + 0.58533 * rr ** 2 - 0.40758 * rr ** 3 + 0.093324 * rr ** 4
    return eta0 + (poly ** 4 - 1) / xi * 1e-7


def yorizane(ps, y):
    n = len(ps)
    vcij = [[(ps[i].vc ** (1 / 3) + ps[j].vc ** (1 / 3)) ** 3 / 8 for j in range(n)] for i in range(n)]
    vcm = sum(y[i] * y[j] * vcij[i][j] for i in range(n) for j in range(n))
    tcm = sum(y[i] * y[j] * vcij[i][j] * math.sqrt(ps[i].tc * ps[j].tc) for i in range(n) for j in range(n)) / vcm
    zcm = 0.291 - 0.08 * sum(yi * c.w for yi, c in zip(y, ps))
    return tcm, zcm * R * tcm / vcm, vcm, zcm, sum(yi * c.m for yi, c in zip(y, ps))


def mason_saxena(ys, ks, mus, mws):
    n = len(ys)
    out = 0.0
    for i in range(n):
        den = sum(ys[j] * (1 + math.sqrt(mus[i] / mus[j]) * (mws[j] / mws[i]) ** 0.25) ** 2
                  / math.sqrt(8 * (1 + mws[i] / mws[j])) for j in range(n))
        out += ys[i] * ks[i] / den
    return out


def safe(f):
    try:
        v = f()
    except (ValueError, ZeroDivisionError, OverflowError):
        return None
    return v if v is not None and math.isfinite(v) and v > 0 else None


def gas_viscosity(ps, pt):
    t, p, y = pt["T"], pt["P"], pt["x"]
    mus = [c.get("vapor_viscosity", t) for c in ps]
    mws = [c.m for c in ps]
    w = Wilke(y, mus, mws) if len(ps) > 1 else mus[0]
    v = pr_volume(ps, t, p, y)
    mean = sum(yi * mi for yi, mi in zip(y, mus))
    f = lucas_published(ps, t, p, y)
    return {
        "Wilke, no pressure correction": w,
        "Herning-Zipperer, no pressure correction": safe(lambda: Herning_Zipperer(y, mus, mws)) if len(ps) > 1 else w,
        "Wilke times Lucas's factor as published": None if f is None else w * f,
        "Wilke plus Jossi-Stiel-Thodos at every density": jst(w, ps, t, v, y),
        "DWSIM's: the mole-fraction mean plus Jossi-Stiel-Thodos": jst(mean, ps, t, v, y),
        "chemicals' Lohrenz-Bray-Clark": safe(lambda: Lorentz_Bray_Clarke(t, p, v, y, mws, [c.tc for c in ps],
                                                                          [c.pc for c in ps], [c.vc for c in ps])),
    }


def gas_conductivity(ps, pt):
    t, p, y = pt["T"], pt["P"], pt["x"]
    ks = [c.get("vapor_thermal_conductivity", t) for c in ps]
    mus = [c.get("vapor_viscosity", t) for c in ps]
    mws = [c.m for c in ps]
    ms = mason_saxena(y, ks, mus, mws)
    lb = Lindsay_Bromley(t, y, ks, mus, [c.tb for c in ps], mws) if len(ps) > 1 else ks[0]
    v = pr_volume(ps, t, p, y)
    tcm, pcm, vcm, zcm, mm = yorizane(ps, y)
    return {
        "Mason-Saxena, no pressure correction": ms,
        "Lindsay-Bromley plus Stiel-Thodos": safe(lambda: Stiel_Thodos_dense(t, mm, tcm, pcm, vcm, zcm, v, lb)),
        "DWSIM's: the mole-fraction mean, no pressure correction": sum(yi * ki for yi, ki in zip(y, ks)),
    }


def liquid_viscosity(ps, pt):
    t, x = pt["T"], pt["x"]
    mus = [c.get("liquid_viscosity", t) for c in ps]
    mass = [xi * c.m for xi, c in zip(x, ps)]
    w = [m / sum(mass) for m in mass]
    return {
        "ln mu by mole fractions, no pressure correction": math.exp(sum(xi * math.log(m) for xi, m in zip(x, mus))),
        "ln mu by mass fractions, no pressure correction": math.exp(sum(wi * math.log(m) for wi, m in zip(w, mus))),
        "DWSIM's: the mole-fraction mean, no pressure correction": sum(xi * m for xi, m in zip(x, mus)),
    }


def liquid_conductivity(ps, pt):
    t, p, x = pt["T"], pt["P"], pt["x"]
    ks = [c.get("liquid_thermal_conductivity", t) for c in ps]
    kp = [Missenard(t, p, c.tc, c.pc, kk) for c, kk in zip(ps, ks)]
    vms = [c.vl(t) for c in ps]
    mass = [xi * c.m for xi, c in zip(x, ps)]
    w = [m / sum(mass) for m in mass]
    out = {
        "Li, no pressure correction": DIPPR9I(x, vms, ks),
        "DIPPR 9H plus Missenard": DIPPR9H(w, kp),
        "DWSIM's: Li on critical volumes, no pressure correction": DIPPR9I(x, [c.vc for c in ps], ks),
    }
    if len(ps) == 2:
        out["Filippov plus Missenard"] = Filippov(w, kp)
    return out


def surface_tension(ps, pt):
    t, x = pt["T"], pt["x"]
    if any(t >= c.tc for c in ps):
        return {"Diguilio-Teja": None, "DWSIM's: the mole-fraction mean": None}
    sig = [c.get("surface_tension", t) for c in ps]
    sig_tb = [c.get("surface_tension", c.tb) for c in ps]
    return {
        "Diguilio-Teja": safe(lambda: Diguilio_Teja(t, x, sig_tb, [c.tb for c in ps], [c.tc for c in ps])),
        "DWSIM's: the mole-fraction mean": sum(xi * s for xi, s in zip(x, sig)),
    }


ALTERNATIVES = {("viscosity", "vapour"): gas_viscosity, ("thermal_conductivity", "vapour"): gas_conductivity,
                ("viscosity", "liquid"): liquid_viscosity, ("thermal_conductivity", "liquid"): liquid_conductivity,
                ("surface_tension", "liquid"): surface_tension}


def figures(x, n=4):
    """x to n significant figures, so the file's numbers are short and the
    same from one run to the next."""
    return float(f"{x:.{n}g}")


def stats(method, pairs):
    devs = [calc / meas - 1 for calc, meas in pairs if calc is not None]
    out = {"method": method, "points": len(devs), "refused": len(pairs) - len(devs)}
    if devs:
        out["aad"] = figures(sum(abs(d) for d in devs) / len(devs))
        out["bias"] = figures(sum(devs) / len(devs))
        out["max"] = figures(max(devs, key=abs))
    return out


def main():
    sets = []
    for key, doi, k, prop, phase in SETS:
        root = fetch(doi)
        comps, cas, points, given = data_set(doi, k)
        ps = [Pure(c) for c in cas]
        res = {}
        for pt in points:
            for m, v in ALTERNATIVES[(prop, phase)](ps, pt).items():
                res.setdefault(m, []).append((v, pt["value"]))
        sets.append({
            "key": key,
            "doi": doi,
            "set": k,
            "citation": citation(root),
            "property": prop,
            "phase": phase,
            "compounds": [{"name": n, "cas": c} for n, c in zip(comps, cas)],
            "pressure_given": given,
            "points": points,
            "alternatives": [stats(m, pairs) for m, pairs in res.items()],
        })
    out = {
        "source": "NIST TRC ThermoML Archive, doi:10.18434/mds2-2422, each file by its paper's DOI from "
                  "https://trc.nist.gov/ThermoML/; points as printed, P 101.325 kPa where a set gives none, "
                  "a mass fraction to mole fractions by the component files' molar masses; u the expanded "
                  "uncertainty as the file gives it",
        "generator": "tools/gen_m5b_published.py (thermo 0.6.1, chemicals 1.5.2)",
        "accessed": datetime.date.today().isoformat(),
        "alternatives": "each on ChemSep's pure values from data/components, a gas's density by thermo's "
                        "PRMIX with kij 0 (its lower-Gibbs root); aad, bias and max as fractions of the "
                        "measured values; refused, the points it gives no value at",
        "sets": sets,
    }
    path = ROOT / "testdata" / "m5b_published.json"
    path.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path}: {len(sets)} sets, {sum(len(s['points']) for s in sets)} points")


if __name__ == "__main__":
    main()
