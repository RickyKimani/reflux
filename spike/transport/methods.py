# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2"]
# ///
"""M5b step 0: candidate mixing rules and pressure corrections against
measured mixtures from NIST's ThermoML archive, each compound's own value
from reflux's component file (ChemSep's correlations), as M5b would take it.

For each data set, every candidate's average and largest deviation, in %:

- gas viscosity (CO2 + CH4 to 30 MPa, CH4 + C3H8 and CH4 + C4H10 to 30
  MPa, H2 + CO2 at 300 kPa to 1100 K): Wilke, Herning-Zipperer and DWSIM's
  mole-fraction mean at low pressure; on Wilke, Lucas's high-pressure
  factor Y Fp Fq, as published or against its own value at zero
  pressure, or Jossi-Stiel-Thodos' residual on Peng-Robinson's
  density, at every density (as DWSIM) or from rho_r 0.1 (as published);
  DWSIM's whole method; chemicals' Lohrenz-Bray-Clark as it is;
- gas conductivity (CO2 + CH4 to 12 MPa): Wassiljewa with Mason and
  Saxena's A_ij (Wilke's Phi of the viscosities), Lindsay-Bromley and
  the mole-fraction mean (DWSIM), each with and without Stiel and Thodos'
  dense residual on Peng-Robinson's density and Yorizane's pseudo-critical
  rules;
- liquid viscosity (heptane + octane to 10 MPa and 513 K, methanol +
  water, diethylene glycol + water, diethylethanolamine + water, at 1
  atm): ln mu by mole fractions (Grunberg-Nissan with G = 0), with and
  without Lucas's compressed-liquid factor, by mass fractions, and DWSIM's
  mole mean;
- liquid conductivity (2-butoxyethanol + water to 150 MPa): Li's rule
  (DIPPR 9I) on the correlations' molar volumes, DIPPR 9H, Filippov, and
  DWSIM's Li on critical volumes, each with and without Missenard's
  pressure factor;
- surface tension (CH4 + C2H6, N2 + CH4, N2 + O2, diethylethanolamine +
  water): Winterfeld-Scriven-Davis on the correlations' molar volumes,
  Diguilio-Teja, and DWSIM's mole mean.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/transport/methods.py <cache dir>

The ThermoML files are fetched into <cache dir> once. Not part of the
simulator: it measures the candidates for the plan's decisions.
"""
import math
import pathlib
import sys
import tomllib
import urllib.request
import xml.etree.ElementTree as ET

from chemicals.interface import Diguilio_Teja, Winterfeld_Scriven_Davis
from chemicals.thermal_conductivity import DIPPR9H, DIPPR9I, Filippov, Lindsay_Bromley, Missenard, Stiel_Thodos_dense
from chemicals.viscosity import Herning_Zipperer, Lorentz_Bray_Clarke, Lucas, Wilke
from thermo.eos_mix import PRMIX

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
CACHE = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
R = 8.31446261815324
NAMES = {"methane": "74-82-8", "carbon dioxide": "124-38-9", "propane": "74-98-6", "butane": "106-97-8",
         "hydrogen": "1333-74-0", "methanol": "67-56-1", "water": "7732-18-5", "diethylene glycol": "111-46-6",
         "nitrogen": "7727-37-9", "oxygen": "7782-44-7", "ethane": "74-84-0", "2-butoxyethan-1-ol": "111-76-2",
         "(diethylamino)ethanol": "100-37-8", "heptane": "142-82-5", "octane": "111-65-9"}


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
        self.m, self.tb = k["mw"] * 1e3, k.get("tb")  # g/mol
        self.dipole = k.get("dipole", 0.0) / 3.33564e-30  # debye

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
    f = CACHE / (doi.replace("/", "_") + ".xml")
    if not f.exists():
        req = urllib.request.Request(f"https://trc.nist.gov/ThermoML/{doi}.xml", headers={"User-Agent": "curl/8.0"})
        with urllib.request.urlopen(req, timeout=120) as r:
            f.write_bytes(r.read())
    return ET.parse(f).getroot()


def rows(doi, k):
    """Set k's points as dicts: T K, P Pa, x (mole fractions in the set's
    component order), value SI; and its components' CAS numbers."""
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
    out = []
    for nv in child(d, "NumValues"):
        p = dict(fixed)
        for x in child(nv, "VariableValue"):
            p[cols[text(x, ["nVarNumber"])]] = float(text(x, ["nVarValue"]))
        val = float(text(child(nv, "PropertyValue")[0], ["nPropValue"]))
        t = p[("Temperature, K", None)]
        pr = 1e3 * p.get(("Pressure, kPa", None), 101.325)
        x = composition(p, comps, cas)
        out.append({"T": t, "P": pr, "x": x, "value": val})
    return cas, out


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
            mw = [Pure(s).m for s in cas]
            w = [0.0] * n
            w[i], w[1 - i] = v, 1 - v
            moles = [w[j] / mw[j] for j in range(n)]
            return [m / sum(moles) for m in moles]
    raise ValueError(p)


# --- methods ---

def pr_volume(ps, t, p, y):
    eos = PRMIX(Tcs=[c.tc for c in ps], Pcs=[c.pc for c in ps], omegas=[c.w for c in ps], zs=y, T=t, P=p,
                kijs=[[0.0] * len(ps) for _ in ps])
    if hasattr(eos, "V_l") and hasattr(eos, "V_g"):
        return eos.V_l if eos.G_dep_l < eos.G_dep_g else eos.V_g
    return eos.V_l if hasattr(eos, "V_l") else eos.V_g


def lucas_factor(ps, t, p, y):
    """Lucas's high-pressure ratio eta/eta0 = Y Fp Fq (Poling et al., 5th
    ed., sec. 9-6), with Lucas's mixture rules."""
    tcm = sum(yi * c.tc for yi, c in zip(y, ps))
    pcm = R * tcm * sum(yi * c.zc for yi, c in zip(y, ps)) / sum(yi * c.vc for yi, c in zip(y, ps))
    tr, pr = t / tcm, p / pcm
    fp0 = 0.0
    fq0 = 0.0
    for yi, c in zip(y, ps):
        mur = 52.46 * c.dipole ** 2 * (c.pc / 1e5) / c.tc ** 2
        tri = t / c.tc
        if mur < 0.022:
            f = 1.0
        elif mur < 0.075:
            f = 1 + 30.55 * (0.292 - c.zc) ** 1.72
        else:
            f = 1 + 30.55 * (0.292 - c.zc) ** 1.72 * abs(0.96 + 0.1 * (tri - 0.7))
        fp0 += yi * f
        q = {"Hydrogen": 0.76}.get(c.name)
        if q:
            s = 1.0 if tri > 12 else -1.0
            fq0 += yi * 1.22 * q ** 0.15 * (1 + 0.00385 * ((tri - 12) ** 2) ** (1 / c.m) * s)
        else:
            fq0 += yi
    z1 = (0.807 * tr ** 0.618 - 0.357 * math.exp(-0.449 * tr) + 0.340 * math.exp(-4.058 * tr) + 0.018) * fp0 * fq0
    if tr <= 1.0:
        if pr > 1.0:  # the form is for a vapour below its saturation pressure
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


def lucas_normalised(ps, t, p, y):
    """Lucas's factor against its own value at zero pressure, so that it's 1
    there below Tr 1 too, where its Z2 tends to 0.6 Tr, not Z1."""
    f = lucas_factor(ps, t, p, y)
    f0 = lucas_factor(ps, t, 1e-3, y)
    return None if f is None else f / f0


def jst(eta0, ps, t, v, y, cutoff):
    """Jossi, Stiel and Thodos' residual on eta0, Pa s, with the mixture's
    mole-fraction-mean Tc, Pc, Vc and M (as DWSIM and Lohrenz-Bray-Clark)."""
    tc = sum(yi * c.tc for yi, c in zip(y, ps))
    pc = sum(yi * c.pc for yi, c in zip(y, ps)) / 101325
    vc = sum(yi * c.vc for yi, c in zip(y, ps))
    m = sum(yi * c.m for yi, c in zip(y, ps))
    rr = vc / v
    if rr < cutoff:
        return eta0
    xi = tc ** (1 / 6) * m ** -0.5 * pc ** (-2 / 3)
    poly = 1.0230 + 0.23364 * rr + 0.58533 * rr ** 2 - 0.40758 * rr ** 3 + 0.093324 * rr ** 4
    return eta0 + (poly ** 4 - 1) / xi * 1e-7


def yorizane(ps, y):
    n = len(ps)
    vcij = [[(ps[i].vc ** (1 / 3) + ps[j].vc ** (1 / 3)) ** 3 / 8 for j in range(n)] for i in range(n)]
    vcm = sum(y[i] * y[j] * vcij[i][j] for i in range(n) for j in range(n))
    tcm = sum(y[i] * y[j] * vcij[i][j] * math.sqrt(ps[i].tc * ps[j].tc) for i in range(n) for j in range(n)) / vcm
    wm = sum(yi * c.w for yi, c in zip(y, ps))
    zcm = 0.291 - 0.08 * wm
    return tcm, zcm * R * tcm / vcm, vcm, zcm, sum(yi * c.m for yi, c in zip(y, ps))


def mason_saxena(ys, ks, mus, mws):
    n = len(ys)
    out = 0.0
    for i in range(n):
        den = sum(ys[j] * (1 + math.sqrt(mus[i] / mus[j]) * (mws[j] / mws[i]) ** 0.25) ** 2
                  / math.sqrt(8 * (1 + mws[i] / mws[j])) for j in range(n))
        out += ys[i] * ks[i] / den
    return out


def stats(pairs):
    devs = [100 * (calc / meas - 1) for calc, meas in pairs if calc is not None and math.isfinite(calc)]
    if not devs:
        return "no value"
    return f"AAD {sum(abs(d) for d in devs) / len(devs):5.2f}%  bias {sum(devs) / len(devs):+6.2f}%  " \
           f"max {max(devs, key=abs):+7.2f}%  ({len(devs)} pts)"


def report(title, results):
    print(f"\n== {title}")
    for name, pairs in results.items():
        print(f"   {name:42} {stats(pairs)}")


def gas_viscosity(doi, k, title):
    cas, pts = rows(doi, k)
    ps = [Pure(c) for c in cas]
    res = {m: [] for m in ["Wilke (no P)", "Herning-Zipperer (no P)", "mole mean (no P)", "Wilke x Lucas",
                           "Wilke x Lucas, against zero pressure",
                           "Wilke + JST, every density", "Wilke + JST from rho_r 0.1", "DWSIM: mole mean + JST",
                           "chemicals' Lohrenz-Bray-Clark"]}
    for pt in pts:
        t, p, y = pt["T"], pt["P"], pt["x"]
        mus = [c.get("vapor_viscosity", t) for c in ps]
        mws = [c.m for c in ps]
        w = Wilke(y, mus, mws) if len(ps) > 1 else mus[0]
        v = pr_volume(ps, t, p, y)
        mean = sum(yi * mi for yi, mi in zip(y, mus))
        res["Wilke (no P)"].append((w, pt["value"]))
        res["Herning-Zipperer (no P)"].append((Herning_Zipperer(y, mus, mws) if len(ps) > 1 else w, pt["value"]))
        res["mole mean (no P)"].append((mean, pt["value"]))
        f = lucas_factor(ps, t, p, y)
        res["Wilke x Lucas"].append((None if f is None else w * f, pt["value"]))
        fn = lucas_normalised(ps, t, p, y)
        res["Wilke x Lucas, against zero pressure"].append((None if fn is None else w * fn, pt["value"]))
        res["Wilke + JST, every density"].append((jst(w, ps, t, v, y, 0.0), pt["value"]))
        res["Wilke + JST from rho_r 0.1"].append((jst(w, ps, t, v, y, 0.1), pt["value"]))
        res["DWSIM: mole mean + JST"].append((jst(mean, ps, t, v, y, 0.0), pt["value"]))
        lbc = Lorentz_Bray_Clarke(t, p, v, y, mws, [c.tc for c in ps], [c.pc for c in ps], [c.vc for c in ps])
        res["chemicals' Lohrenz-Bray-Clark"].append((lbc, pt["value"]))
    t0, t1 = min(q["T"] for q in pts), max(q["T"] for q in pts)
    p0, p1 = min(q["P"] for q in pts), max(q["P"] for q in pts)
    report(f"{title}: {t0:g}-{t1:g} K, {p0 / 1e6:g}-{p1 / 1e6:g} MPa ({doi})", res)


def gas_conductivity(doi, k, title):
    cas, pts = rows(doi, k)
    ps = [Pure(c) for c in cas]
    names = ["Mason-Saxena", "Mason-Saxena + Stiel-Thodos", "Lindsay-Bromley", "Lindsay-Bromley + Stiel-Thodos",
             "DWSIM: mole mean"]
    res = {m: [] for m in names}
    for pt in pts:
        t, p, y = pt["T"], pt["P"], pt["x"]
        ks = [c.get("vapor_thermal_conductivity", t) for c in ps]
        mus = [c.get("vapor_viscosity", t) for c in ps]
        mws = [c.m for c in ps]
        ms = mason_saxena(y, ks, mus, mws)
        lb = Lindsay_Bromley(t, y, ks, mus, [c.tb for c in ps], mws) if len(ps) > 1 else ks[0]
        v = pr_volume(ps, t, p, y)
        tcm, pcm, vcm, zcm, mm = yorizane(ps, y)
        res["Mason-Saxena"].append((ms, pt["value"]))
        res["Mason-Saxena + Stiel-Thodos"].append((Stiel_Thodos_dense(t, mm, tcm, pcm, vcm, zcm, v, ms), pt["value"]))
        res["Lindsay-Bromley"].append((lb, pt["value"]))
        res["Lindsay-Bromley + Stiel-Thodos"].append((Stiel_Thodos_dense(t, mm, tcm, pcm, vcm, zcm, v, lb),
                                                      pt["value"]))
        res["DWSIM: mole mean"].append((sum(yi * ki for yi, ki in zip(y, ks)), pt["value"]))
    p1 = max(q["P"] for q in pts)
    report(f"{title}: to {p1 / 1e6:g} MPa ({doi})", res)


def liquid_viscosity(doi, k, title):
    cas, pts = rows(doi, k)
    ps = [Pure(c) for c in cas]
    res = {m: [] for m in ["ln mu by mole fractions", "ln mu by mole fractions, each x Lucas",
                           "ln mu by mass fractions", "DWSIM: mole mean"]}
    for pt in pts:
        t, x = pt["T"], pt["x"]
        mus = [c.get("liquid_viscosity", t) for c in ps]
        lucas = [Lucas(t, pt["P"], c.tc, c.pc, c.w, c.get("vapor_pressure", t), m) for c, m in zip(ps, mus)]
        res["ln mu by mole fractions, each x Lucas"].append(
            (math.exp(sum(xi * math.log(m) for xi, m in zip(x, lucas))), pt["value"]))
        mass = [xi * c.m for xi, c in zip(x, ps)]
        w = [m / sum(mass) for m in mass]
        res["ln mu by mole fractions"].append((math.exp(sum(xi * math.log(m) for xi, m in zip(x, mus))), pt["value"]))
        res["ln mu by mass fractions"].append((math.exp(sum(wi * math.log(m) for wi, m in zip(w, mus))), pt["value"]))
        res["DWSIM: mole mean"].append((sum(xi * m for xi, m in zip(x, mus)), pt["value"]))
    report(f"{title} ({doi})", res)


def liquid_conductivity(doi, k, title):
    cas, pts = rows(doi, k)
    ps = [Pure(c) for c in cas]
    names = ["Li (DIPPR 9I)", "Li + Missenard", "DIPPR 9H", "DIPPR 9H + Missenard", "Filippov + Missenard",
             "DWSIM: Li on critical volumes"]
    res = {m: [] for m in names}
    for pt in pts:
        t, p, x = pt["T"], pt["P"], pt["x"]
        ks = [c.get("liquid_thermal_conductivity", t) for c in ps]
        kp = [Missenard(t, p, c.tc, c.pc, kk) for c, kk in zip(ps, ks)]
        vms = [c.vl(t) for c in ps]
        mass = [xi * c.m for xi, c in zip(x, ps)]
        w = [m / sum(mass) for m in mass]
        res["Li (DIPPR 9I)"].append((DIPPR9I(x, vms, ks), pt["value"]))
        res["Li + Missenard"].append((DIPPR9I(x, vms, kp), pt["value"]))
        res["DIPPR 9H"].append((DIPPR9H(w, ks), pt["value"]))
        res["DIPPR 9H + Missenard"].append((DIPPR9H(w, kp), pt["value"]))
        res["Filippov + Missenard"].append((Filippov(w, kp), pt["value"]))
        res["DWSIM: Li on critical volumes"].append((DIPPR9I(x, [c.vc for c in ps], ks), pt["value"]))
    p1 = max(q["P"] for q in pts)
    report(f"{title}: to {p1 / 1e6:g} MPa ({doi})", res)


def surface_tension(doi, k, title):
    try:
        cas, pts = rows(doi, k)
    except ValueError as e:  # a set whose composition isn't a variable or constraint
        print(f"\n== {title} ({doi}): no composition in the set: {e}")
        return
    ps = [Pure(c) for c in cas]
    res = {m: [] for m in ["Winterfeld-Scriven-Davis", "Diguilio-Teja", "DWSIM: mole mean"]}
    for pt in pts:
        t, x = pt["T"], pt["x"]
        if any(t >= c.tc for c in ps):
            continue
        sig = [c.get("surface_tension", t) for c in ps]
        res["Winterfeld-Scriven-Davis"].append((Winterfeld_Scriven_Davis(x, sig, [1 / c.vl(t) for c in ps]),
                                                pt["value"]))
        sig_tb = [c.get("surface_tension", c.tb) for c in ps]
        try:
            dt = Diguilio_Teja(t, x, sig_tb, [c.tb for c in ps], [c.tc for c in ps])
        except Exception:  # noqa: BLE001
            dt = None
        res["Diguilio-Teja"].append((dt, pt["value"]))
        res["DWSIM: mole mean"].append((sum(xi * s for xi, s in zip(x, sig)), pt["value"]))
    report(f"{title} ({doi})", res)


def main():
    gas_viscosity("10.1021/je0342419", 2, "gas viscosity, H2 + CO2 at 300 kPa")
    gas_viscosity("10.1016/j.jct.2015.03.007", 0, "gas viscosity, CO2 + CH4 (x 0.5174)")
    gas_viscosity("10.1021/je500857w", 0, "gas viscosity, CH4 + C3H8 (x 0.949)")
    gas_viscosity("10.1021/acs.jced.5b00635", 0, "gas viscosity, CH4 + C4H10 (x 0.9484)")
    gas_conductivity("10.1007/s10765-005-5566-6", 0, "gas conductivity, pure CO2")
    gas_conductivity("10.1007/s10765-005-5566-6", 1, "gas conductivity, CO2 + CH4")
    liquid_viscosity("10.1016/j.jct.2006.01.012", 4, "liquid viscosity, heptane + octane, to 10 MPa")
    liquid_viscosity("10.1016/j.tca.2015.08.005", 2, "liquid viscosity, heptane + octane, 3 to 5 MPa")
    liquid_viscosity("10.1016/j.jct.2019.05.013", 2, "liquid viscosity, methanol + water")
    liquid_viscosity("10.1021/je1009976", 2, "liquid viscosity, diethylene glycol + water")
    liquid_viscosity("10.1021/je700350b", 3, "liquid viscosity, diethylethanolamine + water")
    liquid_conductivity("10.1021/je1007197", 1, "liquid conductivity, 2-butoxyethanol + water")
    surface_tension("10.1016/j.fluid.2013.07.008", 0, "surface tension, CH4 + C2H6")
    surface_tension("10.1016/j.fluid.2010.11.013", 0, "surface tension, N2 + CH4")
    surface_tension("10.1016/j.fluid.2008.06.016", 2, "surface tension, N2 + O2")
    surface_tension("10.1021/je700350b", 1, "surface tension, diethylethanolamine + water")


if __name__ == "__main__":
    main()
