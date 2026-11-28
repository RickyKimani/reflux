# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["chemicals==1.5.2", "thermo==0.6.1", "mpmath==1.3.0"]
# ///
"""Write M5b's transport references (bundle A, steps 1 to 4; D18):
testdata/m5b_transport.json.

Every number is from the component files (data/components/), evaluated
independently of reflux:

- pure: every transport correlation of every compound (the liquid's and
  vapour's viscosity and conductivity, the surface tension, the liquid
  density, by their index in keys) at three temperatures over its range, a
  saturated liquid's below Tc, by chemicals' DIPPR forms (EQ100 to EQ116,
  Chemsep_16) and ChemSep's equation 10 written out, in SI: criterion 1.
  Each row is [CAS, key, T, value], with a fifth, the size of its terms,
  where that isn't the value's own: what rounding is relative to, in a
  correlation that cancels (equation 16's A + e^(...) near a zero, a
  polynomial) or an exponential form whose exponent's terms are large
  (DIPPR 101's reach 100 for a liquid's viscosity).
- lucas_dilute: Lucas's own low-pressure gas viscosity by chemicals'
  Lucas_gas, which holds the polar and quantum factors Fp° and Fq°.
- lucas_factor: Lucas's high-pressure factor Y·Fp·Fq with Lucas's mixture
  rules, Y taken against its value at zero pressure (M5b's D5), by mpmath at
  50 digits from Poling, Prausnitz and O'Connell's equations (5th ed.,
  §9-6), written again here.
- mason_saxena: a gas's conductivity at low pressure, Wassiljewa's form
  with Mason and Saxena's A_ij = Wilke's Phi_ij of the viscosities, by
  mpmath on chemicals' pure values.
- stiel_thodos: chemicals' Stiel_Thodos_dense residual on Yorizane's
  pseudo-critical constants (mpmath), at Peng-Robinson's molar volume by
  thermo (ChemSep's Tc, Pc, omega, no kij).
- lucas_liquid, missenard: chemicals' Lucas and Missenard factors.
- li, winterfeld, grunberg_nissan: chemicals' DIPPR9I and
  Winterfeld_Scriven_Davis, and mpmath's Grunberg-Nissan, on given numbers.
- phases: a phase's transport put together from the pieces above, as M5b's
  D4 to D9 have it, for transport.Mixture.phase end to end.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m5b_transport.py

Its output is the same, byte for byte, for the same versions.
"""
import json
import math
import pathlib
import tomllib

import chemicals
import mpmath as mp
import thermo
from chemicals.dippr import EQ100, EQ101, EQ102, EQ105, EQ106, EQ115, EQ116
from chemicals.interface import Winterfeld_Scriven_Davis
from chemicals.thermal_conductivity import DIPPR9I, Chemsep_16, Missenard, Stiel_Thodos_dense
from chemicals.viscosity import Lucas, Lucas_gas, Wilke
from thermo.eos_mix import PRMIX

mp.mp.dps = 50
ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "components"
OUT = ROOT / "testdata" / "m5b_transport.json"
R = 8.31446261815324
DEBYE = 3.33564e-30
QUANTUM = {"7440-59-7": 1.38, "1333-74-0": 0.76, "7782-39-0": 0.52}
KEYS = ["liquid_viscosity", "vapor_viscosity", "liquid_thermal_conductivity", "vapor_thermal_conductivity",
        "surface_tension", "liquid_density"]
SATURATED = {"liquid_viscosity", "liquid_thermal_conductivity", "surface_tension", "liquid_density",
             "vapor_pressure"}
SCALE = {"Pa.s": 1.0, "W/m/K": 1.0, "N/m": 1.0, "kmol/m3": 1e3, "Pa": 1.0}

CH4, C2H6, C3H8, N2, CO2 = "74-82-8", "74-84-0", "74-98-6", "7727-37-9", "124-38-9"
H2, H2O, CO, MEOH, NH3 = "1333-74-0", "7732-18-5", "630-08-0", "67-56-1", "7664-41-7"
BENZENE, TOLUENE, OCTANE, O2 = "71-43-2", "108-88-3", "111-65-9", "7782-44-7"


def record(cas):
    return tomllib.loads((DATA / f"{cas}.toml").read_text(encoding="utf-8"))


class Pure:
    def __init__(self, cas):
        r = record(cas)
        k = r["constants"]
        self.r, self.cas, self.name = r, cas, r["name"]
        self.tc, self.pc, self.vc, self.zc, self.w = k["tc"], k["pc"], k["vc"], k["zc"], k["omega"]
        self.mw = k["mw"]
        self.dipole = k["dipole"] / DEBYE if "dipole" in k else None
        self.q = QUANTUM.get(cas, 0.0)

    def get(self, key, t):
        return correlation(self.r, key, t)[0]


def correlation(r, key, t):
    """(value in SI, the size of its terms in SI), by chemicals' forms. The
    size is what rounding is relative to: a polynomial's terms, equation
    16's A beside its exponential, or for an exponential form the value
    times its exponent's terms, whose rounding the exponential carries."""
    d = r["correlations"][key]
    c = (d["coeffs"] + [0.0] * 5)[:5]
    a, b, cc, dd, e = c
    tc = d.get("tc", r["constants"].get("tc"))
    eq = d["eq"]
    s = SCALE[d["units"]]
    if eq in (1, 2, 3, 4, 5, 100):
        y = EQ100(t, a, b, cc, dd, e)
        terms = sum(abs(x * t ** i) for i, x in enumerate(c))
    elif eq == 10:
        y = math.exp(a - b / (t + cc))
        terms = abs(y) * max(1.0, abs(a) + abs(b / (t + cc)))
    elif eq == 16:
        y = Chemsep_16(t, a, b, cc, dd, e)
        terms = abs(a) + abs(y - a) * max(1.0, abs(b / t) + abs(cc) + abs(dd * t) + abs(e * t * t))
    elif eq == 101:
        y = EQ101(t, a, b, cc, dd, e)
        terms = abs(y) * max(1.0, abs(a) + abs(b / t) + abs(cc * math.log(t)) + abs(dd * t ** e))
    elif eq == 102:
        y = EQ102(t, a, b, cc, dd)
        terms = abs(y) * max(1.0, abs(b * math.log(t)) + abs(cc / t) + abs(dd / (t * t)))
    elif eq == 105:
        y = EQ105(t, a, b, cc, dd)
        terms = abs(y) * max(1.0, abs((1 + (1 - t / cc) ** dd) * math.log(b)))
    elif eq == 106:
        y = EQ106(t, tc, a, b, cc, dd, e)
        tr = t / tc
        h = b + tr * (cc + tr * (dd + tr * e))
        terms = abs(y) * max(1.0, abs(h * math.log(1 - tr)))
    elif eq == 115:
        y = EQ115(t, a, b, cc, dd, e)
        terms = abs(y) * max(1.0, abs(a) + abs(b / t) + abs(cc * math.log(t)) + abs(dd * t * t) + abs(e / (t * t)))
    elif eq == 116:
        y = EQ116(t, tc, a, b, cc, dd, e)
        terms = abs(a) + abs(b) + abs(cc) + abs(dd) + abs(e)
    else:
        raise ValueError(eq)
    return s * y, s * terms


def finite(x):
    return x if isinstance(x, float) and math.isfinite(x) else None


def pure_rows():
    rows = []
    for f in sorted(DATA.glob("[0-9]*.toml")):
        r = tomllib.loads(f.read_text(encoding="utf-8"))
        tc = r["constants"].get("tc")
        for key in KEYS:
            d = r.get("correlations", {}).get(key)
            if d is None or not d["coeffs"] or d["eq"] == 118:
                continue
            lo, hi = d["t_range"]
            if key in SATURATED and tc is not None:
                hi = min(hi, 0.99 * tc)
            if not lo < hi:
                continue
            for t in (lo, (lo + hi) / 2, hi):
                t = float(t)
                try:
                    y, sc = correlation(r, key, t)
                except (OverflowError, ValueError, ZeroDivisionError):
                    y, sc = float("nan"), float("nan")
                y, sc = finite(y), finite(sc)
                row = [r["ids"]["cas"], KEYS.index(key), t, y]
                if y is None or sc != abs(y):
                    row.append(sc)
                rows.append(row)
    return rows


# --- Lucas, by mpmath from the published equations, written again ---

def m_fp0(c, tr):
    mur = mp.mpf(52.46) * mp.mpf(c.dipole) ** 2 * (mp.mpf(c.pc) * mp.mpf("1e-5")) / mp.mpf(c.tc) ** 2
    if mur < mp.mpf("0.022"):
        return mp.mpf(1)
    f = mp.mpf("30.55") * (mp.mpf("0.292") - mp.mpf(c.zc)) ** mp.mpf("1.72")
    if mur < mp.mpf("0.075"):
        return 1 + f
    return 1 + f * abs(mp.mpf("0.96") + mp.mpf("0.1") * (tr - mp.mpf("0.7")))


def m_fq0(c, tr):
    if c.q == 0:
        return mp.mpf(1)
    x = tr - 12
    sign = 1 if x > 0 else -1
    return mp.mpf("1.22") * mp.mpf(c.q) ** mp.mpf("0.15") * (
        1 + mp.mpf("0.00385") * (x * x) ** (1 / (mp.mpf(c.mw) * 1000)) * sign)


def m_z1(tr, fp0, fq0):
    return (mp.mpf("0.807") * tr ** mp.mpf("0.618") - mp.mpf("0.357") * mp.exp(mp.mpf("-0.449") * tr)
            + mp.mpf("0.340") * mp.exp(mp.mpf("-4.058") * tr) + mp.mpf("0.018")) * fp0 * fq0


def m_z2_low(tr, pr):
    alpha = mp.mpf("3.262") + mp.mpf("14.98") * pr ** mp.mpf("5.508")
    beta = mp.mpf("1.390") + mp.mpf("5.746") * pr
    return mp.mpf("0.600") + mp.mpf("0.760") * pr ** alpha + (mp.mpf("6.990") * pr ** beta - mp.mpf("0.6")) * (1 - tr)


def lucas_factor(cs, y, t, p):
    """Y Fp Fq (M5b's D5), or None where the form doesn't hold."""
    y = [mp.mpf(v) for v in y]
    t, p = mp.mpf(t), mp.mpf(p)
    tcm = sum(yi * c.tc for yi, c in zip(y, cs))
    zs = sum(yi * c.zc for yi, c in zip(y, cs))
    vs = sum(yi * c.vc for yi, c in zip(y, cs))
    fp0 = sum(yi * m_fp0(c, t / c.tc) for yi, c in zip(y, cs))
    fq0 = sum(yi * m_fq0(c, t / c.tc) for yi, c in zip(y, cs))
    heavy = max(range(len(cs)), key=lambda i: cs[i].mw)
    light = min(range(len(cs)), key=lambda i: cs[i].mw)
    ratio = mp.mpf(cs[heavy].mw) / mp.mpf(cs[light].mw)
    # A's bounds compared as doubles, as reflux does: a gas of exactly 0.05
    # CO2, the reforming gas, is on the boundary, and mpmath's decimal 0.05
    # is below the double nearest it.
    if float(ratio) > 9 and 0.05 < float(y[heavy]) < 0.7:
        fq0 *= 1 - mp.mpf("0.01") * ratio ** mp.mpf("0.87")
    pcm = mp.mpf(R) * tcm * zs / vs
    tr, pr = t / tcm, p / pcm
    if tr <= 1:
        if not pr < 1:
            return None
        yy = m_z2_low(tr, pr) / m_z2_low(tr, mp.mpf(0))
    else:
        if not (tr < 40 and pr <= 100):
            return None
        z1 = m_z1(tr, fp0, fq0)
        a = mp.mpf("1.245e-3") / tr * mp.exp(mp.mpf("5.1726") * tr ** mp.mpf("-0.3286"))
        b = a * (mp.mpf("1.6553") * tr - mp.mpf("1.2723"))
        c = mp.mpf("0.4489") / tr * mp.exp(mp.mpf("3.0578") * tr ** mp.mpf("-37.7332"))
        d = mp.mpf("1.7368") / tr * mp.exp(mp.mpf("2.2310") * tr ** mp.mpf("-7.6351"))
        f = mp.mpf("0.9425") * mp.exp(mp.mpf("-0.1853") * tr ** mp.mpf("0.4489"))
        z2 = z1 * (1 + a * pr ** mp.mpf("1.3088") / (b * pr ** f + 1 / (1 + c * pr ** d)))
        yy = z2 / z1
    fp = (1 + (fp0 - 1) * yy ** -3) / fp0
    fq = (1 + (fq0 - 1) * (1 / yy - mp.mpf("0.007") * mp.log(yy) ** 4)) / fq0
    return float(yy * fp * fq)


def mason_saxena(y, k, mu, mw):
    y, k, mu, mw = ([mp.mpf(v) for v in a] for a in (y, k, mu, mw))
    out = mp.mpf(0)
    for i in range(len(y)):
        if y[i] == 0:
            continue
        den = sum(y[j] * (1 + mp.sqrt(mu[i] / mu[j]) * (mw[j] / mw[i]) ** mp.mpf("0.25")) ** 2
                  / mp.sqrt(8 * (1 + mw[i] / mw[j])) for j in range(len(y)))
        out += y[i] * k[i] / den
    return float(out)


def yorizane(cs, y):
    y = [mp.mpf(v) for v in y]
    n = len(cs)
    vcij = [[(mp.cbrt(cs[i].vc) + mp.cbrt(cs[j].vc)) ** 3 / 8 for j in range(n)] for i in range(n)]
    vcm = sum(y[i] * y[j] * vcij[i][j] for i in range(n) for j in range(n))
    tcm = sum(y[i] * y[j] * vcij[i][j] * mp.sqrt(mp.mpf(cs[i].tc) * cs[j].tc) for i in range(n) for j in range(n)) / vcm
    zcm = mp.mpf("0.291") - mp.mpf("0.08") * sum(yi * c.w for yi, c in zip(y, cs))
    return {"tc": float(tcm), "pc": float(zcm * R * tcm / vcm), "vc": float(vcm), "zc": float(zcm),
            "mw": float(sum(yi * c.mw for yi, c in zip(y, cs)))}


def stiel_thodos(cs, y, v):
    ps = yorizane(cs, y)
    if ps["vc"] / v > 2.8:
        return ps, None
    # chemicals' MW in g/mol; with kg = 0 it returns the residual alone
    return ps, Stiel_Thodos_dense(1.0, ps["mw"] * 1e3, ps["tc"], ps["pc"], ps["vc"], ps["zc"], v, 0.0)


def grunberg_nissan(x, mu, g):
    n = len(x)
    x, mu = [mp.mpf(v) for v in x], [mp.mpf(v) for v in mu]
    s = sum(x[i] * mp.log(mu[i]) for i in range(n) if x[i] != 0)
    e = sum(x[i] * x[j] * mp.mpf(g[i][j]) for i in range(n) for j in range(n))
    return float(mp.exp(s + e / 2))


def pr_volume(cs, t, p, y, root):
    eos = PRMIX(Tcs=[c.tc for c in cs], Pcs=[c.pc for c in cs], omegas=[c.w for c in cs], zs=y, T=t, P=p,
                kijs=[[0.0] * len(cs) for _ in cs])
    if root == "vapour":
        return eos.V_g if hasattr(eos, "V_g") else eos.V_l
    return eos.V_l if hasattr(eos, "V_l") else eos.V_g


def norm(a):
    s = sum(a)
    return [v / s for v in a]


# --- the cases ---

GASES = [
    ("natural gas", [CH4, C2H6, C3H8, N2, CO2], [0.85, 0.07, 0.03, 0.03, 0.02]),
    ("the reforming gas mid-reactor", [CH4, H2, H2O, CO2, CO], [0.15, 0.5, 0.2, 0.05, 0.1]),
    ("air", [N2, O2], [0.79, 0.21]),
    ("methane and propane", [CH4, C3H8], [0.949, 0.051]),
    ("hydrogen and carbon dioxide", [H2, CO2], [0.6, 0.4]),
    ("methanol and water", [MEOH, H2O], [0.7, 0.3]),
    ("benzene and toluene", [BENZENE, TOLUENE], [0.56, 0.44]),
    ("carbon dioxide and methane", [CO2, CH4], [0.5174, 0.4826]),
    ("ammonia", [NH3], [1.0]),
    ("steam", [H2O], [1.0]),
    ("hydrogen", [H2], [1.0]),
]
GAS_STATES = {
    "natural gas": [(300.0, 1e5), (300.0, 1e6), (300.0, 1e7), (300.0, 3e7), (400.0, 1e7)],
    "the reforming gas mid-reactor": [(1000.0, 3e6), (1000.0, 1e5), (700.0, 3e6)],
    "air": [(300.0, 1e5), (300.0, 1e7), (500.0, 5e7)],
    "methane and propane": [(320.0, 1e7), (420.0, 3e7)],
    "hydrogen and carbon dioxide": [(500.0, 3e5), (500.0, 5e6), (1100.0, 3e5)],
    "methanol and water": [(350.0, 101325.0), (350.0, 1e3), (500.0, 1e6)],
    "benzene and toluene": [(370.0, 100.0), (370.0, 101325.0)],
    "carbon dioxide and methane": [(350.0, 1e7), (300.0, 2e7)],
    "ammonia": [(420.0, 3e7), (300.0, 5e5)],
    "steam": [(600.0, 5e6), (800.0, 3e7)],
    "hydrogen": [(300.0, 1e7), (50.0, 1e5)],
}

LIQUID_STATES = [(BENZENE, 300.0, 1e5), (BENZENE, 300.0, 1e7), (BENZENE, 300.0, 5e7), (H2O, 350.0, 1e5),
                 (H2O, 350.0, 1e7), (MEOH, 350.0, 1e5), (OCTANE, 450.0, 1e7), (C3H8, 360.0, 5e6),
                 (TOLUENE, 500.0, 2e7)]


def main():
    out = {
        "source": ("M5b's transport references, from the component files (ChemSep 8.50): chemicals "
                   + chemicals.__version__ + "' DIPPR forms, Lucas_gas, Lucas, Missenard, "
                   "Stiel_Thodos_dense, DIPPR9I, Winterfeld_Scriven_Davis and Wilke; mpmath "
                   + mp.__version__ + " at 50 digits for Lucas's high-pressure factor (Y against its "
                   "value at zero pressure, M5b's D5), Mason-Saxena, Grunberg-Nissan and Yorizane's "
                   "rules; thermo " + thermo.__version__ + "'s PRMIX for gas volumes. "
                   "tools/gen_m5b_transport.py"),
    }
    out["keys"] = KEYS
    out["pure"] = pure_rows()

    out["lucas_dilute"] = []
    for cas in [CH4, H2, H2O, MEOH, NH3, CO2, BENZENE, "7440-59-7"]:
        c = Pure(cas)
        for t in (50.0, 300.0, 1000.0):
            mu = Lucas_gas(t, c.tc, c.pc, c.zc, c.mw * 1e3, c.dipole, cas)
            out["lucas_dilute"].append([cas, t, mu])

    out["lucas_factor"], out["mason_saxena"], out["stiel_thodos"] = [], [], []
    for name, cas, y in GASES:
        cs = [Pure(c) for c in cas]
        y = norm(y)
        for t, p in GAS_STATES[name]:
            f = lucas_factor(cs, y, t, p)
            out["lucas_factor"].append({"name": f"{name} at {t:g} K and {p:g} Pa", "cas": cas, "y": y, "t": t,
                                        "p": p, "factor": f})
            mu = [c.get("vapor_viscosity", t) for c in cs]
            k = [c.get("vapor_thermal_conductivity", t) for c in cs]
            out["mason_saxena"].append({"name": f"{name} at {t:g} K", "cas": cas, "y": y, "t": t,
                                        "k": mason_saxena(y, k, mu, [c.mw for c in cs])})
            v = pr_volume(cs, t, p, y, "vapour")
            ps, res = stiel_thodos(cs, y, v)
            out["stiel_thodos"].append({"name": f"{name} at {t:g} K and {p:g} Pa", "cas": cas, "y": y, "v": v,
                                        "pseudo": ps, "residual": res})

    out["lucas_liquid"] = []
    for cas, t, p in LIQUID_STATES:
        c = Pure(cas)
        psat = c.get("vapor_pressure", t)
        out["lucas_liquid"].append([cas, t, p, psat, Lucas(t, p, c.tc, c.pc, c.w, psat, 1.0)])

    out["missenard"] = []
    for tr in (0.4, 0.5, 0.55, 0.65, 0.7, 0.75, 0.8, 0.9):
        for pr in (0.0, 0.5, 1.0, 3.0, 20.0, 75.0, 150.0, 200.0, 300.0):
            tc, pc = 562.05, 4.895e6
            out["missenard"].append([tr * tc, pr * pc, tc, pc, Missenard(tr * tc, pr * pc, tc, pc, 1.0)])

    out["li"], out["winterfeld"], out["grunberg_nissan"] = [], [], []
    mixes = [
        ([0.3, 0.7], [8.9e-5, 1.8e-5], [0.15, 0.6], [0.022, 0.072], [5.5e-4, 8.9e-4], [[0, 0], [0, 0]]),
        ([0.3, 0.7], [8.9e-5, 1.8e-5], [0.15, 0.6], [0.022, 0.072], [5.5e-4, 8.9e-4], [[0, 1.8], [1.8, 0]]),
        ([0.2, 0.5, 0.3], [9e-5, 1.1e-4, 4e-5], [0.13, 0.12, 0.2], [0.028, 0.026, 0.023], [6e-4, 5.5e-4, 3e-4],
         [[0, 0.1, -0.3], [0.1, 0, 0.05], [-0.3, 0.05, 0]]),
        ([0.25, 0.25, 0.25, 0.25], [1e-4, 2e-4, 3e-5, 6e-5], [0.1, 0.14, 0.5, 0.2], [0.02, 0.03, 0.07, 0.025],
         [3e-4, 1e-3, 9e-4, 4e-4], [[0, 0, 0, 0], [0, 0, 0.7, 0], [0, 0.7, 0, 0], [0, 0, 0, 0]]),
        ([1.0, 0.0], [9e-5, 1.8e-5], [0.15, 0.6], [0.022, 0.072], [5.5e-4, 8.9e-4], [[0, 1.8], [1.8, 0]]),
    ]
    for x, v, k, sig, mu, g in mixes:
        out["li"].append({"x": x, "v": v, "k": k, "value": DIPPR9I(x, v, k)})
        out["winterfeld"].append({"x": x, "v": v, "sigma": sig,
                                  "value": Winterfeld_Scriven_Davis(x, sig, [1 / vi for vi in v])})
        out["grunberg_nissan"].append({"x": x, "mu": mu, "g": g, "value": grunberg_nissan(x, mu, g)})

    out["phases"] = phases()
    OUT.write_bytes(dump(out).encode("utf-8"))
    print(f"{len(out['pure'])} pure points, {len(out['phases'])} phases")


def liquid_values(cs, x, t, p):
    mu = [c.get("liquid_viscosity", t) * Lucas(t, p, c.tc, c.pc, c.w, c.get("vapor_pressure", t), 1.0)
          for c in cs]
    k = [c.get("liquid_thermal_conductivity", t) * Missenard(t, p, c.tc, c.pc, 1.0) for c in cs]
    sig = [c.get("surface_tension", t) for c in cs]
    v = [1 / c.get("liquid_density", t) for c in cs]
    if len(cs) == 1:
        return mu[0], k[0], sig[0]
    n = len(cs)
    return (grunberg_nissan(x, mu, [[0] * n for _ in range(n)]), DIPPR9I(x, v, k),
            Winterfeld_Scriven_Davis(x, sig, [1 / vi for vi in v]))


def gas_values(cs, y, t, p, v):
    mu = [c.get("vapor_viscosity", t) for c in cs]
    k = [c.get("vapor_thermal_conductivity", t) for c in cs]
    mws = [c.mw for c in cs]
    mu0 = Wilke(y, mu, mws) if len(cs) > 1 else mu[0]
    f = lucas_factor(cs, y, t, p)
    _, res = stiel_thodos(cs, y, v)
    return (None if f is None else mu0 * f, mason_saxena(y, k, mu, mws) + res, None)


def phases():
    rows = []
    for name, cas, x, t, p in [
        ("benzene and toluene, a liquid at 300 K and 1 atm", [BENZENE, TOLUENE], [0.4, 0.6], 300.0, 101325.0),
        ("benzene and toluene, a liquid at 300 K and 100 bar", [BENZENE, TOLUENE], [0.4, 0.6], 300.0, 1e7),
        ("methanol and water, a liquid at 300 K and 1 atm", [MEOH, H2O], [0.3, 0.7], 300.0, 101325.0),
        ("water alone, a liquid at 350 K and 10 bar", [H2O], [1.0], 350.0, 1e6),
    ]:
        cs = [Pure(c) for c in cas]
        mu, k, sig = liquid_values(cs, x, t, p)
        rows.append({"name": name, "cas": cas, "label": "liquid", "x": x, "t": t, "p": p, "v": 1e-4,
                     "mu": mu, "k": k, "sigma": sig})
    # hydrogen above its Tc is left out of the liquid (D9): water's own values
    cs = [Pure(H2O)]
    mu, k, sig = liquid_values(cs, [1.0], 320.0, 3e6)
    rows.append({"name": "water with 1% hydrogen, a liquid at 320 K and 30 bar", "cas": [H2O, H2],
                 "label": "liquid", "x": [0.99, 0.01], "t": 320.0, "p": 3e6, "v": 1e-4, "mu": mu, "k": k,
                 "sigma": sig})
    for name, cas, y, t, p, label in [
        ("natural gas at 300 K and 100 bar", [CH4, C2H6, C3H8, N2, CO2], [0.85, 0.07, 0.03, 0.03, 0.02], 300.0,
         1e7, "vapour"),
        ("the reforming gas at 1000 K and 30 bar", [CH4, H2, H2O, CO2, CO], [0.15, 0.5, 0.2, 0.05, 0.1], 1000.0,
         3e6, "vapour"),
        ("methanol and water, a vapour at 350 K and 1 atm", [MEOH, H2O], [0.7, 0.3], 350.0, 101325.0, "vapour"),
        ("carbon dioxide labelled liquid at 320 K and 150 bar", [CO2], [1.0], 320.0, 1.5e7, "liquid"),
    ]:
        cs = [Pure(c) for c in cas]
        y = norm(y)
        v = pr_volume(cs, t, p, y, "vapour")
        mu, k, sig = gas_values(cs, y, t, p, v)
        rows.append({"name": name, "cas": cas, "label": label, "x": y, "t": t, "p": p, "v": v, "mu": mu, "k": k,
                     "sigma": sig})
    return rows


def dump(o):
    """One line for each row of the long lists, the rest indented."""
    items = []
    for k, v in o.items():
        if isinstance(v, list):
            body = ",\n".join("  " + json.dumps(r, separators=(",", ":"), ensure_ascii=False) for r in v)
            items.append(f" {json.dumps(k)}: [\n{body}\n ]")
        else:
            items.append(f" {json.dumps(k)}: {json.dumps(v, ensure_ascii=False)}")
    return "{\n" + ",\n".join(items) + "\n}\n"


if __name__ == "__main__":
    main()
