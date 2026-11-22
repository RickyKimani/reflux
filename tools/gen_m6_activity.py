# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2", "mpmath==1.3.0"]
# ///
"""Write M6's activity-model references (steps 1 to 3, criteria 1 and 2):

- testdata/m6_activity.json: Wilson, NRTL, UNIQUAC and the three UNIFACs at
  seeded states of two to five components, on ChemSep's own banks and
  tables (data/vendor/chemsep/) and each compound's component file. For
  every state, thermo 0.6.1's ln gamma, its first two T-derivatives,
  n*d(ln gamma_i)/dn_j, G^E, H^E and Cp^E; for NRTL, UNIQUAC and Wilson the
  same by mpmath at 40 digits from the closed forms, differentiated by
  mpmath.diff.
- testdata/m6_unifac_tables.json: thermo's own copies of the published
  UNIFAC tables (the original, LLE and Dortmund's 2016), which criterion 2
  compares ChemSep's with, parameter by parameter.

It reads ChemSep's tables as shipped, without the repairs and corrections
reflux's reader makes (docs/data-quality.md): none of the systems here uses
a group they touch, so the two read the same numbers, and a test that
didn't would fail on the difference.

Each bank's parameters are taken as A in cal/mol over DECHEMA's R, 1.98721
cal/(mol K): NRTL's tau = A/(1.98721 T), UNIQUAC's tau = exp(-A/(1.98721 T)),
Wilson's Lambda = (Vj/Vi) exp(-A/(1.98721 T)). M6's step 0 confirmed that for
NRTL's bank alone. thermo and mpmath are given the same forms, so this
fixture tests the models' arithmetic and not the banks' conventions.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_activity.py

It reproduces its output byte for byte.
"""
import datetime
import json
import math
import pathlib
import random
import tomllib

import mpmath as mp
import thermo
from thermo import unifac
from thermo.nrtl import NRTL
from thermo.uniquac import UNIQUAC
from thermo.wilson import Wilson
from thermo.unifac import UNIFAC, UNIFAC_subgroup

ROOT = pathlib.Path(__file__).resolve().parents[1]
VENDOR = ROOT / "data" / "vendor" / "chemsep"
COMPONENTS = ROOT / "data" / "components"
R_DECHEMA = 1.98721
# reflux's R, thermo.r_gas and activity.r_gas: N_A*k_B rounded, 1.84e-11 short
# of the exact 8.31446261815324, which thermo (the Python package) uses.
# mpmath's G^E, H^E and Cp^E are on reflux's; thermo's on its own, which
# the tests scale by the ratio.
R = 8.314462618
R_THERMO = 8.31446261815324
SEED = 6

SYSTEMS = [
    ["67-56-1", "7732-18-5"],  # methanol, water
    ["64-17-5", "7732-18-5"],  # ethanol, water
    ["7732-18-5", "71-36-3"],  # water, 1-butanol
    ["67-64-1", "67-56-1", "7732-18-5"],  # acetone, methanol, water
    ["67-56-1", "64-17-5", "7732-18-5", "67-64-1"],  # methanol, ethanol, water, acetone
    ["110-54-3", "71-43-2", "64-17-5", "67-64-1", "7732-18-5"],  # hexane, benzene, ethanol, acetone, water
    ["71-43-2", "110-82-7"],  # benzene, cyclohexane
]


def record(cas):
    with open(COMPONENTS / f"{cas}.toml", "rb") as f:
        return tomllib.load(f)


def ipd_lines(name):
    """A bank's data lines, in file order: (cas1, cas2, numbers)."""
    text = (VENDOR / name).read_text(encoding="latin-1")
    out, inside = [], False
    for line in text.splitlines():
        t = line.strip()
        if t == "[IPD]":
            inside = True
            continue
        if not inside or not t or t.startswith("#") or t.startswith("Comment=") or t.startswith("Units="):
            continue
        out.append(t.split())
    return out


def bank(name, columns):
    """Each pair's first listing, by (cas1, cas2) as the line has it."""
    pairs = {}
    for f in ipd_lines(name):
        if len(f) < 2 + columns:
            continue
        try:
            vals = [float(v) for v in f[2:2 + columns]]
        except ValueError:
            continue
        key = frozenset((f[0], f[1]))
        if key not in pairs and f[0] != f[1]:
            pairs[key] = (f[0], f[1], vals)
    return pairs


def oriented(b, ci, cj):
    """A pair's values turned to i-j: A_ij first."""
    c1, c2, v = b[frozenset((ci, cj))]
    return v if c1 == ci else [v[1], v[0]] + v[2:]


def covered(b, cas):
    return all(frozenset((cas[i], cas[j])) in b for i in range(len(cas)) for j in range(i + 1, len(cas)))


def unifac_table(main_files, sub_files):
    """ChemSep's table by its own numbers: subgroups, and parameters m -> n -> a or (a, b, c)."""
    params = {}
    for coef, pair in enumerate(main_files):
        for name in pair:
            for i, line in enumerate((VENDOR / name).read_text(encoding="latin-1").splitlines()):
                if i == 0:
                    continue
                f = line.split("#")[0].split()
                if len(f) != 6:
                    continue
                m, n = int(f[0]), int(f[2])
                for a, b2, v in [(m, n, float(f[4])), (n, m, float(f[5]))]:
                    cur = params.setdefault(a, {}).setdefault(b2, [None, 0.0, 0.0])
                    if cur[coef] is None or (coef > 0 and cur[coef] == 0.0):
                        cur[coef] = v
    subs = {}
    for name in sub_files:
        for i, line in enumerate((VENDOR / name).read_text(encoding="latin-1").splitlines()):
            if i == 0:
                continue
            f = line.split("#")[0].split()
            if len(f) < 7:
                continue
            sid, main = int(f[0]), int(f[2])
            if sid not in subs:
                subs[sid] = UNIFAC_subgroup(sid, f[1], main, f[3], float(f[4]), float(f[5]))
    return subs, params


TABLES = {
    "UNIFAC": ([["unifacvl.gct", "unifacvl_ext.gct"]], ["unifacrq.gct", "unifacrq_ext.gct"], 0, "unifac_vle"),
    "UNIFAC-LLE": ([["unifacll.gct", "unifacll_ext.gct"]], ["llefacrq.gct", "llefacrq_ext.gct"], 0, "unifac_lle"),
    "UNIFAC-Dortmund": ([["modfac_a.gct", "modfac_a_ext.gct"], ["modfac_b.gct", "modfac_b_ext.gct"],
                         ["modfac_c.gct", "modfac_c_ext.gct"]], ["modfacrq.gct", "modfacrq_ext.gct"], 1,
                        "modified_unifac"),
}


def states(rng, nc, count):
    out = []
    for k in range(count):
        t = 280.0 + 120.0 * rng.random()
        w = [-math.log(1.0 - rng.random()) for _ in range(nc)]
        if k == count - 1 and nc > 2:
            w[rng.randrange(nc)] = 0.0  # one component absent: its gamma at infinite dilution
        s = sum(w)
        out.append((round(t, 6), [round(v / s, 12) for v in w]))
    # The last fraction makes the sum exact in floating point as far as it can.
    return [(t, x[:-1] + [1.0 - sum(x[:-1])]) for t, x in out]


# An answer a source can't give: no values, and ok false, rather than null,
# which the V side's json2 would need an Option for.
NONE = {"ok": False, "ln_gamma": [], "dt": [], "dtt": [], "dn": [], "ge": 0.0, "he": 0.0, "cp": 0.0}


def thermo_answer(model, t, x):
    """thermo's answer, or None where it can't give one: its UNIQUAC divides
    by each x_i, so it fails with a component absent."""
    try:
        return thermo_values(model, t, x)
    except ZeroDivisionError:
        return NONE


def thermo_values(model, t, x):
    nc = len(x)
    lng = [math.log(g) for g in model.gammas()]
    hess = model.d2nGE_dninjs()
    dn = [hess[i][j] / (R_THERMO * t) for i in range(nc) for j in range(nc)]
    # thermo's UNIFAC has the second T-derivative, and its NRTL, UNIQUAC and
    # Wilson don't: theirs are mpmath's alone.
    dtt = list(model.d2lngammas_dT2()) if hasattr(model, "d2lngammas_dT2") else []
    return {"ok": True, "ln_gamma": lng, "dt": list(model.dlngammas_dT()), "dtt": dtt, "dn": dn,
            "ge": model.GE(), "he": model.HE(), "cp": model.CpE()}


# The closed forms, at mpmath's precision: ln gamma as a function of T and
# the mole numbers (x = n/sum n), each parameter a function of T.

def mp_nrtl(tau, alpha):
    def lng(i, t, n):
        s = mp.fsum(n)
        x = [v / s for v in n]
        nc = len(x)
        tt = [[tau[a][b](t) for b in range(nc)] for a in range(nc)]
        g = [[mp.exp(-alpha[a][b] * tt[a][b]) for b in range(nc)] for a in range(nc)]
        bsum = [mp.fsum(x[k] * g[k][j] for k in range(nc)) for j in range(nc)]
        asum = [mp.fsum(x[k] * tt[k][j] * g[k][j] for k in range(nc)) for j in range(nc)]
        return asum[i] / bsum[i] + mp.fsum(x[j] * g[i][j] / bsum[j] * (tt[i][j] - asum[j] / bsum[j]) for j in range(nc))
    return lng


def mp_uniquac(tau, r, q):
    def lng(i, t, n):
        s = mp.fsum(n)
        x = [v / s for v in n]
        nc = len(x)
        tt = [[tau[a][b](t) for b in range(nc)] for a in range(nc)]
        rbar = mp.fsum(x[j] * r[j] for j in range(nc))
        qbar = mp.fsum(x[j] * q[j] for j in range(nc))
        lj = [5 * (r[j] - q[j]) - (r[j] - 1) for j in range(nc)]
        th = [x[j] * q[j] / qbar for j in range(nc)]
        sj = [mp.fsum(th[k] * tt[k][j] for k in range(nc)) for j in range(nc)]
        v = r[i] / rbar
        f = q[i] / qbar
        return (mp.log(v) + 5 * q[i] * mp.log(f / v) + lj[i] - v * mp.fsum(x[j] * lj[j] for j in range(nc))
                + q[i] * (1 - mp.log(sj[i]) - mp.fsum(th[j] * tt[i][j] / sj[j] for j in range(nc))))
    return lng


def mp_wilson(lam):
    def lng(i, t, n):
        s = mp.fsum(n)
        x = [v / s for v in n]
        nc = len(x)
        ll = [[lam[a][b](t) for b in range(nc)] for a in range(nc)]
        sk = [mp.fsum(x[j] * ll[k][j] for j in range(nc)) for k in range(nc)]
        return 1 - mp.log(sk[i]) - mp.fsum(x[k] * ll[k][i] / sk[k] for k in range(nc))
    return lng


def mp_answer(lng, t, x):
    nc = len(x)
    mt = mp.mpf(t)
    mx = [mp.mpf(v) for v in x]
    out = {"ok": True, "ln_gamma": [], "dt": [], "dtt": [], "dn": []}
    for i in range(nc):
        out["ln_gamma"].append(float(lng(i, mt, mx)))
        out["dt"].append(float(mp.diff(lambda tv: lng(i, tv, mx), mt, 1)))
        out["dtt"].append(float(mp.diff(lambda tv: lng(i, tv, mx), mt, 2)))
    for i in range(nc):
        for j in range(nc):
            order = tuple(1 if k == j else 0 for k in range(nc))
            # n*d/dn_j at n = x, whose total is 1.
            out["dn"].append(float(mp.diff(lambda *ns: lng(i, mt, list(ns)), mx, order)))
    # The excess properties from these: G^E = RT sum x ln gamma, H^E = -RT^2
    # sum x d ln gamma/dT, Cp^E its T-derivative.
    g = mp.fsum(mx[i] * lng(i, mt, mx) for i in range(nc))
    d1 = mp.fsum(mx[i] * mp.diff(lambda tv: lng(i, tv, mx), mt, 1) for i in range(nc))
    d2 = mp.fsum(mx[i] * mp.diff(lambda tv: lng(i, tv, mx), mt, 2) for i in range(nc))
    out["ge"] = float(R * mt * g)
    out["he"] = float(-R * mt * mt * d1)
    out["cp"] = float(-R * mt * (2 * d1 + mt * d2))
    return out


def param(f):
    return lambda t: f(t)


def binary_cases(rng, name, b, cas_systems, per):
    cases = []
    for cas in cas_systems:
        if not covered(b, cas):
            continue
        recs = [record(c) for c in cas]
        nc = len(cas)
        p = {(i, j): oriented(b, cas[i], cas[j]) for i in range(nc) for j in range(nc) if i != j}
        for t, x in states(rng, nc, per):
            if name == "NRTL":
                bs = [[0.0 if i == j else p[i, j][0] / R_DECHEMA for j in range(nc)] for i in range(nc)]
                al = [[0.0 if i == j else p[i, j][2] for j in range(nc)] for i in range(nc)]
                model = NRTL(xs=x, T=t, tau_bs=bs, alpha_cs=al)
                taus = [[(lambda t, a=0.0 if i == j else p[i, j][0]: mp.mpf(a) / (mp.mpf(R_DECHEMA) * t))
                         for j in range(nc)] for i in range(nc)]
                ref = mp_answer(mp_nrtl(taus, [[mp.mpf(al[i][j]) for j in range(nc)] for i in range(nc)]), t, x)
            elif name == "UNIQUAC":
                r = [rec["parameters"]["uniquac_r"] for rec in recs]
                q = [rec["parameters"]["uniquac_q"] for rec in recs]
                bs = [[0.0 if i == j else -p[i, j][0] / R_DECHEMA for j in range(nc)] for i in range(nc)]
                model = UNIQUAC(xs=x, T=t, rs=r, qs=q, tau_bs=bs)
                taus = [[(lambda t, a=0.0 if i == j else p[i, j][0]: mp.exp(-mp.mpf(a) / (mp.mpf(R_DECHEMA) * t)))
                         for j in range(nc)] for i in range(nc)]
                ref = mp_answer(mp_uniquac(taus, [mp.mpf(v) for v in r], [mp.mpf(v) for v in q]), t, x)
            else:
                v = [rec["parameters"]["wilson_volume"] for rec in recs]
                la = [[0.0 if i == j else math.log(v[j] / v[i]) for j in range(nc)] for i in range(nc)]
                bs = [[0.0 if i == j else -p[i, j][0] / R_DECHEMA for j in range(nc)] for i in range(nc)]
                model = Wilson(xs=x, T=t, lambda_as=la, lambda_bs=bs)
                lams = [[(lambda t, a=0.0 if i == j else p[i, j][0], i=i, j=j:
                          mp.mpf(1) if i == j else mp.mpf(v[j]) / mp.mpf(v[i]) * mp.exp(-mp.mpf(a) / (mp.mpf(R_DECHEMA) * t)))
                         for j in range(nc)] for i in range(nc)]
                ref = mp_answer(mp_wilson(lams), t, x)
            cases.append({"model": name, "cas": cas, "T": t, "x": x, "thermo": thermo_answer(model, t, x),
                          "mpmath": ref})
    return cases


def unifac_cases(rng, per):
    cases = []
    for name, (mains, subs_files, version, key) in TABLES.items():
        subs, params = unifac_table(mains, subs_files)
        interaction = {m: {n: (tuple(v) if version == 1 else v[0]) for n, v in row.items()} for m, row in params.items()}
        for cas in SYSTEMS:
            groups = [{int(k): v for k, v in record(c).get("groups", {}).get(key, {}).items()} for c in cas]
            if any(not g or any(k not in subs for k in g) for g in groups):
                continue
            mains_present = {subs[k].main_group_id for g in groups for k in g}
            if any(n not in interaction.get(m, {}) for m in mains_present for n in mains_present if m != n):
                continue
            for t, x in states(rng, len(cas), per):
                model = UNIFAC.from_subgroups(T=t, xs=x, chemgroups=groups, subgroups=subs,
                                              interaction_data=interaction, version=version)
                cases.append({"model": name, "cas": cas, "T": t, "x": x, "thermo": thermo_answer(model, t, x),
                              "mpmath": NONE})
    return cases


def tables():
    def subgroups(t):
        return [{"id": k, "name": s.group, "main": s.main_group_id, "main_name": s.main_group, "R": s.R, "Q": s.Q}
                for k, s in sorted(t.items())]

    def interactions(t):
        out = []
        for m in sorted(t):
            for n in sorted(t[m]):
                v = t[m][n]
                out.append({"m": m, "n": n, "abc": list(v) if isinstance(v, tuple) else [v, 0.0, 0.0]})
        return out

    return {
        "UNIFAC": {"subgroups": subgroups(unifac.UFSG), "interactions": interactions(unifac.UFIP)},
        "UNIFAC-LLE": {"subgroups": subgroups(unifac.LLEUFSG), "interactions": interactions(unifac.LLEUFIP)},
        "UNIFAC-Dortmund": {"subgroups": subgroups(unifac.DOUFSG), "interactions": interactions(unifac.DOUFIP2016)},
        "UNIFAC-Dortmund-2006": {"subgroups": [], "interactions": interactions(unifac.DOUFIP2006)},
    }


def main():
    mp.mp.dps = 40
    rng = random.Random(SEED)
    cases = []
    cases += binary_cases(rng, "NRTL", bank("nrtl.ipd", 3), SYSTEMS, 4)
    cases += binary_cases(rng, "UNIQUAC", bank("uniquac.ipd", 2), SYSTEMS, 4)
    cases += binary_cases(rng, "Wilson", bank("wilson.ipd", 2), SYSTEMS, 4)
    cases += unifac_cases(rng, 4)
    source = {"generator": "tools/gen_m6_activity.py", "thermo": thermo.__version__, "mpmath": mp.__version__,
              "generated": datetime.date.today().isoformat(), "seed": SEED, "R_thermo": R_THERMO,
              "R_mpmath": R}
    out = {"source": source,
           "note": "Activity models on ChemSep's banks and tables: thermo's answers, and mpmath's from the closed "
                   "forms at 40 digits for NRTL, UNIQUAC and Wilson. dn is n*d(ln gamma_i)/dn_j, row-major; ge and "
                   "he are J/mol, cp J/(mol K).",
           "cases": cases}
    (ROOT / "testdata" / "m6_activity.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8",
                                                        newline="\n")
    tab = {"source": source, "note": "thermo's copies of the published UNIFAC tables, by their own numbers: the "
                                     "original (UFSG, UFIP), the LLE table (LLEUFSG, LLEUFIP), Dortmund's 2016 "
                                     "(DOUFSG, DOUFIP2016), and Dortmund's 2006 parameters (DOUFIP2006, with "
                                     "DOUFSG's subgroups). abc is a, b, c; the original's and LLE's b and c are 0.",
           **tables()}
    (ROOT / "testdata" / "m6_unifac_tables.json").write_text(json.dumps(tab, indent=1) + "\n", encoding="utf-8",
                                                             newline="\n")
    counts = {}
    for c in cases:
        counts[c["model"]] = counts.get(c["model"], 0) + 1
    print("wrote testdata/m6_activity.json:", counts)
    print("wrote testdata/m6_unifac_tables.json")


if __name__ == "__main__":
    main()
