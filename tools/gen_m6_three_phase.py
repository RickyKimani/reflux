# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["thermo==0.6.1", "chemicals==1.5.2", "mpmath==1.3.0", "scipy==1.14.1"]
# ///
"""Write M6's three-phase references (step 7, criterion 6):
testdata/m6_three_phase.json.

thermo 0.6.1's FlashVLN, with two liquids and a gas, at each state:
- NRTL (ChemSep's bank, the 'Poynting' basis, the ideal gas) on water and
  1-butanol, two liquids, and on water, 1-butanol and methanol, two liquids
  and, at 1 atm near 364 K, a vapour and two liquids;
- Peng-Robinson on water, n-decane and methane (flash/pt_test.v's constants
  and k_ij), a vapour and two liquids, and on water and n-decane, two liquids.

Each answer is recorded as FlashVLN gives it ("flashvln"), and polished on
thermo's own phases ("phases"): successive substitution from FlashVLN's
answer, each step's fractions from Michelsen's multiphase Rachford-Rice
(numerical minimisation of Q over beta >= 0), until every |ln f_ij - ln f_i0|
is below 1e-14. The phases are in order of increasing mass density, each
with thermo's kind (gas or liquid).

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m6_three_phase.py

It reproduces its output byte for byte on the same day.
"""
import datetime
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m6_activity as ga  # noqa: E402
import gen_m6_gamma_phi as gp  # noqa: E402

from thermo import (ChemicalConstantsPackage, PropertyCorrelationsPackage, GibbsExcessLiquid,  # noqa: E402
                    IdealGas, FlashVLN, CEOSGas, CEOSLiquid, PRMIX, HeatCapacityGas)

ROOT = ga.ROOT

# flash/pt_test.v's water, n-decane and methane.
PR_COMPONENTS = [
    {"name": "water", "tc": 647.1, "pc": 220.55e5, "omega": 0.345, "mw": 0.018015},
    {"name": "n-decane", "tc": 617.7, "pc": 21.10e5, "omega": 0.492, "mw": 0.142285},
    {"name": "methane", "tc": 190.6, "pc": 45.99e5, "omega": 0.012, "mw": 0.016043},
]
PR_KIJ = [[0.0, 0.5, 0.5], [0.5, 0.0, 0.04], [0.5, 0.04, 0.0]]


def gamma_phi_flasher(cas):
    """FlashVLN on the gamma-phi liquid, twice, and the ideal gas."""
    flasher, liq = gp.parts("NRTL", cas)
    gas = flasher.gas
    return FlashVLN(flasher.constants, flasher.correlations, liquids=[liq, liq], gas=gas)


def pr_flasher(comps):
    n = len(comps)
    tcs = [c["tc"] for c in comps]
    pcs = [c["pc"] for c in comps]
    omegas = [c["omega"] for c in comps]
    kijs = [row[:n] for row in PR_KIJ[:n]]
    consts = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=omegas, MWs=[c["mw"] * 1e3 for c in comps],
                                      names=[c["name"] for c in comps])
    cps = []
    for _ in comps:
        cp = HeatCapacityGas()
        cp.add_correlation(name="C", model="constant", Tmin=10.0, Tmax=5000.0, value=40.0)
        cp.method = "C"
        cps.append(cp)
    corr = PropertyCorrelationsPackage(constants=consts, HeatCapacityGases=cps, skip_missing=True)
    kw = dict(Tcs=tcs, Pcs=pcs, omegas=omegas, kijs=kijs)
    x0 = [1.0 / n] * n
    gas = CEOSGas(PRMIX, eos_kwargs=kw, HeatCapacityGases=cps, T=300.0, P=1e5, zs=x0)
    liq = CEOSLiquid(PRMIX, eos_kwargs=kw, HeatCapacityGases=cps, T=300.0, P=1e5, zs=x0)
    return FlashVLN(consts, corr, liquids=[liq, liq], gas=gas)


def multiphase_rr(z, lnk, beta):
    """Michelsen's Q(beta) = sum beta - sum z ln E, E_i = sum_j beta_j K_ij,
    minimised over beta >= 0 by Newton's method on the free phases."""
    n, npp = len(z), len(lnk)
    k = []
    for j in range(npp):
        k.append([0.0] * n)
    for i in range(n):
        top = max(lnk[j][i] for j in range(npp))
        for j in range(npp):
            k[j][i] = math.exp(lnk[j][i] - top)
    beta = [max(b, 0.0) for b in beta]
    for _ in range(200):
        e = [sum(beta[j] * k[j][i] for j in range(npp)) for i in range(n)]
        g = [1.0 - sum(z[i] * k[j][i] / e[i] for i in range(n)) for j in range(npp)]
        free = [j for j in range(npp) if beta[j] > 0.0 or g[j] < 0.0]
        if all(abs(g[j]) < 1e-15 for j in free) and all(g[j] > -1e-15 for j in range(npp) if j not in free):
            break
        while True:
            h = [[sum(z[i] * k[a][i] * k[b][i] / e[i] ** 2 for i in range(n)) for b in free] for a in free]
            d = solve(h, [-g[j] for j in free])
            held = [free[a] for a in range(len(free)) if beta[free[a]] == 0.0 and d[a] < 0.0]
            if not held:
                break
            free = [j for j in free if j not in held]
        alpha = 1.0
        for a, j in enumerate(free):
            if d[a] < 0.0:
                alpha = min(alpha, -beta[j] / d[a])
        for a, j in enumerate(free):
            beta[j] = max(0.0, beta[j] + alpha * d[a])
    e = [sum(beta[j] * k[j][i] for j in range(npp)) for i in range(n)]
    x = []
    for j in range(npp):
        xj = [z[i] * k[j][i] / e[i] for i in range(n)]
        s = sum(xj)
        x.append([v / s for v in xj])
    total = sum(beta)
    return [b / total for b in beta], x


def solve(a, b):
    """Gaussian elimination with partial pivoting, for a small system."""
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        for r in range(c + 1, n):
            f = m[r][c] / m[c][c]
            for q in range(c, n + 1):
                m[r][q] -= f * m[c][q]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (m[r][n] - sum(m[r][q] * x[q] for q in range(r + 1, n))) / m[r][r]
    return x


def polish(res, t, p, z, mws):
    """Successive substitution on thermo's phases from FlashVLN's answer."""
    phases = list(res.phases)
    x = [list(ph.zs) for ph in phases]
    beta = list(res.betas)
    for _ in range(2000):
        cur = [ph.to(T=t, P=p, zs=xj) for ph, xj in zip(phases, x)]
        lnphis = [c.lnphis() for c in cur]
        worst = max(abs(math.log(x[j][i]) + lnphis[j][i] - math.log(x[0][i]) - lnphis[0][i])
                    for j in range(1, len(x)) for i in range(len(z)) if x[j][i] > 1e-300 and x[0][i] > 1e-300)
        if worst < 1e-14:
            break
        beta, x = multiphase_rr(z, [[-v for v in lp] for lp in lnphis], beta)
    cur = [ph.to(T=t, P=p, zs=xj) for ph, xj in zip(phases, x)]
    out = []
    for c, b, xj in zip(cur, beta, x):
        out.append({"kind": "gas" if isinstance(c, (IdealGas, CEOSGas)) else "liquid", "fraction": b, "x": xj,
                    "rho_mass": sum(xi * mw for xi, mw in zip(xj, mws)) / c.V()})
    out.sort(key=lambda o: o["rho_mass"])
    return out, worst


def case(model, system, flasher, t, p, z):
    mws = [mw / 1e3 for mw in flasher.constants.MWs]  # kg/mol
    res = flasher.flash(T=t, P=p, zs=z)
    raw = {"phase_count": res.phase_count, "betas": list(res.betas), "zs": [list(ph.zs) for ph in res.phases]}
    phases, worst = polish(res, t, p, z, mws) if res.phase_count > 1 else ([{
        "kind": "gas" if res.gas is not None else "liquid", "fraction": 1.0, "x": list(z),
        "rho_mass": sum(xi * mw for xi, mw in zip(z, mws)) / res.phases[0].V()}], 0.0)
    return {"model": model, "system": system, "T": t, "P": p, "z": z, "flashvln": raw, "phases": phases,
            "residual": worst}


def main():
    out = {
        "source": {"generator": "tools/gen_m6_three_phase.py", "thermo": "0.6.1",
                   "generated": datetime.date.today().isoformat()},
        "note": ("thermo's FlashVLN (two liquids and a gas) as given, and polished on thermo's phases by "
                 "substitution on Michelsen's multiphase Rachford-Rice to 1e-14 in ln f; phases in order of "
                 "increasing mass density."),
        "pr_components": PR_COMPONENTS,
        "pr_kij": PR_KIJ,
        "cases": [],
    }
    wb = ["7732-18-5", "71-36-3"]
    f = gamma_phi_flasher(wb)
    for t in (300.0, 320.0, 340.0):
        for z in ([0.6, 0.4], [0.8, 0.2]):
            out["cases"].append(case("NRTL", wb, f, t, 101325.0, z))
    wbm = ["7732-18-5", "71-36-3", "67-56-1"]
    f = gamma_phi_flasher(wbm)
    for t in (300.0, 330.0, 363.6, 364.2, 364.8):
        out["cases"].append(case("NRTL", wbm, f, t, 101325.0, [0.75, 0.2, 0.05]))
    f = pr_flasher(PR_COMPONENTS)
    for t, p in ((300.0, 10e5), (350.0, 5e5), (400.0, 20e5), (330.0, 30e5)):
        out["cases"].append(case("PR", ["water", "n-decane", "methane"], f, t, p, [0.5, 0.3, 0.2]))
    f = pr_flasher(PR_COMPONENTS[:2])
    for t, p in ((300.0, 10e5), (400.0, 50e5)):
        out["cases"].append(case("PR", ["water", "n-decane"], f, t, p, [0.6, 0.4]))
    (ROOT / "testdata" / "m6_three_phase.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8",
                                                          newline="\n")
    for c in out["cases"]:
        print(c["model"], c["system"], c["T"], c["P"], "FlashVLN", c["flashvln"]["phase_count"], "phases;",
              [(ph["kind"], round(ph["fraction"], 6)) for ph in c["phases"]], "residual", c["residual"])


if __name__ == "__main__":
    main()
