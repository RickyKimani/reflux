# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["mpmath==1.3.0", "chemicals==1.5.2", "thermo==0.6.1", "rdkit==2024.9.6"]
# ///
"""Generate M9's fixture for chemical equilibrium (step 7).

Run from the repository root:

    uv run --python 3.13 tools/gen_m9_equilibrium.py

It writes testdata/m9_equilibrium.json: equilibria of one gas phase, ideal
or Peng–Robinson, each solved here at 50 digits, independently of reflux:

- **The equations:** every species' amount in log space, xⱼ = ln nⱼ, and
  for a basis of the reactions, each one's Σⱼ νⱼ μⱼ/RT = 0, that is ln Q =
  ln K, with what the reactions conserve held at the feed's. A Gibbs
  minimisation is the same over a basis of every reaction the elements
  allow, the null space of the element matrix, found with exact fractions;
  an equilibrium over chosen reactions is over those alone. μⱼ/RT =
  μ°ⱼ/RT + ln(P/p°) + ln yⱼ + ln φⱼ, with μ°ⱼ = H°ⱼ − TS°ⱼ from ChemSep's
  data by mpmath (tools/gen_m9_reactions.py's), R exact, p° = 1 bar.
- **The solver:** Newton on x with a central-difference Jacobian at 50
  digits, halving its step until the residual falls, to a residual below
  1e-40. That's not reflux's method (RAND on element potentials), so the
  two meet only at the answer.
- **Peng–Robinson,** written out here: Ωa and Ωb from the critical-point
  cubic (Michelsen and Mollerup), Soave's α with Peng and Robinson's m(ω),
  van der Waals mixing with every kᵢⱼ zero, the root of lower Gibbs energy,
  and the textbook ln φᵢ; R for the cubic 8.314462618, thermo's.

The output is the same, byte for byte, for the same versions.
"""

import json
import pathlib
import sys
import tomllib
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_m9_reactions as rx  # noqa: E402

mp.mp.dps = 50

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "components"
OUT = ROOT / "testdata" / "m9_equilibrium.json"

R = rx.R
R_PR = mp.mpf("8.314462618")
P0 = mp.mpf(100000)

CH4, H2O, H2, CO, CO2 = "74-82-8", "7732-18-5", "1333-74-0", "630-08-0", "124-38-9"
MEOH, N2, NH3, O2 = "67-56-1", "7727-37-9", "7664-41-7", "7782-44-7"

FORMULAS = {
    CH4: {"C": 1, "H": 4}, H2O: {"H": 2, "O": 1}, H2: {"H": 2}, CO: {"C": 1, "O": 1},
    CO2: {"C": 1, "O": 2}, MEOH: {"C": 1, "H": 4, "O": 1}, N2: {"N": 2}, NH3: {"N": 1, "H": 3},
    O2: {"O": 2},
}
SYMBOLS = {CH4: "CH4", H2O: "H2O", H2: "H2", CO: "CO", CO2: "CO2", MEOH: "CH3OH", N2: "N2", NH3: "NH3",
           O2: "O2"}

REFORMING = [CH4, H2O, H2, CO, CO2]
METHANOL = [CO2, H2, MEOH, H2O, CO]
AMMONIA = [N2, H2, NH3]
COMBUSTION = [CH4, O2, CO2, H2O, CO, H2, N2]

SR = ("CH4 + H2O = CO + 3 H2", {CH4: -1, H2O: -1, CO: 1, H2: 3})
CO2_TO_METHANOL = ("CO2 + 3 H2 = CH3OH + H2O", {CO2: -1, H2: -3, MEOH: 1, H2O: 1})


def case(name, species, model, t, p, n0, reactions=None, inert=()):
    return {"name": name, "species": species, "model": model, "t": t, "p": p, "n0": n0,
            "reactions": reactions, "inert": list(inert)}


LABEL = {"ideal": "ideal", "pr": "PR"}

CASES = [
    *[case(f"reforming, ideal, {t} K, {p} bar", REFORMING, "ideal", t, p, [1, 2, 0, 0, 0])
      for t in (700, 900, 1100) for p in (1, 20)],
    *[case(f"reforming, PR, 900 K, {p} bar", REFORMING, "pr", 900, p, [1, 2, 0, 0, 0]) for p in (20, 30)],
    case("reforming alone, ideal, 900 K, 1 bar", REFORMING, "ideal", 900, 1, [1, 2, 0, 0, 0], [SR]),
    case("reforming alone, PR, 900 K, 20 bar", REFORMING, "pr", 900, 20, [1, 2, 0.05, 0, 0.1], [SR]),
    case("reforming, CO2 inert, ideal, 900 K, 1 bar", REFORMING, "ideal", 900, 1, [1, 2, 0, 0, 0.5],
         inert=[4]),
    *[case(f"methanol, {LABEL[m]}, 500 K, {p} bar", METHANOL, m, 500, p, [1, 3, 0, 0, 0.1])
      for m in ("ideal", "pr") for p in (50, 80)],
    case("methanol from CO2 alone, PR, 500 K, 80 bar", METHANOL, "pr", 500, 80, [1, 3, 0, 0, 0.1],
         [CO2_TO_METHANOL]),
    *[case(f"ammonia, {LABEL[m]}, 700 K, {p} bar", AMMONIA, m, 700, p, [1, 3, 0]) for m in ("ideal", "pr")
      for p in (100, 300)],
    case("lean methane combustion, ideal, 1500 K, 1 bar", COMBUSTION, "ideal", 1500, 1,
         [1, 3, 0, 0, 0, 0, 11.28]),
]


def mu0(cas, t):
    s = rx.chemsep(cas)
    h = s["dhf"] + mp.quad(s["cp"], [rx.T0, t])
    e = s["s0"] + mp.quad(lambda x: s["cp"](x) / x, [rx.T0, t])
    return (h - t * e) / (R * t)


def constants(cas):
    d = tomllib.loads((DATA / f"{cas}.toml").read_text(encoding="utf-8"))["constants"]
    return mp.mpf(float(d["tc"])), mp.mpf(float(d["pc"])), mp.mpf(float(d["omega"]))


def pr_omegas():
    """Ωa and Ωb of Peng–Robinson from its critical point, as Michelsen and
    Mollerup derive them: δ₁ = 1 + √2, δ₂ = 1 − √2."""
    s, m = mp.mpf(2), mp.mpf(-1)

    def at(b):
        zc = (1 - (s - 1) * b) / 3
        a = 3 * zc * zc - m * b * b + s * b * (b + 1)
        return m * b * b * (b + 1) + a * b - zc ** 3

    lo, hi = mp.mpf("1e-6"), mp.mpf("0.3")
    for _ in range(300):
        mid = (lo + hi) / 2
        if (at(lo) < 0) == (at(mid) < 0):
            lo = mid
        else:
            hi = mid
    b = (lo + hi) / 2
    zc = (1 - (s - 1) * b) / 3
    return 3 * zc * zc - m * b * b + s * b * (b + 1), b


OMEGA_A, OMEGA_B = pr_omegas()
D1, D2 = 1 + mp.sqrt(2), 1 - mp.sqrt(2)


def pr_ln_phi(species, t, p, y):
    """ln φᵢ by Peng–Robinson, every kᵢⱼ zero, at the root of lower Gibbs
    energy."""
    ac, bs = [], []
    for cas in species:
        tc, pc, w = constants(cas)
        m = mp.mpf("0.37464") + mp.mpf("1.54226") * w - mp.mpf("0.26992") * w * w
        alpha = (1 + m * (1 - mp.sqrt(t / tc))) ** 2
        ac.append(OMEGA_A * R_PR ** 2 * tc ** 2 / pc * alpha)
        bs.append(OMEGA_B * R_PR * tc / pc)
    n = len(species)
    aij = [[mp.sqrt(ac[i] * ac[j]) for j in range(n)] for i in range(n)]
    a = sum(y[i] * y[j] * aij[i][j] for i in range(n) for j in range(n))
    b = sum(y[i] * bs[i] for i in range(n))
    A = a * p / (R_PR * t) ** 2
    B = b * p / (R_PR * t)
    roots = mp.polyroots([1, -(1 - B), A - 3 * B * B - 2 * B, -(A * B - B * B - B ** 3)], maxsteps=200,
                         extraprec=200)
    zs = [mp.re(z) for z in roots if abs(mp.im(z)) < mp.mpf("1e-30") and mp.re(z) > B]
    best = None
    for z in zs:
        lp = [bs[i] / b * (z - 1) - mp.log(z - B)
              - A / (B * (D1 - D2)) * (2 * sum(y[j] * aij[i][j] for j in range(n)) / a - bs[i] / b)
              * mp.log((z + D1 * B) / (z + D2 * B)) for i in range(n)]
        g = sum(y[i] * lp[i] for i in range(n))
        if best is None or g < best[0]:
            best = (g, lp)
    return best[1]


def null_space(rows, n):
    """A basis of the vectors v with row·v = 0 for every row, exactly."""
    rows = [list(r) for r in rows]
    pivots = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        lead = rows[r][c]
        rows[r] = [x / lead for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
    basis = []
    for f in range(n):
        if f in pivots:
            continue
        v = [Fraction(0)] * n
        v[f] = Fraction(1)
        for i, p in enumerate(pivots):
            v[p] = -rows[i][f]
        basis.append(v)
    return basis


def solve(c):
    species, t, p = c["species"], mp.mpf(c["t"]), mp.mpf(c["p"]) * P0
    n0 = [mp.mpf(repr(float(x))) for x in c["n0"]]
    s = len(species)
    if c["reactions"]:
        nu = [[Fraction(r[1].get(cas, 0)) for cas in species] for r in c["reactions"]]
    else:
        elements = sorted({e for cas in species for e in FORMULAS[cas]})
        rows = [[Fraction(FORMULAS[cas].get(e, 0)) if j not in c["inert"] else Fraction(0)
                 for j, cas in enumerate(species)] for e in elements]
        rows += [[Fraction(1 if j == i else 0) for j in range(s)] for i in c["inert"]]
        nu = null_space(rows, s)
    # A species no reaction touches, which the feed hasn't, stays at zero.
    touched = [any(r[j] != 0 for r in nu) for j in range(s)]
    active = [j for j in range(s) if touched[j] or n0[j] > 0]
    nu_a = [[r[j] for j in active] for r in nu]
    cons = null_space(nu_a, len(active))
    g = {j: mu0(species[j], t) for j in active if touched[j]}
    lnp = mp.log(p / P0)

    def residual(x):
        n = [mp.exp(v) for v in x]
        total = sum(n)
        y = [v / total for v in n]
        phi = [mp.mpf(0)] * len(active)
        if c["model"] == "pr":
            phi = pr_ln_phi([species[j] for j in active], t, p, y)
        mu = [(g.get(j, 0) + lnp + mp.log(y[q]) + phi[q]) for q, j in enumerate(active)]
        f = [sum(r[q] * mu[q] for q in range(len(active))) for r in nu_a]
        for row in cons:
            want = sum(row[q] * n0[j] for q, j in enumerate(active))
            scale = max(abs(row[q] * n0[j]) for q, j in enumerate(active)) or 1
            f.append((sum(row[q] * n[q] for q in range(len(active))) - want) / scale)
        return f

    x = [mp.log(sum(n0) / len(active))] * len(active)
    if c.get("start"):
        x = c["start"]
    fx = residual(x)
    norm = max(abs(v) for v in fx)
    h = mp.mpf("1e-20")
    for _ in range(400):
        if norm < mp.mpf("1e-40"):
            break
        jac = mp.matrix(len(fx), len(x))
        for i in range(len(x)):
            up = list(x)
            dn = list(x)
            up[i] += h
            dn[i] -= h
            fu, fd = residual(up), residual(dn)
            for k in range(len(fx)):
                jac[k, i] = (fu[k] - fd[k]) / (2 * h)
        dx = mp.lu_solve(jac, mp.matrix([-v for v in fx]))
        step = max(abs(v) for v in dx)
        lam = min(mp.mpf(1), 5 / step) if step > 0 else mp.mpf(1)
        while True:
            trial = [x[i] + lam * dx[i] for i in range(len(x))]
            ft = residual(trial)
            nt = max(abs(v) for v in ft)
            if nt < norm or lam < mp.mpf("1e-12"):
                break
            lam /= 2
        x, fx, norm = trial, ft, nt
    if norm >= mp.mpf("1e-40"):
        raise SystemExit(f"{c['name']}: no convergence, residual {norm}")
    n = [mp.mpf(0)] * s
    for q, j in enumerate(active):
        n[j] = mp.exp(x[q])
    total = sum(n)
    return n, [v / total for v in n], x


def main():
    cases = []
    solved = {}
    for c in CASES:
        if c["model"] == "pr":
            # Start from the ideal gas's answer for the same case.
            ideal = dict(c, model="ideal")
            _, _, x = solve(ideal)
            c = dict(c, start=x)
        n, y, _ = solve(c)
        solved[c["name"]] = y
        cases.append({
            "name": c["name"],
            "species": c["species"],
            "symbols": [SYMBOLS[cas] for cas in c["species"]],
            "model": c["model"],
            "reactions": [r[0] for r in c["reactions"]] if c["reactions"] else [],
            "inert": c["inert"],
            "t": float(c["t"]),
            "p": float(mp.mpf(c["p"]) * P0),
            "n0": [float(v) for v in c["n0"]],
            "n": [float(v) for v in n],
            "ln_y": [float(mp.log(v)) if v > 0 else None for v in y],
        })
        print(f"{c['name']}: y = {[mp.nstr(v, 8) for v in y]}")
    data = {
        "source": {
            "by": f"tools/gen_m9_equilibrium.py, with mpmath {mp.__version__} at 50 digits",
            "method": "each case solved for every species' ln n by Newton: each reaction of a basis at ln Q = ln K "
                      "(a basis of every reaction the elements allow, for a Gibbs minimisation; the chosen reactions "
                      "alone, otherwise) and what they conserve at the feed's, to a residual below 1e-40",
            "data": "mu°/RT = (H° - TS°)/RT from ChemSep's data in data/components/, by mpmath; R = "
                    "8.31446261815324, p° = 1 bar; Peng-Robinson written out with every kij zero, Soave's alpha "
                    "with Peng and Robinson's m(omega), R 8.314462618 for the cubic",
            "units": "T in K, P in Pa, amounts in mol; ln_y is each species' ln y (null where it's zero); inert "
                     "lists species held at their feed",
        },
        "cases": cases,
    }
    OUT.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {OUT.relative_to(ROOT)}: {len(cases)} cases")


if __name__ == "__main__":
    main()
