# /// script
# requires-python = ">=3.11"
# dependencies = ["chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Generate the M3 reference fixtures for reflux.

Run from the repository root:

    uv run tools/gen_m3_fixtures.py

It writes to testdata/ the files below. To write only some of them, name
them: `uv run tools/gen_m3_fixtures.py thermo stability` (the others are
pip and rachford_rice).

m3_rachford_rice.json
    Rachford-Rice splits for fixed K-values, read by
    numeric/rachford_rice_test.v. The cases are seeded and cover ordinary
    mixtures of 2 to 20 components, extreme K-values (1e-30 to 1e25), a
    trace liquid or vapour (a phase fraction down to 1e-15), K-values within
    1e-8 of 1, negative flashes (beta outside [0, 1]), and absent
    components.

    The reference is written independently of reflux and of chemicals: the
    root is bracketed by bisection between the two poles nearest zero, at
    100 significant digits, from the exact binary values of z and K. So it
    is the exact answer for the floats given, rounded once. Each case also
    records what chemicals' own Rachford-Rice solvers return, for
    comparison. The script stops if chemicals' mpmath solver disagrees with
    the reference.

m3_pip.json
    Phase identification, read by thermo/phase_test.v. Caleb Bell's `thermo`
    computes PR and SRK states, the same models and constants as reflux,
    for pure components and mixtures on a grid of reduced T and P that
    crosses the saturation curve and the critical point. Each root carries
    its molar volume and thermo's analytic phase identification parameter
    (PIP). Pure components' exact critical points are left out, since the
    PIP is 0/0 there. No state in the grid is past any component's Soave
    zero, where reflux and thermo differ on purpose (see thermo/alpha.v).

    It also holds pure CO2 through the point where Soave's alpha reaches
    zero, T0 = Tc(1 + 1/m)^2, for PR and SRK, with the properties M1's
    fixtures check. For a pure component thermo computes alpha from
    k = 1 + m(1 - sqrt(Tr)) without dividing by it, so it's a valid
    reference right through T0.

m3_stability.json
    Stability, read by flash/stability_test.v: PR (and SRK for one binary)
    states of binaries, a ternary and a 10-component gas condensate, over
    grids of T, P and z, with extra states near the methane/n-butane
    critical locus. For each binary a brute-force scan of the tangent-plane
    distance over the whole composition range gives the global minimum and
    the verdict, independent of any iteration. Every state also records
    thermo's FlashVL result (its stability test plus split), and for a split
    the lowest tm of its phases against the feed, which must be negative.

m3_thermo.json
    The flash references, read by flash/thermo_test.v now and by the tests
    of the PH, PS, TV, UV and vapour-fraction flashes to come. thermo's PR
    (and SRK for one binary), with the stability components plus isobutane,
    and a shared ideal-gas Cp: a quadratic in T, which thermo takes as a
    poly_fit and reflux as a Polynomial, so total H and S compare as exactly
    as the departures do.

    - pt: PT flashes of binaries, the ternary, the gas condensate, states
      near the methane/n-butane critical locus, pure components and a
      narrow-boiling n-butane/isobutane mixture. Each has the split (phases
      sorted by mass density, the vapour first) with its fugacity residual,
      and molar V, H, S and U. These are the specifications the other
      flashes are given, and the PT state is the answer they must recover.
      Each case also records what thermo's own PH, PS and TV flashes return
      from those values, or the error they raise: thermo 0.6.1's FlashVL
      fails on PS and TV and on vapour fractions strictly between 0 and 1
      (an unbound local in flash_vl.py), so the PT state is the reference.
    - saturation: pure components at a given P or T and vapour fraction,
      from FlashPureVLS, which covers what no PT flash can: a pure
      component's two-phase states.
    - boundaries: bubble and dew points of the mixtures at a given P or T,
      with the incipient phase.
"""

import datetime
import json
import math
import pathlib
import platform
import random
import sys

import chemicals
import chemicals.rachford_rice as rr
import mpmath as mp
import thermo
from thermo import (PR, PRMIX, SRK, SRKMIX, CEOSGas, CEOSLiquid, ChemicalConstantsPackage, FlashPureVLS, FlashVL,
                    HeatCapacityGas, PropertyCorrelationsPackage)
from thermo.phases import Phase

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT_RR = ROOT / "testdata" / "m3_rachford_rice.json"
OUT_PIP = ROOT / "testdata" / "m3_pip.json"
OUT_STAB = ROOT / "testdata" / "m3_stability.json"
OUT_THERMO = ROOT / "testdata" / "m3_thermo.json"
SEED = 20261027

# As in tools/gen_m1_fixtures.py (Smith, Van Ness & Abbott, Table B.1),
# plus n-decane from the same table.
COMPONENTS = [
    {"name": "methane", "tc": 190.6, "pc": 45.99e5, "omega": 0.012, "mw": 0.016043},
    {"name": "n-butane", "tc": 425.1, "pc": 37.96e5, "omega": 0.200, "mw": 0.058123},
    {"name": "carbon dioxide", "tc": 304.2, "pc": 73.83e5, "omega": 0.224, "mw": 0.044010},
    {"name": "n-decane", "tc": 617.7, "pc": 21.10e5, "omega": 0.492, "mw": 0.142285},
]
KIJ = [[0.0, 0.02, 0.10, 0.04], [0.02, 0.0, 0.13, 0.0], [0.10, 0.13, 0.0, 0.11], [0.04, 0.0, 0.11, 0.0]]

mp.mp.dps = 100


def reference(zs, ks):
    """The Rachford-Rice root for these exact floats, by bisection at 100
    digits between c1 = 1/(1 - Kmax) and cn = 1/(1 - Kmin), over the
    components present. Returns beta, 1 - beta, x and y, each rounded once."""
    present = [i for i, z in enumerate(zs) if z > 0]
    z = [mp.mpf(v) for v in zs]
    k = [mp.mpf(v) for v in ks]
    kmax = max(k[i] for i in present)
    kmin = min(k[i] for i in present)
    lo, hi = 1 / (1 - kmax), 1 / (1 - kmin)

    def f(beta):
        return mp.fsum(z[i] * (k[i] - 1) / (1 + beta * (k[i] - 1)) for i in present)

    # Each halving gains a bit; 100 digits need about 330 below the window's width.
    width_bits = int(mp.log(hi - lo, 2)) + 1
    for _ in range(width_bits + 340):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    beta = (lo + hi) / 2
    x = [z[i] / (1 + beta * (k[i] - 1)) if z[i] > 0 else mp.mpf(0) for i in range(len(zs))]
    y = [k[i] * x[i] for i in range(len(zs))]
    return {
        "beta": float(beta),
        "liquid": float(1 - beta),
        "x": [float(v) for v in x],
        "y": [float(v) for v in y],
    }


def chemicals_solvers(zs, ks):
    """What each of chemicals' two-phase solvers returns, or the error it raises."""
    out = []
    solvers = {
        "flash_inner_loop": lambda: rr.flash_inner_loop(zs, ks),
        "Rachford_Rice_solution": lambda: rr.Rachford_Rice_solution(zs, ks),
        # This one returns (L/F, V/F, x, y).
        "Leibovici_Neoschil": lambda: rr.Rachford_Rice_solution_Leibovici_Neoschil(zs, ks)[1:],
        "LN2": lambda: rr.Rachford_Rice_solution_LN2(zs, ks),
    }
    for name, solve in solvers.items():
        try:
            beta, x, y = solve()
            if not all(math.isfinite(v) for v in [beta, *x, *y]):
                raise ValueError("non-finite result")
            out.append({"method": name, "beta": float(beta), "x": [float(v) for v in x],
                        "y": [float(v) for v in y]})
        except Exception as e:  # recorded, not fatal: a solver refusing a case is information
            out.append({"method": name, "error": f"{type(e).__name__}: {e}"[:200]})
    return out


def check_against_chemicals_mpmath(case):
    """chemicals has its own mpmath solver. It must agree with the reference,
    or one of them is wrong and no fixture is written."""
    try:
        liquid, beta, x, y = rr.Rachford_Rice_solution_mpmath(case["z"], case["k"])
    except Exception as e:
        case["chemicals"].append({"method": "mpmath", "error": f"{type(e).__name__}: {e}"[:200]})
        return
    case["chemicals"].append({"method": "mpmath", "beta": float(beta), "x": list(x), "y": list(y)})
    ref = case["reference"]
    scale = max(1.0, abs(ref["beta"]))
    assert abs(beta - ref["beta"]) <= 1e-15 * scale, (case["name"], beta, ref["beta"])
    for a, b in zip(x + y, ref["x"] + ref["y"]):
        assert abs(a - b) <= 1e-14 * abs(b), (case["name"], a, b)


def straddles(zs, ks):
    present = [ks[i] for i in range(len(zs)) if zs[i] > 0]
    return max(present) > 1 and min(present) < 1


def sign_at(zs, ks, beta):
    return sum(z * (k - 1) / (1 + beta * (k - 1)) for z, k in zip(zs, ks) if z > 0)


def normalised(raw):
    s = math.fsum(raw)
    return [v / s for v in raw]


def cases():
    rng = random.Random(SEED)
    out = []

    def add(name, zs, ks):
        assert straddles(zs, ks), name
        out.append({"name": name, "z": zs, "k": ks})

    # chemicals' own documented examples.
    add("chemicals docstring, ternary", [0.5, 0.3, 0.2], [1.685, 0.742, 0.532])
    add("chemicals docstring, trace liquid", [0.999999999999, 1e-12], [2.0, 1e-12])

    # Ordinary mixtures: K log-uniform over 1e-3 to 1e3.
    for n in (2, 3, 5, 10, 20):
        made = 0
        while made < 8:
            ks = [10 ** rng.uniform(-3, 3) for _ in range(n)]
            zs = normalised([rng.uniform(0.01, 1) for _ in range(n)])
            if straddles(zs, ks):
                made += 1
                add(f"random, {n} components, #{made}", zs, ks)

    # Extreme K: heavies down to 1e-30 and lights up to 1e25.
    for m in range(10):
        n = rng.choice((4, 6, 10))
        ks = [10 ** rng.uniform(-2, 2) for _ in range(n)]
        ks[0] = 10 ** rng.uniform(15, 25)
        ks[-1] = 10 ** rng.uniform(-30, -15)
        if m % 2:
            ks[-2] = 10 ** rng.uniform(-30, -15)
        zs = normalised([rng.uniform(0.01, 1) for _ in range(n)])
        add(f"extreme K, #{m + 1}", zs, ks)

    # A trace of liquid: a very heavy component present in a tiny amount,
    # in an otherwise light feed, so 1 - beta is about that amount.
    for m, amount in enumerate((1e-6, 1e-9, 1e-12, 1e-15)):
        n = 5
        ks = [10 ** rng.uniform(0.3, 2) for _ in range(n - 1)] + [10 ** rng.uniform(-25, -15)]
        zs = normalised([rng.uniform(0.1, 1) for _ in range(n - 1)])
        zs = [v * (1 - amount) for v in zs] + [amount]
        add(f"trace liquid {amount:g}", zs, ks)
    # And a trace of vapour: a very light component, tiny, in a heavy feed.
    for m, amount in enumerate((1e-6, 1e-9, 1e-12, 1e-15)):
        n = 5
        ks = [10 ** rng.uniform(15, 25)] + [10 ** rng.uniform(-2, -0.3) for _ in range(n - 1)]
        zs = normalised([rng.uniform(0.1, 1) for _ in range(n - 1)])
        zs = [amount] + [v * (1 - amount) for v in zs]
        add(f"trace vapour {amount:g}", zs, ks)

    # Near-critical: every K within 1e-8 to 1e-4 of 1.
    for m in range(8):
        n = rng.choice((2, 3, 5))
        while True:
            ks = [1 + rng.choice((-1, 1)) * 10 ** rng.uniform(-8, -4) for _ in range(n)]
            zs = normalised([rng.uniform(0.1, 1) for _ in range(n)])
            if straddles(zs, ks):
                break
        add(f"near-critical, #{m + 1}", zs, ks)

    # Negative flashes: K-values of a single-phase feed, so beta < 0 or > 1.
    low = high = 0
    while low < 8 or high < 8:
        n = rng.choice((2, 3, 5, 8))
        ks = [10 ** rng.uniform(-1.5, 1.5) for _ in range(n)]
        zs = normalised([rng.uniform(0.01, 1) for _ in range(n)])
        if not straddles(zs, ks):
            continue
        if sign_at(zs, ks, 0.0) < 0 and low < 8:
            low += 1
            add(f"negative flash, beta < 0, #{low}", zs, ks)
        elif sign_at(zs, ks, 1.0) > 0 and high < 8:
            high += 1
            add(f"negative flash, beta > 1, #{high}", zs, ks)

    # Absent components (z = 0, with K-values outside the others' range) and
    # a K of exactly 1.
    add("absent components", [0.0, 0.4, 0.35, 0.25, 0.0], [50.0, 2.5, 0.8, 0.3, 1e-3])
    add("a K of exactly 1", [0.3, 0.4, 0.3], [2.0, 1.0, 0.4])
    add("two components share K max", [0.25, 0.25, 0.5], [3.0, 3.0, 0.2])
    return out


def soave_m(kind, omega):
    if kind == "pr":
        return 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    return 0.480 + 1.574 * omega - 0.176 * omega * omega


def eos_at(kind, idx, zs, t, p):
    comps = [COMPONENTS[i] for i in idx]
    if len(comps) == 1:
        c = comps[0]
        cls = PR if kind == "pr" else SRK
        return cls(Tc=c["tc"], Pc=c["pc"], omega=c["omega"], T=t, P=p), [0.0]
    cls = PRMIX if kind == "pr" else SRKMIX
    kij = [[KIJ[i][j] for j in idx] for i in idx]
    eos = cls(Tcs=[c["tc"] for c in comps], Pcs=[c["pc"] for c in comps],
              omegas=[c["omega"] for c in comps], zs=zs, kijs=kij, T=t, P=p)
    return eos, [x for row in kij for x in row]


def roots(eos, full):
    """Each root thermo found, with its PIP and, if full, what M1 checks."""
    out = []
    for suffix in ("l", "g"):
        if not hasattr(eos, "V_" + suffix):
            continue
        pip = getattr(eos, "PIP_" + suffix)
        r = {"root": "liquid" if suffix == "l" else "vapor", "v": getattr(eos, "V_" + suffix), "pip": pip,
             "label": "liquid" if pip > 1 else "vapor"}
        if full:
            mixture = hasattr(eos, "zs") and len(eos.zs) > 1
            ln_phi = getattr(eos, "lnphis_" + suffix) if mixture else [getattr(eos, "lnphi_" + suffix)]
            r.update({"z": getattr(eos, "Z_" + suffix), "ln_phi": list(ln_phi),
                      "h_dep": getattr(eos, "H_dep_" + suffix), "s_dep": getattr(eos, "S_dep_" + suffix),
                      "cp_dep": getattr(eos, "Cp_dep_" + suffix)})
        out.append(r)
    return out


def pip_cases():
    systems = [([0], [1.0]), ([1], [1.0]), ([2], [1.0]), ([3], [1.0]),
               ([0, 1], [0.6, 0.4]), ([0, 1, 2], [0.2, 0.5, 0.3]), ([0, 3], [0.7, 0.3])]
    trs = [0.6, 0.8, 0.95, 0.99, 1.0, 1.01, 1.05, 1.2, 1.6, 2.5]
    prs = [0.01, 0.3, 0.9, 0.99, 1.0, 1.01, 1.1, 2.0, 5.0]
    out = []
    for kind in ("pr", "srk"):
        for idx, zs in systems:
            tpc = sum(z * COMPONENTS[i]["tc"] for z, i in zip(zs, idx))
            ppc = sum(z * COMPONENTS[i]["pc"] for z, i in zip(zs, idx))
            for tr in trs:
                t = tr * tpc
                # No state past a Soave zero: there reflux and thermo differ on purpose.
                for i in idx:
                    c = COMPONENTS[i]
                    t0 = c["tc"] * (1 + 1 / soave_m(kind, c["omega"])) ** 2
                    assert t < 0.9 * t0, (kind, c["name"], t, t0)
                for pr in prs:
                    if len(idx) == 1 and tr == 1.0 and pr == 1.0:
                        # A pure component's critical point: dP/dV and d2P/dV2 are both
                        # zero, so the PIP is 0/0 there, and any label is arbitrary.
                        continue
                    eos, kij = eos_at(kind, idx, zs, t, pr * ppc)
                    out.append({"eos": kind, "components": idx, "kij": kij, "n": zs, "t": t, "p": pr * ppc,
                                "roots": roots(eos, False)})
        # Liquids squeezed to within 0.02% to 2.5% of the covolume (V/b from
        # 1.0002 to 1.025). Not physical pressures, but they test that the
        # PIP's volume step scales with V - b, not V.
        for idx in ([1], [3]):
            t = 0.6 * COMPONENTS[idx[0]]["tc"]
            for p in (1e9, 1e10, 1e11):
                eos, kij = eos_at(kind, idx, [1.0], t, p)
                out.append({"eos": kind, "components": idx, "kij": kij, "n": [1.0], "t": t, "p": p,
                            "roots": roots(eos, False)})
    return out


def soave_zero_cases():
    out = []
    co2 = COMPONENTS[2]
    for kind in ("pr", "srk"):
        t0 = co2["tc"] * (1 + 1 / soave_m(kind, co2["omega"])) ** 2
        for ratio in (0.9, 0.99, 0.999, 1.0, 1.001, 1.01, 1.1, 1.5):
            for p in (1e5, 1e6, 1e7):
                eos, _ = eos_at(kind, [2], [1.0], ratio * t0, p)
                out.append({"eos": kind, "components": [2], "kij": [0.0], "n": [1.0], "t": ratio * t0, "p": p,
                            "t0": t0, "roots": roots(eos, True)})
    return out


def pip_fixture():
    return {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_fixtures.py",
            "python": platform.python_version(),
            "thermo": thermo.__version__,
        },
        "note": "thermo's PR and SRK with these constants and kij. Each root: molar volume v (m3/mol), thermo's "
        "analytic PIP, and its label (liquid when PIP > 1). soave_zero is pure CO2 around T0 = Tc(1 + 1/m)^2, "
        "where Soave's alpha is zero, with Z, ln phi and the departures (J/mol, J/(mol K)) as in M1.",
        "components": COMPONENTS,
        "cases": pip_cases(),
        "soave_zero": soave_zero_cases(),
    }


def write_pip():
    pip = pip_fixture()
    OUT_PIP.write_text(json.dumps(pip, indent=1) + "\n", encoding="utf-8", newline="\n")
    labels = [r["label"] for c in pip["cases"] for r in c["roots"]]
    near = sum(1 for c in pip["cases"] for r in c["roots"] if abs(r["pip"] - 1) < 1e-6)
    print(f"wrote {OUT_PIP.relative_to(ROOT)}: {len(pip['cases'])} states, {len(labels)} roots "
          f"({labels.count('liquid')} liquid, {labels.count('vapor')} vapour, {near} within 1e-6 of PIP = 1); "
          f"{len(pip['soave_zero'])} CO2 states around the Soave zero")


# Stability (m3_stability.json). Smith, Van Ness & Abbott, Table B.1, as above.
STAB_COMPONENTS = [
    {"name": "methane", "tc": 190.6, "pc": 45.99e5, "omega": 0.012, "mw": 0.016043},
    {"name": "ethane", "tc": 305.3, "pc": 48.72e5, "omega": 0.100, "mw": 0.030070},
    {"name": "propane", "tc": 369.8, "pc": 42.48e5, "omega": 0.152, "mw": 0.044097},
    {"name": "n-butane", "tc": 425.1, "pc": 37.96e5, "omega": 0.200, "mw": 0.058123},
    {"name": "n-pentane", "tc": 469.7, "pc": 33.70e5, "omega": 0.252, "mw": 0.072150},
    {"name": "n-hexane", "tc": 507.6, "pc": 30.25e5, "omega": 0.301, "mw": 0.086177},
    {"name": "n-heptane", "tc": 540.2, "pc": 27.40e5, "omega": 0.350, "mw": 0.100204},
    {"name": "n-decane", "tc": 617.7, "pc": 21.10e5, "omega": 0.492, "mw": 0.142285},
    {"name": "carbon dioxide", "tc": 304.2, "pc": 73.83e5, "omega": 0.224, "mw": 0.044010},
    {"name": "nitrogen", "tc": 126.2, "pc": 34.00e5, "omega": 0.038, "mw": 0.028014},
]
C1, C2, C3, NC4, NC5, NC6, NC7, NC10, CO2, N2 = range(10)


def stab_kij(i, j):
    """Illustrative kij: CO2 and N2 with hydrocarbons, and methane with the
    heavier ones. The tests check the algorithm, not the data."""
    pair = {i, j}
    if i == j:
        return 0.0
    if pair == {CO2, N2}:
        return -0.02
    if CO2 in pair:
        return 0.12
    if N2 in pair:
        return 0.03
    if C1 in pair:
        other = (pair - {C1}).pop()
        return {NC4: 0.02, NC10: 0.04}.get(other, 0.0)
    return 0.0


def stab_eos(kind, idx, zs, t, p):
    comps = [STAB_COMPONENTS[i] for i in idx]
    cls = PRMIX if kind == "pr" else SRKMIX
    kij = [[stab_kij(i, j) for j in idx] for i in idx]
    return cls(Tcs=[c["tc"] for c in comps], Pcs=[c["pc"] for c in comps], omegas=[c["omega"] for c in comps],
               zs=zs, kijs=kij, T=t, P=p), kij


def lnphis_stable(eos):
    """ln phi at the root of lower Gibbs energy, as stability and reflux use."""
    if hasattr(eos, "V_l") and hasattr(eos, "V_g"):
        return list(eos.lnphis_l) if eos.G_dep_l < eos.G_dep_g else list(eos.lnphis_g)
    return list(eos.lnphis_l) if hasattr(eos, "V_l") else list(eos.lnphis_g)


def tpd(kind, idx, t, p, d, w):
    """The normalised tangent-plane distance of composition w against a feed
    whose d_i = ln z_i + ln phi_i(z)."""
    eos, _ = stab_eos(kind, idx, w, t, p)
    lp = lnphis_stable(eos)
    return sum(wi * (math.log(wi) + lpi - di) for wi, lpi, di in zip(w, lp, d) if wi > 0)


def brute_binary(kind, idx, t, p, z):
    """The global minimum of tm over the whole binary composition range: a
    scan of 4,002 compositions, linear and logistic in w1, then each local
    minimum refined by golden-section search. It depends on no iteration, so
    it settles stable or unstable outright. The trivial minimum at w = z
    (tm = 0) is left out."""
    eos, _ = stab_eos(kind, idx, z, t, p)
    d = [math.log(zi) + lpi for zi, lpi in zip(z, lnphis_stable(eos))]

    def tm_of(s):  # s is the logistic coordinate: w1 = 1/(1 + exp(-s))
        w1 = 1.0 / (1.0 + math.exp(-s))
        return tpd(kind, idx, t, p, d, [w1, 1.0 - w1])

    grid = sorted({math.log(w / (1 - w)) for w in [(k + 0.5) / 2001 for k in range(2001)]}
                  | {-40.0 + 80.0 * k / 2000 for k in range(2001)})
    vals = [tm_of(s) for s in grid]
    s_feed = math.log(z[0] / z[1])
    best = (0.0, None)
    for k in range(1, len(grid) - 1):
        if vals[k] <= vals[k - 1] and vals[k] <= vals[k + 1]:
            a, b = grid[k - 1], grid[k + 1]
            g = (math.sqrt(5) - 1) / 2
            c, e = b - g * (b - a), a + g * (b - a)
            fc, fe = tm_of(c), tm_of(e)
            for _ in range(120):
                if fc < fe:
                    b, e, fe = e, c, fc
                    c = b - g * (b - a)
                    fc = tm_of(c)
                else:
                    a, c, fc = c, e, fe
                    e = a + g * (b - a)
                    fe = tm_of(e)
                if b - a < 1e-13 * max(1.0, abs(a)):
                    break
            s = 0.5 * (a + b)
            if abs(s - s_feed) < 1e-4:
                continue  # the trivial minimum, w = z
            val = tm_of(s)
            if val < best[0]:
                w1 = 1.0 / (1.0 + math.exp(-s))
                best = (val, [w1, 1.0 - w1])
    return {"tm_min": best[0], "w": best[1]}


def thermo_flash(kind, idx, t, p, z):
    """thermo's FlashVL at the tolerance step 0's spike found it needs."""
    comps = [STAB_COMPONENTS[i] for i in idx]
    kij = [[stab_kij(i, j) for j in idx] for i in idx]
    tcs, pcs, ws = [c["tc"] for c in comps], [c["pc"] for c in comps], [c["omega"] for c in comps]
    constants = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=ws, MWs=[c["mw"] * 1000 for c in comps],
                                         CASs=[None] * len(comps))
    # A PT flash doesn't use Cp; FlashVL only needs the objects.
    cps = [HeatCapacityGas(poly_fit=(1.0, 2000.0, [0.0, 30.0])) for _ in comps]
    correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=cps, skip_missing=True)
    eos_cls = PRMIX if kind == "pr" else SRKMIX
    args = dict(Tcs=tcs, Pcs=pcs, omegas=ws, kijs=kij)
    flasher = FlashVL(constants, correlations, liquid=CEOSLiquid(eos_cls, args, HeatCapacityGases=cps),
                      gas=CEOSGas(eos_cls, args, HeatCapacityGases=cps))
    # At 1e-30 its substitution can end in a rounding cycle (OscillationError),
    # so step the tolerance down until one converges; the spike found 1e-25
    # still leaves the fugacities equal to 1e-13.
    res, tol = None, None
    for tol in (1e-30, 1e-25, 1e-20, 1e-13):
        flasher.PT_SS_TOL = tol
        try:
            res = flasher.flash(T=t, P=p, zs=z)
            break
        except Exception as e:
            err = f"{type(e).__name__}: {e}"[:200]
    if res is None:
        return {"phases": 0, "error": err}
    if res.phase_count == 1:
        return {"phases": 1, "tol": tol}
    phases = [res.gas] + list(res.liquids) if res.gas is not None else list(res.liquids)
    return {"phases": res.phase_count, "tol": tol, "betas": list(res.betas),
            "compositions": [list(ph.zs) for ph in phases]}


def stability_cases():
    out = []

    def add(kind, idx, z, t, p, brute):
        eos, kij = stab_eos(kind, idx, z, t, p)
        case = {"eos": kind, "components": idx, "kij": [x for row in kij for x in row], "z": z, "t": t, "p": p,
                "thermo": thermo_flash(kind, idx, t, p, z)}
        d = [math.log(zi) + lpi for zi, lpi in zip(z, lnphis_stable(eos))]
        if case["thermo"]["phases"] > 1:
            # A split with lower G is a witness: at least one of its phases
            # has tm < 0 against the feed's tangent plane.
            case["witness_tm"] = min(tpd(kind, idx, t, p, d, w) for w in case["thermo"]["compositions"])
        if brute:
            case["brute"] = brute_binary(kind, idx, t, p, z)
            case["stable"] = case["brute"]["tm_min"] > -1e-10
        elif case["thermo"]["phases"] == 0:
            print("  left out, no reference:", kind, idx, z, t, p, case["thermo"]["error"])
            return
        else:
            case["stable"] = case["thermo"]["phases"] == 1
        out.append(case)

    for kind in ("pr", "srk"):
        for t in (200.0, 250.0, 300.0, 350.0, 380.0):
            for p in (5e5, 20e5, 50e5, 80e5, 110e5, 130e5):
                for z1 in (0.1, 0.3, 0.5, 0.7, 0.9):
                    add(kind, [C1, NC4], [z1, 1 - z1], t, p, True)
    for t in (300.0, 400.0, 500.0):
        for p in (20e5, 100e5, 200e5, 300e5):
            for z1 in (0.3, 0.6, 0.9):
                add("pr", [C1, NC10], [z1, 1 - z1], t, p, True)
    for t in (280.0, 320.0, 360.0, 400.0):
        for p in (20e5, 50e5, 80e5, 100e5):
            for z1 in (0.2, 0.5, 0.8):
                add("pr", [CO2, NC4], [z1, 1 - z1], t, p, True)
    # Near the methane/n-butane critical locus, where the two phases become
    # alike and substitution crawls.
    for t in (340.0, 350.0, 360.0):
        for k in range(9):
            add("pr", [C1, NC4], [0.6, 0.4], t, 116e5 + k * 2e5, True)
    # Closer still: where substitution alone needs 100 to 465 steps to reach
    # the trivial solution, and the second-order finish has to do the work.
    # (z1, T, the first P in units of 0.1 bar); 11 pressures 0.1 bar apart
    # from each. The pressures are built from integers so that they're the
    # exact doubles 9.23e6 and so on: a few of these states are sensitive
    # enough that an ulp decides whether the hard path is taken.
    for z1, t, p0 in ((0.5, 380.0, 915), (0.6, 340.0, 1245), (0.7, 330.0, 1315)):
        for k in range(11):
            add("pr", [C1, NC4], [z1, 1 - z1], t, (p0 + k) * 1e4, True)
    for z1, t, p in ((0.6, 370.0, 963), (0.8, 340.0, 938)):
        add("pr", [C1, NC4], [z1, 1 - z1], t, p * 1e4, True)
    for t in (250.0, 300.0, 350.0):
        for p in (20e5, 50e5, 80e5):
            for z in ([0.5, 0.3, 0.2], [0.2, 0.5, 0.3], [0.3, 0.3, 0.4]):
                add("pr", [C1, NC4, CO2], z, t, p, False)
    condensate = [0.70, 0.08, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01]
    for t in (250.0, 300.0, 350.0, 400.0, 450.0):
        for p in (10e5, 50e5, 100e5, 150e5, 200e5, 250e5, 300e5, 350e5):
            add("pr", list(range(10)), condensate, t, p, False)
    return out


def stability_fixture():
    return {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_fixtures.py",
            "python": platform.python_version(),
            "thermo": thermo.__version__,
        },
        "note": "stable is the reference verdict: the brute-force scan's where there is one (tm_min > -1e-10), "
        "otherwise thermo's FlashVL (one phase). thermo holds FlashVL's result at PT_SS_TOL 1e-30; witness_tm "
        "is the lowest tm of its phases against the feed. tm is normalised: sum w_i(ln w_i + ln phi_i(w) - d_i).",
        "components": STAB_COMPONENTS,
        "cases": stability_cases(),
    }


def write_stability():
    data = stability_fixture()
    OUT_STAB.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    cases = data["cases"]
    unstable = sum(1 for c in cases if not c["stable"])
    disagree = [c for c in cases if "brute" in c and c["thermo"]["phases"] > 0
                and (c["thermo"]["phases"] == 1) != c["stable"]]
    failed = sum(1 for c in cases if c["thermo"]["phases"] == 0)
    weak = [c for c in cases if "witness_tm" in c and c["witness_tm"] >= 0]
    print(f"wrote {OUT_STAB.relative_to(ROOT)}: {len(cases)} states, {unstable} unstable; "
          f"thermo failed on {failed}, disagrees with the brute-force scan on {len(disagree)}; "
          f"{len(weak)} thermo splits without a witness")
    for c in disagree:
        print("  disagreement:", c["eos"], c["components"], c["z"], c["t"], c["p"], c["brute"]["tm_min"], c["thermo"])


# Flash references (m3_thermo.json): the stability components plus
# isobutane (Smith, Van Ness & Abbott, Table B.1), each with an ideal-gas
# Cp/R = A + B T + C T^2 after the same book's Table C.1. N2 and CO2 drop
# their D/T^2 term, and n-decane, which the table lacks, is n-heptane plus
# three CH2 increments (n-heptane minus n-hexane). As with the kij, the tests
# check the algorithm, not the data. The fixture stores Cp itself, in
# J/(mol K), so both programs read the same floats whichever R made them.
R_GAS = 8.314462618
CP_OVER_R = {
    "methane": (1.702, 9.081e-3, -2.164e-6),
    "ethane": (1.131, 19.225e-3, -5.561e-6),
    "propane": (1.213, 28.785e-3, -8.824e-6),
    "n-butane": (1.935, 36.915e-3, -11.402e-6),
    "n-pentane": (2.464, 45.351e-3, -14.111e-6),
    "n-hexane": (3.025, 53.722e-3, -16.791e-6),
    "n-heptane": (3.570, 62.127e-3, -19.486e-6),
    "n-decane": (5.205, 87.342e-3, -27.571e-6),
    "carbon dioxide": (5.457, 1.045e-3, 0.0),
    "nitrogen": (3.280, 0.593e-3, 0.0),
    "isobutane": (1.677, 37.853e-3, -11.945e-6),
}
FLASH_COMPONENTS = [
    dict(c, cp=[R_GAS * a for a in CP_OVER_R[c["name"]]])
    for c in STAB_COMPONENTS + [{"name": "isobutane", "tc": 408.1, "pc": 36.48e5, "omega": 0.181, "mw": 0.058123}]
]
IC4 = 10
CONDENSATE = [0.70, 0.08, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01]


def flash_kij(i, j):
    """stab_kij, with isobutane taking n-butane's values, and 0 with n-butane."""
    a, b = (NC4 if i == IC4 else i), (NC4 if j == IC4 else j)
    return 0.0 if a == b else stab_kij(a, b)


def flasher_for(kind, idx):
    """thermo's flasher for these components: FlashPureVLS for one, FlashVL
    for a mixture. Returns it with the kij matrix, flattened."""
    comps = [FLASH_COMPONENTS[i] for i in idx]
    kij = [[flash_kij(i, j) for j in idx] for i in idx]
    # poly_fit takes the coefficients highest power first. 50 to 3000 K
    # covers every state, so thermo never extrapolates.
    cps = [HeatCapacityGas(poly_fit=(50.0, 3000.0, list(reversed(c["cp"])))) for c in comps]
    return make_flasher(PRMIX if kind == "pr" else SRKMIX, comps, kij, cps), [x for row in kij for x in row]


def make_flasher(eos_cls, comps, kij, cps=None):
    """thermo's flasher with this equation of state, constants (Tc, Pc,
    omega and M in kg/mol) and kij matrix: FlashPureVLS for one component,
    FlashVL for a mixture. Without cps it gets a constant Cp, which PT,
    bubble and dew flashes never use."""
    tcs, pcs, ws = [c["tc"] for c in comps], [c["pc"] for c in comps], [c["omega"] for c in comps]
    constants = ChemicalConstantsPackage(Tcs=tcs, Pcs=pcs, omegas=ws, MWs=[c["mw"] * 1000 for c in comps],
                                         CASs=[None] * len(comps))
    if cps is None:
        cps = [HeatCapacityGas(poly_fit=(1.0, 3000.0, [0.0, 30.0])) for _ in comps]
    correlations = PropertyCorrelationsPackage(constants, HeatCapacityGases=cps, skip_missing=True)
    args = dict(Tcs=tcs, Pcs=pcs, omegas=ws, kijs=kij)
    liquid = CEOSLiquid(eos_cls, args, HeatCapacityGases=cps)
    gas = CEOSGas(eos_cls, args, HeatCapacityGases=cps)
    if len(comps) == 1:
        f = FlashPureVLS(constants, correlations, gas=gas, liquids=[liquid], solids=[])
        f.PVF_xtol = f.TVF_xtol = f.TPV_HSGUA_xtol = 1e-14
    else:
        f = FlashVL(constants, correlations, liquid=liquid, gas=gas)
        f.DEW_BUBBLE_QUASI_NEWTON_XTOL = 1e-14
        f.DEW_BUBBLE_NEWTON_XTOL = 1e-14
        # Every other flash type runs PT flashes inside it, at this tolerance.
        f.PT_SS_TOL = 1e-30
    return f


def thermo_pt(f, t, p, z):
    """A PT flash, stepping FlashVL's tolerance down from 1e-30 as in
    thermo_flash. Returns the result and the tolerance, or the error. The
    flasher is left at 1e-30, so that what runs next doesn't depend on the
    order of the cases."""
    if isinstance(f, FlashPureVLS):
        return f.flash(T=t, P=p, zs=z), None
    err = None
    try:
        for tol in (1e-30, 1e-25, 1e-20, 1e-13):
            f.PT_SS_TOL = tol
            try:
                return f.flash(T=t, P=p, zs=z), tol
            except Exception as e:
                err = f"{type(e).__name__}: {e}"[:200]
    finally:
        f.PT_SS_TOL = 1e-30
    raise RuntimeError(err)


def by_density(res):
    """A result's phases, the least dense by mass first (reflux's D6), each
    with its fraction. Molar volume can order them the other way: for
    methane with n-decane at 300 K and 200 bar, the methane-rich phase has
    the smaller molar volume but a third of the mass density."""
    mws = res.constants.MWs
    pairs = list(zip(res.betas, res.phases))
    return sorted(pairs, key=lambda bp: sum(x * mw for x, mw in zip(bp[1].zs, mws)) / bp[1].V())


def fugacity_residual(a, b):
    """max |ln f_i(a) - ln f_i(b)| over the components in both phases."""
    return max(abs(math.log(x) + lx - math.log(y) - ly)
               for x, lx, y, ly in zip(a.zs, a.lnphis(), b.zs, b.lnphis()) if x > 0 and y > 0)


def describe(res):
    """A thermo result as the fixture stores it: each phase (fraction,
    composition, molar V, H and S, and thermo's label for it), the
    mixture's molar V, H, S and U, and the fugacity residual of a split."""
    pairs = by_density(res)
    out = {
        "phases": [{"fraction": beta, "z": list(ph.zs), "v": ph.V(), "h": ph.H(), "s": ph.S(),
                    "label": "vapor" if ph is res.gas else "liquid"} for beta, ph in pairs],
        "v": res.V(), "h": res.H(), "s": res.S(), "u": res.U(),
    }
    if len(pairs) == 2:
        out["residual"] = fugacity_residual(pairs[0][1], pairs[1][1])
    return out


def landing(f, z, **spec):
    """Where thermo's own flash lands for this specification, or its error."""
    try:
        res = f.flash(zs=z, **spec)
    except Exception as e:
        return {"error": f"{type(e).__name__}: {e}"[:120]}
    pairs = by_density(res)
    return {"t": res.T, "p": res.P, "phases": len(pairs), "beta": pairs[0][0]}


def flash_pt_cases():
    out = []

    def add(kind, idx, z, t, p):
        f, kij = flasher_for(kind, idx)
        res, tol = thermo_pt(f, t, p, z)
        state = describe(res)
        case = {"eos": kind, "components": idx, "kij": kij, "z": z, "t": t, "p": p, "tol": tol, **state}
        case["thermo"] = {"ph": landing(f, z, P=p, H=state["h"]), "ps": landing(f, z, P=p, S=state["s"]),
                          "tv": landing(f, z, T=t, V=state["v"])}
        out.append(case)

    for t in (250.0, 300.0, 350.0):
        for p in (5e5, 30e5, 60e5, 100e5):
            for z1 in (0.3, 0.6):
                add("pr", [C1, NC4], [z1, 1 - z1], t, p)
    for t, p in ((300.0, 30e5), (300.0, 60e5), (350.0, 30e5)):
        add("srk", [C1, NC4], [0.5, 0.5], t, p)
    for t in (300.0, 360.0):
        for p in (30e5, 60e5):
            add("pr", [CO2, NC4], [0.5, 0.5], t, p)
    for t in (300.0, 450.0):
        for p in (50e5, 200e5):
            add("pr", [C1, NC10], [0.6, 0.4], t, p)
    for t, p in ((120.0, 10e5), (150.0, 20e5), (150.0, 40e5)):
        add("pr", [N2, C1], [0.3, 0.7], t, p)
    for t, p in ((250.0, 20e5), (300.0, 50e5), (350.0, 80e5)):
        add("pr", [C1, NC4, CO2], [0.5, 0.3, 0.2], t, p)
    for t in (300.0, 400.0):
        for p in (50e5, 100e5, 200e5, 350e5):
            add("pr", list(range(10)), CONDENSATE, t, p)
    # Near the methane/n-butane critical locus, from the stability fixture's
    # states 0.1 bar apart (exact doubles, as there).
    for k in (0, 3, 6, 9):
        add("pr", [C1, NC4], [0.6, 0.4], 340.0, (1245 + k) * 1e4)
    # Pure components: liquid, vapour and supercritical states.
    for t in (250.0, 300.0, 350.0, 450.0):
        for p in (1e5, 5e5, 20e5, 50e5):
            add("pr", [NC4], [1.0], t, p)
    for t in (280.0, 310.0):
        for p in (50e5, 80e5):
            add("pr", [CO2], [1.0], t, p)
    # Narrow-boiling: n-butane with isobutane, whose two-phase range at
    # 5 bar is under 4 K wide (and under 0.4 K at 5% isobutane).
    for z1 in (0.5, 0.95):
        f, _ = flasher_for("pr", [NC4, IC4])
        bubble = f.flash(P=5e5, VF=0.0, zs=[z1, 1 - z1]).T
        dew = f.flash(P=5e5, VF=1.0, zs=[z1, 1 - z1]).T
        for frac in (-0.5, 0.25, 0.5, 0.75, 1.5):
            add("pr", [NC4, IC4], [z1, 1 - z1], bubble + frac * (dew - bubble), 5e5)
    return out


def saturation_cases():
    """Pure components on their saturation curve, at a given P or T and
    vapour fraction, from FlashPureVLS."""
    out = []
    for idx, ps, ts in (([NC4], (1e5, 5e5, 20e5), (280.0, 350.0, 410.0)),
                        ([CO2], (10e5, 40e5, 70e5), (230.0, 280.0, 300.0))):
        f, kij = flasher_for("pr", idx)
        for spec, values in (("p", ps), ("t", ts)):
            for value in values:
                for vf in (0.0, 0.3, 1.0):
                    res = f.flash(zs=[1.0], VF=vf, **{spec.upper(): value})
                    out.append({"eos": "pr", "components": idx, "kij": kij, "z": [1.0], "spec": spec, "vf": vf,
                                "t": res.T, "p": res.P, **describe(res)})
    return out


def polish_boundary(f, z, vf, spec, t, p, w):
    """Newton from thermo's bubble or dew point, on thermo's own ln phi.

    The unknowns are u_i = ln(w_i / z_i), for the incipient composition w,
    and ln T (P given) or ln P (T given). The equations are
    u_i + ln phi_i(w) - ln phi_i(z) = 0, each phase on thermo's root for it
    (the feed's own, the liquid for a bubble point), and sum(w) = 1. The
    Jacobian is by central differences, so each step gains about seven
    digits; three steps take thermo's 1e-7 to round-off."""
    import numpy as np

    feed_phase, new_phase = (f.liquid, f.gas) if vf == 0.0 else (f.gas, f.liquid)
    n = len(z)

    def tp(x):
        return (math.exp(x), p) if spec == "p" else (t, math.exp(x))

    def g(v):
        tt, pp = tp(v[-1])
        ws = [z[i] * math.exp(v[i]) for i in range(n)]
        total = math.fsum(ws)
        lf = feed_phase.to(T=tt, P=pp, zs=z).lnphis()
        ln = new_phase.to(T=tt, P=pp, zs=[x / total for x in ws]).lnphis()
        return np.array([v[i] + ln[i] - lf[i] for i in range(n)] + [total - 1.0])

    v = np.array([math.log(w[i] / z[i]) for i in range(n)] + [math.log(t if spec == "p" else p)])
    for _ in range(6):
        r = g(v)
        if max(abs(r)) < 1e-15:
            break
        jac = np.empty((n + 1, n + 1))
        for k in range(n + 1):
            h = 1e-6 * max(1.0, abs(v[k]))
            up, down = v.copy(), v.copy()
            up[k] += h
            down[k] -= h
            jac[:, k] = (g(up) - g(down)) / (2 * h)
        v = v - np.linalg.solve(jac, r)
    tt, pp = tp(v[-1])
    ws = [z[i] * math.exp(v[i]) for i in range(n)]
    total = math.fsum(ws)
    return tt, pp, [x / total for x in ws]


def boundary_problem(f, z, spec, t, p, w, residual, shift):
    """Why a polished bubble or dew point isn't one, or None. The polish
    must have converged without moving far, and PT flashes 1e-7 to either
    side (relatively, in the unknown T or P) must show a phase boundary:
    one side a single phase, the other a split in which the phase nearer
    the incipient composition is under 0.1% of the feed. (The step is
    small because n-butane with isobutane splits over only 1.2 K.)"""
    if not residual < 1e-12 or not abs(shift) < 1e-6:
        return f"the polish didn't converge (residual {residual:.2g}, shift {shift:.2g})"
    sides = []
    for d in (1e-7, -1e-7):
        tt, pp = (t * (1 + d), p) if spec == "p" else (t, p * (1 + d))
        res, _ = thermo_pt(f, tt, pp, z)
        if res.phase_count == 1:
            sides.append(0.0)
            continue
        near = min(zip(res.betas, res.phases), key=lambda bp: sum(abs(a - b) for a, b in zip(bp[1].zs, w)))
        sides.append(near[0])
    if sorted(sides)[0] != 0.0 or not 0.0 < sorted(sides)[1] < 1e-3:
        return f"PT flashes either side give incipient-phase fractions {sides}"
    return None


def boundary_cases():
    """Bubble and dew points of mixtures, at a given P or T. incipient is
    the composition of the phase that is just appearing. thermo's own
    answer stops at a fugacity residual up to 2e-7 (thermo_residual), so
    each is polished by polish_boundary; shift is how far that moved the
    unknown T or P, relatively. An answer that boundary_problem finds isn't
    a bubble or dew point is recorded as rejected, with the reason."""
    out = []
    systems = [
        ("pr", [C1, NC4], [0.4, 0.6], (10e5, 30e5, 60e5), (250.0, 300.0)),
        ("srk", [C1, NC4], [0.4, 0.6], (30e5,), (300.0,)),
        ("pr", [CO2, NC4], [0.5, 0.5], (30e5,), (300.0,)),
        ("pr", [C1, NC10], [0.6, 0.4], (50e5,), (450.0,)),
        ("pr", [C1, NC4, CO2], [0.5, 0.3, 0.2], (30e5,), (250.0,)),
        ("pr", [NC4, IC4], [0.5, 0.5], (1e5, 5e5), (300.0,)),
        ("pr", [NC4, IC4], [0.95, 0.05], (5e5,), (300.0,)),
        ("pr", list(range(10)), CONDENSATE, (50e5,), (300.0, 400.0)),
    ]
    for kind, idx, z, ps, ts in systems:
        f, kij = flasher_for(kind, idx)
        for spec, values in (("p", ps), ("t", ts)):
            for value in values:
                for vf in (0.0, 1.0):
                    case = {"eos": kind, "components": idx, "kij": kij, "z": z, "spec": spec, "vf": vf}
                    try:
                        res = f.flash(zs=z, VF=vf, **{spec.upper(): value})
                    except Exception as e:
                        case[spec] = value
                        case["error"] = f"{type(e).__name__}: {e}"[:120]
                        out.append(case)
                        continue
                    feed, new = (res.liquid0, res.gas) if vf == 0.0 else (res.gas, res.liquid0)
                    raw = fugacity_residual(feed, new)
                    t, p, w = polish_boundary(f, z, vf, spec, res.T, res.P, list(new.zs))
                    feed_phase, new_phase = (f.liquid, f.gas) if vf == 0.0 else (f.gas, f.liquid)
                    feed, new = feed_phase.to(T=t, P=p, zs=z), new_phase.to(T=t, P=p, zs=w)
                    shift = (t - res.T) / res.T if spec == "p" else (p - res.P) / res.P
                    residual = fugacity_residual(feed, new)
                    problem = boundary_problem(f, z, spec, t, p, w, residual, shift)
                    if problem:
                        case.update({"t": res.T, "p": res.P, "thermo_residual": raw,
                                     "rejected": "thermo's answer isn't a boundary: " + problem})
                        out.append(case)
                        continue
                    case.update({"t": t, "p": p, "incipient": w, "residual": residual,
                                 "thermo_residual": raw, "shift": shift,
                                 "v": feed.V(), "h": feed.H(), "s": feed.S()})
                    out.append(case)
    return out


def write_thermo():
    pt = flash_pt_cases()
    saturation = saturation_cases()
    boundaries = boundary_cases()
    data = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_fixtures.py",
            "python": platform.python_version(),
            "thermo": thermo.__version__,
        },
        "note": "thermo's PR (SRK where eos says so) with these constants, kij and ideal-gas Cp = cp[0] + "
        "cp[1] T + cp[2] T^2 in J/(mol K). Molar v (m3/mol), h (J/mol), s (J/(mol K)) and u (J/mol); h is "
        "relative to the ideal gas at t_ref, and s to the ideal gas at t_ref and p_ref_s. Phases are sorted "
        "by mass density, the vapour first; label is thermo's (by PIP), residual is max |ln fV - ln fL|. tol is "
        "FlashVL's PT_SS_TOL. thermo holds where thermo's own PH, PS and TV flashes land from h, s and v, or "
        "the error they raise.",
        "t_ref": Phase.T_REF_IG,
        "p_ref_s": Phase.P_REF_IG,
        "components": FLASH_COMPONENTS,
        "pt": pt,
        "saturation": saturation,
        "boundaries": boundaries,
    }
    OUT_THERMO.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8", newline="\n")
    splits = [c for c in pt if len(c["phases"]) == 2]
    good = [c for c in boundaries if "incipient" in c]
    print(f"wrote {OUT_THERMO.relative_to(ROOT)}: {len(pt)} PT states ({len(splits)} split, worst residual "
          f"{max(c['residual'] for c in splits):.2g}), {len(saturation)} saturation states (worst residual "
          f"{max(c.get('residual', 0.0) for c in saturation):.2g}), {len(boundaries)} bubble and dew points "
          f"({sum(1 for c in boundaries if 'error' in c)} refused by thermo, "
          f"{sum(1 for c in boundaries if 'rejected' in c)} rejected; worst residual "
          f"{max(c['residual'] for c in good):.2g} after polishing thermo's "
          f"{max(c['thermo_residual'] for c in good):.2g}, which moved T or P by up to "
          f"{max(abs(c['shift']) for c in good):.2g}; incipient phase at least "
          f"{min(max(abs(a - b) for a, b in zip(c['incipient'], c['z'])) for c in good):.2g} from the feed)")
    for c in boundaries:
        if "rejected" in c or "error" in c:
            print(f"  {c['components']} {c['spec']} {c['vf']}: {c.get('rejected') or c.get('error')}")
    for name in ("ph", "ps", "tv"):
        landed = [c for c in pt if "error" not in c["thermo"][name]]
        worst = max((max(abs(c["thermo"][name]["t"] - c["t"]) / c["t"], abs(c["thermo"][name]["p"] - c["p"]) / c["p"])
                     for c in landed), default=0.0)
        print(f"  thermo's {name.upper()} flash: {len(landed)} of {len(pt)} ran, worst |dT/T|, |dP/P| {worst:.2g}")


def write_rachford_rice():
    data = cases()
    for case in data:
        case["reference"] = reference(case["z"], case["k"])
        case["chemicals"] = chemicals_solvers(case["z"], case["k"])
        check_against_chemicals_mpmath(case)
    fixture = {
        "source": {
            "generated": datetime.date.today().isoformat(),
            "generator": "tools/gen_m3_fixtures.py",
            "python": platform.python_version(),
            "chemicals": chemicals.__version__,
            "mpmath": mp.__version__,
            "seed": SEED,
        },
        "note": "reference is the exact root for these floats (bisection at 100 digits), rounded once; "
        "liquid is 1 - beta at full precision. chemicals lists each of chemicals' solvers' result, or the "
        "error it raised; its mpmath solver is checked against the reference before writing.",
        "cases": data,
    }
    OUT_RR.write_text(json.dumps(fixture, indent=1) + "\n", encoding="utf-8", newline="\n")
    refused = {}
    for case in data:
        for r in case["chemicals"]:
            if "error" in r:
                refused[r["method"]] = refused.get(r["method"], 0) + 1
    print(f"wrote {OUT_RR.relative_to(ROOT)}: {len(data)} cases; chemicals refused: {refused or 'none'}")


def main():
    writers = {"stability": write_stability, "pip": write_pip, "rachford_rice": write_rachford_rice,
               "thermo": write_thermo}
    which = sys.argv[1:] or list(writers)
    unknown = [w for w in which if w not in writers]
    if unknown:
        sys.exit(f"no fixture called {', '.join(unknown)}; the fixtures are {', '.join(writers)}")
    for name in which:
        writers[name]()


if __name__ == "__main__":
    main()
