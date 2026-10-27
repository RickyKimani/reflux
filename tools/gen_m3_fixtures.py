# /// script
# requires-python = ">=3.11"
# dependencies = ["chemicals==1.5.2", "mpmath==1.3.0", "thermo==0.6.1"]
# ///
"""Generate the M3 reference fixtures for reflux.

Run from the repository root:

    uv run tools/gen_m3_fixtures.py

It writes to testdata/:

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
"""

import datetime
import json
import math
import pathlib
import platform
import random

import chemicals
import chemicals.rachford_rice as rr
import mpmath as mp
import thermo
from thermo import PR, PRMIX, SRK, SRKMIX

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT_RR = ROOT / "testdata" / "m3_rachford_rice.json"
OUT_PIP = ROOT / "testdata" / "m3_pip.json"
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


def main():
    write_pip()
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


if __name__ == "__main__":
    main()
