# /// script
# requires-python = ">=3.11"
# dependencies = ["chemicals==1.5.2", "mpmath==1.3.0"]
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

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT_RR = ROOT / "testdata" / "m3_rachford_rice.json"
SEED = 20261027

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


def main():
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
