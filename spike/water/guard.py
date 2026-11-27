# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "mpmath==1.3.0"]
# ///
"""M5a plan spike: a branch guard for thermo's volume solvers, and a
bracketed start for saturation_pressure, on IAPWS-95 near Tc.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/water/guard.py [d7 | d7b | d7b-near]

The first guard (no argument, and `d7` on the whole grid): on its own
branch, P(rho) is concave for the vapour and P(V) convex for the liquid, so
Newton from each solver's end never crosses p: the vapour's iterates keep
P <= p and the liquid's P >= p. An iterate on the wrong side of p (beyond
rounding) has jumped the spinodal gap onto the other branch, and the solver
returns none, as it does at a spinodal. It works near Tc, but the grid
(`d7`) found IAPWS-95's liquid branch isn't convex at 251 to 260 K, where
it refused real liquid roots.

The second (`d7b` on the grid, `d7b-near` near Tc), D7 as the plan has it:
below Tc a vapour iterate denser than rho_c, or a liquid iterate lighter
than it, has left its branch; the vapour's start is bounded by the
auxiliary saturated-vapour density; a liquid walk that never reaches p
above vmin returns none.

The start: saturation_pressure's x1.1 and x0.9 steps, but once a pressure
with no vapour root and one with no liquid root are both known, the next is
their geometric mean, which finds a window of metastability however narrow.
"""
import importlib.util
import math
import pathlib

spec = importlib.util.spec_from_file_location(
    "s", pathlib.Path(__file__).with_name("iapws95.py"))
s = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)
from iapws import IAPWS95  # noqa: E402

EPS = s.EPS


def liquid_volume_g(fx, p):
    vmin = fx.vmin
    v = 1.05 * vmin
    for _ in range(60):
        if s.p_dpdv(fx, v)[0] > p:
            break
        v = vmin + 0.5 * (v - vmin)
    for _ in range(s.VOLUME_MAX_ITER):
        sp, sd, noise = s.p_dpdv(fx, v)
        if sd >= 0.0:
            return None
        if sp < p - noise:
            return None  # the guard: jumped onto the vapour branch
        dv = -(sp - p) / sd
        if abs(sp - p) <= noise and abs(dv) > s.VOLUME_RTOL * v:
            return v
        dv = min(max(dv, -0.5 * (v - vmin)), 0.5 * (v - vmin))
        v += dv
        if abs(dv) <= s.VOLUME_RTOL * v:
            return v
    return None


def vapor_volume_g(fx, p):
    ntot = fx.ntot
    rho_max = ntot / fx.vmin
    rho = min(p / (s.R * fx.t), 0.5 * rho_max)
    for _ in range(s.VOLUME_MAX_ITER):
        v = ntot / rho
        sp, sd, noise = s.p_dpdv(fx, v)
        dp_drho = -sd * v * v / ntot
        if dp_drho <= 0.0:
            return None
        if sp > p + noise:
            return None  # the guard: jumped onto the liquid branch
        drho = -(sp - p) / dp_drho
        if abs(sp - p) <= noise and abs(drho) > s.VOLUME_RTOL * rho:
            return v
        drho = min(max(drho, -0.5 * rho), min(0.5 * rho, 0.5 * (rho_max - rho)))
        rho += drho
        if abs(drho) <= s.VOLUME_RTOL * rho:
            if rho_max - rho <= 1e-9 * rho_max:
                return None
            return ntot / rho
    return None


def psat(t, guarded, bracket):
    vmin = 1.0 / (5.0 * s.RHOC)
    fx = s.Fixed(t, vmin)
    lv = liquid_volume_g if guarded else s.liquid_volume
    vv_ = vapor_volume_g if guarded else s.vapor_volume
    p = 22.064e6 * 10 ** (7 / 3 * (1 + 0.3443) * (1 - s.TC / t))
    lo = hi = None  # a P with no liquid root (too low), one with no vapour root
    for it in range(1, 501):
        vl = lv(fx, p)
        if vl is None:
            lo = p
            p = math.sqrt(lo * hi) if bracket and hi else p * 1.1
            continue
        vv = vv_(fx, p)
        if vv is None:
            hi = p
            p = math.sqrt(lo * hi) if bracket and lo else p * 0.9
            continue
        if abs(vl - vv) <= 1e-8 * vv:
            return None, it, "lost the two-phase region", fx.calls
        diff = s.gibbs_residual(fx, p, vl) - s.gibbs_residual(fx, p, vv)
        p *= math.exp(diff)
        if abs(diff) < 1e-10:
            return p, it, "converged", fx.calls
    return None, 500, "no convergence", fx.calls


def main():
    print("saturation_pressure near Tc, IAPWS-95 (iapws package's Psat for reference)")
    print("  T K        | unguarded           | guard alone         | guard and bracket")
    for t in [600.0, 625.0, 640.0, 645.0, 646.0, 646.5, 647.0, 647.05, 647.09, 647.095,
              647.0959]:
        ref = IAPWS95(T=t, x=0).P * 1e6
        row = []
        for g, b in [(False, False), (True, False), (True, True)]:
            p, it, st, calls = psat(t, g, b)
            row.append(f"{st[:9]:9} {it:3} it" + (f" {abs(p / ref - 1):.0e}" if p else "      "))
        print(f"  {t:10.4f} | " + " | ".join(row))
    # Branch roots: below Tc, at pressures across the metastable window and
    # beyond, the guarded solvers' roots against the scan's branch roots.
    print("branch roots near Tc: each solver's root against the branches' (none where its"
          " branch doesn't reach p)")
    vmin = 1.0 / (5.0 * s.RHOC)
    wrong_u = wrong_g = total = 0
    for t in [600.0, 640.0, 645.0, 646.0, 647.0, 647.09]:
        ps = IAPWS95(T=t, x=0).P * 1e6
        import numpy as np
        grid = np.exp(np.linspace(math.log(1e-6 * s.RHOC), math.log(1 / vmin), 200001))
        P = np.array([s.pressure(t, r) for r in grid])
        ev = int(np.argmax(np.diff(P) <= 0))
        sl = len(P) - 1 - int(np.argmax(np.diff(P)[::-1] <= 0))
        for f in np.linspace(0.97, 1.03, 61):
            p = ps * f
            fx = s.Fixed(t, vmin)
            branch_v = s.bisect_root(t, p, grid[np.searchsorted(P[:ev + 1], p) - 1],
                                     grid[np.searchsorted(P[:ev + 1], p)]) \
                if P[0] <= p <= P[ev] else None
            branch_l = None
            if P[sl] <= p <= P[-1]:
                j = sl + int(np.searchsorted(P[sl:], p))
                branch_l = s.bisect_root(t, p, grid[j - 1], grid[j])
            for name, solver, want in [("liquid", s.liquid_volume, branch_l),
                                       ("vapour", s.vapor_volume, branch_v),
                                       ("liquid guarded", liquid_volume_g, branch_l),
                                       ("vapour guarded", vapor_volume_g, branch_v)]:
                v = solver(fx, p)
                got = None if v is None else 1 / v
                ok = (got is None and want is None) or (
                    got is not None and want is not None and abs(got / want - 1) < 1e-8)
                if "guarded" in name:
                    wrong_g += not ok
                else:
                    wrong_u += not ok
                    total += 1
    print(f"  {total} (T, p) states x 2 solvers: wrong branch root unguarded {wrong_u}, "
          f"guarded {wrong_g}")




# --- D7 as proposed, on the volume grid ------------------------------------

def vapor_start_cap(t):
    """Wagner and Pruss's auxiliary saturated-vapour density (SR1-86), mol/m3,
    below Tc; none above it."""
    if t >= s.TC:
        return None
    return IAPWS95._Vapor_Density(t) / s.M


def liquid_volume_d7(fx, p):
    vmin = fx.vmin
    v = 1.05 * vmin
    reached = False
    for _ in range(60):
        if s.p_dpdv(fx, v)[0] > p:
            reached = True
            break
        v = vmin + 0.5 * (v - vmin)
    if not reached:
        return None  # the branch never reaches p above vmin: refused, not vmin
    last = None
    for _ in range(s.VOLUME_MAX_ITER):
        sp, sd, noise = s.p_dpdv(fx, v)
        if sd >= 0.0:
            return None
        if last is not None and last > p and sp < p - noise:
            return None  # the guard: a step up in V carried it past p
        last = sp
        dv = -(sp - p) / sd
        if abs(sp - p) <= noise and abs(dv) > s.VOLUME_RTOL * v:
            return v
        dv = min(max(dv, -0.5 * (v - vmin)), 0.5 * (v - vmin))
        v += dv
        if abs(dv) <= s.VOLUME_RTOL * v:
            return v
    return None


def vapor_volume_d7(fx, p):
    ntot = fx.ntot
    rho_max = ntot / fx.vmin
    rho = min(p / (s.R * fx.t), 0.5 * rho_max)
    cap = vapor_start_cap(fx.t)
    if cap is not None:
        rho = min(rho, cap)
    last = None
    for _ in range(s.VOLUME_MAX_ITER):
        v = ntot / rho
        sp, sd, noise = s.p_dpdv(fx, v)
        dp_drho = -sd * v * v / ntot
        if dp_drho <= 0.0:
            return None
        if last is not None and last < p and sp > p + noise:
            return None  # the guard: a step up in rho carried it past p
        last = sp
        drho = -(sp - p) / dp_drho
        if abs(sp - p) <= noise and abs(drho) > s.VOLUME_RTOL * rho:
            return v
        drho = min(max(drho, -0.5 * rho), min(0.5 * rho, 0.5 * (rho_max - rho)))
        rho += drho
        if abs(drho) <= s.VOLUME_RTOL * rho:
            if rho_max - rho <= 1e-9 * rho_max:
                return None
            return ntot / rho
    return None


def d7_volume():
    """The volume grid of iapws95.py's `volume`, with D7's first solvers.
    Found: the liquid guard refuses real liquid roots at 251 to 260 K, where
    IAPWS-95's liquid branch isn't convex between vmin and the root."""
    s.liquid_volume = liquid_volume_d7
    s.vapor_volume = vapor_volume_d7
    s.part_volume()


# --- D7 as amended: the branches split at rho_c below Tc ---------------------
#
# Below Tc the spinodal curve has its vapour side below rho_c and its liquid
# side above it, so a vapour iterate denser than rho_c, or a liquid iterate
# lighter than it, has left its branch. That needs no convexity, and only a
# model that gives its split density (and its vapour start) takes it: a
# model without them takes the solvers' paths unchanged.

def split_density(t):
    return s.RHOC if t < s.TC else None


def liquid_volume_d7b(fx, p):
    vmin = fx.vmin
    v = 1.05 * vmin
    reached = False
    for _ in range(60):
        if s.p_dpdv(fx, v)[0] > p:
            reached = True
            break
        v = vmin + 0.5 * (v - vmin)
    if not reached:
        return None  # the branch never reaches p above vmin: refused, not vmin
    split = split_density(fx.t)
    for _ in range(s.VOLUME_MAX_ITER):
        if split is not None and fx.ntot / v < split:
            return None  # lighter than rho_c below Tc: off the liquid's branch
        sp, sd, noise = s.p_dpdv(fx, v)
        if sd >= 0.0:
            return None
        dv = -(sp - p) / sd
        if abs(sp - p) <= noise and abs(dv) > s.VOLUME_RTOL * v:
            return v
        dv = min(max(dv, -0.5 * (v - vmin)), 0.5 * (v - vmin))
        v += dv
        if abs(dv) <= s.VOLUME_RTOL * v:
            if split is not None and fx.ntot / v < split:
                return None
            return v
    return None


def vapor_volume_d7b(fx, p):
    ntot = fx.ntot
    rho_max = ntot / fx.vmin
    rho = min(p / (s.R * fx.t), 0.5 * rho_max)
    cap = vapor_start_cap(fx.t)
    if cap is not None:
        rho = min(rho, cap)
    split = split_density(fx.t)
    for _ in range(s.VOLUME_MAX_ITER):
        if split is not None and rho > split:
            return None  # denser than rho_c below Tc: off the vapour's branch
        v = ntot / rho
        sp, sd, noise = s.p_dpdv(fx, v)
        dp_drho = -sd * v * v / ntot
        if dp_drho <= 0.0:
            return None
        drho = -(sp - p) / dp_drho
        if abs(sp - p) <= noise and abs(drho) > s.VOLUME_RTOL * rho:
            return v
        drho = min(max(drho, -0.5 * rho), min(0.5 * rho, 0.5 * (rho_max - rho)))
        rho += drho
        if abs(drho) <= s.VOLUME_RTOL * rho:
            if rho_max - rho <= 1e-9 * rho_max:
                return None
            if split is not None and rho > split:
                return None
            return ntot / rho
    return None


def d7b_volume():
    """The volume grid with D7's amended solvers."""
    s.liquid_volume = liquid_volume_d7b
    s.vapor_volume = vapor_volume_d7b
    s.part_volume()


def d7b_near():
    """main()'s near-Tc checks, with D7's amended solvers in place of the
    guarded ones: Psat with the bracket, and the branch roots."""
    global liquid_volume_g, vapor_volume_g
    liquid_volume_g, vapor_volume_g = liquid_volume_d7b, vapor_volume_d7b
    main()


if __name__ == "__main__":
    import sys

    # `guard.py` runs the near-Tc comparison; `d7` the volume grid with D7's
    # first solvers; `d7b` the grid with the amended ones, `d7b-near` the
    # near-Tc comparison with them (in the "guard" columns).
    {"d7": d7_volume, "d7b": d7b_volume, "d7b-near": d7b_near}.get(
        (sys.argv[1:] or ["main"])[0], main)()
