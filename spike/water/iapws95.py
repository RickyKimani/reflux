# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = ["iapws==1.5.4", "mpmath==1.3.0"]
# ///
"""M5a plan spike: IAPWS-95 under reflux's generic volume solver, and the
non-analytic terms near the critical point.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> spike/water/iapws95.py [part ...]

Parts (all by default):
  check    my transcription of Table 5 against the iapws package's, and
           Table 6's printed values;
  critical the six derivatives of the residual part near (tau, delta) =
           (1, 1), float64 against mpmath's numerical derivatives at 40 digits;
  volume   thermo.solve_tp's three volume solvers (liquid_volume,
           vapor_volume, stable_volume in thermo/model.v), ported line for
           line, on a (T, P) grid, each answer against a brute-force scan of
           every root of P(rho) = p and the one of least Gibbs energy;
  sat      thermo.saturation_pressure (successive substitution from
           Edmister's start), ported, against Table 8 and the iapws package;
  gas      the gas constant: ln phi of pure water at low pressure when F
           carries IAPWS-95's R against reflux's (option a of the plan).

The port follows reflux's model contract: F(T, V, n) = A^r/(R T) with R
reflux's 8.314462618, so F = n (R95/R) phi^r + n (R95/R - 1) ln delta, which
makes P(T, rho) IAPWS-95's exactly. One mole throughout: V is molar.
"""
import math
import sys
import time

import mpmath
from iapws import IAPWS95

# IAPWS-95's constants, as the release gives them.
TC = 647.096  # K
RHOC_MASS = 322.0  # kg/m3
M = 0.018015268  # kg/mol (IAPWS-95's 18.015268 g/mol)
R95 = 0.46151805e3 * M  # J/(mol K): 8.314371357587
R = 8.314462618  # reflux's thermo.r_gas
RHOC = RHOC_MASS / M  # mol/m3

C = IAPWS95._constants
EPS = 2.220446049250313e-16


# --- The residual part, Table 5, in float64 -------------------------------

def phir(tau, delta):
    """phi^r and its six derivatives: f, f_d, f_dd, f_t, f_tt, f_dt."""
    f = fd = fdd = ft = ftt = fdt = 0.0
    for n, d, t in zip(C["nr1"], C["d1"], C["t1"]):
        dd, tt = delta ** d, tau ** t
        f += n * dd * tt
        fd += n * d * delta ** (d - 1) * tt
        fdd += n * d * (d - 1) * delta ** (d - 2) * tt
        ft += n * t * dd * tau ** (t - 1)
        ftt += n * t * (t - 1) * dd * tau ** (t - 2)
        fdt += n * d * t * delta ** (d - 1) * tau ** (t - 1)
    for n, d, t, c in zip(C["nr2"], C["d2"], C["t2"], C["c2"]):
        dc = delta ** c
        e = math.exp(-dc)
        base = n * delta ** d * tau ** t * e
        a = d - c * dc
        f += base
        fd += base * a / delta
        # (d - 1 - c dc), not (a - 1): at a low density a - 1 cancels.
        fdd += base * (a * (d - 1 - c * dc) - c * c * dc) / delta ** 2
        ft += base * t / tau
        ftt += base * t * (t - 1) / tau ** 2
        fdt += base * a * t / (delta * tau)
    for n, d, t, al, be, ga, ep in zip(C["nr3"], C["d3"], C["t3"], C["alfa3"], C["beta3"],
                                       C["gamma3"], C["epsilon3"]):
        e = math.exp(-al * (delta - ep) ** 2 - be * (tau - ga) ** 2)
        base = n * delta ** d * tau ** t * e
        xd = d / delta - 2 * al * (delta - ep)
        xt = t / tau - 2 * be * (tau - ga)
        f += base
        fd += base * xd
        fdd += base * (xd * xd - d / delta ** 2 - 2 * al)
        ft += base * xt
        ftt += base * (xt * xt - t / tau ** 2 - 2 * be)
        fdt += base * xd * xt
    for n, a, b, A, B, Cc, D, bt in zip(C["nr4"], C["a4"], C["b4"], C["A"], C["B"], C["C"],
                                         C["D"], C["beta4"]):
        dm = delta - 1.0
        d2 = dm * dm
        theta = (1 - tau) + A * d2 ** (0.5 / bt)
        Delta = theta * theta + B * d2 ** a
        psi = math.exp(-Cc * d2 - D * (tau - 1) ** 2)
        psid = -2 * Cc * dm * psi
        psidd = (2 * Cc * d2 - 1) * 2 * Cc * psi
        psit = -2 * D * (tau - 1) * psi
        psitt = (2 * D * (tau - 1) ** 2 - 1) * 2 * D * psi
        psidt = 4 * Cc * D * dm * (tau - 1) * psi
        # dDelta/ddelta over (delta - 1): finite, 0 at delta = 1, so the
        # release's (1/(delta-1)) dDelta/ddelta in d2Delta/ddelta2 is
        # written without the 0/0.
        g = A * theta * 2 / bt * (d2 ** (0.5 / bt - 1) if dm != 0 else 0.0) \
            + 2 * B * a * d2 ** (a - 1)
        Deltad = dm * g
        Deltadd = g + 4 * B * a * (a - 1) * d2 ** (a - 1) \
            + 2 * A * A / bt ** 2 * d2 ** (1 / bt - 1) \
            + A * theta * 4 / bt * (0.5 / bt - 1) * d2 ** (0.5 / bt - 1)
        Db = Delta ** b
        Db1 = Delta ** (b - 1) if Delta != 0 else math.inf
        Db2 = Delta ** (b - 2) if Delta != 0 else math.inf
        Dbd = b * Db1 * Deltad
        Dbdd = b * (Db1 * Deltadd + (b - 1) * Db2 * Deltad ** 2)
        Dbt = -2 * theta * b * Db1
        Dbtt = 2 * b * Db1 + 4 * theta * theta * b * (b - 1) * Db2
        Dbdt = -A * b * 2 / bt * Db1 * dm * (d2 ** (0.5 / bt - 1) if dm != 0 else 0.0) \
            - 2 * theta * b * (b - 1) * Db2 * Deltad
        f += n * Db * delta * psi
        fd += n * (Db * (psi + delta * psid) + Dbd * delta * psi)
        fdd += n * (Db * (2 * psid + delta * psidd) + 2 * Dbd * (psi + delta * psid)
                    + Dbdd * delta * psi)
        ft += n * delta * (Dbt * psi + Db * psit)
        ftt += n * delta * (Dbtt * psi + 2 * Dbt * psit + Db * psitt)
        fdt += n * (Db * (psit + delta * psidt) + delta * Dbd * psit + Dbt * (psi + delta * psid)
                    + Dbdt * delta * psi)
    return f, fd, fdd, ft, ftt, fdt


def phir_mp(tau, delta):
    """phi^r alone, in mpmath, for its derivatives by mpmath.diff."""
    mp = mpmath
    f = mp.mpf(0)
    for n, d, t in zip(C["nr1"], C["d1"], C["t1"]):
        f += mp.mpf(n) * delta ** d * tau ** mp.mpf(t)
    for n, d, t, c in zip(C["nr2"], C["d2"], C["t2"], C["c2"]):
        f += mp.mpf(n) * delta ** d * tau ** t * mp.exp(-delta ** c)
    for n, d, t, al, be, ga, ep in zip(C["nr3"], C["d3"], C["t3"], C["alfa3"], C["beta3"],
                                       C["gamma3"], C["epsilon3"]):
        f += mp.mpf(n) * delta ** d * tau ** t * mp.exp(
            -al * (delta - ep) ** 2 - be * (tau - mp.mpf(ga)) ** 2)
    for n, a, b, A, B, Cc, D, bt in zip(C["nr4"], C["a4"], C["b4"], C["A"], C["B"], C["C"],
                                         C["D"], C["beta4"]):
        d2 = (delta - 1) ** 2
        theta = (1 - tau) + mp.mpf(A) * d2 ** (mp.mpf(1) / (2 * mp.mpf(bt)))
        Delta = theta ** 2 + mp.mpf(B) * d2 ** mp.mpf(a)
        f += mp.mpf(n) * Delta ** mp.mpf(b) * delta * mp.exp(-Cc * d2 - D * (tau - 1) ** 2)
    return f


# --- check ----------------------------------------------------------------

def part_check():
    w = IAPWS95()
    worst = 0.0
    import random
    rnd = random.Random(1)
    for _ in range(2000):
        tau = rnd.uniform(0.25, 2.6)
        delta = rnd.choice([rnd.uniform(1e-6, 4.0), 10 ** rnd.uniform(-8, 0.6)])
        mine = phir(tau, delta)
        ref = w._phir(tau, delta)
        theirs = [ref["fir"], ref["fird"], ref["firdd"], ref["firt"], ref["firtt"], ref["firdt"]]
        for x, y in zip(mine, theirs):
            worst = max(worst, abs(x - y) / max(abs(y), 1e-300))
    print(f"check: my Table 5 against the iapws package's, 2000 random states: "
          f"worst relative difference {worst:.2e}")
    tau, delta = TC / 500.0, 838.025 / RHOC_MASS
    f, fd, fdd, ft, ftt, fdt = phir(tau, delta)
    table6 = dict(f=-0.342693206e1, fd=-0.364366650, fdd=0.856063701, ft=-0.581403435e1,
                  ftt=-0.223440737e1, fdt=-0.112176915e1)
    for k, v in zip(["f", "fd", "fdd", "ft", "ftt", "fdt"], [f, fd, fdd, ft, ftt, fdt]):
        print(f"  Table 6 phi^r_{k:3}: {v: .9e}  printed {table6[k]: .9e}  "
              f"rel {abs(v - table6[k]) / abs(table6[k]):.1e}")


# --- critical -------------------------------------------------------------

def part_critical():
    mpmath.mp.dps = 40
    names = ["f", "f_d", "f_dd", "f_t", "f_tt", "f_dt"]
    orders = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0), (1, 1)]  # (tau, delta)
    print("critical: float64 Table 5 against mpmath.diff of phi^r at 40 digits")
    print("  direction (dtau, ddelta), distance s: worst relative error of the six; "
          "f_tt's value")
    dirs = [(0.0, 1.0), (1.0, 0.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 1.0), (-1.0, 0.5)]
    worst_all = 0.0
    for dt, dd in dirs:
        line = []
        for k in range(1, 13, 1):
            s = 10.0 ** (-k)
            tau, delta = 1.0 + dt * s, 1.0 + dd * s
            got = phir(tau, delta)
            tm, dm_ = mpmath.mpf(tau), mpmath.mpf(delta)
            errs = []
            for (ot, od), g in zip(orders, got):
                ref = mpmath.diff(lambda a, b: phir_mp(a, b), (tm, dm_), (ot, od))
                errs.append(float(abs((mpmath.mpf(g) - ref) / ref)) if ref != 0 else float(abs(g)))
            e = max(errs)
            worst_all = max(worst_all, e if k <= 10 else 0)
            line.append(f"1e-{k}:{e:.0e}")
        print(f"  ({dt:+.0f},{dd:+.1f}) " + " ".join(line))
    # Along the critical isochore and isotherm, how f_tt and f_dd grow.
    print("  growth on the critical isochore (delta = 1, tau = 1 + s) and isotherm "
          "(tau = 1, delta = 1 + s):")
    for k in [2, 4, 6, 8, 10, 12, 14]:
        s = 10.0 ** (-k)
        a = phir(1 + s, 1.0)
        b = phir(1.0, 1 + s)
        print(f"    s = 1e-{k:<2}  isochore f_tt {a[4]: .6e} f_dd {a[2]: .6e} | "
              f"isotherm f_tt {b[4]: .6e} f_dd {b[2]: .6e}")
    at = phir(1.0, 1.0)
    print(f"  at the critical point itself: {dict(zip(names, at))}")


# --- volume ---------------------------------------------------------------

class Fixed:
    """thermo.Fixed: the model at fixed T (n = 1 mol)."""

    def __init__(self, t, vmin):
        self.t = t
        self.tau = TC / t
        self.vmin = vmin
        self.ntot = 1.0
        self.calls = 0

    def at(self, v):
        """F, F_V, F_VV at molar volume v."""
        self.calls += 1
        delta = 1.0 / (v * RHOC)
        f, fd, fdd, _, _, _ = phir(self.tau, delta)
        k = R95 / R
        F = k * f + (k - 1) * math.log(delta)
        Fd = k * fd + (k - 1) / delta
        Fdd = k * fdd - (k - 1) / delta ** 2
        Fv = -delta / v * Fd
        Fvv = (delta * delta * Fdd + 2 * delta * Fd) / (v * v)
        return F, Fv, Fvv


VOLUME_MAX_ITER = 100
VOLUME_RTOL = 1e-12


def p_dpdv(fx, v):
    _, fv, fvv = fx.at(v)
    rt = R * fx.t
    return (rt * (fx.ntot / v - fv), -rt * (fvv + fx.ntot / (v * v)),
            16.0 * EPS * rt * (fx.ntot / v + abs(fv)))


def liquid_volume(fx, p):
    vmin = fx.vmin
    v = 1.05 * vmin
    for _ in range(60):
        if p_dpdv(fx, v)[0] > p:
            break
        v = vmin + 0.5 * (v - vmin)
    for _ in range(VOLUME_MAX_ITER):
        sp, sd, noise = p_dpdv(fx, v)
        if sd >= 0.0:
            return None
        dv = -(sp - p) / sd
        if abs(sp - p) <= noise and abs(dv) > VOLUME_RTOL * v:
            return v
        dv = min(max(dv, -0.5 * (v - vmin)), 0.5 * (v - vmin))
        v += dv
        if abs(dv) <= VOLUME_RTOL * v:
            return v
    return None


def vapor_volume(fx, p):
    ntot = fx.ntot
    rho_max = ntot / fx.vmin
    rho = min(p / (R * fx.t), 0.5 * rho_max)
    for _ in range(VOLUME_MAX_ITER):
        v = ntot / rho
        sp, sd, noise = p_dpdv(fx, v)
        dp_drho = -sd * v * v / ntot
        if dp_drho <= 0.0:
            return None
        drho = -(sp - p) / dp_drho
        if abs(sp - p) <= noise and abs(drho) > VOLUME_RTOL * rho:
            return v
        drho = min(max(drho, -0.5 * rho), min(0.5 * rho, 0.5 * (rho_max - rho)))
        rho += drho
        if abs(drho) <= VOLUME_RTOL * rho:
            if rho_max - rho <= 1e-9 * rho_max:
                return None
            return ntot / rho
    return None


def gibbs_residual(fx, p, v):
    z = p * v / (fx.ntot * R * fx.t)
    return fx.at(v)[0] + fx.ntot * (z - 1.0) - fx.ntot * math.log(z)


def stable_volume(fx, p):
    vl = liquid_volume(fx, p)
    if vl is None:
        return vapor_volume(fx, p), "vapor only"
    vv = vapor_volume(fx, p)
    if vv is None:
        return vl, "liquid only"
    return (vl, "liquid") if gibbs_residual(fx, p, vl) < gibbs_residual(fx, p, vv) else (vv, "vapor")


def pressure(t, rho):
    """P from IAPWS-95 at molar density rho, Pa."""
    f, fd, *_ = phir(TC / t, rho / RHOC)
    return rho * R95 * t * (1 + rho / RHOC * fd)


def brute_roots(t, p, grid, pgrid):
    """Every root of P(rho) = p on the scan, refined by bisection, with
    dP/drho there and the residual Gibbs energy."""
    roots = []
    for i in range(len(grid) - 1):
        a, b = pgrid[i] - p, pgrid[i + 1] - p
        if a == 0 or (a < 0) != (b < 0):
            lo, hi = grid[i], grid[i + 1]
            flo = a
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                fm = pressure(t, mid) - p
                if (fm < 0) == (flo < 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
                if hi - lo <= 4 * EPS * hi:
                    break
            rho = 0.5 * (lo + hi)
            h = 1e-7 * rho
            slope = (pressure(t, rho + h) - pressure(t, rho - h)) / (2 * h)
            roots.append((rho, slope))
    return roots


def bisect_root(t, p, lo, hi):
    flo = pressure(t, lo) - p
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = pressure(t, mid) - p
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo <= 4 * EPS * hi:
            break
    return 0.5 * (lo + hi)


def part_volume():
    temps = sorted(set([251.165, 260, 273.16, 280, 300, 325, 350, 375, 400, 450, 500, 550, 600,
                        625, 640, 645, 646.5, 647.0, 647.09, 647.096, 647.1, 647.5, 650, 660, 700,
                        800, 900, 1073.15, 1273, 1500, 2000, 2273.15, 5000]))
    pressures = [10 ** (k / 8) for k in range(0, 8 * 9 + 1)]  # 1 Pa to 1 GPa
    pressures += [2e9, 1e10, 1e11]
    for dmax in [5.0]:
        vmin = 1.0 / (dmax * RHOC)
        print(f"volume: vmin at delta = {dmax} ({dmax * RHOC_MASS:.0f} kg/m3); "
              f"{len(temps)} temperatures x {len(pressures)} pressures")
        for t in [251.165, 300, 647.096, 1273, 5000]:
            print(f"  P at delta 5, 8, 10, T = {t:8.3f} K: "
                  + ", ".join(f"{pressure(t, d * RHOC):.3e}" for d in (5, 8, 10)) + " Pa")
        n_ok = n_bad = n_none = n_trap = 0
        bad, traps = [], []
        calls = 0
        t0 = time.perf_counter()
        for t in temps:
            # Scan rho from 1e-9 of critical to vmin's, log-spaced, fine.
            npts = 60000
            lo, hi = math.log(1e-9 * RHOC), math.log(1 / vmin)
            grid = [math.exp(lo + (hi - lo) * i / (npts - 1)) for i in range(npts)]
            pgrid = [pressure(t, r) for r in grid]
            # The physical branches: the vapour's, rising from rho -> 0 until P
            # first stops rising, and the liquid's, rising to vmin from where
            # P last stopped rising. Between them (when they differ) lie the
            # spinodal region and IAPWS-95's spurious loops.
            ev = 0
            while ev + 1 < npts and pgrid[ev + 1] > pgrid[ev]:
                ev += 1
            sl = npts - 1
            while sl > 0 and pgrid[sl - 1] < pgrid[sl]:
                sl -= 1
            for p in pressures:
                fx = Fixed(t, vmin)
                try:
                    v, which = stable_volume(fx, p)
                except (OverflowError, ValueError, ZeroDivisionError) as e:
                    v, which = None, f"exception {e}"
                calls += fx.calls
                cands = []
                for a, b in [(0, ev), (sl, npts - 1)]:
                    if pgrid[a] <= p <= pgrid[b]:
                        for i in range(a, b):
                            if pgrid[i] <= p <= pgrid[i + 1]:
                                cands.append(bisect_root(t, p, grid[i], grid[i + 1]))
                                break
                cands = sorted(set(cands))
                if not cands:
                    if v is None:
                        n_none += 1
                    else:
                        n_bad += 1
                        bad.append((t, p, "no root on either branch, solver gave", 1 / v * M,
                                    which))
                    continue
                fxg = Fixed(t, vmin)
                gbest, rbest = min((gibbs_residual(fxg, p, 1 / r), r) for r in cands)
                # Spurious roots: mechanically stable, off both branches, and
                # lower in Gibbs energy than the physical answer.
                if ev < sl:
                    for i in range(ev, sl):
                        if (pgrid[i] - p) * (pgrid[i + 1] - p) < 0 and pgrid[i + 1] > pgrid[i]:
                            r = bisect_root(t, p, grid[i], grid[i + 1])
                            g = gibbs_residual(fxg, p, 1 / r)
                            if g < gbest:
                                n_trap += 1
                                traps.append((t, p, r * M, g - gbest))
                                break
                if v is None:
                    n_bad += 1
                    bad.append((t, p, "solver found none; the branches' answer at rho =",
                                rbest * M, f"{len(cands)} branch roots"))
                    continue
                rho = 1 / v
                if abs(rho - rbest) <= 1e-7 * rbest:
                    n_ok += 1
                else:
                    g_solver = gibbs_residual(fxg, p, v)
                    if abs(g_solver - gbest) <= 1e-9:
                        n_ok += 1  # on the saturation line, within rounding
                    else:
                        n_bad += 1
                        bad.append((t, p, f"solver {which} rho = {rho * M:.9g}, branches'",
                                    rbest * M, f"dG/RT = {g_solver - gbest:.3e}; branch roots "
                                    f"(kg/m3) {[round(r * M, 6) for r in cands]}"))
        dt = time.perf_counter() - t0
        print(f"  agree {n_ok}, disagree {n_bad}, neither has a root {n_none}; "
              f"{calls} model calls ({dt:.0f} s with the scans)")
        for b in bad[:40]:
            print("   ", b)
        kinds = {}
        for b in bad:
            k = b[2].split(" rho")[0].split(";")[0]
            kinds[k] = kinds.get(k, 0) + 1
        print(f"  disagreements by kind: {kinds}")
        print(f"  (again) agree {n_ok}, disagree {n_bad}, neither has a root {n_none}")
        print(f"  states with a spurious, mechanically stable root of lower G: {n_trap}")
        tt = sorted(set(x[0] for x in traps))
        print(f"  at temperatures {tt}")
        for x in traps[:3] + traps[-3:]:
            print(f"    T {x[0]} K, P {x[1]:.3g} Pa: rho {x[2]:.4f} kg/m3, G below by {x[3]:.3e} RT")


def part_sat():
    w = IAPWS95()
    vmin = 1.0 / (5.0 * RHOC)
    print("sat: thermo.saturation_pressure ported (P <- P phiL/phiV from Edmister's start)")
    table8 = {275.0: 0.698451167e-3, 450.0: 0.932203564, 625.0: 0.169082693e2}
    omega, pc = 0.3443, 22.064e6
    for t in [273.16, 275.0, 300, 350, 400, 450, 500, 550, 600, 625.0, 640, 645, 646, 646.5,
              647.0, 647.05, 647.09, 647.095]:
        fx = Fixed(t, vmin)
        p = pc * 10 ** (7 / 3 * (1 + omega) * (1 - TC / t))
        it = 0
        status = "no convergence in 500"
        for it in range(1, 501):
            vl = liquid_volume(fx, p)
            if vl is None:
                p *= 1.1
                continue
            vv = vapor_volume(fx, p)
            if vv is None:
                p *= 0.9
                continue
            if abs(vl - vv) <= 1e-8 * vv:
                status = "lost the two-phase region"
                break
            diff = gibbs_residual(fx, p, vl) - gibbs_residual(fx, p, vv)
            p *= math.exp(diff)
            if abs(diff) < 1e-10:
                status = "converged"
                break
        ref = w._Vapor_Pressure(t) if hasattr(w, "_Vapor_Pressure") else float("nan")
        extra = ""
        if t in table8:
            extra = f" Table 8 {table8[t] * 1e6:.9e} rel {abs(p / (table8[t] * 1e6) - 1):.1e}"
        print(f"  T = {t:8.3f} K: {status:26} {it:3} iterations, {fx.calls:5} calls, "
              f"P = {p:.9e} Pa (Wagner-Pruss {ref * 1e6:.6e}){extra}")


def part_gas():
    print("gas: pure water's ln phi as P -> 0 with F carrying R95 (option a)")
    print(f"  R95 = {R95:.12f} J/(mol K), R95/R - 1 = {R95 / R - 1:.6e}")
    for t in [300.0, 373.15, 500.0, 1000.0]:
        for p in [1.0, 1e3, 1e5]:
            fx = Fixed(t, 1.0 / (5.0 * RHOC))
            v = vapor_volume(fx, p)
            if v is None:
                continue
            delta = 1 / (v * RHOC)
            f, fd, *_ = phir(TC / t, delta)
            k = R95 / R
            fn = k * (f + delta * fd) + (k - 1) * (math.log(delta) + 1)
            z = p * v / (R * t)
            # The same state with reflux's R throughout (option b): phi^r alone.
            z_b = 1 + delta * fd
            lnphi_b = f + delta * fd - math.log(z_b)
            print(f"  T = {t:7.2f} K, P = {p:8.0e} Pa: ln phi (a) {fn - math.log(z): .4e}, "
                  f"(b) {lnphi_b: .4e}, difference {fn - math.log(z) - lnphi_b: .3e}")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["check", "critical", "volume", "sat", "gas"]
    for part in parts:
        globals()["part_" + part]()
        print()
