# /// script
# requires-python = ">=3.11,<3.14"
# dependencies = ["sympy==1.13.3", "mpmath==1.3.0"]
# ///
"""M8 step 9's reference (criterion 7): Hock and Schittkowski's problems,
their published starts and optima, and each optimum again at 40 digits.

W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming Codes,
Lecture Notes in Economics and Mathematical Systems 187, Springer, 1981.

Each problem is written here in SymPy, as the book states it:
- minimise f(x),
- equalities c(x) = 0 first, then inequalities c(x) >= 0,
- bounds lower <= x <= upper, None where there's none.

Its optimum is found again from the published one: the first-order
conditions on the active set the book's optimum has, solved by mpmath's
Newton at 40 digits with SymPy's exact gradients. Each answer is checked
before it's written:
- within 1e-6 of the published x* and f* (the book prints 7 to 10 figures),
  which catches a problem mistyped here;
- every inactive inequality and bound strictly satisfied, every active
  inequality's and bound's multiplier strictly positive.

    PYTHONIOENCODING=utf-8 uv run --python <python 3.13> tools/gen_m8_hock_schittkowski.py

writes testdata/m8_hock_schittkowski.json. reflux's test writes each
problem again in V.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 40


def rosenbrock(x):
    return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2


# Each problem: n, f, equalities, inequalities, bounds, the start, the
# published optimum (x*, f*), and the active set there: the inequalities by
# index from 0, and the bounds as (variable, 'lower' or 'upper').
def problems():
    out = []

    def add(number, n, f, eq, ineq, lower, upper, x0, x_star, f_star, active, bounds):
        out.append(dict(number=number, n=n, f=f, eq=eq, ineq=ineq, lower=lower, upper=upper,
                        x0=x0, x_star=x_star, f_star=f_star, active=active, bounds=bounds))

    x = sp.symbols('x1:11')
    pi, sqrt = sp.pi, sp.sqrt
    add(1, 2, rosenbrock(x), [], [], [None, -1.5], [None, None], [-2, 1], [1, 1], 0, [], [])
    add(2, 2, rosenbrock(x), [], [], [None, 1.5], [None, None], [-2, 1],
        ['1.224370749', '1.5'], '0.05042618790', [], [(1, 'lower')])
    add(3, 2, x[1] + sp.Rational(1, 100000) * (x[1] - x[0]) ** 2, [], [], [None, 0], [None, None],
        [10, 1], [0, 0], 0, [], [(1, 'lower')])
    add(4, 2, (x[0] + 1) ** 3 / 3 + x[1], [], [], [1, 0], [None, None], ['1.125', '0.125'],
        [1, 0], sp.Rational(8, 3), [], [(0, 'lower'), (1, 'lower')])
    add(5, 2, sp.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - sp.Rational(3, 2) * x[0]
        + sp.Rational(5, 2) * x[1] + 1, [], [], [-1.5, -3], [4, 3], [0, 0],
        [-pi / 3 + sp.Rational(1, 2), -pi / 3 - sp.Rational(1, 2)], -sqrt(3) / 2 - pi / 3, [], [])
    add(6, 2, (1 - x[0]) ** 2, [10 * (x[1] - x[0] ** 2)], [], [None] * 2, [None] * 2,
        ['-1.2', 1], [1, 1], 0, [], [])
    add(7, 2, sp.log(1 + x[0] ** 2) - x[1], [(1 + x[0] ** 2) ** 2 + x[1] ** 2 - 4], [],
        [None] * 2, [None] * 2, [2, 2], [0, sqrt(3)], -sqrt(3), [], [])
    add(9, 2, sp.sin(pi * x[0] / 12) * sp.cos(pi * x[1] / 16), [4 * x[0] - 3 * x[1]], [],
        [None] * 2, [None] * 2, [0, 0], [-3, -4], sp.Rational(-1, 2), [], [])
    add(10, 2, x[0] - x[1], [], [-3 * x[0] ** 2 + 2 * x[0] * x[1] - x[1] ** 2 + 1], [None] * 2,
        [None] * 2, [-10, 10], [0, 1], -1, [0], [])
    add(11, 2, (x[0] - 5) ** 2 + x[1] ** 2 - 25, [], [-x[0] ** 2 + x[1]], [None] * 2, [None] * 2,
        ['4.9', '0.1'], ['1.234772825', '1.524663929'], '-8.498464223', [0], [])
    add(12, 2, x[0] ** 2 / 2 + x[1] ** 2 - x[0] * x[1] - 7 * x[0] - 7 * x[1], [],
        [25 - 4 * x[0] ** 2 - x[1] ** 2], [None] * 2, [None] * 2, [0, 0], [2, 3], -30, [0], [])
    add(14, 2, (x[0] - 2) ** 2 + (x[1] - 1) ** 2, [x[0] - 2 * x[1] + 1],
        [-x[0] ** 2 / 4 - x[1] ** 2 + 1], [None] * 2, [None] * 2, [2, 2],
        [(sqrt(7) - 1) / 2, (sqrt(7) + 1) / 4], 9 - sp.Rational(23, 8) * sqrt(7), [0], [])
    add(15, 2, rosenbrock(x), [], [x[0] * x[1] - 1, x[0] + x[1] ** 2], [None, None],
        ['0.5', None], [-2, 1], ['0.5', 2], '306.5', [0], [(0, 'upper')])
    add(21, 2, x[0] ** 2 / 100 + x[1] ** 2 - 100, [], [10 * x[0] - x[1] - 10], [2, -50],
        [50, 50], [-1, -1], [2, 0], '-99.96', [], [(0, 'lower')])
    add(28, 3, (x[0] + x[1]) ** 2 + (x[1] + x[2]) ** 2, [x[0] + 2 * x[1] + 3 * x[2] - 1], [],
        [None] * 3, [None] * 3, [-4, 1, 1], ['0.5', '-0.5', '0.5'], 0, [], [])
    add(35, 3, 9 - 8 * x[0] - 6 * x[1] - 4 * x[2] + 2 * x[0] ** 2 + 2 * x[1] ** 2 + x[2] ** 2
        + 2 * x[0] * x[1] + 2 * x[0] * x[2], [], [3 - x[0] - x[1] - 2 * x[2]], [0, 0, 0],
        [None] * 3, ['0.5', '0.5', '0.5'],
        [sp.Rational(4, 3), sp.Rational(7, 9), sp.Rational(4, 9)], sp.Rational(1, 9), [0], [])
    add(39, 4, -x[0], [x[1] - x[0] ** 3 - x[2] ** 2, x[0] ** 2 - x[1] - x[3] ** 2], [],
        [None] * 4, [None] * 4, [2, 2, 2, 2], [1, 1, 0, 0], -1, [], [])
    add(43, 4, x[0] ** 2 + x[1] ** 2 + 2 * x[2] ** 2 + x[3] ** 2 - 5 * x[0] - 5 * x[1]
        - 21 * x[2] + 7 * x[3], [],
        [8 - x[0] ** 2 - x[1] ** 2 - x[2] ** 2 - x[3] ** 2 - x[0] + x[1] - x[2] + x[3],
         10 - x[0] ** 2 - 2 * x[1] ** 2 - x[2] ** 2 - 2 * x[3] ** 2 + x[0] + x[3],
         5 - 2 * x[0] ** 2 - x[1] ** 2 - x[2] ** 2 - 2 * x[0] + x[1] + x[3]],
        [None] * 4, [None] * 4, [0, 0, 0, 0], [0, 1, 2, -1], -44, [0, 2], [])
    add(48, 5, (x[0] - 1) ** 2 + (x[1] - x[2]) ** 2 + (x[3] - x[4]) ** 2,
        [x[0] + x[1] + x[2] + x[3] + x[4] - 5, x[2] - 2 * (x[3] + x[4]) + 3], [], [None] * 5,
        [None] * 5, [3, 5, -3, 2, -2], [1, 1, 1, 1, 1], 0, [], [])
    add(71, 4, x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2],
        [x[0] ** 2 + x[1] ** 2 + x[2] ** 2 + x[3] ** 2 - 40], [x[0] * x[1] * x[2] * x[3] - 25],
        [1] * 4, [5] * 4, [1, 5, 5, 1], [1, '4.7429994', '3.8211503', '1.3794082'],
        '17.0140173', [0], [(0, 'lower')])
    add(100, 7, (x[0] - 10) ** 2 + 5 * (x[1] - 12) ** 2 + x[2] ** 4 + 3 * (x[3] - 11) ** 2
        + 10 * x[4] ** 6 + 7 * x[5] ** 2 + x[6] ** 4 - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6], [],
        [127 - 2 * x[0] ** 2 - 3 * x[1] ** 4 - x[2] - 4 * x[3] ** 2 - 5 * x[4],
         282 - 7 * x[0] - 3 * x[1] - 10 * x[2] ** 2 - x[3] + x[4],
         196 - 23 * x[0] - x[1] ** 2 - 6 * x[5] ** 2 + 8 * x[6],
         -4 * x[0] ** 2 - x[1] ** 2 + 3 * x[0] * x[1] - 2 * x[2] ** 2 - 5 * x[5] + 11 * x[6]],
        [None] * 7, [None] * 7, [1, 2, 0, 4, 0, 1, 1],
        ['2.330499', '1.951372', '-0.4775414', '4.365726', '-0.6244870', '1.038131', '1.594227'],
        '680.6300573', [0, 3], [])
    return out


def mpf(v):
    return mp.mpf(str(sp.N(v, 60))) if isinstance(v, sp.Basic) else mp.mpf(v)


def solve(p):
    n = p['n']
    xs = sp.symbols('x1:11')[:n]
    f, eq, ineq = p['f'], p['eq'], p['ineq']
    active = [eq[i] for i in range(len(eq))] + [ineq[j] for j in p['active']]
    grad = lambda e: [sp.diff(e, v) for v in xs]
    gf, ga = grad(f), [grad(c) for c in active]
    nb = len(p['bounds'])
    na = len(active)
    # Unknowns: x, then each active constraint's multiplier, then each
    # active bound's.
    def residual(*z):
        xv = dict(zip(xs, z[:n]))
        lam = z[n:n + na]
        mu = z[n + na:]
        out = []
        for i in range(n):
            r = mpf(gf[i].subs(xv))
            for k in range(na):
                r -= lam[k] * mpf(ga[k][i].subs(xv))
            for k, (var, side) in enumerate(p['bounds']):
                if var == i:
                    r += -mu[k] if side == 'lower' else mu[k]
            out.append(r)
        for c in active:
            out.append(mpf(c.subs(xv)))
        for (var, side) in p['bounds']:
            bound = p['lower'][var] if side == 'lower' else p['upper'][var]
            out.append(z[var] - mpf(bound))
        return out

    x_pub = [mpf(sp.sympify(v)) for v in p['x_star']]
    # The multipliers' first values: least squares at the published x*.
    xv = dict(zip(xs, x_pub))
    cols = [[mpf(ga[k][i].subs(xv)) for i in range(n)] for k in range(na)]
    for (var, side) in p['bounds']:
        cols.append([(mp.mpf(1) if side == 'lower' else mp.mpf(-1)) if i == var else mp.mpf(0)
                     for i in range(n)])
    if cols:
        a = mp.matrix(n, len(cols))
        for k, col in enumerate(cols):
            for i in range(n):
                a[i, k] = col[i]
        g = mp.matrix([mpf(gf[i].subs(xv)) for i in range(n)])
        m0 = list(mp.lu_solve(a.T * a, a.T * g))
    else:
        m0 = []
    z = mp.findroot(residual, x_pub + m0, tol=mp.mpf(10) ** -35, maxsteps=200)
    z = [z[i] for i in range(len(z))] if n + na + nb > 1 else [z]
    x_opt, lam, mu = z[:n], z[n:n + na], z[n + na:]
    xv = dict(zip(xs, x_opt))
    f_opt = mpf(f.subs(xv))
    # Checks: the published optimum, then feasibility and the multipliers.
    scale = lambda v: max(1, abs(v))
    for i in range(n):
        assert abs(x_opt[i] - x_pub[i]) <= 1e-6 * scale(x_pub[i]), (p['number'], i, x_opt[i])
    f_pub = mpf(sp.sympify(p['f_star']))
    assert abs(f_opt - f_pub) <= 1e-6 * scale(f_pub), (p['number'], f_opt, f_pub)
    for i, c in enumerate(eq):
        assert abs(mpf(c.subs(xv))) < mp.mpf(10) ** -30
    lam_ineq = [mp.mpf(0)] * len(ineq)
    for k, j in enumerate(p['active']):
        lam_ineq[j] = lam[len(eq) + k]
        assert lam_ineq[j] > 0, (p['number'], 'multiplier', j, lam_ineq[j])
    for j, c in enumerate(ineq):
        if j not in p['active']:
            assert mpf(c.subs(xv)) > 1e-8, (p['number'], 'inactive', j)
    mu_lower, mu_upper = [mp.mpf(0)] * n, [mp.mpf(0)] * n
    for k, (var, side) in enumerate(p['bounds']):
        assert mu[k] > 0, (p['number'], 'bound multiplier', var, mu[k])
        (mu_lower if side == 'lower' else mu_upper)[var] = mu[k]
    for i in range(n):
        if p['lower'][i] is not None and (i, 'lower') not in p['bounds']:
            assert x_opt[i] - mpf(p['lower'][i]) > 1e-8, (p['number'], 'lower', i)
        if p['upper'][i] is not None and (i, 'upper') not in p['bounds']:
            assert mpf(p['upper'][i]) - x_opt[i] > 1e-8, (p['number'], 'upper', i)
    num = lambda v: None if v is None else float(mpf(v))
    return {
        'number': p['number'],
        'n': n,
        'equalities': len(eq),
        'inequalities': len(ineq),
        'f': sp.sstr(f),
        'c': [sp.sstr(c) for c in eq + ineq],
        'lower': [num(v) for v in p['lower']],
        'upper': [num(v) for v in p['upper']],
        'x0': [float(mpf(sp.sympify(v))) for v in p['x0']],
        'published_x': [str(v) for v in p['x_star']],
        'published_f': str(p['f_star']),
        'x': [float(v) for v in x_opt],
        'f_value': float(f_opt),
        'multipliers': [float(v) for v in lam[:len(eq)]] + [float(v) for v in lam_ineq],
        'lower_multipliers': [float(v) for v in mu_lower],
        'upper_multipliers': [float(v) for v in mu_upper],
    }


def main():
    cases = [solve(p) for p in problems()]
    out = {
        'provenance': 'tools/gen_m8_hock_schittkowski.py: Hock and Schittkowski (1981), each '
        'problem in SymPy 1.13.3, its optimum solved again by mpmath 1.3.0 at 40 digits on the '
        "published optimum's active set, within 1e-6 of the published x* and f*, every multiplier "
        'of an active inequality or bound positive; c lists the equalities, then the '
        'inequalities, c >= 0',
        'problems': cases,
    }
    path = Path(__file__).resolve().parent.parent / 'testdata' / 'm8_hock_schittkowski.json'
    path.write_bytes((json.dumps(out, indent=1) + '\n').encode())
    print(f'wrote {len(cases)} problems to {path}')


if __name__ == '__main__':
    main()
