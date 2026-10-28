module numeric

import math

// A closure wrapped here is never freed, which a test can afford (see
// functions.v).
// SystemFn is a system as a function value, and Sys wraps one as a System.
type SystemFn = fn (x []f64, mut f []f64, mut jac []f64) !

struct Sys {
	f SystemFn @[required]
}

fn (w Sys) eval(x []f64, mut f []f64, mut jac []f64) ! {
	w.f(x, mut f, mut jac)!
}

// Problem is one test system: its size, standard starting point, and the
// root where one is known in closed form.
struct Problem {
	name string
	x0   []f64
	root []f64 // empty when only the residual can be checked
	// How close a converged x must be to root. |F| ≤ 1e-10 puts a simple
	// root within about 1e-10/‖J‖.
	root_tol f64      = 1e-8
	f        SystemFn = unsafe { nil }
}

// mgh returns the Moré, Garbow & Hillstrom (1981) problems for systems of
// nonlinear equations, in the form MINPACK's hybrj tests use, each with its
// analytic Jacobian and standard starting point. Watson's and Chebyquad
// are left out.
fn mgh() []Problem {
	n10 := 10
	h := 1.0 / f64(n10 + 1)
	return [
		Problem{
			name: 'Rosenbrock'
			x0:   [-1.2, 1.0]
			root: [1.0, 1.0]
			f:    fn (x []f64, mut f []f64, mut j []f64) ! {
				f[0] = 1.0 - x[0]
				f[1] = 10.0 * (x[1] - x[0] * x[0])
				j[0], j[1] = -1.0, 0.0
				j[2], j[3] = -20.0 * x[0], 10.0
			}
		},
		Problem{
			name:     'Powell singular'
			x0:       [3.0, -1.0, 0.0, 1.0]
			root:     [0.0, 0.0, 0.0, 0.0]
			// J is singular at the root and F is quadratic there, so |F| ≤ 1e-10
			// only puts x within about √1e-10.
			root_tol: 1e-4
			f:        fn (x []f64, mut f []f64, mut j []f64) ! {
				s5, s10 := math.sqrt(5.0), math.sqrt(10.0)
				f[0] = x[0] + 10.0 * x[1]
				f[1] = s5 * (x[2] - x[3])
				f[2] = (x[1] - 2.0 * x[2]) * (x[1] - 2.0 * x[2])
				f[3] = s10 * (x[0] - x[3]) * (x[0] - x[3])
				for i in 0 .. 16 {
					j[i] = 0.0
				}
				j[0], j[1] = 1.0, 10.0
				j[6], j[7] = s5, -s5
				j[9], j[10] = 2.0 * (x[1] - 2.0 * x[2]), -4.0 * (x[1] - 2.0 * x[2])
				j[12], j[15] = 2.0 * s10 * (x[0] - x[3]), -2.0 * s10 * (x[0] - x[3])
			}
		},
		Problem{
			name: 'Powell badly scaled'
			x0:   [0.0, 1.0]
			f:    fn (x []f64, mut f []f64, mut j []f64) ! {
				f[0] = 1e4 * x[0] * x[1] - 1.0
				f[1] = math.exp(-x[0]) + math.exp(-x[1]) - 1.0001
				j[0], j[1] = 1e4 * x[1], 1e4 * x[0]
				j[2], j[3] = -math.exp(-x[0]), -math.exp(-x[1])
			}
		},
		Problem{
			// The gradient of Wood's function. It has a second zero besides
			// (1, 1, 1, 1), the function's saddle near (−0.968, 0.947, −0.970,
			// 0.951), and Newton may find either, so only the residual is checked.
			name: 'Wood'
			x0:   [-3.0, -1.0, -3.0, -1.0]
			f:    fn (x []f64, mut f []f64, mut j []f64) ! {
				t1 := x[1] - x[0] * x[0]
				t2 := x[3] - x[2] * x[2]
				f[0] = -200.0 * x[0] * t1 - (1.0 - x[0])
				f[1] = 200.0 * t1 + 20.2 * (x[1] - 1.0) + 19.8 * (x[3] - 1.0)
				f[2] = -180.0 * x[2] * t2 - (1.0 - x[2])
				f[3] = 180.0 * t2 + 20.2 * (x[3] - 1.0) + 19.8 * (x[1] - 1.0)
				for i in 0 .. 16 {
					j[i] = 0.0
				}
				j[0], j[1] = -200.0 * (t1 - 2.0 * x[0] * x[0]) + 1.0, -200.0 * x[0]
				j[4], j[5], j[7] = -400.0 * x[0], 220.2, 19.8
				j[10], j[11] = -180.0 * (t2 - 2.0 * x[2] * x[2]) + 1.0, -180.0 * x[2]
				j[13], j[14], j[15] = 19.8, -360.0 * x[2], 200.2
			}
		},
		Problem{
			name: 'helical valley'
			x0:   [-1.0, 0.0, 0.0]
			root: [1.0, 0.0, 0.0]
			f:    fn (x []f64, mut f []f64, mut j []f64) ! {
				tpi := 2.0 * math.pi
				mut theta := math.atan(x[1] / x[0]) / tpi
				if x[0] < 0.0 {
					theta += 0.5
				}
				r2 := x[0] * x[0] + x[1] * x[1]
				r := math.sqrt(r2)
				f[0] = 10.0 * (x[2] - 10.0 * theta)
				f[1] = 10.0 * (r - 1.0)
				f[2] = x[2]
				j[0], j[1], j[2] = 100.0 * x[1] / (tpi * r2), -100.0 * x[0] / (tpi * r2), 10.0
				j[3], j[4], j[5] = 10.0 * x[0] / r, 10.0 * x[1] / r, 0.0
				j[6], j[7], j[8] = 0.0, 0.0, 1.0
			}
		},
		Problem{
			name: 'Brown almost-linear'
			x0:   []f64{len: n10, init: 0.5}
			f:    fn [n10] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				s := sum(x)
				mut prod := 1.0
				for v in x {
					prod *= v
				}
				for i in 0 .. n - 1 {
					f[i] = x[i] + s - f64(n + 1)
					for k in 0 .. n {
						j[i * n + k] = if k == i { 2.0 } else { 1.0 }
					}
				}
				f[n - 1] = prod - 1.0
				for k in 0 .. n {
					mut p := 1.0
					for m in 0 .. n {
						if m != k {
							p *= x[m]
						}
					}
					j[(n - 1) * n + k] = p
				}
			}
		},
		Problem{
			name: 'discrete boundary value'
			x0:   []f64{len: n10, init: f64(index + 1) * h * (f64(index + 1) * h - 1.0)}
			f:    fn [n10, h] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				for i in 0 .. n * n {
					j[i] = 0.0
				}
				for i in 0 .. n {
					t := f64(i + 1) * h
					c := x[i] + t + 1.0
					left := if i > 0 { x[i - 1] } else { 0.0 }
					right := if i < n - 1 { x[i + 1] } else { 0.0 }
					f[i] = 2.0 * x[i] - left - right + h * h * c * c * c / 2.0
					j[i * n + i] = 2.0 + 1.5 * h * h * c * c
					if i > 0 {
						j[i * n + i - 1] = -1.0
					}
					if i < n - 1 {
						j[i * n + i + 1] = -1.0
					}
				}
			}
		},
		Problem{
			name: 'discrete integral equation'
			x0:   []f64{len: n10, init: f64(index + 1) * h * (f64(index + 1) * h - 1.0)}
			f:    fn [n10, h] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				for i in 0 .. n {
					ti := f64(i + 1) * h
					mut lower, mut upper := 0.0, 0.0
					for k in 0 .. n {
						tk := f64(k + 1) * h
						c := x[k] + tk + 1.0
						if k <= i {
							lower += tk * c * c * c
							j[i * n + k] = 1.5 * h * (1.0 - ti) * tk * c * c
						} else {
							upper += (1.0 - tk) * c * c * c
							j[i * n + k] = 1.5 * h * ti * (1.0 - tk) * c * c
						}
					}
					f[i] = x[i] + h * ((1.0 - ti) * lower + ti * upper) / 2.0
					j[i * n + i] += 1.0
				}
			}
		},
		Problem{
			name: 'trigonometric'
			x0:   []f64{len: n10, init: 1.0 / f64(n10)}
			f:    fn [n10] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				mut c := 0.0
				for v in x {
					c += math.cos(v)
				}
				for i in 0 .. n {
					f[i] = f64(n) - c + f64(i + 1) * (1.0 - math.cos(x[i])) - math.sin(x[i])
					for k in 0 .. n {
						j[i * n + k] = math.sin(x[k])
					}
					j[i * n + i] += f64(i + 1) * math.sin(x[i]) - math.cos(x[i])
				}
			}
		},
		Problem{
			name: 'variably dimensioned'
			x0:   []f64{len: n10, init: 1.0 - f64(index + 1) / f64(n10)}
			root: []f64{len: n10, init: 1.0}
			f:    fn [n10] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				mut s := 0.0
				for k in 0 .. n {
					s += f64(k + 1) * (x[k] - 1.0)
				}
				for i in 0 .. n {
					f[i] = x[i] - 1.0 + f64(i + 1) * s * (1.0 + 2.0 * s * s)
					for k in 0 .. n {
						j[i * n + k] = f64(i + 1) * f64(k + 1) * (1.0 + 6.0 * s * s)
					}
					j[i * n + i] += 1.0
				}
			}
		},
		Problem{
			name: 'Broyden tridiagonal'
			x0:   []f64{len: n10, init: -1.0}
			f:    fn [n10] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				for i in 0 .. n * n {
					j[i] = 0.0
				}
				for i in 0 .. n {
					left := if i > 0 { x[i - 1] } else { 0.0 }
					right := if i < n - 1 { x[i + 1] } else { 0.0 }
					f[i] = (3.0 - 2.0 * x[i]) * x[i] - left - 2.0 * right + 1.0
					j[i * n + i] = 3.0 - 4.0 * x[i]
					if i > 0 {
						j[i * n + i - 1] = -1.0
					}
					if i < n - 1 {
						j[i * n + i + 1] = -2.0
					}
				}
			}
		},
		Problem{
			name: 'Broyden banded'
			x0:   []f64{len: n10, init: -1.0}
			f:    fn [n10] (x []f64, mut f []f64, mut j []f64) ! {
				n := n10
				for i in 0 .. n {
					f[i] = x[i] * (2.0 + 5.0 * x[i] * x[i]) + 1.0
					for k in 0 .. n {
						j[i * n + k] = 0.0
					}
					j[i * n + i] = 2.0 + 15.0 * x[i] * x[i]
					for k in math.max(0, i - 5) .. math.min(n, i + 2) {
						if k != i {
							f[i] -= x[k] * (1.0 + x[k])
							j[i * n + k] = -(1.0 + 2.0 * x[k])
						}
					}
				}
			}
		},
	]
}

// freudenstein_roth has a root at (5, 4), but from the standard start
// (0.5, −2) methods without luck run into a local minimum of ‖F‖ near
// (11.41, −0.8968), where F ≠ 0 and J is singular.
fn freudenstein_roth() Problem {
	return Problem{
		name: 'Freudenstein–Roth'
		x0:   [0.5, -2.0]
		root: [5.0, 4.0]
		f:    fn (x []f64, mut f []f64, mut j []f64) ! {
			f[0] = -13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1]
			f[1] = -29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1]
			j[0], j[1] = 1.0, 10.0 * x[1] - 3.0 * x[1] * x[1] - 2.0
			j[2], j[3] = 1.0, 3.0 * x[1] * x[1] + 2.0 * x[1] - 14.0
		}
	}
}

// residual_at evaluates p's F at x afresh, independently of the solver.
fn residual_at(p Problem, x []f64) f64 {
	n := x.len
	mut f := []f64{len: n}
	mut j := []f64{len: n * n}
	p.f(x, mut f, mut j) or { return math.inf(1) }
	return max_abs(f)
}

fn test_the_test_jacobians_match_finite_differences() {
	// A wrong analytic Jacobian in a test problem would look like a solver
	// bug, so check each against central differences first, at x0 and at a
	// second, less symmetric point.
	mut problems := mgh()
	problems << freudenstein_roth()
	for p in problems {
		n := p.x0.len
		for shift in [0.0, 0.137] {
			x := []f64{len: n, init: p.x0[index] + shift * f64(index + 1)}
			mut f := []f64{len: n}
			mut j := []f64{len: n * n}
			p.f(x, mut f, mut j)!
			for k in 0 .. n {
				step := 1e-6 * math.max(math.abs(x[k]), 1.0)
				mut xp := x.clone()
				mut xm := x.clone()
				xp[k] += step
				xm[k] -= step
				mut fp := []f64{len: n}
				mut fm := []f64{len: n}
				mut scratch := []f64{len: n * n}
				p.f(xp, mut fp, mut scratch)!
				p.f(xm, mut fm, mut scratch)!
				for i in 0 .. n {
					fd := (fp[i] - fm[i]) / (2.0 * step)
					assert math.abs(fd - j[i * n + k]) <= 1e-5 * math.max(math.abs(fd), 1.0), '${p.name}: ∂F${i}/∂x${k} is ${j[i * n + k]}, finite differences give ${fd}'
				}
			}
		}
	}
}

fn test_newton_solves_the_mgh_problems_or_says_it_did_not() {
	// Each problem from 1× and 10× its standard start (MINPACK's convention:
	// a zero in x0 becomes the factor), with both methods. A run must either
	// converge, which is checked here by evaluating F afresh, or fail with a
	// NewtonError: never claim a point that isn't a root. The runs that fail
	// are recorded, and the list can't grow. Both failures stop at a local
	// minimum of ‖F‖ that isn't a root (see the next test); MINPACK's hybrj
	// fails the trigonometric problem from 10×x0 too.
	expected_failures := ['trigonometric 10.0x0 line_search', 'trigonometric 10.0x0 dogleg']
	mut failures := []string{}
	for p in mgh() {
		for factor in [1.0, 10.0] {
			x0 := p.x0.map(if it == 0.0 && factor != 1.0 { factor } else { it * factor })
			for method in [Globalization.line_search, .dogleg] {
				label := '${p.name} ${factor}x0 ${method}'
				sol := newton(Sys{p.f}, x0, method: method) or {
					assert err is NewtonError, '${label}: ${err}'
					failures << label
					continue
				}
				assert sol.residual <= 1e-10, label
				assert residual_at(p, sol.x) == sol.residual, label
				if p.root.len > 0 {
					for i, v in p.root {
						assert math.abs(sol.x[i] - v) <= p.root_tol, '${label}: x[${i}] = ${sol.x[i]}, root ${v}'
					}
				}
			}
		}
	}
	for f in failures {
		assert f in expected_failures, 'a new failure: ${f}'
	}
}

fn test_a_reported_local_minimum_is_one() {
	// The dogleg's failure on the trigonometric problem says it found a
	// local minimum of ‖F‖ that isn't a root. Check that's true: at the
	// point it returns, the gradient JᵀF of ½‖F‖² vanishes while F doesn't.
	p := mgh().filter(it.name == 'trigonometric')[0]
	newton(Sys{p.f}, p.x0.map(it * 10.0), method: .dogleg) or {
		assert err is NewtonError
		e := err as NewtonError
		assert e.reason.contains('local minimum')
		n := e.x.len
		mut f := []f64{len: n}
		mut j := []f64{len: n * n}
		p.f(e.x, mut f, mut j)!
		assert max_abs(f) == e.residual
		mut g := []f64{len: n}
		for i in 0 .. n {
			for k in 0 .. n {
				g[k] += j[i * n + k] * f[i]
			}
		}
		assert e.residual > 1e-4
		assert max_abs(g) <= 1e-8
		return
	}
	assert false, 'expected the trigonometric problem from 10×x0 to fail'
}

fn test_freudenstein_roth_is_solved_or_refused() {
	// From the standard start, Newton heads for the local minimum near
	// (11.41, −0.8968). Whatever happens, a success must be the real root.
	p := freudenstein_roth()
	for method in [Globalization.line_search, .dogleg] {
		sol := newton(Sys{p.f}, p.x0, method: method) or {
			assert err is NewtonError, '${method}: ${err}'
			continue
		}
		assert residual_at(p, sol.x) <= 1e-10, '${method}'
		assert math.abs(sol.x[0] - 5.0) <= 1e-8 && math.abs(sol.x[1] - 4.0) <= 1e-8, '${method}: ${sol.x}'
	}
}

fn test_a_system_without_a_root_fails() {
	f := fn (x []f64, mut f []f64, mut j []f64) ! {
		f[0] = x[0] * x[0] + 1.0
		j[0] = 2.0 * x[0]
	}
	for method in [Globalization.line_search, .dogleg] {
		if sol := newton(Sys{f}, [3.0], method: method) {
			assert false, '${method} claimed a root of x² + 1 at ${sol.x}'
		} else {
			assert err is NewtonError, '${method}: ${err}'
			e := err as NewtonError
			// The best point is the minimum of x² + 1, at 0.
			assert math.abs(e.x[0]) <= 1e-3, '${method}: ${e.x}'
		}
	}
}

fn test_newton_steps_back_into_the_domain() {
	// ln x = 1 from x = 10. The first Newton step lands at x ≈ −3, where
	// the system refuses to evaluate, so the solver must take a shorter
	// step instead of failing.
	f := fn (x []f64, mut f []f64, mut j []f64) ! {
		if x[0] <= 0.0 {
			return error('x must be positive')
		}
		f[0] = math.log(x[0]) - 1.0
		j[0] = 1.0 / x[0]
	}
	for method in [Globalization.line_search, .dogleg] {
		sol := newton(Sys{f}, [10.0], method: method)!
		assert math.abs(sol.x[0] - math.e) <= 1e-9, '${method}: ${sol.x}'
		assert sol.residual <= 1e-10
	}
}

fn test_newton_gets_past_a_singular_jacobian() {
	// x² + y = 2, x = y, from (−0.5, 0): J = [[2x, 1], [1, −1]] is singular
	// there, so the first step is the regularised one. The roots are (1, 1)
	// and (−2, −2).
	f := fn (x []f64, mut f []f64, mut j []f64) ! {
		f[0] = x[0] * x[0] + x[1] - 2.0
		f[1] = x[0] - x[1]
		j[0], j[1] = 2.0 * x[0], 1.0
		j[2], j[3] = 1.0, -1.0
	}
	for method in [Globalization.line_search, .dogleg] {
		sol := newton(Sys{f}, [-0.5, 0.0], method: method)!
		assert sol.residual <= 1e-10
		assert math.abs(sol.x[0] - 1.0) <= 1e-9 || math.abs(sol.x[0] + 2.0) <= 1e-9, '${method}: ${sol.x}'
	}
}

fn test_newton_refuses_to_start_from_a_stationary_point() {
	// x² + y² = 4, x = y, from (0, 0): J = [[0, 0], [1, −1]] and F = (−4, 0),
	// so the gradient JᵀF of ‖F‖² is exactly zero. No method can move, and
	// the error has to say why.
	f := fn (x []f64, mut f []f64, mut j []f64) ! {
		f[0] = x[0] * x[0] + x[1] * x[1] - 4.0
		f[1] = x[0] - x[1]
		j[0], j[1] = 2.0 * x[0], 2.0 * x[1]
		j[2], j[3] = 1.0, -1.0
	}
	for method in [Globalization.line_search, .dogleg] {
		if _ := newton(Sys{f}, [0.0, 0.0], method: method) {
			assert false, '${method} moved from a stationary point'
		} else {
			e := err as NewtonError
			assert e.reason.contains('stationary point'), '${method}: ${e.reason}'
			assert e.iterations == 0 && e.residual == 4.0
		}
	}
}

fn test_newton_stops_on_xtol_when_ftol_is_out_of_reach() {
	// 1e9·(x² − 2): at the double nearest √2, F is 4.4e-7 of pure rounding,
	// so ftol = 1e-10 can't be met. The Newton step becomes shorter than
	// xtol instead, which is convergence, and the residual is reported as
	// it is.
	f := fn (x []f64, mut f []f64, mut j []f64) ! {
		f[0] = 1e9 * (x[0] * x[0] - 2.0)
		j[0] = 2e9 * x[0]
	}
	for method in [Globalization.line_search, .dogleg] {
		sol := newton(Sys{f}, [1.0], method: method)!
		assert math.abs(sol.x[0] - math.sqrt2) <= 4.0 * f64_epsilon, '${method}: ${sol.x}'
		assert sol.residual > 1e-10 && sol.residual < 1e-5, '${method}: ${sol.residual}'
	}
}

fn test_newton_counts_and_limits() {
	p := mgh()[0] // Rosenbrock
	// Already at the root: nothing to do, one evaluation.
	sol := newton(Sys{p.f}, [1.0, 1.0])!
	assert sol.iterations == 0 && sol.evaluations == 1 && sol.residual == 0.0
	// Out of steps: an error that says so, carrying the point reached.
	if _ := newton(Sys{mgh()[3].f}, mgh()[3].x0, max_iter: 2) {
		assert false, 'Wood in two steps'
	} else {
		e := err as NewtonError
		assert e.reason == 'no convergence in 2 steps'
		assert e.iterations == 2 && e.x.len == 4
	}
}

fn test_newton_refuses_bad_input() {
	p := mgh()[0]
	for x0 in [[]f64{}, [math.nan(), 1.0], [math.inf(1), 1.0]] {
		if _ := newton(Sys{p.f}, x0) {
			assert false, 'x0 = ${x0} should be refused'
		}
	}
	if _ := newton(Sys{p.f}, p.x0, ftol: 0.0) {
		assert false, 'ftol = 0 should be refused'
	}
	if _ := newton(Sys{p.f}, p.x0, max_iter: 0) {
		assert false, 'max_iter = 0 should be refused'
	}
	nan_at_start := fn (x []f64, mut f []f64, mut j []f64) ! {
		f[0] = math.nan()
		j[0] = 1.0
	}
	if _ := newton(Sys{nan_at_start}, [1.0]) {
		assert false, 'a system that is NaN at x0 should be refused'
	}
}
