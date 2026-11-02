// Broyden's method against what's known of it exactly: in one dimension it
// is the secant method, step for step; on a linear system it finishes within
// 2n steps (Gay, SIAM J. Numer. Anal. 16, 1979); and on nonlinear systems it
// reaches the root Newton's method finds.
import math
import numeric

// Cosine is cos x − x = 0, as a fixed point x = cos x.
struct Cosine {
mut:
	calls int
}

fn (mut c Cosine) eval(x []f64, mut f []f64) ! {
	c.calls++
	f[0] = math.cos(x[0]) - x[0]
}

fn test_one_dimension_is_the_secant_method() {
	// Its first step is x + F(x), and each one after is the secant step
	// through the last two points. Checked by writing the secant out.
	mut sys := Cosine{}
	r := numeric.broyden(mut sys, [0.5], ftol: 1e-14)!
	mut a := 0.5
	mut fa := math.cos(a) - a
	mut b := a + fa
	mut fb := math.cos(b) - b
	mut steps := 1
	for math.abs(fb) > 1e-14 {
		c := b - fb * (b - a) / (fb - fa)
		a, fa = b, fb
		b = c
		fb = math.cos(b) - b
		steps++
	}
	assert r.x[0] == b, '${r.x[0]} vs ${b}'
	assert r.iterations == steps
	assert r.evaluations == sys.calls
	assert r.evaluations == steps + 1
	// The Dottie number, to 20 figures, and its residual as returned.
	assert math.abs(r.x[0] - 0.73908513321516064166) <= 1e-14
	assert r.residual == math.abs(r.f[0])
	assert r.residual <= 1e-14
}

// Linear is A·x − b = 0.
struct Linear {
	a [][]f64
	b []f64
}

fn (mut l Linear) eval(x []f64, mut f []f64) ! {
	for i, row in l.a {
		f[i] = -l.b[i]
		for j, v in row {
			f[i] += v * x[j]
		}
	}
}

fn test_a_linear_system_finishes_within_2n_steps() {
	// Gay's theorem: exact in 2n steps, in exact arithmetic. Each A is a fixed
	// diagonally dominant matrix built from integers, so the root is known by
	// construction, x = 1, 2, 3, ….
	for n in [1, 2, 3, 5, 8, 13] {
		mut a := [][]f64{len: n, init: []f64{len: n}}
		for i in 0 .. n {
			for j in 0 .. n {
				a[i][j] = f64(((i + 1) * (j + 2)) % 5) / 10.0
			}
			a[i][i] = -3.0 - f64(i % 3)
		}
		root := []f64{len: n, init: f64(index + 1)}
		mut b := []f64{len: n}
		for i in 0 .. n {
			for j in 0 .. n {
				b[i] += a[i][j] * root[j]
			}
		}
		mut sys := Linear{a, b}
		r := numeric.broyden(mut sys, []f64{len: n}, ftol: 1e-10)!
		// One step more than Gay's count, for the residual to show it.
		assert r.iterations <= 2 * n + 1, 'n = ${n}: ${r.iterations} steps'
		for i, v in r.x {
			assert math.abs(v - root[i]) <= 1e-9, 'n = ${n}: x[${i}] = ${v}'
		}
	}
}

// Tridiagonal is Broyden's own tridiagonal function, from the 1965 paper:
// (3 − 2xᵢ)xᵢ − xᵢ₋₁ − 2xᵢ₊₁ + 1 = 0, with x₀ = xₙ₊₁ = 0.
struct Tridiagonal {}

fn (mut t Tridiagonal) eval(x []f64, mut f []f64) ! {
	n := x.len
	for i in 0 .. n {
		left := if i > 0 { x[i - 1] } else { 0.0 }
		right := if i < n - 1 { x[i + 1] } else { 0.0 }
		f[i] = (3.0 - 2.0 * x[i]) * x[i] - left - 2.0 * right + 1.0
	}
}

// Newtonian is the same function with its Jacobian, for numeric.newton.
struct Newtonian {}

fn (s Newtonian) eval(x []f64, mut f []f64, mut jac []f64) ! {
	n := x.len
	for i in 0 .. n {
		left := if i > 0 { x[i - 1] } else { 0.0 }
		right := if i < n - 1 { x[i + 1] } else { 0.0 }
		f[i] = (3.0 - 2.0 * x[i]) * x[i] - left - 2.0 * right + 1.0
		for j in 0 .. n {
			jac[i * n + j] = 0.0
		}
		jac[i * n + i] = 3.0 - 4.0 * x[i]
		if i > 0 {
			jac[i * n + i - 1] = -1.0
		}
		if i < n - 1 {
			jac[i * n + i + 1] = -2.0
		}
	}
}

fn test_the_tridiagonal_system_reaches_newtons_root() {
	// Broyden's own test, from his starting point, x = −1, with the Jacobian's
	// diagonal for the first guess. Newton's root is the reference.
	for n in [5, 10, 20] {
		x0 := []f64{len: n, init: -1.0}
		want := numeric.newton(Newtonian{}, x0, ftol: 1e-13)!
		mut sys := Tridiagonal{}
		r := numeric.broyden(mut sys, x0, ftol: 1e-11, scale: 6.0, max_iter: 200)!
		for i, v in r.x {
			assert math.abs(v - want.x[i]) <= 1e-9, 'n = ${n}: x[${i}]'
		}
		// Far fewer evaluations than Newton's Jacobians would take by differences.
		assert r.evaluations <= 60, '${r.evaluations}'
	}
}

// Fenced is 3(x − 1) = 0, with x below −0.5 out of its domain.
struct Fenced {
mut:
	refused int
}

fn (mut s Fenced) eval(x []f64, mut f []f64) ! {
	if x[0] < -0.5 {
		s.refused++
		return error('out of the domain')
	}
	f[0] = 3.0 * (x[0] - 1.0)
}

fn test_a_step_outside_the_domain_is_halved() {
	// With scale 1, the first step from 2 is 2 − 3 = −1, refused, and it's
	// halved to 0.5.
	mut s := Fenced{}
	r := numeric.broyden(mut s, [2.0], scale: 1.0)!
	assert math.abs(r.x[0] - 1.0) <= 1e-10
	assert s.refused == 1
}

// Nowhere can't be evaluated anywhere but where it starts.
struct Nowhere {}

fn (mut s Nowhere) eval(x []f64, mut f []f64) ! {
	if x[0] != 2.0 {
		return error('out of the domain')
	}
	f[0] = 1.0
}

// Bad returns a value that isn't finite away from the start.
struct Bad {}

fn (mut s Bad) eval(x []f64, mut f []f64) ! {
	f[0] = if x[0] == 2.0 { 1.0 } else { math.nan() }
}

// Unsolvable is 1 + x² = 0, with no real root. It keeps the smallest value it
// returned.
struct Unsolvable {
mut:
	least f64 = math.inf(1)
}

fn (mut s Unsolvable) eval(x []f64, mut f []f64) ! {
	f[0] = 1.0 + x[0] * x[0]
	s.least = math.min(s.least, f[0])
}

// Flat doesn't depend on x, so the update's denominator is zero.
struct Flat {}

fn (mut s Flat) eval(_ []f64, mut f []f64) ! {
	f[0] = 1.0
}

// Hopeless isn't finite even where it starts.
struct Hopeless {}

fn (mut s Hopeless) eval(_ []f64, mut f []f64) ! {
	f[0] = math.nan()
}

fn test_a_failure_is_an_error_with_the_best_point() {
	mut a := Nowhere{}
	numeric.broyden(mut a, [2.0]) or {
		assert err is numeric.BroydenError, err.msg()
		e := err as numeric.BroydenError
		assert e.reason.contains('no step'), err.msg()
		assert e.x == [2.0] && e.residual == 1.0
		assert err.msg().starts_with('broyden: no step'), err.msg()
		mut b := Bad{}
		numeric.broyden(mut b, [2.0]) or {
			assert err.msg().contains('no step'), err.msg()
			mut c := Unsolvable{}
			numeric.broyden(mut c, [2.0], max_iter: 5) or {
				assert err is numeric.BroydenError, err.msg()
				f := err as numeric.BroydenError
				assert f.reason == 'ran out of iterations' && f.iterations == 5, err.msg()
				// The best point really is the best it saw: 1 + x² there.
				assert math.abs(f.residual - (1.0 + f.x[0] * f.x[0])) <= 1e-12
				assert f.residual == c.least, '${f.residual} vs ${c.least}'
				// With nothing to learn from, H is kept, and the iterations just run out.
				mut d := Flat{}
				numeric.broyden(mut d, [2.0], max_iter: 5) or {
					assert err is numeric.BroydenError, err.msg()
					assert (err as numeric.BroydenError).reason == 'ran out of iterations', err.msg()
					return
				}
				assert false, 'a root of 1'
				return
			}
			assert false, 'a root of 1 + x²'
			return
		}
		assert false, 'a root from nothing but nan'
		return
	}
	assert false, 'a root from a point it could not leave'
}

// Stuck fails at its start.
struct Stuck {}

fn (mut s Stuck) eval(_ []f64, mut _ []f64) ! {
	return error('stuck')
}

fn test_a_bad_call_is_refused() {
	mut c := Cosine{}
	for message in ['x0 is empty', 'x0[0] = nan', 'must be positive'] {
		x := if message == 'x0 is empty' {
			[]f64{}
		} else if message == 'x0[0] = nan' {
			[math.nan()]
		} else {
			[1.0]
		}
		numeric.broyden(mut c, x, scale: if message == 'must be positive' { 0.0 } else { -1.0 }) or {
			assert err.msg().contains(message), err.msg()
			continue
		}
		assert false, message
	}
	// A start that isn't finite isn't a start.
	mut h := Hopeless{}
	numeric.broyden(mut h, [1.0]) or {
		assert err.msg().contains('not finite at x0'), err.msg()
		// A failure at x0 is the system's own.
		mut s := Stuck{}
		numeric.broyden(mut s, [1.0]) or {
			assert err.msg() == 'stuck', err.msg()
			return
		}
		assert false, 'a root of nothing'
		return
	}
	assert false, 'a root from nan'
}

fn test_the_last_iteration_can_converge() {
	// Allowed exactly the steps it needs, it succeeds.
	mut a := Cosine{}
	r := numeric.broyden(mut a, [0.5], ftol: 1e-14)!
	mut b := Cosine{}
	s := numeric.broyden(mut b, [0.5], ftol: 1e-14, max_iter: r.iterations)!
	assert s.x == r.x && s.iterations == r.iterations
}

// Plateau is flat up to 3, then 4 − x: its first steps teach it nothing.
struct Plateau {}

fn (mut s Plateau) eval(x []f64, mut f []f64) ! {
	f[0] = if x[0] < 3.0 { 1.0 } else { 4.0 - x[0] }
}

fn test_a_step_that_changes_nothing_keeps_h() {
	// On the plateau the update's denominator is zero, so H stays as it was and
	// the steps go on, off the end of it and to the root.
	mut s := Plateau{}
	r := numeric.broyden(mut s, [0.0])!
	assert math.abs(r.x[0] - 4.0) <= 1e-10, '${r.x[0]}'
}
