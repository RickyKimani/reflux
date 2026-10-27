module numeric

import math

fn rosenbrock(x []f64, mut g []f64, mut h []f64) !f64 {
	a, b := x[0], x[1]
	g[0] = -400.0 * a * (b - a * a) - 2.0 * (1.0 - a)
	g[1] = 200.0 * (b - a * a)
	h[0] = 1200.0 * a * a - 400.0 * b + 2.0
	h[1] = -400.0 * a
	h[2] = h[1]
	h[3] = 200.0
	return 100.0 * (b - a * a) * (b - a * a) + (1.0 - a) * (1.0 - a)
}

// powell_singular has its minimum 0 at the origin, where the Hessian is
// singular: the quartic terms make the minimum degenerate, as the trivial
// solution of a stability test is next to a critical point.
fn powell_singular(x []f64, mut g []f64, mut h []f64) !f64 {
	t1, t2, t3, t4 := x[0] + 10.0 * x[1], x[2] - x[3], x[1] - 2.0 * x[2], x[0] - x[3]
	for i in 0 .. 16 {
		h[i] = 0.0
	}
	g[0] = 2.0 * t1 + 40.0 * t4 * t4 * t4
	g[1] = 20.0 * t1 + 4.0 * t3 * t3 * t3
	g[2] = 10.0 * t2 - 8.0 * t3 * t3 * t3
	g[3] = -10.0 * t2 - 40.0 * t4 * t4 * t4
	q3, q4 := 12.0 * t3 * t3, 120.0 * t4 * t4
	h[0] = 2.0 + q4
	h[1], h[4] = 20.0, 20.0
	h[3], h[12] = -q4, -q4
	h[5] = 200.0 + q3
	h[6], h[9] = -2.0 * q3, -2.0 * q3
	h[10] = 10.0 + 4.0 * q3
	h[11], h[14] = -10.0, -10.0
	h[15] = 10.0 + q4
	return t1 * t1 + 5.0 * t2 * t2 + t3 * t3 * t3 * t3 + 10.0 * t4 * t4 * t4 * t4
}

// wood has its minimum 0 at (1, 1, 1, 1), and a saddle near
// (−0.97, 0.95, −0.97, 0.95) that a gradient-zero search can stop at.
fn wood(x []f64, mut g []f64, mut h []f64) !f64 {
	a, b, c, d := x[0], x[1], x[2], x[3]
	for i in 0 .. 16 {
		h[i] = 0.0
	}
	g[0] = -400.0 * a * (b - a * a) - 2.0 * (1.0 - a)
	g[1] = 200.0 * (b - a * a) + 20.2 * (b - 1.0) + 19.8 * (d - 1.0)
	g[2] = -360.0 * c * (d - c * c) - 2.0 * (1.0 - c)
	g[3] = 180.0 * (d - c * c) + 20.2 * (d - 1.0) + 19.8 * (b - 1.0)
	h[0] = 1200.0 * a * a - 400.0 * b + 2.0
	h[1], h[4] = -400.0 * a, -400.0 * a
	h[5] = 220.2
	h[7], h[13] = 19.8, 19.8
	h[10] = 1080.0 * c * c - 360.0 * d + 2.0
	h[11], h[14] = -360.0 * c, -360.0 * c
	h[15] = 200.2
	return 100.0 * (b - a * a) * (b - a * a) + (1.0 - a) * (1.0 - a) + 90.0 * (d - c * c) * (d - c * c) +
		(1.0 - c) * (1.0 - c) + 10.1 * ((b - 1.0) * (b - 1.0) + (d - 1.0) * (d - 1.0)) +
		19.8 * (b - 1.0) * (d - 1.0)
}

fn check_derivatives(name string, f ObjectiveFn, x []f64) ! {
	// The tests' own gradients and Hessians, against central differences.
	n := x.len
	mut g := []f64{len: n}
	mut h := []f64{len: n * n}
	f(x, mut g, mut h)!
	for k in 0 .. n {
		step := 1e-6 * math.max(1.0, math.abs(x[k]))
		mut xp := x.clone()
		mut xm := x.clone()
		xp[k] += step
		xm[k] -= step
		mut gp := []f64{len: n}
		mut gm := []f64{len: n}
		mut scratch := []f64{len: n * n}
		fp := f(xp, mut gp, mut scratch)!
		fm := f(xm, mut gm, mut scratch)!
		assert math.abs((fp - fm) / (2.0 * step) - g[k]) <= 1e-5 * math.max(1.0, math.abs(g[k])), '${name}: g[${k}]'
		for i in 0 .. n {
			fd := (gp[i] - gm[i]) / (2.0 * step)
			assert math.abs(fd - h[i * n + k]) <= 1e-5 * math.max(1.0, math.abs(fd)), '${name}: H[${i}][${k}]'
		}
	}
}

fn test_the_test_functions_derivatives() {
	check_derivatives('Rosenbrock', rosenbrock, [-1.2, 1.0])!
	check_derivatives('Powell singular', powell_singular, [3.0, -1.0, 0.0, 1.0])!
	check_derivatives('Wood', wood, [-3.0, -1.0, -3.0, -1.0])!
}

fn test_minimize_finds_the_standard_minima() {
	r := minimize(rosenbrock, [-1.2, 1.0])!
	assert math.abs(r.x[0] - 1.0) <= 1e-8 && math.abs(r.x[1] - 1.0) <= 1e-8 && r.minimum
	// From (0, 1) the Hessian is indefinite, so the first step needs a shift.
	r2 := minimize(rosenbrock, [0.0, 1.0])!
	assert math.abs(r2.x[0] - 1.0) <= 1e-8 && r2.minimum
	// Wood's function, past its saddle, to the true minimum.
	w := minimize(wood, [-3.0, -1.0, -3.0, -1.0])!
	for v in w.x {
		assert math.abs(v - 1.0) <= 1e-8, '${w.x}'
	}
	assert w.minimum && w.f <= 1e-16
	// The degenerate minimum: Newton slows to linear convergence there, but
	// gets to it, and the Hessian at the end is singular.
	p := minimize(powell_singular, [3.0, -1.0, 0.0, 1.0])!
	assert p.f <= 1e-12 && max_abs(p.g) <= 1e-10
	assert max_abs(p.x) <= 1e-2
}

fn test_minimize_takes_a_quadratic_in_one_step() {
	// f = ½xᵀAx − bᵀx with A positive definite: one Newton step, then the
	// gradient is zero to rounding.
	n := 6
	b := []f64{len: n, init: f64(index + 1)}
	f := fn [n, b] (x []f64, mut g []f64, mut h []f64) !f64 {
		mut v := 0.0
		for i in 0 .. n {
			mut ax := 0.0
			for j in 0 .. n {
				a := if i == j { 4.0 } else { 1.0 / f64(1 + i + j) }
				h[i * n + j] = a
				ax += a * x[j]
			}
			g[i] = ax - b[i]
			v += 0.5 * x[i] * ax - b[i] * x[i]
		}
		return v
	}
	s := minimize(f, []f64{len: n})!
	assert s.iterations == 1 && max_abs(s.g) <= 1e-12 && s.minimum
}

fn test_minimize_finishes_on_the_gradient_where_f_is_rounding() {
	// f = (10⁶ + q(x)) − 10⁶ with q = Σ (xᵢ − 1)² + (xᵢ − 1)⁴: f carries a
	// rounding of about 10⁶·ε ≈ 2e-10, so once q is that small no step can be
	// seen to lower it. The gradient is exact, though, and the last Newton
	// steps are judged by it, down to gtol.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		mut q := 0.0
		for i, v in x {
			d := v - 1.0
			q += d * d + d * d * d * d
			g[i] = 2.0 * d + 4.0 * d * d * d
			for j in 0 .. x.len {
				h[i * x.len + j] = if i == j { 2.0 + 12.0 * d * d } else { 0.0 }
			}
		}
		return (1e6 + q) - 1e6
	}
	s := minimize(f, [3.0, -2.0, 0.5], gtol: 1e-14, fscale: 1e6)!
	assert max_abs(s.g) <= 1e-14
	for v in s.x {
		assert math.abs(v - 1.0) <= 1e-14
	}
}

fn test_minimize_rejects_a_step_that_makes_f_worse() {
	// f = √(1 + x²) is convex, but Newton's step from x is to −x³, so from
	// x = 2 plain Newton runs off to infinity. Each such step raises f and
	// has to be cut back by the shift.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		r := math.sqrt(1.0 + x[0] * x[0])
		g[0] = x[0] / r
		h[0] = 1.0 / (r * r * r)
		return r
	}
	s := minimize(f, [2.0])!
	assert math.abs(s.x[0]) <= 1e-9 && s.minimum
}

fn test_minimize_stops_honestly_when_the_gradient_is_noise() {
	// The gradient carries noise of 1e-10 that never vanishes, as a computed
	// gradient's rounding wouldn't, so gtol = 1e-14 is out of reach. Once
	// neither f nor a Newton step can show progress, minimize returns the
	// point it has, close to 1, rather than wander until it runs out of
	// steps.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		d := x[0] - 1.0
		g[0] = 2.0 * d + 1e-10 * math.sin(1e12 * x[0])
		h[0] = 2.0
		return (1e6 + d * d) - 1e6
	}
	s := minimize(f, [3.0], gtol: 1e-14, fscale: 1e6)!
	assert math.abs(s.x[0] - 1.0) <= 1e-9 && max_abs(s.g) <= 1e-9
}

fn test_minimize_never_claims_a_minimum_against_a_wall() {
	// (x − 10)² where only x ≤ 0 can be evaluated: the lowest point there
	// is the wall, x = 0, where g = −20. Every step past it is refused, the
	// shift grows and the steps shrink, but a short step forced by the
	// shift says nothing about convergence. It must fail.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		if x[0] > 0.0 {
			return error('outside the domain')
		}
		g[0] = 2.0 * (x[0] - 10.0)
		h[0] = 2.0
		return (x[0] - 10.0) * (x[0] - 10.0)
	}
	if s := minimize(f, [-1.0]) {
		assert false, 'claimed a minimum at ${s.x} with g = ${s.g}'
	} else {
		assert err is NewtonError
	}
}

fn test_minimize_steps_back_into_the_domain() {
	// f = x − ln x has its minimum at 1. From 5 the Newton step overshoots
	// to x < 0, where f is undefined.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		if x[0] <= 0.0 {
			return error('x must be positive')
		}
		g[0] = 1.0 - 1.0 / x[0]
		h[0] = 1.0 / (x[0] * x[0])
		return x[0] - math.log(x[0])
	}
	s := minimize(f, [5.0])!
	assert math.abs(s.x[0] - 1.0) <= 1e-9
}

fn test_minimize_reports_a_saddle_for_what_it_is() {
	// f = x² − y² + y⁴ from (1, 0): the gradient never has a y component,
	// so the iteration ends at the saddle (0, 0). It says so: minimum is
	// false.
	f := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		g[0] = 2.0 * x[0]
		g[1] = -2.0 * x[1] + 4.0 * x[1] * x[1] * x[1]
		h[0], h[1], h[2] = 2.0, 0.0, 0.0
		h[3] = -2.0 + 12.0 * x[1] * x[1]
		return x[0] * x[0] - x[1] * x[1] + x[1] * x[1] * x[1] * x[1]
	}
	s := minimize(f, [1.0, 0.0])!
	assert max_abs(s.x) <= 1e-9 && !s.minimum
	// Off the axis the shift steers away from the saddle, to a minimum at
	// y = ±1/√2.
	t := minimize(f, [1.0, 0.1])!
	assert math.abs(math.abs(t.x[1]) - math.sqrt(0.5)) <= 1e-9 && t.minimum
}

fn test_minimize_fails_honestly() {
	// Unbounded below: f = −x² has no minimum. Never a success.
	unbounded := fn (x []f64, mut g []f64, mut h []f64) !f64 {
		g[0] = -2.0 * x[0]
		h[0] = -2.0
		return -x[0] * x[0]
	}
	if s := minimize(unbounded, [1.0]) {
		assert false, 'claimed a minimum of −x² at ${s.x}'
	} else {
		assert err is NewtonError
	}
	if _ := minimize(rosenbrock, []f64{}) {
		assert false
	}
	if _ := minimize(rosenbrock, [math.nan(), 1.0]) {
		assert false
	}
	if _ := minimize(rosenbrock, [-1.2, 1.0], max_iter: 2) {
		assert false, 'Rosenbrock in two steps'
	} else {
		assert (err as NewtonError).reason == 'no convergence in 2 steps'
	}
}
