module numeric

import math

// Line fits a polynomial of degree n − 1 to points: linear in θ, so its
// minimum is the normal equations' solution, which the test finds by LU.
struct Line {
	t []f64
	y []f64
	n int
}

fn (mut p Line) eval(theta []f64, mut r []f64, mut j []f64) ! {
	for i, t in p.t {
		mut v := 0.0
		for k in 0 .. p.n {
			v += theta[k] * math.pow(t, k)
			j[i * p.n + k] = math.pow(t, k)
		}
		r[i] = v - p.y[i]
	}
}

// A polynomial fit against the normal equations, (AᵀA)·θ = Aᵀy and
// (AᵀA)⁻¹, by LU. The fit stops at the gradient test, |Jₖᵀr| ≤
// gtol·‖Jₖ‖·‖r‖, which leaves θ within Σₗ |(JᵀJ)⁻¹ₖₗ|·gtol·‖Jₗ‖·‖r‖ of the
// minimum: the bound each θₖ is held to, with 1e-12 for the LU's own
// rounding. The inverse is held to 1e-10.
fn test_a_linear_fit_against_the_normal_equations() {
	t := []f64{len: 15, init: 2.0 * f64(index) / 14.0}
	y := t.map(1.0 - 2.0 * it + 0.5 * it * it * it + 0.01 * math.sin(7.0 * it))
	mut p := Line{t, y, 4}
	s := least_squares(mut p, [0.0, 0, 0, 0], t.len)!
	mut ata := []f64{len: 16}
	mut aty := []f64{len: 4}
	for i, ti in t {
		for a in 0 .. 4 {
			aty[a] += math.pow(ti, a) * y[i]
			for b in 0 .. 4 {
				ata[a * 4 + b] += math.pow(ti, a) * math.pow(ti, b)
			}
		}
	}
	f := lu_factor(ata, 4)!
	want := f.solve(aty)!
	norms := column_norms(s.j, t.len, 4)
	rnorm := math.sqrt(s.ssr)
	for k in 0 .. 4 {
		mut bound := 1e-12 * math.max(1.0, math.abs(want[k]))
		for l in 0 .. 4 {
			bound += math.abs(s.inverse[k * 4 + l]) * 1e-10 * norms[l] * rnorm
		}
		assert math.abs(s.theta[k] - want[k]) <= bound, '${k}: ${s.theta[k]} against ${want[k]}, bound ${bound}'
	}
	for c in 0 .. 4 {
		mut e := []f64{len: 4}
		e[c] = 1.0
		col := f.solve(e)!
		for k in 0 .. 4 {
			assert math.abs(s.inverse[k * 4 + c] - col[k]) <= 1e-10 * math.max(1.0, math.abs(col[k])), 'inverse ${k}, ${c}'
		}
	}
	assert s.gradient <= 1e-10
	assert math.abs(s.ssr - dot(s.r, s.r)) == 0.0
}

// Decay fits a·exp(−b·t), nonlinear in b, to points with no noise, so its
// minimum is the exact a and b with nothing left over; b must stay above
// zero, and an evaluation below it fails.
struct Decay {
	t []f64
	y []f64
}

fn (mut p Decay) eval(theta []f64, mut r []f64, mut j []f64) ! {
	if theta[1] <= 0.0 {
		return error('b = ${theta[1]} is not above 0')
	}
	for i, t in p.t {
		e := math.exp(-theta[1] * t)
		r[i] = theta[0] * e - p.y[i]
		j[i * 2] = e
		j[i * 2 + 1] = -theta[0] * t * e
	}
}

// An exponential with exact data, from far off, through a step b refuses.
fn test_an_exact_exponential() {
	t := []f64{len: 12, init: 0.5 * f64(index)}
	y := t.map(3.0 * math.exp(-0.7 * it))
	mut p := Decay{t, y}
	s := least_squares(mut p, [1.0, 0.05], t.len)!
	assert math.abs(s.theta[0] - 3.0) <= 1e-12 && math.abs(s.theta[1] - 0.7) <= 1e-12, '${s.theta}'
	assert s.ssr <= 1e-24
	// Measured: 9 iterations, which can't grow.
	assert s.iterations <= 9, '${s.iterations}'
}

// Rosen is Rosenbrock's function as residuals, 10(x₂ − x₁²) and 1 − x₁,
// whose minimum, (1, 1), has none.
struct Rosen {}

fn (mut p Rosen) eval(x []f64, mut r []f64, mut j []f64) ! {
	r[0] = 10.0 * (x[1] - x[0] * x[0])
	r[1] = 1.0 - x[0]
	j[0], j[1], j[2], j[3] = -20.0 * x[0], 10.0, -1.0, 0.0
}

fn test_rosenbrocks_residuals() {
	mut p := Rosen{}
	s := least_squares(mut p, [-1.2, 1.0], 2)!
	assert math.abs(s.theta[0] - 1.0) <= 1e-12 && math.abs(s.theta[1] - 1.0) <= 1e-12, '${s.theta}'
	// Measured: 22 iterations, which can't grow.
	assert s.iterations <= 22, '${s.iterations}'
}

fn refused(mut p Residuals, theta0 []f64, m int, says string, params LeastSquaresParams) {
	if s := least_squares(mut p, theta0, m, params) {
		assert false, 'fitted ${s.theta}'
	} else {
		assert err is LeastSquaresError, err.msg()
		assert (err as LeastSquaresError).reason.contains(says), err.msg()
	}
}

// What least_squares refuses, each with its reason: a parameter the data
// don't determine, a line's two coefficients given the same column; fewer
// residuals than parameters; a gtol of 0; a start that isn't finite; one
// where the problem fails; and the budget spent.
fn test_refusals() {
	mut same := Line{[0.0, 1, 2, 3], [1.0, 2, 3, 5], 2}
	mut twin := Twin{same}
	refused(mut twin, [0.0, 0], 4, 'parameter 1 is not determined by the data')
	refused(mut same, [0.0, 0], 1, 'at least as many residuals (1) as parameters (2)')
	refused(mut same, [0.0, 0], 4, 'gtol', gtol: 0.0)
	refused(mut same, [math.nan(), 0], 4, 'not finite')
	mut decay := Decay{[0.0, 1, 2], [3.0, 1.5, 0.75]}
	refused(mut decay, [1.0, -1.0], 3, 'at the start: b = -1.0 is not above 0')
	refused(mut decay, [1.0, 0.01], 3, 'not met in 1 iterations', max_iter: 1)
}

// Twin is a Line whose second column repeats its first.
struct Twin {
mut:
	line Line
}

fn (mut p Twin) eval(theta []f64, mut r []f64, mut j []f64) ! {
	for i, _ in p.line.t {
		r[i] = theta[0] + theta[1] - p.line.y[i]
		j[i * 2] = 1.0
		j[i * 2 + 1] = 1.0
	}
}

// Near is a line whose second column is its first but for 1e-13 in one
// row: dependent to rounding, which the refusal must see, not only an
// exact zero on R's diagonal.
struct Near {}

fn (mut p Near) eval(theta []f64, mut r []f64, mut j []f64) ! {
	y := [1.0, 2, 3, 5]
	for i in 0 .. 4 {
		second := if i == 3 { 1.0 + 1e-13 } else { 1.0 }
		r[i] = theta[0] + second * theta[1] - y[i]
		j[i * 2] = 1.0
		j[i * 2 + 1] = second
	}
}

fn test_a_column_dependent_to_rounding() {
	mut p := Near{}
	refused(mut p, [0.0, 0], 4, 'parameter 1 is not determined by the data')
}
