module fit

import math
import numeric

// The point test of a binary's VLE data for Gibbs–Duhem consistency (M6's
// criterion 8), as Fredenslund, Gmehling and Rasmussen give Van Ness, Byer
// and Gibbs' test (Vapor–Liquid Equilibria Using UNIFAC, 1977, ch. 3): Gᴱ
// fitted to the bubble pressures alone by Barker's method, then the vapour
// it predicts against the vapour measured. Data whose y and P (or T) agree
// with one Gᴱ are consistent; y that the P–x data can't give aren't.

// The point test's verdict: consistent when the mean |Δy| is at most this
// (Fredenslund et al.'s criterion).
pub const point_test_limit = 0.01

// The Legendre series' terms: Fredenslund et al. use up to five, and four
// fit every binary of criterion 8 to within its data's scatter.
pub const point_test_terms = 4

// PointTest is a data set's point test: Gᴱ/RT's Legendre coefficients, each
// point's Δy = y_model − y_measured and ΔP/P, and the verdict.
pub struct PointTest {
pub:
	coefficients []f64
	dy           []f64
	dp           []f64 // (P_model − P)/P at each point
	mean_dy      f64
	max_dy       f64
	consistent   bool
	iterations   int
}

// PointData is a binary's bubble points for the point test: each point's
// liquid and vapour mole fractions of the first component, its pressure
// (Pa), and each component's vapour pressure at its T (Pa).
pub struct PointData {
pub:
	x1    []f64
	y1    []f64
	p     []f64
	psat1 []f64
	psat2 []f64
	terms int = point_test_terms
}

// legendre is L₀ … L_{n−1} at z and their derivatives.
fn legendre(z f64, n int) ([]f64, []f64) {
	mut l := []f64{len: n}
	mut d := []f64{len: n}
	if n > 0 {
		l[0] = 1.0
	}
	if n > 1 {
		l[1] = z
		d[1] = 1.0
	}
	for k in 2 .. n {
		// Bonnet's recursion, and its derivative.
		kf := f64(k)
		l[k] = ((2.0 * kf - 1.0) * z * l[k - 1] - (kf - 1.0) * l[k - 2]) / kf
		d[k] = d[k - 2] + (2.0 * kf - 1.0) * l[k - 1]
	}
	return l, d
}

// ln_gamma_terms are each coefficient's part of ln γ₁ and ln γ₂ at x₁, for
// Gᴱ/RT = x₁x₂ Σₖ aₖ Lₖ(x₁ − x₂): with g = Gᴱ/RT, ln γ₁ = g + x₂ g′ and
// ln γ₂ = g − x₁ g′, g′ = dg/dx₁, all linear in the aₖ.
fn ln_gamma_terms(x1 f64, n int) ([]f64, []f64) {
	x2 := 1.0 - x1
	l, d := legendre(x1 - x2, n)
	mut c1 := []f64{len: n}
	mut c2 := []f64{len: n}
	for k in 0 .. n {
		g := x1 * x2 * l[k]
		dg := (x2 - x1) * l[k] + 2.0 * x1 * x2 * d[k]
		c1[k] = g + x2 * dg
		c2[k] = g - x1 * dg
	}
	return c1, c2
}

// PointProblem is Barker's method as least squares' residuals: each bubble
// point's (P_model − P)/P, the vapour the ideal gas, with its exact
// Jacobian.
struct PointProblem {
	d PointData
}

fn (mut p PointProblem) eval(theta []f64, mut r []f64, mut j []f64) ! {
	n := theta.len
	for i, x1 in p.d.x1 {
		c1, c2 := ln_gamma_terms(x1, n)
		mut l1, mut l2 := 0.0, 0.0
		for k in 0 .. n {
			l1 += theta[k] * c1[k]
			l2 += theta[k] * c2[k]
		}
		a := x1 * math.exp(l1) * p.d.psat1[i]
		b := (1.0 - x1) * math.exp(l2) * p.d.psat2[i]
		if !math.is_finite(a + b) {
			return error('the model pressure at x₁ = ${x1} is ${a + b}')
		}
		r[i] = (a + b - p.d.p[i]) / p.d.p[i]
		for k in 0 .. n {
			j[i * n + k] = (a * c1[k] + b * c2[k]) / p.d.p[i]
		}
	}
}

// LinearProblem is A·θ − b as least squares' residuals, A m×n row-major.
struct LinearProblem {
	a []f64
	b []f64
}

fn (mut p LinearProblem) eval(theta []f64, mut r []f64, mut j []f64) ! {
	n := theta.len
	for i in 0 .. p.b.len {
		mut v := -p.b[i]
		for k in 0 .. n {
			v += p.a[i * n + k] * theta[k]
			j[i * n + k] = p.a[i * n + k]
		}
		r[i] = v
	}
}

// measured_start is the Gᴱ the measured vapour implies, the modified
// Raoult's law's γᵢ = yᵢP/(xᵢPᵢˢᵃᵗ) at each point, its coefficients by
// linear least squares on gᵢ/(x₁x₂) = Σₖ aₖ Lₖ(x₁ − x₂): a start for
// Barker's fit, which then sees the pressures alone. None where a y is 0
// or 1.
fn measured_start(d PointData) ?[]f64 {
	m := d.x1.len
	mut a := []f64{len: m * d.terms}
	mut b := []f64{len: m}
	for i, x1 in d.x1 {
		x2 := 1.0 - x1
		if !(d.y1[i] > 0.0 && d.y1[i] < 1.0) {
			return none
		}
		g1 := math.log(d.y1[i] * d.p[i] / (x1 * d.psat1[i]))
		g2 := math.log((1.0 - d.y1[i]) * d.p[i] / (x2 * d.psat2[i]))
		b[i] = (x1 * g1 + x2 * g2) / (x1 * x2)
		l, _ := legendre(x1 - x2, d.terms)
		for k in 0 .. d.terms {
			a[i * d.terms + k] = l[k]
		}
	}
	mut lin := LinearProblem{a, b}
	s := numeric.least_squares(mut lin, []f64{len: d.terms}, m) or { return none }
	return s.theta
}

// point_test runs the point test on a binary's bubble points, each strictly
// inside (0, 1): Gᴱ/RT's Legendre coefficients fitted to the pressures, the
// vapour the ideal gas, then y₁ = x₁γ₁P₁ˢᵃᵗ/P_model at each point against
// the y₁ measured. For isobaric data each point is at its own T, which the
// vapour pressures carry; Gᴱ's T-dependence over the set is neglected, as
// in Fredenslund et al.'s test.
//
// Barker's fit can have more than one minimum where the data cover part of
// the range (ethanol–water's azeotropic region alone, say), so it starts
// twice, from no excess energy and from the Gᴱ the measured vapour implies,
// and keeps the lower Σr²; the first where they tie.
pub fn point_test(d PointData) !PointTest {
	m := d.x1.len
	if d.y1.len != m || d.p.len != m || d.psat1.len != m || d.psat2.len != m {
		return error('the point test needs x₁, y₁, P and both vapour pressures at every point')
	}
	if d.terms < 1 || m < d.terms + 1 {
		return error('the point test needs more points (${m}) than Legendre terms (${d.terms})')
	}
	for i in 0 .. m {
		if !(d.x1[i] > 0.0 && d.x1[i] < 1.0) || !(d.p[i] > 0.0) || !(d.psat1[i] > 0.0)
			|| !(d.psat2[i] > 0.0) {
			return error('point ${i + 1}: x₁ ${d.x1[i]} must be inside (0, 1), and P and the vapour pressures above 0')
		}
	}
	mut problem := PointProblem{d}
	mut s := numeric.least_squares(mut problem, []f64{len: d.terms}, m)!
	if start := measured_start(d) {
		if other := numeric.least_squares(mut problem, start, m) {
			if other.ssr < s.ssr {
				s = other
			}
		}
	}
	mut dy := []f64{len: m}
	mut dp := []f64{len: m}
	mut sum, mut worst := 0.0, 0.0
	for i, x1 in d.x1 {
		c1, c2 := ln_gamma_terms(x1, d.terms)
		mut l1, mut l2 := 0.0, 0.0
		for k in 0 .. d.terms {
			l1 += s.theta[k] * c1[k]
			l2 += s.theta[k] * c2[k]
		}
		a := x1 * math.exp(l1) * d.psat1[i]
		b := (1.0 - x1) * math.exp(l2) * d.psat2[i]
		dy[i] = a / (a + b) - d.y1[i]
		dp[i] = (a + b - d.p[i]) / d.p[i]
		sum += math.abs(dy[i])
		worst = math.max(worst, math.abs(dy[i]))
	}
	mean := sum / f64(m)
	return PointTest{
		coefficients: s.theta
		dy:           dy
		dp:           dp
		mean_dy:      mean
		max_dy:       worst
		consistent:   mean <= point_test_limit
		iterations:   s.iterations
	}
}
