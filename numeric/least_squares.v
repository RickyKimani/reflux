module numeric

import math

// Residuals is a nonlinear least-squares problem for least_squares (M8's
// D14): eval writes the residuals r(θ) into r and their Jacobian
// ∂rᵢ/∂θⱼ into j, m×n row-major, or returns an error where θ is outside the
// problem's domain, which shortens the step. It takes a `mut` receiver, as
// a problem that keeps a workspace needs.
pub interface Residuals {
mut:
	eval(theta []f64, mut r []f64, mut j []f64) !
}

// LeastSquaresParams are least_squares' options. gtol is the gradient
// test's: the largest cosine between r and a column of J, which is zero at
// a minimum, whatever the residuals' and the parameters' scales.
@[params]
pub struct LeastSquaresParams {
pub:
	gtol     f64 = 1e-10
	max_iter int = 200
}

// LeastSquaresSolution is a minimum of ½‖r‖²: θ, the residuals and their
// Jacobian there, Σ rᵢ², (JᵀJ)⁻¹ (n×n, which the residuals' variance
// scales into θ's covariance), and the gradient test's value.
pub struct LeastSquaresSolution {
pub:
	theta       []f64
	r           []f64
	j           []f64
	ssr         f64
	inverse     []f64
	gradient    f64
	iterations  int
	evaluations int
}

// LeastSquaresError is least_squares' failure: the reason, and the best θ
// reached.
pub struct LeastSquaresError {
	Error
pub:
	reason      string
	theta       []f64
	ssr         f64
	iterations  int
	evaluations int
	// undetermined is the parameter the data don't determine, by index, or
	// -1.
	undetermined int = -1
}

// msg is the reason, with the iterations taken.
pub fn (e LeastSquaresError) msg() string {
	return 'least squares: ${e.reason} (${e.iterations} iterations, ${e.evaluations} evaluations)'
}

// least_squares minimises ½‖r(θ)‖² from θ0 by Levenberg and Marquardt's
// method as Moré gives it (Lecture Notes in Mathematics 630, 1978), with
// Nielsen's update of the damping:
//
// - each step δ minimises ‖J·δ + r‖² + μ·‖D·δ‖², by Householder QR on the
//   stacked [J; √μ·D] and [−r; 0], so that JᵀJ is never formed;
// - D's entries are each column of J's norm, the largest seen so far, which
//   makes the steps independent of the parameters' scales;
// - a step is taken if the residuals fall, ρ, their fall over the fall the
//   linear model predicted, above zero; then μ ← μ·max(1/3, 1 − (2ρ − 1)³),
//   and otherwise μ ← μ·ν, ν doubling each time;
// - where the model predicts a fall below 1e-10 of Σr², which its rounding
//   can't resolve, ρ is noise, and a step is taken if the gradient test
//   falls, as ρ = 1: near the minimum it falls where Σr² can no longer tell.
//
// θ is returned only where the gradient test holds, max |Jⱼᵀr|/(‖Jⱼ‖·‖r‖)
// ≤ gtol, or r is zero; with (JᵀJ)⁻¹ from J's QR, its columns scaled to
// unit length. A Jacobian whose scaled R has a diagonal entry below 1e-10 of
// its largest has a parameter the data don't determine, and is refused.
pub fn least_squares(mut problem Residuals, theta0 []f64, m int, params LeastSquaresParams) !LeastSquaresSolution {
	n := theta0.len
	if n == 0 || m < n {
		return ls_fail('it needs at least one parameter, and at least as many residuals (${m}) as parameters (${n})',
			theta0, math.nan(), 0, 0)
	}
	if !(params.gtol > 0.0) || params.max_iter < 1 {
		return ls_fail('gtol and max_iter must be positive', theta0, math.nan(), 0, 0)
	}
	for i, v in theta0 {
		if !math.is_finite(v) {
			return ls_fail('θ0[${i}] = ${v} is not finite', theta0, math.nan(), 0, 0)
		}
	}
	mut theta := theta0.clone()
	mut r := []f64{len: m}
	mut j := []f64{len: m * n}
	problem.eval(theta, mut r, mut j) or {
		return ls_fail('at the start: ${err.msg()}', theta, math.nan(), 0, 1)
	}
	if r.any(!math.is_finite(it)) || j.any(!math.is_finite(it)) {
		return ls_fail('at the start, r or J is not finite', theta, math.nan(), 0, 1)
	}
	mut evaluations := 1
	mut ssr := dot(r, r)
	mut d := column_norms(j, m, n)
	mut mu := 1e-3
	mut nu := 2.0
	mut rn := []f64{len: m}
	mut jn := []f64{len: m * n}
	for iteration in 0 .. params.max_iter {
		g := gradient_test(j, r, m, n)
		if ssr == 0.0 || g <= params.gtol {
			inverse := scaled_inverse(j, m, n) or {
				return LeastSquaresError{
					reason:       err.msg()
					theta:        theta.clone()
					ssr:          ssr
					iterations:   iteration
					evaluations:  evaluations
					undetermined: err.code()
				}
			}
			return LeastSquaresSolution{
				theta:       theta
				r:           r
				j:           j
				ssr:         ssr
				inverse:     inverse
				gradient:    g
				iterations:  iteration
				evaluations: evaluations
			}
		}
		// The step: the least-squares solution of [J; √μ·D]·δ = [−r; 0].
		mut a := []f64{len: (m + n) * n}
		mut b := []f64{len: m + n}
		for i in 0 .. m {
			for k in 0 .. n {
				a[i * n + k] = j[i * n + k]
			}
			b[i] = -r[i]
		}
		for k in 0 .. n {
			a[(m + k) * n + k] = math.sqrt(mu) * d[k]
		}
		delta := householder_solve(a, m + n, n, b) or {
			return ls_fail('the step: ${err.msg()}', theta, ssr, iteration, evaluations)
		}
		mut tn := []f64{len: n}
		mut small := true
		for k in 0 .. n {
			tn[k] = theta[k] + delta[k]
			if tn[k] != theta[k] {
				small = false
			}
		}
		if small {
			return ls_fail('the steps fell below rounding in θ with the gradient test at ${g}, above gtol: the residuals may be too noisy for gtol, or J wrong',
				theta, ssr, iteration, evaluations)
		}
		// The fall the linear model predicts, ‖r‖² − ‖r + J·δ‖².
		mut predicted := 0.0
		for i in 0 .. m {
			mut jd := 0.0
			for k in 0 .. n {
				jd += j[i * n + k] * delta[k]
			}
			predicted -= jd * (2.0 * r[i] + jd)
		}
		evaluations++
		problem.eval(tn, mut rn, mut jn) or {
			mu *= nu
			nu *= 2.0
			continue
		}
		if rn.any(!math.is_finite(it)) || jn.any(!math.is_finite(it)) {
			mu *= nu
			nu *= 2.0
			continue
		}
		ssr_n := dot(rn, rn)
		// Near the minimum the fall the model predicts is below Σr²'s own
		// rounding, tens of ε of it, and ρ is noise; there a step is judged
		// by the gradient test, which J and r give to full precision,
		// falling. 1e-10 of Σr² leaves a thousandfold margin.
		within := predicted <= 1e-10 * ssr
		rho := if within {
			if gradient_test(jn, rn, m, n) < g { 1.0 } else { -1.0 }
		} else if predicted > 0.0 {
			(ssr - ssr_n) / predicted
		} else {
			-1.0
		}
		if rho > 0.0 {
			theta = tn.clone()
			r = rn.clone()
			j = jn.clone()
			ssr = ssr_n
			norms := column_norms(j, m, n)
			for k in 0 .. n {
				d[k] = math.max(d[k], norms[k])
			}
			mu *= math.max(1.0 / 3.0, 1.0 - math.pow(2.0 * rho - 1.0, 3))
			nu = 2.0
		} else {
			mu *= nu
			nu *= 2.0
		}
	}
	return ls_fail('the gradient test not met in ${params.max_iter} iterations', theta, ssr,
		params.max_iter, evaluations)
}

fn ls_fail(reason string, theta []f64, ssr f64, iterations int, evaluations int) IError {
	return LeastSquaresError{
		reason:      reason
		theta:       theta.clone()
		ssr:         ssr
		iterations:  iterations
		evaluations: evaluations
	}
}

// column_norms is each column of the m×n matrix j's norm, 1 for a column of
// zeros, which has no scale of its own.
fn column_norms(j []f64, m int, n int) []f64 {
	mut out := []f64{len: n}
	for k in 0 .. n {
		mut s := 0.0
		for i in 0 .. m {
			s += j[i * n + k] * j[i * n + k]
		}
		out[k] = if s > 0.0 { math.sqrt(s) } else { 1.0 }
	}
	return out
}

// gradient_test is max |Jₖᵀr|/(‖Jₖ‖·‖r‖) over the columns k: the largest
// cosine between r and a column, 0 at a minimum (Moré's gtol test).
fn gradient_test(j []f64, r []f64, m int, n int) f64 {
	rnorm := math.sqrt(dot(r, r))
	if rnorm == 0.0 {
		return 0.0
	}
	norms := column_norms(j, m, n)
	mut worst := 0.0
	for k in 0 .. n {
		mut s := 0.0
		for i in 0 .. m {
			s += j[i * n + k] * r[i]
		}
		worst = math.max(worst, math.abs(s) / (norms[k] * rnorm))
	}
	return worst
}

// householder_solve is the least-squares solution of a·x = b, a m×n with
// m ≥ n, by Householder QR: Qᵀ applied to b as R is formed, then R·x = Qᵀb
// by back substitution. It fails where R has a zero on its diagonal.
fn householder_solve(a0 []f64, m int, n int, b0 []f64) ![]f64 {
	mut a := a0.clone()
	mut b := b0.clone()
	householder(mut a, m, n, mut b)
	mut x := []f64{len: n}
	for k := n - 1; k >= 0; k-- {
		mut s := b[k]
		for c in k + 1 .. n {
			s -= a[k * n + c] * x[c]
		}
		if a[k * n + k] == 0.0 {
			return error('R is singular in column ${k}')
		}
		x[k] = s / a[k * n + k]
	}
	return x
}

// householder overwrites a's upper triangle with R from a = Q·R, and b with
// Qᵀ·b.
fn householder(mut a []f64, m int, n int, mut b []f64) {
	for k in 0 .. n {
		mut norm := 0.0
		for i in k .. m {
			norm += a[i * n + k] * a[i * n + k]
		}
		norm = math.sqrt(norm)
		if norm == 0.0 {
			continue
		}
		alpha := if a[k * n + k] > 0.0 { -norm } else { norm }
		// v = a[k..m, k] − α·e₁, kept in place below the diagonal.
		mut v := []f64{len: m - k}
		for i in k .. m {
			v[i - k] = a[i * n + k]
		}
		v[0] -= alpha
		vv := dot(v, v)
		if vv == 0.0 {
			continue
		}
		for c in k .. n {
			mut s := 0.0
			for i in k .. m {
				s += v[i - k] * a[i * n + c]
			}
			f := 2.0 * s / vv
			for i in k .. m {
				a[i * n + c] -= f * v[i - k]
			}
		}
		mut s := 0.0
		for i in k .. m {
			s += v[i - k] * b[i]
		}
		f := 2.0 * s / vv
		for i in k .. m {
			b[i] -= f * v[i - k]
		}
	}
}

// scaled_inverse is (JᵀJ)⁻¹ from the QR of J with its columns scaled to unit
// length, J·S⁻¹ = Q·R: (JᵀJ)⁻¹ = S⁻¹·R⁻¹·R⁻ᵀ·S⁻¹. It refuses a J whose R has a
// diagonal entry below 1e-10 of its largest: that parameter isn't
// determined by the data.
fn scaled_inverse(j []f64, m int, n int) ![]f64 {
	s := column_norms(j, m, n)
	mut a := []f64{len: m * n}
	for i in 0 .. m {
		for k in 0 .. n {
			a[i * n + k] = j[i * n + k] / s[k]
		}
	}
	mut b := []f64{len: m}
	householder(mut a, m, n, mut b)
	mut largest := 0.0
	for k in 0 .. n {
		largest = math.max(largest, math.abs(a[k * n + k]))
	}
	for k in 0 .. n {
		if !(math.abs(a[k * n + k]) > 1e-10 * largest) {
			return error_with_code('parameter ${k} is not determined by the data: its column of J is a combination of the others',
				k)
		}
	}
	// R⁻¹, upper triangular, by back substitution on each unit vector.
	mut ri := []f64{len: n * n}
	for c in 0 .. n {
		for k := c; k >= 0; k-- {
			mut v := if k == c { 1.0 } else { 0.0 }
			for q in k + 1 .. c + 1 {
				v -= a[k * n + q] * ri[q * n + c]
			}
			ri[k * n + c] = v / a[k * n + k]
		}
	}
	mut out := []f64{len: n * n}
	for p in 0 .. n {
		for q in 0 .. n {
			mut v := 0.0
			for k in math.max(p, q) .. n {
				v += ri[p * n + k] * ri[q * n + k]
			}
			out[p * n + q] = v / (s[p] * s[q])
		}
	}
	return out
}
