module numeric

import math

// BroydenParams are broyden's options.
@[params]
pub struct BroydenParams {
pub:
	// Converged when every |Fᵢ| ≤ ftol, so the system should be scaled for
	// that to mean something.
	ftol     f64 = 1e-10
	max_iter int = 100
	// The first Jacobian is scale·I. For F(x) = g(x) − x, −1 makes the first
	// step x ← g(x), which is direct substitution.
	scale f64 = -1.0
}

// BroydenSolution is a converged solve.
pub struct BroydenSolution {
pub:
	x           []f64
	f           []f64 // F(x)
	residual    f64   // max |Fᵢ(x)|
	iterations  int   // steps taken
	evaluations int   // calls to the system
}

// BroydenError is broyden's failure: the reason, and the best point reached.
pub struct BroydenError {
	Error
pub:
	reason      string
	x           []f64 // the point with the smallest max |Fᵢ| found
	residual    f64   // max |Fᵢ| there
	iterations  int
	evaluations int
}

// msg is the reason, with the steps taken and the residual left.
pub fn (e BroydenError) msg() string {
	return 'broyden: ${e.reason} (${e.iterations} steps, ${e.evaluations} evaluations, max |F| = ${e.residual})'
}

// How many times a step is halved when the system can't be evaluated there.
const broyden_halvings = 30

// broyden solves F(x) = 0 from x0 by Broyden's "good" method (Math. Comp. 19,
// 1965), with no Jacobian: it's for a system where each evaluation is costly
// and a Jacobian would cost an evaluation per unknown.
//
// It keeps an approximate inverse Jacobian H, and steps by p = −H·F. After a
// step s with a change y in F, the update
//
//     H ← H + (s − H·y)·(sᵀ·H) / (sᵀ·H·y)
//
// makes H·y = s, by the Sherman–Morrison formula, and changes H as little as
// that allows. Where the denominator vanishes, H is kept.
//
// There's no line search, since a fixed point's residual needn't fall at
// every step. If the system can't be evaluated at a step, or returns
// something that isn't finite, the step is halved. If it fails at x0, that's
// the system's own error.
//
// It succeeds only at a point where every |Fᵢ| ≤ ftol. Anything else is a
// BroydenError.
pub fn broyden(mut system Equations, x0 []f64, params BroydenParams) !BroydenSolution {
	n := x0.len
	if n == 0 {
		return error('broyden: x0 is empty')
	}
	for i, v in x0 {
		if !math.is_finite(v) {
			return error('broyden: x0[${i}] = ${v} is not finite')
		}
	}
	if !(params.ftol > 0.0) || params.max_iter < 1 || params.scale == 0.0
		|| !math.is_finite(params.scale) {
		return error('broyden: ftol and max_iter must be positive and scale finite and not zero')
	}
	mut x := x0.clone()
	mut f := []f64{len: n}
	system.eval(x, mut f)!
	mut evaluations := 1
	if !f.all(math.is_finite(it)) {
		return error('broyden: the system is not finite at x0')
	}
	mut best := x.clone()
	mut best_residual := max_abs(f)
	mut h := []f64{len: n * n}
	for i in 0 .. n {
		h[i * n + i] = 1.0 / params.scale
	}
	mut p := []f64{len: n}
	mut trial := []f64{len: n}
	mut f_trial := []f64{len: n}
	mut s := []f64{len: n}
	mut y := []f64{len: n}
	mut hy := []f64{len: n}
	mut sh := []f64{len: n}
	for k in 0 .. params.max_iter {
		if max_abs(f) <= params.ftol {
			return BroydenSolution{
				x:           x
				f:           f
				residual:    max_abs(f)
				iterations:  k
				evaluations: evaluations
			}
		}
		for i in 0 .. n {
			p[i] = 0.0
			for j in 0 .. n {
				p[i] -= h[i * n + j] * f[j]
			}
		}
		mut t := 1.0
		mut stepped := false
		for _ in 0 .. broyden_halvings {
			for i in 0 .. n {
				trial[i] = x[i] + t * p[i]
			}
			evaluations++
			if _ := system.eval(trial, mut f_trial) {
				if f_trial.all(math.is_finite(it)) {
					stepped = true
					break
				}
			}
			t *= 0.5
		}
		if !stepped {
			return BroydenError{
				reason:      'no step from here can be evaluated'
				x:           best
				residual:    best_residual
				iterations:  k
				evaluations: evaluations
			}
		}
		for i in 0 .. n {
			s[i] = trial[i] - x[i]
			y[i] = f_trial[i] - f[i]
		}
		mut d := 0.0
		for i in 0 .. n {
			hy[i] = 0.0
			sh[i] = 0.0
			for j in 0 .. n {
				hy[i] += h[i * n + j] * y[j]
				sh[i] += s[j] * h[j * n + i]
			}
		}
		for i in 0 .. n {
			d += s[i] * hy[i]
		}
		if d != 0.0 && math.is_finite(d) {
			for i in 0 .. n {
				for j in 0 .. n {
					h[i * n + j] += (s[i] - hy[i]) * sh[j] / d
				}
			}
		}
		for i in 0 .. n {
			x[i] = trial[i]
			f[i] = f_trial[i]
		}
		r := max_abs(f)
		if r < best_residual {
			best_residual = r
			for i in 0 .. n {
				best[i] = x[i]
			}
		}
	}
	if max_abs(f) <= params.ftol {
		return BroydenSolution{
			x:           x
			f:           f
			residual:    max_abs(f)
			iterations:  params.max_iter
			evaluations: evaluations
		}
	}
	return BroydenError{
		reason:      'ran out of iterations'
		x:           best
		residual:    best_residual
		iterations:  params.max_iter
		evaluations: evaluations
	}
}
