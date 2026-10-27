module numeric

import math

// ObjectiveFn evaluates a smooth function at x. It returns f(x) and writes
// the gradient into g and the Hessian into h, n×n row-major. It returns an
// error when x is outside the function's domain; minimize then takes a
// shorter step.
pub type ObjectiveFn = fn (x []f64, mut g []f64, mut h []f64) !f64

// MinimizeParams are minimize's options.
@[params]
pub struct MinimizeParams {
pub:
	// Converged when every |gᵢ| ≤ gtol,
	gtol f64 = 1e-10
	// or when an unshifted Newton step moves no xⱼ by more than
	// xtol·max(|xⱼ|, 1).
	xtol f64 = 1e-13
	// fscale is the size of the terms f is computed from. f's rounding is
	// about ε·max(|f|, fscale), and a predicted fall below that is no fall
	// at all. Leave it at 0 when f is computed without cancellation; set it
	// when f is a small difference of large terms, as a tangent-plane
	// distance 1 + Σ… near 0 is.
	fscale   f64
	max_iter int = 100
}

// MinimizeSolution is a converged minimisation.
pub struct MinimizeSolution {
pub:
	x []f64
	f f64
	g []f64 // the gradient at x
	// minimum is true when the Hessian at x is positive definite, so x is a
	// strict local minimum. A stationary point reached along a direction of
	// zero gradient can be a saddle, and then this is false.
	minimum     bool
	iterations  int
	evaluations int
}

// minimize finds a local minimum of a smooth function from x0, by Newton's
// method with a Levenberg–Marquardt shift (Nocedal & Wright, §3.4 and
// chapter 4). Each step solves (H + μI)·p = −g. The shift μ ≥ 0 is raised
// until H + μI is positive definite, so that p goes downhill, and until the
// step lowers f by at least a small fraction of what the quadratic model
// predicts. It falls again while the model predicts well, and the steps
// become plain Newton steps near a minimum. Raising μ shortens and turns the
// step towards −g, which is a trust region in all but name.
//
// Unlike newton, which looks for zeros of a gradient and can stop wherever
// ‖g‖ stops falling, this only ever accepts a step that lowers f. Along a
// long, nearly flat valley that's the difference between getting to the end
// and stopping halfway.
//
// It succeeds when |gᵢ| ≤ gtol for every i, when an unshifted Newton step is
// shorter than xtol, or when the model predicts a fall in f smaller than f's
// own rounding (see fscale), so that no step could be seen to help.
// Anything else is a NewtonError, whose residual is then max |gᵢ|.
pub fn minimize(objective ObjectiveFn, x0 []f64, params MinimizeParams) !MinimizeSolution {
	n := x0.len
	if n == 0 {
		return error('minimize: x0 is empty')
	}
	for i, v in x0 {
		if !math.is_finite(v) {
			return error('minimize: x0[${i}] = ${v} is not finite')
		}
	}
	if !(params.gtol > 0.0) || !(params.xtol > 0.0) || params.max_iter < 1 {
		return error('minimize: gtol, xtol and max_iter must all be positive')
	}
	mut s := Minimizer{
		objective: objective
		n:         n
		params:    params
		x:         x0.clone()
		g:         []f64{len: n}
		h:         []f64{len: n * n}
	}
	s.f = s.evaluate(s.x, mut s.g, mut s.h) or {
		return error('minimize: the function cannot be evaluated at x0')
	}
	mut xt := []f64{len: n}
	mut gt := []f64{len: n}
	mut ht := []f64{len: n * n}
	mut mu := 0.0
	for s.iterations < params.max_iter {
		if max_abs(s.g) <= params.gtol {
			return s.solution()
		}
		mut scale := 1.0
		for i in 0 .. n {
			scale = math.max(scale, math.abs(s.h[i * n + i]))
		}
		mu_min := 1e-10 * scale
		for {
			if mu > 1e20 * scale {
				return s.fail('no step lowers f: a stationary point the shifts cannot leave, or a function unbounded below')
			}
			mut shifted := s.h.clone()
			for i in 0 .. n {
				shifted[i * n + i] += mu
			}
			chol := cholesky(shifted, n) or {
				mu = math.max(4.0 * mu, mu_min)
				continue
			}
			sol := chol.solve(s.g) or {
				mu = math.max(4.0 * mu, mu_min)
				continue
			}
			p := sol.map(-it)
			if mu == 0.0 && s.small_step(p) {
				// An unshifted Newton step this short: x is as accurate as
				// the problem allows.
				return s.solution()
			}
			// The fall in f the shifted quadratic model predicts for p:
			// with (H + μI)p = −g it is −½gᵀp, which is positive.
			predicted := -0.5 * dot(s.g, p)
			noise := 4.0 * f64_epsilon * math.max(math.abs(s.f), params.fscale)
			if predicted <= noise {
				// f can't see any more progress, but the gradient still can.
				// Finish with plain Newton steps for as long as H is positive
				// definite and each one lowers max |g|.
				if s.polish() {
					mu = 0.0
					break
				}
				// That's convergence only if the unshifted model agrees there
				// is nothing left. A large shift makes any step look useless:
				// against a domain wall, say, where g is far from zero.
				if s.unshifted_gain() <= noise {
					return s.solution() // nothing measurable left to gain
				}
				mu = math.max(4.0 * mu, mu_min)
				continue
			}
			for i in 0 .. n {
				xt[i] = s.x[i] + p[i]
			}
			ft := s.evaluate(xt, mut gt, mut ht) or {
				mu = math.max(4.0 * mu, mu_min)
				continue
			}
			ratio := (s.f - ft) / predicted
			if ratio < 1e-4 {
				mu = math.max(4.0 * mu, mu_min)
				continue
			}
			if ratio > 0.75 {
				mu = if mu / 4.0 < mu_min { 0.0 } else { mu / 4.0 }
			} else if ratio < 0.25 {
				mu = math.max(2.0 * mu, mu_min)
			}
			s.x = xt.clone()
			s.g = gt.clone()
			s.h = ht.clone()
			s.f = ft
			s.iterations++
			break
		}
	}
	return s.fail('no convergence in ${params.max_iter} steps')
}

// Minimizer is one minimisation in progress.
struct Minimizer {
	objective ObjectiveFn @[required]
	n         int
	params    MinimizeParams
mut:
	x           []f64
	f           f64
	g           []f64
	h           []f64
	iterations  int
	evaluations int
}

// evaluate calls the function at x and returns f, or none when it failed
// or gave anything that isn't finite.
fn (mut s Minimizer) evaluate(x []f64, mut g []f64, mut h []f64) ?f64 {
	s.evaluations++
	v := s.objective(x, mut g, mut h) or { return none }
	if !math.is_finite(v) {
		return none
	}
	for e in g {
		if !math.is_finite(e) {
			return none
		}
	}
	for e in h {
		if !math.is_finite(e) {
			return none
		}
	}
	return v
}

// polish takes one unshifted Newton step, and keeps it if H is positive
// definite and the step lowers max |g|. It's for the last steps to a
// minimum, where f is flat to rounding but g is still accurate.
fn (mut s Minimizer) polish() bool {
	chol := cholesky(s.h, s.n) or { return false }
	sol := chol.solve(s.g) or { return false }
	n := s.n
	xt := []f64{len: n, init: s.x[index] - sol[index]}
	mut gt := []f64{len: n}
	mut ht := []f64{len: n * n}
	ft := s.evaluate(xt, mut gt, mut ht) or { return false }
	if max_abs(gt) >= max_abs(s.g) {
		return false
	}
	s.x = xt
	s.g = gt
	s.h = ht
	s.f = ft
	s.iterations++
	return true
}

// unshifted_gain is the fall in f the plain Newton model predicts,
// ½gᵀH⁻¹g, or infinity when H isn't positive definite and the model
// predicts no minimum at all.
fn (s &Minimizer) unshifted_gain() f64 {
	chol := cholesky(s.h, s.n) or { return math.inf(1) }
	sol := chol.solve(s.g) or { return math.inf(1) }
	return 0.5 * dot(s.g, sol)
}

fn (s &Minimizer) small_step(p []f64) bool {
	for j in 0 .. s.n {
		if math.abs(p[j]) > s.params.xtol * math.max(math.abs(s.x[j]), 1.0) {
			return false
		}
	}
	return true
}

fn (s &Minimizer) solution() MinimizeSolution {
	return MinimizeSolution{
		x:           s.x
		f:           s.f
		g:           s.g
		minimum:     if _ := cholesky(s.h, s.n) { true } else { false }
		iterations:  s.iterations
		evaluations: s.evaluations
	}
}

fn (s &Minimizer) fail(reason string) NewtonError {
	return NewtonError{
		reason:      reason
		x:           s.x
		residual:    max_abs(s.g)
		iterations:  s.iterations
		evaluations: s.evaluations
	}
}
