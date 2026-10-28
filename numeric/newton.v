module numeric

import math

// Globalization is how newton keeps a step from making things worse when
// the full Newton step is too long or points the wrong way.
pub enum Globalization {
	// line_search backtracks along the Newton direction until ‖F‖ falls
	// enough (Dennis & Schnabel's A6.3.1).
	line_search
	// dogleg keeps a trust region and steps along Powell's dogleg path
	// between the steepest-descent and Newton steps (Dennis & Schnabel's
	// A6.4.3–A6.4.5).
	dogleg
}

// NewtonParams are newton's options. `@[params]` lets a caller name just
// the ones it changes, as in `newton(f, x0, method: .dogleg)`; the rest
// keep these defaults. It's the closest V gets to keyword arguments.
@[params]
pub struct NewtonParams {
pub:
	method Globalization = .line_search
	// Converged when every |Fᵢ| ≤ ftol, so the system should be scaled for
	// that to mean something.
	ftol f64 = 1e-10
	// Also converged when a full Newton step moves no xⱼ by more than
	// xtol·max(|xⱼ|, 1): x is then that accurate.
	xtol     f64 = 1e-13
	max_iter int = 200
	// The longest step allowed, as a multiple of max(‖x₀‖, 1). It keeps a
	// near-singular Jacobian from throwing x far away. Dennis and Schnabel
	// use 1000.
	max_step f64 = 1000.0
}

// NewtonSolution is a converged solve.
pub struct NewtonSolution {
pub:
	x           []f64
	f           []f64 // F(x)
	residual    f64   // max |Fᵢ(x)|
	iterations  int   // steps taken
	evaluations int   // calls to the system
}

// NewtonError is newton's failure. It carries the reason and the best point
// reached, so a caller can try something else from there. V errors are any
// type that embeds `Error` and has a msg method, and a caller can test for
// this one with `if err is NewtonError`.
pub struct NewtonError {
	Error
pub:
	reason      string
	x           []f64 // the point with the smallest ‖F‖ found
	residual    f64   // max |Fᵢ| there
	iterations  int
	evaluations int
}

// msg is the reason, with the steps taken and the residual left.
pub fn (e NewtonError) msg() string {
	return 'newton: ${e.reason} (${e.iterations} steps, ${e.evaluations} evaluations, max |F| = ${e.residual})'
}

// The sufficient-decrease constant in the Armijo condition, and the
// smallest relative step the line search tries before giving up (ε^(2/3),
// as Dennis and Schnabel suggest).
const armijo = 1e-4
const min_relative_step = 3.7e-11

// newton solves F(x) = 0 from x0 by Newton's method, made globally
// convergent by a line search or a trust region (params.method). Each step
// solves J·p = −F by LU. Where J is singular to working precision it uses
// the regularised step (JᵀJ + μI)p = −JᵀF instead, which still goes
// downhill on ‖F‖².
//
// It succeeds only at a point where every |Fᵢ| ≤ ftol, or where a full
// Newton step was shorter than xtol. Anything else is a NewtonError, never a
// point that merely stopped moving: a local minimum of ‖F‖ that isn't a
// root, a Jacobian that stays singular, or max_iter running out.
pub fn newton(system System, x0 []f64, params NewtonParams) !NewtonSolution {
	n := x0.len
	if n == 0 {
		return error('newton: x0 is empty')
	}
	for i, v in x0 {
		if !math.is_finite(v) {
			return error('newton: x0[${i}] = ${v} is not finite')
		}
	}
	if !(params.ftol > 0.0) || !(params.xtol > 0.0) || params.max_iter < 1 || !(params.max_step > 0.0) {
		return error('newton: ftol, xtol, max_iter and max_step must all be positive')
	}
	mut s := Newton{
		system:   system
		n:        n
		params:   params
		max_step: params.max_step * math.max(norm2(x0), 1.0)
		x:        x0.clone()
		f:        []f64{len: n}
		jac:      []f64{len: n * n}
	}
	if !s.evaluate(s.x, mut s.f, mut s.jac) {
		return error('newton: the system cannot be evaluated at x0')
	}
	if max_abs(s.f) <= params.ftol {
		return s.solution()
	}
	return match params.method {
		.line_search { s.line_search()! }
		.dogleg { s.dogleg()! }
	}
}

// Newton is one solve in progress: the current point, F and J there, and
// the counts.
struct Newton {
	system   System @[required]
	n        int
	params   NewtonParams
	max_step f64
mut:
	x           []f64
	f           []f64
	jac         []f64
	iterations  int
	evaluations int
}

// evaluate calls the system at x, and reports whether it gave finite values.
fn (mut s Newton) evaluate(x []f64, mut f []f64, mut jac []f64) bool {
	s.evaluations++
	s.system.eval(x, mut f, mut jac) or { return false }
	for v in f {
		if !math.is_finite(v) {
			return false
		}
	}
	for v in jac {
		if !math.is_finite(v) {
			return false
		}
	}
	return true
}

fn (s &Newton) solution() NewtonSolution {
	return NewtonSolution{
		x:           s.x
		f:           s.f
		residual:    max_abs(s.f)
		iterations:  s.iterations
		evaluations: s.evaluations
	}
}

fn (s &Newton) fail(reason string) NewtonError {
	return NewtonError{
		reason:      reason
		x:           s.x
		residual:    max_abs(s.f)
		iterations:  s.iterations
		evaluations: s.evaluations
	}
}

// direction returns the step to try from the current point, and whether it
// is the true Newton step −J⁻¹F. When J is singular to working precision it
// returns the step of (JᵀJ + μI)p = −JᵀF instead, with
// μ = √(nε)·‖JᵀJ‖₁ (Dennis & Schnabel §6.5): a Levenberg–Marquardt step,
// which goes downhill on ‖F‖² and stays bounded.
fn (s &Newton) direction() ([]f64, bool) {
	n := s.n
	if lu := lu_factor(s.jac, n) {
		if p := lu.solve(s.f) {
			return p.map(-it), true
		}
	}
	mut jtj := []f64{len: n * n}
	for i in 0 .. n {
		for j in 0 .. n {
			mut acc := 0.0
			for k in 0 .. n {
				acc += s.jac[k * n + i] * s.jac[k * n + j]
			}
			jtj[i * n + j] = acc
		}
	}
	mut norm1 := 0.0
	for j in 0 .. n {
		mut col := 0.0
		for i in 0 .. n {
			col += math.abs(jtj[i * n + j])
		}
		norm1 = math.max(norm1, col)
	}
	mu := math.sqrt(f64(n) * f64_epsilon) * norm1
	for i in 0 .. n {
		jtj[i * n + i] += mu
	}
	g := s.gradient()
	if lu := lu_factor(jtj, n) {
		if p := lu.solve(g) {
			return p.map(-it), false
		}
	}
	return []f64{len: n}, false // J is zero: no direction at all
}

// gradient returns Jᵀ·F, the gradient of ½‖F‖².
fn (s &Newton) gradient() []f64 {
	n := s.n
	mut g := []f64{len: n}
	for i in 0 .. n {
		for j in 0 .. n {
			g[j] += s.jac[i * n + j] * s.f[i]
		}
	}
	return g
}

// small_step reports whether step p moves no xⱼ by more than tol·max(|xⱼ|, 1).
fn (s &Newton) small_step(p []f64, tol f64) bool {
	for j in 0 .. s.n {
		if math.abs(p[j]) > tol * math.max(math.abs(s.x[j]), 1.0) {
			return false
		}
	}
	return true
}

// relative_length is the largest |pⱼ|/max(|xⱼ|, 1).
fn (s &Newton) relative_length(p []f64) f64 {
	mut r := 0.0
	for j in 0 .. s.n {
		r = math.max(r, math.abs(p[j]) / math.max(math.abs(s.x[j]), 1.0))
	}
	return r
}

// accept moves to xt, where F and J are ft and jt, and reports whether that
// point is converged: |F| within ftol, or a full Newton step shorter than xtol.
fn (mut s Newton) accept(xt []f64, ft []f64, jt []f64, p []f64, full_newton bool) bool {
	s.iterations++
	short := full_newton && s.small_step(p, s.params.xtol)
	s.x = xt
	s.f = ft
	s.jac = jt
	return max_abs(s.f) <= s.params.ftol || short
}

// line_search is Dennis & Schnabel's A6.3.1: try the Newton step, and while
// ‖F‖² hasn't fallen by the Armijo condition, backtrack to the minimum of a
// quadratic (then cubic) model of ‖F‖² along the step.
fn (mut s Newton) line_search() !NewtonSolution {
	n := s.n
	mut xt := []f64{len: n}
	mut ft := []f64{len: n}
	mut jt := []f64{len: n * n}
	for s.iterations < s.params.max_iter {
		mut p, newton_step := s.direction()
		mut full := newton_step
		length := norm2(p)
		if length > s.max_step {
			for j in 0 .. n {
				p[j] *= s.max_step / length
			}
			full = false
		}
		if full && s.small_step(p, s.params.xtol) {
			// A Newton step this short is taken without the decrease test:
			// x is as accurate as F can be evaluated, and ‖F‖ there is
			// rounding noise that needn't fall.
			for j in 0 .. n {
				xt[j] = s.x[j] + p[j]
			}
			if s.evaluate(xt, mut ft, mut jt) {
				_ = s.accept(xt.clone(), ft.clone(), jt.clone(), p, true)
				return s.solution()
			}
		}
		phi := 0.5 * dot(s.f, s.f)
		slope := dot(s.gradient(), p) // d(½‖F‖²)/dλ at λ = 0
		if !(slope < 0.0) {
			return s.fail('no downhill direction: a stationary point of |F| that is not a root')
		}
		lambda_min := min_relative_step / s.relative_length(p)
		mut lambda := 1.0
		mut lambda_prev, mut phi_prev := 0.0, 0.0
		for {
			// The full step is always tried: close to a root it can be far
			// shorter than lambda_min allows a backtracked one to be.
			if lambda < 1.0 && lambda < lambda_min {
				return s.fail('the line search cannot reduce |F|: a local minimum of |F| that is not a root, or a singular Jacobian')
			}
			for j in 0 .. n {
				xt[j] = s.x[j] + lambda * p[j]
			}
			if !s.evaluate(xt, mut ft, mut jt) {
				// Outside the domain, or not finite there: step back hard.
				lambda_prev = 0.0
				lambda *= 0.1
				continue
			}
			phi_t := 0.5 * dot(ft, ft)
			if phi_t <= phi + armijo * lambda * slope {
				break
			}
			mut next := 0.0
			if lambda_prev == 0.0 {
				// The minimum of the quadratic through φ(0), φ'(0) and φ(λ).
				next = -slope * lambda * lambda / (2.0 * (phi_t - phi - slope * lambda))
			} else {
				// The minimum of the cubic through φ(0), φ'(0) and the last two trials.
				r1 := phi_t - phi - lambda * slope
				r2 := phi_prev - phi - lambda_prev * slope
				d := lambda - lambda_prev
				a := (r1 / (lambda * lambda) - r2 / (lambda_prev * lambda_prev)) / d
				b := (-lambda_prev * r1 / (lambda * lambda) + lambda * r2 / (lambda_prev * lambda_prev)) / d
				next = if a == 0.0 {
					-slope / (2.0 * b)
				} else {
					disc := b * b - 3.0 * a * slope
					if disc < 0.0 { 0.5 * lambda } else { (-b + math.sqrt(disc)) / (3.0 * a) }
				}
			}
			lambda_prev, phi_prev = lambda, phi_t
			lambda = if math.is_finite(next) {
				math.max(0.1 * lambda, math.min(next, 0.5 * lambda))
			} else {
				0.5 * lambda
			}
		}
		for j in 0 .. n {
			p[j] *= lambda
		}
		if s.accept(xt.clone(), ft.clone(), jt.clone(), p, full && lambda == 1.0) {
			return s.solution()
		}
	}
	return s.fail('no convergence in ${s.params.max_iter} steps')
}

// dogleg is a trust-region method (Dennis & Schnabel's A6.4.3–A6.4.5, and
// Nocedal & Wright's algorithm 4.1). Each step minimises the linear model
// ½‖F + J·p‖² within a radius Δ along the dogleg path: from x to the
// Cauchy point, the model's minimum along −∇, then on towards the Newton
// step. The radius grows when the model predicts well and shrinks when it
// doesn't.
fn (mut s Newton) dogleg() !NewtonSolution {
	n := s.n
	mut xt := []f64{len: n}
	mut ft := []f64{len: n}
	mut jt := []f64{len: n * n}
	mut radius := -1.0 // set from the first Newton step
	for s.iterations < s.params.max_iter {
		pn, newton_step := s.direction()
		g := s.gradient()
		if norm2(g) == 0.0 {
			return s.fail('a stationary point of |F| that is not a root')
		}
		if radius < 0.0 {
			radius = math.min(norm2(pn), s.max_step)
		}
		// The Cauchy point: the model's minimum along −g, at −(‖g‖²/‖Jg‖²)·g.
		jg := mat_vec(s.jac, g, n)
		cauchy := g.map(-it * dot(g, g) / dot(jg, jg))
		phi := 0.5 * dot(s.f, s.f)
		for {
			if radius < min_relative_step * math.max(norm2(s.x), 1.0) {
				return s.fail('the trust region collapsed: a local minimum of |F| that is not a root, or a singular Jacobian')
			}
			p, full := dogleg_step(pn, cauchy, g, radius)
			for j in 0 .. n {
				xt[j] = s.x[j] + p[j]
			}
			length := norm2(p)
			if !s.evaluate(xt, mut ft, mut jt) {
				radius = 0.1 * length
				continue
			}
			if full && newton_step && s.small_step(p, s.params.xtol) {
				// As in line_search: so short a Newton step is taken whatever
				// the model's ratio, which is rounding noise by now.
				_ = s.accept(xt.clone(), ft.clone(), jt.clone(), p, true)
				return s.solution()
			}
			// The model's predicted fall in ½‖F‖², against the actual one.
			jp := mat_vec(s.jac, p, n)
			mut model := 0.0
			for i in 0 .. n {
				r := s.f[i] + jp[i]
				model += r * r
			}
			predicted := phi - 0.5 * model
			actual := phi - 0.5 * dot(ft, ft)
			ratio := if predicted > 0.0 { actual / predicted } else { -1.0 }
			if ratio < 0.25 {
				radius = 0.25 * length
			} else if ratio > 0.75 && length >= 0.99 * radius {
				radius = math.min(2.0 * radius, s.max_step)
			}
			if ratio > armijo {
				if s.accept(xt.clone(), ft.clone(), jt.clone(), p, full && newton_step) {
					return s.solution()
				}
				break
			}
		}
	}
	return s.fail('no convergence in ${s.params.max_iter} steps')
}

// dogleg_step returns the point where the dogleg path leaves the trust
// region of the given radius, and whether that's the whole Newton step.
fn dogleg_step(pn []f64, cauchy []f64, g []f64, radius f64) ([]f64, bool) {
	if norm2(pn) <= radius {
		return pn, true
	}
	cn := norm2(cauchy)
	if cn >= radius {
		gn := norm2(g)
		return g.map(-it * radius / gn), false
	}
	// cauchy + τ(pn − cauchy) with ‖·‖ = radius, 0 ≤ τ ≤ 1.
	d := []f64{len: pn.len, init: pn[index] - cauchy[index]}
	a := dot(d, d)
	b := 2.0 * dot(cauchy, d)
	c := cn * cn - radius * radius
	tau := (-b + math.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
	return []f64{len: pn.len, init: cauchy[index] + tau * d[index]}, false
}

// mat_vec returns A·v for the n×n row-major matrix a.
fn mat_vec(a []f64, v []f64, n int) []f64 {
	mut out := []f64{len: n}
	for i in 0 .. n {
		mut acc := 0.0
		for j in 0 .. n {
			acc += a[i * n + j] * v[j]
		}
		out[i] = acc
	}
	return out
}

// norm2 returns the Euclidean norm ‖x‖₂.
pub fn norm2(x []f64) f64 {
	return math.sqrt(dot(x, x))
}

// max_abs returns max |xᵢ|, the infinity norm.
pub fn max_abs(x []f64) f64 {
	mut m := 0.0
	for v in x {
		m = math.max(m, math.abs(v))
	}
	return m
}
