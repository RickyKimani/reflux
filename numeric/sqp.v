module numeric

import math

// Constrained is a smooth constrained problem for sqp (M8's D12): minimise f(x)
// subject to c(x), its equalities c = 0 first and then its inequalities
// c ≥ 0, within the bounds sqp is given.
//
// values writes c(x) into c and returns f(x), or an error where x is outside
// the problem's domain, which shortens sqp's step. gradients writes ∇f into g
// and each ∇cᵢ into a's row i, m×n row-major, at the x of the last call to
// values, so that a problem that solves something at x, a flowsheet, can
// start its differences from that solve. Both take a `mut` receiver for
// that.
pub interface Constrained {
mut:
	values(x []f64, mut c []f64) !f64
	gradients(x []f64, mut g []f64, mut a []f64) !
}

// SqpParams are sqp's options. The bounds may be empty, for none, or one a
// variable, ±inf where there's none. tol is absolute on c and relative to
// max(1, ‖∇f‖∞) on the Lagrangian's gradient, so the caller scales its
// problem's constraints to be of order one where they matter.
@[params]
pub struct SqpParams {
pub:
	lower        []f64
	upper        []f64
	equalities   int
	inequalities int
	tol          f64 = 1e-8
	max_iter     int = 100
}

// SqpIterate is one of sqp's iterates: x, f, the largest violation of a
// constraint, the Lagrangian's gradient's largest entry relative to
// max(1, ‖∇f‖∞), and the step the line search took from it (0 at the last).
pub struct SqpIterate {
pub:
	x            []f64
	f            f64
	violation    f64
	stationarity f64
	step         f64
}

// SqpSolution is a point where the first-order conditions hold to tol:
// - c met, equalities to tol and inequalities to within −tol;
// - ∇f = Σ λᵢ·∇cᵢ + μ_lower − μ_upper, to tol;
// - each inequality's λᵢ and each bound's μ at least zero;
// - each λᵢ·cᵢ and μ times its bound's distance zero, to tol.
pub struct SqpSolution {
pub:
	x                 []f64
	f                 f64
	c                 []f64
	multipliers       []f64 // λ, one a constraint
	lower_multipliers []f64 // μ_lower, one a variable, 0 where there's no bound
	upper_multipliers []f64
	iterations        int // QPs solved, the last at the answer
	evaluations       int // calls to values
	gradients         int // calls to gradients
	history           []SqpIterate
}

// SqpFailure is what kind of failure an OptimizeError is.
pub enum SqpFailure {
	bad_input           // sizes that don't agree, bounds crossed, or a value not finite
	evaluation          // the problem failed at the start, or its gradients did
	infeasible          // the constraints' linearisation has no point that lowers their violation
	line_search         // no step along the QP's direction lowered the merit function
	too_many_iterations // max_iter QPs without the conditions met
}

// OptimizeError is sqp's failure (D12): its reason, and where it got to.
pub struct OptimizeError {
	Error
pub:
	kind        SqpFailure
	reason      string
	x           []f64
	iterations  int
	evaluations int
	history     []SqpIterate // the iterates before it failed
}

// msg is the reason, with the iterations taken.
pub fn (e OptimizeError) msg() string {
	return 'sqp: ${e.reason} (${e.iterations} iterations, ${e.evaluations} evaluations)'
}

// The line search's sufficient decrease, the least fraction of the merit's
// directional derivative a step must achieve (Nocedal and Wright, ch. 18).
const sqp_armijo = 1e-4

// How many times the line search shortens its step before giving up.
const sqp_backtracks = 40

// sqp minimises f(x) subject to c(x) and the bounds, from x0 moved into the
// bounds, by Han and Powell's sequential quadratic programming, as Nocedal
// and Wright's chapter 18 gives it:
//
// - At each iterate x, the QP with the Lagrangian's Hessian B, f's gradient
//   and c's linearisation, minimise ½·dᵀ·B·d + ∇fᵀ·d subject to
//   c + ∇c·d = 0 or ≥ 0 and the bounds on x + d, gives a step d and
//   multipliers λ̂.
// - If the linearisation has no point, the QP is relaxed as Powell's and
//   Schittkowski's: a variable δ in [0, 1] scales each equality and violated
//   inequality, c·(1 − δ), at a cost of ½·ρ·δ². δ = 1 always has d = 0, so the
//   relaxed QP has a point; a δ that stays at 1 means nothing nearby lowers
//   the violation, and the problem is refused as infeasible.
// - The first-order conditions are checked at x with λ̂ (SqpSolution), and x
//   returned once they hold.
// - Otherwise the step αd is taken by backtracking on the ℓ1 merit function
//   f + Σ νᵢ·|violationᵢ|, Powell's weights νᵢ = max(|λ̂ᵢ|, (νᵢ + |λ̂ᵢ|)/2) making
//   d a descent direction.
// - B is updated by BFGS on the Lagrangian's gradients at λ̂, with Powell's
//   damping, which keeps B positive definite; the first update scales the
//   identity it starts from by yᵀy/sᵀy.
//
// Every iterate is within the bounds, so the problem is never evaluated
// outside them.
pub fn sqp(mut problem Constrained, x0 []f64, params SqpParams) !SqpSolution {
	n := x0.len
	m := params.equalities + params.inequalities
	mut lower := params.lower.clone()
	mut upper := params.upper.clone()
	if lower.len == 0 {
		lower = []f64{len: n, init: math.inf(-1)}
	}
	if upper.len == 0 {
		upper = []f64{len: n, init: math.inf(1)}
	}
	if n == 0 || lower.len != n || upper.len != n || params.equalities < 0
		|| params.inequalities < 0 {
		return sqp_fail(.bad_input, 'x0 needs a value a variable, and the bounds one a variable',
			x0, 0, 0)
	}
	if !(params.tol > 0.0) || params.max_iter < 1 {
		return sqp_fail(.bad_input, 'tol and max_iter must be positive', x0, 0, 0)
	}
	mut x := x0.clone()
	for i in 0 .. n {
		if !math.is_finite(x0[i]) {
			return sqp_fail(.bad_input, 'x0[${i}] = ${x0[i]} is not finite', x0, 0, 0)
		}
		if math.is_nan(lower[i]) || math.is_nan(upper[i]) || !(lower[i] < upper[i])
			|| math.is_inf(lower[i], 1) || math.is_inf(upper[i], -1) {
			return sqp_fail(.bad_input, 'variable ${i} has bounds ${lower[i]} and ${upper[i]}',
				x0, 0, 0)
		}
		x[i] = math.min(math.max(x[i], lower[i]), upper[i])
	}
	mut c := []f64{len: m}
	mut f := problem.values(x, mut c) or {
		return sqp_fail(.evaluation, 'at the start: ${err.msg()}', x, 0, 1)
	}
	if !math.is_finite(f) || c.any(!math.is_finite(it)) {
		return sqp_fail(.evaluation, 'at the start, f or c is not finite', x, 0, 1)
	}
	mut evaluations := 1
	mut gradients := 1
	mut g := []f64{len: n}
	mut a := []f64{len: m * n}
	problem.gradients(x, mut g, mut a) or {
		return sqp_fail(.evaluation, 'the gradients at the start: ${err.msg()}', x, 0, 1)
	}
	if g.any(!math.is_finite(it)) || a.any(!math.is_finite(it)) {
		return sqp_fail(.evaluation, 'the gradients at the start are not finite', x, 0, 1)
	}
	mut b := []f64{len: n * n}
	for i in 0 .. n {
		b[i * n + i] = 1.0
	}
	mut nu := []f64{len: m}
	mut history := []SqpIterate{}
	mut scaled := false
	for iteration in 1 .. params.max_iter + 1 {
		step := sqp_step(b, g, a, c, x, lower, upper, params.equalities) or {
			return sqp_failed(.infeasible, err.msg(), x, iteration, evaluations, history)
		}
		lam := step.lambda
		viol := violation(c, params.equalities)
		gscale := math.max(1.0, max_abs(g))
		mut lagrangian := g.clone()
		for i in 0 .. m {
			for j in 0 .. n {
				lagrangian[j] -= lam[i] * a[i * n + j]
			}
		}
		for j in 0 .. n {
			lagrangian[j] -= step.mu_lower[j] - step.mu_upper[j]
		}
		stationarity := max_abs(lagrangian) / gscale
		if viol <= params.tol && stationarity <= params.tol
			&& complementary(c, lam, step, x, lower, upper, params.equalities, params.tol * gscale) {
			history << SqpIterate{
				x:            x.clone()
				f:            f
				violation:    viol
				stationarity: stationarity
			}
			return SqpSolution{
				x:                 x
				f:                 f
				c:                 c
				multipliers:       lam
				lower_multipliers: step.mu_lower
				upper_multipliers: step.mu_upper
				iterations:        iteration
				evaluations:       evaluations
				gradients:         gradients
				history:           history
			}
		}
		if step.delta > 1.0 - 1e-10 {
			return sqp_failed(.infeasible, "no step lowers the constraints' violation, ${viol}: the linearised constraints have no point in common, and the problem may have no feasible point near x",
				x, iteration, evaluations, history)
		}
		for i in 0 .. m {
			nu[i] = math.max(math.abs(lam[i]), 0.5 * (nu[i] + math.abs(lam[i])))
		}
		phi := f + weighted(c, nu, params.equalities)
		d := step.d
		slope := dot(g, d) - (1.0 - step.delta) * weighted(c, nu, params.equalities)
		if !(slope < 0.0) {
			return sqp_failed(.line_search, "the QP's step isn't a descent direction for the merit function (its slope is ${slope}): the gradients may be too inaccurate for tol",
				x, iteration, evaluations, history)
		}
		mut alpha := 1.0
		mut xn := []f64{len: n}
		mut cn := []f64{len: m}
		mut fn_ := 0.0
		mut accepted := false
		mut reason := ''
		for _ in 0 .. sqp_backtracks {
			for i in 0 .. n {
				xn[i] = math.min(math.max(x[i] + alpha * d[i], lower[i]), upper[i])
			}
			evaluations++
			value := problem.values(xn, mut cn) or {
				reason = err.msg()
				alpha *= 0.5
				continue
			}
			if !math.is_finite(value) || cn.any(!math.is_finite(it)) {
				reason = 'f or c is not finite'
				alpha *= 0.5
				continue
			}
			fn_ = value
			phin := fn_ + weighted(cn, nu, params.equalities)
			if phin <= phi + sqp_armijo * alpha * slope {
				accepted = true
				break
			}
			// The minimum of the quadratic through φ(0), its slope and φ(α),
			// kept within [α/10, α/2].
			q := -slope * alpha * alpha / (2.0 * (phin - phi - slope * alpha))
			alpha = math.min(math.max(q, 0.1 * alpha), 0.5 * alpha)
		}
		if !accepted {
			why := if reason == '' { '' } else { ': ${reason}' }
			return sqp_failed(.line_search, "no step along the QP's direction lowered the merit function${why}",
				x, iteration, evaluations, history)
		}
		history << SqpIterate{
			x:            x.clone()
			f:            f
			violation:    viol
			stationarity: stationarity
			step:         alpha
		}
		mut gn := []f64{len: n}
		mut an := []f64{len: m * n}
		gradients++
		problem.gradients(xn, mut gn, mut an) or {
			return sqp_failed(.evaluation, 'the gradients: ${err.msg()}', xn, iteration,
				evaluations, history)
		}
		if gn.any(!math.is_finite(it)) || an.any(!math.is_finite(it)) {
			return sqp_failed(.evaluation, 'the gradients are not finite', xn, iteration,
				evaluations, history)
		}
		// BFGS on the Lagrangian's gradients at λ̂; the bounds' part is
		// linear, and cancels.
		mut s := []f64{len: n}
		mut y := []f64{len: n}
		for j in 0 .. n {
			s[j] = xn[j] - x[j]
			y[j] = gn[j] - g[j]
			for i in 0 .. m {
				y[j] -= lam[i] * (an[i * n + j] - a[i * n + j])
			}
		}
		if !scaled && dot(s, y) > 0.0 {
			gamma := dot(y, y) / dot(s, y)
			for i in 0 .. n * n {
				b[i] *= gamma
			}
			scaled = true
		}
		bfgs_damped(mut b, s, y)
		x, c, f, g, a = xn.clone(), cn.clone(), fn_, gn.clone(), an.clone()
	}
	return sqp_failed(.too_many_iterations, 'the first-order conditions not met in ${params.max_iter} iterations',
		x, params.max_iter, evaluations, history)
}

fn sqp_fail(kind SqpFailure, reason string, x []f64, iterations int, evaluations int) IError {
	return sqp_failed(kind, reason, x, iterations, evaluations, [])
}

fn sqp_failed(kind SqpFailure, reason string, x []f64, iterations int, evaluations int, history []SqpIterate) IError {
	return OptimizeError{
		kind:        kind
		reason:      reason
		x:           x.clone()
		iterations:  iterations
		evaluations: evaluations
		history:     history
	}
}

// violation is the largest by which c misses: |cᵢ| for an equality, −cᵢ for
// an inequality below zero.
fn violation(c []f64, equalities int) f64 {
	mut v := 0.0
	for i, ci in c {
		v = math.max(v, if i < equalities { math.abs(ci) } else { -ci })
	}
	return v
}

// weighted is Σ νᵢ·violationᵢ, the merit function's penalty.
fn weighted(c []f64, nu []f64, equalities int) f64 {
	mut s := 0.0
	for i, ci in c {
		s += nu[i] * (if i < equalities { math.abs(ci) } else { math.max(0.0, -ci) })
	}
	return s
}

// complementary is whether each inequality's and each bound's multiplier is
// at least zero and, times how far from its bound it is, zero, to tol.
fn complementary(c []f64, lam []f64, step SqpStep, x []f64, lower []f64, upper []f64, equalities int, tol f64) bool {
	for i in equalities .. c.len {
		if lam[i] < -tol || math.abs(lam[i] * c[i]) > tol {
			return false
		}
	}
	for j in 0 .. x.len {
		scale := math.max(1.0, math.abs(x[j]))
		if step.mu_lower[j] < -tol || step.mu_upper[j] < -tol {
			return false
		}
		if step.mu_lower[j] > 0.0 && step.mu_lower[j] * (x[j] - lower[j]) > tol * scale {
			return false
		}
		if step.mu_upper[j] > 0.0 && step.mu_upper[j] * (upper[j] - x[j]) > tol * scale {
			return false
		}
	}
	return true
}

// SqpStep is the QP's answer at an iterate: the step, the multipliers of the
// constraints and the bounds, and δ, 0 unless the QP had to be relaxed.
struct SqpStep {
	d        []f64
	lambda   []f64
	mu_lower []f64
	mu_upper []f64
	delta    f64
}

// sqp_step solves the QP at x, and its relaxed form if the QP has no point.
fn sqp_step(b []f64, g []f64, a []f64, c []f64, x []f64, lower []f64, upper []f64, equalities int) !SqpStep {
	if s := sqp_qp(b, g, a, c, x, lower, upper, equalities, false) {
		return s
	} else {
		if err is QpError && (err as QpError).kind != .infeasible {
			return err
		}
	}
	return sqp_qp(b, g, a, c, x, lower, upper, equalities, true) or {
		return error('the linearised constraints have no point in common, and the relaxed QP, which always has one, failed too (${err.msg()}): the problem may have no feasible point near x')
	}
}

// sqp_qp builds and solves the QP at x: its rows are the constraints, then
// each finite bound on d, lower then upper. Relaxed, it has one more
// variable, δ, as sqp says.
fn sqp_qp(b []f64, g []f64, a []f64, c []f64, x []f64, lower []f64, upper []f64, equalities int, relaxed bool) !SqpStep {
	n := x.len
	m := c.len
	size := if relaxed { n + 1 } else { n }
	mut h := []f64{len: size * size}
	mut big := 1.0
	for i in 0 .. n {
		for j in 0 .. n {
			h[i * size + j] = b[i * n + j]
		}
		big = math.max(big, b[i * n + i])
	}
	mut lin := []f64{len: size}
	for i in 0 .. n {
		lin[i] = g[i]
	}
	if relaxed {
		h[n * size + n] = 1e4 * big
	}
	mut rows := [][]f64{}
	mut rhs := []f64{}
	for i in 0 .. m {
		mut row := []f64{len: size}
		for j in 0 .. n {
			row[j] = a[i * n + j]
		}
		if relaxed && (i < equalities || c[i] < 0.0) {
			row[n] = -c[i]
		}
		rows << row
		rhs << -c[i]
	}
	// Which variable and side each bound's row is.
	mut which := []int{}
	for j in 0 .. n {
		if !math.is_inf(lower[j], -1) {
			mut row := []f64{len: size}
			row[j] = 1.0
			rows << row
			rhs << lower[j] - x[j]
			which << j
		}
	}
	for j in 0 .. n {
		if !math.is_inf(upper[j], 1) {
			mut row := []f64{len: size}
			row[j] = -1.0
			rows << row
			rhs << x[j] - upper[j]
			which << -1 - j
		}
	}
	if relaxed {
		mut lo := []f64{len: size}
		lo[n] = 1.0
		rows << lo
		rhs << 0.0
		mut hi := []f64{len: size}
		hi[n] = -1.0
		rows << hi
		rhs << -1.0
	}
	s := qp(Qp{
		h:    h
		g:    lin
		rows: rows
		b:    rhs
		eq:   equalities
	})!
	mut mu_lower := []f64{len: n}
	mut mu_upper := []f64{len: n}
	for k, j in which {
		if j >= 0 {
			mu_lower[j] = s.u[m + k]
		} else {
			mu_upper[-1 - j] = s.u[m + k]
		}
	}
	return SqpStep{
		d:        s.x[..n].clone()
		lambda:   s.u[..m].clone()
		mu_lower: mu_lower
		mu_upper: mu_upper
		delta:    if relaxed { s.x[n] } else { 0.0 }
	}
}

// bfgs_damped updates B with the step s and the gradients' change y by
// BFGS, with Powell's damping: where sᵀy < 0.2·sᵀBs, y is replaced by
// θ·y + (1 − θ)·Bs, θ = 0.8·sᵀBs/(sᵀBs − sᵀy), which keeps B positive
// definite (Nocedal and Wright, procedure 18.2).
fn bfgs_damped(mut b []f64, s []f64, y []f64) {
	n := s.len
	mut bs := []f64{len: n}
	for i in 0 .. n {
		for j in 0 .. n {
			bs[i] += b[i * n + j] * s[j]
		}
	}
	sbs := dot(s, bs)
	if !(sbs > 0.0) {
		return
	}
	sy := dot(s, y)
	theta := if sy >= 0.2 * sbs { 1.0 } else { 0.8 * sbs / (sbs - sy) }
	mut r := []f64{len: n}
	for i in 0 .. n {
		r[i] = theta * y[i] + (1.0 - theta) * bs[i]
	}
	sr := dot(s, r)
	if !(sr > 0.0) {
		return
	}
	for i in 0 .. n {
		for j in 0 .. n {
			b[i * n + j] += -bs[i] * bs[j] / sbs + r[i] * r[j] / sr
		}
	}
	for i in 0 .. n {
		for j in 0 .. i {
			v := 0.5 * (b[i * n + j] + b[j * n + i])
			b[i * n + j] = v
			b[j * n + i] = v
		}
	}
}
