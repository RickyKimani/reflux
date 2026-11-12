module numeric

import math

// RODAS4, Hairer and Wanner's Rosenbrock method of order 4 (Hairer and
// Wanner II, §IV.7): six stages, stiffly accurate, L-stable, with an
// embedded solution of order 3 for the error estimate. A Rosenbrock method
// solves no nonlinear equations: each step factors E = I/(hγ) − J once and
// solves six linear systems with it, J the Jacobian at the step's start
// (D8). The coefficients are Hairer's (ROCOE, METH = 1), in his transformed
// form, and the step is his:
//
//     E·kᵢ = f(t + cᵢh, y + Σⱼ aᵢⱼkⱼ) + Σⱼ (cᵢⱼ/h)kⱼ + h·dᵢ·∂f/∂t,
//
// with y₅ = y + Σⱼ a₅ⱼkⱼ, y₆ = y₅ + k₅ and y_new = y₆ + k₆: the last two
// stages refine the solution, so k₆ is the error estimate.

const ro_gamma = 0.25
const ro_c2 = 0.386
const ro_c3 = 0.21
const ro_c4 = 0.63
const ro_d1 = 0.25
const ro_d2 = -0.1043
const ro_d3 = 0.1035
const ro_d4 = -0.3620000000000023e-01
const ro_a21 = 0.1544e+01
const ro_a31 = 0.9466785280815826
const ro_a32 = 0.2557011698983284
const ro_a41 = 0.3314825187068521e+01
const ro_a42 = 0.2896124015972201e+01
const ro_a43 = 0.9986419139977817
const ro_a51 = 0.1221224509226641e+01
const ro_a52 = 0.6019134481288629e+01
const ro_a53 = 0.1253708332932087e+02
const ro_a54 = -0.6878860361058950
const ro_c21 = -0.56688e+01
const ro_c31 = -0.2430093356833875e+01
const ro_c32 = -0.2063599157091915
const ro_c41 = -0.1073529058151375
const ro_c42 = -0.9594562251023355e+01
const ro_c43 = -0.2047028614809616e+02
const ro_c51 = 0.7496443313967647e+01
const ro_c52 = -0.1024680431464352e+02
const ro_c53 = -0.3399990352819905e+02
const ro_c54 = 0.1170890893206160e+02
const ro_c61 = 0.8083246795921522e+01
const ro_c62 = -0.7981132988064893e+01
const ro_c63 = -0.3152159432874371e+02
const ro_c64 = 0.1631930543123136e+02
const ro_c65 = -0.6058818238834054e+01

// RODAS's step control: steps growing at most 6 times and shrinking at most
// 5, with Gustafsson's predictive controller (Hairer and Wanner II,
// §IV.8), and the step halved where E is singular, five times in all at
// most, as RODAS counts them (NSING).
const ro_grow = 6.0
const ro_shrink = 5.0
const ro_singular = 5

// Rodas4 is one RODAS4 step's workspace: f, the Jacobian and ∂f/∂t at the
// step's start, E, the stages and the new point, made once for an
// integration.
struct Rodas4 {
mut:
	dy1  []f64 // f at the step's start
	jac  []f64 // ∂f/∂y there, row-major
	dfdt []f64 // ∂f/∂t there
	e    []f64
	dy   []f64
	rhs  []f64
	k1   []f64
	k2   []f64
	k3   []f64
	k4   []f64
	k5   []f64
	k6   []f64 // after a step, the error estimate: y_new − y₆
	ynew []f64 // after a step, the new point
}

fn new_rodas4(n int) Rodas4 {
	return Rodas4{
		dy1:  []f64{len: n}
		jac:  []f64{len: n * n}
		dfdt: []f64{len: n}
		e:    []f64{len: n * n}
		dy:   []f64{len: n}
		rhs:  []f64{len: n}
		k1:   []f64{len: n}
		k2:   []f64{len: n}
		k3:   []f64{len: n}
		k4:   []f64{len: n}
		k5:   []f64{len: n}
		k6:   []f64{len: n}
		ynew: []f64{len: n}
	}
}

// RoFailure is why a RODAS4 step couldn't be taken.
enum RoFailure {
	none
	singular   // E is singular
	evaluation // f failed, or wasn't finite
}

// attempt takes a step of h from (t, y), with f, J and ∂f/∂t there in dy1,
// jac and dfdt, leaving the new point in ynew and the error estimate in k6.
// It says why where it can't, with f's reason.
fn (mut r Rodas4) attempt(mut ode StiffOde, t f64, y []f64, h f64) (RoFailure, string) {
	n := y.len
	// E = I/(hγ) − J, factored once for the six stages.
	fac := 1.0 / (h * ro_gamma)
	for i in 0 .. n * n {
		r.e[i] = -r.jac[i]
	}
	for i in 0 .. n {
		r.e[i * n + i] += fac
	}
	lu := lu_factor(r.e, n) or { return RoFailure.singular, err.msg() }
	// Stage 1.
	for i in 0 .. n {
		r.rhs[i] = r.dy1[i] + h * ro_d1 * r.dfdt[i]
	}
	lu.solve_into(r.rhs, mut r.k1) or { return RoFailure.singular, err.msg() }
	// Stage 2.
	for i in 0 .. n {
		r.ynew[i] = y[i] + ro_a21 * r.k1[i]
	}
	mut why := ode_stage(mut ode, t + ro_c2 * h, r.ynew, mut r.dy)
	if why != '' {
		return RoFailure.evaluation, why
	}
	for i in 0 .. n {
		r.rhs[i] = r.dy[i] + ro_c21 / h * r.k1[i] + h * ro_d2 * r.dfdt[i]
	}
	lu.solve_into(r.rhs, mut r.k2) or { return RoFailure.singular, err.msg() }
	// Stage 3.
	for i in 0 .. n {
		r.ynew[i] = y[i] + ro_a31 * r.k1[i] + ro_a32 * r.k2[i]
	}
	why = ode_stage(mut ode, t + ro_c3 * h, r.ynew, mut r.dy)
	if why != '' {
		return RoFailure.evaluation, why
	}
	for i in 0 .. n {
		r.rhs[i] = r.dy[i] + (ro_c31 * r.k1[i] + ro_c32 * r.k2[i]) / h + h * ro_d3 * r.dfdt[i]
	}
	lu.solve_into(r.rhs, mut r.k3) or { return RoFailure.singular, err.msg() }
	// Stage 4.
	for i in 0 .. n {
		r.ynew[i] = y[i] + ro_a41 * r.k1[i] + ro_a42 * r.k2[i] + ro_a43 * r.k3[i]
	}
	why = ode_stage(mut ode, t + ro_c4 * h, r.ynew, mut r.dy)
	if why != '' {
		return RoFailure.evaluation, why
	}
	for i in 0 .. n {
		r.rhs[i] = r.dy[i] + (ro_c41 * r.k1[i] + ro_c42 * r.k2[i] + ro_c43 * r.k3[i]) / h +
			h * ro_d4 * r.dfdt[i]
	}
	lu.solve_into(r.rhs, mut r.k4) or { return RoFailure.singular, err.msg() }
	// Stage 5, at y₅.
	for i in 0 .. n {
		r.ynew[i] = y[i] + ro_a51 * r.k1[i] + ro_a52 * r.k2[i] + ro_a53 * r.k3[i] +
			ro_a54 * r.k4[i]
	}
	why = ode_stage(mut ode, t + h, r.ynew, mut r.dy)
	if why != '' {
		return RoFailure.evaluation, why
	}
	for i in 0 .. n {
		r.rhs[i] = r.dy[i] + (ro_c51 * r.k1[i] + ro_c52 * r.k2[i] + ro_c53 * r.k3[i] +
			ro_c54 * r.k4[i]) / h
	}
	lu.solve_into(r.rhs, mut r.k5) or { return RoFailure.singular, err.msg() }
	// Stage 6, at y₆ = y₅ + k₅, the embedded solution.
	for i in 0 .. n {
		r.ynew[i] += r.k5[i]
	}
	why = ode_stage(mut ode, t + h, r.ynew, mut r.dy)
	if why != '' {
		return RoFailure.evaluation, why
	}
	for i in 0 .. n {
		r.rhs[i] = r.dy[i] + (ro_c61 * r.k1[i] + ro_c62 * r.k2[i] + ro_c63 * r.k3[i] +
			ro_c64 * r.k4[i] + ro_c65 * r.k5[i]) / h
	}
	lu.solve_into(r.rhs, mut r.k6) or { return RoFailure.singular, err.msg() }
	for i in 0 .. n {
		r.ynew[i] += r.k6[i]
	}
	return RoFailure.none, ''
}

// prepare evaluates f, the Jacobian and ∂f/∂t at an accepted point, for the
// steps from it.
fn (mut r Rodas4) prepare(mut ode StiffOde, t f64, y []f64) ! {
	ode.eval(t, y, mut r.dy1)!
	if !finite(r.dy1) {
		return error("f isn't finite")
	}
	ode.jacobian(t, y, mut r.jac, mut r.dfdt)!
	if !finite(r.jac) || !finite(r.dfdt) {
		return error("the Jacobian isn't finite")
	}
}

// rodas4 integrates y' = f(t, y) from (t0, y0) to t1, either way, by RODAS4
// on the problem's own Jacobian. It returns y at t1 only once every step has
// met the tolerance. A step whose evaluation fails or isn't finite is
// rejected and halved, and so is one whose E is singular, five times in
// all at most.
pub fn rodas4(mut ode StiffOde, t0 f64, y0 []f64, t1 f64, params OdeParams) !OdeSolution {
	w := weights(t0, y0, t1, params)!
	n := y0.len
	dir := if t1 > t0 { 1.0 } else { -1.0 }
	hmax := if params.hmax > 0 { params.hmax } else { math.abs(t1 - t0) }
	mut r := new_rodas4(n)
	mut y := y0.clone()
	mut rec := Record{
		on: params.record
	}
	mut t := t0
	mut steps, mut rejected, mut evaluations, mut jacobians, mut factorisations := 0, 0, 1, 1, 0
	r.prepare(mut ode, t, y) or { return ode_fail(.bad_input, '${err.msg()} at the start', t, y, 0) }
	rec.add(t, y)
	mut h := if params.h0 != 0 {
		dir * math.min(math.abs(params.h0), hmax)
	} else {
		evaluations++
		hinit(mut ode, t, y, r.dy1, dir, 4, hmax, w)
	}
	mut reject, mut last := false, false
	mut singular := 0
	mut hacc, mut erracc := 0.0, 0.0
	mut failure := ''
	for {
		if steps + rejected >= params.max_steps {
			return ode_fail(.too_many_steps, 'more than ${params.max_steps} steps', t, y,
				steps)
		}
		if 0.1 * math.abs(h) <= math.abs(t) * uround {
			why := if failure != '' { ': ${failure}' } else { '' }
			return ode_fail(.step_too_small, 'the step fell to ${math.abs(h)}${why}', t, y,
				steps)
		}
		// Hairer's end-of-interval rule, on a new step only: a rejected one
		// only shrinks.
		if !reject && (t + 1.0001 * h - t1) * dir >= 0 {
			h = t1 - t
			last = true
		}
		factorisations++
		kind, why := r.attempt(mut ode, t, y, h)
		if kind == .singular {
			singular++
			if singular >= ro_singular {
				return ode_fail(.singular, 'I/(hγ) − J is singular, after halving the step ${ro_singular} times',
					t, y, steps)
			}
		}
		if kind == .none {
			evaluations += 5
		}
		size := if kind == .none { w.norm(r.k6, y, r.ynew) } else { math.inf(1) }
		failure = if kind != .none {
			why
		} else if !math.is_finite(size) {
			"the error estimate isn't finite"
		} else {
			''
		}
		if failure != '' {
			// Halve the step and try again from the same point.
			rejected++
			h *= 0.5
			reject, last = true, false
			continue
		}
		mut f := math.max(1.0 / ro_grow, math.min(ro_shrink, math.pow(size, 0.25) / ode_safety))
		mut hnew := h / f
		if size > 1 {
			rejected++
			reject, last = true, false
			h = hnew
			continue
		}
		steps++
		if steps > 1 {
			// Gustafsson's predictive controller.
			mut gus := (hacc / h) * math.pow(size * size / erracc, 0.25) / ode_safety
			gus = math.max(1.0 / ro_grow, math.min(ro_shrink, gus))
			f = math.max(f, gus)
			hnew = h / f
		}
		hacc = h
		erracc = math.max(1e-2, size)
		for i in 0 .. n {
			y[i] = r.ynew[i]
		}
		t = if last { t1 } else { t + h }
		rec.add(t, y)
		if last {
			return OdeSolution{
				t:              t
				y:              y
				steps:          steps
				rejected:       rejected
				evaluations:    evaluations
				jacobians:      jacobians
				factorisations: factorisations
				ts:             rec.ts
				ys:             rec.ys
			}
		}
		evaluations++
		jacobians++
		r.prepare(mut ode, t, y) or {
			return ode_fail(.bad_input, '${err.msg()} at an accepted point', t, y, steps)
		}
		if math.abs(hnew) > hmax {
			hnew = dir * hmax
		}
		if reject {
			hnew = dir * math.min(math.abs(hnew), math.abs(h))
		}
		reject = false
		h = hnew
	}
	return ode_fail(.too_many_steps, 'unreachable', t, y, steps)
}
