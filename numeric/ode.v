module numeric

import math

// Initial-value problems y' = f(t, y), for M9's reactors (D8): dopri5 for
// what isn't stiff, rodas4 for what is. Both follow Hairer's codes (DOPRI5
// and RODAS, from Hairer, Nørsett and Wanner, "Solving Ordinary Differential
// Equations" I and II): the same steps, the same error norm and the same
// step-size control, so their behaviour can be held to his.

// Ode is y' = f(t, y), for dopri5. eval writes f(t, y) into dy. It returns
// an error where y is outside the problem's domain, a negative mole number
// say, and the integrator then rejects the step and tries a shorter one. Its
// receiver is mut, unlike System's: a reactor's evaluation fills buffers it
// keeps, and counts.
pub interface Ode {
mut:
	eval(t f64, y []f64, mut dy []f64) !
}

// StiffOde is an Ode with its Jacobian, for rodas4. jacobian writes
// ∂fᵢ/∂yⱼ into jac, row-major (jac[i*n + j]), the layout lu_factor takes,
// and ∂f/∂t into dfdt, zero for a problem whose f has no t in it.
pub interface StiffOde {
	Ode
mut:
	jacobian(t f64, y []f64, mut jac []f64, mut dfdt []f64) !
}

// OdeParams are an integration's options. A step is accepted when its error
// estimate, measured against atol + rtol·|y| component by component, has a
// root mean square of at most 1 (Hairer's norm).
@[params]
pub struct OdeParams {
pub:
	rtol f64 = 1e-6
	atol f64 = 1e-6
	// atols are absolute tolerances component by component, in place of
	// atol, where the components' scales differ.
	atols []f64
	// h0 is the first step's size; 0 chooses one from f at the start
	// (Hairer's HINIT).
	h0 f64
	// hmax bounds a step's size; 0 is the whole interval.
	hmax      f64
	max_steps int = 100000
	// record keeps t and y at the start and at the end of every accepted
	// step, for a profile.
	record bool
}

// OdeSolution is y at the interval's end, and what it took.
pub struct OdeSolution {
pub:
	t              f64
	y              []f64
	steps          int // accepted
	rejected       int
	evaluations    int // of f
	jacobians      int
	factorisations int
	// ts and ys are the start and every accepted step's end, if recorded.
	ts []f64
	ys [][]f64
}

// OdeFailure is what kind of failure an OdeError is.
pub enum OdeFailure {
	bad_input      // a tolerance, an interval or a starting value that can't be
	too_many_steps // more than max_steps
	step_too_small // the step fell to rounding in t
	singular       // rodas4's matrix singular, after halving the step five times
	stiff          // dopri5 found the problem stiff: a job for rodas4
}

// OdeError is every integration's failure, with where it got to: the last
// accepted point and the steps taken.
pub struct OdeError {
	Error
pub:
	kind   OdeFailure
	reason string
	t      f64
	y      []f64
	steps  int
}

// msg is the reason, where, and after how many steps.
pub fn (e OdeError) msg() string {
	return 'ode: ${e.reason} at t = ${e.t} (${e.steps} steps)'
}

fn ode_fail(kind OdeFailure, reason string, t f64, y []f64, steps int) IError {
	return OdeError{
		kind:   kind
		reason: reason
		t:      t
		y:      y.clone()
		steps:  steps
	}
}

// uround is the unit roundoff Hairer's codes use, about 2⁻⁵².
const uround = 2.3e-16

// Hairer's step-size control: a new step is at most fac_grow times the
// last and at least 1/fac_shrink of it, times a safety factor. dopri5 uses
// 10 and 5, rodas4 6 and 5.
const ode_safety = 0.9

// Weights are an integration's tolerances, checked.
struct Weights {
	rtol  f64
	atols []f64
}

// weights checks the tolerances and the problem's start: Hairer's codes
// refuse an rtol within 10 units of roundoff, and an atol that isn't above
// zero.
fn weights(t0 f64, y0 []f64, t1 f64, p OdeParams) !Weights {
	n := y0.len
	if n == 0 {
		return ode_fail(.bad_input, 'no equations', t0, y0, 0)
	}
	if !math.is_finite(t0) || !math.is_finite(t1) || t0 == t1 {
		return ode_fail(.bad_input, 'no interval from ${t0} to ${t1}', t0, y0, 0)
	}
	for i, x in y0 {
		if !math.is_finite(x) {
			return ode_fail(.bad_input, 'y[${i}] is ${x} at the start', t0, y0, 0)
		}
	}
	if !(p.rtol > 10 * uround) || math.is_inf(p.rtol, 1) {
		return ode_fail(.bad_input, 'rtol ${p.rtol} is within rounding, or not a tolerance', t0,
			y0, 0)
	}
	if p.atols.len != 0 && p.atols.len != n {
		return ode_fail(.bad_input, '${p.atols.len} absolute tolerances for ${n} equations',
			t0, y0, 0)
	}
	atols := if p.atols.len == n { p.atols.clone() } else { []f64{len: n, init: p.atol} }
	for i, a in atols {
		if !(a > 0) || math.is_inf(a, 1) {
			return ode_fail(.bad_input, 'atol ${a} for y[${i}] is not above zero', t0, y0, 0)
		}
	}
	if p.max_steps <= 0 || p.hmax < 0 || !math.is_finite(p.hmax) || !math.is_finite(p.h0) {
		return ode_fail(.bad_input, 'max_steps, h0 or hmax is not one', t0, y0, 0)
	}
	return Weights{p.rtol, atols}
}

// norm is Hairer's error norm: the root mean square of errᵢ/sᵢ with
// sᵢ = atolᵢ + rtol·max(|y0ᵢ|, |y1ᵢ|).
fn (w Weights) norm(err []f64, y0 []f64, y1 []f64) f64 {
	mut s := 0.0
	for i, e in err {
		sk := w.atols[i] + w.rtol * math.max(math.abs(y0[i]), math.abs(y1[i]))
		s += (e / sk) * (e / sk)
	}
	return math.sqrt(s / f64(err.len))
}

// finite says whether every value is finite.
fn finite(x []f64) bool {
	for v in x {
		if !math.is_finite(v) {
			return false
		}
	}
	return true
}

// hinit is Hairer's first step (HINIT, Hairer, Nørsett and Wanner I, §II.4):
// a step of 1% of y over y', scaled, then one with 1% of the error a step
// of that order would make from f's second derivative, estimated by one
// explicit Euler step. f0 is f at the start. It costs one evaluation, for
// the caller to count; one that fails leaves the first guess.
fn hinit(mut ode Ode, t f64, y []f64, f0 []f64, dir f64, order int, hmax f64, w Weights) f64 {
	n := y.len
	mut dnf, mut dny := 0.0, 0.0
	for i in 0 .. n {
		sk := w.atols[i] + w.rtol * math.abs(y[i])
		dnf += (f0[i] / sk) * (f0[i] / sk)
		dny += (y[i] / sk) * (y[i] / sk)
	}
	mut h := if dnf <= 1e-10 || dny <= 1e-10 { 1e-6 } else { math.sqrt(dny / dnf) * 0.01 }
	h = math.min(h, hmax)
	y1 := []f64{len: n, init: y[index] + dir * h * f0[index]}
	mut f1 := []f64{len: n}
	ode.eval(t + dir * h, y1, mut f1) or { return dir * h }
	if !finite(f1) {
		return dir * h
	}
	mut der2 := 0.0
	for i in 0 .. n {
		sk := w.atols[i] + w.rtol * math.abs(y[i])
		d := (f1[i] - f0[i]) / sk
		der2 += d * d
	}
	der2 = math.sqrt(der2) / h
	der12 := math.max(math.abs(der2), math.sqrt(dnf))
	h1 := if der12 <= 1e-15 {
		math.max(1e-6, h * 1e-3)
	} else {
		math.pow(0.01 / der12, 1.0 / f64(order))
	}
	return dir * math.min(math.min(100 * h, h1), hmax)
}

// Record keeps a profile's points, when asked to.
struct Record {
	on bool
mut:
	ts []f64
	ys [][]f64
}

fn (mut r Record) add(t f64, y []f64) {
	if r.on {
		r.ts << t
		r.ys << y.clone()
	}
}
