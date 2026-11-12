module numeric

import math

// Dormand and Prince's explicit Runge–Kutta pair of orders 5 and 4 (J.
// Comput. Appl. Math. 6, 19, 1980), as Hairer's DOPRI5 runs it: the fifth-
// order solution carried on, the fourth-order one for the error estimate,
// and the last stage's f the next step's first (FSAL). The coefficients
// are Hairer's (CDOPRI), which SciPy's RK45 has too, as exact fractions.

const dp_c2 = 1.0 / 5
const dp_c3 = 3.0 / 10
const dp_c4 = 4.0 / 5
const dp_c5 = 8.0 / 9
const dp_a21 = 1.0 / 5
const dp_a31 = 3.0 / 40
const dp_a32 = 9.0 / 40
const dp_a41 = 44.0 / 45
const dp_a42 = -56.0 / 15
const dp_a43 = 32.0 / 9
const dp_a51 = 19372.0 / 6561
const dp_a52 = -25360.0 / 2187
const dp_a53 = 64448.0 / 6561
const dp_a54 = -212.0 / 729
const dp_a61 = 9017.0 / 3168
const dp_a62 = -355.0 / 33
const dp_a63 = 46732.0 / 5247
const dp_a64 = 49.0 / 176
const dp_a65 = -5103.0 / 18656
const dp_a71 = 35.0 / 384
const dp_a73 = 500.0 / 1113
const dp_a74 = 125.0 / 192
const dp_a75 = -2187.0 / 6784
const dp_a76 = 11.0 / 84
// The error estimate's weights: the fifth-order solution's less the
// fourth-order one's, the seventh stage's included.
const dp_e1 = 71.0 / 57600
const dp_e3 = -71.0 / 16695
const dp_e4 = 71.0 / 1920
const dp_e5 = -17253.0 / 339200
const dp_e6 = 22.0 / 525
const dp_e7 = -1.0 / 40

// DOPRI5's step control: Hairer's PI controller with β = 0.04 (Hairer,
// Nørsett and Wanner I, §IV.2), steps growing at most 10 times and
// shrinking at most 5.
const dp_beta = 0.04
const dp_grow = 10.0
const dp_shrink = 5.0

// Stiffness detection, as DOPRI5's: every dp_nstiff accepted steps, and on
// each step while it suspects stiffness, h·λ is estimated from the last two
// stages. 15 suspect steps in a row is stiff; 6 that aren't clear it.
const dp_nstiff = 1000
const dp_stiff_hl = 3.25

// Dopri5 is one Dormand–Prince step's workspace: the stages, the new point
// and the error estimate, made once for an integration.
struct Dopri5 {
mut:
	k1   []f64 // f at the step's start
	k2   []f64 // after a step, f at its end: the seventh stage
	k3   []f64
	k4   []f64
	k5   []f64
	k6   []f64
	y1   []f64 // the fifth-order solution at the step's end
	ysti []f64 // the sixth stage's point, for stiffness detection
	est  []f64 // the error estimate: fifth order less fourth
}

fn new_dopri5(n int) Dopri5 {
	return Dopri5{
		k1:   []f64{len: n}
		k2:   []f64{len: n}
		k3:   []f64{len: n}
		k4:   []f64{len: n}
		k5:   []f64{len: n}
		k6:   []f64{len: n}
		y1:   []f64{len: n}
		ysti: []f64{len: n}
		est:  []f64{len: n}
	}
}

// attempt takes a step of h from (t, y), with f(t, y) in k1, leaving the new
// point in y1, f there in k2 and the error estimate in est. It returns why
// where an evaluation fails or isn't finite, '' where it's done.
fn (mut d Dopri5) attempt(mut ode Ode, t f64, y []f64, h f64) string {
	n := y.len
	for i in 0 .. n {
		d.y1[i] = y[i] + h * dp_a21 * d.k1[i]
	}
	mut why := ode_stage(mut ode, t + dp_c2 * h, d.y1, mut d.k2)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.y1[i] = y[i] + h * (dp_a31 * d.k1[i] + dp_a32 * d.k2[i])
	}
	why = ode_stage(mut ode, t + dp_c3 * h, d.y1, mut d.k3)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.y1[i] = y[i] + h * (dp_a41 * d.k1[i] + dp_a42 * d.k2[i] + dp_a43 * d.k3[i])
	}
	why = ode_stage(mut ode, t + dp_c4 * h, d.y1, mut d.k4)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.y1[i] = y[i] + h * (dp_a51 * d.k1[i] + dp_a52 * d.k2[i] + dp_a53 * d.k3[i] +
			dp_a54 * d.k4[i])
	}
	why = ode_stage(mut ode, t + dp_c5 * h, d.y1, mut d.k5)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.ysti[i] = y[i] + h * (dp_a61 * d.k1[i] + dp_a62 * d.k2[i] + dp_a63 * d.k3[i] +
			dp_a64 * d.k4[i] + dp_a65 * d.k5[i])
	}
	why = ode_stage(mut ode, t + h, d.ysti, mut d.k6)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.y1[i] = y[i] + h * (dp_a71 * d.k1[i] + dp_a73 * d.k3[i] + dp_a74 * d.k4[i] +
			dp_a75 * d.k5[i] + dp_a76 * d.k6[i])
	}
	// k2 is free again: it takes the seventh stage, f at the new point.
	why = ode_stage(mut ode, t + h, d.y1, mut d.k2)
	if why != '' {
		return why
	}
	for i in 0 .. n {
		d.est[i] = h * (dp_e1 * d.k1[i] + dp_e3 * d.k3[i] + dp_e4 * d.k4[i] + dp_e5 * d.k5[i] +
			dp_e6 * d.k6[i] + dp_e7 * d.k2[i])
	}
	return ''
}

// stiffness is DOPRI5's estimate of |h·λ| from the step just taken: the
// sixth and seventh stages' difference over their points'. Zero where the
// points coincide.
fn (d Dopri5) stiffness(h f64) f64 {
	mut num, mut den := 0.0, 0.0
	for i in 0 .. d.y1.len {
		num += (d.k2[i] - d.k6[i]) * (d.k2[i] - d.k6[i])
		den += (d.y1[i] - d.ysti[i]) * (d.y1[i] - d.ysti[i])
	}
	return if den > 0 { math.abs(h) * math.sqrt(num / den) } else { 0.0 }
}

// dopri5 integrates y' = f(t, y) from (t0, y0) to t1, either way, with
// Dormand and Prince's pair and Hairer's step control. It returns y at t1
// only once every step has met the tolerance. A step whose evaluation fails
// or isn't finite is rejected and halved. A problem it finds stiff is
// refused (`.stiff`), for rodas4.
pub fn dopri5(mut ode Ode, t0 f64, y0 []f64, t1 f64, params OdeParams) !OdeSolution {
	w := weights(t0, y0, t1, params)!
	n := y0.len
	dir := if t1 > t0 { 1.0 } else { -1.0 }
	hmax := if params.hmax > 0 { params.hmax } else { math.abs(t1 - t0) }
	mut d := new_dopri5(n)
	mut y := y0.clone()
	mut rec := Record{
		on: params.record
	}
	mut t := t0
	mut steps, mut rejected, mut evaluations := 0, 0, 1
	ode.eval(t, y, mut d.k1) or { return ode_fail(.bad_input, 'f fails at the start: ${err.msg()}', t, y, 0) }
	if !finite(d.k1) {
		return ode_fail(.bad_input, "f isn't finite at the start", t, y, 0)
	}
	rec.add(t, y)
	mut h := if params.h0 != 0 {
		dir * math.min(math.abs(params.h0), hmax)
	} else {
		evaluations++
		hinit(mut ode, t, y, d.k1, dir, 5, hmax, w)
	}
	mut facold := 1e-4
	expo := 0.2 - dp_beta * 0.75
	mut reject, mut last := false, false
	mut suspect, mut clear := 0, 0
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
		if (t + 1.01 * h - t1) * dir > 0 {
			h = t1 - t
			last = true
		}
		failure = d.attempt(mut ode, t, y, h)
		evaluations += 6
		e := if failure == '' { w.norm(d.est, y, d.y1) } else { math.inf(1) }
		if failure == '' && !math.is_finite(e) {
			failure = "the error estimate isn't finite"
		}
		if failure != '' {
			// An evaluation failed: halve the step and try again.
			rejected++
			h *= 0.5
			reject, last = true, false
			continue
		}
		fac11 := math.pow(e, expo)
		fac := math.max(1.0 / dp_grow, math.min(dp_shrink, fac11 / math.pow(facold, dp_beta) / ode_safety))
		mut hnew := h / fac
		if e > 1 {
			rejected++
			reject, last = true, false
			h /= math.min(dp_shrink, fac11 / ode_safety)
			continue
		}
		facold = math.max(e, 1e-4)
		steps++
		if steps % dp_nstiff == 0 || suspect > 0 {
			if d.stiffness(h) > dp_stiff_hl {
				clear = 0
				suspect++
				if suspect == 15 {
					return ode_fail(.stiff, 'the problem seems stiff (|h·λ| above ${dp_stiff_hl} for 15 steps in a row)',
						t, y, steps)
				}
			} else {
				clear++
				if clear == 6 {
					suspect = 0
				}
			}
		}
		for i in 0 .. n {
			d.k1[i] = d.k2[i]
			y[i] = d.y1[i]
		}
		t = if last { t1 } else { t + h }
		rec.add(t, y)
		if last {
			return OdeSolution{
				t:           t
				y:           y
				steps:       steps
				rejected:    rejected
				evaluations: evaluations
				ts:          rec.ts
				ys:          rec.ys
			}
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

// ode_stage evaluates one stage into k, and says why where f fails or
// isn't finite: '' where it's fine. A V parameter can't be a mut string, so
// the reason comes back as the result, as it would in Go.
fn ode_stage(mut ode Ode, t f64, y []f64, mut k []f64) string {
	ode.eval(t, y, mut k) or { return err.msg() }
	if !finite(k) {
		return "f isn't finite"
	}
	return ''
}
