module study

import expr
import flowsheet
import math

// The tolerance the solves of a derivative are converged to (D10): a
// difference over a step of 1e-4 of the input is good to about 1e-8 of the
// output's scale then, which an optimiser's gradient needs.
pub const sensitivity_tol = 1e-12

// The step, as a fraction of the input's value: the truncation error of
// the difference over h, h²·f‴/6, and its rounding, ε·f/h, balance near
// ε^(1/3), 1e-4 for ε = 1e-12.
const sensitivity_step = 1e-4

// Derivative is an output's derivative with respect to an input at a point
// (D10), by central differences: its value and its error's estimate, or
// NaN for both and why there's none.
pub struct Derivative {
pub:
	value  f64
	error  f64
	reason string
}

// derivatives are the outputs' derivatives at a solved point (D10): for
// each input in turn, for each output, `plan.of` within `plan.to`. Each
// input is moved to v − h, v − h/2, v, v + h/2 and v + h, h = 1e-4 of its
// value, and the file solved at each from the point's solution, at a
// tolerance of 1e-12, specs and all. With D(s) = (f(v + s) − f(v − s))/(2s),
// the derivative is Richardson's D(h/2) + (D(h/2) − D(h))/3. Its error's
// estimate is two parts:
// - |D(h/2) − D(h)|, three times D(h/2)'s truncation error, which the
//   extrapolation's is far below;
// - |Δ⁴|/h, from the fourth difference over the five points,
//   f(v − h) − 4f(v − h/2) + 6f(v) − 4f(v + h/2) + f(v + h): its truncation
//   is O(h⁴), so it measures the solves' own noise, about √70 times it, and
//   the extrapolation carries about 1.9 times that noise over h.
// A unit's own tolerance, a PFR's rtol, isn't tightened: its noise shows in
// the second part.
pub fn derivatives(l flowsheet.Loaded, sol flowsheet.Solution, plan flowsheet.StudyPlan) []Derivative {
	mut out := []Derivative{}
	tight := flowsheet.Loaded{
		...l
		params: flowsheet.SolveParams{
			...l.params
			tol: math.min(l.params.tol, sensitivity_tol)
		}
	}
	for input in plan.to {
		out << slopes(tight, sol, plan, input)
	}
	return out
}

// missing is n derivatives missing, for one reason.
fn missing(n int, reason string) []Derivative {
	return []Derivative{len: n, init: Derivative{
		value:  math.nan()
		error:  math.nan()
		reason: reason
	}}
}

// slopes are a study's outputs' derivatives with respect to one input, by
// central differences over a step of 1e-4 of its value (D10).
fn slopes(l flowsheet.Loaded, sol flowsheet.Solution, plan flowsheet.StudyPlan, key flowsheet.Input) []Derivative {
	// The input as this point has it: a study's value where it varies it.
	input := l.input(key.key) or { return missing(plan.of.len, err.msg()) }
	v := input.value
	if v == 0.0 {
		return missing(plan.of.len, '${input.key} is 0 here, which gives its step no scale')
	}
	return differences(l, sol, plan.of.map(plan.outputs[it].expr), plan.of.map(plan.outputs[it].name),
		input, sensitivity_step * math.abs(v), .both)
}

// Side is where an input's five solves go: about its value, as a study's
// do, or to one side of it, as an optimiser's do where a bound is within the
// step (D13).
enum Side {
	both
	above
	below
}

// differences are the derivatives of exprs, named by names, with respect to
// one input at a solved point, l at the input's value, over a step h.
//
// About v, the steps are v ± h and v ± h/2, as derivatives says. To one side,
// they're v + σ·h·k/4 for k = 0 to 4, σ = ±1, and with D(s) = (f(v + σs) −
// f(v))/(σs):
// - R₁ = 2·D(h/2) − D(h) truncates to −h²·f‴/12, and the derivative,
//   R₂ = 2·D(h/4) − D(h/2), to a quarter of that;
// - so 4/3·|R₁ − R₂| is R₁'s truncation, four times R₂'s, and |Δ⁴|/h the
//   solves' noise as before. R₂ carries about twice R₁'s noise, 10ε/h.
fn differences(l flowsheet.Loaded, sol flowsheet.Solution, exprs []expr.Expr, names []string, input flowsheet.Input, h f64, side Side) []Derivative {
	v := input.value
	mut why := []string{len: exprs.len}
	// An output with no value at the point itself has no derivative there,
	// for the reason its own value gives, not for a step's.
	for k, e in exprs {
		_ := l.value(e, sol) or {
			why[k] = '${names[k]} has no value here'
			0.0
		}
	}
	steps := match side {
		.both { [v - h, v - 0.5 * h, v, v + 0.5 * h, v + h] }
		.above { [v, v + 0.25 * h, v + 0.5 * h, v + 0.75 * h, v + h] }
		.below { [v, v - 0.25 * h, v - 0.5 * h, v - 0.75 * h, v - h] }
	}
	// f[s][k]: output k at step s.
	mut f := [][]f64{len: 5, init: []f64{len: exprs.len}}
	for s, x in steps {
		if why.all(it != '') {
			break
		}
		at := l.with([input], [x]) or {
			return missing(exprs.len, 'at ${input.key} = ${x}: ${err.msg()}')
		}
		at_x := at.solve_from(flowsheet.Start{
			last: sol
		}) or { return missing(exprs.len, 'at ${input.key} = ${x}: ${err.msg()}') }
		for k, e in exprs {
			f[s][k] = at.value(e, at_x) or {
				if why[k] == '' {
					why[k] = 'at ${input.key} = ${x}: ${reason_of(err)}'
				}
				math.nan()
			}
		}
	}
	mut out := []Derivative{}
	for k in 0 .. exprs.len {
		if why[k] != '' {
			out << Derivative{
				value:  math.nan()
				error:  math.nan()
				reason: why[k]
			}
			continue
		}
		noise := f[0][k] - 4.0 * f[1][k] + 6.0 * f[2][k] - 4.0 * f[3][k] + f[4][k]
		if side == .both {
			// Each difference over the steps the doubles actually are.
			wide := (f[4][k] - f[0][k]) / (steps[4] - steps[0])
			narrow := (f[3][k] - f[1][k]) / (steps[3] - steps[1])
			out << Derivative{
				value: narrow + (narrow - wide) / 3.0
				error: math.abs(narrow - wide) + math.abs(noise) / h
			}
			continue
		}
		whole := (f[4][k] - f[0][k]) / (steps[4] - steps[0])
		half := (f[2][k] - f[0][k]) / (steps[2] - steps[0])
		quarter := (f[1][k] - f[0][k]) / (steps[1] - steps[0])
		r1 := 2.0 * half - whole
		r2 := 2.0 * quarter - half
		out << Derivative{
			value: r2
			error: 4.0 / 3.0 * math.abs(r1 - r2) + math.abs(noise) / h
		}
	}
	return out
}
