module study

import expr
import flowsheet
import math
import numeric

// The tolerance an optimum's first-order conditions are met to (D12, D13),
// on the problem as the optimiser sees it, scaled: each input over its
// bounds, the objective over its largest derivative at the start, each
// constraint over its sides' size there. The gradients, differences of
// solves at 1e-12, are good to about 1e-8 of their scale, so 1e-6 is well
// above their noise.
pub const optimize_tol = 1e-6

// OptimizeParams are an optimisation's options: workers to take the
// gradient's inputs at once (D8), the SQP's iterations and its tolerance.
@[params]
pub struct OptimizeParams {
pub:
	workers  int = 1
	max_iter int = 100
	tol      f64 = optimize_tol
}

// Iteration is an iterate of an optimisation: its inputs and objective, in
// SI, its worst constraint's violation and the Lagrangian's gradient, both
// scaled, and the step the line search took from it (0 at the last).
pub struct Iteration {
pub:
	inputs       []f64
	objective    f64
	violation    f64
	stationarity f64
	step         f64
}

// ConstraintValue is a constraint at the optimum: its two sides, in SI,
// whether it's active, and its multiplier: how fast the optimum's objective
// improves as the constraint is relaxed, per unit of its sides, or for an
// equality as its right side grows.
pub struct ConstraintValue {
pub:
	lhs        f64
	rhs        f64
	active     bool
	multiplier f64
}

// Optimum is an optimisation's outcome (D13): the optimum, shown to be one
// by its first-order conditions, with its multipliers, or why there's none
// and the last iterate. Every number is SI.
pub struct Optimum {
pub:
	plan        flowsheet.OptimizePlan
	ok          bool
	error       string // why there's no optimum, else ''
	inputs      []f64
	objective   f64
	constraints []ConstraintValue // in the file's order
	// Each input's bounds' multipliers: how fast the objective improves as
	// that bound is relaxed, per unit of the input; zero where it's inside.
	lower       []f64
	upper       []f64
	specs       []f64 // the file's specs' values at the optimum
	iterations  []Iteration
	evaluations int // solves at points the optimiser chose
	solves      int // every solve, the gradients' included
	solution    ?flowsheet.Solution
}

// Sheet is a flowsheet as the SQP's problem (numeric.Constrained): each input
// scaled over its bounds to [0, 1], the objective and the constraints as
// Optimum says, the equalities first. Each evaluation is a converged solve
// at a tolerance of 1e-12, specs and all, started from the last (D6); the
// gradients are D10's differences, one-sided where a bound is within the
// step.
struct Sheet {
	l          flowsheet.Loaded
	plan       flowsheet.OptimizePlan
	order      []int // the constraints, equalities first, by index in the plan
	equalities int
	sign       f64 // 1 to minimise, −1 to maximise
	workers    int
mut:
	f_scale     f64
	c_scale     []f64
	here        ?flowsheet.Solution // the solution at the last evaluation's x
	at          flowsheet.Loaded    // the file at that x
	first       []f64               // the start, whose gradient was taken to scale
	first_g     []f64               // its gradient, raw, the objective's and each side's
	solves      int
	evaluations int
}

// x_of is the inputs at a scaled point, on a bound exactly where u is 0 or
// 1.
fn (s Sheet) x_of(u []f64) []f64 {
	mut x := []f64{len: u.len}
	for k, v in s.plan.vary {
		x[k] = if u[k] <= 0.0 {
			v.lo
		} else if u[k] >= 1.0 {
			v.hi
		} else {
			math.min(math.max(v.lo + u[k] * (v.hi - v.lo), v.lo), v.hi)
		}
	}
	return x
}

// exprs are what the gradient differentiates: the objective, then each
// constraint's two sides, in the SQP's order.
fn (s Sheet) exprs() ([]expr.Expr, []string) {
	mut e := [s.plan.objective.expr]
	mut names := [s.plan.objective.name]
	for i in s.order {
		q := s.plan.constraints[i]
		e << q.lhs
		e << q.rhs
		names << q.lhs.plain()
		names << q.rhs.plain()
	}
	return e, names
}

// solve solves the file at x, from the last solution where there is one.
fn (mut s Sheet) solve(x []f64) !(flowsheet.Loaded, flowsheet.Solution) {
	at := s.l.with(s.plan.vary.map(it.input), x)!
	s.solves++
	s.evaluations++
	sol := if last := s.here {
		at.solve_from(flowsheet.Start{
			last: last
		})!
	} else {
		at.solve()!
	}
	return at, sol
}

// sides are the objective and each constraint's two sides at a solution.
fn (s Sheet) sides(at flowsheet.Loaded, sol flowsheet.Solution) ![]f64 {
	mut out := [at.value(s.plan.objective.expr, sol) or {
		return error('the objective, ${s.plan.objective.expr.plain()}: ${reason_of(err)}')
	}]
	for i in s.order {
		q := s.plan.constraints[i]
		for e in [q.lhs, q.rhs] {
			out << at.value(e, sol) or { return error('${q.plain()}: ${reason_of(err)}') }
		}
	}
	return out
}

// scaled is the objective and the constraints as the SQP sees them, from
// their sides.
fn (s Sheet) scaled(v []f64, mut c []f64) f64 {
	for k, i in s.order {
		lhs, rhs := v[1 + 2 * k], v[2 + 2 * k]
		c[k] = (if s.plan.constraints[i].relation == .le { rhs - lhs } else { lhs - rhs }) / s.c_scale[k]
	}
	return s.sign * v[0] / s.f_scale
}

fn (mut s Sheet) values(u []f64, mut c []f64) !f64 {
	at, sol := s.solve(s.x_of(u))!
	v := s.sides(at, sol)!
	s.here, s.at = sol, at
	return s.scaled(v, mut c)
}

// raw is the gradient of the objective and of each constraint's sides, by
// each input in turn, at the last evaluation: g[k][j] is the j-th's by
// input k. Each input's five solves are on a worker of their own, where
// there are workers; each starts from the point's solution, so the bits are
// the same.
fn (mut s Sheet) raw(x []f64) ![][]f64 {
	sol := s.here or { return error('no solution to differentiate at') }
	exprs, names := s.exprs()
	n := x.len
	mut out := [][]f64{len: n}
	if s.workers > 1 && n > 1 {
		mut threads := []thread []Derivative{}
		for k in 0 .. n {
			threads << spawn s.slope(sol, exprs, names, k, x[k])
		}
		results := threads.wait()
		for k in 0 .. n {
			out[k] = s.checked(results[k], k)!
		}
	} else {
		for k in 0 .. n {
			out[k] = s.checked(s.slope(sol, exprs, names, k, x[k]), k)!
		}
	}
	s.solves += 5 * n
	return out
}

// slope is one input's derivatives: by central differences over a step of
// 1e-4 of the largest value its bounds allow, so that one at 0 has a step,
// or to one side where a bound is within it.
fn (s Sheet) slope(sol flowsheet.Solution, exprs []expr.Expr, names []string, k int, x f64) []Derivative {
	v := s.plan.vary[k]
	input := s.at.input(v.input.key) or { return missing(exprs.len, err.msg()) }
	h, side := step_of(x, v.lo, v.hi)
	return differences(s.at, sol, exprs, names, input, h, side)
}

// step_of is an input's step at x within [lo, hi], and the side its
// differences go: 1e-4 of the largest value the bounds allow, at most a
// quarter of the range, about x unless a bound is within the step.
fn step_of(x f64, lo f64, hi f64) (f64, Side) {
	scale := math.max(math.abs(x), math.max(math.abs(lo), math.abs(hi)))
	h := math.min(sensitivity_step * scale, 0.25 * (hi - lo))
	side := if x - h < lo {
		Side.above
	} else if x + h > hi {
		Side.below
	} else {
		Side.both
	}
	return h, side
}

// checked is an input's derivatives' values, or why one has none.
fn (s Sheet) checked(d []Derivative, k int) ![]f64 {
	for e in d {
		if e.reason != '' {
			return error('the gradient by ${s.plan.vary[k].input.key}: ${e.reason}')
		}
	}
	return d.map(it.value)
}

fn (mut s Sheet) gradients(u []f64, mut g []f64, mut a []f64) ! {
	n := u.len
	raw := if u == s.first && s.first_g.len > 0 { s.first_g_rows(n) } else { s.raw(s.x_of(u))! }
	for k in 0 .. n {
		width := s.plan.vary[k].hi - s.plan.vary[k].lo
		g[k] = s.sign * raw[k][0] / s.f_scale * width
		for j, i in s.order {
			d := raw[k][1 + 2 * j] - raw[k][2 + 2 * j]
			sense := if s.plan.constraints[i].relation == .le { -1.0 } else { 1.0 }
			a[j * n + k] = sense * d / s.c_scale[j] * width
		}
	}
}

// first_g_rows is the start's gradient, kept from scaling, as raw gives it.
fn (s Sheet) first_g_rows(n int) [][]f64 {
	w := s.first_g.len / n
	return [][]f64{len: n, init: s.first_g[index * w..(index + 1) * w]}
}

// optimize runs one of the file's optimisations (D13): the SQP (D12) on the
// flowsheet from the plan's starts, every evaluation a converged solve.
// It's scaled at the start: the objective over its largest derivative there,
// in the inputs scaled over their bounds, so that its gradient starts at
// about 1, and each constraint over the larger of its two sides there, 1 in
// SI where both are 0.
pub fn optimize(l flowsheet.Loaded, plan flowsheet.OptimizePlan, params OptimizeParams) Optimum {
	tight := flowsheet.Loaded{
		...l
		params: flowsheet.SolveParams{
			...l.params
			tol: math.min(l.params.tol, sensitivity_tol)
		}
	}
	mut sequence := []int{}
	for i, q in plan.constraints {
		if q.relation == .eq {
			sequence << i
		}
	}
	equalities := sequence.len
	for i, q in plan.constraints {
		if q.relation != .eq {
			sequence << i
		}
	}
	mut s := Sheet{
		l:          tight
		plan:       plan
		order:      sequence
		equalities: equalities
		sign:       if plan.maximize { -1.0 } else { 1.0 }
		workers:    params.workers
		at:         tight
	}
	n := plan.vary.len
	start := plan.vary.map(it.start)
	mut u := []f64{len: n}
	for k, v in plan.vary {
		u[k] = (start[k] - v.lo) / (v.hi - v.lo)
	}
	// The scales, from the start's solution and its gradient, which the
	// SQP's first gradient then reuses.
	at, sol := s.solve(start) or { return s.failed(start, 'at the start: ${err.msg()}', []) }
	v := s.sides(at, sol) or { return s.failed(start, 'at the start: ${err.msg()}', []) }
	s.here, s.at = sol, at
	raw := s.raw(start) or { return s.failed(start, 'at the start: ${err.msg()}', []) }
	mut largest := 0.0
	for k in 0 .. n {
		largest = math.max(largest, math.abs(raw[k][0]) * (plan.vary[k].hi - plan.vary[k].lo))
	}
	s.f_scale = if largest > 0.0 {
		largest
	} else if v[0] != 0.0 {
		math.abs(v[0])
	} else {
		1.0
	}
	for k in 0 .. sequence.len {
		size := math.max(math.abs(v[1 + 2 * k]), math.abs(v[2 + 2 * k]))
		s.c_scale << if size > 0.0 { size } else { 1.0 }
	}
	s.first = u.clone()
	for row in raw {
		s.first_g << row
	}
	result := numeric.sqp(mut s, u,
		lower:        []f64{len: n}
		upper:        []f64{len: n, init: 1.0}
		equalities:   equalities
		inequalities: sequence.len - equalities
		tol:          params.tol
		max_iter:     params.max_iter
	) or {
		if err is numeric.OptimizeError {
			e := err as numeric.OptimizeError
			return s.failed(s.x_of(e.x), e.reason, e.history)
		}
		return s.failed(start, err.msg(), [])
	}
	return s.optimum(result)
}

// iterations are the SQP's iterates in SI.
fn (s Sheet) iterations(h []numeric.SqpIterate) []Iteration {
	return h.map(Iteration{
		inputs:       s.x_of(it.x)
		objective:    s.sign * it.f * s.f_scale
		violation:    it.violation
		stationarity: it.stationarity
		step:         it.step
	})
}

// failed is an optimisation with no optimum, and why.
fn (s Sheet) failed(x []f64, reason string, h []numeric.SqpIterate) Optimum {
	return Optimum{
		plan:        s.plan
		error:       reason
		inputs:      x
		objective:   math.nan()
		iterations:  s.iterations(h)
		evaluations: s.evaluations
		solves:      s.solves
	}
}

// optimum is the SQP's answer in SI: the objective and the constraints'
// sides from the optimum's own solution, and the multipliers unscaled as
// Optimum says.
fn (s Sheet) optimum(r numeric.SqpSolution) Optimum {
	x := s.x_of(r.x)
	sol := s.here or { return s.failed(x, 'no solution at the optimum', r.history) }
	v := s.sides(s.at, sol) or { return s.failed(x, err.msg(), r.history) }
	mut constraints := []ConstraintValue{len: s.plan.constraints.len}
	for k, i in s.order {
		q := s.plan.constraints[i]
		lam := r.multipliers[k]
		multiplier := (if q.relation == .eq { -lam } else { lam }) * s.f_scale / s.c_scale[k]
		constraints[i] = ConstraintValue{
			lhs:        v[1 + 2 * k]
			rhs:        v[2 + 2 * k]
			active:     q.relation == .eq || lam > 0.0 || math.abs(r.c[k]) <= optimize_tol
			multiplier: multiplier
		}
	}
	mut lower := []f64{}
	mut upper := []f64{}
	for k, vary in s.plan.vary {
		width := vary.hi - vary.lo
		lower << r.lower_multipliers[k] * s.f_scale / width
		upper << r.upper_multipliers[k] * s.f_scale / width
	}
	return Optimum{
		plan:        s.plan
		ok:          true
		inputs:      x
		objective:   v[0]
		constraints: constraints
		lower:       lower
		upper:       upper
		specs:       sol.specs.map(it.value)
		iterations:  s.iterations(r.history)
		evaluations: s.evaluations
		solves:      s.solves
		solution:    sol
	}
}
