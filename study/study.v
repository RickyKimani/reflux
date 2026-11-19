// Module study runs a flowsheet file's studies (M8's D2, D7 and D9): one or
// two inputs over their values, every output at every point, each point
// with its status and, where something has no value, why.
module study

import expr
import flowsheet
import math

// Status is a point's outcome (D9).
pub enum Status {
	solved   // solved, and every output has a value
	no_value // solved, and an output has none, for the reason given
	failed   // not solved, for the reason given; its outputs have none
}

// Point is one point of a study: its inputs' values in SI, its status, each
// output's value and, where it has none, why, and what solving it took.
pub struct Point {
pub:
	index   int   // in the grid's order: the first input's values outer
	inputs  []f64 // SI
	status  Status
	outputs []f64    // SI; NaN where an output has no value
	reasons []string // each output's reason for having none, else ''
	error   string   // why a failed point failed
	specs   []f64    // each spec's value there, SI
	// What it took: the passes of blocks with tears, the specs'
	// evaluations, the integrators' accepted steps, and whether it started
	// from its neighbour (D6) or cold.
	passes      int
	evaluations int
	steps       int
	warm        bool
	// derivatives are its outputs' with respect to the plan's inputs, for
	// each input each output (D10), none where the plan asks for none.
	derivatives []Derivative
}

// Result is a study run: its plan, and every point in the grid's order.
pub struct Result {
pub:
	plan   flowsheet.StudyPlan
	points []Point
}

// RunParams are run's options.
@[params]
pub struct RunParams {
pub:
	// warm starts each point but a chunk's first from its neighbour (D6);
	// without it every point is cold.
	warm bool = true
	// workers is how many threads take chunks (D8): 1 runs them in order on
	// the calling thread.
	workers int = 1
}

// order is the grid's points in the order they're solved (D7): along the
// first input's values, and for two inputs along each row of the second's,
// alternating direction, so that each point's neighbour is the last solved.
// Each is its index in the grid's order.
pub fn order(plan flowsheet.StudyPlan) []int {
	n0 := plan.axes[0].values.len
	if plan.axes.len == 1 {
		return []int{len: n0, init: index}
	}
	n1 := plan.axes[1].values.len
	mut out := []int{cap: n0 * n1}
	for i in 0 .. n0 {
		for k in 0 .. n1 {
			j := if i % 2 == 0 { k } else { n1 - 1 - k }
			out << i * n1 + j
		}
	}
	return out
}

// values are a grid point's inputs, by its index in the grid's order.
fn values(plan flowsheet.StudyPlan, index int) []f64 {
	if plan.axes.len == 1 {
		return [plan.axes[0].values[index]]
	}
	n1 := plan.axes[1].values.len
	return [plan.axes[0].values[index / n1], plan.axes[1].values[index % n1]]
}

// ratio is the predictor's ratio for c after b after a (flowsheet.Start):
// this step over the last along the one input that moved in both, or none
// where the two steps aren't along one input in one direction.
fn ratio(a []f64, b []f64, c []f64) ?f64 {
	mut axis := -1
	for k in 0 .. a.len {
		if a[k] != b[k] || b[k] != c[k] {
			if axis >= 0 || a[k] == b[k] || b[k] == c[k] {
				return none
			}
			axis = k
		}
	}
	if axis < 0 {
		return none
	}
	r := (c[axis] - b[axis]) / (b[axis] - a[axis])
	return if r > 0 && math.is_finite(r) { r } else { none }
}

// run runs a study (D7, D8): the points in order, in chunks of the plan's
// size, each chunk's first point cold and the rest from their neighbours. A
// point that fails doesn't stop the study, and the next starts cold. With
// more than one worker, threads share the file and take whole chunks as
// they come free; since every chunk starts cold, its points are the same
// bits whichever thread solves them, and whenever.
pub fn run(l flowsheet.Loaded, plan flowsheet.StudyPlan, params RunParams) Result {
	seq := order(plan)
	mut chunks := [][]int{}
	for start := 0; start < seq.len; start += plan.chunk {
		chunks << seq[start..math.min(seq.len, start + plan.chunk)]
	}
	mut points := []Point{len: seq.len}
	workers := math.min(math.max(params.workers, 1), chunks.len)
	if workers <= 1 {
		for chunk in chunks {
			for p in run_chunk(l, plan, chunk, params) {
				points[p.index] = p
			}
		}
	} else {
		queue := chan int{cap: chunks.len}
		for k in 0 .. chunks.len {
			queue <- k
		}
		queue.close()
		mut threads := []thread []Point{}
		for _ in 0 .. workers {
			threads << spawn work(l, plan, chunks, queue, params)
		}
		for part in threads.wait() {
			for p in part {
				points[p.index] = p
			}
		}
	}
	return Result{
		plan:   plan
		points: points
	}
}

// work is one worker: chunks from the queue until it's empty.
fn work(l flowsheet.Loaded, plan flowsheet.StudyPlan, chunks [][]int, queue chan int, params RunParams) []Point {
	mut out := []Point{}
	for {
		k := <-queue or { break }
		out << run_chunk(l, plan, chunks[k], params)
	}
	return out
}

// run_chunk solves a chunk's points in order: the first cold, each other
// from its neighbour, and after a failure cold again.
fn run_chunk(l flowsheet.Loaded, plan flowsheet.StudyPlan, chunk []int, params RunParams) []Point {
	keys := plan.axes.map(it.input)
	mut out := []Point{cap: chunk.len}
	mut last := ?flowsheet.Solution(none)
	mut before := ?flowsheet.Solution(none)
	mut last_x, mut before_x := []f64{}, []f64{}
	for index in chunk {
		x := values(plan, index)
		solve := l.with(keys, x) or {
			out << failed(index, x, plan, err)
			last, before = none, none
			continue
		}
		mut warm := false
		mut from := flowsheet.Start{}
		if params.warm {
			if b := last {
				warm = true
				mut r := 1.0
				mut a := ?flowsheet.Solution(none)
				if s := before {
					if q := ratio(before_x, last_x, x) {
						a, r = s, q
					}
				}
				from = flowsheet.Start{
					last:   b
					before: a
					ratio:  r
				}
			}
		}
		sol := solved(solve, warm, from) or {
			out << Point{
				...failed(index, x, plan, err)
				warm: warm
			}
			last, before = none, none
			continue
		}
		p := measured(solve, sol, index, x, plan, warm)
		out << if plan.to.len > 0 {
			Point{
				...p
				derivatives: derivatives(solve, sol, plan)
			}
		} else {
			p
		}
		before, before_x = last, last_x.clone()
		last, last_x = sol, x.clone()
	}
	return out
}

// solved is a point solved from its neighbour, or cold.
fn solved(l flowsheet.Loaded, warm bool, from flowsheet.Start) !flowsheet.Solution {
	if warm {
		return l.solve_from(from)
	}
	return l.solve()
}

// measured is a solved point: each output's value, or why it has none.
fn measured(l flowsheet.Loaded, sol flowsheet.Solution, index int, x []f64, plan flowsheet.StudyPlan, warm bool) Point {
	mut outputs := []f64{len: plan.outputs.len}
	mut reasons := []string{len: plan.outputs.len}
	mut status := Status.solved
	for k, o in plan.outputs {
		outputs[k] = l.value(o.expr, sol) or {
			reasons[k] = reason_of(err)
			status = .no_value
			math.nan()
		}
	}
	mut passes, mut evaluations, mut steps := 0, 0, 0
	for h in sol.history {
		passes += h.residuals.len
		evaluations += h.spec_residuals.len
	}
	for u in sol.units {
		if p := u.profile {
			steps += p.points.len - 1
		}
	}
	return Point{
		index:       index
		inputs:      x
		status:      status
		outputs:     outputs
		reasons:     reasons
		specs:       sol.specs.map(it.value)
		passes:      passes
		evaluations: evaluations
		steps:       steps
		warm:        warm
	}
}

// failed is a point that didn't solve, and why.
fn failed(index int, x []f64, plan flowsheet.StudyPlan, err IError) Point {
	return Point{
		index:   index
		inputs:  x
		status:  .failed
		outputs: []f64{len: plan.outputs.len, init: math.nan()}
		reasons: []string{len: plan.outputs.len}
		error:   err.msg()
	}
}

// reason_of is an expression's reason, without its place in the text.
fn reason_of(err IError) string {
	if err is expr.ExprError {
		return err.reason
	}
	return err.msg()
}
