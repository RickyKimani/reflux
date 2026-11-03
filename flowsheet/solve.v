module flowsheet

import flash
import math
import numeric
import thermo

// Method is how a block's tears are converged (D7).
pub enum Method {
	direct   // direct substitution: each pass starts from the last one's answer
	wegstein // bounded Wegstein, every other iteration
	broyden  // Broyden's quasi-Newton method, on all a block's tear variables
}

// SolveParams are solve's options.
@[params]
pub struct SolveParams {
pub:
	method Method = .wegstein
	// tol is how far a tear variable may move in a pass, converged: a component
	// flow as a fraction of the stream's total, P of itself, and H of RT.
	tol            f64 = 1e-9
	max_iterations int = 100
	// tears are streams to tear, as for plan.
	tears []string
	// estimates are first values for streams, each of which is then torn.
	estimates map[string]Stream
	// work is the budget for the tear search, as for plan.
	work i64 = 200_000_000
}

// UnitResult is a unit as it was last calculated: its balance is exact to
// these inlets, which for a tear's consumer are a pass behind the tear's value
// in Solution.streams.
pub struct UnitResult {
pub:
	name    string
	inlets  []Stream
	outlets []Stream
	energy  []EnergyStream
}

// TearHistory is how a block's tears converged.
pub struct TearHistory {
pub:
	tears []string
	// residuals has a row for each pass, with each tear's largest change as a
	// fraction of its scale.
	residuals [][]f64
}

// Balance is how far a balance is from closing: the worst component's
// imbalance as a fraction of the flow, and the energy's as a fraction of F·RT.
pub struct Balance {
pub:
	component f64
	energy    f64
}

// Solution is a flowsheet solved: every tear converged and every balance
// checked (D8).
pub struct Solution {
pub:
	streams  map[string]Stream
	units    []UnitResult // in the order of Flowsheet.units
	products []string
	plan     Plan
	history  []TearHistory
	balance  Balance // of the whole flowsheet
	worst    Balance // of the worst unit, for each balance
}

// How far a unit's own balances may be from closing (criterion 2), and how
// far rounding may take the whole, as a fraction of its flows.
const unit_tolerance = 1e-10
const round_off = 1e-12

// solve solves a flowsheet: each block once, in the plan's order, and a block
// with tears until they converge (D7). A tear starts from its estimate, or
// from what enters its block. Where its unit has another inlet it starts empty,
// so that the first pass is once-through. Otherwise it starts as a mixture of
// what enters the block, at its lowest pressure.
//
// Failure is a FlowsheetError, with the history so far: a unit's failure as it
// happened, .no_convergence, or .unbalanced.
pub fn (f Flowsheet) solve(fl thermo.Fluid, params SolveParams) !Solution {
	if !(params.tol > 0.0) || params.max_iterations < 1 {
		return bad_spec('', 'tol and max_iterations must be above zero')
	}
	mut named := params.tears.clone()
	for name, _ in params.estimates {
		named << name
	}
	plan := f.plan(tears: named, work: params.work)!
	mut s := Solver{
		sheet:  f
		fl:     fl
		params: params
		nc:     fl.model.components().len
		runs:   []UnitResult{len: f.units.len}
	}
	for feed in f.feeds {
		s.streams[feed.name] = feed.stream
	}
	for b in plan.blocks {
		if b.tears.len == 0 {
			s.pass(b) or { return s.with_history(err) }
		} else {
			s.converge(b)!
		}
	}
	return s.solution(plan)
}

// Solver is one solve. It's its own numeric.Equations, for Broyden.
struct Solver {
	sheet  Flowsheet
	fl     thermo.Fluid
	params SolveParams
	nc     int
mut:
	streams map[string]Stream
	runs    []UnitResult
	history []TearHistory
	// The block in hand.
	block Block
	flow  f64 // what enters it, mol/s
	rows  [][]f64
	scale []f64 // Broyden's unknowns are the tear variables over these
	fresh bool  // the tears hold their first values, not yet used
	// What the last pass started from and came to, and so what the tears of the
	// blocks so far have left open: each component's flow, and the energy.
	x_last     []f64
	g_last     []f64
	gap_flow   f64
	gap_energy f64
}

// pass calculates a block's units once, in order.
fn (mut s Solver) pass(b Block) ! {
	for i in b.units {
		p := s.sheet.units[i]
		inlets := p.inlets.map(s.streams[it])
		o := p.unit.calculate(s.fl, inlets)!
		if o.outlets.len != p.outlets.len {
			return bad_spec(p.unit.name, 'it gave ${o.outlets.len} outlets for ${p.outlets.len} ports')
		}
		for k, name in p.outlets {
			s.streams[name] = o.outlets[k]
		}
		s.runs[i] = UnitResult{
			name:    p.unit.name
			inlets:  inlets
			outlets: o.outlets
			energy:  o.energy
		}
	}
}

// begin gives a block's tears their first values.
fn (mut s Solver) begin(b Block) ! {
	mut entering := []Stream{}
	for i in b.units {
		for name in s.sheet.units[i].inlets {
			if name in s.streams {
				entering << s.streams[name]
			}
		}
	}
	mix := blend('', s.fl, entering, entering.len)!
	r := flash.ph(s.fl, mix.p, mix.enthalpy / mix.total, mix.flows, estimate: mix.t) or {
		return flash_failed('', err)
	}
	full := of_flash(s.fl, mix.flows, r)
	empty := Stream{
		...full
		flows: []f64{len: s.nc}
	}
	s.block = b
	s.flow = mix.total
	s.rows = [][]f64{}
	for name in b.tears {
		s.streams[name] = s.params.estimates[name] or {
			inlets := s.sheet.units.filter(name in it.inlets)[0].inlets
			if inlets.any(it !in b.tears) { empty } else { full }
		}
	}
}

// converge iterates a block until its tears settle.
fn (mut s Solver) converge(b Block) ! {
	s.begin(b)!
	if s.params.method == .broyden {
		s.by_broyden() or { return s.failed(err) }
	} else {
		s.by_substitution() or { return s.failed(err) }
	}
	s.keep_history()
	s.add_gap()
}

// add_gap counts what a converged block's tears left open: the balance of the
// whole can't close better than the difference between the last pass's start
// and its end, and it should close to that.
fn (mut s Solver) add_gap() {
	n := s.nc + 2
	for t in 0 .. s.block.tears.len {
		mut flow_x, mut flow_g, mut worst := 0.0, 0.0, 0.0
		for i in 0 .. s.nc {
			flow_x += s.x_last[t * n + i]
			flow_g += s.g_last[t * n + i]
			worst = math.max(worst, math.abs(s.g_last[t * n + i] - s.x_last[t * n + i]))
		}
		s.gap_flow += worst
		s.gap_energy += math.abs(flow_g * s.g_last[t * n + s.nc + 1] - flow_x * s.x_last[t * n +
			s.nc + 1])
	}
}

// failed adds the history to a failure in a block with tears, the block's own
// passes included.
fn (mut s Solver) failed(err IError) IError {
	s.keep_history()
	return s.with_history(err)
}

// with_history adds the history so far to a failure.
fn (s Solver) with_history(err IError) IError {
	if err is FlowsheetError {
		return FlowsheetError{
			kind:    err.kind
			unit:    err.unit
			reason:  err.reason
			cause:   err.cause
			history: s.history
		}
	}
	return err
}

fn (mut s Solver) keep_history() {
	s.history << TearHistory{
		tears:     s.block.tears
		residuals: s.rows
	}
}

// tear_values is the tears' variables, for each its flows, P and molar H.
fn (s Solver) tear_values() []f64 {
	mut x := []f64{cap: s.block.tears.len * (s.nc + 2)}
	for name in s.block.tears {
		t := s.streams[name]
		x << t.flows
		x << t.p
		x << t.h
	}
	return x
}

// compare writes into f how far each variable moved in a pass, as a fraction
// of its scale (D7), keeps each tear's largest in the history, and returns the
// largest of all.
fn (mut s Solver) compare(x []f64, g []f64, mut f []f64) f64 {
	n := s.nc + 2
	s.x_last = x
	s.g_last = g
	mut row := []f64{len: s.block.tears.len}
	for t, name in s.block.tears {
		c := s.streams[name]
		total := if c.total() > 0.0 { c.total() } else { s.flow }
		for i in 0 .. n {
			scale := if i < s.nc {
				total
			} else if i == s.nc {
				c.p
			} else {
				thermo.r_gas * c.t
			}
			k := t * n + i
			f[k] = (g[k] - x[k]) / scale
			row[t] = math.max(row[t], math.abs(f[k]))
		}
	}
	s.rows << row
	mut worst := 0.0
	for r in row {
		worst = math.max(worst, r)
	}
	return worst
}

// set_tears rebuilds the tears from their variables by PH flashes, each from
// its last T (D6). An acceleration can overshoot, so a flow is held at 0. With
// no flow, a tear keeps its last state.
fn (mut s Solver) set_tears(x []f64) ! {
	n := s.nc + 2
	for t, name in s.block.tears {
		last := s.streams[name]
		v := x[t * n..(t + 1) * n]
		flows := v[..s.nc].map(math.max(it, 0.0))
		if numeric.sum(flows) == 0.0 {
			s.streams[name] = Stream{
				...last
				flows: flows
			}
			continue
		}
		r := flash.ph(s.fl, v[s.nc], v[s.nc + 1], flows, estimate: last.t) or {
			return FlowsheetError{
				kind:   .flash_failed
				reason: 'tear ${name}: ${err.msg()}'
				cause:  if err is flash.FlashError { err.kind } else { none }
			}
		}
		s.streams[name] = of_flash(s.fl, flows, r)
	}
}

// wegstein is the next x from the last two passes, variable by variable
// (Wegstein, Commun. ACM 1, 1958). With s the slope of g, q = s/(s − 1) gives
// x ← q·x + (1 − q)·g(x), held to q in [−5, 0] so that it only ever
// accelerates, and by no more than 5 times the step (D7).
fn wegstein(x0 []f64, g0 []f64, x []f64, g []f64) []f64 {
	mut next := g.clone()
	for i in 0 .. x.len {
		dx := x[i] - x0[i]
		if dx == 0.0 {
			continue
		}
		slope := (g[i] - g0[i]) / dx
		q := math.max(-5.0, math.min(0.0, slope / (slope - 1.0)))
		next[i] = q * x[i] + (1.0 - q) * g[i]
	}
	return next
}

// by_substitution is direct substitution, with Wegstein's step every other
// iteration if that's the method, from the third. The first pass starts from a
// guess, an empty stream as like as not, and a slope through that says nothing:
// taking one almost doubled the passes on a nonlinear loop. A direct step needs
// no flash: the pass has left the tear's new value in place.
fn (mut s Solver) by_substitution() ! {
	mut f := []f64{len: s.block.tears.len * (s.nc + 2)}
	mut x0, mut g0 := []f64{}, []f64{}
	for k in 1 .. s.params.max_iterations + 1 {
		x := s.tear_values()
		s.pass(s.block)!
		g := s.tear_values()
		if s.compare(x, g, mut f) <= s.params.tol {
			return
		}
		if s.params.method == .wegstein && k >= 3 && k % 2 == 1 {
			s.set_tears(wegstein(x0, g0, x, g))!
		}
		x0 = x.clone()
		g0 = g.clone()
	}
	return s.not_converged('the tears are still moving')
}

fn (s Solver) not_converged(reason string) IError {
	last := s.rows.last()
	return FlowsheetError{
		kind:   .no_convergence
		reason: 'tears ${s.block.tears.join(', ')}: ${reason} after ${s.rows.len} passes, by ${last}'
	}
}

// by_broyden runs numeric.broyden on all the block's tear variables together,
// after one pass, for the same reason as Wegstein: its first secant through a
// guess would stay in its Jacobian.
fn (mut s Solver) by_broyden() ! {
	start := s.tear_values()
	s.pass(s.block)!
	mut change := []f64{len: start.len}
	if s.compare(start, s.tear_values(), mut change) <= s.params.tol {
		return
	}
	x := s.tear_values()
	n := s.nc + 2
	s.scale = []f64{len: x.len}
	for t in 0 .. s.block.tears.len {
		for i in 0 .. s.nc {
			s.scale[t * n + i] = s.flow
		}
		first := s.streams[s.block.tears[t]]
		s.scale[t * n + s.nc] = first.p
		s.scale[t * n + s.nc + 1] = thermo.r_gas * first.t
	}
	mut y := []f64{len: x.len}
	for i, v in x {
		y[i] = v / s.scale[i]
	}
	s.fresh = true
	numeric.broyden(mut s, y, ftol: s.params.tol, max_iter: s.params.max_iterations - 2) or {
		if err is numeric.BroydenError {
			return s.not_converged(err.reason)
		}
		return err
	}
}

// eval is one pass, for numeric.broyden: the tears' changes, at y.
fn (mut s Solver) eval(y []f64, mut f []f64) ! {
	if s.fresh {
		s.fresh = false
	} else {
		mut x := []f64{len: y.len}
		for i, v in y {
			x[i] = v * s.scale[i]
		}
		s.set_tears(x)!
	}
	x := s.tear_values()
	s.pass(s.block)!
	s.compare(x, s.tear_values(), mut f)
}

// products is the streams that leave the flowsheet.
fn (f Flowsheet) products() []string {
	mut used := map[string]bool{}
	for p in f.units {
		for name in p.inlets {
			used[name] = true
		}
	}
	mut out := []string{}
	for p in f.units {
		for name in p.outlets {
			if name !in used {
				out << name
			}
		}
	}
	return out
}

// imbalance is the worst component's imbalance, and the energy's, from the
// streams in and out and the energy added.
fn imbalance(inlets []Stream, outlets []Stream, added f64, nc int) (f64, f64) {
	mut comp := []f64{len: nc}
	mut energy := added
	for s in inlets {
		energy += s.total() * s.h
		for i, n in s.flows {
			comp[i] += n
		}
	}
	for s in outlets {
		energy -= s.total() * s.h
		for i, n in s.flows {
			comp[i] -= n
		}
	}
	mut worst := 0.0
	for c in comp {
		worst = math.max(worst, math.abs(c))
	}
	return worst, math.abs(energy)
}

// added is the sum of some energy streams.
fn added(energy []EnergyStream) f64 {
	mut q := 0.0
	for e in energy {
		q += e.rate
	}
	return q
}

// balance is how far a unit closes on the inlets it was calculated from: its
// worst component against its flow, and its energy against F·RT at its hottest
// inlet. With nothing coming in there's nothing to measure against, and solve's
// check of the whole flowsheet is what would catch such a unit making material.
pub fn (r UnitResult) balance() Balance {
	mut flow, mut t := 0.0, 0.0
	for x in r.inlets {
		flow += x.total()
		t = math.max(t, x.t)
	}
	if flow == 0.0 {
		return Balance{}
	}
	comp, energy := imbalance(r.inlets, r.outlets, added(r.energy), r.inlets[0].flows.len)
	return Balance{comp / flow, energy / (flow * thermo.r_gas * t)}
}

// solution checks every balance and gives the answer (D8). A unit closes to
// 1e-10 on the inlets it was calculated from. The whole closes to what the
// tears left open, within what rounding and the units' own tolerance allow.
fn (s Solver) solution(plan Plan) !Solution {
	mut worst_comp, mut worst_energy := 0.0, 0.0
	for r in s.runs {
		b := r.balance()
		worst_comp = math.max(worst_comp, b.component)
		worst_energy = math.max(worst_energy, b.energy)
		if b.component > unit_tolerance || b.energy > unit_tolerance {
			return FlowsheetError{
				kind:    .unbalanced
				unit:    r.name
				reason:  'it closes to ${b.component} of its flow and ${b.energy} of F·RT'
				history: s.history
			}
		}
	}
	feeds := s.sheet.feeds.map(it.stream)
	products := s.sheet.products()
	outs := products.map(s.streams[it])
	mut flow := 0.0
	mut t_flow := 0.0
	for x in feeds {
		flow += x.total()
		t_flow += x.total() * x.t
	}
	rt := thermo.r_gas * t_flow / flow
	open_flow := s.gap_flow + round_off * flow * f64(s.runs.len)
	open_energy := s.gap_energy + unit_tolerance * flow * rt * f64(s.runs.len)
	mut q := 0.0
	for r in s.runs {
		q += added(r.energy)
	}
	comp, energy := imbalance(feeds, outs, q, s.nc)
	if comp > open_flow || energy > open_energy {
		return FlowsheetError{
			kind:    .unbalanced
			reason:  'the flowsheet closes to ${comp} mol/s and ${energy} W, with ${open_flow} mol/s and ${open_energy} W left open by its tears'
			history: s.history
		}
	}
	return Solution{
		streams:  s.streams
		units:    s.runs
		products: products
		plan:     plan
		history:  s.history
		balance:  Balance{comp / flow, energy / (flow * rt)}
		worst:    Balance{worst_comp, worst_energy}
	}
}
