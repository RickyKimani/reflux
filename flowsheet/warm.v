module flowsheet

import math
import thermo

// Warm is a unit that can start from its last answer at a neighbouring
// point (M8's D6): calculate_from gives what calculate would, to its own
// solver's tolerance, starting from last, the unit as it was calculated
// there. A unit without it starts cold.
pub interface Warm {
	calculate_from(fl thermo.Fluid, inlets []Stream, last UnitResult) !Outcome
}

// Start is where a solve starts from (M8's D6): the solution of a
// neighbouring point, and the one before it for a secant predictor. Each
// part can be left out, for measuring what it saves.
pub struct Start {
pub:
	last   Solution
	before ?Solution
	tears  bool = true // a block's tears from the last answer
	// predict starts the tears and the specs' values from the secant through
	// the last two answers, b + ratio·(b − a), where there are two: ratio is
	// this point's step over the last one's.
	predict bool = true
	ratio   f64  = 1.0
	units   bool = true // each Warm unit from its last answer
	specs   bool = true // each spec from its last value
}

// last_unit is a unit's result at the last point, if the solve starts from
// one and the unit was calculated there.
fn (s Solver) last_unit(i int) ?UnitResult {
	start := s.start or { return none }
	if !start.units {
		return none
	}
	return start.last.units[i] or { return none }
}

// calculated is a unit's outcome on these inlets: from its last answer
// where it can start from one, else cold.
fn (s Solver) calculated(i int, u Unit, inlets []Stream) !Outcome {
	if last := s.last_unit(i) {
		if u is Warm {
			return u.calculate_from(s.fl, inlets, last)
		}
	}
	return u.calculate(s.fl, inlets)
}

// warm_tear is a tear's first value from the last point's answer, or from
// the secant through the last two, b + r·(b − a) in its flows, P and molar
// H, a flow held at 0, by a PH flash from b's T. A predictor the flash
// refuses falls back to b.
fn (s Solver) warm_tear(name string) ?Stream {
	start := s.start or { return none }
	// Measured in M8 step 5: only direct substitution, and Broyden with
	// specs, take fewer passes from a neighbour's tears; Wegstein and Broyden
	// alone take more than from their cold start.
	if !start.tears || !(s.params.method == .direct || s.block.specs.len > 0) {
		return none
	}
	b := start.last.streams[name] or { return none }
	if !start.predict {
		return b
	}
	before := start.before or { return b }
	a := before.streams[name] or { return b }
	mut flows := []f64{len: b.flows.len}
	r := start.ratio
	for i in 0 .. flows.len {
		flows[i] = math.max(b.flows[i] + r * (b.flows[i] - a.flows[i]), 0)
	}
	p := b.p + r * (b.p - a.p)
	h := b.h + r * (b.h - a.h)
	if !(p > 0) || flows.all(it == 0) {
		return b
	}
	return feed_ph_from(s.fl, flows, p, h, b.t) or { b }
}

// warm_value is a spec's first value: the last point's, or the secant
// through the last two, held within its bounds.
fn (s Solver) warm_value(k int, lo f64, hi f64) ?f64 {
	start := s.start or { return none }
	if !start.specs {
		return none
	}
	b := start.last.specs[k] or { return none }
	if !start.predict {
		return b.value
	}
	before := start.before or { return b.value }
	a := before.specs[k] or { return b.value }
	return math.max(lo, math.min(hi, b.value + start.ratio * (b.value - a.value)))
}
