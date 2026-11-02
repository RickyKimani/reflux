module flowsheet

import math
import numeric
import thermo

// Splitter divides its inlet by fractions of the total (D4). Each outlet
// has the inlet's intensive state, so nothing is flashed.
pub struct Splitter {
pub:
	name      string
	fractions []f64 // one per outlet, summing to 1
}

// new_splitter is a splitter with a fraction for each outlet, at least 2.
// The fractions must sum to 1 within 1e-9, a sum no typing slip can reach.
pub fn new_splitter(name string, fractions []f64) !Splitter {
	if fractions.len < 2 {
		return bad_spec(name, 'a splitter needs at least 2 outlets, not ${fractions.len}')
	}
	for f in fractions {
		if !(f >= 0.0) {
			return bad_spec(name, 'a split fraction of ${f}')
		}
	}
	sum := numeric.sum(fractions)
	if math.abs(sum - 1.0) > 1e-9 {
		return bad_spec(name, 'the split fractions sum to ${sum}, not 1')
	}
	return Splitter{
		name:      name
		fractions: fractions.map(it / sum)
	}
}

// inlet_ports is in.
pub fn (s Splitter) inlet_ports() []string {
	return ['in']
}

// outlet_ports are out1, out2, and so on.
pub fn (s Splitter) outlet_ports() []string {
	return []string{len: s.fractions.len, init: 'out${index + 1}'}
}

// calculate splits the inlet.
pub fn (s Splitter) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(s.name, fl, inlets, 1)!
	inlet := inlets[0]
	mut outlets := []Stream{cap: s.fractions.len}
	for f in s.fractions {
		outlets << Stream{
			...inlet
			flows: inlet.flows.map(it * f)
		}
	}
	return Outcome{
		outlets: outlets
	}
}
