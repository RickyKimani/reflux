module flowsheet

import flash
import math
import numeric
import thermo

// ComponentSeparator sends a given fraction of each component to its first
// outlet and the rest to its second, both at a given T and P (D4). It's a
// black box, so its duty is what the balance needs.
pub struct ComponentSeparator {
pub:
	name      string
	fractions []f64 // of each component, to the first outlet
	t         f64   // K
	p         f64   // Pa
}

// SeparatorParams are a component separator's specification.
@[params]
pub struct SeparatorParams {
pub:
	fractions []f64 @[required]
	t         f64   @[required]
	p         f64   @[required]
}

// new_component_separator is a component separator with its specification.
pub fn new_component_separator(name string, params SeparatorParams) !ComponentSeparator {
	for f in params.fractions {
		if !(f >= 0.0 && f <= 1.0) {
			return bad_spec(name, 'a fraction of ${f} to the first outlet')
		}
	}
	for x in [params.t, params.p] {
		if !(x > 0.0) || !math.is_finite(x) {
			return bad_spec(name, 'an outlet at ${params.t} K and ${params.p} Pa')
		}
	}
	return ComponentSeparator{
		name:      name
		fractions: params.fractions.clone()
		t:         params.t
		p:         params.p
	}
}

// inlet_ports is in.
pub fn (s ComponentSeparator) inlet_ports() []string {
	return ['in']
}

// outlet_ports are out1 and out2.
pub fn (s ComponentSeparator) outlet_ports() []string {
	return ['out1', 'out2']
}

// calculate separates the inlet.
pub fn (s ComponentSeparator) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(s.name, fl, inlets, 1)!
	inlet := inlets[0]
	if s.fractions.len != inlet.flows.len {
		return bad_spec(s.name, 'got ${s.fractions.len} fractions for ${inlet.flows.len} components')
	}
	first := []f64{len: inlet.flows.len, init: inlet.flows[index] * s.fractions[index]}
	second := []f64{len: inlet.flows.len, init: inlet.flows[index] - first[index]}
	a := s.outlet(fl, first, inlet)!
	b := s.outlet(fl, second, inlet)!
	out := a.total() * a.h + b.total() * b.h
	return Outcome{
		outlets: [a, b]
		energy:  [EnergyStream{'duty', out - inlet.total() * inlet.h}]
	}
}

// outlet flashes one outlet at the given T and P. An empty one has the
// inlet's composition.
fn (s ComponentSeparator) outlet(fl thermo.Fluid, flows []f64, inlet Stream) !Stream {
	z := if numeric.sum(flows) > 0.0 { flows } else { inlet.result.z }
	r := flash.pt(fl.model, s.t, s.p, z) or { return flash_failed(s.name, err) }
	return of_flash(fl, flows, r)
}
