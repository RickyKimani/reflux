module reactors

import flowsheet
import kinetics
import math
import numeric
import thermo

// Conversion is one step of a conversion reactor: a reaction, the reactant
// it's measured by, and the share of what enters the step that reacts.
pub struct Conversion {
pub:
	reaction   kinetics.Reaction
	key        int // the species, by index
	conversion f64 // 0 to 1
}

// ConversionReactor runs its reactions in order (D10), each converting its
// share of its key reactant as it reaches that step, so a later reaction
// sees what the earlier ones left. The outlet is at a temperature, or where
// the energy balance puts it.
pub struct ConversionReactor {
pub:
	name    string
	species []kinetics.Species
	steps   []Conversion
	thermal Thermal
	t       f64 // K, if isothermal: 0 for the inlet's
	duty    f64 // W, for .duty
	drop    f64 // Pa
}

// ConversionParams are a conversion reactor's specification.
@[params]
pub struct ConversionParams {
pub:
	species []kinetics.Species @[required]
	steps   []Conversion       @[required]
	thermal Thermal
	t       f64
	duty    f64
	drop    f64
}

// new_conversion_reactor is a conversion reactor with its specification:
// each step's key a reactant of its reaction, and its conversion from 0 to 1.
pub fn new_conversion_reactor(name string, p ConversionParams) !ConversionReactor {
	check_thermal(name, p.thermal, p.t, p.duty, p.drop)!
	if p.steps.len == 0 {
		return spec_error(name, 'no reactions')
	}
	for s in p.steps {
		r := s.reaction
		k := r.species.index(s.key)
		if k < 0 || r.nu[k] >= 0 {
			what := if s.key >= 0 && s.key < p.species.len {
				p.species[s.key].name
			} else {
				'${s.key}'
			}
			return spec_error(name, '${what} is not a reactant of ${r.name}')
		}
		if !(s.conversion >= 0 && s.conversion <= 1) {
			return spec_error(name, "${r.name}'s conversion is ${s.conversion}, not from 0 to 1")
		}
	}
	return ConversionReactor{
		name:    name
		species: p.species
		steps:   p.steps
		thermal: p.thermal
		t:       p.t
		duty:    p.duty
		drop:    p.drop
	}
}

// inlet_ports is in.
pub fn (c ConversionReactor) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (c ConversionReactor) outlet_ports() []string {
	return ['out']
}

// flows_out runs the steps on the inlet's flows: each step's extent is its
// conversion of its key as it enters the step, over the key's coefficient.
// A step that would take a reactant below zero beyond rounding, 1e-12 of
// the flow, is refused, with how far; within it, the reactant is zero.
fn (c ConversionReactor) flows_out(flows []f64) ![]f64 {
	total := numeric.sum(flows)
	mut n := flows.clone()
	for s in c.steps {
		r := s.reaction
		before := n[s.key]
		extent := s.conversion * before / -r.nu[r.species.index(s.key)]
		for q, j in r.species {
			n[j] += r.nu[q] * extent
		}
		for j in r.species {
			if n[j] >= 0 {
				continue
			}
			if n[j] >= -1e-12 * total {
				n[j] = 0
				continue
			}
			return spec_error(c.name, '${r.name} needs more ${c.species[j].name} than there is: converting ${s.conversion} of its ${c.species[s.key].name} would leave ${n[j]} mol/s')
		}
	}
	return n
}

// calculate runs the steps, then finds the outlet at T or by the energy
// balance.
pub fn (c ConversionReactor) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	inlet := one_inlet(c.name, fl, inlets, c.species.len)!
	flows := c.flows_out(inlet.flows)!
	at := ConversionAt{c, fl, flows, pressure_out(c.name, inlet, c.drop)!}
	return thermal_outcome(c.name, at, inlet, c.thermal, c.t, c.duty, element_rows(c.species))
}

// ConversionAt is a conversion reactor's outlet at a temperature: its flows
// fixed, only its state moving.
struct ConversionAt {
	c     ConversionReactor
	fl    thermo.Fluid
	flows []f64
	p     f64
}

fn (a ConversionAt) outlet_at(t f64) !flowsheet.Stream {
	if !(math.is_finite(t) && t > 0) {
		return failed(a.c.name, 'no outlet at ${t} K')
	}
	return at_tp(a.c.name, a.fl, a.flows, t, a.p)
}
