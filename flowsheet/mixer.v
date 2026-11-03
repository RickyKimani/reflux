module flowsheet

import flash
import math
import thermo

// Mixer joins its inlets into one outlet (D4): the flows add, the enthalpy
// flows add, and a PH flash at the outlet pressure finds the state.
pub struct Mixer {
pub:
	name   string
	inlets int
	// outlet_p is the outlet pressure, Pa, or 0 for the lowest inlet's.
	outlet_p f64
}

// MixerParams are a mixer's options.
@[params]
pub struct MixerParams {
pub:
	// outlet_p is the outlet pressure, Pa, at most the lowest inlet's. Left
	// at 0, it's the lowest inlet's.
	outlet_p f64
}

// new_mixer is a mixer with this many inlets, at least 2.
pub fn new_mixer(name string, inlets int, params MixerParams) !Mixer {
	if inlets < 2 {
		return bad_spec(name, 'a mixer needs at least 2 inlets, not ${inlets}')
	}
	if !math.is_finite(params.outlet_p) || params.outlet_p < 0.0 {
		return bad_spec(name, 'the outlet pressure is ${params.outlet_p} Pa')
	}
	return Mixer{
		name:     name
		inlets:   inlets
		outlet_p: params.outlet_p
	}
}

// inlet_ports are in1, in2, and so on.
pub fn (m Mixer) inlet_ports() []string {
	return []string{len: m.inlets, init: 'in${index + 1}'}
}

// outlet_ports is out.
pub fn (m Mixer) outlet_ports() []string {
	return ['out']
}

// Blend is inlets added together, before any flash.
struct Blend {
	flows    []f64 // mol/s
	enthalpy f64   // W
	total    f64   // mol/s
	t        f64   // K, the flow-weighted mean: a start for the flash
	p        f64   // Pa, the lowest inlet's
}

// blend adds up a unit's inlets. An empty inlet adds nothing, and doesn't
// set the pressure.
fn blend(unit string, fl thermo.Fluid, inlets []Stream, want int) !Blend {
	check_inlets(unit, fl, inlets, want)!
	mut flows := []f64{len: fl.model.components().len}
	mut enthalpy := 0.0 // W
	mut t_flow := 0.0 // K·mol/s
	mut total := 0.0 // mol/s
	mut lowest := math.inf(1)
	for s in inlets {
		f := s.total()
		if f == 0.0 {
			continue
		}
		for i, n in s.flows {
			flows[i] += n
		}
		enthalpy += f * s.h
		t_flow += f * s.t
		total += f
		lowest = math.min(lowest, s.p)
	}
	if total == 0.0 {
		return bad_spec(unit, 'no inlet has any flow')
	}
	return Blend{
		flows:    flows
		enthalpy: enthalpy
		total:    total
		t:        t_flow / total
		p:        lowest
	}
}

// outlet_pressure is a given outlet pressure, or 0 for the lowest inlet's.
fn (b Blend) outlet_pressure(unit string, given f64) !f64 {
	p := if given > 0.0 { given } else { b.p }
	if p > b.p {
		return bad_spec(unit, 'an outlet at ${p} Pa is above the lowest inlet, ${b.p} Pa')
	}
	return p
}

// calculate mixes the inlets, by a PH flash from their mean T.
pub fn (m Mixer) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	b := blend(m.name, fl, inlets, m.inlets)!
	p := b.outlet_pressure(m.name, m.outlet_p)!
	r := flash.ph(fl, p, b.enthalpy / b.total, b.flows, estimate: b.t) or {
		return flash_failed(m.name, err)
	}
	return Outcome{
		outlets: [of_flash(fl, b.flows, r)]
	}
}
