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

// calculate mixes the inlets. An empty inlet adds nothing, and doesn't set
// the pressure.
pub fn (m Mixer) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(m.name, fl, inlets, m.inlets)!
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
		return bad_spec(m.name, 'no inlet has any flow')
	}
	p := if m.outlet_p > 0.0 { m.outlet_p } else { lowest }
	if p > lowest {
		return bad_spec(m.name, 'an outlet at ${p} Pa is above the lowest inlet, ${lowest} Pa')
	}
	r := flash.ph(fl, p, enthalpy / total, flows, estimate: t_flow / total) or {
		return flash_failed(m.name, err)
	}
	return Outcome{
		outlets: [of_flash(fl, flows, r)]
	}
}
