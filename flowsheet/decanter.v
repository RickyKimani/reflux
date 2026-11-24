module flowsheet

import thermo

// Decanter is a three-phase separator (M6's D11): it mixes its inlets and
// flashes them at a T or a duty and a P, as the flash drum does, and each
// phase leaves by its own outlet, the vapour by vapour, the lighter liquid
// by light and the heavier by heavy. Any of them may be empty. A single
// liquid is D9's liquid, and leaves by light.
pub struct Decanter {
pub:
	name string
	spec FlashDrum // its specification, as a flash drum's
}

// new_decanter is a decanter with its specification, which is a flash
// drum's.
pub fn new_decanter(name string, params DrumParams) !Decanter {
	return Decanter{
		name: name
		spec: new_flash_drum(name, params)!
	}
}

// inlet_ports are in1, in2, and so on.
pub fn (d Decanter) inlet_ports() []string {
	return d.spec.inlet_ports()
}

// outlet_ports are vapour, light and heavy.
pub fn (d Decanter) outlet_ports() []string {
	return ['vapour', 'light', 'heavy']
}

// energy_ports names its energy stream, for an expression (M8's D3).
pub fn (d Decanter) energy_ports() []string {
	return ['duty']
}

// calculate mixes the inlets and flashes them, as a flash drum does, and
// sends each phase out by its label: vapour by vapour, liquid by light and
// liquid2 by heavy. A port with no phase gets an empty stream, at the
// flash's T and P.
pub fn (d Decanter) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	b := blend(d.name, fl, inlets, d.spec.inlets)!
	p := b.outlet_pressure(d.name, d.spec.outlet_p)!
	r := d.spec.flash_mixture(fl, b, p) or { return flash_failed(d.name, err) }
	mut outlets := []Stream{}
	mut out := 0.0
	for label in [thermo.PhaseLabel.vapor, .liquid, .liquid2] {
		k := r.phases.map(it.label).index(label)
		s := if k < 0 {
			of_flash(fl, []f64{len: b.flows.len}, r)
		} else if r.phases.len == 1 {
			of_flash(fl, b.flows, r)
		} else {
			of_phases(fl, b.total, r, [k])
		}
		out += s.total() * s.h
		outlets << s
	}
	return Outcome{
		outlets: outlets
		energy:  [EnergyStream{'duty', out - b.enthalpy}]
	}
}
