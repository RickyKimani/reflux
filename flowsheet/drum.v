module flowsheet

import flash
import math
import thermo

// DrumMode says what a flash drum's value is.
pub enum DrumMode {
	outlet_t // the drum's temperature, K
	duty     // the heat added, W: 0 is adiabatic
}

// FlashDrum mixes its inlets and flashes the mixture, with the vapour and
// the liquid leaving by separate outlets (D4).
pub struct FlashDrum {
pub:
	name   string
	inlets int
	mode   DrumMode
	value  f64
	// outlet_p is the drum's pressure, Pa, or 0 for the lowest inlet's.
	outlet_p f64
}

// DrumParams are a flash drum's specification.
@[params]
pub struct DrumParams {
pub:
	inlets int = 1
	mode   DrumMode @[required]
	value  f64      @[required]
	// outlet_p is the drum's pressure, Pa, at most the lowest inlet's. Left
	// at 0, it's the lowest inlet's.
	outlet_p f64
}

// new_flash_drum is a flash drum with its specification.
pub fn new_flash_drum(name string, params DrumParams) !FlashDrum {
	if params.inlets < 1 {
		return bad_spec(name, 'a flash drum needs an inlet, and has ${params.inlets}')
	}
	ok := match params.mode {
		.outlet_t { params.value > 0.0 }
		.duty { true }
	}
	if !ok || !math.is_finite(params.value) {
		return bad_spec(name, 'a flash drum with ${params.mode} ${params.value}')
	}
	if !math.is_finite(params.outlet_p) || params.outlet_p < 0.0 {
		return bad_spec(name, 'the outlet pressure is ${params.outlet_p} Pa')
	}
	return FlashDrum{
		name:     name
		inlets:   params.inlets
		mode:     params.mode
		value:    params.value
		outlet_p: params.outlet_p
	}
}

// inlet_ports are in1, in2, and so on.
pub fn (d FlashDrum) inlet_ports() []string {
	return []string{len: d.inlets, init: 'in${index + 1}'}
}

// outlet_ports are vapour and liquid.
pub fn (d FlashDrum) outlet_ports() []string {
	return ['vapour', 'liquid']
}

// of_phase is phase k of a two-phase flash as a stream of its own: saturated
// at the flash's T and P, so it needs no flash. It has the form a bubble or
// dew point's answer has, both phases with the other at no moles.
fn of_phase(fl thermo.Fluid, total f64, r flash.FlashResult, k int) Stream {
	phase := r.phases[k]
	mut phases := []flash.Phase{cap: r.phases.len}
	for i, other in r.phases {
		phases << flash.Phase{
			...other
			fraction: if i == k { 1.0 } else { 0.0 }
		}
	}
	saturated := flash.FlashResult{
		...r
		z:      phase.state.n
		beta:   if k == 0 { 1.0 } else { 0.0 }
		phases: phases
	}
	return of_flash(fl, phase.state.n.map(it * total * phase.fraction), saturated)
}

// separated is the flash's vapour and liquid, or for a single phase, that
// phase and an empty stream.
fn separated(fl thermo.Fluid, flows []f64, total f64, r flash.FlashResult) (Stream, Stream) {
	if r.phases.len == 2 {
		return of_phase(fl, total, r, 0), of_phase(fl, total, r, 1)
	}
	full := of_flash(fl, flows, r)
	empty := of_flash(fl, []f64{len: flows.len}, r)
	return if r.beta == 1.0 { full, empty } else { empty, full }
}

// calculate mixes the inlets and flashes them. An empty inlet adds nothing
// and doesn't set the pressure, as in a mixer. The duty is what the balance
// needs, and with a duty given that's the duty the flash reached.
pub fn (d FlashDrum) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	b := blend(d.name, fl, inlets, d.inlets)!
	p := b.outlet_pressure(d.name, d.outlet_p)!
	r := d.flash_mixture(fl, b, p) or { return flash_failed(d.name, err) }
	vapour, liquid := separated(fl, b.flows, b.total, r)
	out := vapour.total() * vapour.h + liquid.total() * liquid.h
	return Outcome{
		outlets: [vapour, liquid]
		energy:  [EnergyStream{'duty', out - b.enthalpy}]
	}
}

// flash_mixture is the flash that the mode calls for.
fn (d FlashDrum) flash_mixture(fl thermo.Fluid, b Blend, p f64) !flash.FlashResult {
	return match d.mode {
		.outlet_t { flash.pt(fl.model, d.value, p, b.flows)! }
		.duty { flash.ph(fl, p, (b.enthalpy + d.value) / b.total, b.flows, estimate: b.t)! }
	}
}
