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

// energy_ports names its energy stream, for an expression (M8's D3).
pub fn (d FlashDrum) energy_ports() []string {
	return ['duty']
}

// of_phases is some of a flash's phases as a stream of their own: at the
// flash's T and P, so it needs no flash. It keeps every phase, those not
// taken at no moles, as a bubble or dew point's answer keeps the other
// phase, and those taken in their shares; its β is the vapour's share.
fn of_phases(fl thermo.Fluid, total f64, r flash.FlashResult, take []int) Stream {
	mut share := 0.0
	for k in take {
		share += r.phases[k].fraction
	}
	nc := r.z.len
	mut z := []f64{len: nc}
	mut phases := []flash.Phase{cap: r.phases.len}
	mut beta := 0.0
	for k, ph in r.phases {
		f := if k in take { ph.fraction / share } else { 0.0 }
		for i in 0 .. nc {
			z[i] += f * ph.state.n[i]
		}
		if ph.label == .vapor {
			beta = f
		}
		phases << flash.Phase{
			...ph
			fraction: f
		}
	}
	saturated := flash.FlashResult{
		...r
		z:      z
		beta:   beta
		phases: phases
	}
	return of_flash(fl, z.map(it * total * share), saturated)
}

// separated is the flash's vapour and its liquids, by their labels (M6's
// D9 and D11): the vapour, if any, out of the vapour port, and every liquid,
// liquid2 too, out of the liquid port. A port with no phase gets an empty
// stream.
fn separated(fl thermo.Fluid, flows []f64, total f64, r flash.FlashResult) (Stream, Stream) {
	vapour := []int{len: r.phases.len, init: index}.filter(r.phases[it].label == .vapor)
	liquids := []int{len: r.phases.len, init: index}.filter(r.phases[it].label != .vapor)
	if vapour.len == 0 || liquids.len == 0 {
		full := of_flash(fl, flows, r)
		empty := of_flash(fl, []f64{len: flows.len}, r)
		return if vapour.len > 0 { full, empty } else { empty, full }
	}
	return of_phases(fl, total, r, vapour), of_phases(fl, total, r, liquids)
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
