module flowsheet

import flash
import math
import thermo

// HeaterMode says what a heater's value is.
pub enum HeaterMode {
	outlet_t               // the outlet temperature, K
	duty                   // the heat added, W: below 0 it's a cooler
	outlet_vapour_fraction // the outlet's β, 0 to 1
	delta_t                // the temperature rise, K: below 0 it's a cooler
}

// Heater adds or removes heat, by whichever of its modes it's given, with
// a pressure drop (D3, D4). It's the heater and the cooler both.
pub struct Heater {
pub:
	name          string
	mode          HeaterMode
	value         f64
	pressure_drop f64 // Pa
}

// HeaterParams are a heater's specification.
@[params]
pub struct HeaterParams {
pub:
	mode          HeaterMode @[required]
	value         f64        @[required]
	pressure_drop f64
}

// new_heater is a heater with its specification.
pub fn new_heater(name string, params HeaterParams) !Heater {
	ok := match params.mode {
		.outlet_t { params.value > 0.0 }
		.duty, .delta_t { true }
		.outlet_vapour_fraction { params.value >= 0.0 && params.value <= 1.0 }
	}
	if !ok || !math.is_finite(params.value) {
		return bad_spec(name, 'a heater with ${params.mode} ${params.value}')
	}
	if !(params.pressure_drop >= 0.0) || !math.is_finite(params.pressure_drop) {
		return bad_spec(name, 'the pressure drop is ${params.pressure_drop} Pa')
	}
	return Heater{
		name:          name
		mode:          params.mode
		value:         params.value
		pressure_drop: params.pressure_drop
	}
}

// inlet_ports is in.
pub fn (h Heater) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (h Heater) outlet_ports() []string {
	return ['out']
}

// calculate finds the outlet and the duty, F·(H_out − H_in). With a duty
// given, that's the duty the flash reached, within its energy tolerance.
pub fn (h Heater) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(h.name, fl, inlets, 1)!
	inlet := inlets[0]
	p := inlet.p - h.pressure_drop
	if !(p > 0.0) {
		return bad_spec(h.name, 'a drop of ${h.pressure_drop} Pa from ${inlet.p} Pa leaves no pressure')
	}
	total := inlet.total()
	if h.mode == .duty && total == 0.0 {
		return bad_spec(h.name, 'a duty on an empty stream')
	}
	r := h.flash_outlet(fl, inlet, p, total) or { return flash_failed(h.name, err) }
	outlet := of_flash(fl, inlet.flows, r)
	return Outcome{
		outlets: [outlet]
		energy:  [EnergyStream{'duty', total * (outlet.h - inlet.h)}]
	}
}

// flash_outlet is the flash that the mode calls for, at the outlet's P.
fn (h Heater) flash_outlet(fl thermo.Fluid, inlet Stream, p f64, total f64) !flash.FlashResult {
	z := inlet.result.z
	return match h.mode {
		.outlet_t { flash.pt(fl.model, h.value, p, z)! }
		.delta_t { flash.pt(fl.model, inlet.t + h.value, p, z)! }
		.duty { flash.ph(fl, p, inlet.h + h.value / total, z, estimate: inlet.t)! }
		.outlet_vapour_fraction { flash.pvf(fl.model, p, h.value, z)! }
	}
}
