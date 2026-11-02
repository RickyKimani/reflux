module flowsheet

import flash
import math
import thermo

// PressureMode says what a pressure-changing unit's value is.
pub enum PressureMode {
	outlet_p // the outlet pressure, Pa
	delta_p  // the size of the change, Pa, at least 0
}

// MachineParams are a pump's, compressor's or expander's specification.
@[params]
pub struct MachineParams {
pub:
	mode       PressureMode @[required]
	value      f64          @[required]
	efficiency f64          @[required] // above 0, at most 1
}

// check_pressure checks a pressure specification.
fn check_pressure(unit string, mode PressureMode, value f64) ! {
	ok := match mode {
		.outlet_p { value > 0.0 }
		.delta_p { value >= 0.0 }
	}
	if !ok || !math.is_finite(value) {
		return bad_spec(unit, 'a pressure specification of ${mode} ${value} Pa')
	}
}

// check_machine checks a machine's specification.
fn check_machine(unit string, params MachineParams) ! {
	check_pressure(unit, params.mode, params.value)!
	if !(params.efficiency > 0.0 && params.efficiency <= 1.0) {
		return bad_spec(unit, 'an efficiency of ${params.efficiency}')
	}
}

// pressure_up is the outlet pressure of a unit that raises it.
fn pressure_up(unit string, mode PressureMode, value f64, inlet f64) !f64 {
	p := match mode {
		.outlet_p { value }
		.delta_p { inlet + value }
	}
	if p < inlet {
		return bad_spec(unit, "it can't lower ${inlet} Pa to ${p} Pa")
	}
	return p
}

// pressure_down is the outlet pressure of a unit that lowers it.
fn pressure_down(unit string, mode PressureMode, value f64, inlet f64) !f64 {
	p := match mode {
		.outlet_p { value }
		.delta_p { inlet - value }
	}
	if !(p > 0.0) {
		return bad_spec(unit, 'a drop of ${value} Pa from ${inlet} Pa leaves no pressure')
	}
	if p > inlet {
		return bad_spec(unit, "it can't raise ${inlet} Pa to ${p} Pa")
	}
	return p
}

// worked is a machine's outcome: the inlet taken to P and molar H by a PH
// flash from T, with the work that took, F·(H_out − H_in).
fn worked(unit string, fl thermo.Fluid, inlet Stream, p f64, h f64, t f64) !Outcome {
	r := flash.ph(fl, p, h, inlet.result.z, estimate: t) or { return flash_failed(unit, err) }
	outlet := of_flash(fl, inlet.flows, r)
	return Outcome{
		outlets: [outlet]
		energy:  [EnergyStream{'work', inlet.total() * (outlet.h - inlet.h)}]
	}
}
