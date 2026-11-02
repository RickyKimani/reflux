module flowsheet

import flash
import math
import thermo

// ValveMode says what a valve's value is.
pub enum ValveMode {
	outlet_p      // the outlet pressure, Pa
	pressure_drop // the drop across it, Pa, at least 0
}

// Valve expands its inlet at constant enthalpy (D4).
pub struct Valve {
pub:
	name  string
	mode  ValveMode
	value f64
}

// ValveParams are a valve's specification.
@[params]
pub struct ValveParams {
pub:
	mode  ValveMode @[required]
	value f64       @[required]
}

// new_valve is a valve with its specification.
pub fn new_valve(name string, params ValveParams) !Valve {
	ok := match params.mode {
		.outlet_p { params.value > 0.0 }
		.pressure_drop { params.value >= 0.0 }
	}
	if !ok || !math.is_finite(params.value) {
		return bad_spec(name, 'a valve with ${params.mode} ${params.value} Pa')
	}
	return Valve{
		name:  name
		mode:  params.mode
		value: params.value
	}
}

// inlet_ports is in.
pub fn (v Valve) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (v Valve) outlet_ports() []string {
	return ['out']
}

// calculate expands the inlet to the outlet pressure, by a PH flash from
// the inlet's T.
pub fn (v Valve) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(v.name, fl, inlets, 1)!
	inlet := inlets[0]
	p := match v.mode {
		.outlet_p { v.value }
		.pressure_drop { inlet.p - v.value }
	}
	if !(p > 0.0) {
		return bad_spec(v.name, 'a drop of ${v.value} Pa from ${inlet.p} Pa leaves no pressure')
	}
	if p > inlet.p {
		return bad_spec(v.name, "a valve can't raise ${inlet.p} Pa to ${p} Pa")
	}
	r := flash.ph(fl, p, inlet.h, inlet.result.z, estimate: inlet.t) or {
		return flash_failed(v.name, err)
	}
	return Outcome{
		outlets: [of_flash(fl, inlet.flows, r)]
	}
}
