module flowsheet

import flash
import thermo

// Valve expands its inlet at constant enthalpy (D4).
pub struct Valve {
pub:
	name  string
	mode  PressureMode
	value f64
}

// ValveParams are a valve's specification.
@[params]
pub struct ValveParams {
pub:
	mode  PressureMode @[required]
	value f64          @[required]
}

// new_valve is a valve with its specification.
pub fn new_valve(name string, params ValveParams) !Valve {
	check_pressure(name, params.mode, params.value)!
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
	p := pressure_down(v.name, v.mode, v.value, inlet.p)!
	r := flash.ph(fl, p, inlet.h, inlet.result.z, estimate: inlet.t) or {
		return flash_failed(v.name, err)
	}
	return Outcome{
		outlets: [of_flash(fl, inlet.flows, r)]
	}
}
