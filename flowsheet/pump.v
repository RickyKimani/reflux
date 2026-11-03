module flowsheet

import thermo

// Pump raises a liquid's pressure (D4): ΔH = v·ΔP/η, with v the inlet's
// molar volume. It takes no vapour.
pub struct Pump {
pub:
	name       string
	mode       PressureMode
	value      f64
	efficiency f64
}

// new_pump is a pump with its specification.
pub fn new_pump(name string, params MachineParams) !Pump {
	check_machine(name, params)!
	return Pump{
		name:       name
		mode:       params.mode
		value:      params.value
		efficiency: params.efficiency
	}
}

// inlet_ports is in.
pub fn (u Pump) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (u Pump) outlet_ports() []string {
	return ['out']
}

// calculate pumps the inlet. An empty inlet has no vapour to refuse.
pub fn (u Pump) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(u.name, fl, inlets, 1)!
	inlet := inlets[0]
	if inlet.beta > 0.0 && inlet.total() > 0.0 {
		return bad_spec(u.name, 'a pump takes liquid, and this inlet has β = ${inlet.beta}')
	}
	p := pressure_up(u.name, u.mode, u.value, inlet.p)!
	h := inlet.h + inlet.v * (p - inlet.p) / u.efficiency
	return worked(u.name, fl, inlet, p, h, inlet.t)
}
