module flowsheet

import flash
import thermo

// Compressor raises a gas's pressure (D4). A PS flash finds the isentropic
// outlet, and the real one takes ΔH_s/η.
pub struct Compressor {
pub:
	name       string
	mode       PressureMode
	value      f64
	efficiency f64
}

// Expander lowers a gas's pressure and takes work from it (D4). A PS flash
// finds the isentropic outlet, and the real one gives up η·ΔH_s, so its
// outlet can be two-phase.
pub struct Expander {
pub:
	name       string
	mode       PressureMode
	value      f64
	efficiency f64
}

// new_compressor is a compressor with its specification.
pub fn new_compressor(name string, params MachineParams) !Compressor {
	check_machine(name, params)!
	return Compressor{
		name:       name
		mode:       params.mode
		value:      params.value
		efficiency: params.efficiency
	}
}

// new_expander is an expander with its specification.
pub fn new_expander(name string, params MachineParams) !Expander {
	check_machine(name, params)!
	return Expander{
		name:       name
		mode:       params.mode
		value:      params.value
		efficiency: params.efficiency
	}
}

// inlet_ports is in.
pub fn (u Compressor) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (u Compressor) outlet_ports() []string {
	return ['out']
}

// inlet_ports is in.
pub fn (u Expander) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (u Expander) outlet_ports() []string {
	return ['out']
}

// isentropic is the inlet taken to P at constant S.
fn isentropic(unit string, fl thermo.Fluid, inlet Stream, p f64) !Stream {
	r := flash.ps(fl, p, inlet.s, inlet.result.z, estimate: inlet.t) or {
		return flash_failed(unit, err)
	}
	return of_flash(fl, inlet.flows, r)
}

// calculate compresses the inlet.
pub fn (u Compressor) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(u.name, fl, inlets, 1)!
	inlet := inlets[0]
	p := pressure_up(u.name, u.mode, u.value, inlet.p)!
	ideal := isentropic(u.name, fl, inlet, p)!
	return worked(u.name, fl, inlet, p, inlet.h + (ideal.h - inlet.h) / u.efficiency, ideal.t)
}

// calculate expands the inlet.
pub fn (u Expander) calculate(fl thermo.Fluid, inlets []Stream) !Outcome {
	check_inlets(u.name, fl, inlets, 1)!
	inlet := inlets[0]
	p := pressure_down(u.name, u.mode, u.value, inlet.p)!
	ideal := isentropic(u.name, fl, inlet, p)!
	return worked(u.name, fl, inlet, p, inlet.h + u.efficiency * (ideal.h - inlet.h), ideal.t)
}
