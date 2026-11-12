module flowsheet

import thermo

// EnergyStream is a rate of heat or work (D2), reported by the unit that
// takes it. A rate above zero goes into the unit's material, so each unit
// balances as ΣH_in + Σrate = ΣH_out.
pub struct EnergyStream {
pub:
	name string // 'duty' for heat
	rate f64    // W
}

// Outcome is what a unit calculated.
pub struct Outcome {
pub:
	outlets []Stream // in the order of outlet_ports
	energy  []EnergyStream
	// conserved is what a unit that reacts conserves, as rows over the
	// components (its elements), on which its balance is checked: empty for
	// a unit that conserves every component, as M4's do.
	conserved [][]f64
	// profile is a reactor's along its volume, if it has one.
	profile ?Profile
}

// Profile is a reactor's state along it: the inlet, then each point its
// integrator stepped to (M9's report format 2).
pub struct Profile {
pub:
	points []ProfilePoint
}

// ProfilePoint is a reactor's state at one place along it.
pub struct ProfilePoint {
pub:
	volume f64   // m³ from the inlet
	length f64   // m from the inlet, for a reactor of known length; else 0
	t      f64   // K
	p      f64   // Pa
	flows  []f64 // mol/s of each component
}

// Unit is one unit operation. It holds its specification and nothing about
// a solve, so calculating it twice from the same inlets gives the same
// answer. Later unit families (columns, reactors) implement this from their
// own modules (D1).
pub interface Unit {
	name string
	// inlet_ports names its inlets, in the order calculate takes them.
	inlet_ports() []string
	// outlet_ports names its outlets, in the order calculate returns them.
	outlet_ports() []string
	// calculate finds the outlets from the inlets. Failure is a FlowsheetError.
	calculate(fl thermo.Fluid, inlets []Stream) !Outcome
}

// check_inlets checks that a unit was given the streams it has ports for,
// each with a flow per component of the fluid.
fn check_inlets(unit string, fl thermo.Fluid, inlets []Stream, want int) ! {
	if inlets.len != want {
		return bad_spec(unit, 'got ${inlets.len} inlets for ${want} ports')
	}
	nc := fl.model.components().len
	for s in inlets {
		if s.flows.len != nc {
			return bad_spec(unit, 'got an inlet of ${s.flows.len} components for a fluid of ${nc}')
		}
	}
}
