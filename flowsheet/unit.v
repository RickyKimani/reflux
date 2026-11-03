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
