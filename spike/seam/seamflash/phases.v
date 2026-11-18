module seamflash

import thermo

// Phases is M6's D3 in its least form: what the PT flash and the stability
// test ask of a model, a phase at (T, P, n) with its ln φ, their derivatives
// and its label. Not part of the simulator.
pub interface Phases {
	components() []thermo.Component
	phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State
	derivs(st thermo.State) thermo.LnPhiDerivs
	label(st thermo.State) thermo.PhaseLabel
}

// Eos is a Helmholtz model as a Phases, by its volume roots: exactly what
// flash calls on a thermo.Model, one call further away.
pub struct Eos {
pub:
	m thermo.Model
}

pub fn (e Eos) components() []thermo.Component {
	return e.m.components()
}

pub fn (e Eos) phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State {
	return thermo.solve_tp(e.m, t, p, n, root)
}

pub fn (e Eos) derivs(st thermo.State) thermo.LnPhiDerivs {
	return thermo.ln_phi_derivs(e.m, st)
}

pub fn (e Eos) label(st thermo.State) thermo.PhaseLabel {
	return thermo.phase_label(e.m, st)
}
