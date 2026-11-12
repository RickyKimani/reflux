module reactors

import equilibrium
import flowsheet
import kinetics
import thermo

// EquilibriumReactor brings chosen reactions to equilibrium (D10, as
// amended in step 7): each one's Q its K, K from the Gibbs energies (D7),
// with every species no reaction touches passing through. It's restricted
// equilibrium, for where only some reactions run, and the reactions must be
// independent.
pub struct EquilibriumReactor {
pub:
	name      string
	species   []kinetics.Species
	reactions []kinetics.Reaction
	thermal   Thermal
	t         f64 // K, if isothermal: 0 for the inlet's
	duty      f64 // W, for .duty
	drop      f64 // Pa
}

// EquilibriumParams are an equilibrium reactor's specification.
@[params]
pub struct EquilibriumParams {
pub:
	species   []kinetics.Species  @[required]
	reactions []kinetics.Reaction @[required]
	thermal   Thermal
	t         f64
	duty      f64
	drop      f64
}

// new_equilibrium_reactor is an equilibrium reactor with its specification.
pub fn new_equilibrium_reactor(name string, p EquilibriumParams) !EquilibriumReactor {
	check_thermal(name, p.thermal, p.t, p.duty, p.drop)!
	if p.reactions.len == 0 {
		return spec_error(name, 'no reactions to bring to equilibrium')
	}
	d := kinetics.dependence(p.reactions, p.species.len) or { return spec_error(name, err.msg()) }
	if d.combinations.len > 0 {
		return spec_error(name, 'its reactions must be independent, but ${d.combinations.map(it.describe(p.reactions)).join('; ')}')
	}
	return EquilibriumReactor{
		name:      name
		species:   p.species
		reactions: p.reactions
		thermal:   p.thermal
		t:         p.t
		duty:      p.duty
		drop:      p.drop
	}
}

// inlet_ports is in.
pub fn (e EquilibriumReactor) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (e EquilibriumReactor) outlet_ports() []string {
	return ['out']
}

// calculate finds the equilibrium over the reactions, at T or by the energy
// balance.
pub fn (e EquilibriumReactor) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	inlet := one_inlet(e.name, fl, inlets, e.species.len)!
	at := EquilibriumAt{e, fl, inlet.flows, pressure_out(e.name, inlet, e.drop)!}
	return thermal_outcome(e.name, at, inlet, e.thermal, e.t, e.duty, element_rows(e.species))
}

// EquilibriumAt is an equilibrium reactor's outlet at a temperature.
struct EquilibriumAt {
	e     EquilibriumReactor
	fl    thermo.Fluid
	flows []f64
	p     f64
}

fn (a EquilibriumAt) outlet_at(t f64) !flowsheet.Stream {
	m := equilibrium.Mixture{
		species: a.e.species
		model:   a.fl.model
	}
	ln_k := equilibrium.ln_k_gibbs(a.e.reactions, a.e.species, t) or {
		return failed(a.e.name, err.msg())
	}
	r := equilibrium.over_reactions(m, a.e.reactions, ln_k, t, a.p, a.flows) or {
		return failed(a.e.name, err.msg())
	}
	return at_tp(a.e.name, a.fl, r.n, t, a.p)
}
