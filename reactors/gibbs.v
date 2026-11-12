module reactors

import equilibrium
import flowsheet
import kinetics
import math
import thermo

// Gibbs is the Gibbs reactor (D9, D10): its outlet the equilibrium over every
// species of the file that its feed's elements allow, in one phase, at a
// temperature or where its energy balance puts it. Species held inert pass
// through.
pub struct Gibbs {
pub:
	name    string
	species []kinetics.Species
	inert   []int
	thermal Thermal
	t       f64 // K, if isothermal: 0 for the inlet's
	duty    f64 // W, for .duty
	drop    f64 // Pa
}

// GibbsParams are a Gibbs reactor's specification.
@[params]
pub struct GibbsParams {
pub:
	species []kinetics.Species @[required]
	inert   []int
	thermal Thermal
	t       f64
	duty    f64
	drop    f64
}

// new_gibbs is a Gibbs reactor with its specification.
pub fn new_gibbs(name string, p GibbsParams) !Gibbs {
	check_thermal(name, p.thermal, p.t, p.duty, p.drop)!
	if p.species.len == 0 {
		return spec_error(name, 'no species')
	}
	for i in p.inert {
		if i < 0 || i >= p.species.len {
			return spec_error(name, 'no species ${i} to hold inert')
		}
	}
	return Gibbs{
		name:    name
		species: p.species
		inert:   p.inert
		thermal: p.thermal
		t:       p.t
		duty:    p.duty
		drop:    p.drop
	}
}

// check_thermal checks a reactor's temperature, duty and pressure drop.
fn check_thermal(name string, thermal Thermal, t f64, duty f64, drop f64) ! {
	if thermal == .isothermal && (!(t >= 0) || math.is_inf(t, 1)) {
		return spec_error(name, 'a temperature of ${t} K')
	}
	if !math.is_finite(duty) {
		return spec_error(name, 'a duty of ${duty} W')
	}
	if !(drop >= 0) || math.is_inf(drop, 1) {
		return spec_error(name, 'a pressure drop of ${drop} Pa')
	}
}

// inlet_ports is in.
pub fn (g Gibbs) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (g Gibbs) outlet_ports() []string {
	return ['out']
}

// calculate finds the equilibrium, at T or by the energy balance, and the
// duty an isothermal reactor takes, F·H out less F·H in.
pub fn (g Gibbs) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	inlet := one_inlet(g.name, fl, inlets, g.species.len)!
	at := GibbsAt{g, fl, inlet.flows, pressure_out(g.name, inlet, g.drop)!}
	return thermal_outcome(g.name, at, inlet, g.thermal, g.t, g.duty, element_rows(g.species))
}

// GibbsAt is a Gibbs reactor's outlet at a temperature, for the energy
// balance.
struct GibbsAt {
	g     Gibbs
	fl    thermo.Fluid
	flows []f64
	p     f64
}

fn (a GibbsAt) outlet_at(t f64) !flowsheet.Stream {
	m := equilibrium.Mixture{
		species: a.g.species
		model:   a.fl.model
	}
	r := equilibrium.gibbs(m, t, a.p, a.flows, a.g.inert) or { return failed(a.g.name, err.msg()) }
	return at_tp(a.g.name, a.fl, r.n, t, a.p)
}

// one_inlet checks a reactor's one inlet: a flow per species, and something
// flowing.
fn one_inlet(name string, fl thermo.Fluid, inlets []flowsheet.Stream, n int) !flowsheet.Stream {
	if inlets.len != 1 {
		return spec_error(name, 'got ${inlets.len} inlets for 1 port')
	}
	if fl.model.components().len != n || inlets[0].flows.len != n {
		return spec_error(name, 'got an inlet of ${inlets[0].flows.len} components for ${n} species')
	}
	if !(inlets[0].total() > 0) {
		return spec_error(name, 'nothing flows in')
	}
	return inlets[0]
}

// thermal_outcome is a reactor's outlet by its thermal mode, from what it
// makes at each T: at T (the inlet's if 0), or where its enthalpy is the
// inlet's (adiabatic) or the inlet's plus the duty. An isothermal reactor's
// duty is F·H out less F·H in; an adiabatic one has none, its balance
// closing to Brent's tolerance.
fn thermal_outcome(name string, at Outlet, inlet flowsheet.Stream, thermal Thermal, t f64, duty f64, rows [][]f64) !flowsheet.Outcome {
	outlet := match thermal {
		.isothermal { at.outlet_at(if t > 0 { t } else { inlet.t })! }
		.adiabatic { solve_t(name, at, enthalpy_flow(inlet), inlet.t)! }
		.duty { solve_t(name, at, enthalpy_flow(inlet) + duty, inlet.t)! }
	}
	one_phase(name, outlet)!
	energy := match thermal {
		.isothermal {
			[flowsheet.EnergyStream{'duty', enthalpy_flow(outlet) - enthalpy_flow(inlet)}]
		}
		.adiabatic { []flowsheet.EnergyStream{} }
		.duty { [flowsheet.EnergyStream{'duty', duty}] }
	}
	return flowsheet.Outcome{
		outlets:   [outlet]
		energy:    energy
		conserved: rows
	}
}
