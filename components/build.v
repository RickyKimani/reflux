module components

import tdep
import thermo

// Turning a record into what the rest of reflux uses: a thermo.Component,
// a thermo.IdealGas, and property correlations in SI per mol.
//
// The importer kept ChemSep's defects verbatim and listed them in
// data/components/CHANGES.md (M2 step 3). This is where they're decided:
// - a correlation with no coefficients, or an impossible range, is refused;
// - a vapour pressure that can't be used falls back to the Antoine record,
//   which every ChemSep compound has, and says so;
// - a compound without Tc, Pc, ω or M can't become a thermo.Component, and
//   the error says which it lacks.

struct Unit {
	scale f64    // SI value = scale × the source's value
	si    string // the SI unit
}

// The output units ChemSep's correlations use, and their SI per mol.
const si_units = {
	'Pa':       Unit{1.0, 'Pa'}
	'J/kmol/K': Unit{1e-3, 'J/(mol·K)'}
	'J/kmol':   Unit{1e-3, 'J/mol'}
	'kmol/m3':  Unit{1e3, 'mol/m³'}
	'm3/kmol':  Unit{1e-3, 'm³/mol'}
	'Pa.s':     Unit{1.0, 'Pa·s'}
	'W/m/K':    Unit{1.0, 'W/(m·K)'}
	'N/m':      Unit{1.0, 'N/m'}
	'_':        Unit{1.0, ''}
}

// Property is one of a compound's correlations, ready to evaluate in SI per
// mol: value = scale × c.eval(T), in `unit`.
pub struct Property {
pub:
	key   string // the correlation used, e.g. 'antoine_vapor_pressure'
	note  string // why, when it isn't the one asked for
	c     tdep.Correlation
	scale f64
	unit  string
}

// Value is a property's value and whether T was inside the correlation's
// fitted range. Outside it the value is an extrapolation (M2 decision D7):
// a flash's trial temperatures stray outside briefly, so it's computed
// anyway, and the caller decides whether to trust it.
pub struct Value {
pub:
	value    f64
	in_range bool
}

// at evaluates the property at T.
pub fn (p Property) at(t f64) Value {
	return Value{
		value:    p.scale * p.c.eval(t)
		in_range: p.c.in_range(t)
	}
}

// property returns correlation `key`, such as 'vapor_pressure' or
// 'liquid_density', converted to SI per mol.
pub fn (r Record) property(key string) !Property {
	return r.usable(key) or {
		if key != 'vapor_pressure' {
			return err
		}
		why := err.msg()
		p := r.usable('antoine_vapor_pressure') or { return error('${why}, and ${err}') }
		Property{
			...p
			note: 'the Antoine vapour pressure, because ${why}'
		}
	}
}

fn (r Record) usable(key string) !Property {
	what := '${r.name} (${r.ids.cas}) ${key}'
	d := r.correlations.get(key) or { return error('${r.name} (${r.ids.cas}) has ${err}') }
	if d.coeffs.len == 0 {
		return error('${what} has no coefficients (see data/components/CHANGES.md)')
	}
	if d.tmin > d.tmax || d.tmin < 0.0 || d.tmax <= 0.0 {
		return error('${what} has an impossible range [${d.tmin}, ${d.tmax}] K (see data/components/CHANGES.md)')
	}
	if !tdep.is_supported(d.eq) {
		return error("${what} is equation ${d.eq}, which reflux can't evaluate")
	}
	u := si_units[d.units] or {
		return error("${what} is in `${d.units}`, which reflux can't convert to SI")
	}
	// The reduced-temperature forms take the correlation's own Tc if it has
	// one, the compound's otherwise; the rest ignore it.
	tc := d.tc or { r.constants.tc or { 0.0 } }
	c := tdep.new(d.eq, d.coeffs, d.tmin, d.tmax, tc) or { return error('${what}: ${err}') }
	return Property{
		key:   key
		c:     c
		scale: u.scale
		unit:  u.si
	}
}

// component returns the constants every thermo model needs.
pub fn (r Record) component() !thermo.Component {
	c := r.constants
	mut missing := []string{}
	if c.tc == none {
		missing << 'Tc'
	}
	if c.pc == none {
		missing << 'Pc'
	}
	if c.omega == none {
		missing << 'ω'
	}
	if c.mw == none {
		missing << 'M'
	}
	if missing.len > 0 {
		return error('${r.name} (${r.ids.cas}) has no ${missing.join(', ')}, so no equation of state can use it')
	}
	return thermo.Component{
		name:  r.name
		tc:    c.tc or { 0.0 }
		pc:    c.pc or { 0.0 }
		omega: c.omega or { 0.0 }
		mw:    c.mw or { 0.0 }
	}
}

// ideal_gas returns the compound's ideal-gas Cp, from its ideal_gas_cp
// correlation, with its ΔH°f and S° when it has them.
pub fn (r Record) ideal_gas() !thermo.IdealGas {
	p := r.property('ideal_gas_cp')!
	if p.unit != 'J/(mol·K)' {
		return error('${r.name} (${r.ids.cas}) ideal_gas_cp is in ${p.unit}, not J/(mol·K)')
	}
	return thermo.IdealGas{
		cp:  thermo.Correlated{
			c:     p.c
			scale: p.scale
		}
		dhf: r.constants.dhf
		s0:  r.constants.s0
	}
}
