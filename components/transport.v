module components

import math
import transport
import water

// A compound's transport data, and a fluid's transport model (M5b's D2,
// D3, D10): the transport module knows nothing of component files, so they
// become its compounds here, as activity's models are built here.

// Lucas's quantum parameter Q for the quantum gases (Poling et al., §9-4).
const lucas_quantum = {
	'7440-59-7': 1.38 // helium-4
	'1333-74-0': 0.76 // hydrogen
	'7782-39-0': 0.52 // deuterium
}

// The SI unit each quantity's correlation must come in.
fn transport_unit(q transport.Quantity) string {
	return match q {
		.liquid_viscosity, .vapour_viscosity { 'Pa·s' }
		.liquid_conductivity, .vapour_conductivity { 'W/(m·K)' }
		.surface_tension { 'N/m' }
		.liquid_density { 'mol/m³' }
		.vapour_pressure { 'Pa' }
	}
}

// transport_pure is the compound as transport needs it (D3): its Tc, Pc, ω
// and M, which every fluid's compound has, its Vc, Zc and dipole moment
// (NaN where the file has none), and each correlation, or why it has none.
// A vapour pressure that can't be used falls back to the Antoine record, as
// every property does (Record.property).
pub fn (r Record) transport_pure() !transport.Pure {
	c := r.component()!
	mut sources := []transport.Source{}
	for q in transport.quantities {
		p := r.property(q.key()) or {
			sources << transport.Source{
				why: err.msg()
			}
			continue
		}
		want := transport_unit(q)
		if p.unit != want {
			sources << transport.Source{
				why: '${r.name} (${r.ids.cas}) ${q.key()} is in ${p.unit}, not ${want}'
			}
			continue
		}
		sources << transport.Source{
			has: true
			c:   transport.Correlated{
				key:   p.key
				c:     p.c
				scale: p.scale
			}
		}
	}
	k := r.constants
	dipole := if d := k.dipole { d / transport.debye } else { math.nan() }
	return transport.Pure{
		name:    r.name
		cas:     r.ids.cas
		tc:      c.tc
		pc:      c.pc
		vc:      k.vc or { math.nan() }
		zc:      k.zc or { math.nan() }
		omega:   c.omega
		mw:      c.mw
		dipole:  dipole
		quantum: lucas_quantum[r.ids.cas] or { 0.0 }
		sources: sources
	}
}

// mixture is the compounds' transport by their correlations (D3 to D9),
// with no liquid pair's Gᵢⱼ; a file's pairs come with format 5 (D15).
pub fn mixture(records []Record) !transport.Mixture {
	mut pures := []transport.Pure{}
	for r in records {
		pures << r.transport_pure()!
	}
	return transport.Mixture{
		pures: pures
	}
}

// transport is the named compounds' transport by their correlations.
pub fn (db Database) transport(names []string) !transport.Mixture {
	mut records := []Record{}
	for name in names {
		records << db.find(name)!
	}
	return mixture(records)
}

// transport is the fluid's transport model (D10): IAPWS's own for water by
// IAPWS-95, their industrial forms for water by IF97, the compounds'
// correlations for any other equation.
pub fn (s Setup) transport() !transport.Model {
	if s.equation == 'IAPWS-95' {
		return water.Transport95{}
	}
	if s.equation == 'IF97' {
		return water.TransportIf97{}
	}
	return mixture(s.records)!
}

// wilke is a gas mixture's viscosity at low pressure by Wilke's rule, for a
// PFR's Ergun drop (M9's D10, moved by M5b's D2): each compound's
// vapour-viscosity correlation and molar mass.
pub fn wilke(records []Record) !transport.Wilke {
	mut pure := []transport.Correlated{}
	mut mw := []f64{}
	for r in records {
		p := r.property('vapor_viscosity')!
		pure << transport.Correlated{
			key:   p.key
			c:     p.c
			scale: p.scale
		}
		m := r.constants.mw or { return error('${r.name} has no molar mass') }
		if !(m > 0) {
			return error("${r.name}'s molar mass, ${m} kg/mol, isn't one")
		}
		mw << m
	}
	return transport.new_wilke(pure, mw)!
}
