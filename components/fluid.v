module components

import thermo

// AlphaChoice picks each component's α function.
pub enum AlphaChoice {
	equation        // the equation's own: Soave's for SRK and PR
	mathias_copeman // ChemSep's Mathias–Copeman constants for PR (mc-PR.lib)
}

@[params]
pub struct FluidOptions {
pub:
	alpha AlphaChoice
	// kij takes binary interaction parameters from ChemSep's bank for the
	// equation (PR and SRK; there's none for vdW and RK). Off, every pair is
	// 0 except the overrides.
	kij       bool = true
	overrides []Kij // pairs set by hand
}

// Setup is a fluid built from the database, with a record of how: the
// equation, the records used, the kᵢⱼ matrix and where each value came from,
// and how each α was chosen.
pub struct Setup {
pub:
	kind    thermo.CubicKind
	fluid   thermo.Fluid
	records []Record
	kij     KijMatrix
	alphas  []string // per component, e.g. 'Mathias–Copeman, mc-PR.lib line 99'
}

// fluid builds a thermo.Fluid for a cubic equation of state from compounds
// named by CAS number, name, alias or InChIKey:
//
//     fl := db.fluid(.pr, ['methane', 'n-butane'])!
//
// See setup for what it does, and for a record of how it was done.
pub fn (db Database) fluid(kind thermo.CubicKind, names []string, opts FluidOptions) !thermo.Fluid {
	return db.setup(kind, names, opts)!.fluid
}

// setup builds the fluid:
// - each compound's Tc, Pc, ω and M, and its ideal-gas Cp, ΔH°f and S°;
// - kᵢⱼ from ChemSep's bank for the equation, in the database's root, with
//   any overrides (the first listing of a repeated pair, M2 decision D12);
// - α as chosen. Mathias–Copeman is for PR only, since that's what mc-PR.lib
//   was fitted for. A compound without a complete set gets the equation's
//   own α, and `alphas` says so.
pub fn (db Database) setup(kind thermo.CubicKind, names []string, opts FluidOptions) !Setup {
	if names.len == 0 {
		return error('a fluid needs at least one compound')
	}
	if opts.alpha == .mathias_copeman && kind != .pr {
		return error("ChemSep's Mathias–Copeman constants (mc-PR.lib) are fitted for PR, not ${kind}")
	}
	mut records := []Record{}
	mut comps := []thermo.Component{}
	mut ideal := []thermo.IdealGas{}
	for name in names {
		r := db.find(name)!
		comps << r.component()!
		ideal << r.ideal_gas()!
		records << r
	}
	cas := records.map(it.ids.cas)
	bank := if opts.kij && kind in [.pr, .srk] {
		load_kij(db.vendor, kind)!
	} else {
		KijBank{
			file: 'no kij file'
		}
	}
	kij := bank.matrix(cas, overrides: opts.overrides)!
	mut alphas := []thermo.Alpha{}
	mut how := []string{}
	if opts.alpha == .mathias_copeman {
		lib := load_mathias_copeman_pr(db.vendor)!
		for i, r in records {
			own := kind.default_alpha(comps[i].omega)
			e := lib.find(r.ids.cas) or {
				alphas << own
				how << "the equation's own: mc-PR.lib has no ${r.name}"
				continue
			}
			a := e.alpha() or {
				alphas << own
				how << "the equation's own: ${err}"
				continue
			}
			alphas << a
			how << 'Mathias–Copeman, mc-PR.lib line ${e.line}'
		}
	} else {
		for c in comps {
			alphas << kind.default_alpha(c.omega)
			how << "the equation's own"
		}
	}
	eos := thermo.new_cubic(kind, comps, kij: kij.k, alphas: alphas)!
	return Setup{
		kind:    kind
		fluid:   thermo.new_fluid(eos, ideal)!
		records: records
		kij:     kij
		alphas:  how
	}
}

// to_toml writes how the fluid was built, for a caller to read back (the C
// API's reflux_fluid_setup): the equation and the kᵢⱼ bank; each component's
// CAS number, name, Tc (K), Pc (Pa), ω, molar mass (kg/mol), ΔH°f (J/mol)
// and S° (J/(mol·K)) where it has them, and α; and each pair i < j's kᵢⱼ with
// where it came from. Floats are V's shortest
// round-trip form, so they read back exactly.
pub fn (s Setup) to_toml() !string {
	mut out := ['equation = ${quote(equation_name(s.kind))!}', 'kij_bank = ${quote(s.kij.bank)!}']
	comps := s.fluid.model.components()
	for i, r in s.records {
		c := comps[i]
		out << ''
		out << '[[components]]'
		out << 'cas = ${quote(r.ids.cas)!}'
		out << 'name = ${quote(r.name)!}'
		out << 'tc = ${float(c.tc)!}'
		out << 'pc = ${float(c.pc)!}'
		out << 'omega = ${float(c.omega)!}'
		out << 'mw = ${float(c.mw)!}'
		ig := s.fluid.ideal[i]
		if dhf := ig.dhf {
			out << 'dhf = ${float(dhf)!}'
		}
		if s0 := ig.s0 {
			out << 's0 = ${float(s0)!}'
		}
		out << 'alpha = ${quote(s.alphas[i])!}'
	}
	for p in s.kij.pairs {
		out << ''
		out << '[[pairs]]'
		out << 'i = ${p.i}'
		out << 'j = ${p.j}'
		out << 'k = ${float(p.k)!}'
		out << 'source = ${quote(p.source.str())!}'
		if used := p.used {
			out << 'line = ${used.line}'
			out << 'comment = ${quote(used.comment)!}'
		}
		if p.others.len > 0 {
			out << 'others = [${p.others.map(it.line.str()).join(', ')}]'
		}
	}
	return out.join('
') + '
'
}

// equation_name is an equation's name as a flowsheet file writes it.
pub fn equation_name(kind thermo.CubicKind) string {
	return match kind {
		.vdw { 'vdW' }
		.rk { 'RK' }
		.srk { 'SRK' }
		.pr { 'PR' }
	}
}
