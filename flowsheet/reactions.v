module flowsheet

import chem
import components
import expr
import kinetics
import os
import quantity
import thermo
import toml
import toml.ast

// A file's reactions (M9, D11): [definitions], [site.NAME] and
// [reaction.NAME], in a format-2 flowsheet file or a file of reactions
// alone, read into kinetics.Kinetics with every value checked and every
// failure at its line. docs/flowsheet-format.md describes the sections.

// The format whose files may hold reactions.
const reactions_format = 2

// The flowsheet's own sections, which a file's reactions are read beside
// and left to the flowsheet's reader.
const flowsheet_sections = ['model', 'solver', 'feed', 'estimate', 'unit']

// KineticsFile is a file's reactions, read and built.
pub struct KineticsFile {
pub:
	kinetics kinetics.Kinetics
	names    []string // the components, as the file names them
}

// load_kinetics reads a file's reactions. See load_kinetics_text.
pub fn load_kinetics(path string, params LoadParams) !KineticsFile {
	text := os.read_file(path) or { return file_error(0, path, "it can't be read: ${err.msg()}", []) }
	return load_kinetics_text(text, LoadParams{
		...params
		dir: os.dir(path)
	})
}

// load_kinetics_text reads a file's components, definitions, sites and
// reactions, and builds them: every equation balanced, every name bound and
// every rate's units checked. The flowsheet's own sections, if it has them,
// are left to load.
pub fn load_kinetics_text(text string, params LoadParams) !KineticsFile {
	check_text(text)!
	doc := toml.parse_text(text) or {
		line, reason := syntax_error(err.msg())
		return file_error(line, '', reason, [])
	}
	root := doc.ast.table as map[string]ast.Value
	mut top := Section{
		path:  ''
		table: root
		text:  text
	}
	given := top.integer('format')!
	if given != reactions_format {
		return top.fail('format', 'reactions need format ${reactions_format}, not ${given}', [])
	}
	mut names := top.child('components')!
	component_names := names.strings('names')!
	if component_names.len == 0 {
		return names.fail('names', 'a file needs at least one component', [])
	}
	data := if names.has('data') { names.strings('data')! } else { []string{} }
	names.finish()!
	db := database(params, data) or { return top.fail('components', err.msg(), []) }
	mut species := []kinetics.Species{}
	for i, name in component_names {
		if name in component_names[..i] {
			return names.fail('names', '${name} is listed twice', [])
		}
		r := db.find(name) or {
			if err is components.LookupError {
				return file_error(names.line('names'), 'components.names', err.msg(), err.suggestions)
			}
			return names.fail('names', err.msg(), [])
		}
		species << species_of(r)
	}
	for s in flowsheet_sections {
		if top.has(s) {
			top.used << s
		}
	}
	k := read_reactions(mut top, species)!
	top.finish()!
	return KineticsFile{
		kinetics: k
		names:    component_names
	}
}

// read_reactions reads a file's [definitions], [site.NAME] and
// [reaction.NAME] over its species, and builds them: every equation
// balanced, every name bound and every rate's units checked, each failure
// at its line. The caller finishes the file's top table.
fn read_reactions(mut top Section, species []kinetics.Species) !kinetics.Kinetics {
	mut own := []expr.Source{}
	if top.has('definitions') {
		mut d := top.child('definitions')!
		for key in d.table.keys() {
			own << expr.Source{key, d.text(key)!}
		}
	}
	mut sites := []kinetics.Site{}
	mut site_sections := map[string]Section{}
	if top.has('site') {
		mut all := top.child('site')!
		for name in all.table.keys() {
			mut s := all.child(name)!
			sites << read_site(mut s, name)!
			s.finish()!
			site_sections[name] = s
		}
	}
	mut reactions := []kinetics.Reaction{}
	mut laws := []kinetics.Rated{}
	mut reaction_sections := []Section{}
	mut all := top.child('reaction')!
	for name in all.table.keys() {
		mut s := all.child(name)!
		equation := s.text('equation')!
		r := kinetics.parse_equation(name, equation, species) or {
			return s.fail('equation', err.msg().all_after('${name}: '), [])
		}
		if s.has('law') {
			laws << kinetics.Rated{reactions.len, read_law(mut s)!}
		} else {
			no_law(s)!
		}
		s.finish()!
		reactions << r
		reaction_sections << s
	}
	if reactions.len == 0 {
		return all.fail('', 'a file needs at least one reaction', [])
	}
	return kinetics.build(species, reactions, sites, laws, own) or {
		if err is kinetics.KineticsError {
			if err.reaction >= 0 {
				return reaction_sections[err.reaction].fail('', err.msg(), [])
			}
			if err.site != '' {
				if s := site_sections[err.site] {
					return s.fail('', err.msg(), [])
				}
			}
			if err.definition != '' && top.has('definitions') {
				mut d := top.child('definitions')!
				if d.has(err.definition) {
					return d.fail(err.definition, err.msg(), [])
				}
			}
		}
		return top.fail('', err.msg(), [])
	}
}

// species_of is a compound as reactions know it, from its component file.
// A mixture's pseudo-compound has no formula, and a compound without
// ideal-gas data no K from the Gibbs energies: only a reaction that needs
// them refuses.
pub fn species_of(r components.Record) kinetics.Species {
	// Each Option is made apart: V 0.5.2 accepts `r.formula() or { none }`
	// for an Option field and writes C that no compiler takes.
	mut formula := ?chem.Formula(none)
	if f := r.formula() {
		formula = f
	}
	mut gas := ?thermo.IdealGas(none)
	if g := r.ideal_gas() {
		gas = g
	}
	return kinetics.Species{
		name:    r.name
		cas:     r.ids.cas
		aliases: r.ids.aliases
		formula: formula
		charge:  r.constants.charge or { 0 }
		gas:     gas
		dgf:     r.constants.dgf
	}
}

// read_site reads [site.NAME]: its symbol, if it names its term, and its
// terms.
fn read_site(mut s Section, name string) !kinetics.Site {
	symbol := if s.has('symbol') { s.text('symbol')! } else { '' }
	terms := read_terms(mut s, 'terms')!
	return kinetics.Site{
		name:   name
		symbol: symbol
		terms:  terms
	}
}

// read_terms reads a list of terms, each { name, K, a }.
fn read_terms(mut s Section, key string) ![]kinetics.Monomial {
	mut out := []kinetics.Monomial{}
	mut tables := s.tables(key)!
	for mut t in tables {
		name := if t.has('name') { t.text('name')! } else { '' }
		k := read_constant(mut t, 'K', 'dH', 'K_ref')!
		powers := read_powers(mut t, 'a')!
		if powers.len == 0 {
			return t.fail('a', 'a term needs at least one species', [])
		}
		t.finish()!
		out << kinetics.Monomial{name, k, powers}
	}
	return out
}

// read_constant reads a constant: an expression, as a string, or a table:
// { A, Ea } or { k_ref, Ea, T_ref } for a rate constant, { A, dH } or
// { K_ref, dH, T_ref } for an adsorption constant, each a number with its
// unit.
fn read_constant(mut s Section, key string, energy string, centred string) !kinetics.TempFn {
	v := s.require(key)!
	if v is ast.Quoted {
		return kinetics.Written{v.text}
	}
	mut t := s.child(key)!
	reference := if t.has(centred) { centred } else { 'A' }
	a := t.text(reference)!
	e := t.text(energy)!
	t_ref := if reference == centred { t.text('T_ref')! } else { '' }
	t.finish()!
	if energy == 'Ea' {
		return kinetics.Arrhenius{a, e, t_ref}
	}
	return kinetics.VantHoff{a, e, t_ref}
}

// read_powers reads a table of species and their exponents.
fn read_powers(mut s Section, key string) ![]kinetics.Power {
	v := s.require(key)!
	if v is map[string]ast.Value {
		mut out := []kinetics.Power{}
		for species, x in v {
			if x is ast.Number {
				out << kinetics.Power{species, x.f64()}
			} else {
				return file_error(line_of(x), s.full(key), "${species}'s exponent should be a number",
					[])
			}
		}
		return out
	}
	return s.fail(key, 'it should be a table of species and their exponents', [])
}

// choice reads one of a key's words, or its default.
fn choice(mut s Section, key string, words []string, default string) !string {
	if !s.has(key) {
		return default
	}
	w := s.text(key)!
	if w !in words {
		return s.fail(key, 'no ${key} ${w}: one of ${words.join(', ')}', quantity.closest(w,
			words))
	}
	return w
}

// no_law refuses a law's keys on a reaction that has no law.
fn no_law(s Section) ! {
	for key in ['k', 'k_name', 'orders', 'K', 'K_name', 'basis', 'per', 'driving', 'site', 'site_power',
		'terms', 'symbol', 'adsorbed', 'gas', 'rate'] {
		if s.has(key) {
			return s.fail(key, 'a reaction without a law has no ${key}: give its law', [])
		}
	}
}

// read_law reads a reaction's rate law.
fn read_law(mut s Section) !kinetics.RateLaw {
	law := choice(mut s, 'law', ['power_law', 'lhhw', 'eley_rideal', 'expression'], '')!
	basis := match choice(mut s, 'basis', ['partial_pressure', 'fugacity', 'concentration',
		'mole_fraction'], 'partial_pressure')! {
		'fugacity' { kinetics.Basis.fugacity }
		'concentration' { kinetics.Basis.concentration }
		'mole_fraction' { kinetics.Basis.fraction }
		else { kinetics.Basis.pressure }
	}
	per := if choice(mut s, 'per', ['catalyst_mass', 'volume'], 'catalyst_mass')! == 'volume' {
		kinetics.Per.volume
	} else {
		kinetics.Per.mass
	}
	if law == 'expression' {
		return kinetics.Expression{
			text:  s.text('rate')!
			basis: basis
			per:   per
		}
	}
	k := read_constant(mut s, 'k', 'Ea', 'k_ref')!
	k_name := if s.has('k_name') { s.text('k_name')! } else { '' }
	k_text := s.text('K')!
	eq := match k_text {
		'gibbs' {
			kinetics.Equilibrium{
				kind: .thermo
			}
		}
		'none' {
			kinetics.Equilibrium{}
		}
		else {
			kinetics.Equilibrium{
				kind:        .correlation
				correlation: kinetics.Written{k_text}
				name:        if s.has('K_name') { s.text('K_name')! } else { '' }
			}
		}
	}
	if law == 'eley_rideal' {
		return kinetics.EleyRideal{
			k:           k
			k_name:      k_name
			adsorbed:    s.text('adsorbed')!
			gas:         s.text('gas')!
			equilibrium: eq
			basis:       basis
			per:         per
			site:        s.text('site')!
		}
	}
	orders := read_powers(mut s, 'orders')!
	form := if choice(mut s, 'driving', ['approach', 'difference'], 'approach')! == 'difference' {
		kinetics.Form.difference
	} else {
		kinetics.Form.approach
	}
	if law == 'power_law' {
		return kinetics.PowerLaw{
			k:           k
			k_name:      k_name
			orders:      orders
			equilibrium: eq
			form:        form
			basis:       basis
			per:         per
		}
	}
	terms := if s.has('terms') { ?[]kinetics.Monomial(read_terms(mut s, 'terms')!) } else { none }
	return kinetics.Lhhw{
		k:           k
		k_name:      k_name
		orders:      orders
		equilibrium: eq
		form:        form
		basis:       basis
		per:         per
		site:        s.text('site')!
		power:       if s.has('site_power') { s.number('site_power')! } else { 1.0 }
		terms:       terms
		symbol:      if s.has('symbol') { s.text('symbol')! } else { '' }
	}
}
