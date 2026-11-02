module flowsheet

import components
import math
import numeric
import os
import quantity
import thermo
import toml
import toml.ast

// The format this reflux reads (D9). A later one comes with a migration.
const file_format = 1

// How deep arrays and inline tables may nest. No flowsheet goes past 2, and
// V's parser recurses, so it's what stops a file overflowing the stack.
const max_nesting = 32

// LoadParams are load's options.
@[params]
pub struct LoadParams {
pub:
	// db is the component data, if the caller has it: otherwise it's read, with
	// the file's own data directories after the built-in ones.
	db ?components.Database
	// dir is what a file's data directories are relative to.
	dir string
	// makers are unit types beyond the built-in ones (D1).
	makers map[string]Maker
}

// Loaded is a flowsheet file, read and checked: ready to solve, or to write
// back.
pub struct Loaded {
pub:
	fluid  thermo.Fluid
	sheet  Flowsheet
	params SolveParams
	root   map[string]ast.Value // the file, for write
}

// solve solves the flowsheet with the file's own solver settings.
pub fn (l Loaded) solve() !Solution {
	return l.sheet.solve(l.fluid, l.params)
}

// load reads a flowsheet file. See load_text.
pub fn load(path string, params LoadParams) !Loaded {
	text := os.read_file(path) or { return file_error(0, path, "it can't be read: ${err.msg()}", []) }
	return load_text(text, LoadParams{
		...params
		dir: os.dir(path)
	})
}

// load_text reads a flowsheet (D9): docs/flowsheet-format.md describes the
// format. Everything is checked here, each value for its own reason, and the
// streams for joining up, so that a file that loads can only fail in solving.
// Failure is a FlowsheetError of kind .bad_file with the line, or the failure
// of whatever it was building.
pub fn load_text(text string, params LoadParams) !Loaded {
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
	version := top.integer('format')!
	if version != file_format {
		return top.fail('format', 'this reflux reads format ${file_format}, not ${version}', [])
	}
	mut names := top.child('components')!
	component_names := names.strings('names')!
	if component_names.len == 0 {
		return names.fail('names', 'a flowsheet needs at least one component', [])
	}
	for i, name in component_names {
		if name in component_names[..i] {
			return names.fail('names', '${name} is listed twice', [])
		}
	}
	data := if names.has('data') { names.strings('data')! } else { []string{} }
	names.finish()!
	db := database(params, data) or { return top.fail('components', err.msg(), []) }

	mut model := top.child('model')!
	fl := read_model(mut model, db, component_names, names.line('names'))!

	mut solver := if top.has('solver') {
		top.child('solver')!
	} else {
		Section{ path: 'solver', text: text }
	}
	settings := read_solver(mut solver)!

	mut sheet := Flowsheet{}
	mut feeds := top.child('feed')!
	for name in feeds.table.keys() {
		mut s := feeds.child(name)!
		s.components = component_names
		stream := read_stream(mut s, fl) or { return s.located(err) }
		sheet.add_feed(name, stream) or { return s.located(err) }
	}

	mut estimates := map[string]Stream{}
	if top.has('estimate') {
		mut es := top.child('estimate')!
		for name in es.table.keys() {
			mut s := es.child(name)!
			s.components = component_names
			estimates[name] = read_stream(mut s, fl) or { return s.located(err) }
		}
	}

	mut makers := built_in_makers()
	for kind, maker in params.makers {
		makers[kind] = maker
	}
	mut units := top.child('unit')!
	for name in units.table.keys() {
		mut s := units.child(name)!
		s.components = component_names
		read_unit(mut s, name, makers, mut sheet)!
	}
	top.finish()!

	solve_params := SolveParams{
		...settings
		estimates: estimates
	}
	mut torn := solve_params.tears.clone()
	torn << estimates.keys()
	sheet.plan(tears: torn) or {
		if err is FlowsheetError && err.unit != '' {
			mut s := units.child(err.unit)!
			return s.located(err)
		}
		return top.located(err)
	}
	return Loaded{
		fluid:  fl
		sheet:  sheet
		params: solve_params
		root:   root
	}
}

// check_text refuses what would bring V's parser down, which nothing can catch:
// a \u escape, whose check runs off the end of its string when it's the right
// distance from it (a bug in V 0.5.2's toml.checker), and nesting deep enough
// to overflow the stack. It walks strings and comments as TOML has them, so
// that a Windows path in a 'literal string' is left alone.
fn check_text(text string) ! {
	mut line := 1
	mut depth := 0
	mut i := 0
	for i < text.len {
		c := text[i]
		match c {
			`\n` {
				line++
				i++
			}
			`#` {
				for i < text.len && text[i] != `\n` {
					i++
				}
			}
			`'` {
				i, line = skip_literal(text, i, line)
			}
			`"` {
				i, line = skip_basic(text, i, line)!
			}
			`[`, `{` {
				depth++
				if depth > max_nesting {
					return file_error(line, '', 'it nests more than ${max_nesting} deep', [])
				}
				i++
			}
			`]`, `}` {
				depth--
				i++
			}
			else {
				i++
			}
		}
	}
}

// skip_literal is past a literal string, single or triple quoted.
fn skip_literal(text string, start int, first int) (int, int) {
	mut line := first
	triple := text[start..].starts_with("'''")
	mut i := start + if triple { 3 } else { 1 }
	for i < text.len {
		if text[i] == `\n` {
			line++
			if !triple {
				return i, line - 1
			}
		}
		if triple && text[i..].starts_with("'''") {
			return i + 3, line
		}
		if !triple && text[i] == `'` {
			return i + 1, line
		}
		i++
	}
	return i, line
}

// skip_basic is past a basic string, single or triple quoted.
fn skip_basic(text string, start int, first int) !(int, int) {
	mut line := first
	triple := text[start..].starts_with('"""')
	mut i := start + if triple { 3 } else { 1 }
	for i < text.len {
		c := text[i]
		if c == `\\` {
			if i + 1 < text.len && text[i + 1] in [`u`, `U`] {
				return file_error(line, '', r'a \u escape: write the character itself', [])
			}
			i += 2
			continue
		}
		if c == `\n` {
			line++
			if !triple {
				return i, line - 1
			}
		}
		if triple && text[i..].starts_with('"""') {
			return i + 3, line
		}
		if !triple && c == `"` {
			return i + 1, line
		}
		i++
	}
	return i, line
}

// syntax_error is the line and reason in an error from V's parser, which looks
// like `parser.Parser.value did not expect "x" ... (3,14) in ...text...`.
fn syntax_error(message string) (int, string) {
	reason := message.all_after(' ')
	open := reason.last_index('(') or { return 0, reason }
	close := reason.index_after(')', open) or { return 0, reason }
	line := reason[open + 1..close].all_before(',').int()
	if line > 0 {
		return line, reason[..open].trim_space()
	}
	return 0, reason
}

// database is the component data for a file.
fn database(params LoadParams, data []string) !components.Database {
	if db := params.db {
		if data.len == 0 {
			return db
		}
	}
	mut dirs := [os.getenv_opt('REFLUX_DATA') or { components.default_dir }]
	for d in data {
		dirs << if os.is_abs_path(d) { d } else { os.join_path(params.dir, d) }
	}
	return components.load(dirs)!
}

// read_model builds the fluid from [model].
fn read_model(mut s Section, db components.Database, names []string, names_line int) !thermo.Fluid {
	equation := s.text('equation')!
	kinds := {
		'PR':  thermo.CubicKind.pr
		'SRK': .srk
		'RK':  .rk
		'vdW': .vdw
	}
	kind := kinds[equation] or {
		return s.fail('equation', 'no equation ${equation}', quantity.closest(equation, kinds.keys()))
	}
	mut opts := components.FluidOptions{
		kij:   if s.has('kij') { s.boolean('kij')! } else { true }
		alpha: read_alpha(mut s)!
	}
	mut pairs := []components.Kij{}
	if s.has('kij_overrides') {
		for mut t in s.tables('kij_overrides')! {
			a := t.text('a')!
			b := t.text('b')!
			k := t.number('k')!
			t.finish()!
			for name in [a, b] {
				if name !in names {
					return t.fail('a', '${name} is not one of the components', quantity.closest(name,
						names))
				}
			}
			pairs << components.Kij{
				cas1: (db.find(a) or { return t.fail('a', err.msg(), []) }).ids.cas
				cas2: (db.find(b) or { return t.fail('b', err.msg(), []) }).ids.cas
				k:    k
			}
		}
	}
	s.finish()!
	fl := db.fluid(kind, names,
		kij:       opts.kij
		alpha:     opts.alpha
		overrides: pairs
	) or {
		mut hint := []string{}
		if err is components.LookupError {
			hint = err.suggestions.clone()
			return file_error(names_line, 'components.names', 'no compound ${err.query}', hint)
		}
		return s.fail('equation', err.msg(), [])
	}
	return fl
}

fn read_alpha(mut s Section) !components.AlphaChoice {
	if !s.has('alpha') {
		return .equation
	}
	alpha := s.text('alpha')!
	return match alpha {
		'equation' {
			components.AlphaChoice.equation
		}
		'mathias_copeman' {
			.mathias_copeman
		}
		else {
			s.fail('alpha', 'no alpha ${alpha}', quantity.closest(alpha, ['equation', 'mathias_copeman']))
		}
	}
}

// read_solver reads [solver].
fn read_solver(mut s Section) !SolveParams {
	mut method := Method.wegstein
	if s.has('method') {
		name := s.text('method')!
		method = match name {
			'direct' {
				Method.direct
			}
			'wegstein' {
				.wegstein
			}
			'broyden' {
				.broyden
			}
			else {
				return s.fail('method', 'no method ${name}', quantity.closest(name, [
					'direct',
					'wegstein',
					'broyden',
				]))
			}
		}
	}
	defaults := SolveParams{}
	tol := if s.has('tol') { s.number('tol')! } else { defaults.tol }
	if !(tol > 0.0) || !math.is_finite(tol) {
		return s.fail('tol', 'a tolerance of ${tol}', [])
	}
	passes := if s.has('max_iterations') { s.integer('max_iterations')! } else { 100 }
	if passes < 1 || passes > 1000000 {
		return s.fail('max_iterations', '${passes} iterations', [])
	}
	tears := if s.has('tears') { s.strings('tears')! } else { []string{} }
	s.finish()!
	return SolveParams{
		method:         method
		tol:            tol
		max_iterations: int(passes)
		tears:          tears
	}
}

// read_stream reads a feed or an estimate: its state (D2), and its flows.
fn read_stream(mut s Section, fl thermo.Fluid) !Stream {
	flows := read_flows(mut s, fl)!
	has_t, has_p := s.has('T'), s.has('P')
	has_h, has_b := s.has('H'), s.has('vapour_fraction')
	given := [has_t, has_p, has_h, has_b].filter(it).len
	if given != 2 || (has_t && has_h) || (!has_p && !has_t) {
		return s.fail('', 'a feed needs T and P, P and H, T and vapour_fraction, or P and vapour_fraction',
			[])
	}
	mut out := Stream{}
	if has_t && has_p {
		t, p := s.quantity('T', .temperature)!, s.quantity('P', .pressure)!
		s.finish()!
		out = feed_tp(fl, flows, t, p)!
	} else if has_p && has_h {
		p, h := s.quantity('P', .pressure)!, s.quantity('H', .molar_enthalpy)!
		s.finish()!
		out = feed_ph(fl, flows, p, h)!
	} else if has_t {
		t, b := s.quantity('T', .temperature)!, s.fraction('vapour_fraction')!
		s.finish()!
		out = feed_tvf(fl, flows, t, b)!
	} else {
		p, b := s.quantity('P', .pressure)!, s.fraction('vapour_fraction')!
		s.finish()!
		out = feed_pvf(fl, flows, p, b)!
	}
	return out
}

// read_flows is a stream's molar flows, by whichever basis it's given (D9).
fn read_flows(mut s Section, fl thermo.Fluid) ![]f64 {
	comps := fl.model.components()
	bases := ['mole_fractions', 'mass_fractions', 'molar_flows', 'mass_flows']
	basis := s.one_of(bases)!
	mut flows := []f64{len: comps.len}
	match basis {
		'molar_flows' {
			for i, x in s.by_component(basis, .molar_flow)! {
				flows[i] = x
			}
		}
		'mass_flows' {
			for i, x in s.by_component(basis, .mass_flow)! {
				flows[i] = x / comps[i].mw
			}
		}
		else {
			x := s.by_component(basis, .fraction)!
			sum := numeric.sum(x)
			if math.abs(sum - 1.0) > 1e-6 {
				return s.fail(basis, 'they sum to ${sum}, not 1', [])
			}
			// Moles of each component per unit of the mixture.
			mut per := []f64{len: comps.len}
			mut mass := 0.0 // kg/mol of mixture
			mut moles := 0.0 // mol/kg of mixture
			for i, v in x {
				per[i] = v / sum
				mass += per[i] * comps[i].mw
				moles += per[i] / comps[i].mw
			}
			total := s.quantity_either('flow', .molar_flow, .mass_flow)!
			for i in 0 .. x.len {
				flows[i] = match true {
					basis == 'mole_fractions' && total.dim == .molar_flow {
						per[i] * total.value
					}
					basis == 'mole_fractions' {
						per[i] * total.value / mass
					}
					total.dim == .mass_flow {
						per[i] / comps[i].mw * total.value
					}
					else {
						per[i] / comps[i].mw / moles * total.value
					}
				}
			}
		}
	}
	return flows
}

// read_unit builds one unit and adds it to the flowsheet.
fn read_unit(mut s Section, name string, makers map[string]Maker, mut sheet Flowsheet) ! {
	kind := s.text('type')!
	maker := makers[kind] or {
		return s.fail('type', 'no unit type ${kind}', quantity.closest(kind, makers.keys()))
	}
	s.name = name
	u := maker(mut s) or { return s.located(err) }
	ports := s.ports(u)!
	s.finish()!
	sheet.add_unit(u, ports) or { return s.located(err) }
}
