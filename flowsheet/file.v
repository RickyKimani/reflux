module flowsheet

import activity
import components
import kinetics
import math
import numeric
import os
import quantity
import thermo
import toml
import toml.ast

// The formats this reflux reads (D9): 1, M4's; 2, which adds M9's
// reactions and reactors (D11); 3, which adds M8's specs (D16); and 4, which
// adds M6's γ–φ fluids and the decanter (D12). A later one comes with a
// migration.
const file_formats = [1, 2, 3, 4]

// How deep arrays and inline tables may nest. No flowsheet goes past 2, and
// V's parser recurses, so it's what stops a file overflowing the stack.
const max_nesting = 32

// LoadParams are load's options.
@[params]
pub struct LoadParams {
pub:
	// db is the component data, if the caller has it: otherwise the data root's
	// (components.data_root()) is read. A file's own data directories go after
	// its root's compounds.
	db ?components.Database
	// dir is what a file's data directories are relative to.
	dir string
	// makers are unit types beyond the built-in ones (D1).
	makers map[string]Maker
}

// Context is what a unit's maker may need beyond its own table: the file's
// format, its components as reactions know them, with their records, and
// its reactions, built (format 2's, D11).
pub struct Context {
pub:
	format   int
	species  []kinetics.Species
	records  []components.Record
	kinetics ?kinetics.Kinetics
}

// Loaded is a flowsheet file, read and checked: ready to solve, or to write
// back.
pub struct Loaded {
pub:
	fluid  thermo.Fluid
	model  ModelReport // how the fluid was built, for the report (M6's D13)
	sheet  Flowsheet
	params SolveParams
	root   map[string]ast.Value // the file, for write
	// For the report: the components as the file names them, each unit's type in
	// the order of sheet.units, and the plan solve will follow.
	names []string
	types []string
	plan  Plan
	// What remaking a feed or unit needs (M8's D4): the file's text, for its
	// lines, the makers and what they're given.
	text    string
	makers  map[string]Maker
	context Context
	// What an expression over the flowsheet may name (M8's D3), and the
	// file's specs (D5).
	variables Variables
	specs     []Spec
	studies   []StudyPlan // the file's [study.NAME], for the study module (D16)
	// The file's [optimize.NAME], for the study module (D13).
	optimizations []OptimizePlan
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
	given := top.integer('format')!
	if given !in file_formats.map(i64(it)) {
		return top.fail('format', 'this reflux reads format 1, 2, 3 or 4, not ${given}', [])
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
	modelled := read_model(mut model, db, component_names, names.line('names'), int(given))!
	fl := modelled.fluid

	// The file's components as reactions know them, for a reactor's maker,
	// and in format 2 its reactions.
	mut records := []components.Record{}
	mut species := []kinetics.Species{}
	for name in component_names {
		r := db.find(name) or { return names.fail('names', err.msg(), []) }
		records << r
		species << species_of(r)
	}
	mut built := ?kinetics.Kinetics(none)
	if given >= 2 && (top.has('reaction') || top.has('site') || top.has('definitions')) {
		built = read_reactions(mut top, species)!
	}
	context := Context{
		format:   int(given)
		species:  species
		records:  records
		kinetics: built
	}

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
	mut types := []string{}
	for name in units.table.keys() {
		mut s := units.child(name)!
		s.components = component_names
		s.context = context
		types << read_unit(mut s, name, makers, mut sheet)!
	}
	solve_params := SolveParams{
		...settings
		estimates: estimates
	}
	members, keys := member_names(component_names, species)
	mut loaded := Loaded{
		fluid:     fl
		model:     modelled.report
		sheet:     sheet
		params:    solve_params
		root:      root
		names:     component_names
		types:     types
		text:      text
		makers:    makers
		context:   context
		variables: sheet.variables(members, keys)
	}
	if top.has('spec') {
		if given < 3 {
			return top.fail('spec', 'specs need format 3', [])
		}
		mut specs := top.child('spec')!
		loaded = Loaded{
			...loaded
			specs: read_specs(mut specs, loaded)!
		}
	}
	if top.has('study') {
		if given < 3 {
			return top.fail('study', 'studies need format 3', [])
		}
		mut studies := top.child('study')!
		loaded = Loaded{
			...loaded
			studies: read_studies(mut studies, loaded)!
		}
	}
	if top.has('optimize') {
		if given < 3 {
			return top.fail('optimize', 'optimisations need format 3', [])
		}
		mut tables := top.child('optimize')!
		loaded = Loaded{
			...loaded
			optimizations: read_optimizations(mut tables, loaded)!
		}
	}
	top.finish()!

	mut torn := solve_params.tears.clone()
	torn << estimates.keys()
	plan := sheet.plan(tears: torn, links: loaded.links()!) or {
		if err is FlowsheetError && err.unit != '' {
			mut s := units.child(err.unit)!
			return s.located(err)
		}
		if err is FlowsheetError && err.key.starts_with('spec.') {
			return err
		}
		return top.located(err)
	}
	return Loaded{
		...loaded
		plan: plan
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

// database is the component data for a file: the caller's, or the data root's
// (components.data_root()), with the file's own data directories after its
// compounds.
fn database(params LoadParams, data []string) !components.Database {
	mut root := components.data_root()
	if db := params.db {
		if data.len == 0 {
			return db
		}
		root = db.root
	}
	mut extra := []string{}
	for d in data {
		extra << if os.is_abs_path(d) { d } else { os.join_path(params.dir, d) }
	}
	return components.open(root, extra)!
}

// Modelled is a file's fluid and how it was built: its report (M6's D13),
// and a γ–φ fluid's setup, whose pairs a fit may vary (D16).
struct Modelled {
	fluid     thermo.Fluid
	report    ModelReport
	gamma_phi ?components.GammaPhiSetup
}

// The equations of state a file names.
const equations = {
	'PR':  thermo.CubicKind.pr
	'SRK': .srk
	'RK':  .rk
	'vdW': .vdw
}

// The activity models a file names (M6's D12).
const liquid_models = {
	'Wilson':          components.LiquidModel.wilson
	'NRTL':            .nrtl
	'UNIQUAC':         .uniquac
	'UNIFAC':          .unifac
	'UNIFAC-LLE':      .unifac_lle
	'UNIFAC-Dortmund': .unifac_dortmund
}

// read_model builds the fluid from [model]: an equation of state, or in
// format 4 a γ–φ fluid, whose liquid is an activity model (D12).
fn read_model(mut s Section, db components.Database, names []string, names_line int, format int) !Modelled {
	if s.has('liquid') {
		if format < 4 {
			return s.fail('liquid', 'an activity model needs format 4 (M6), and this file is format ${format}',
				[])
		}
		return read_gamma_phi(mut s, db, names, names_line)
	}
	if format >= 4 && !s.has('equation') {
		return s.fail('', 'a model needs equation, for an equation of state, or liquid, for an activity model',
			[])
	}
	for key in ['vapour', 'poynting', 'henry', 'missing', 'pairs'] {
		if s.has(key) && format >= 4 {
			return s.fail(key, "${key} is a γ–φ fluid's, with an activity model as its liquid", [])
		}
	}
	equation := s.text('equation')!
	if equation == 'IAPWS-95' {
		return read_iapws95(mut s, db, names, names_line, format)
	}
	kind := equations[equation] or {
		mut known := equations.keys()
		known << 'IAPWS-95'
		return s.fail('equation', 'no equation ${equation}', quantity.closest(equation, known))
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
	setup := db.setup(kind, names,
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
	return Modelled{
		fluid:  setup.fluid
		report: cubic_report(equation, setup, names)
	}
}

// read_iapws95 builds water's fluid by IAPWS-95 (M5a's D12): format 4's,
// for water alone, with none of a cubic's keys.
fn read_iapws95(mut s Section, db components.Database, names []string, names_line int, format int) !Modelled {
	if format < 4 {
		return s.fail('equation', 'IAPWS-95 needs format 4 (M5a), and this file is format ${format}',
			[])
	}
	for key in ['alpha', 'kij', 'kij_overrides'] {
		if s.has(key) {
			return s.fail(key, "${key} is a cubic's, and IAPWS-95 has none", [])
		}
	}
	s.finish()!
	if names.len != 1 {
		return file_error(names_line, 'components.names', "IAPWS-95 is water's alone; water with other compounds needs a mixture model, which reflux doesn't have yet",
			[])
	}
	setup := db.iapws95(names[0]) or {
		if err is components.LookupError {
			return file_error(names_line, 'components.names', 'no compound ${err.query}', err.suggestions)
		}
		return file_error(names_line, 'components.names', err.msg(), [])
	}
	return Modelled{
		fluid:  setup.fluid
		report: ModelReport{
			equation: 'IAPWS-95'
			notes:    ["IAPWS R6-95(2018), with its own gas constant; H and S on reflux's reference"]
		}
	}
}

// cubic_report says how an equation of state's fluid was built: each pair's
// kᵢⱼ with its source, and each α not the equation's own.
fn cubic_report(equation string, setup components.Setup, names []string) ModelReport {
	mut pairs := []PairReport{}
	for p in setup.kij.pairs {
		pairs << PairReport{
			a:          names[p.i]
			b:          names[p.j]
			parameters: [ParameterReport{'kij', p.k.str()}]
			source:     match p.source {
				.file {
					if used := p.used {
						'${setup.kij.bank} line ${used.line}: ${used.comment}'
					} else {
						setup.kij.bank
					}
				}
				.override {
					'the file'
				}
				.missing {
					'none listed: 0'
				}
			}
		}
	}
	mut notes := []string{}
	for i, how in setup.alphas {
		if how != "the equation's own" {
			notes << '${names[i]}: α ${how}'
		}
	}
	return ModelReport{
		equation: equation
		pairs:    pairs
		notes:    notes
	}
}

// The parameters of a binary model's pair, as a file names them (D12).
fn pair_keys(liquid components.LiquidModel) []string {
	return match liquid {
		.nrtl { ['tau_ab', 'tau_ba', 'alpha'] }
		.uniquac { ['tau_ab', 'tau_ba'] }
		.wilson { ['lambda_ab', 'lambda_ba'] }
		else { []string{} }
	}
}

// read_gamma_phi builds a γ–φ fluid from format 4's [model] (D12): liquid,
// the activity model; vapour, "ideal gas" or an equation; poynting; henry,
// solutes beyond those above their Tc; missing, "refuse" or "ideal"; and
// [[model.pairs]], a binary model's pairs in place of the bank's.
fn read_gamma_phi(mut s Section, db components.Database, names []string, names_line int) !Modelled {
	if s.has('equation') {
		return s.fail('equation', 'a γ–φ fluid has its liquid by liquid and its vapour by vapour, not equation',
			[])
	}
	for key in ['alpha', 'kij_overrides'] {
		if s.has(key) {
			return s.fail(key, "${key} is an equation of state's: a γ–φ fluid's vapour takes its equation's own α and the bank's kᵢⱼ",
				[])
		}
	}
	liquid_name := s.text('liquid')!
	liquid := liquid_models[liquid_name] or {
		return s.fail('liquid', 'no activity model ${liquid_name}', quantity.closest(liquid_name,
			liquid_models.keys()))
	}
	mut vapour := ?thermo.CubicKind(none)
	mut vapour_name := 'ideal gas'
	if s.has('vapour') {
		vapour_name = s.text('vapour')!
		if vapour_name != 'ideal gas' {
			vapour = equations[vapour_name] or {
				mut all := ['ideal gas']
				all << equations.keys()
				return s.fail('vapour', 'no vapour ${vapour_name}', quantity.closest(vapour_name,
					all))
			}
		}
	}
	if s.has('kij') && vapour == none {
		return s.fail('kij', "kij is the vapour's equation's, and the vapour is the ideal gas", [])
	}
	kij := if s.has('kij') { s.boolean('kij')! } else { true }
	poynting := if s.has('poynting') { s.boolean('poynting')! } else { true }
	henry := if s.has('henry') { s.strings('henry')! } else { []string{} }
	for name in henry {
		if name !in names {
			return s.fail('henry', '${name} is not one of the components', quantity.closest(name,
				names))
		}
	}
	missing := if s.has('missing') { s.text('missing')! } else { 'refuse' }
	if missing !in ['refuse', 'ideal'] {
		return s.fail('missing', 'missing is "refuse" or "ideal", not ${missing}', quantity.closest(missing,
			['refuse', 'ideal']))
	}
	keys := pair_keys(liquid)
	mut overrides := []components.PairOverride{}
	if s.has('pairs') {
		if keys.len == 0 {
			return s.fail('pairs', "UNIFAC's parameters are its groups', from its table: pairs are NRTL's, UNIQUAC's and Wilson's",
				[])
		}
		for mut t in s.tables('pairs')! {
			a := t.text('a')!
			b := t.text('b')!
			for k, name in [a, b] {
				if name !in names {
					return t.fail(['a', 'b'][k], '${name} is not one of the components',
						quantity.closest(name, names))
				}
			}
			if a == b {
				return t.fail('b', 'a pair is of two components, and this one is ${a} twice', [])
			}
			mut texts := []string{}
			for key in keys {
				text := t.text(key)!
				what := match key {
					'alpha' { 'α' }
					'lambda_ab', 'lambda_ba' { 'Λ' }
					else { 'τ' }
				}
				activity.plain(text, what) or { return t.fail(key, err.msg(), []) }
				texts << text
			}
			t.finish()!
			overrides << components.PairOverride{
				a:     (db.find(a) or { return t.fail('a', err.msg(), []) }).ids.cas
				b:     (db.find(b) or { return t.fail('b', err.msg(), []) }).ids.cas
				texts: texts
			}
		}
	}
	s.finish()!
	setup := db.gamma_phi(names,
		liquid:   liquid
		vapour:   vapour
		kij:      kij
		poynting: poynting
		henry:    henry
		ideal:    missing == 'ideal'
		pairs:    overrides
	) or {
		if err is components.LookupError {
			return file_error(names_line, 'components.names', 'no compound ${err.query}', err.suggestions)
		}
		if err.msg().starts_with('a pair ') || err.msg().starts_with('the pair ') {
			return s.fail('pairs', err.msg(), [])
		}
		return s.fail('liquid', err.msg(), [])
	}
	return Modelled{
		fluid:     setup.fluid
		report:    gamma_phi_report(liquid_name, vapour_name, poynting, setup, names, keys)
		gamma_phi: setup
	}
}

// gamma_phi_report says how a γ–φ fluid was built: its solutes, each pair's
// parameters with their source, its notes, and a warning for each pair given
// none (D13).
fn gamma_phi_report(liquid string, vapour string, poynting bool, setup components.GammaPhiSetup, names []string, keys []string) ModelReport {
	cas := setup.records.map(it.ids.cas)
	bank := match liquid {
		'NRTL' { 'nrtl.ipd' }
		'UNIQUAC' { 'uniquac.ipd' }
		else { 'wilson.ipd' }
	}
	mut pairs := []PairReport{}
	mut warnings := []string{}
	for u in setup.pairs {
		a := names[setup.model.solvents[u.i]]
		b := names[setup.model.solvents[u.j]]
		source := match u.source {
			.bank {
				if l := u.listing {
					'${bank} line ${l.line}: ${l.comment}'
				} else {
					bank
				}
			}
			.file {
				'the file'
			}
			.ideal {
				'none listed: an ideal pair'
			}
			.no_residual {
				'none listed: no residual interaction'
			}
		}
		if u.source in [.ideal, .no_residual] {
			warnings << '${a} with ${b}: ${bank} has no pair, so ${source.all_after(': ')} (missing = "ideal")'
		}
		mut parameters := []ParameterReport{}
		for k, key in keys {
			parameters << ParameterReport{key, u.texts[k] or { '' }}
		}
		pairs << PairReport{
			a:          a
			b:          b
			parameters: parameters
			source:     source
		}
	}
	return ModelReport{
		liquid:   liquid
		vapour:   vapour
		poynting: poynting
		solutes:  setup.solutes.map(names[cas.index(it)])
		pairs:    pairs
		notes:    setup.notes
		warnings: warnings
	}
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

// read_unit builds one unit and adds it to the flowsheet, and returns its type.
fn read_unit(mut s Section, name string, makers map[string]Maker, mut sheet Flowsheet) !string {
	u, kind, ports := make_unit(mut s, name, makers)!
	sheet.add_unit(u, ports) or { return s.located(err) }
	return kind
}

// make_unit builds a unit from its table by its type's maker: the unit, its
// type, and the stream on each port.
fn make_unit(mut s Section, name string, makers map[string]Maker) !(Unit, string, map[string]string) {
	kind := s.text('type')!
	maker := makers[kind] or {
		return s.fail('type', 'no unit type ${kind}', quantity.closest(kind, makers.keys()))
	}
	s.name = name
	u := maker(mut s) or { return s.located(err) }
	ports := s.ports(u)!
	s.finish()!
	return u, kind, ports
}
