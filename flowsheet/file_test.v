// The file format (exit criterion 6, and the rest of criterion 1): a file builds
// the same flowsheet as code does, loading, writing and loading again changes
// nothing, every wrong value is refused with its line and reason and with
// suggestions for a misspelt name, and fuzzed files never panic.
import components
import flowsheet { Flowsheet, FlowsheetError, Loaded, Section, Solution, Unit }
import math
import os
import thermo

const db = components.load_default() or { panic(err) }

fn read_flowsheet(name string) string {
	path := os.join_path(@VMODROOT, 'testdata', 'flowsheets', name)
	return os.read_file(path) or { panic(err) }
}

fn reference_text() string {
	return read_flowsheet('reference.toml')
}

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

// same_numbers asserts that two solutions are the same to the last bit.
fn same_numbers(a Solution, b Solution, at string) {
	assert a.streams.len == b.streams.len, at
	for name, x in a.streams {
		y := b.streams[name] or { panic('${at}: no ${name}') }
		assert x.t == y.t && x.p == y.p && x.h == y.h && x.flows == y.flows, '${at}: ${name}'
	}
}

fn test_a_file_builds_what_code_builds() {
	// The same flowsheet in V, unit by unit: every stream is the same to the last
	// bit, since nothing about a file should change a number.
	l := load(reference_text())!
	fl := db.fluid(.pr, ['methane', 'ethane', 'propane', 'n-butane', 'n-pentane'])!
	mut f := Flowsheet{}
	z := [0.7, 0.1, 0.08, 0.07, 0.05]
	mut flows := []f64{}
	for x in z {
		flows << x * 100.0
	}
	f.add_feed('feed', flowsheet.feed_tp(fl, flows, 300.0, 30e5)!)!
	f.add_unit(flowsheet.new_mixer('mixer', 2)!, {
		'in1': 'feed'
		'in2': 'recycle'
		'out': 'mixed'
	})!
	f.add_unit(flowsheet.new_heater('cooler', mode: .outlet_t, value: 230.0)!, {
		'in':  'mixed'
		'out': 'cold'
	})!
	f.add_unit(flowsheet.new_flash_drum('drum', mode: .duty, value: 0.0)!, {
		'in1':    'cold'
		'vapour': 'gas'
		'liquid': 'liquid'
	})!
	f.add_unit(flowsheet.new_splitter('split', [0.5, 0.5])!, {
		'in':   'liquid'
		'out1': 'back'
		'out2': 'on'
	})!
	f.add_unit(flowsheet.new_pump('pump', mode: .outlet_p, value: 31e5, efficiency: 0.75)!, {
		'in':  'back'
		'out': 'recycle'
	})!
	f.add_unit(flowsheet.new_valve('valve', mode: .outlet_p, value: 5e5)!, {
		'in':  'on'
		'out': 'let down'
	})!
	f.add_unit(flowsheet.new_flash_drum('second', mode: .duty, value: 0.0)!, {
		'in1':    'let down'
		'vapour': 'flash gas'
		'liquid': 'product'
	})!
	same_numbers(l.solve()!, f.solve(fl)!, 'file and code')
}

// every_unit has one of each unit type, each specified a different way from
// the reference flowsheet.
const every_unit = read_flowsheet('every_unit.toml')

fn test_every_unit_type_loads_and_solves() {
	l := load(every_unit)!
	assert l.sheet.units.len == 9
	assert l.params.method == .direct && l.params.tol == 1e-8
	assert l.params.max_iterations == 200 && l.params.tears == ['back']
	assert l.params.estimates['back'].t == 300.0
	// The model is what the file says: SRK, with its one kᵢⱼ set by hand.
	cas := ['74-82-8', '74-98-6', '109-66-0']
	mine := db.fluid(.srk, cas,
		kij:       false
		overrides: [components.Kij{
			cas1: cas[0]
			cas2: cas[1]
			k:    0.01
		}]
	)!
	assert flowsheet.feed_tp(mine, [1.0, 1.0, 1.0], 250.0, 20e5)!.h == flowsheet.feed_tp(l.fluid, [
		1.0,
		1.0,
		1.0,
	], 250.0, 20e5)!.h
	s := l.solve()!
	// The duty is what the file says.
	warm := s.units.filter(it.name == 'warm')[0].energy[0].rate
	assert math.abs(warm + 10e3) <= 1e-5, '${warm}'
	assert s.plan.blocks.filter(it.tears.len > 0)[0].tears == ['back']
	// The kij override reached the model, and the units their specifications.
	assert s.streams['eased'].p == 4.5e5
	assert s.streams['top'].p == 2e5 && s.streams['top'].t == 290.0
	assert s.streams['lifted'].p == 12e5
	assert math.abs(s.streams['warmed'].p - 19.5e5) <= 1e-6
	assert math.abs(s.streams['turned'].p - 7.5e5) <= 1e-6
	// 25% of what the separator's rich outlet carries goes back.
	assert math.abs(s.streams['back'].total() - 0.25 * s.streams['eased'].total()) <= 1e-7
	// In the second feed, 1800 kg/h is half a kg/s, half of it each compound.
	mw := l.fluid.model.components()
	oil := l.sheet.feeds[1].stream
	assert math.abs(oil.flows[1] - 0.25 / mw[1].mw) <= 1e-10
	assert math.abs(oil.flows[2] - 0.25 / mw[2].mw) <= 1e-10
	assert oil.flows[0] == 0.0 && oil.t == 300.15 && oil.beta == 0.0
}

// stream_flows loads a feed of four bases into one flowsheet.
fn with_feed(lines string) string {
	base := reference_text()
	head := base.all_before('[feed.feed]')
	tail := '[unit.mixer]' + base.all_after('[unit.mixer]')
	return '${head}[feed.feed]\nT = "300 K"\nP = "30 bar"\n${lines}\n\n${tail}'
}

fn test_each_basis_gives_the_same_flows() {
	// One mixture stated four ways (D9), with the molar masses from the file's own
	// fluid: methane and n-pentane, as moles of 0.75 and 0.25, say.
	l := load(reference_text())!
	mw := l.fluid.model.components().map(it.mw)
	mole := [0.75, 0.0, 0.0, 0.0, 0.25]
	f := 8.0 // mol/s
	mut mass := 0.0
	for i, x in mole {
		mass += x * mw[i]
	}
	kg := f * mass
	m0, m4 := 0.75 * mw[0] / mass, 0.25 * mw[4] / mass
	texts := [
		'flow = "8 mol/s"\nmole_fractions = { methane = 0.75, n-pentane = 0.25 }',
		'flow = "${kg} kg/s"\nmole_fractions = { methane = 0.75, n-pentane = 0.25 }',
		'flow = "8 mol/s"\nmass_fractions = { methane = ${m0}, n-pentane = ${m4} }',
		'flow = "${kg} kg/s"\nmass_fractions = { methane = ${m0}, n-pentane = ${m4} }',
		'molar_flows = { methane = "6 mol/s", n-pentane = "2 mol/s" }',
		'mass_flows = { methane = "${6.0 * mw[0]} kg/s", n-pentane = "${2.0 * mw[4]} kg/s" }',
	]
	for t in texts {
		s := load(with_feed(t))!.sheet.feeds[0].stream
		for i, want in [6.0, 0.0, 0.0, 0.0, 2.0] {
			assert math.abs(s.flows[i] - want) <= 1e-12, '${t}: ${i}: ${s.flows[i]}'
		}
	}
	// Fractions a little off 1 are taken as they stand, shared out of the flow.
	off := load(with_feed('flow = "8 mol/s"\nmole_fractions = { methane = 0.7500004, n-pentane = 0.25 }'))!
	assert math.abs(off.sheet.feeds[0].stream.total() - 8.0) <= 1e-12
	// Each state a feed can be given by (D2).
	for t in [
		'[feed.feed]\nP = "30 bar"\nH = "-90 kJ/mol"\n',
		'[feed.feed]\nT = "230 K"\nvapour_fraction = 0.5\n',
		'[feed.feed]\nP = "30 bar"\nvapour_fraction = "50 %"\n',
	] {
		text := reference_text().replace('[feed.feed]\nT = "300 K"\nP = "30 bar"\n', t)
		assert text != reference_text()
		load(text)!
	}
}

fn test_load_write_load_is_the_same_flowsheet() {
	// Criterion 6. Writing is a fixed point after the first, and the flowsheet
	// solves to the same numbers exactly, in the units the file used.
	for text in [reference_text(), every_unit] {
		a := load(text)!
		written := a.write()
		b := load(written)!
		assert b.write() == written
		same_numbers(a.solve()!, b.solve()!, 'written')
		assert b.params.method == a.params.method && b.params.tol == a.params.tol
		assert b.params.tears == a.params.tears
		assert b.params.max_iterations == a.params.max_iterations
		assert b.params.estimates.len == a.params.estimates.len
		assert b.sheet.units.len == a.sheet.units.len
		assert written.starts_with('format = 1\n')
	}
	written := load(every_unit)!.write()
	assert written.contains('[unit.share]')
	assert written.contains('fractions = ["25 %", "75 %"]')
	assert written.contains('flow = "1800 kg/h"')
	assert written.contains('T = "27 degC"')
}

// refused asserts that a file fails as a bad file, with this line and reason.
fn refused(text string, line int, key string, want string, suggestion string) {
	load(text) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		assert e.key == key, 'key ${e.key}, not ${key}: ${err.msg()}'
		assert e.reason.contains(want), err.msg()
		if suggestion != '' {
			assert suggestion in e.suggestions, err.msg()
			assert err.msg().contains('did you mean'), err.msg()
		}
		return
	}
	assert false, 'accepted: ${want}'
}

// lined is the 1-based line of the first occurrence of a text in the reference.
fn lined(text string) int {
	return reference_text().all_before(text).count('\n') + 1
}

fn sub(a string, b string) string {
	t := reference_text()
	assert t.count(a) == 1, a
	return t.replace(a, b)
}

fn test_a_wrong_field_is_refused_with_its_line_and_reason() {
	feed := '[feed.feed]\n'
	// The format.
	refused(sub('format = 1', 'format = 2'), lined('format = 1'), 'format', 'reads format 1',
		'')
	refused(sub('format = 1', 'format = "1"'), lined('format = 1'), 'format', 'a whole number',
		'')
	refused(sub('format = 1', 'format = 1.5'), lined('format = 1'), 'format', 'a whole number',
		'')
	refused(sub('format = 1\n', ''), 0, '', 'format is missing', '')
	// Components.
	names := lined('names =')
	refused(sub('"n-pentane"]', '"n-pentan"]'), names, 'components.names', 'no compound n-pentan',
		'N-pentane (109-66-0)')
	refused(sub('["methane", "ethane", "propane", "n-butane", "n-pentane"]', '[]'), lined('[components]'),
		'components.names', 'at least one component', '')
	refused(sub('"ethane", ', '"methane", '), names, 'components.names', 'listed twice', '')
	refused(sub('names = [', 'nmaes = ['), lined('[components]'), 'components', 'names is missing',
		'')
	refused(sub('[model]\n', 'color = "blue"\n\n[model]\n'), names + 2, 'components.color', 'no key color',
		'')
	// The model.
	eq := lined('equation =')
	refused(sub('equation = "PR"', 'equation = "Pr"'), eq, 'model.equation', 'no equation Pr',
		'PR')
	refused(sub('equation = "PR"', 'equation = 1'), eq, 'model.equation', 'should be a string',
		'')
	refused(sub('equation = "PR"', 'equation = "PR"\nalpha = "soave"'), eq + 1, 'model.alpha',
		'no alpha soave', '')
	refused(sub('equation = "PR"', 'equation = "PR"\nkij = "yes"'), eq + 1, 'model.kij', 'true or false',
		'')
	refused(sub('equation = "PR"', 'equation = "SRK"\nalpha = "mathias_copeman"'), eq, 'model.equation',
		'fitted for PR', '')
	refused(sub('equation = "PR"', 'equation = "PR"\nkij_overrides = [{ a = "methane", b = "ethan", k = 0.1 }]'),
		eq + 1, 'model.kij_overrides[0].a', 'ethan is not one of the components', 'ethane')
	refused(sub('equation = "PR"', 'equation = "PR"\nequaton = "SRK"'), eq + 1, 'model.equaton',
		'no key equaton', 'equation')
	// The solver.
	m := lined('method =')
	refused(sub('method = "wegstein"', 'method = "wegstien"'), m, 'solver.method', 'no method wegstien',
		'wegstein')
	refused(sub('method = "wegstein"', 'tol = -1.0'), m, 'solver.tol', 'a tolerance of', '')
	refused(sub('method = "wegstein"', 'max_iterations = 0'), m, 'solver.max_iterations', '0 iterations',
		'')
	refused(sub('method = "wegstein"', 'max_iterations = 2.5'), m, 'solver.max_iterations', 'a whole number',
		'')
	refused(sub('method = "wegstein"', 'tears = ["nothing"]'), 0, '', 'tear nothing is not a stream',
		'')
	// A feed.
	t := lined('T = "300 K"')
	refused(sub('T = "300 K"', 'T = 300'), t, 'feed.feed.T', 'needs a unit', '')
	refused(sub('T = "300 K"', 'T = "300"'), t, 'feed.feed.T', 'needs a unit', '')
	refused(sub('T = "300 K"', 'T = "300 k"'), t, 'feed.feed.T', 'no unit k', 'K')
	refused(sub('T = "300 K"', 'T = "30 bar"'), t, 'feed.feed.T', 'is a unit of pressure', '')
	refused(sub('T = "300 K"', 'T = "-300 degC"'), t, 'feed.feed.T', 'absolute zero', '')
	refused(sub('T = "300 K"', 'Temp = "300 K"'), lined(feed), 'feed.feed', 'a feed needs T and P',
		'')
	refused(sub('T = "300 K"', 'T = "300 K"\nH = "1 kJ/mol"'), lined(feed), 'feed.feed', 'a feed needs T and P',
		'')
	refused(sub('P = "30 bar"', 'H = "1 kJ/mol"'), lined(feed), 'feed.feed', 'a feed needs T and P',
		'')
	refused(sub('T = "300 K"\nP = "30 bar"', 'H = "1 kJ/mol"\nvapour_fraction = 0.5'), lined(feed),
		'feed.feed', 'a feed needs T and P', '')
	refused(sub('T = "300 K"', 'T = "300 K"\nnotes = "hot"'), t + 1, 'feed.feed.notes', 'no key notes',
		'')
	refused(sub('flow = "100 mol/s"', 'flow = "100 W"'), lined('flow ='), 'feed.feed.flow', 'is a unit of power',
		'')
	refused(sub('flow = "100 mol/s"\n', ''), lined(feed), 'feed.feed', 'flow is missing', '')
	fractions := lined('mole_fractions')
	refused(sub('methane = 0.7', 'methane = 0.8'), fractions, 'feed.feed.mole_fractions', 'they sum to',
		'')
	refused(sub('methane = 0.7', 'methan = 0.7'), fractions, 'feed.feed.mole_fractions', 'methan is not one of the components',
		'methane')
	refused(sub('mole_fractions', 'mole_fraction'), lined(feed), 'feed.feed', 'a specification is missing',
		'')
	refused(sub('flow = "100 mol/s"', 'flow = "100 mol/s"\nmolar_flows = { methane = "1 mol/s" }'),
		fractions, 'feed.feed.molar_flows', 'give only one of', '')
	refused(sub('methane = 0.7', 'methane = "x"'), fractions, 'feed.feed.mole_fractions', 'not a number',
		'')
	refused(sub('P = "30 bar"', 'P = "-5 bar"'), lined('P = "30 bar"'), 'feed.feed.P', 'not above zero',
		'')
	// A unit.
	u := lined('[unit.cooler]')
	refused(sub('type = "heater"', 'type = "heatr"'), u + 1, 'unit.cooler.type', 'no unit type heatr',
		'heater')
	refused(sub('type = "heater"', 'type = 1'), u + 1, 'unit.cooler.type', 'should be a string',
		'')
	refused(sub('type = "heater"\n', ''), u, 'unit.cooler', 'type is missing', '')
	refused(sub('outlet_T = "230 K"', 'outlet_t = "230 K"'), u, 'unit.cooler', 'a specification is missing',
		'')
	refused(sub('outlet_T = "230 K"', 'outlet_T = "230 K"\nduty = "1 W"'), u + 5, 'unit.cooler.duty',
		'give only one of', '')
	refused(sub('outlet_T = "230 K"', 'outlet_T = "230 K"\ncolor = "red"'), u + 5, 'unit.cooler.color',
		'no key color', '')
	refused(sub('out = "cold"', 'out = ["cold", "other"]'), u + 3, 'unit.cooler.out', 'needs 1 streams',
		'')
	refused(sub('in = "mixed"\n', ''), u, 'unit.cooler', 'in is missing', '')
	refused(sub('out = "cold"\n', ''), u, 'unit.cooler', 'out is missing', '')
	refused(sub('efficiency = 0.75', 'efficiency = 1.5'), lined('[unit.pump]'), 'unit.pump', 'an efficiency of',
		'')
	refused(sub('efficiency = 0.75\n', ''), lined('[unit.pump]'), 'unit.pump', 'efficiency is missing',
		'')
	refused(sub('fractions = [0.5, 0.5]', 'fractions = [0.5, 0.4]'), lined('[unit.split]'), 'unit.split',
		'sum to 0.9', '')
	refused(sub('fractions = [0.5, 0.5]', 'fractions = 0.5'), lined('fractions = ['), 'unit.split.fractions',
		'a list of fractions', '')
	// A temperature rise isn't a temperature, and a cooler's is below zero.
	cooled := load(sub('outlet_T = "230 K"', 'delta_T = "-70 K"'))!
	solved := cooled.solve()!
	assert math.abs(solved.streams['mixed'].t - solved.streams['cold'].t - 70.0) <= 1e-9
	// A header with a comment after it still has its line.
	refused(sub('[unit.cooler]\n', '[unit.cooler] # the cooler\n').replace('outlet_T = "230 K"',
		'outlet_t = "230 K"'), u, 'unit.cooler', 'a specification is missing', '')
	// An estimate tears its stream, and that stream has to be in a loop.
	refused(sub('[unit.mixer]', '[estimate.on]\nT = "300 K"\nP = "1 bar"\nmolar_flows = { methane = "1 mol/s" }\n\n[unit.mixer]'),
		0, '', 'tear on is in no loop', '')
	// A flowsheet that doesn't join up, from the unit at fault.
	refused(sub('out = "cold"', 'out = "mixed"'), lined('[unit.cooler]'), 'unit.cooler', 'two sources',
		'')
	refused(sub('in = "on"', 'in = "back"'), lined('[unit.valve]'), 'unit.valve', 'goes to two inlets',
		'')
	refused(sub('in = "on"', 'in = "nothing"'), lined('[unit.valve]'), 'unit.valve', 'comes from no feed',
		'')
}

fn test_the_top_level_is_checked_too() {
	// A table where none should be, misspelt.
	refused(sub('[model]\n', '[modle]\nequation = "PR"\n\n[model]\n'), lined('[model]') + 1, 'modle',
		'no key modle', 'model')
	// A value where a table should be.
	refused(sub('[feed.feed]\nT = "300 K"\nP = "30 bar"\nflow = "100 mol/s"\nmole_fractions = { methane = 0.7, ethane = 0.1, propane = 0.08, n-butane = 0.07, n-pentane = 0.05 }\n',
		'[feed]\nfeed = 1\n'), lined('[feed.feed]') + 1, 'feed.feed', 'should be a table', '')
}

fn test_an_empty_file_is_refused() {
	refused('', 0, '', 'format is missing', '')
	refused('# nothing but a comment\n', 0, '', 'format is missing', '')
}

fn test_names_with_quotes_and_spaces_survive_a_round_trip() {
	// A name is whatever TOML allows, and writing mustn't break it.
	// The stream's name has a backslash, a tab and quotes in it.
	name := r'"c:\\temp\t \"x\""'
	text := sub('[unit.cooler]', '[unit."the cold one"]').replace('out = "cold"', 'out = ${name}').replace('in = "cold"',
		'in = ${name}')
	a := load(text)!
	written := a.write()
	b := load(written)!
	assert b.write() == written
	assert a.sheet.units.map(it.unit.name) == b.sheet.units.map(it.unit.name)
	same_numbers(a.solve()!, b.solve()!, 'odd names')
	assert a.solve()!.streams['c:\\temp\t "x"'].t == 230.0
	// The same for a literal string, whose backslashes are its own.
	lit := sub('out = "mixed"', "out = 'c:\\users'").replace('in = "mixed"', "in = 'c:\\users'")
	c := load(lit)!
	assert c.write().contains("out = 'c:\\users'")
	d := load(c.write())!
	same_numbers(c.solve()!, d.solve()!, 'a literal string')
}

fn test_a_syntax_error_has_its_line() {
	load('format = 1\n\n[components\n') or {
		e := err as FlowsheetError
		assert e.kind == .bad_file, err.msg()
		return
	}
	assert false, 'accepted'
}

fn test_a_key_given_twice_is_refused() {
	// Exactly one of each (D3): the parser refuses the same key twice.
	load(sub('outlet_T = "230 K"', 'outlet_T = "230 K"\noutlet_T = "240 K"')) or {
		assert err is FlowsheetError
		assert (err as FlowsheetError).kind == .bad_file
		return
	}
	assert false, 'accepted'
}

fn test_what_would_bring_the_parser_down_is_refused_first() {
	// V 0.5.2's toml checker runs off the end of a string with a \u escape the
	// right distance from its end. It's a panic, which nothing can catch, so
	// check_text refuses every \u before the parser sees it, valid or not.
	refused('format = 1\na = "' + r'\u12recycle' + '"\n', 2, '', 'escape', '')
	refused('format = 1\na = "' + r'\u00e9abcdefg' + '"\n', 2, '', 'escape', '')
	refused('format = 1\n\n\na = """\n\n' + r'\U0001F600' + '"""\n', 6, '', 'escape', '')
	// Comments and literal strings don't count, so a Windows path is fine.
	load(reference_text() + '\n# ' + r'"C:\users\me"' + '\n')!
	load(sub('out = "mixed"', "out = '" + r'"\users' + "'").replace('in = "mixed"', "in = '" +
		r'"\users' + "'"))!
	load(sub('[model]\n', "[model]\n# 'a comment's quote\n"))!
	text := sub('names = [', "data = ['" + r'.' + "']\nnames = [")
	load(text)!
	// And a depth no flowsheet comes near, which would overflow the stack.
	deep := 'format = 1\na = ' + '['.repeat(100000)
	refused(deep, 2, '', 'nests more than', '')
	refused('format = 1\na = ' + '{ a = '.repeat(50000), 2, '', 'nests more than', '')
}

fn test_the_database_given_is_the_one_used() {
	// With an empty one, no compound is found.
	dir := os.join_path(os.vtmp_dir(), 'reflux_empty_db_${os.getpid()}')
	os.mkdir_all(dir)!
	defer {
		os.rmdir_all(dir) or {}
	}
	empty := components.load([dir])!
	flowsheet.load_text(reference_text(), db: empty) or {
		assert err.msg().contains('no compound methane'), err.msg()
		return
	}
	assert false, 'found methane in an empty database'
}

fn test_your_own_component_data_overrides_the_database() {
	// data is relative to the file, and read after the built-in data.
	dir := os.join_path(os.vtmp_dir(), 'reflux_file_test_${os.getpid()}')
	os.mkdir_all(os.join_path(dir, 'mine'))!
	defer {
		os.rmdir_all(dir) or {}
	}
	methane := os.join_path(@VMODROOT, 'data', 'components', '74-82-8.toml')
	original := os.read_file(methane)!
	line := original.split_into_lines().filter(it.starts_with('critical_temperature'))
	assert line.len <= 1
	os.write_file(os.join_path(dir, 'mine', '74-82-8.toml'), original)!
	path := os.join_path(dir, 'plant.toml')
	os.write_file(path, sub('names = [', 'data = ["mine"]\nnames = ['))!
	l := flowsheet.load(path)!
	assert l.fluid.model.components().len == 5
	if _ := flowsheet.load(os.join_path(dir, 'no such file.toml')) {
		assert false, 'read a file that is not there'
	} else {
		assert err.msg().contains("can't be read"), err.msg()
	}
}

// Doubler is a unit from outside this module, as a later family's will be: it
// sends its inlet straight out.
struct Doubler {
	name string
}

fn (d Doubler) inlet_ports() []string {
	return ['in']
}

fn (d Doubler) outlet_ports() []string {
	return ['out']
}

fn (d Doubler) calculate(_ thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	return flowsheet.Outcome{
		outlets: [inlets[0]]
	}
}

fn make_doubler(mut s Section) !Unit {
	return Doubler{s.name}
}

fn test_a_later_unit_family_plugs_in() {
	// D1: the loader takes a table of makers, and needn't import a family's module.
	text := sub('type = "valve"\nin = "on"\nout = "let down"\noutlet_P = "5 bar"', 'type = "passes"\nin = "on"\nout = "let down"')
	refused(text, lined('[unit.valve]') + 1, 'unit.valve.type', 'no unit type passes', '')
	mut makers := map[string]flowsheet.Maker{}
	makers['passes'] = make_doubler
	l := flowsheet.load_text(text, db: db, makers: makers)!
	assert l.sheet.units.map(it.unit.name).contains('valve')
	assert l.solve()!.streams['let down'].p == 30e5
}

// Rng is a small seeded generator for the fuzzing: splitmix64.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) below(n int) int {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return int((z ^ (z >> 31)) % u64(n))
}

const pieces = ['"', "'", '[', ']', '[[', ']]', '{', '}', '=', ',', '.', '#', '\n', ' ', '\t',
	r'\ ', r'\u12', '"""', "'''", '1e', '0x', '-', '+', '_', '1979-05-27', 'T', ':', 'true', 'nan',
	'inf', 'é', '\r\n', '\x7f', 'a.b', '[a]', '[a.b]', 'x = 1', '[unit.x]', 'type = "mixer"',
	'[feed.x]', '"0 K"', '"-1 bar"', '1e400', '-1', '0', '[]', '{}', '""', 'in = "feed"', 'out = "gas"']

fn mutate(mut r Rng, text string) string {
	mut s := text
	for _ in 0 .. 1 + r.below(4) {
		i := r.below(s.len)
		j := i + r.below(12)
		e := if j < s.len { j } else { s.len }
		p := pieces[r.below(pieces.len)]
		s = match r.below(4) {
			0 { s[..i] + p + s[i..] }
			1 { s[..i] + s[e..] }
			2 { s[..i] + p + s[e..] }
			else { s[..i] + s[i..e] + s[i..] }
		}
		if s == '' {
			s = 'x'
		}
	}
	return s
}

fn test_fuzzed_files_never_panic() {
	// 4,000 seeded mutations of the reference flowsheet and of every_unit: any
	// answer is fine as long as there is one, and a failure is typed. A panic
	// would end the test. Whatever loads is written and loaded again.
	mut r := Rng{20261103}
	mut accepted, mut refusals := 0, 0
	bases := [reference_text(), every_unit]
	for _ in 0 .. 4000 {
		t := mutate(mut r, bases[r.below(2)])
		l := load(t) or {
			assert err is FlowsheetError, '${t}\n${err.msg()}'
			refusals++
			continue
		}
		accepted++
		written := l.write()
		again := load(written) or {
			assert false, 'it loaded, and what it wrote did not: ${err.msg()}\n${written}'
			continue
		}
		assert again.write() == written
	}
	println('fuzz: ${accepted} accepted, ${refusals} refused')
	assert accepted > 0 && refusals > 0
}
