// Format 4 (M6 step 9, criterion 10): a γ–φ fluid's [model], its pairs as
// expressions, the decanter, and report format 4's model; each refusal at
// its line. testdata/flowsheets/decanter.toml is the file most of it edits.
import components
import flash
import flowsheet { FlowsheetError, Loaded }
import math
import os
import thermo

const db = components.load_default() or { panic(err) }

fn decanter_text() string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'decanter.toml')) or {
		panic(err)
	}
}

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

// with is the decanter file with a text replaced, which must be there once.
fn with(a string, b string) string {
	t := decanter_text()
	assert t.count(a) == 1, a
	return t.replace(a, b)
}

// lined is the 1-based line of a text's first occurrence in the decanter file.
fn lined(text string) int {
	return decanter_text().all_before(text).count('\n') + 1
}

fn refused(text string, line int, key string, want string, suggestion string) {
	load(text) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		assert e.key == key, 'key ${e.key}, not ${key}: ${err.msg()}'
		assert e.reason.contains(want), err.msg()
		if suggestion != '' {
			assert suggestion in e.suggestions, err.msg()
		}
		return
	}
	assert false, 'accepted: ${want}'
}

// The decanter's outlets are the PT flash's phases at the heater's outlet,
// each by its label out of its own port, in the fluid the database builds
// for the same names: the file changes no number. The decanter is
// adiabatic, so its PH flash lands on the heater's T, to its tolerance.
fn test_a_decanter_sends_each_phase_out_of_its_port() {
	l := load(decanter_text())!
	sol := l.solve()!
	fl := db.gamma_phi(['water', '1-butanol', 'methanol'])!.fluid
	r := flash.pt(fl.model, 364.2, 101325.0, [0.75, 0.2, 0.05])!
	assert r.phases.map(it.label) == [thermo.PhaseLabel.vapor, .liquid, .liquid2]
	hot := sol.streams['hot'] or { panic('no hot') }
	assert hot.result.phases.len == 3 && hot.beta == r.beta
	for k, name in ['vapour', 'organic', 'aqueous'] {
		s := sol.streams[name] or { panic('no ${name}') }
		assert math.abs(s.t - 364.2) <= 1e-7, '${name}: ${s.t}'
		assert math.abs(s.total() - 10.0 * r.phases[k].fraction) <= 1e-8, name
		for i in 0 .. 3 {
			assert math.abs(s.flows[i] / s.total() - r.phases[k].state.n[i]) <= 1e-8, name
		}
		// One phase present, the others kept at no moles, as a drum's are.
		assert s.result.phases.filter(it.fraction > 0.0).map(it.label) == [r.phases[k].label]
	}
	// The report (format 4's model, now format 5) has three phases in the hot
	// stream and the model's every pair from ChemSep's bank.
	rep := l.report(sol)
	assert rep.format == 5
	streams := rep.streams.filter(it.name == 'hot')
	assert streams[0].phases.map(it.phase) == ['vapour', 'liquid', 'liquid2']
	m := rep.model or { panic('no model') }
	assert m.liquid == 'NRTL' && m.vapour == 'ideal gas' && m.poynting && m.equation == ''
	assert m.pairs.map('${it.a} ${it.b}') == ['water 1-butanol', 'water methanol', '1-butanol methanol']
	assert m.pairs.all(it.source.starts_with('nrtl.ipd line '))
	assert m.pairs[1].parameters.map(it.name) == ['tau_ab', 'tau_ba', 'alpha']
	assert m.pairs[1].parameters.map(it.text) == ['792.8020 K / (1.98721*T)',
		'-189.0469 K / (1.98721*T)', '0.2999']
	assert m.warnings.len == 0
	assert rep.table().contains('The liquid by NRTL, the vapour the ideal gas, with the Poynting correction.')
}

// A port with no phase gets an empty stream: at 300 K the feed is two
// liquids, so no vapour; and nearly pure water is one liquid, D9's liquid,
// out of light (NRTL's water holds 0.0070 butanol at 300 K).
fn test_a_decanter_with_fewer_phases() {
	cold := with('outlet_T = "364.2 K"', 'outlet_T = "300 K"')
	sol := load(cold)!.solve()!
	v := sol.streams['vapour'] or { panic('no vapour') }
	assert v.total() == 0.0
	for name in ['organic', 'aqueous'] {
		s := sol.streams[name] or { panic('no ${name}') }
		assert s.total() > 1.0, name
	}
	wet := with('mole_fractions = { water = 0.75, 1-butanol = 0.2, methanol = 0.05 }', 'mole_fractions = { water = 0.995, 1-butanol = 0.004, methanol = 0.001 }')
	dry := load(wet.replace('outlet_T = "364.2 K"', 'outlet_T = "300 K"'))!.solve()!
	light := dry.streams['organic'] or { panic('no organic') }
	heavy := dry.streams['aqueous'] or { panic('no aqueous') }
	hot := dry.streams['hot'] or { panic('no hot') }
	// The inlet's own flows, not a phase's share of them.
	assert light.flows == hot.flows && heavy.total() == 0.0
	assert light.result.phases.map(it.label) == [thermo.PhaseLabel.liquid]
}

// A flash drum sends both liquids out of its liquid port, and says so
// (D11).
fn test_a_drum_with_two_liquids_says_so() {
	text := with('type = "decanter"', 'type = "flash_drum"').replace('light = "organic"\nheavy = "aqueous"',
		'liquid = "liquids"').replace('outlet_T = "364.2 K"', 'outlet_T = "300 K"')
	l := load(text)!
	sol := l.solve()!
	liquid := sol.streams['liquids'] or { panic('no liquids') }
	assert liquid.result.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2]
	assert math.abs(liquid.total() - 10.0) <= 1e-12
	rep := l.report(sol)
	drum := rep.units.filter(it.name == 'decanter')[0]
	assert drum.warnings == ['its liquid outlet has two liquid phases, liquid and liquid2, which a decanter would separate']
}

// A file's pair is used in place of the bank's: the bank's own numbers as
// the file's give the same bits, either way round; another value moves the
// answer; and the report says the file gave it.
fn test_a_files_pair_replaces_the_banks() {
	plain := load(decanter_text())!.solve()!
	same := '[model]\nliquid = "NRTL"\n\n[[model.pairs]]\na = "water"\nb = "methanol"\ntau_ab = "792.8020 K / (1.98721*T)"\ntau_ba = "-189.0469 K / (1.98721*T)"\nalpha = "0.2999"\n'
	turned := '[model]\nliquid = "NRTL"\n\n[[model.pairs]]\na = "methanol"\nb = "water"\ntau_ab = "-189.0469 K / (1.98721*T)"\ntau_ba = "792.8020 K / (1.98721*T)"\nalpha = "0.2999"\n'
	for text in [same, turned] {
		l := load(with('[model]\nliquid = "NRTL"\n', text))!
		sol := l.solve()!
		for name in ['vapour', 'organic', 'aqueous'] {
			a := plain.streams[name] or { panic(name) }
			b := sol.streams[name] or { panic(name) }
			assert a.flows == b.flows && a.t == b.t && a.h == b.h, name
		}
		m := l.report(sol).model or { panic('no model') }
		assert m.pairs[1].source == 'the file'
		assert m.pairs[1].parameters.map(it.text) == ['792.8020 K / (1.98721*T)',
			'-189.0469 K / (1.98721*T)', '0.2999']
		assert m.pairs[0].source.starts_with('nrtl.ipd line ')
	}
	// In SI, with R: the same τ to 5e-7 (DECHEMA's R against the SI's over
	// 4.184 J), a different answer.
	si := '[model]\nliquid = "NRTL"\n\n[[model.pairs]]\na = "water"\nb = "methanol"\ntau_ab = "3317.08 J/mol / (R*T)"\ntau_ba = "-790.97 J/mol / (R*T)"\nalpha = "0.3"\n'
	moved := load(with('[model]\nliquid = "NRTL"\n', si))!.solve()!
	a := plain.streams['vapour'] or { panic('vapour') }
	b := moved.streams['vapour'] or { panic('vapour') }
	assert a.flows != b.flows
	assert math.abs(a.total() - b.total()) / a.total() < 0.05
}

// A pair the bank lacks is refused by default, naming it; missing = "ideal"
// makes it ideal, and the report warns of it.
fn test_a_missing_pair_is_refused_or_ideal() {
	hexane := with('names = ["water", "1-butanol", "methanol"]', 'names = ["water", "1-butanol", "n-hexane"]').replace('methanol = 0.05',
		'n-hexane = 0.05')
	refused(hexane, lined('liquid ='), 'model.liquid', "nrtl.ipd doesn't list Water with N-hexane",
		'')
	ideal := hexane.replace('liquid = "NRTL"', 'liquid = "NRTL"\nmissing = "ideal"')
	l := load(ideal)!
	m := l.model
	pair := m.pairs.filter(it.a == 'water' && it.b == 'n-hexane')[0]
	assert pair.source == 'none listed: an ideal pair'
	assert pair.parameters.map(it.text) == ['0', '0', '0.3']
	assert m.warnings == [
		'water with n-hexane: nrtl.ipd has no pair, so an ideal pair (missing = "ideal")',
		'1-butanol with n-hexane: nrtl.ipd has no pair, so an ideal pair (missing = "ideal")',
	]
	uniquac := ideal.replace('liquid = "NRTL"', 'liquid = "UNIQUAC"')
	mu := load(uniquac)!.model
	assert mu.pairs.filter(it.a == 'water' && it.b == 'n-hexane')[0].source == 'none listed: no residual interaction'
	sol := l.solve()!
	assert sol.streams.len == 5
}

// UNIFAC, Wilson and UNIQUAC files load, each with its own keys; a cubic
// vapour takes the bank's kᵢⱼ, or none; and an equation of state's report
// gives each kᵢⱼ with its source.
fn test_each_model_and_vapour_loads() {
	for liquid in ['UNIFAC', 'UNIFAC-LLE', 'UNIFAC-Dortmund', 'UNIQUAC', 'Wilson'] {
		l := load(with('liquid = "NRTL"', 'liquid = "${liquid}"'))!
		assert l.model.liquid == liquid
		assert l.model.pairs.len == if liquid.starts_with('UNIFAC') { 0 } else { 3 }, liquid
		l.solve()!
	}
	wilson := '[model]\nliquid = "Wilson"\n\n[[model.pairs]]\na = "water"\nb = "methanol"\nlambda_ab = "0.5"\nlambda_ba = "0.9"\n'
	w := load(with('[model]\nliquid = "NRTL"\n', wilson))!
	assert w.model.pairs[1].parameters.map(it.name) == ['lambda_ab', 'lambda_ba']
	pr := load(with('liquid = "NRTL"', 'liquid = "NRTL"\nvapour = "PR"\nkij = false\npoynting = false'))!
	assert pr.model.vapour == 'PR' && !pr.model.poynting
	pr.solve()!
	// The vapour's kᵢⱼ: none with kij = false, the bank's by default.
	for kij, any in {
		'false': false
		'true':  true
	} {
		v := load(with('liquid = "NRTL"', 'liquid = "NRTL"\nvapour = "PR"\nkij = ${kij}'))!
		if v.fluid.model is thermo.GammaPhi {
			g := v.fluid.model as thermo.GammaPhi
			c := g.vapour or { panic('no vapour') }
			if c is thermo.Cubic {
				assert c.kij.any(it != 0.0) == any, kij
			} else {
				assert false, 'not a cubic'
			}
		} else {
			assert false, 'not γ–φ'
		}
	}
	cubic := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml'))!
	eos := load(cubic)!.model
	assert eos.equation == 'PR' && eos.liquid == '' && eos.pairs.len == 10
	assert eos.pairs[0].parameters[0].name == 'kij'
	assert eos.pairs.all(it.source.starts_with('pr.ipd line ') || it.source == 'none listed: 0')
	assert eos.notes.len == 0
	// No bank: every kᵢⱼ none listed, 0; and an α not the equation's own,
	// each a note.
	unlisted := load(cubic.replace('equation = "PR"', 'equation = "PR"\nkij = false'))!.model
	assert unlisted.pairs.all(it.source == 'none listed: 0' && it.parameters[0].text == '0.0')
	mc := load(cubic.replace('equation = "PR"', 'equation = "PR"\nalpha = "mathias_copeman"'))!.model
	assert mc.notes.len == 5 && mc.notes[0].starts_with('methane: α Mathias–Copeman, mc-PR.lib line '), '${mc.notes}'
}

// Every key is checked, at its line (criterion 10).
fn test_a_wrong_model_is_refused_at_its_line() {
	model := lined('liquid =')
	refused(with('format = 4', 'format = 3'), model, 'model.liquid', 'an activity model needs format 4',
		'')
	refused(with('liquid = "NRTL"', 'liquid = "NTRL"'), model, 'model.liquid', 'no activity model NTRL',
		'NRTL')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nequation = "PR"'), model + 1, 'model.equation',
		'not equation', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nalpha = "equation"'), model + 1, 'model.alpha',
		"an equation of state's", '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nvapour = "Pr"'), model + 1, 'model.vapour',
		'no vapour Pr', 'PR')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nkij = true'), model + 1, 'model.kij', 'the vapour is the ideal gas',
		'')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\npoynting = "yes"'), model + 1, 'model.poynting',
		'true or false', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nhenry = ["hydrogen"]'), model + 1, 'model.henry',
		'hydrogen is not one of the components', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"\nmissing = "skip"'), model + 1, 'model.missing',
		'"refuse" or "ideal"', '')
	refused(with('liquid = "NRTL"', 'liquid = "UNIFAC"\npairs = [{ a = "water", b = "methanol", tau_ab = "1", tau_ba = "1" }]'),
		model + 1, 'model.pairs', "UNIFAC's parameters are its groups'", '')
	refused(with('liquid = "NRTL"', 'colour = "blue"'), model - 1, 'model', 'a model needs equation, for an equation of state, or liquid',
		'')
	refused(with('liquid = "NRTL"', 'equation = "PR"\npoynting = true'), model + 1, 'model.poynting',
		"a γ–φ fluid's", '')
	// A pair.
	pair := '\n\n[[model.pairs]]\na = "water"\nb = "methanol"\ntau_ab = "1"\ntau_ba = "2"\nalpha = "0.3"'
	at := model + 2 // the pair's [[model.pairs]] line
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair.replace('"methanol"', '"methanl"')),
		at + 2, 'model.pairs[0].b', 'methanl is not one of the components', 'methanol')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair.replace('"methanol"', '"water"')),
		at + 2, 'model.pairs[0].b', 'water twice', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair.replace('tau_ab = "1"', 'tau_ab = "100 K"')),
		at + 3, 'model.pairs[0].tau_ab', 'τ', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair.replace('tau_ab = "1"', 'tau_ab = "1 +"')),
		at + 3, 'model.pairs[0].tau_ab', '', '')
	// A key missing from an inline table is at its first key, the table's
	// own line in TOML's array of tables being only its header.
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair.replace('\nalpha = "0.3"', '')),
		at + 1, 'model.pairs[0]', 'alpha is missing', '')
	refused(with('liquid = "NRTL"', 'liquid = "UNIQUAC"' + pair), at + 5, 'model.pairs[0].alpha',
		'no key alpha', '')
	refused(with('liquid = "NRTL"', 'liquid = "NRTL"' + pair + pair.replace('a = "water"\nb = "methanol"',
		'a = "methanol"\nb = "water"')), at + 1, 'model.pairs', 'given twice', '')
	// The decanter needs format 4, and its ports.
	unit := lined('type = "decanter"')
	refused(with('format = 4', 'format = 3').replace('liquid = "NRTL"', 'equation = "PR"'), unit,
		'unit.decanter.type', 'a decanter needs format 4', '')
	refused(with('heavy = "aqueous"\n', ''), unit - 1, 'unit.decanter', 'out is missing', '')
}

// A Henry's-law solute the file names is one, and a pair can't name it:
// the activity model is the solvents'. And the database refuses pairs for
// UNIFAC itself, as the file does.
fn test_henrys_law_solutes_and_pairs() {
	co2 := with('names = ["water", "1-butanol", "methanol"]', 'names = ["water", "methanol", "carbon dioxide"]').replace('mole_fractions = { water = 0.75, 1-butanol = 0.2, methanol = 0.05 }',
		'mole_fractions = { water = 0.75, methanol = 0.2, "carbon dioxide" = 0.05 }')
	plain := load(co2.replace('liquid = "NRTL"', 'liquid = "NRTL"
henry = ["carbon dioxide"]'))!
	assert plain.model.solutes == ['carbon dioxide']
	assert plain.model.pairs.map('${it.a} ${it.b}') == ['water methanol']
	pair := '

[[model.pairs]]
a = "water"
b = "carbon dioxide"
tau_ab = "1"
tau_ba = "1"
alpha = "0.3"'
	refused(co2.replace('liquid = "NRTL"', 'liquid = "NRTL"
henry = ["carbon dioxide"]' + pair),
		lined('liquid =') + 4, 'model.pairs', "a Henry's-law solute", '')
	cas := ['7732-18-5', '67-56-1']
	if _ := db.gamma_phi(cas,
		liquid: .unifac
		pairs:  [components.PairOverride{cas[0], cas[1], ['1', '1', '0.3']}]
	) {
		assert false, 'UNIFAC took pairs'
	} else {
		assert err.msg().contains("UNIFAC's parameters are its groups'"), err.msg()
	}
	// A pair of the wrong count, which the file can't give.
	if _ := db.gamma_phi(cas, pairs: [components.PairOverride{cas[0], cas[1], ['1', '1']}]) {
		assert false, 'NRTL took two parameters'
	} else {
		assert err.msg().contains('has 2 parameters, where the model takes 3'), err.msg()
	}
}

// Load, write and load again: the same file, and the same numbers.
fn test_a_format_4_file_writes_back() {
	pair := '[model]\nliquid = "NRTL"\nvapour = "PR"\n\n[[model.pairs]]\na = "water"\nb = "methanol"\ntau_ab = "3317.08 J/mol / (R*T)"\ntau_ba = "-790.97 J/mol / (R*T)"\nalpha = "0.3"\n'
	l := load(with('[model]\nliquid = "NRTL"\n', pair))!
	again := load(l.write())!
	assert again.write() == l.write()
	a := l.solve()!
	b := again.solve()!
	for name, s in a.streams {
		t := b.streams[name] or { panic(name) }
		assert s.flows == t.flows && s.h == t.h, name
	}
}
