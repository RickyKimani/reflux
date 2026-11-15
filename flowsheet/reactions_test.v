// A file's reactions (M9 step 6, D11): the three kinetics files in
// testdata/kinetics/ read into the same rates as the papers' and DWSIM's,
// from testdata/m9_kinetics.json, or as kinetics' own structured forms;
// their checks and factor tables; and every way a file is refused, at its
// line.
module flowsheet

import components
import json2
import kinetics
import math
import os

const kinetics_dir = os.join_path(@VMODROOT, 'testdata', 'kinetics')

fn kfile(name string) KineticsFile {
	return load_kinetics(os.join_path(kinetics_dir, name)) or { panic(err) }
}

struct KCase {
	t      f64
	p      []f64
	rates  []f64
	scales []f64
}

struct KModel {
	name    string
	species []string
	cases   []KCase
}

struct KFixture {
	models []KModel
}

fn kfixture() KFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_kinetics.json')
	return json2.decode[KFixture](os.read_file(path) or { panic(err) }) or { panic(err) }
}

// point_at is a fixture's state as a point: T, and the partial pressures
// exactly, in the fixture's species order, which is the file's; P their
// sum.
fn point_at(k &kinetics.Kinetics, c KCase) kinetics.Point {
	mut pt := k.point()
	pt.x[0] = c.t
	mut total := 0.0
	for i, p in c.p {
		pt.x[k.ctx.slot(0, i)] = p
		total += p
	}
	pt.x[1] = total
	return pt
}

// A file's rates are the fixture's, to 1e-13 of each rate's largest term,
// as kinetics' own structured forms are (kinetics/law_test.v).
fn rates_match(k &kinetics.Kinetics, m KModel) {
	mut buf := k.buffer()
	mut out := []f64{len: k.rates.len}
	for c in m.cases {
		pt := point_at(k, c)
		k.rates(&pt, mut buf, mut out)
		for i, r in out {
			assert math.abs(r - c.rates[i]) <= 1e-13 * c.scales[i], '${m.name} rate ${i} at ${c.t} K: ${r} vs ${c.rates[i]}'
		}
	}
}

fn test_vanden_bussche_froment_from_its_file() {
	f := kfile('vanden_bussche_froment.toml')
	assert f.names == ['carbon dioxide', 'hydrogen', 'methanol', 'water', 'carbon monoxide']
	rates_match(&f.kinetics, kfixture().models[1])
}

// The owner's DWSIM case's rates, its three denominators and typed 8.314
// included, are DWSIM's own expressions'.
fn test_dwsim_from_its_file() {
	f := kfile('dwsim_reforming.toml')
	fx := kfixture()
	assert fx.models[2].name == 'DWSIM'
	rates_match(&f.kinetics, fx.models[2])
}

// Xu and Froment's file, K from the Gibbs energies, gives the bits of the
// same kinetics built from kinetics' structured forms: the file is read
// into them key for key.
fn test_xu_froment_from_its_file() {
	f := kfile('xu_froment.toml')
	k := f.kinetics
	db := components.load_default()!
	sp := ['methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide'].map(species_of(db.find(it) or {
		panic(err)
	}))
	rs := [
		kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!,
		kinetics.parse_equation('shift', 'CO + H2O = CO2 + H2', sp)!,
		kinetics.parse_equation('overall', 'CH4 + 2 H2O = CO2 + 4 H2', sp)!,
	]
	site := kinetics.Site{
		name:   'Ni'
		symbol: 'DEN'
		terms:  [
			kinetics.Monomial{'K_CO', kinetics.VantHoff{'8.23e-5 bar^-1', '-70.65 kJ/mol', ''}, [
				kinetics.Power{'CO', 1},
			]},
			kinetics.Monomial{'K_H2', kinetics.VantHoff{'6.12e-9 bar^-1', '-82.90 kJ/mol', ''}, [
				kinetics.Power{'H2', 1},
			]},
			kinetics.Monomial{'K_CH4', kinetics.VantHoff{'6.65e-4 bar^-1', '-38.28 kJ/mol', ''}, [
				kinetics.Power{'CH4', 1},
			]},
			kinetics.Monomial{'K_H2O', kinetics.VantHoff{'1.77e5', '88.68 kJ/mol', ''}, [
				kinetics.Power{'H2O', 1},
				kinetics.Power{'H2', -1},
			]},
		]
	}
	ks := [
		kinetics.Arrhenius{'4.225e15 kmol*bar^0.5/(kg*h)', '240.1 kJ/mol', ''},
		kinetics.Arrhenius{'1.955e6 kmol/(kg*h*bar)', '67.13 kJ/mol', ''},
		kinetics.Arrhenius{'1.020e15 kmol*bar^0.5/(kg*h)', '243.9 kJ/mol', ''},
	]
	orders := [
		[kinetics.Power{'CH4', 1}, kinetics.Power{'H2O', 1}, kinetics.Power{'H2', -2.5}],
		[kinetics.Power{'CO', 1}, kinetics.Power{'H2O', 1}, kinetics.Power{'H2', -1}],
		[kinetics.Power{'CH4', 1}, kinetics.Power{'H2O', 2}, kinetics.Power{'H2', -3.5}],
	]
	mut laws := []kinetics.Rated{}
	for i in 0 .. 3 {
		laws << kinetics.Rated{i, kinetics.Lhhw{
			k:           ks[i]
			orders:      orders[i]
			equilibrium: kinetics.Equilibrium{
				kind: .thermo
			}
			form:        .difference
			site:        'Ni'
			power:       2
		}}
	}
	want := kinetics.build(sp, rs, [site], laws, [])!
	assert k.rates.map(it.plain()) == want.rates.map(it.plain())
	s := kinetics.State{
		t: 900
		p: 15e5
		y: [0.2, 0.45, 0.2, 0.1, 0.05]
	}
	mut p1 := k.point()
	mut p2 := want.point()
	k.fill(s, mut p1)!
	want.fill(s, mut p2)!
	mut b1 := k.buffer()
	mut b2 := want.buffer()
	mut o1 := []f64{len: 3}
	mut o2 := []f64{len: 3}
	k.rates(&p1, mut b1, mut o1)
	want.rates(&p2, mut b2, mut o2)
	assert o1.map(math.f64_bits(it)) == o2.map(math.f64_bits(it))
}

fn messages(k kinetics.Kinetics, s ?kinetics.State) []string {
	return k.check(s).map(it.str())
}

const reforming_state = kinetics.State{
	t: 900
	p: 15e5
	y: [0.2, 0.45, 0.2, 0.1, 0.05]
}

// What `reflux kinetics check` finds in each file.
fn test_the_checks() {
	xf := messages(kfile('xu_froment.toml').kinetics, none)
	assert xf == [
		'note: overall = reforming + shift: fine for kinetics, each having its own rate, but an equilibrium over these reactions needs them independent (2 of 3 are)',
		'warning: K_H2O rises with T: its ΔH, 88.68 kJ/mol, is above zero, where adsorption is exothermic (unless the term lumps more than an adsorption)',
	]
	vbf := messages(kfile('vanden_bussche_froment.toml').kinetics, kinetics.State{
		t: 513.15
		p: 50e5
		y: [0.08, 0.65, 0.04, 0.03, 0.08]
	})
	assert vbf.len == 3
	assert vbf[0].starts_with("warning: methanol's rate constant has a negative activation energy, -36.696 kJ/mol")
	assert vbf[1].starts_with("note: methanol's fitted K, K1, is a factor of 0.76987 of the Gibbs energies' at 513.15 K")
	assert vbf[2].starts_with("note: RWGS's fitted K, K2, is a factor of 1.0452 of the Gibbs energies' at 513.15 K")
	dw := messages(kfile('dwsim_reforming.toml').kinetics, reforming_state)
	assert dw.len == 13
	assert dw.filter(it.contains('has its own adsorption term')) == [
		"note: shift has its own adsorption term in place of site Ni's, without K_CH4: a site's term belongs to every reaction on it",
		"note: overall has its own adsorption term in place of site Ni's, without K_CO: a site's term belongs to every reaction on it",
	]
	assert dw.filter(it.starts_with("warning: overall's fitted K")) == [
		"warning: overall's fitted K, K3, is a factor of 0.20173 of the Gibbs energies' at 900 K (ln K -1.6008); from 300 to 1500 K its ln K is -24.687 to +2.0051 from theirs",
	]
	assert dw.filter(it.contains('8.314')).len == 7
	assert dw.contains("warning: k1 has 8.314 J/(mol*K), 5.564e-5 below R, the molar gas constant, N_A·k_B: write R unless the paper fitted with this number; at 900 K it moves reforming's rate by -0.17838% from what R gives")
}

// The factor table: each part's value at the state, its site's term the
// sum of its terms, and the rate the tape's.
fn test_the_factor_table() {
	k := kfile('xu_froment.toml').kinetics
	ex := k.explain(reforming_state)!
	assert ex.len == 3
	e := ex[0]
	assert e.law == "LHHW over site Ni's adsorption term squared on partial pressures, per mass of catalyst, reversible, K from the Gibbs energies"
	assert e.unicode == 'k₁/p_H₂^2.5·(p_CH₄·p_H₂O − p_CO·p_H₂³/(Keq_reforming·p°²))/DEN²'
	assert e.uses.map(it.label) == ['K_CO', 'K_H2', 'K_CH4', 'K_H2O', 'DEN', 'k1']
	labels := e.parts.map(it.label)
	assert labels[0] == 'its rate constant'
	assert labels.filter(it == 'a term of it').len == 4
	den := e.parts.filter(it.label == "its site's adsorption term")[0].value
	terms := e.parts.filter(it.label == 'a term of it').map(it.value)
	assert math.abs(den - (1 + terms[0] + terms[1] + terms[2] + terms[3])) <= 1e-15 * den
	q := e.parts.filter(it.label.starts_with('Q, '))[0].value
	p := 15e5
	assert math.abs(q - 0.1 * p * math.pow(0.2 * p, 3) / (0.2 * p * (0.45 * p))) <= 1e-14 * q
	assert e.parts.filter(it.label.starts_with('its K'))[0].unit == 'Pa^2'
	assert e.parts[0].unit == 'mol·Pa^0.5/(kg·s)'
	mut pt := k.point()
	k.fill(reforming_state, mut pt)!
	mut buf := k.buffer()
	mut out := []f64{len: 3}
	k.rates(&pt, mut buf, mut out)
	assert ex.map(it.rate) == out
	text := e.str()
	assert text.starts_with('reforming: CH4 + H2O = CO + 3 H2\n  LHHW over')
	assert text.ends_with('  r = ${kinetics.sig(out[0])} mol/(kg·s)')
}

fn refused_text(text string) string {
	f := load_kinetics_text(text) or { return err.msg() }
	return 'loaded ${f.kinetics.rates.len} rates'
}

const xu_froment_text = os.read_file(os.join_path(kinetics_dir, 'xu_froment.toml')) or {
	panic(err)
}

fn with(old string, new string) string {
	assert xu_froment_text.contains(old), old
	return xu_froment_text.replace_once(old, new)
}

// Every way a file is refused, with its line and key.
fn test_refused() {
	assert refused_text(xu_froment_text) == 'loaded 3 rates'
	cases := [
		['format = 2', 'format = 1', 'flowsheet: line 6: format: reactions need format 2, not 1'],
		['equation = "CH4 + H2O = CO + 3 H2"', 'equation = "CH4 + H2O = CO + H2"',
			"flowsheet: line 22: reaction.reforming.equation: CH4 + H2O = CO + H2 doesn't balance: H: 6 on the left, 2 on the right; balanced, it would be CH4 + H2O = CO + 3 H2"],
		['equation = "CH4 + H2O = CO + 3 H2"', 'equation = "CH4 + H2O = CO + 3 H3"',
			'flowsheet: line 22: reaction.reforming.equation: column 18: no species here is H3 (did you mean H2, CH4, H2O?)'],
		['law = "lhhw"', 'law = "lhw"',
			'flowsheet: line 23: reaction.reforming.law: no law lhw: one of power_law, lhhw, eley_rideal, expression (did you mean lhhw?)'],
		['K = "gibbs"', 'Kx = "gibbs"', 'flowsheet: line 21: reaction.reforming: K is missing'],
		['driving = "difference"', 'driving = "different"',
			'flowsheet: line 27: reaction.reforming.driving: no driving different: one of approach, difference (did you mean difference?)'],
		['site = "Ni"', 'site = "Pt"', 'flowsheet: line 21: reaction.reforming: reforming: no site Pt'],
		['Ea = "240.1 kJ/mol"', 'Ea = "240.1 kJ"',
			"flowsheet: line 21: reaction.reforming: definition k1: reforming's rate constant: k1: line 1, column 32: the argument of exp must be a pure number, but it is an amount of substance, mol"],
		['orders = { CH4 = 1, H2O = 1, H2 = -2.5 }', 'orders = { CH4 = "1" }',
			"flowsheet: line 25: reaction.reforming.orders: CH4's exponent should be a number"],
		['[reaction.reforming]', '[reaction.reforming]\nspeed = 3',
			'flowsheet: line 22: reaction.reforming.speed: no key speed here'],
		['law = "lhhw"\nk = { A = "4.225e15', 'k = { A = "4.225e15',
			'flowsheet: line 23: reaction.reforming.k: a reaction without a law has no k: give its law'],
		['a = { CO = 1 } },', 'a = {} },',
			'flowsheet: line 15: site.Ni.terms[0].a: a term needs at least one species'],
		['symbol = "DEN"', 'symbol = "DEN"\nwhat = 1',
			'flowsheet: line 14: site.Ni.what: no key what here'],
		['[site.Ni]', '[definitions]\nDEN = "2"\n\n[site.Ni]',
			"flowsheet: line 15: site.Ni: site Ni's adsorption term: DEN already names a definition of your own: give this one a name of its own"],
	]
	for c in cases {
		assert refused_text(with(c[0], c[1])) == c[2], c[1]
	}
	assert refused_text('format = 2\n[components]\nnames = ["methane"]\n[reaction]\n') == 'flowsheet: line 4: reaction: a file needs at least one reaction'
	// A whole flowsheet file's other sections are left to the flowsheet's
	// reader.
	assert refused_text(with('[site.Ni]', '[model]\nequation = "PR"\n\n[unit.r]\ntype = "pfr"\n\n[site.Ni]')) == 'loaded 3 rates'
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

// 500 seeded files, each Xu and Froment's with a value or two replaced by
// another of any kind: each loads, or is refused with a FlowsheetError
// that has its reason, and none panics.
fn test_random_values_never_panic() {
	values := ['"lhhw"', '"power_law"', '"eley_rideal"', '"expression"', '"gibbs"', '"none"', '1',
		'-2.5', '0', '"CH4"', '"H5"', '"3 bar"', '"240.1 kJ/mol"', '"exp(1/T)"', '"p[CH4]"', '{}',
		'[]', '{ A = "1" }', '{ A = "1 bar", dH = "1 J/mol" }', '{ CO = 1 }', '{ CO = "x" }', 'true',
		'"Ni"', '"Pt"', '"difference"', '"x = y"', '"1 mol/(kg*s)"', '[{ K = "1", a = { CO = 1 } }]',
		'"ü"', '""']
	lines := xu_froment_text.split_into_lines()
	mut rng := Rng{20261110}
	mut loaded := 0
	db := components.load_default()!
	for _ in 0 .. 500 {
		mut out := lines.clone()
		for _ in 0 .. 1 + rng.below(2) {
			i := rng.below(out.len)
			if out[i].contains(' = ') && !out[i].starts_with('#') && !out[i].starts_with(' ') {
				out[i] = out[i].all_before(' = ') + ' = ' + values[rng.below(values.len)]
			}
		}
		text := out.join('\n')
		f := load_kinetics_text(text, db: db) or {
			assert err is FlowsheetError, '${err.msg()}\n${text}'
			e := err as FlowsheetError
			assert e.reason != '', text
			continue
		}
		assert f.kinetics.rates.len <= 3
		loaded++
	}
	assert loaded > 50
}

// Every law, basis and kind of constant read from a file is the same
// kinetics as built from kinetics' own structures: the same rates and
// definitions, written the same, and the same law said by explain.
fn test_every_law_and_basis_from_a_file() {
	text := 'format = 2

[components]
names = ["methane", "water", "hydrogen", "carbon monoxide", "carbon dioxide"]

[site.S]
terms = [{ K = { K_ref = "2 bar^-1", dH = "-50 kJ/mol", T_ref = "800 K" }, a = { CO = 1 } }]

[reaction.a]
equation = "CH4 + H2O = CO + 3 H2"
law = "power_law"
k = { k_ref = "1e-3 m^3/(mol*s)", Ea = "100 kJ/mol", T_ref = "900 K" }
orders = { CH4 = 1, H2O = 1 }
K = "gibbs"
basis = "concentration"
per = "volume"

[reaction.b]
equation = "CO + H2O = CO2 + H2"
law = "eley_rideal"
k = "1e-9 mol/(kg*s*Pa)"
adsorbed = "CO"
gas = "H2O"
K = "none"
site = "S"
basis = "fugacity"

[reaction.c]
equation = "CH4 + 2 H2O = CO2 + 4 H2"
law = "expression"
rate = "1e-3 mol/(kg*s) * y[CH4]"
basis = "mole_fraction"
'
	f := load_kinetics_text(text)!
	db := components.load_default()!
	sp := ['methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide'].map(species_of(db.find(it) or {
		panic(err)
	}))
	rs := [
		kinetics.parse_equation('a', 'CH4 + H2O = CO + 3 H2', sp)!,
		kinetics.parse_equation('b', 'CO + H2O = CO2 + H2', sp)!,
		kinetics.parse_equation('c', 'CH4 + 2 H2O = CO2 + 4 H2', sp)!,
	]
	site := kinetics.Site{
		name:  'S'
		terms: [
			kinetics.Monomial{'', kinetics.VantHoff{'2 bar^-1', '-50 kJ/mol', '800 K'}, [
				kinetics.Power{'CO', 1},
			]},
		]
	}
	laws := [
		kinetics.Rated{0, kinetics.PowerLaw{
			k:           kinetics.Arrhenius{'1e-3 m^3/(mol*s)', '100 kJ/mol', '900 K'}
			orders:      [kinetics.Power{'CH4', 1}, kinetics.Power{'H2O', 1}]
			equilibrium: kinetics.Equilibrium{
				kind: .thermo
			}
			basis:       .concentration
			per:         .volume
		}},
		kinetics.Rated{1, kinetics.EleyRideal{
			k:        kinetics.Written{'1e-9 mol/(kg*s*Pa)'}
			adsorbed: 'CO'
			gas:      'H2O'
			site:     'S'
			basis:    .fugacity
		}},
		kinetics.Rated{2, kinetics.Expression{
			text:  '1e-3 mol/(kg*s) * y[CH4]'
			basis: .fraction
		}},
	]
	want := kinetics.build(sp, rs, [site], laws, [])!
	k := f.kinetics
	assert k.rates.map(it.plain()) == want.rates.map(it.plain())
	assert k.defs.names == want.defs.names
	assert k.defs.exprs.map(it.plain()) == want.defs.exprs.map(it.plain())
	s := kinetics.State{
		t: 900
		p: 2e5
		y: [0.2, 0.4, 0.2, 0.1, 0.1]
	}
	got := k.explain(s)!.map(it.law)
	assert got == want.explain(s)!.map(it.law)
	assert got == [
		'a power law on concentrations, per volume, reversible, K from the Gibbs energies',
		'Eley–Rideal on site S on fugacities, per mass of catalyst, irreversible',
		'an expression on mole fractions, per mass of catalyst',
	]
	// A format other than 2 is refused, whichever way it's off.
	for n in ['1', '3'] {
		g := load_kinetics_text(text.replace_once('format = 2', 'format = ${n}')) or {
			assert err.msg() == 'flowsheet: line 1: format: reactions need format 2, not ${n}'
			continue
		}
		assert false, '${g.names}'
	}
}
