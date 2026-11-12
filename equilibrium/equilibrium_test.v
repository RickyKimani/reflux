// Chemical equilibrium (M9 step 7, criterion 5): Gibbs minimisation and the
// equilibrium over chosen reactions, ideal gas and Peng–Robinson, against
// mpmath's at 50 digits by another method (testdata/m9_equilibrium.json);
// and everything refused.
module equilibrium

import components
import json2
import kinetics
import math
import os

const db = components.load_default() or { panic(err) }

fn species(cas []string) []kinetics.Species {
	mut out := []kinetics.Species{}
	for c in cas {
		r := db.find(c) or { panic(err) }
		out << kinetics.Species{
			name:    r.name
			cas:     r.ids.cas
			formula: r.formula() or { panic(err) }
			gas:     r.ideal_gas() or { panic(err) }
			dgf:     r.constants.dgf
		}
	}
	return out
}

struct ECase {
	name      string
	species   []string
	model     string
	reactions []string
	inert     []int
	t         f64
	p         f64
	n0        []f64
	n         []f64
	ln_y      []?f64
}

struct EFixture {
	cases []ECase
}

fn fixture() EFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_equilibrium.json')
	return json2.decode[EFixture](os.read_file(path) or { panic(err) }) or { panic(err) }
}

// mixture is a case's species, with Peng–Robinson (every kᵢⱼ zero) where
// the case has it.
fn mixture(c ECase) Mixture {
	sp := species(c.species)
	if c.model == 'pr' {
		names := c.species.map(db.find(it) or { panic(err) }.name)
		fl := db.fluid(.pr, names, kij: false) or { panic(err) }
		return Mixture{
			species: sp
			model:   fl.model
		}
	}
	return Mixture{
		species: sp
	}
}

fn run(c ECase) !Result {
	m := mixture(c)
	if c.reactions.len == 0 {
		return gibbs(m, c.t, c.p, c.n0, c.inert)
	}
	mut rs := []kinetics.Reaction{}
	for i, e in c.reactions {
		rs << kinetics.parse_equation('r${i + 1}', e, m.species)!
	}
	return over_reactions(m, rs, ln_k_gibbs(rs, m.species, c.t)!, c.t, c.p, c.n0)
}

// Every case's ln y is mpmath's to 1e-11: μ°/RT comes from Cp integrals good
// to 1e-13 of themselves, in values of order 100 (1.6e-12 at worst). A
// species the case has none of has none in reflux's either. The RAND steps
// and the φ updates are recorded, and can't grow.
fn test_against_mpmath() {
	f := fixture()
	assert f.cases.len == 21
	mut worst := 0.0
	mut most_steps, mut most_outer := 0, 0
	mut longest, mut most_updated := '', ''
	for c in f.cases {
		r := run(c) or { panic('${c.name}: ${err}') }
		for j, want in c.ln_y {
			if w := want {
				d := math.abs(math.log(r.y[j]) - w)
				worst = math.max(worst, d)
				assert d <= 1e-11, '${c.name}, species ${j}: ln y ${math.log(r.y[j])} vs ${w}'
			} else {
				assert r.n[j] == 0, '${c.name}, species ${j}: ${r.n[j]}'
			}
		}
		if r.steps > most_steps {
			most_steps, longest = r.steps, c.name
		}
		if r.outer > most_outer {
			most_outer, most_updated = r.outer, c.name
		}
	}
	println('equilibria: ln y within ${worst}; at most ${most_steps} RAND steps (${longest}) and ${most_outer} φ updates (${most_updated})')
	assert most_steps <= 36
	assert most_outer <= 13
}

// Held inert, CO₂ leaves reforming alone among the other four: the same
// answer as the equilibrium over reforming by itself.
fn test_inert_is_one_reaction_fewer() {
	sp := species(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	m := Mixture{
		species: sp
	}
	n0 := [1.0, 2, 0.1, 0, 0.5]
	a := gibbs(m, 900, 1e5, n0, [4])!
	sr := kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!
	b := over_reactions(m, [sr], ln_k_gibbs([sr], sp, 900)!, 900, 1e5, n0)!
	assert a.n[4] == 0.5 && b.n[4] == 0.5
	for j in 0 .. 5 {
		assert math.abs(math.log(a.y[j]) - math.log(b.y[j])) <= 1e-12
	}
}

// From a feed without carbon, the equilibrium over reforming and the shift
// makes none of carbon's compounds: the elements say so, where what the
// reactions conserve, a basis that mixes the elements, doesn't. The
// consistency check found this, as a singular matrix.
fn test_an_element_the_feed_has_not() {
	sp := species(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	m := Mixture{
		species: sp
	}
	rs := [
		kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!,
		kinetics.parse_equation('shift', 'CO + H2O = CO2 + H2', sp)!,
	]
	n0 := [0.0, 2, 1, 0, 0]
	r := over_reactions(m, rs, ln_k_gibbs(rs, sp, 900)!, 900, 15e5, n0)!
	assert r.n[0] == 0 && r.n[3] == 0 && r.n[4] == 0, r.n.str()
	assert math.abs(r.n[1] - 2) <= 1e-15 * 2 && math.abs(r.n[2] - 1) <= 1e-15, r.n.str()
}

// The checks an answer passes before it's returned refuse one a little
// wrong, and pass the right one: a balance off by 1e-9 of the feed, and a
// mole fraction 1e-6 off, which moves μ/RT and ln Q by 1e-6.
fn test_the_checks_refuse_a_wrong_answer() {
	sp := species(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	m := Mixture{
		species: sp
	}
	t, p := 900.0, 1e5
	n0 := [1.0, 2, 0, 0, 0]
	r := gibbs(m, t, p, n0, [])!
	rows := element_rows(m, [])!
	cons := active_rows(rows, [], n0)
	g := standard(m, t, cons.active)!
	rs := [
		kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!,
		kinetics.parse_equation('shift', 'CO + H2O = CO2 + H2', sp)!,
	]
	ln_k := ln_k_gibbs(rs, sp, t)!
	mut nu := [][]f64{}
	for x in rs {
		mut row := []f64{len: sp.len}
		for q, j in x.species {
			row[j] = x.nu[q]
		}
		nu << row
	}
	verify_balance(rows, n0, r)!
	verify_stationary(m, cons, g, p, r, [])!
	verify_k(rs, nu, ln_k, p, r)!
	mut n := r.n.clone()
	n[0] += 1e-9
	mut y := r.y.clone()
	y[0] *= 1 + 1e-6
	off := Result{
		...r
		y: y
	}
	mut refused := 0
	verify_balance(rows, n0, Result{ ...r, n: n }) or { refused++ }
	verify_stationary(m, cons, g, p, off, []) or { refused++ }
	verify_k(rs, nu, ln_k, p, off) or { refused++ }
	assert refused == 3
}

fn refusal(f fn () !Result) string {
	r := f() or { return err.msg() }
	return 'answered ${r.y}'
}

fn test_refused() {
	sp := species(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	m := Mixture{
		species: sp
	}
	n0 := [1.0, 2, 0, 0, 0]
	assert refusal(fn [m, n0] () !Result {
		return gibbs(m, 0, 1e5, n0, [])
	}) == 'equilibrium: no equilibrium at T = 0.0 K and P = 100000 Pa'
	assert refusal(fn [m] () !Result {
		return gibbs(m, 900, 1e5, [1.0, 2], [])
	}) == 'equilibrium: 2 amounts for 5 species'
	assert refusal(fn [m] () !Result {
		return gibbs(m, 900, 1e5, [1.0, -2, 0, 0, 0], [])
	}) == "equilibrium: Water's amount, -2 mol, isn't one"
	assert refusal(fn [m] () !Result {
		return gibbs(m, 900, 1e5, [0.0, 0, 0, 0, 0], [])
	}) == 'equilibrium: nothing to bring to equilibrium'
	// A dependent set: the overall reaction is reforming plus the shift.
	rs := [
		kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!,
		kinetics.parse_equation('shift', 'CO + H2O = CO2 + H2', sp)!,
		kinetics.parse_equation('overall', 'CH4 + 2 H2O = CO2 + 4 H2', sp)!,
	]
	assert refusal(fn [m, n0, rs] () !Result {
		return over_reactions(m, rs, [1.0, 1, 2], 900, 1e5, n0)
	}) == 'equilibrium: the reactions must be independent, but overall = reforming + shift'
	assert refusal(fn [m, n0, rs] () !Result {
		return over_reactions(m, rs[..2], [1.0], 900, 1e5, n0)
	}) == 'equilibrium: 1 equilibrium constants for 2 reactions'
	// A species without the data its μ° needs.
	mut bare := sp.clone()
	bare[3] = kinetics.Species{
		...bare[3]
		gas: none
	}
	assert refusal(fn [bare, n0] () !Result {
		return gibbs(Mixture{ species: bare }, 900, 1e5, n0, [])
	}) == 'equilibrium: Carbon monoxide has no ideal-gas data'
	// Methanol and water at 300 K and 100 bar would condense: one phase only.
	ms := species(['124-38-9', '1333-74-0', '67-56-1', '7732-18-5', '630-08-0'])
	names := ms.map(it.name)
	fl := db.fluid(.pr, names, kij: false)!
	pr := Mixture{
		species: ms
		model:   fl.model
	}
	got := refusal(fn [pr] () !Result {
		return gibbs(pr, 300, 100e5, [1.0, 3, 0, 0, 0], [])
	})
	assert got == 'equilibrium: the equilibrium mixture at 300 K and 1e7 Pa would split into two phases, which needs M6', got
}
