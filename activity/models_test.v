// M6 steps 2 and 3, criterion 1: Wilson, NRTL, UNIQUAC and the three
// UNIFACs against thermo 0.6.1 and, for the binary models, mpmath at 40
// digits, on ChemSep's own banks and tables (testdata/m6_activity.json, from
// tools/gen_m6_activity.py); and each model's Gibbs–Duhem, its ln γ against
// its own Gibbs energy, and a pure component's γ, checked from scratch.
import activity
import components
import json2
import math
import os

struct Answer {
	ok       bool
	ln_gamma []f64
	dt       []f64
	dtt      []f64
	dn       []f64
	ge       f64
	he       f64
	cp       f64
}

struct Case {
	model  string
	cas    []string
	t      f64 @[json: 'T']
	x      []f64
	thermo Answer
	mpmath Answer
}

struct Source {
	r_thermo f64 @[json: 'R_thermo']
	r_mpmath f64 @[json: 'R_mpmath']
}

struct Fixture {
	source Source
	cases  []Case
}

fn fixture() !Fixture {
	return json2.decode[Fixture](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_activity.json'))!)!
}

fn model_of(db components.Database, name string, cas []string) !activity.Model {
	return match name {
		'NRTL' { activity.Model(db.nrtl(cas)!.model) }
		'UNIQUAC' { activity.Model(db.uniquac(cas)!.model) }
		'Wilson' { activity.Model(db.wilson(cas)!.model) }
		'UNIFAC' { activity.Model(db.unifac(.original, cas)!) }
		'UNIFAC-LLE' { activity.Model(db.unifac(.lle, cas)!) }
		else { activity.Model(db.unifac(.dortmund, cas)!) }
	}
}

// gap is |a − b| over max(|b|, floor).
fn gap(a f64, b f64, floor f64) f64 {
	return math.abs(a - b) / math.max(math.abs(b), floor)
}

struct Worst {
mut:
	m map[string]f64
}

fn (mut w Worst) note(key string, v f64) {
	if v > (w.m[key] or { 0.0 }) {
		w.m[key] = v
	}
}

fn test_the_models_against_thermo_and_mpmath() {
	f := fixture()!
	db := components.load_default()!
	mut worst := Worst{}
	mut counted := map[string]int{}
	for c in f.cases {
		m := model_of(db, c.model, c.cas)!
		mut w := m.work()
		activity.eval(m, c.t, c.x, .t | .n, mut w)!
		e := activity.excess(c.t, c.x, &w)
		nc := c.x.len
		counted[c.model] = (counted[c.model] or { 0 }) + 1
		// Each row of dn by its own size, since an entry can be the small
		// difference of large terms.
		mut row := []f64{len: nc}
		for i in 0 .. nc {
			for j in 0 .. nc {
				row[i] = math.max(row[i], math.abs(w.dn[i * nc + j]))
			}
		}
		for src in [c.mpmath, c.thermo] {
			if !src.ok {
				continue
			}
			tag := if src == c.mpmath { 'mpmath' } else { 'thermo' }
			// thermo's energies are on its R, N_A·k_B to every digit; reflux's
			// r_gas is that rounded, 1.84e-11 below.
			scale := if tag == 'thermo' { activity.r_gas / f.source.r_thermo } else { 1.0 }
			for i in 0 .. nc {
				worst.note('${tag} ln γ', gap(w.ln_gamma[i], src.ln_gamma[i], 1.0))
				worst.note('${tag} dT', gap(w.dt[i], src.dt[i], 1e-4))
				if src.dtt.len > 0 {
					worst.note('${tag} dT²', gap(w.dtt[i], src.dtt[i], 1e-6))
				}
				for j in 0 .. nc {
					worst.note('${tag} dn', math.abs(w.dn[i * nc + j] - src.dn[i * nc + j]) / math.max(row[i],
						1e-3))
				}
			}
			worst.note('${tag} Gᴱ', gap(e.g, src.ge * scale, 10.0))
			worst.note('${tag} Hᴱ', gap(e.h, src.he * scale, 10.0))
			worst.note('${tag} Cpᴱ', gap(e.cp, src.cp * scale, 1.0))
		}
	}
	println('cases: ${counted}')
	for k, v in worst.m {
		println('  ${k}: ${v:.2e}')
	}
	// Measured, under tcc, against mpmath: ln γ 1.9e-15, its T-derivatives
	// 4.4e-15 and 2.1e-15, dn 4.9e-14 of its row, the energies 1.5e-14; a sum
	// of a few terms each rounded, so a few times ε, and dn's a difference
	// (∂ᵢⱼ − Σₖ xₖ∂ᵢₖ) of terms its row's size. Against thermo, whose own
	// floats are no closer: dn 2.5e-13, its float formula cancelling where an
	// entry is small (4e-13 absolute, against mpmath, in the generator).
	limits := {
		'mpmath ln γ': 1e-14
		'mpmath dT':   2e-14
		'mpmath dT²':  2e-14
		'mpmath dn':   2e-13
		'mpmath Gᴱ':   1e-13
		'mpmath Hᴱ':   1e-13
		'mpmath Cpᴱ':  1e-13
		'thermo ln γ': 1e-14
		'thermo dT':   1e-13
		'thermo dT²':  1e-13
		'thermo dn':   1e-12
		'thermo Gᴱ':   2e-13
		'thermo Hᴱ':   2e-13
		'thermo Cpᴱ':  1e-13
	}
	for k, v in worst.m {
		assert v <= limits[k] or { 0.0 }, '${k}: ${v}'
	}
	assert worst.m.len == limits.len
	assert counted == {
		'NRTL':            20
		'UNIQUAC':         20
		'Wilson':          20
		'UNIFAC':          28
		'UNIFAC-LLE':      28
		'UNIFAC-Dortmund': 28
	}
}

// Each model's own consistency, from scratch at every fixture state: dn
// symmetric and Gibbs–Duhem's Σᵢ xᵢ·dn[i][j] = 0, each to rounding; and ln γᵢ
// as ∂(nGᴱ/RT)/∂nᵢ by central differences of nGᴱ/RT = Σₖ nₖ ln γₖ, with
// Richardson's extrapolation, which uses none of the jets the derivatives
// come from.
fn test_each_model_is_consistent() {
	f := fixture()!
	db := components.load_default()!
	mut worst_gd, mut worst_sym, mut worst_g := 0.0, 0.0, 0.0
	for c in f.cases {
		m := model_of(db, c.model, c.cas)!
		mut w := m.work()
		activity.eval(m, c.t, c.x, .n, mut w)!
		nc := c.x.len
		for j in 0 .. nc {
			mut s, mut size := 0.0, 0.0
			for i in 0 .. nc {
				s += c.x[i] * w.dn[i * nc + j]
				size += math.abs(c.x[i] * w.dn[i * nc + j])
				worst_sym = math.max(worst_sym, math.abs(w.dn[i * nc + j] - w.dn[j * nc + i]) / math.max(math.abs(w.dn[i * nc + j]),
					1e-3))
			}
			worst_gd = math.max(worst_gd, math.abs(s) / math.max(size, 1e-300))
		}
		lng := w.ln_gamma.clone()
		for i in 0 .. nc {
			if c.x[i] == 0.0 {
				continue
			}
			d1 := ng_slope(m, c.t, c.x, i, 1e-4)!
			d2 := ng_slope(m, c.t, c.x, i, 5e-5)!
			worst_g = math.max(worst_g, math.abs((4.0 * d2 - d1) / 3.0 - lng[i]) / math.max(math.abs(lng[i]),
				1e-2))
		}
	}
	println('Gibbs–Duhem ${worst_gd:.1e}, symmetry ${worst_sym:.1e}, ln γ against nGᴱ ${worst_g:.1e}')
	// Measured: 2.4e-14 of the sum's terms, 3.1e-13 and 1.6e-9. dn's entries
	// are differences of the formula's raw derivatives, which grow as 1/xᵢ
	// for a scarce component, so their rounding is relative to those, not to
	// dn. The last is the differences' own error: h⁴ truncation left by
	// Richardson, and ε·nGᴱ/h rounding.
	assert worst_gd <= 1e-13
	assert worst_sym <= 1e-12
	assert worst_g <= 1e-8
}

// Sᴱ as −∂Gᴱ/∂T, by central differences of Gᴱ, at every fixture state; and
// a Work used at one T, then another, giving a fresh Work's bits there, so
// that what it keeps until T moves is let go when T does.
fn test_the_entropy_and_what_is_kept() {
	f := fixture()!
	db := components.load_default()!
	mut worst := 0.0
	for c in f.cases {
		m := model_of(db, c.model, c.cas)!
		mut w := m.work()
		activity.eval(m, c.t, c.x, .t, mut w)!
		e := activity.excess(c.t, c.x, &w)
		mut g := []f64{len: 2}
		for k, dt in [0.01, -0.01] {
			mut wk := m.work()
			activity.eval(m, c.t + dt, c.x, .t, mut wk)!
			g[k] = activity.excess(c.t + dt, c.x, &wk).g
		}
		s := -(g[0] - g[1]) / 0.02
		worst = math.max(worst, math.abs(e.s - s) / math.max(math.abs(e.s), 0.1))
		// The same Work at a second T.
		activity.eval(m, c.t + 7.0, c.x, .t | .n, mut w)!
		mut fresh := m.work()
		activity.eval(m, c.t + 7.0, c.x, .t | .n, mut fresh)!
		assert w.ln_gamma == fresh.ln_gamma && w.dt == fresh.dt && w.dtt == fresh.dtt && w.dn == fresh.dn, c.model
	}
	println('Sᴱ against −∂Gᴱ/∂T: ${worst:.1e}')
	// Measured: 2.7e-8, the central difference's own error at a step of
	// 0.01 K: truncation h²·Gᴱ‴/6 and rounding ε·Gᴱ/h.
	assert worst <= 1e-7
}

// A T a model refuses leaves nothing kept (found after step 3): params
// writes a Work's parameters as it goes, so a refusal part of the way
// through leaves some of them at the T refused. Each model here refuses
// 150 K on its second pair of parameters, the first written by then, and
// must then give a fresh Work's bits at the T it had, and refuse 150 K again.
fn test_a_refused_temperature_leaves_nothing_kept() {
	names := ['a', 'b']
	first := activity.plain('(100 K)/T', 'first')!
	below := activity.plain('sqrt(T/(300 K) - 1)', 'second')!
	negative := activity.plain('T/(300 K) - 1', 'second')!
	alpha := activity.plain('0.3', 'α')!
	table := activity.Table{
		variant:      .original
		subgroups:    {
			1: activity.Subgroup{1, 'CH3', 1, 0.9011, 0.848, 15.035}
			2: activity.Subgroup{2, 'H2O', 7, 0.92, 1.4, 18.015}
		}
		mains:        {
			1: 'CH2'
			7: 'H2O'
		}
		interactions: {
			activity.interaction_key(1, 7): activity.Interaction{1000.0, 0.0, 0.0}
			// ψ is exp(−1000) at 150 K, which is 0 in a double.
			activity.interaction_key(7, 1): activity.Interaction{150000.0, 0.0, 0.0}
		}
	}
	models := [
		activity.Model(activity.new_nrtl(names, [
			activity.NrtlPair{0, 1, first, below, alpha},
		])!),
		activity.Model(activity.new_uniquac(names, [
			0.92,
			1.43,
		], [
			1.4,
			1.43,
		], [
			activity.UniquacPair{0, 1, first, below},
		])!),
		activity.Model(activity.new_wilson(names, [
			activity.WilsonPair{0, 1, first, negative},
		])!),
		activity.Model(activity.new_unifac(table, names, [
			{
				1: 1
			},
			{
				2: 1
			},
		])!),
	]
	mut wrong := []string{}
	for m in models {
		mut fresh := m.work()
		activity.eval(m, 350.0, [0.4, 0.6], .t | .n, mut fresh)!
		// A Work that held 350 K, and one that held nothing.
		for held in [true, false] {
			mut w := m.work()
			if held {
				activity.eval(m, 350.0, [0.4, 0.6], .t | .n, mut w)!
			}
			for _ in 0 .. 2 {
				mut refused := false
				activity.eval(m, 150.0, [0.4, 0.6], .t | .n, mut w) or { refused = true }
				assert refused, '${m.name()} accepted 150 K'
			}
			activity.eval(m, 350.0, [0.4, 0.6], .t | .n, mut w)!
			if w.ln_gamma != fresh.ln_gamma || w.dt != fresh.dt || w.dtt != fresh.dtt
				|| w.dn != fresh.dn {
				wrong << '${m.name()}: ${w.ln_gamma}, not ${fresh.ln_gamma}'
			}
		}
	}
	assert wrong == [], wrong.join('; ')
}

// A pair the bank doesn't list, by the caller's leave (D7): ideal in NRTL
// and Wilson, and in UNIQUAC only without a residual term, its ln γ then
// the combinatorial term, written here in UNIFAC's way rather than
// UNIQUAC's.
fn test_an_unlisted_pair() {
	db := components.load_default()!
	cas := ['7732-18-5', '110-54-3'] // water and n-hexane, in none of the banks
	x := [0.4, 0.6]
	nrtl := db.nrtl(cas, ideal: true)!
	wilson := db.wilson(cas, ideal: true)!
	assert nrtl.pairs[0].source == .ideal && wilson.pairs[0].source == .ideal
	for m in [activity.Model(nrtl.model), activity.Model(wilson.model)] {
		mut w := m.work()
		activity.eval(m, 298.15, x, .t, mut w)!
		for i in 0 .. 2 {
			// Σ xⱼ is 1 to a unit in the last place.
			assert math.abs(w.ln_gamma[i]) <= 4.5e-16 && w.dt[i] == 0.0, m.name()
		}
	}
	u := db.uniquac(cas, ideal: true)!
	assert u.pairs[0].source == .no_residual
	r, q := u.model.r, u.model.q
	mut w := u.model.work()
	activity.eval(u.model, 298.15, x, .t, mut w)!
	rbar, qbar := x[0] * r[0] + x[1] * r[1], x[0] * q[0] + x[1] * q[1]
	for i in 0 .. 2 {
		v, f := r[i] / rbar, q[i] / qbar
		want := 1.0 - v + math.log(v) - 5.0 * q[i] * (1.0 - v / f + math.log(v / f))
		// Terms of up to 20, each rounded, cancel to ln γ.
		assert math.abs(w.ln_gamma[i] - want) <= 2e-14, '${i}: ${w.ln_gamma[i]}, not ${want}'
		if i == 0 {
			assert want > 0.1, "water's combinatorial term here"
		}
		assert w.dt[i] == 0.0
	}
}

// ng_slope is ∂(nGᴱ/RT)/∂nᵢ at x, one mole in all, by a central difference
// of step h in nᵢ.
fn ng_slope(m activity.Model, t f64, x []f64, i int, h f64) !f64 {
	mut w := m.work()
	mut g := []f64{len: 2}
	for side, sign in [1.0, -1.0] {
		mut n := x.clone()
		n[i] += sign * h
		total := 1.0 + sign * h
		xs := n.map(it / total)
		activity.eval(m, t, xs, activity.ln_gamma_only, mut w)!
		mut s := 0.0
		for k, v in n {
			s += v * w.ln_gamma[k]
		}
		g[side] = s
	}
	return (g[0] - g[1]) / (2.0 * h)
}

// A pure component's γ is exactly 1, and its T-derivatives 0, in every
// model: each formula's terms cancel to the bit there (UNIFAC's by taking
// each pure group's ln Γ by the mixture's own arithmetic).
fn test_a_pure_component_is_ideal() {
	db := components.load_default()!
	cas := ['67-64-1', '67-56-1', '7732-18-5']
	for name in ['NRTL', 'UNIQUAC', 'Wilson', 'UNIFAC', 'UNIFAC-LLE', 'UNIFAC-Dortmund'] {
		m := model_of(db, name, cas)!
		mut w := m.work()
		for i in 0 .. 3 {
			x := []f64{len: 3, init: if index == i { 1.0 } else { 0.0 }}
			activity.eval(m, 333.15, x, .t, mut w)!
			assert w.ln_gamma[i] == 0.0 && w.dt[i] == 0.0 && w.dtt[i] == 0.0, '${name}, ${cas[i]}: ${w.ln_gamma[i]}'
		}
	}
}

fn test_refusals() {
	db := components.load_default()!
	m := db.nrtl(['67-56-1', '7732-18-5'])!.model
	mut w := m.work()
	for case in [
		Refusal{0.0, [0.5, 0.5], 'a temperature of 0'},
		Refusal{300.0, [-0.1, 1.1], 'a mole fraction of -0.1'},
		Refusal{300.0, [0.5, 0.25, 0.25], '3 mole fractions for 2 components'},
		Refusal{300.0, [math.inf(1), 0.5], 'a mole fraction of +inf'},
	] {
		activity.eval(m, case.t, case.x, .t, mut w) or {
			assert err.msg().contains(case.why), err.msg()
			continue
		}
		assert false, 'accepted ${case}'
	}
	// A pair the bank doesn't list: refused by name, unless asked for ideal.
	missing := ['67-56-1', '7732-18-5', '7440-37-1']
	if _ := db.nrtl(missing) {
		assert false, 'built NRTL with an unlisted pair'
	} else {
		assert err.msg().contains("doesn't list"), err.msg()
	}
	ideal := db.nrtl(missing, ideal: true)!
	assert ideal.pairs.filter(it.source == .ideal).len == 2
	// A τ that can't be evaluated at a T, named with its text.
	bad := activity.new_nrtl(['a', 'b'], [
		activity.NrtlPair{0, 1, activity.plain('ln((T - 300 K) / (1 K))', 'τ')!, activity.plain('0', 'τ')!, activity.plain('0.3', 'α')!},
	])!
	mut wb := bad.work()
	mut refused := false
	activity.eval(bad, 250.0, [0.5, 0.5], .t, mut wb) or {
		assert err.msg().contains('ln((T - 300 K) / (1 K))'), err.msg()
		refused = true
	}
	assert refused, 'a τ of ln of a negative number accepted'

	// A τ with a dimension.
	if _ := activity.plain('300 K', 'τ') {
		assert false, 'a τ in kelvins accepted'
	}
	// UNIFAC: a subgroup the table hasn't, and two main groups with no
	// parameters between them.
	table := activity.Table{
		variant:      .original
		subgroups:    {
			1: activity.Subgroup{1, 'CH3', 1, 0.9011, 0.848, 15.035}
			2: activity.Subgroup{2, 'H2O', 7, 0.92, 1.4, 18.015}
		}
		mains:        {
			1: 'CH2'
			7: 'H2O'
		}
		interactions: {
			activity.interaction_key(1, 7): activity.Interaction{1318.0, 0.0, 0.0}
		}
	}
	if _ := activity.new_unifac(table, ['a', 'b'], [{
		1: 1
	}, {
		3: 1
	}]) {
		assert false, 'a subgroup the table has not accepted'
	} else {
		assert err.msg().contains("subgroup 3 isn't in the table"), err.msg()
	}
	if _ := activity.new_unifac(table, ['a', 'b'], [{
		1: 2
	}, {
		2: 1
	}]) {
		assert false, 'a missing parameter accepted'
	} else {
		assert err.msg().contains('no parameters for main group 7 (H2O) with 1 (CH2)'), err.msg()
	}
	// UNIFAC: a compound with no groups for the variant.
	if _ := db.unifac(.original, ['7440-37-1', '7732-18-5']) {
		assert false, 'UNIFAC built for argon'
	} else {
		assert err.msg().contains('no UNIFAC groups'), err.msg()
	}
}

struct Refusal {
	t   f64
	x   []f64
	why string
}
