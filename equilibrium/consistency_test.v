// The consistency check (M9 step 7, criterion 4): the papers' rates stop at
// the equilibrium the Gibbs energies give and run towards it; a fitted K
// stops where its own Q is K; a rate on pressures in a real gas stops where
// Q on pressures is K; and typed mistakes are found.
module equilibrium

import components
import flowsheet
import kinetics
import math
import os
import rand
import thermo

const cdb = components.load_default() or { panic(err) }

const kinetics_dir = os.join_path(@VMODROOT, 'testdata', 'kinetics')

fn file_kinetics(name string) kinetics.Kinetics {
	f := flowsheet.load_kinetics(os.join_path(kinetics_dir, name), db: cdb) or { panic(err) }
	return f.kinetics
}

fn gases(names ...string) []kinetics.Species {
	mut out := []kinetics.Species{}
	for n in names {
		r := cdb.find(n) or { panic(err) }
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

// at is a state at T and P from mole fractions by species name.
fn at(k kinetics.Kinetics, t f64, p f64, ys map[string]f64) kinetics.State {
	mut y := []f64{len: k.species.len}
	for name, x in ys {
		y[kinetics.find(k.species, name) or { panic(err) }] = x
	}
	return kinetics.State{
		t: t
		p: p
		y: y
	}
}

// affinity_of is reaction i's ΔG/RT at amounts n, from scratch: ln K from
// its standard properties, ln φ from the model's own state at n.
fn affinity_of(k kinetics.Kinetics, model ?thermo.Model, i int, t f64, p f64, n []f64) f64 {
	r := k.reactions[i]
	st := r.standard(k.species, t) or { panic(err) }
	mut total := 0.0
	for x in n {
		total += x
	}
	mut ln_phi := []f64{len: n.len}
	if m := model {
		ln_phi = (thermo.solve_tp(m, t, p, n, .stable) or { panic(err) }).ln_phi
	}
	mut a := -st.ln_k
	for q, j in r.species {
		a += r.nu[q] * (math.log(n[j] / total) + math.log(p / 1e5) + ln_phi[j])
	}
	return a
}

// rate_of is rate j at amounts n.
fn rate_of(k kinetics.Kinetics, model ?thermo.Model, j int, t f64, p f64, n []f64) f64 {
	mut total := 0.0
	for x in n {
		total += x
	}
	mut s := kinetics.State{
		t: t
		p: p
		y: n.map(it / total)
	}
	if m := model {
		st := thermo.solve_tp(m, t, p, n, .stable) or { panic(err) }
		s = kinetics.State{
			...s
			phi: st.ln_phi.map(math.exp(it))
			z:   st.z
		}
	}
	mut pt := k.point()
	k.fill(s, mut pt) or { panic(err) }
	mut buf := k.buffer()
	mut out := []f64{len: k.rates.len}
	k.rates(&pt, mut buf, mut out)
	return out[j]
}

// verify holds a judged rate's stop to what it claims: ΔG/RT there is its
// reaction's, from scratch, to 1e-12, and the rate changes sign across it,
// 1e-9 of the amounts either way along its reaction.
fn verify(k kinetics.Kinetics, model ?thermo.Model, t f64, p f64, j Judged) {
	assert j.found, '${k.reactions[j.reaction].name}: no stop'
	a := affinity_of(k, model, j.reaction, t, p, j.at)
	assert math.abs(j.stops - a) <= 1e-12, '${j.stops} vs ${a}'
	r := k.reactions[j.reaction]
	mut total := 0.0
	for x in j.at {
		total += x
	}
	mut before := j.at.clone()
	mut after := j.at.clone()
	for q, s in r.species {
		before[s] -= 1e-9 * total * r.nu[q]
		after[s] += 1e-9 * total * r.nu[q]
	}
	rb := rate_of(k, model, j.rate, t, p, before)
	ra := rate_of(k, model, j.rate, t, p, after)
	assert rb * ra <= 0, '${k.reactions[j.reaction].name}: ${rb} and ${ra} either side of its stop'
}

// Xu and Froment's rates, K from the Gibbs energies, stop at the equilibrium
// over reforming and the shift, overall's included, at every state, and run
// towards it: 60 seeded states from 600 to 1200 K and 1 to 40 bar. Their
// stops are within the check's tolerance of zero in ΔG/RT, 1e-10: measured,
// 2.7e-15 at worst.
fn test_the_papers_rates_stop_at_equilibrium() {
	k := file_kinetics('xu_froment.toml')
	rand.seed([u32(7), 11])
	mut worst := 0.0
	for _ in 0 .. 60 {
		mut y := []f64{len: 5, init: rand.f64_in_range(0.01, 1) or { panic(err) }}
		mut total := 0.0
		for x in y {
			total += x
		}
		y = y.map(it / total)
		s := kinetics.State{
			t: rand.f64_in_range(600, 1200) or { panic(err) }
			p: (rand.f64_in_range(1, 40) or { panic(err) }) * 1e5
			y: y
		}
		c := consistency(&k, none, s) or { panic(err) }
		assert c.over == [0, 1]
		assert c.judged.len == 3
		for j in c.judged {
			assert j.skipped == ''
			verify(k, none, s.t, s.p, j)
			worst = math.max(worst, math.abs(j.stops))
			assert j.consistent(), '${k.reactions[j.reaction].name} at ${s}: stops at ${j.stops}, runs ${j.at_state} where ΔG/RT is ${j.affinity}'
			assert j.at_state * j.affinity <= 0
		}
	}
	println('Xu and Froment: their rates stop within ΔG/RT ${worst} of equilibrium')
	assert worst <= stop_tol
	s := at(k, 900, 15e5, {
		'CH4': 0.2
		'H2O': 0.45
		'H2':  0.2
		'CO':  0.1
		'CO2': 0.05
	})
	assert check(&k, none, s).map(it.str()) == [
		"note: reforming's, shift's and overall's rates are zero at the equilibrium this state reaches, over reforming and shift, and run towards it at the state",
	]
}

// A fitted K stops where its own Q is K: ΔG/RT = ln K_fitted − ln K there,
// in K's activities at 1 bar. Vanden Bussche and Froment's are Graaf's, in
// bar; the owner's DWSIM case's are in atm, the overall's a factor of about
// 0.2 from the Gibbs energies' at 900 K (step 6), which is warned of.
fn test_a_fitted_k_stops_where_its_own_q_is_k() {
	vbf := file_kinetics('vanden_bussche_froment.toml')
	t := 513.15
	s := at(vbf, t, 50e5, {
		'CO2':   0.2
		'H2':    0.65
		'CH3OH': 0.04
		'H2O':   0.03
		'CO':    0.08
	})
	c := consistency(&vbf, none, s)!
	fitted := [
		math.ln10 * (3066 / t - 10.592),
		math.ln10 * (2.029 - 2073 / t),
	]
	assert c.judged.len == 2
	for j in c.judged {
		verify(vbf, none, t, 50e5, j)
		want := fitted[j.reaction] - vbf.reactions[j.reaction].standard(vbf.species, t)!.ln_k
		assert math.abs(j.stops - want) <= 1e-9, '${vbf.reactions[j.reaction].name}: ${j.stops} vs ${want}'
	}
	dw := file_kinetics('dwsim_reforming.toml')
	t2 := 900.0
	s2 := at(dw, t2, 15e5, {
		'CH4': 0.2
		'H2O': 0.45
		'H2':  0.2
		'CO':  0.1
		'CO2': 0.05
	})
	c2 := consistency(&dw, none, s2)!
	atm := 2 * math.log(1.01325)
	fitted2 := [30.42 - 27106 / t2 + atm, -3.798 + 4160 / t2, 34.218 - 31266 / t2 + atm]
	for j in c2.judged {
		verify(dw, none, t2, 15e5, j)
		want := fitted2[j.reaction] - dw.reactions[j.reaction].standard(dw.species, t2)!.ln_k
		assert math.abs(j.stops - want) <= 1e-9, '${dw.reactions[j.reaction].name}: ${j.stops} vs ${want}'
	}
	got := c2.findings(&dw).map(it.str())
	assert got.len == 3, got.str()
	assert got[2].starts_with("warning: overall's rate isn't zero at the equilibrium this state reaches: on its fitted K it stops where Q/K is 0.20173 (ΔG/RT -1.6008)"), got[2]
}

// On fugacities, a rate stops at the equilibrium in Peng–Robinson's gas
// too; on partial pressures it stops where Q on pressures is K, so ΔG/RT
// there is Σν ln φ, 0.052 at 900 K and 30 bar; and on
// concentrations, cᵢ = yᵢP/(ZRT) with K·(p°/RT)^Δn, where Q on pressures is
// K·Z^Δn, so ΔG/RT is Σν ln φ + Δn ln Z.
fn test_a_real_gas_on_pressures() {
	sp := gases('methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide')
	r := kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!
	fl := cdb.fluid(.pr, sp.map(it.name), kij: false)!
	model := ?thermo.Model(fl.model)
	t, p := 900.0, 30e5
	for basis in [kinetics.Basis.fugacity, .pressure, .concentration] {
		unit := if basis == .concentration { 'm^3/(mol*s)' } else { 'mol/(kg*s*Pa^2)' }
		k := kinetics.build(sp, [r], [], [
			kinetics.Rated{0, kinetics.PowerLaw{
				k:           kinetics.Fixed{'1e-12 ${unit}'}
				orders:      [kinetics.Power{'CH4', 1}, kinetics.Power{'H2O', 1}]
				equilibrium: kinetics.Equilibrium{
					kind: .thermo
				}
				basis:       basis
				per:         if basis == .concentration { kinetics.Per.volume } else { .mass }
			}},
		], [])!
		s := at(k, t, p, {
			'CH4': 0.25
			'H2O': 0.6
			'H2':  0.1
			'CO':  0.05
		})
		c := consistency(&k, model, s)!
		j := c.judged[0]
		verify(k, model, t, p, j)
		if basis == .fugacity {
			assert j.consistent(), '${j.stops}'
			continue
		}
		st := thermo.solve_tp(fl.model, t, p, j.at, .stable)!
		mut sum := 0.0
		for q, i in r.species {
			sum += r.nu[q] * st.ln_phi[i]
		}
		if basis == .concentration {
			sum += 2 * math.log(st.z)
		}
		assert math.abs(j.stops - sum) <= 1e-9, '${basis}: ${j.stops} vs ${sum}'
		println('on ${basis} in Peng–Robinson at 900 K and 30 bar, it stops at ΔG/RT ${j.stops}')
		assert math.abs(j.stops) > 1e-3
		assert c.findings(&k).map(it.severity) == [kinetics.Severity.note]
	}
}

// A rate typed as an expression is held to equilibrium as a structured one
// is: one on mole fractions without (p°/P)^Δn stops at ΔG/RT = 2 ln(P/p°);
// one with CO₂ for CO stops elsewhere; one with a factor of 0.99999, a
// constant typed to five figures, stops at ΔG/RT = +1e-5; one with its sign
// reversed is zero at equilibrium but runs away from it; and one without a
// backward term is irreversible as written.
fn test_typed_mistakes_are_found() {
	sp := gases('methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide')
	r := kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!
	t, p := 900.0, 15e5
	cases := {
		'y[CH4]*y[H2O] - y[CO]*y[H2]^3/Keq[reforming]':                       "warning: reforming's rate isn't zero at the equilibrium this state reaches: it stops where Q/K is 225 (ΔG/RT +5.4161), so a reactor on it heads for a different equilibrium"
		'y[CH4]*y[H2O] - y[CO2]*y[H2]^3/(Keq[reforming]*(p_std/P)^2)':        "warning: reforming's rate isn't zero at the equilibrium this state reaches: it stops where Q/K is"
		'y[CO]*y[H2]^3/(Keq[reforming]*(p_std/P)^2) - y[CH4]*y[H2O]':         "warning: reforming's rate runs backwards at the state, where its reaction's ΔG/RT is -"
		'y[CH4]*y[H2O] - 0.99999*y[CO]*y[H2]^3/(Keq[reforming]*(p_std/P)^2)': "note: reforming's rate isn't zero at the equilibrium this state reaches: it stops where ΔG/RT is +1e-5, so a reactor on it heads for a different equilibrium"
		'y[CH4]*y[H2O]':                                                      "note: reforming's rate is never zero along its reaction from the equilibrium this state reaches: it runs one way until a species runs out, irreversible as written"
		'y[CH4]*y[H2O] - y[CO]*y[H2]^3/(Keq[reforming]*(p_std/P)^2) + 0':     "note: reforming's rate is zero at the equilibrium this state reaches, over reforming, and runs towards it at the state"
	}
	for text, want in cases {
		k := kinetics.build(sp, [r], [], [
			kinetics.Rated{0, kinetics.Expression{
				text: '1e-3 mol/(kg*s) * (${text})'
			}},
		], [])!
		s := at(k, t, p, {
			'CH4': 0.25
			'H2O': 0.6
			'H2':  0.1
			'CO':  0.04
			'CO2': 0.01
		})
		got := check(&k, none, s).map(it.str())
		assert got.len == 1 && got[0].starts_with(want), '${text}: ${got}'
		c := consistency(&k, none, s)!
		if c.judged[0].found {
			verify(k, none, t, p, c.judged[0])
		}
	}
}

// Where a rate is zero twice along its reaction, it stops at the nearer, the
// way it runs: here at Q/K = 1.0001, not 1.001, a band the probes must not
// step over.
fn test_the_nearest_zero() {
	sp := gases('methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide')
	r := kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!
	t, p := 900.0, 15e5
	q := 'y[CO]*y[H2]^3/(y[CH4]*y[H2O]*Keq[reforming]*(p_std/P)^2)'
	k := kinetics.build(sp, [r], [], [
		kinetics.Rated{0, kinetics.Expression{
			text: '1e-3 mol/(kg*s) * (1.0001 - ${q}) * (1.001 - ${q})'
		}},
	], [])!
	s := at(k, t, p, {
		'CH4': 0.25
		'H2O': 0.6
		'H2':  0.1
		'CO':  0.04
		'CO2': 0.01
	})
	c := consistency(&k, none, s)!
	j := c.judged[0]
	verify(k, none, t, p, j)
	assert math.abs(j.stops - math.log(1.0001)) <= 1e-9, '${j.stops}'
}

// A fitted K 1.5 times the Gibbs energies' stops at ΔG/RT = ln 1.5, so at a
// state between, at ΔG/RT = 0.2, its rate runs forwards against its ΔG.
// That's its K, already noted, not running away from equilibrium.
fn test_a_fitted_k_between() {
	sp := gases('methane', 'water', 'hydrogen', 'carbon monoxide', 'carbon dioxide')
	r := kinetics.parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp)!
	t, p := 900.0, 15e5
	kt := math.exp(r.standard(sp, t)!.ln_k)
	k := kinetics.build(sp, [r], [], [
		kinetics.Rated{0, kinetics.PowerLaw{
			k:           kinetics.Fixed{'1e-12 mol/(kg*s*Pa^2)'}
			orders:      [kinetics.Power{'CH4', 1}, kinetics.Power{'H2O', 1}]
			equilibrium: kinetics.Equilibrium{
				kind:        .correlation
				correlation: kinetics.Fixed{'${1.5 * kt} * p_std^2'}
			}
		}},
	], [])!
	// y_CO·y_H2³/(y_CH4·y_H2O)·(P/p°)² = K·e^0.2.
	y_ch4, y_h2o, y_co := 0.1, 0.3, 0.05
	y_h2 := math.cbrt(kt * math.exp(0.2) * math.pow(1e5 / p, 2) * y_ch4 * y_h2o / y_co)
	y_co2 := 1 - y_ch4 - y_h2o - y_co - y_h2
	assert y_co2 > 0
	s := kinetics.State{
		t: t
		p: p
		y: [y_ch4, y_h2o, y_h2, y_co, y_co2]
	}
	c := consistency(&k, none, s)!
	j := c.judged[0]
	verify(k, none, t, p, j)
	assert math.abs(j.stops - math.log(1.5)) <= 1e-9, '${j.stops}'
	assert math.abs(j.affinity - 0.2) <= 1e-12, '${j.affinity}'
	assert j.at_state * j.affinity > 0
	assert c.findings(&k).map(it.str()) == [
		"note: reforming's rate isn't zero at the equilibrium this state reaches: on its fitted K it stops where Q/K is 1.5 (ΔG/RT +0.40547), so a reactor on it heads for a different equilibrium",
	]
}

// What can't be judged says why: a state without carbon reaches an
// equilibrium without methane, and a state that's wrong is a note.
fn test_what_cannot_be_judged() {
	k := file_kinetics('xu_froment.toml')
	s := at(k, 900, 15e5, {
		'H2O': 0.5
		'H2':  0.5
	})
	assert check(&k, none, s).map(it.str()) == [
		"note: reforming's rate can't be held to equilibrium: no Methane at the equilibrium this state reaches",
		"note: shift's rate can't be held to equilibrium: no Carbon monoxide at the equilibrium this state reaches",
		"note: overall's rate can't be held to equilibrium: no Methane at the equilibrium this state reaches",
	]
	bad := kinetics.State{
		t: 900
		p: 15e5
		y: [0.5, 0.5]
	}
	assert check(&k, none, bad).map(it.str()) == [
		"note: the rates can't be held to equilibrium from this state: equilibrium: 2 mole fractions for 5 species",
	]
}
