// M5b step 4, criterion 4: a phase's transport, by its label and D9's
// rules, end to end against the pieces put together by
// tools/gen_m5b_transport.py (testdata/m5b_transport.json); what's left
// out of a liquid and why; water's own by IAPWS; and every phase of seeded
// flashes of every fluid kind with its values, or why not.
import components
import flash
import json2
import math
import os
import thermo
import transport
import water

fn fixture() !map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5b_transport.json'))!)!.as_map()
}

fn floats(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

fn label_of(s string) thermo.PhaseLabel {
	return match s {
		'vapour' { thermo.PhaseLabel.vapor }
		'liquid2' { .liquid2 }
		else { .liquid }
	}
}

fn check(name string, what string, got transport.Result, want json2.Any, mut worst []f64) {
	if want is json2.Null {
		assert !got.ok() && got.why != '', '${name}: ${what} ${got.value} where there should be none'
		return
	}
	assert got.ok(), '${name}: ${what}: ${got.why}'
	e := math.abs(got.value / want.f64() - 1)
	worst[0] = math.max(worst[0], e)
	assert e <= 1e-13, '${name}: ${what} ${got.value} vs ${want.f64()}'
}

// Each phase, liquids and gases, against the generator's own putting
// together of chemicals' and mpmath's pieces as D4 to D9 have it, within
// 1e-13: the pieces' own agreements, compounded (measured: 1.4e-14). A
// liquid with 1% hydrogen is water's own; carbon dioxide above its Tc
// labelled liquid takes the gas's methods.
fn test_phases_against_fixture() {
	db := components.load_default()!
	f := fixture()!
	mut worst := [0.0]
	for c in f['phases']!.as_array() {
		m := c.as_map()
		names := m['cas']!.as_array().map(it.str())
		mix := db.transport(names)!
		name := m['name']!.str()
		got := mix.phase(label_of(m['label']!.str()), m['t']!.f64(), m['p']!.f64(), floats(m['x']!),
			m['v']!.f64())
		check(name, 'μ', got.viscosity, m['mu']!, mut worst)
		check(name, 'λ', got.conductivity, m['k']!, mut worst)
		check(name, 'σ', got.surface_tension, m['sigma']!, mut worst)
	}
	println('phases within ${worst[0]} of the fixture')
}

// D9's rules, each said in the phase's notes or its reasons: a compound
// above its Tc in a liquid is left out, and noted with its share; more than
// dissolved_limit of them leave the liquid none; a phase of compounds all
// above their Tc takes the gas's methods, whatever its label; a gas has no
// surface tension.
fn test_dissolved_and_supercritical() {
	db := components.load_default()!
	mix := db.transport(['water', 'hydrogen'])!
	one := mix.phase(.liquid, 320, 3e6, [0.99, 0.01], 1e-4)
	assert one.viscosity.ok() && one.conductivity.ok() && one.surface_tension.ok()
	assert one.notes.any(it.contains("Hydrogen, above its Tc, left out: 0.01 of the liquid's moles (D9)")), '${one.notes}'
	alone := db.transport(['water'])!.phase(.liquid, 320, 3e6, [1.0], 1e-4)
	assert alone.viscosity.value == one.viscosity.value
	assert alone.conductivity.value == one.conductivity.value
	assert alone.surface_tension.value == one.surface_tension.value
	// What's kept is renormalised: methanol and water with 1% hydrogen are
	// methanol and water at 0.3 and 0.7, to the rounding of the division
	// (a planted mistake that kept 0.297 and 0.693 survived until this).
	three := db.transport(['methanol', 'water', 'hydrogen'])!.phase(.liquid, 320, 3e6, [
		0.297,
		0.693,
		0.01,
	], 1e-4)
	two := db.transport(['methanol', 'water'])!.phase(.liquid, 320, 3e6, [0.3, 0.7], 1e-4)
	for r in [[three.viscosity.value, two.viscosity.value], [three.conductivity.value,
		two.conductivity.value], [three.surface_tension.value, two.surface_tension.value]] {
		assert math.abs(r[0] / r[1] - 1) <= 1e-14, '${r[0]} vs ${r[1]}'
	}
	many := mix.phase(.liquid, 320, 3e6, [0.8, 0.2], 1e-4)
	for r in [many.viscosity, many.conductivity, many.surface_tension] {
		assert !r.ok() && math.is_nan(r.value)
		assert r.why.contains("Hydrogen, above its Tc, is 0.2 of the liquid's moles, more than 0.1 (D9)"), r.why
	}
	co2 := db.transport(['carbon dioxide'])!
	v := 6.409559770234321e-05
	as_liquid := co2.phase(.liquid, 320, 1.5e7, [1.0], v)
	as_vapour := co2.phase(.vapor, 320, 1.5e7, [1.0], v)
	assert as_liquid.viscosity.value == as_vapour.viscosity.value
	assert as_liquid.conductivity.value == as_vapour.conductivity.value
	assert as_liquid.notes.any(it.contains("every compound is above its Tc, so the gas's methods (D9)"))
	assert as_vapour.surface_tension.why == 'a gas has no surface tension'
	// A gas at low pressure is its low-pressure value: Lucas's factor is 1
	// there, Stiel and Thodos' residual 0 at zero density.
	gas := db.transport(['methane', 'ethane'])!
	low := gas.phase(.vapor, 300, 1e-300, [0.7, 0.3], 1e300)
	w := components.wilke([db.find('methane')!, db.find('ethane')!])!
	assert low.viscosity.value == w.viscosity(300, 0, [0.7, 0.3])!
}

// One compound alone is its own value, times its compressed-liquid factor,
// to the bit: no mixing rule's rounding.
fn test_one_compound_is_its_own() {
	db := components.load_default()!
	b := db.find('benzene')!.transport_pure()!
	t, p := 320.0, 5e6
	got := db.transport(['benzene', 'toluene'])!.phase(.liquid, t, p, [1.0, 0.0], 1e-4)
	psat := b.value(.vapour_pressure, t)!.value
	assert got.viscosity.value == b.value(.liquid_viscosity, t)!.value * transport.lucas_liquid(t, p, b.tc,
		b.pc, b.omega, psat)
	assert got.conductivity.value == b.value(.liquid_conductivity, t)!.value * transport.missenard(t, p,
		b.tc, b.pc)
	assert got.surface_tension.value == b.value(.surface_tension, t)!.value
}

// What has no value says why, and what's taken beyond its range is noted:
// squalane's liquid viscosity is equation 118; P-diisopropylbenzene has no
// dipole moment for Lucas's polar term; benzene's liquid at 250 K is below
// its correlations' ranges; a wrong composition or state has nothing.
fn test_reasons_and_notes() {
	db := components.load_default()!
	sq := db.transport(['squalane', 'n-hexane'])!.phase(.liquid, 400, 1e5, [0.5, 0.5], 1e-4)
	assert !sq.viscosity.ok() && sq.viscosity.why.contains('is equation 118'), sq.viscosity.why
	assert sq.conductivity.ok() && sq.surface_tension.ok()
	dipb := db.transport(['100-18-5'])!.phase(.vapor, 800, 1e5, [1.0], 0.066)
	assert !dipb.viscosity.ok() && dipb.viscosity.why.contains('dipole moment'), dipb.viscosity.why
	assert dipb.conductivity.ok()
	cold := db.transport(['benzene'])!.phase(.liquid, 250, 1e5, [1.0], 1e-4)
	assert cold.notes.any(it.contains("Benzene's liquid_viscosity at 250.0 K, outside its fitted range, 278.68 K to 562.0 K")), '${cold.notes}'
	mix := db.transport(['benzene', 'toluene'])!
	for bad in [mix.phase(.liquid, 300, 1e5, [1.0], 1e-4),
		mix.phase(.liquid, -1, 1e5, [0.5, 0.5], 1e-4), mix.phase(.liquid, 300, 1e5, [0.0, 0.0], 1e-4),
		mix.phase(.liquid, 300, 1e5, [-0.5, 1.5], 1e-4)] {
		assert !bad.viscosity.ok() && !bad.conductivity.ok() && !bad.surface_tension.ok()
		assert bad.viscosity.why != ''
	}
}

// A stream of IAPWS-95 water reports M5a's R12-08, R15-11 and R1-76 at its
// phase's density, to the bit (D10), and its Setup's model is water's own;
// a cubic's is its compounds'.
fn test_water_by_iapws() {
	db := components.load_default()!
	setup := db.iapws95('water')!
	m := setup.transport()!
	assert m is water.Transport95
	pr := db.setup(.pr, ['methane', 'water'])!.transport()!
	assert pr is transport.Mixture
	// IF97 water's is its own too, the industrial forms (D10); and IF97 is
	// water's alone.
	assert db.if97('water')!.transport()! is water.TransportIf97
	if s := db.if97('methane') {
		assert false, 'an IF97 fluid of ${s.records[0].name}'
	} else {
		assert err.msg().contains("IF97 is water's alone, not Methane"), err.msg()
	}
	for tp in [[300.0, 1e5], [450.0, 1e6], [700.0, 3e7], [373.0, 1e5], [1200.0, 1e7]] {
		r := flash.pt(setup.fluid.model, tp[0], tp[1], [1.0])!
		for ph in r.phases {
			v := ph.state.molar_volume()
			rho := water.molar_mass / v
			got := m.phase(ph.label, r.t, r.p, [1.0], v)
			assert got.viscosity.value == water.viscosity(r.t, rho)!
			assert got.conductivity.value == water.thermal_conductivity(r.t, rho)!
			if ph.label == .vapor {
				assert got.surface_tension.why == 'a gas has no surface tension'
			} else {
				assert got.surface_tension.value == water.surface_tension(r.t)!
			}
		}
	}
	hot := m.phase(.vapor, 1300, 1e8, [1.0], water.molar_mass / 300)
	assert !hot.viscosity.ok() && hot.viscosity.why != ''
}

struct Tally {
mut:
	phases  int
	values  []int    // by property: viscosity, conductivity, surface tension
	reasons []string // each distinct reason
}

fn (mut t Tally) add(r transport.PhaseTransport, label thermo.PhaseLabel) {
	t.phases++
	for i, x in [r.viscosity, r.conductivity, r.surface_tension] {
		if x.ok() {
			assert x.value > 0 && !math.is_inf(x.value, 1), '${x.value}'
			t.values[i]++
		} else {
			assert x.why != '' && math.is_nan(x.value)
			if x.why !in t.reasons {
				t.reasons << x.why
			}
		}
	}
	if label == .vapor {
		assert !r.surface_tension.ok()
	}
}

fn fractions(n []f64) []f64 {
	mut total := 0.0
	for v in n {
		total += v
	}
	return n.map(it / total)
}

fn seeded_states(n int, seed u64) [][]f64 {
	mut s := seed
	mut out := [][]f64{}
	for _ in 0 .. n {
		s = s * 6364136223846793005 + 1442695040888963407
		u := f64(s >> 11) / f64(u64(1) << 53)
		s = s * 6364136223846793005 + 1442695040888963407
		w := f64(s >> 11) / f64(u64(1) << 53)
		out << [u, w]
	}
	return out
}

// Criterion 4: every phase of seeded PT flashes of a cubic's mixtures, a
// γ–φ fluid's with two liquids, and IAPWS-95 water has its values or
// reasons, a value positive and finite, a reason said; how many have each
// is recorded, numbers that can't fall. Of 272 phases, every viscosity but
// one, a dense vapour of the methanol plant's gas just above its
// pseudo-critical pressure below its pseudo-critical temperature, where
// Lucas's form doesn't hold (D5); every conductivity; and 126 surface
// tensions, every liquid's but those with carbon dioxide between 302.6 K,
// where its correlation goes negative, and its Tc of 304.2 K (D3), and
// water labelled liquid above its Tc.
fn test_seeded_flashes() {
	db := components.load_default()!
	mut tally := Tally{
		values: []int{len: 3}
	}
	cubics := [
		['methane', 'ethane', 'propane', 'n-butane', 'nitrogen', 'carbon dioxide'],
		['benzene', 'toluene'],
		['methanol', 'water', 'hydrogen', 'carbon dioxide', 'carbon monoxide'],
		['n-hexane', 'water'],
	]
	feeds := [
		[0.8, 0.07, 0.05, 0.03, 0.03, 0.02],
		[0.4, 0.6],
		[0.3, 0.6, 0.05, 0.03, 0.02],
		[0.5, 0.5],
	]
	for k, names in cubics {
		setup := db.setup(.pr, names)!
		model := setup.transport()!
		for st in seeded_states(40, u64(k + 1)) {
			t := 200.0 + 400.0 * st[0]
			p := math.pow(10, 4 + 3 * st[1])
			r := flash.pt(setup.fluid.model, t, p, feeds[k]) or { continue }
			for ph in r.phases {
				x := fractions(ph.state.n)
				tally.add(model.phase(ph.label, r.t, r.p, x, ph.state.molar_volume()), ph.label)
			}
		}
	}
	gp := db.gamma_phi(['water', '1-butanol'], components.GammaPhiOptions{})!
	gm := components.mixture(gp.records)!
	for t in [290.0, 300.0, 320.0, 340.0, 360.0, 380.0] {
		for z in [0.3, 0.7, 0.95] {
			r := flash.pt(gp.fluid.model, t, 101325, [z, 1 - z]) or { continue }
			for ph in r.phases {
				x := fractions(ph.state.n)
				tally.add(gm.phase(ph.label, r.t, r.p, x, ph.state.molar_volume()), ph.label)
			}
		}
	}
	ws := db.iapws95('water')!
	wm := ws.transport()!
	for st in seeded_states(40, 99) {
		t := 280.0 + 900.0 * st[0]
		p := math.pow(10, 3 + 5 * st[1])
		r := flash.pt(ws.fluid.model, t, p, [1.0]) or { continue }
		for ph in r.phases {
			tally.add(wm.phase(ph.label, r.t, r.p, [1.0], ph.state.molar_volume()), ph.label)
		}
	}
	println('${tally.phases} phases; values ${tally.values}; ${tally.reasons.len} reasons:')
	for why in tally.reasons {
		println('  ${why}')
	}
	assert tally.phases == 272
	for i, want in [271, 272, 126] {
		assert tally.values[i] >= want, '${tally.values}'
	}
	assert tally.reasons.len == 4
}
