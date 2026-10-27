// Checks the PT flash. The reference splits are Caleb Bell's `thermo`
// 0.6.1 FlashVL results stored in testdata/m3_stability.json (written by
// tools/gen_m3_fixtures.py at PT_SS_TOL 1e-30, or 1e-25 where 1e-30 cycled),
// and every split is checked to be an equilibrium independently (exit
// criterion 3).
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct FixtureComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

struct ThermoResult {
	phases       int
	betas        []f64
	compositions [][]f64
}

struct FlashCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	thermo     ThermoResult
	stable     bool
}

struct Fixture {
	components []FixtureComponent
	cases      []FlashCase
}

fn read_fixture() !Fixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_stability.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py`.') }
	return json2.decode[Fixture](text)!
}

fn model_for(f Fixture, c FlashCase) !thermo.Cubic {
	comps := c.components.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if c.eos == 'pr' { CubicKind.pr } else { CubicKind.srk }
	return thermo.new_cubic(kind, comps, kij: c.kij)!
}

// check_equilibrium checks exit criterion 3 on a two-phase result, from
// scratch: fugacities equal, each component's mass balance, Σx = Σy = 1, and
// neither phase unstable.
fn check_equilibrium(m thermo.Model, r flash.FlashResult, at string) ! {
	assert r.phases.len == 2, at
	vap, liq := r.phases[0], r.phases[1]
	assert vap.label == .vapor && liq.label == .liquid, at
	assert vap.state.molar_volume() > liq.state.molar_volume(), '${at}: the vapour is the less dense'
	assert r.beta > 0.0 && r.beta < 1.0 && vap.fraction == r.beta, at
	assert math.abs(vap.fraction + liq.fraction - 1.0) <= 1e-15, at
	y, x := vap.state.n, liq.state.n
	sy := thermo.solve_tp(m, r.t, r.p, y, .stable)!
	sx := thermo.solve_tp(m, r.t, r.p, x, .stable)!
	mut sum_x, mut sum_y := 0.0, 0.0
	for i, z in r.z {
		sum_x += x[i]
		sum_y += y[i]
		assert math.abs(vap.fraction * y[i] + liq.fraction * x[i] - z) <= 1e-12, '${at}: mass balance of component ${i}'
		if z > 0.0 {
			res := math.abs(math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
			assert res <= 1e-10, '${at}: |ln fV − ln fL| of component ${i} is ${res}'
		}
	}
	assert math.abs(sum_x - 1.0) <= 1e-14 && math.abs(sum_y - 1.0) <= 1e-14, at
	for pair in [[x, y], [y, x]] {
		s := flash.stability(m, r.t, r.p, pair[0], trials: [pair[1]], tol: 1e-8)!
		assert s.stable, '${at}: a phase of the split is unstable'
	}
}

fn test_pt_matches_thermos_splits() {
	// The 166 unstable states of the stability fixture: binaries, states
	// near the methane/n-butane critical locus, a ternary and a
	// 10-component gas condensate. β and both compositions agree with
	// thermo to 1e-8 (exit criterion 2): measured, 4e-9 and 2e-10.
	f := read_fixture()!
	mut splits := 0
	for c in f.cases {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		m := model_for(f, c)!
		r := flash.pt(m, c.t, c.p, c.z)!
		if c.stable {
			assert r.phases.len == 1 && r.path == .single_phase, at
			continue
		}
		assert r.phases.len == 2, at
		// thermo lists its phases gas first, but labels them by PIP, so match
		// reflux's vapour to the nearer of its two.
		vap := r.phases[0].state.n
		mut d := [0.0, 0.0]
		for k in 0 .. 2 {
			for i, v in vap {
				d[k] += math.abs(v - c.thermo.compositions[k][i])
			}
		}
		g := if d[0] <= d[1] { 0 } else { 1 }
		assert math.abs(r.beta - c.thermo.betas[g]) <= 1e-8, '${at}: β ${r.beta} vs ${c.thermo.betas[g]}'
		for i in 0 .. vap.len {
			assert math.abs(vap[i] - c.thermo.compositions[g][i]) <= 1e-8, '${at}: y[${i}]'
			assert math.abs(r.phases[1].state.n[i] - c.thermo.compositions[1 - g][i]) <= 1e-8, '${at}: x[${i}]'
		}
		splits++
	}
	assert splits == 166
}

fn test_every_split_is_an_equilibrium() {
	f := read_fixture()!
	mut newton, mut most := 0, 0
	for c in f.cases {
		if c.stable {
			continue
		}
		m := model_for(f, c)!
		r := flash.pt(m, c.t, c.p, c.z)!
		check_equilibrium(m, r, '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa')!
		if r.path == .newton {
			newton++
		}
		most = math.max(most, r.iterations)
	}
	// Near the critical locus successive substitution slows and Newton on G
	// takes over: 44 of the 166 splits. None takes more than 28 steps.
	assert newton > 20
	assert most <= 35, '${most} steps'
}

fn test_the_split_budget() {
	// Recorded, and can't grow by more than 2.4%: 2,246 split steps over
	// the fixture's 166 splits, 44 of them finished by Newton on G. Newton
	// rescues a poor start or a slow iteration, so mistakes there show up
	// as work, not wrong answers: without DEM it's 2,697 steps and 70
	// finishes, and with K = w/z instead of W/z it's 2,372 and 49.
	f := read_fixture()!
	mut total, mut newton := 0, 0
	for c in f.cases.filter(!it.stable) {
		m := model_for(f, c)!
		r := flash.pt(m, c.t, c.p, c.z, verify: false)!
		total += r.iterations
		if r.path == .newton {
			newton++
		}
	}
	assert total <= 2300, '${total} split steps'
	assert newton <= 46, '${newton} Newton finishes'
}

fn test_a_single_phase_is_labelled_by_its_pip() {
	f := read_fixture()!
	for c in f.cases {
		if !c.stable {
			continue
		}
		m := model_for(f, c)!
		r := flash.pt(m, c.t, c.p, c.z)!
		assert r.phases.len == 1 && r.phases[0].fraction == 1.0
		ph := r.phases[0]
		assert ph.label == thermo.phase_label(m, ph.state)
		assert r.beta == if ph.label == .vapor { 1.0 } else { 0.0 }
		for i, z in r.z {
			assert ph.state.n[i] == z
		}
	}
}

fn condensate() !thermo.Cubic {
	f := read_fixture()!
	c := f.cases.filter(it.components.len == 10)[0]
	return model_for(f, c)!
}

const condensate_z = [0.70, 0.08, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01]

fn test_retrograde_condensation() {
	// At 400 K the gas condensate is below its cricondentherm and above its
	// critical temperature. Compressing it condenses liquid and then, past a
	// maximum, vaporises it again until the dew point: retrograde
	// behaviour, with every split an equilibrium on the way.
	m := condensate()!
	mut liquid := []f64{}
	for pbar in [20.0, 50.0, 80.0, 110.0, 140.0, 155.0] {
		r := flash.pt(m, 400.0, pbar * 1e5, condensate_z)!
		check_equilibrium(m, r, 'condensate at 400 K, ${pbar} bar')!
		liquid << 1.0 - r.beta
	}
	assert liquid[0] < liquid[1] && liquid[1] < liquid[2] && liquid[2] < liquid[3], '${liquid}'
	assert liquid[3] > liquid[4] && liquid[4] > liquid[5], '${liquid}'
	// And at 170 bar it is a single phase again, above its dew point.
	assert flash.pt(m, 400.0, 170e5, condensate_z)!.phases.len == 1
}

fn test_light_and_heavy_side_by_side() {
	// Nitrogen at Tr ≈ 2 beside n-decane at Tr ≈ 0.4: the decane's K is
	// about 1e-5 at 1 bar, and the nitrogen dissolves only as a trace. The
	// split is still an equilibrium to criterion 3, traces included.
	n2 := Component{
		name:  'nitrogen'
		tc:    126.2
		pc:    34.00e5
		omega: 0.038
		mw:    0.028014
	}
	decane := Component{
		name:  'n-decane'
		tc:    617.7
		pc:    21.10e5
		omega: 0.492
		mw:    0.142285
	}
	m := thermo.new_cubic(.pr, [n2, decane], kij: [0.0, 0.11, 0.11, 0.0])!
	for c in [[250.0, 1e5], [300.0, 10e5], [300.0, 100e5], [450.0, 200e5]] {
		r := flash.pt(m, c[0], c[1], [0.5, 0.5])!
		at := 'N2/n-decane at ${c[0]} K, ${c[1]} Pa'
		check_equilibrium(m, r, at)!
		assert r.phases[0].state.n[1] < 0.05, at // little decane in the gas
	}
	// At 1 bar the gas is nearly ideal and the liquid nearly pure decane, so
	// Raoult's law, y = x·Psat/P with PR's own Psat, holds to a few percent
	// for the decane, a trace of 4e-5 in the gas.
	low := flash.pt(m, 250.0, 1e5, [0.5, 0.5])!
	psat := thermo.saturation_pressure(thermo.new_cubic(.pr, [decane])!, 250.0)!
	raoult := low.phases[1].state.n[1] * psat / 1e5
	assert math.abs(low.phases[0].state.n[1] / raoult - 1.0) < 0.05, '${low.phases[0].state.n[1]} vs ${raoult}'
	assert low.phases[1].state.n[0] < 0.01 // a little dissolved nitrogen
}

fn test_three_phases_are_refused_not_hidden() {
	// Water, n-decane and methane at 350 K and 5 bar: PR (with a large
	// water–hydrocarbon kᵢⱼ) has a methane-rich vapour, nearly pure liquid
	// water and a decane-rich liquid. A two-phase flash can't represent
	// that. Re-testing each phase of its split finds the third, so the flash
	// reports more_phases rather than a split that only looks converged.
	water := Component{
		name:  'water'
		tc:    647.1
		pc:    220.55e5
		omega: 0.345
		mw:    0.018015
	}
	decane := Component{
		name:  'n-decane'
		tc:    617.7
		pc:    21.10e5
		omega: 0.492
		mw:    0.142285
	}
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	m := thermo.new_cubic(.pr, [water, decane, methane],
		kij: [0.0, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.04, 0.0]
	)!
	z := [0.5, 0.3, 0.2]
	for c in [[300.0, 10e5], [350.0, 5e5], [400.0, 20e5]] {
		if r := flash.pt(m, c[0], c[1], z) {
			assert false, '${c}: a two-phase answer for a three-phase state (β ${r.beta})'
		} else {
			assert err is flash.FlashError
			assert (err as flash.FlashError).kind == .more_phases, err.msg()
		}
	}
	// Without the re-test the flash returns that split, whose fugacities
	// are equal (it's a stationary point of G) but which isn't the
	// equilibrium: which is exactly what verify is for.
	r := flash.pt(m, 350.0, 5e5, z, verify: false)!
	assert r.phases.len == 2
}

fn test_a_split_meets_the_tolerance_it_was_given() {
	// tol = 1e-16 is below what rounding in ln f allows, so neither
	// substitution nor Newton on G can be sure to reach it; Newton stops
	// where G itself stops showing progress. Whatever comes back must meet
	// it or be an error: never a split claimed at a tolerance it didn't
	// reach. The residual recomputed here carries its own rounding, a few
	// ulps of its O(1) terms, so it's held to 1e-15.
	f := read_fixture()!
	mut met, mut refused := 0, 0
	for c in f.cases.filter(!it.stable && it.components == [0, 3] && it.t >= 330.0) {
		m := model_for(f, c)!
		r := flash.pt(m, c.t, c.p, c.z, tol: 1e-16, verify: false) or {
			assert (err as flash.FlashError).kind == .no_convergence, err.msg()
			refused++
			continue
		}
		vap, liq := r.phases[0].state, r.phases[1].state
		for i in 0 .. 2 {
			res := math.abs(math.log(vap.n[i]) + vap.ln_phi[i] - math.log(liq.n[i]) - liq.ln_phi[i])
			assert res <= 1e-15, '${c.z} at ${c.t} K, ${c.p} Pa: residual ${res}'
		}
		met++
	}
	assert met + refused > 20
}

fn test_verify_changes_nothing_when_the_split_is_right() {
	f := read_fixture()!
	for c in f.cases.filter(!it.stable)[..40] {
		m := model_for(f, c)!
		a := flash.pt(m, c.t, c.p, c.z)!
		b := flash.pt(m, c.t, c.p, c.z, verify: false)!
		assert a.beta == b.beta && a.phases[0].state.n == b.phases[0].state.n
	}
}

fn test_pt_refuses_bad_input() {
	f := read_fixture()!
	m := model_for(f, f.cases[0])!
	for args in [[250.0, 30e5, 0.4], [250.0, 30e5, 0.4, -0.6], [0.0, 30e5, 0.4, 0.6],
		[250.0, 0.0, 0.4, 0.6]] {
		flash.pt(m, args[0], args[1], args[2..]) or {
			assert (err as flash.FlashError).kind == .bad_input, err.msg()
			continue
		}
		assert false, '${args} should be refused'
	}
	if _ := flash.pt(m, 250.0, 30e5, [0.4, 0.6], tol: 0.0) {
		assert false, 'tol = 0 should be refused'
	}
}
