// Checks the vapour-fraction flashes (tvf, pvf) and bubble and dew points,
// against testdata/m3_thermo.json (Caleb Bell's `thermo` 0.6.1) and against
// the PT flash: from every PT split, the vapour-fraction flash at the same β
// must come back to the same state (exit criterion 4). DWSIM's bubble and
// dew points are checked in dwsim_test.v.
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct VfComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

struct VfPhase {
	fraction f64
	z        []f64
}

struct VfPtCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	phases     []VfPhase
}

struct VfSaturation {
	components []int
	kij        []f64
	spec       string
	vf         f64
	t          f64
	p          f64
}

struct VfBoundary {
	eos        string
	components []int
	kij        []f64
	z          []f64
	spec       string
	vf         f64
	t          f64
	p          f64
	incipient  []f64
	error      string
	rejected   string
}

struct VfThermo {
	components []VfComponent
	pt         []VfPtCase
	saturation []VfSaturation
	boundaries []VfBoundary
}

fn read_thermo() !VfThermo {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_thermo.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py thermo`.')
	}
	return json2.decode[VfThermo](text)!
}

fn cubic(f VfThermo, eos string, idx []int, kij []f64) !thermo.Cubic {
	cs := idx.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if eos == 'srk' { CubicKind.srk } else { CubicKind.pr }
	return thermo.new_cubic(kind, cs, kij: kij)!
}

// model_of is the fixture's model for the first PT case with these
// components, so that tests use its kᵢⱼ.
fn model_of(f VfThermo, idx []int) !thermo.Cubic {
	for c in f.pt {
		if c.eos == 'pr' && c.components == idx {
			return cubic(f, 'pr', idx, c.kij)!
		}
	}
	return error('no PT case for ${idx}')
}

const methane_butane = [0, 3]
const condensate = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
const condensate_z = [0.70, 0.08, 0.05, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02, 0.01]

// density is a phase's mass density, kg/m³: D6's measure of which phase is
// the vapour.
fn density(m thermo.Model, st thermo.State) f64 {
	comps := m.components()
	mut mass := 0.0
	for i, n in st.n {
		mass += n * comps[i].mw
	}
	return mass / st.v
}

// equilibrium checks exit criterion 3 on a vapour-fraction result, from
// scratch: equal fugacities, each component's mass balance, Σx = Σy = 1,
// the vapour the less dense, and neither phase unstable.
fn equilibrium(m thermo.Model, r flash.FlashResult, at string) ! {
	assert r.phases.len == 2, at
	vap, liq := r.phases[0], r.phases[1]
	assert vap.label == .vapor && liq.label == .liquid, at
	assert density(m, vap.state) < density(m, liq.state), '${at}: the vapour is the less dense'
	assert vap.fraction == r.beta && liq.fraction == 1.0 - r.beta, at
	y, x := vap.state.n, liq.state.n
	sy := thermo.solve_tp(m, r.t, r.p, y, .stable)!
	sx := thermo.solve_tp(m, r.t, r.p, x, .stable)!
	mut sum_x, mut sum_y := 0.0, 0.0
	for i, z in r.z {
		sum_x += x[i]
		sum_y += y[i]
		assert math.abs(r.beta * y[i] + (1.0 - r.beta) * x[i] - z) <= 1e-12, '${at}: mass balance of component ${i}'
		if z > 0.0 {
			res := math.abs(math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
			assert res <= 1e-10, '${at}: |ln fV − ln fL| of component ${i} is ${res}'
		}
	}
	assert math.abs(sum_x - 1.0) <= 1e-14 && math.abs(sum_y - 1.0) <= 1e-14, '${at}: Σx − 1 = ${sum_x - 1.0}, Σy − 1 = ${sum_y - 1.0}'
	for pair in [[x, y], [y, x]] {
		s := flash.stability(m, r.t, r.p, pair[0], trials: [pair[1]], tol: 1e-8)!
		assert s.stable, '${at}: a phase is unstable'
	}
}

// is_boundary checks that a bubble or dew point is one in the PT flash's
// eyes: 1e-7 to one side of it (relatively, in the unknown T or P) a single
// phase, and to the other a split with a trace of the incipient phase.
fn is_boundary(m thermo.Model, r flash.FlashResult, given_t bool, at string) ! {
	incipient := if r.beta == 0.0 { r.phases[0].state.n } else { r.phases[1].state.n }
	mut sides := []f64{}
	for d in [1e-7, -1e-7] {
		t, p := if given_t { r.t, r.p * (1.0 + d) } else { r.t * (1.0 + d), r.p }
		s := flash.pt(m, t, p, r.z)!
		if s.phases.len == 1 {
			sides << 0.0
			continue
		}
		mut best, mut frac := math.inf(1), 0.0
		for ph in s.phases {
			mut dist := 0.0
			for i, w in ph.state.n {
				dist += math.abs(w - incipient[i])
			}
			if dist < best {
				best, frac = dist, ph.fraction
			}
		}
		sides << frac
	}
	sides.sort()
	assert sides[0] == 0.0 && sides[1] > 0.0 && sides[1] < 1e-3, '${at}: incipient fractions either side ${sides}'
}

fn test_bubble_and_dew_points_match_thermo() {
	// thermo's 40 bubble and dew points, polished to a fugacity residual of
	// 1e-13: T or P to 1e-11 relative and the incipient phase to 1e-11
	// (measured 9e-13 and 1.3e-13), every one an equilibrium, in at most
	// 170 Newton steps in all (163 measured with the line search; the dogleg
	// takes 182).
	f := read_thermo()!
	mut worst_tp, mut worst_w := 0.0, 0.0
	mut runs, mut steps := 0, 0
	for c in f.boundaries {
		if c.error != '' || c.rejected != '' {
			continue
		}
		at := '${c.eos} ${c.components} z ${c.z}, ${c.spec} given, vf ${c.vf}'
		m := cubic(f, c.eos, c.components, c.kij)!
		r := if c.spec == 'p' {
			flash.pvf(m, c.p, c.vf, c.z)!
		} else {
			flash.tvf(m, c.t, c.vf, c.z)!
		}
		equilibrium(m, r, at)!
		assert r.path == .newton, at
		worst_tp = math.max(worst_tp, if c.spec == 'p' {
			math.abs(r.t - c.t) / c.t
		} else {
			math.abs(r.p - c.p) / c.p
		})
		w := if c.vf == 0.0 { r.phases[0].state.n } else { r.phases[1].state.n }
		for i, wi in w {
			worst_w = math.max(worst_w, math.abs(wi - c.incipient[i]))
		}
		runs++
		steps += r.iterations
	}
	println('thermo boundaries: T or P ${worst_tp:.2g}, incipient phase ${worst_w:.2g}; ${runs} runs, ${steps} steps')
	assert worst_tp <= 1e-11 && worst_w <= 1e-11
	assert runs == 40 && steps <= 170
}

fn test_pure_saturation_matches_thermo() {
	// A pure component's answer is its saturation point whatever β: T or P
	// to 1e-11 of thermo's FlashPureVLS (9.1e-13 measured), the two roots
	// with equal fugacities.
	f := read_thermo()!
	mut worst := 0.0
	for c in f.saturation {
		at := '${c.components} ${c.spec} given, vf ${c.vf}'
		m := cubic(f, 'pr', c.components, c.kij)!
		r := if c.spec == 'p' {
			flash.pvf(m, c.p, c.vf, [1.0])!
		} else {
			flash.tvf(m, c.t, c.vf, [1.0])!
		}
		worst = math.max(worst, if c.spec == 'p' {
			math.abs(r.t - c.t) / c.t
		} else {
			math.abs(r.p - c.p) / c.p
		})
		assert r.beta == c.vf && r.phases.len == 2, at
		assert r.phases[0].state.v > r.phases[1].state.v, '${at}: two distinct roots'
		assert math.abs(r.phases[0].state.ln_phi[0] - r.phases[1].state.ln_phi[0]) <= 1e-10, at
	}
	println('pure saturation: T or P ${worst:.2g}')
	assert worst <= 1e-11
}

fn test_round_trips_from_pt_splits() {
	// Exit criterion 4: from each of thermo's 39 PT splits, tvf at (T, β)
	// gives back P, and pvf at (P, β) gives back T, with the same split: to
	// 1e-8 relative (measured 2.6e-12 in P, 8.5e-14 in T and 7e-13 in the
	// split).
	//
	// One state has a twin. Methane with n-decane at 200 bar has β = 0.115
	// at 450 K and again at 345.6 K, and pvf from Wilson's start finds the
	// cooler one. That one is checked to be an equilibrium, and pvf started
	// 2% from 450 K finds the original.
	f := read_thermo()!
	mut worst_t, mut worst_p, mut worst_split := 0.0, 0.0, 0.0
	mut splits, mut twins := 0, 0
	for c in f.pt {
		if c.phases.len != 2 {
			continue
		}
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		m := cubic(f, c.eos, c.components, c.kij)!
		beta := c.phases[0].fraction
		rp := flash.tvf(m, c.t, beta, c.z)!
		mut rt := flash.pvf(m, c.p, beta, c.z)!
		equilibrium(m, rp, at)!
		equilibrium(m, rt, at)!
		if math.abs(rt.t - c.t) > 1e-6 * c.t {
			twins++
			rt = flash.pvf(m, c.p, beta, c.z, estimate: 1.02 * c.t)!
			equilibrium(m, rt, at)!
		}
		worst_p = math.max(worst_p, math.abs(rp.p - c.p) / c.p)
		worst_t = math.max(worst_t, math.abs(rt.t - c.t) / c.t)
		for r in [rp, rt] {
			for k in 0 .. 2 {
				for i, x in r.phases[k].state.n {
					worst_split = math.max(worst_split, math.abs(x - c.phases[k].z[i]))
				}
			}
		}
		splits++
	}
	println('round trips: P ${worst_p:.2g}, T ${worst_t:.2g}, split ${worst_split:.2g}')
	assert worst_p <= 1e-8 && worst_t <= 1e-8 && worst_split <= 1e-8
	assert splits == 39 && twins == 1
}

fn test_nothing_above_the_cricondenbar() {
	// Methane/n-butane (40/60) splits at 90 bar from 327 to 379 K, by PT
	// flashes 1 K apart, and nowhere at 100 bar: its cricondenbar is in
	// between. So close to its critical point both ends of the 90-bar
	// isobar are bubble points: Wilson's start finds the cooler, and a start
	// at 375 K the hotter, both boundaries. There's no dew point at 90 bar,
	// and dew_t says it found none from a split it knows is there. At
	// 100 bar neither exists, and both say the point is outside the
	// envelope.
	f := read_thermo()!
	m := model_of(f, methane_butane)!
	z := [0.4, 0.6]
	cool := flash.bubble_t(m, 90e5, z)!
	hot := flash.bubble_t(m, 90e5, z, estimate: 375.0)!
	for r in [cool, hot] {
		equilibrium(m, r, 'bubble point at 90 bar, ${r.t} K')!
		is_boundary(m, r, false, 'bubble point at 90 bar, ${r.t} K')!
	}
	assert cool.t > 320.0 && cool.t < 330.0 && hot.t > 375.0 && hot.t < 383.0, '${cool.t} and ${hot.t}'
	if _ := flash.dew_t(m, 90e5, z) {
		assert false, 'a dew point at 90 bar'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .no_convergence, err.msg()
	}
	for k in 0 .. 61 {
		assert flash.pt(m, 200.0 + 5.0 * f64(k), 100e5, z)!.phases.len == 1
	}
	for beta in [0.0, 1.0] {
		if _ := flash.pvf(m, 100e5, beta, z) {
			assert false, 'a point with β = ${beta} above the cricondenbar'
		} else {
			assert err is flash.FlashError && (err as flash.FlashError).kind == .outside_envelope, err.msg()
		}
	}
}

fn test_no_dew_point_past_the_cricondentherm() {
	// At 393 K the same mixture splits from 58 to 77.5 bar, and both ends
	// are dew points: Wilson's start finds the lower, and a start at 77 bar,
	// inside the split, the upper. At 400 K, past its cricondentherm, no PT
	// flash splits and dew_p says so.
	f := read_thermo()!
	m := model_of(f, methane_butane)!
	z := [0.4, 0.6]
	lower := flash.dew_p(m, 393.0, z)!
	upper := flash.dew_p(m, 393.0, z, estimate: 77e5)!
	for r in [lower, upper] {
		equilibrium(m, r, 'dew point at 393 K, ${r.p} Pa')!
		is_boundary(m, r, true, 'dew point at 393 K, ${r.p} Pa')!
	}
	assert lower.p > 57e5 && lower.p < 59e5 && upper.p > 77e5 && upper.p < 78e5, '${lower.p} and ${upper.p}'
	for k in 0 .. 41 {
		assert flash.pt(m, 400.0, 1e5 * math.pow(10.0, 0.06 * f64(k)), z)!.phases.len == 1
	}
	if _ := flash.dew_p(m, 400.0, z) {
		assert false, 'a dew point past the cricondentherm'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .outside_envelope, err.msg()
	}
}

fn test_two_dew_pressures_of_a_gas_condensate() {
	// The 10-component condensate at 400 K condenses as P rises from its
	// lower dew point and vaporises again up to its upper one. Wilson's
	// start finds the lower, thermo's (13.26 bar, in m3_thermo.json), and a
	// start at 200 bar finds the upper; both are boundaries.
	f := read_thermo()!
	m := model_of(f, condensate)!
	lower := flash.dew_p(m, 400.0, condensate_z)!
	upper := flash.dew_p(m, 400.0, condensate_z, estimate: 200e5)!
	for r in [lower, upper] {
		equilibrium(m, r, 'condensate dew point at 400 K, ${r.p} Pa')!
		is_boundary(m, r, true, 'condensate dew point at 400 K, ${r.p} Pa')!
	}
	mut want := 0.0
	for c in f.boundaries {
		if c.components == condensate && c.spec == 't' && c.vf == 1.0 && c.t == 400.0 {
			want = c.p
		}
	}
	println('condensate dew pressures at 400 K: ${lower.p} and ${upper.p} Pa')
	assert math.abs(lower.p - want) <= 1e-10 * want
	assert upper.p > 100e5 && upper.p < 170e5, '${upper.p}'
}

fn test_three_phases_are_refused() {
	// pt_test's water, n-decane and methane at 350 K and 5 bar has three
	// phases, and a two-phase split there only looks converged. Asked for
	// that split's β at 350 K, tvf finds 5 bar again, and its re-test of
	// the phases reports more_phases; without the re-test it returns the
	// false split.
	m := thermo.new_cubic(.pr, [
		Component{
			name:  'water'
			tc:    647.1
			pc:    220.55e5
			omega: 0.345
			mw:    0.018015
		},
		Component{
			name:  'n-decane'
			tc:    617.7
			pc:    21.10e5
			omega: 0.492
			mw:    0.142285
		},
		Component{
			name:  'methane'
			tc:    190.6
			pc:    45.99e5
			omega: 0.012
			mw:    0.016043
		},
	], kij: [0.0, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.04, 0.0])!
	z := [0.5, 0.3, 0.2]
	false_split := flash.pt(m, 350.0, 5e5, z, verify: false)!
	assert false_split.phases.len == 2
	if _ := flash.tvf(m, 350.0, false_split.beta, z) {
		assert false, 'a two-phase answer for a three-phase state'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
	}
	r := flash.tvf(m, 350.0, false_split.beta, z, verify: false)!
	assert math.abs(r.p - 5e5) <= 1e-8 * 5e5, '${r.p}'
}

fn test_a_pure_component_from_a_poor_start() {
	// Started at ten times its saturation pressure, where n-butane has only
	// a liquid root, or a tenth of it, where it has only a vapour root, the
	// pure-component search still reaches both roots and the same answer;
	// likewise from 100 K too hot and too cold.
	f := read_thermo()!
	m := cubic(f, 'pr', [3], [0.0])!
	psat := thermo.saturation_pressure(m, 300.0)!
	for est in [10.0 * psat, 0.1 * psat] {
		r := flash.tvf(m, 300.0, 0.5, [1.0], estimate: est)!
		assert math.abs(r.p - psat) <= 1e-9 * psat, 'from ${est} Pa: ${r.p}'
	}
	tsat := flash.pvf(m, psat, 0.5, [1.0])!.t
	assert math.abs(tsat - 300.0) <= 1e-9 * 300.0
	for est in [200.0, 400.0] {
		r := flash.pvf(m, psat, 0.5, [1.0], estimate: est)!
		assert math.abs(r.t - 300.0) <= 1e-9 * 300.0, 'from ${est} K: ${r.t}'
	}
}

fn test_the_tolerance_contract() {
	// At any tol the answer meets it: fugacities and each component's mass
	// balance to tol. (An earlier version weighted Rachford–Rice's equation
	// to promise tol/100 in the mass balance, but Newton's last step always
	// overshot the promise, even at 1e-4, so no answer could tell.)
	f := read_thermo()!
	for tol in [1e-4, 1e-8] {
		for c in f.boundaries {
			if c.error != '' || c.rejected != '' {
				continue
			}
			at := '${c.eos} ${c.components} z ${c.z}, ${c.spec} given, vf ${c.vf}, tol ${tol}'
			m := cubic(f, c.eos, c.components, c.kij)!
			r := if c.spec == 'p' {
				flash.pvf(m, c.p, c.vf, c.z, tol: tol)!
			} else {
				flash.tvf(m, c.t, c.vf, c.z, tol: tol)!
			}
			y, x := r.phases[0].state.n, r.phases[1].state.n
			for i, z in r.z {
				assert math.abs(r.beta * y[i] + (1.0 - r.beta) * x[i] - z) <= tol, '${at}: mass balance of ${i}'
				res := math.abs(math.log(y[i]) + r.phases[0].state.ln_phi[i] - math.log(x[i]) - r.phases[1].state.ln_phi[i])
				assert res <= tol, '${at}: fugacity of ${i}'
			}
		}
	}
}

fn test_bad_input() {
	// V can't keep results in an array, so each call is a closure.
	f := read_thermo()!
	m := model_of(f, methane_butane)!
	z := [0.4, 0.6]
	bad := {
		'β < 0':            fn [m, z] () !flash.FlashResult {
			return flash.tvf(m, 250.0, -0.1, z)
		}
		'β > 1':            fn [m, z] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 1.5, z)
		}
		'T = 0':            fn [m, z] () !flash.FlashResult {
			return flash.tvf(m, 0.0, 0.5, z)
		}
		'P < 0':            fn [m, z] () !flash.FlashResult {
			return flash.pvf(m, -1.0, 0.5, z)
		}
		'a negative z':     fn [m] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 0.5, [0.4, -0.6])
		}
		'too few z':        fn [m] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 0.5, [1.0])
		}
		'no moles':         fn [m] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 0.5, [0.0, 0.0])
		}
		'a negative start': fn [m, z] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 0.5, z, estimate: -1.0)
		}
		'tol = 0':          fn [m, z] () !flash.FlashResult {
			return flash.tvf(m, 250.0, 0.5, z, tol: 0.0)
		}
	}
	// Each must be refused for its own reason, not some later one.
	why := {
		'β < 0':            'outside [0, 1]'
		'β > 1':            'outside [0, 1]'
		'T = 0':            'not a temperature'
		'P < 0':            'not a pressure'
		'a negative z':     'not a mole fraction'
		'too few z':        'entries'
		'no moles':         'no moles'
		'a negative start': 'estimate'
		'tol = 0':          'tol and max_steps'
	}
	for name, call in bad {
		call() or {
			assert err is flash.FlashError && (err as flash.FlashError).kind == .bad_input, '${name}: ${err}'
			assert err.msg().contains(why[name]), '${name}: ${err}'
			continue
		}
		assert false, '${name} was accepted'
	}
}
