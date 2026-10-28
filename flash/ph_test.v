// Checks the PH and PS flashes against testdata/m3_thermo.json (Caleb
// Bell's `thermo` 0.6.1, with a shared quadratic ideal-gas Cp): from each
// PT state's H and S, and each pure component's saturation states, the
// flash must come back to the same T and split (exit criterion 4). DWSIM's
// PH and PS cases are in dwsim_test.v.
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct PhComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
	cp    []f64
}

struct PhPhase {
	fraction f64
	z        []f64
}

struct PhCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	phases     []PhPhase
	h          f64
	s          f64
}

struct PhSaturation {
	components []int
	kij        []f64
	vf         f64
	t          f64
	p          f64
	h          f64
	s          f64
}

struct PhThermo {
	p_ref_s    f64
	components []PhComponent
	pt         []PhCase
	saturation []PhSaturation
}

fn read_thermo() !PhThermo {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_thermo.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py thermo`.')
	}
	return json2.decode[PhThermo](text)!
}

// fluid builds the fixture's model and its ideal-gas Cp, cp[0] + cp[1]·T +
// cp[2]·T², with no formation data, so that H and S are on thermo's zero.
fn fluid(f PhThermo, eos string, idx []int, kij []f64) !thermo.Fluid {
	cs := idx.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if eos == 'srk' { CubicKind.srk } else { CubicKind.pr }
	m := thermo.new_cubic(kind, cs, kij: kij)!
	return thermo.new_fluid(m, idx.map(thermo.IdealGas{
		cp: thermo.Polynomial{
			terms: [
				thermo.PowerTerm{f.components[it].cp[0], 0},
				thermo.PowerTerm{f.components[it].cp[1], 1},
				thermo.PowerTerm{f.components[it].cp[2], 2},
			]
		}
	}))!
}

// reflux_s turns thermo's S, relative to the ideal gas at 101325 Pa, to
// reflux's, relative to it at 1 bar.
fn reflux_s(f PhThermo, s f64) f64 {
	return s + thermo.r_gas * math.log(thermo.p_ref / f.p_ref_s)
}

// split_dev is the largest difference in the phase fractions and
// compositions, phases in the same order (the vapour first).
fn split_dev(r flash.FlashResult, phases []PhPhase) f64 {
	mut d := 0.0
	for k, ph in r.phases {
		d = math.max(d, math.abs(ph.fraction - phases[k].fraction))
		for i, x in ph.state.n {
			d = math.max(d, math.abs(x - phases[k].z[i]))
		}
	}
	return d
}

fn test_round_trips_from_pt_states() {
	// Exit criterion 4: from each of the 83 PT states' H and S, the PH and
	// PS flashes come back to T and the split, to 1e-8 (measured 2.2e-11 in
	// T and 2.6e-11 in the split), in at most 1990 PT flashes and Newton
	// steps in all (1949 measured; 2018 without the Q-function Newton).
	// Among them are n-butane with isobutane
	// at 5 bar, whose two-phase range is under 4 K wide (under 0.4 K at 5%
	// isobutane): its 6 two-phase states come back too.
	f := read_thermo()!
	mut worst_t, mut worst_split := 0.0, 0.0
	mut at_t := ''
	mut steps, mut narrow := 0, 0
	for c in f.pt {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		fl := fluid(f, c.eos, c.components, c.kij)!
		rh := flash.ph(fl, c.p, c.h, c.z) or { panic('${at}: ph: ${err}') }
		rs := flash.ps(fl, c.p, reflux_s(f, c.s), c.z) or { panic('${at}: ps: ${err}') }
		for r in [rh, rs] {
			assert r.phases.len == c.phases.len, '${at}: ${r.phases.len} phases at ${r.t} K'
			dt := math.abs(r.t - c.t) / c.t
			if dt > worst_t {
				worst_t, at_t = dt, at
			}
			worst_split = math.max(worst_split, split_dev(r, c.phases))
			steps += r.iterations
		}
		if c.components == [3, 10] && c.phases.len == 2 {
			narrow++
		}
	}
	println('round trips: T ${worst_t:.2g} (${at_t}), split ${worst_split:.2g}; ${steps} steps')
	assert worst_t <= 1e-8 && worst_split <= 1e-8
	assert steps <= 1990 && narrow == 6
}

fn test_pure_saturation() {
	// A pure component's H jumps at Tsat. An H or S from inside the jump,
	// the saturated mixture at vapour fraction vf, must give back Tsat to
	// 1e-9 and that vf to 1e-8 (measured 1.5e-11 and 8e-11), the saturated
	// liquid and vapour (vf 0 and 1) included.
	f := read_thermo()!
	mut worst_t, mut worst_b := 0.0, 0.0
	for c in f.saturation {
		at := '${c.components} at ${c.p} Pa, vf ${c.vf}'
		fl := fluid(f, 'pr', c.components, c.kij)!
		rh := flash.ph(fl, c.p, c.h, [1.0]) or { panic('${at}: ph: ${err}') }
		rs := flash.ps(fl, c.p, reflux_s(f, c.s), [1.0]) or { panic('${at}: ps: ${err}') }
		for r in [rh, rs] {
			worst_t = math.max(worst_t, math.abs(r.t - c.t) / c.t)
			worst_b = math.max(worst_b, math.abs(r.beta - c.vf))
		}
	}
	println('pure saturation: T ${worst_t:.2g}, beta ${worst_b:.2g}')
	assert worst_t <= 1e-9 && worst_b <= 1e-8
}

fn test_three_phases_are_refused() {
	// pt_test's water, n-decane and methane at 350 K and 5 bar has three
	// phases. Given the H of the false two-phase split there, ph's PT
	// flashes near 350 K re-test their phases and report more_phases.
	cs := [
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
	]
	m := thermo.new_cubic(.pr, cs, kij: [0.0, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.04, 0.0])!
	fl := thermo.new_fluid(m, [4.0, 25.0, 4.3].map(thermo.IdealGas{
		cp: thermo.Polynomial{
			terms: [thermo.PowerTerm{thermo.r_gas * it, 0}]
		}
	}))!
	z := [0.5, 0.3, 0.2]
	r := flash.pt(m, 350.0, 5e5, z, verify: false)!
	mut h := 0.0
	for ph in r.phases {
		h += ph.fraction * fl.props(ph.state).h
	}
	if _ := flash.ph(fl, 5e5, h, z) {
		assert false, 'a two-phase answer for a three-phase state'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
	}
}

fn test_the_tolerance_contract() {
	// At any etol the answer meets it, |H − H_spec| ≤ etol·RT and
	// |S − S_spec| ≤ etol·R: at 1e-4, where Newton and Brent stop early, on
	// every other PT state (all of them take 50 s under tcc).
	f := read_thermo()!
	for etol in [1e-4] {
		for k, c in f.pt {
			if k % 2 == 1 {
				continue
			}
			at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa, etol ${etol}'
			fl := fluid(f, c.eos, c.components, c.kij)!
			s_spec := reflux_s(f, c.s)
			rh := flash.ph(fl, c.p, c.h, c.z, etol: etol)!
			rs := flash.ps(fl, c.p, s_spec, c.z, etol: etol)!
			mut h, mut s := 0.0, 0.0
			for ph in rh.phases {
				h += ph.fraction * fl.props(ph.state).h
			}
			for ph in rs.phases {
				s += ph.fraction * fl.props(ph.state).s
			}
			assert math.abs(h - c.h) <= etol * thermo.r_gas * rh.t, '${at}: H'
			assert math.abs(s - s_spec) <= etol * thermo.r_gas, '${at}: S'
		}
	}
}

fn test_an_impossible_tolerance() {
	// At etol 1e-16, below what H and S can be computed to, every answer
	// meets it or is an error: never a result that only looks converged.
	f := read_thermo()!
	mut met, mut refused := 0, 0
	for k, c in f.pt {
		if k % 8 != 0 {
			continue
		}
		fl := fluid(f, c.eos, c.components, c.kij)!
		r := flash.ph(fl, c.p, c.h, c.z, etol: 1e-16) or {
			assert err is flash.FlashError, err.msg()
			refused++
			continue
		}
		mut h := 0.0
		for ph in r.phases {
			h += ph.fraction * fl.props(ph.state).h
		}
		assert math.abs(h - c.h) <= 1e-16 * thermo.r_gas * r.t, '${c.components} at ${c.t} K'
		met++
	}
	println('etol 1e-16: ${met} met, ${refused} refused')
	assert met + refused == 11
}

fn test_an_unreachable_specification() {
	// An H that no temperature up to 10⁴ K reaches, and one below what any
	// temperature down to 10 K gives, are bad input, said so.
	f := read_thermo()!
	c := f.pt[0]
	fl := fluid(f, c.eos, c.components, c.kij)!
	for h in [1e7, -1e7] {
		if _ := flash.ph(fl, c.p, h, c.z) {
			assert false, 'H = ${h} was reached'
		} else {
			assert err is flash.FlashError && (err as flash.FlashError).kind == .bad_input, err.msg()
			assert err.msg().contains('outside the model'), err.msg()
		}
	}
}

fn test_bad_input() {
	f := read_thermo()!
	c := f.pt[0]
	fl := fluid(f, c.eos, c.components, c.kij)!
	bad := {
		'P = 0':            fn [fl, c] () !flash.FlashResult {
			return flash.ph(fl, 0.0, c.h, c.z)
		}
		'H is NaN':         fn [fl, c] () !flash.FlashResult {
			return flash.ph(fl, c.p, math.nan(), c.z)
		}
		'a negative z':     fn [fl, c] () !flash.FlashResult {
			return flash.ps(fl, c.p, c.s, [0.5, -0.5])
		}
		'etol = 0':         fn [fl, c] () !flash.FlashResult {
			return flash.ph(fl, c.p, c.h, c.z, etol: 0.0)
		}
		'a negative start': fn [fl, c] () !flash.FlashResult {
			return flash.ph(fl, c.p, c.h, c.z, estimate: -5.0)
		}
	}
	why := {
		'P = 0':            'not a pressure'
		'H is NaN':         'not a number'
		'a negative z':     'not a mole fraction'
		'etol = 0':         'must be positive'
		'a negative start': 'estimate'
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
