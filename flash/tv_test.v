// Checks the TV and UV flashes against testdata/m3_thermo.json (Caleb
// Bell's `thermo` 0.6.1, with a shared quadratic ideal-gas Cp): from each
// PT state's V and U, and each pure component's saturation states, the
// flash must come back to the same T, P and split (exit criterion 4).
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct TvComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
	cp    []f64
}

struct TvPhase {
	fraction f64
	z        []f64
}

struct TvCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	phases     []TvPhase
	v          f64
	u          f64
}

struct TvSaturation {
	components []int
	kij        []f64
	vf         f64
	t          f64
	p          f64
	v          f64
	u          f64
}

struct TvThermo {
	components []TvComponent
	pt         []TvCase
	saturation []TvSaturation
}

fn read_thermo() !TvThermo {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_thermo.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py thermo`.')
	}
	return json2.decode[TvThermo](text)!
}

// fluid builds the fixture's model and its ideal-gas Cp, cp[0] + cp[1]·T +
// cp[2]·T², with no formation data, so that U is on thermo's zero.
fn fluid(f TvThermo, eos string, idx []int, kij []f64) !thermo.Fluid {
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

// own is reflux's own molar V and U at a PT state, from its PT flash. A
// liquid's P is far more sensitive to V than V is to P (by its bulk modulus
// over P, a thousand here and there), so a round trip starts from reflux's
// V, not thermo's, which differs from it by about 1e-11. That V agrees with
// thermo's is thermo_test's job, and thermo's own TV flash, where it ran,
// lands on the same PT state.
fn own(fl thermo.Fluid, c TvCase) !(f64, f64) {
	r := flash.pt(fl.model, c.t, c.p, c.z)!
	mut v, mut u := 0.0, 0.0
	for ph in r.phases {
		v += ph.fraction * ph.state.molar_volume()
		u += ph.fraction * fl.props(ph.state).u
	}
	return v, u
}

fn split_dev(r flash.FlashResult, phases []TvPhase) f64 {
	mut d := 0.0
	for k, ph in r.phases {
		d = math.max(d, math.abs(ph.fraction - phases[k].fraction))
		for i, x in ph.state.n {
			d = math.max(d, math.abs(x - phases[k].z[i]))
		}
	}
	return d
}

fn test_tv_round_trips_from_pt_states() {
	// Exit criterion 4: from each of the 83 PT states' V, at its T, the TV
	// flash comes back to P and the split, to 1e-8 (measured 2e-12 and
	// 2e-11), in at most 450 PT flashes and Helmholtz steps in all (436).
	f := read_thermo()!
	mut worst_p, mut worst_split := 0.0, 0.0
	mut at_p := ''
	mut steps := 0
	for c in f.pt {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		fl := fluid(f, c.eos, c.components, c.kij)!
		v, _ := own(fl, c)!
		r := flash.tv(fl.model, c.t, v, c.z) or { panic('${at}: tv: ${err}') }
		assert r.phases.len == c.phases.len, '${at}: ${r.phases.len} phases at ${r.p} Pa'
		dp := math.abs(r.p - c.p) / c.p
		if dp > worst_p {
			worst_p, at_p = dp, at
		}
		worst_split = math.max(worst_split, split_dev(r, c.phases))
		steps += r.iterations
	}
	println('TV round trips: P ${worst_p:.2g} (${at_p}), split ${worst_split:.2g}; ${steps} steps')
	assert worst_p <= 1e-8 && worst_split <= 1e-8
	assert steps <= 450
}

fn test_uv_round_trips_from_pt_states() {
	// From each PT state's U and V, the UV flash comes back to T, P and the
	// split, to 1e-8 (measured 4.6e-12, 3.9e-11 and 1.2e-11), in at most
	// 980 TV flashes in all (958; 1122 without the split's Newton finish).
	f := read_thermo()!
	mut worst_t, mut worst_p, mut worst_split := 0.0, 0.0, 0.0
	mut at_t := ''
	mut steps := 0
	for c in f.pt {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		fl := fluid(f, c.eos, c.components, c.kij)!
		v, u := own(fl, c)!
		r := flash.uv(fl, u, v, c.z) or { panic('${at}: uv: ${err}') }
		assert r.phases.len == c.phases.len, '${at}: ${r.phases.len} phases at ${r.t} K, ${r.p} Pa'
		dt := math.abs(r.t - c.t) / c.t
		if dt > worst_t {
			worst_t, at_t = dt, at
		}
		worst_p = math.max(worst_p, math.abs(r.p - c.p) / c.p)
		worst_split = math.max(worst_split, split_dev(r, c.phases))
		steps += r.iterations
	}
	println('UV round trips: T ${worst_t:.2g} (${at_t}), P ${worst_p:.2g}, split ${worst_split:.2g}; ${steps} steps')
	assert worst_t <= 1e-8 && worst_p <= 1e-8 && worst_split <= 1e-8
	assert steps <= 980
}

fn test_pure_saturation() {
	// A pure component's V jumps across its saturation dome. From inside it,
	// TV gives back Psat and the vapour fraction, and UV Tsat, Psat and the
	// vapour fraction: P and T to 1e-9, the vapour fraction to 1e-8
	// (measured 1.2e-10, 1.5e-11 and 1.5e-10). Inside the dome P is Psat(T)
	// whatever V, so thermo's V and U serve here.
	f := read_thermo()!
	mut worst_p, mut worst_t, mut worst_b := 0.0, 0.0, 0.0
	for c in f.saturation {
		at := '${c.components} at ${c.t} K, ${c.p} Pa, vf ${c.vf}'
		fl := fluid(f, 'pr', c.components, c.kij)!
		rt := flash.tv(fl.model, c.t, c.v, [1.0]) or { panic('${at}: tv: ${err}') }
		ru := flash.uv(fl, c.u, c.v, [1.0]) or { panic('${at}: uv: ${err}') }
		worst_p = math.max(worst_p, math.max(math.abs(rt.p - c.p), math.abs(ru.p - c.p)) / c.p)
		worst_t = math.max(worst_t, math.abs(ru.t - c.t) / c.t)
		worst_b = math.max(worst_b, math.max(math.abs(rt.beta - c.vf), math.abs(ru.beta - c.vf)))
	}
	println('pure saturation: P ${worst_p:.2g}, T ${worst_t:.2g}, beta ${worst_b:.2g}')
	assert worst_p <= 1e-9 && worst_t <= 1e-9 && worst_b <= 1e-8
}

fn test_three_phases_are_refused() {
	// pt_test's water, n-decane and methane at 350 K and 5 bar has three
	// phases. Given the V of the false two-phase split there, tv's PT
	// flashes near 5 bar re-test their phases and report more_phases.
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
	z := [0.5, 0.3, 0.2]
	r := flash.pt(m, 350.0, 5e5, z, verify: false)!
	mut v := 0.0
	for ph in r.phases {
		v += ph.fraction * ph.state.molar_volume()
	}
	if _ := flash.tv(m, 350.0, v, z) {
		assert false, 'a two-phase answer for a three-phase state'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
	}
}

fn test_the_tolerance_contract() {
	// At any vtol and etol the answer meets them: at 1e-4, where the
	// searches stop early, on every fourth PT state; and at 1e-16, below
	// what V and U can be computed to, every answer meets it or is an
	// error, on every eighth.
	f := read_thermo()!
	mut met, mut refused := 0, 0
	for k, c in f.pt {
		if k % 4 != 0 {
			continue
		}
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		fl := fluid(f, c.eos, c.components, c.kij)!
		v, u := own(fl, c)!
		for tol in [1e-4, 1e-16] {
			if tol == 1e-16 && k % 8 != 0 {
				continue
			}
			rt := flash.tv(fl.model, c.t, v, c.z, vtol: tol) or {
				assert tol == 1e-16 && err is flash.FlashError, '${at}: ${err}'
				refused++
				continue
			}
			mut got := 0.0
			for ph in rt.phases {
				got += ph.fraction * ph.state.molar_volume()
			}
			assert math.abs(got - v) <= tol * v, '${at}: TV, vtol ${tol}'
			ru := flash.uv(fl, u, v, c.z, vtol: math.max(tol, 1e-12), etol: tol) or {
				assert tol == 1e-16 && err is flash.FlashError, '${at}: ${err}'
				refused++
				continue
			}
			mut gu := 0.0
			for ph in ru.phases {
				gu += ph.fraction * fl.props(ph.state).u
			}
			assert math.abs(gu - u) <= tol * thermo.r_gas * ru.t, '${at}: UV, etol ${tol}'
			if tol == 1e-16 {
				met++
			}
		}
	}
	println('at 1e-16: ${met} met, ${refused} refused')
}

fn test_bad_input() {
	f := read_thermo()!
	c := f.pt[0]
	fl := fluid(f, c.eos, c.components, c.kij)!
	bad := {
		'T = 0':                fn [fl, c] () !flash.FlashResult {
			return flash.tv(fl.model, 0.0, c.v, c.z)
		}
		'V < 0':                fn [fl, c] () !flash.FlashResult {
			return flash.tv(fl.model, c.t, -1e-3, c.z)
		}
		'V below the covolume': fn [fl, c] () !flash.FlashResult {
			return flash.tv(fl.model, c.t, 1e-6, c.z)
		}
		'U is NaN':             fn [fl, c] () !flash.FlashResult {
			return flash.uv(fl, math.nan(), c.v, c.z)
		}
		'a negative z':         fn [fl, c] () !flash.FlashResult {
			return flash.tv(fl.model, c.t, c.v, [0.5, -0.5])
		}
		'vtol = 0':             fn [fl, c] () !flash.FlashResult {
			return flash.tv(fl.model, c.t, c.v, c.z, vtol: 0.0)
		}
		'an unreachable U':     fn [fl, c] () !flash.FlashResult {
			return flash.uv(fl, 1e8, c.v, c.z)
		}
	}
	why := {
		'T = 0':                'not a temperature'
		'V < 0':                'not a volume'
		'V below the covolume': 'smallest volume'
		'U is NaN':             'not a number'
		'a negative z':         'not a mole fraction'
		'vtol = 0':             'must be positive'
		'an unreachable U':     'outside the model'
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
