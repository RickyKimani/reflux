// bench_m5 times reflux for M5a step 10 (criterion 10), as
// tools/gen_m5_bench.py runs it before and after DWSIM's, CoolProp's and
// the iapws package's runs of the same:
//
//     v -prod -cc gcc -o bench_m5.exe tools/bench_m5.v
//     bench_m5.exe --build          # how it was built
//     bench_m5.exe ROUNDS           # every case, ROUNDS rounds each
//
// The cases are water at three states (300 K and 3 MPa, a liquid; 500 K and
// 1 MPa, a vapour; 700 K and 30 MPa, dense above Tc): IAPWS-95's PT flash
// (flash.pt on its model) and PH flash back from each state's h
// (flash.ph), and IF97's state_pt and state_ph. A round is each case's
// three calls, after one untimed round; it prints one JSON object: each
// round's µs per call, and the answers (h in kJ/kg on IAPWS's reference for
// a PT call, T in K for a PH call), for the test to solve again.
module main

import flash
import json2
import os
import thermo
import time
import water
import water.if97

pub const states = [[300.0, 3e6], [500.0, 1e6], [700.0, 3e7]] // T K, P Pa

struct Case {
	model   string
	kind    string
	us      []f64 // each round's µs per call
	answers []f64 // each call's h (kJ/kg) or T (K)
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

// Iapws95 is what an IAPWS-95 call needs: its phase model, its fluid, and
// D6's constants between its h and IAPWS's.
struct Iapws95 {
	model thermo.PhaseModel
	fluid thermo.Fluid
	dh    f64   // J/mol
	h     []f64 // each state's h, J/mol on reflux's reference, for PH
}

fn iapws95() !Iapws95 {
	fl := water.fluid(none, none)!
	dh, _ := water.reference_offsets(fl)
	m := thermo.by_roots(water.iapws95())
	mut h := []f64{}
	for s in states {
		r := flash.pt(m, s[0], s[1], [1.0])!
		h << fl.props(r.phases[0].state).h
	}
	return Iapws95{m, fl, dh, h}
}

// call is one of a case's calls, its answer.
fn call(w Iapws95, model string, kind string, k int) !f64 {
	t, p := states[k][0], states[k][1]
	if model == 'IAPWS-95' {
		if kind == 'pt' {
			r := flash.pt(w.model, t, p, [1.0])!
			return water.props95(t, 1.0 / r.phases[0].state.v).h / water.molar_mass * 1e-3
		}
		return flash.ph(w.fluid, p, w.h[k], [1.0])!.t
	}
	if kind == 'pt' {
		return if97.state_pt(p, t)!.h * 1e-3
	}
	return if97.state_ph(p, if97.state_pt(p, t)!.h)!.t
}

fn run(w Iapws95, model string, kind string, rounds int) !Case {
	mut answers := []f64{}
	for k in 0 .. 3 {
		answers << call(w, model, kind, k)!
	}
	// IF97's PH call above includes its PT call for h; the timed one takes
	// h ready, as the others' do.
	mut h97 := []f64{}
	for s in states {
		h97 << if97.state_pt(s[1], s[0])!.h
	}
	mut us := []f64{}
	for _ in 0 .. rounds {
		t0 := time.sys_mono_now()
		for k in 0 .. 3 {
			if model == 'IF97' && kind == 'ph' {
				_ := if97.state_ph(states[k][1], h97[k])!
			} else {
				_ := call(w, model, kind, k)!
			}
		}
		us << ms_since(t0) * 1e3 / 3.0
	}
	return Case{model, kind, us, answers}
}

fn main() {
	args := os.args[1..]
	if args.len == 1 && args[0] == '--build' {
		println(json2.encode({
			'compiler': $if gcc { 'gcc' } $else { 'other' }
			'prod':     $if prod { 'yes' } $else { 'no' }
		}))
		return
	}
	rounds := if args.len > 0 { args[0].int() } else { 50 }
	w := iapws95() or { panic(err) }
	mut cases := []Case{}
	for model in ['IAPWS-95', 'IF97'] {
		for kind in ['pt', 'ph'] {
			cases << run(w, model, kind, rounds) or { panic('${model} ${kind}: ${err}') }
		}
	}
	println(json2.encode(cases))
}
