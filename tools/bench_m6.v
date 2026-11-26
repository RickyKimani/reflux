// bench_m6 times reflux for M6 step 13 (D18, criterion 12), as
// tools/gen_m6_bench.py runs it before and after DWSIM's and thermo's runs
// of the same:
//
//     v -prod -cc gcc -o bench_m6.exe tools/bench_m6.v
//     bench_m6.exe --build          # how it was built
//     bench_m6.exe ROUNDS           # every case, ROUNDS rounds each
//
// The cases are methanol and water at 1 atm, by NRTL and UNIFAC as γ–φ
// fluids (the vapour the ideal gas, Poynting on) and by Peng–Robinson: PT
// flashes at three states, two-phase by NRTL, and bubble temperatures of
// three liquids. A round is each case's three calls, after one untimed
// round; it prints one JSON object: each round's µs per call, and the
// answers (β, or the bubble T), for the test to solve again.
module main

import components
import flash
import json2
import os
import thermo
import time

pub const p_atm = 101325.0

pub const pt_states = [[350.0, 0.5], [357.0, 0.3], [346.0, 0.7]] // T, methanol's z

pub const bubble_x = [0.1, 0.5, 0.9] // methanol's x

struct Case {
	model   string
	kind    string
	us      []f64 // each round's µs per call
	answers []f64 // each call's β or T
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

fn fluid(db components.Database, model string) !thermo.Fluid {
	names := ['methanol', 'water']
	return match model {
		'NRTL' { db.gamma_phi(names)!.fluid }
		'UNIFAC' { db.gamma_phi(names, liquid: .unifac)!.fluid }
		else { db.fluid(.pr, names)! }
	}
}

// call is one of a case's calls, its answer.
fn call(fl thermo.Fluid, kind string, k int) !f64 {
	if kind == 'pt' {
		s := pt_states[k]
		return flash.pt(fl.model, s[0], p_atm, [s[1], 1.0 - s[1]])!.beta
	}
	x := bubble_x[k]
	return flash.bubble_t(fl.model, p_atm, [x, 1.0 - x])!.t
}

fn run(db components.Database, model string, kind string, rounds int) !Case {
	fl := fluid(db, model)!
	mut answers := []f64{}
	for k in 0 .. 3 {
		answers << call(fl, kind, k)!
	}
	mut us := []f64{}
	for _ in 0 .. rounds {
		t0 := time.sys_mono_now()
		for k in 0 .. 3 {
			_ := call(fl, kind, k)!
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
	db := components.load_default() or { panic(err) }
	mut cases := []Case{}
	for model in ['NRTL', 'UNIFAC', 'PR'] {
		for kind in ['pt', 'bubble'] {
			cases << run(db, model, kind, rounds) or { panic('${model} ${kind}: ${err}') }
		}
	}
	println(json2.encode(cases))
}
