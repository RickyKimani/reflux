// bench_m5b times reflux for M5b step 11 (criterion 10), as
// tools/gen_m5b_bench.py runs it before and after DWSIM's and thermo's runs
// of the same:
//
//     v -prod -cc gcc -o bench_m5b.exe tools/bench_m5b.v
//     bench_m5b.exe --build          # how it was built
//     bench_m5b.exe ROUNDS           # every case, ROUNDS rounds each
//
// - transport: a phase's viscosity, thermal conductivity and surface
//   tension (transport.Mixture.phase, all three at once) on Peng–Robinson's
//   phase at three states: benzene–toluene's liquid at 300 K and 1 atm, the
//   reforming gas at 1000 K and 30 bar, and methane–propane at 320 K and
//   100 bar, a dense gas;
// - water: IF97's PT and PH flashes (flash.pt and flash.ph on its phase
//   model, M5b's D12), IAPWS-95's, and water.if97's state_pt and state_ph,
//   at M5a's three states (300 K and 3 MPa, 500 K and 1 MPa, 700 K and
//   30 MPa);
// - the report: a solve's report (Loaded.report) of the reference
//   flowsheet and of the IF97 Rankine cycle, with each phase's transport
//   and without it (transport: none), D11's cost.
//
// A round is each case's calls, after one untimed round. It prints one JSON
// object: each round's µs per call, and the answers, for the test to solve
// again.
module main

import components
import flash
import flowsheet
import json2
import os
import thermo
import time
import water.if97

// A transport case: compounds by name, mole fractions, T K, P Pa, label.
struct TransportCase {
	name  string
	names []string
	z     []f64
	t     f64
	p     f64
	label string
}

const transport_cases = [
	TransportCase{'benzene-toluene liquid', ['benzene', 'toluene'], [0.4, 0.6], 300.0, 101325.0, 'liquid'},
	TransportCase{'reforming gas', ['methane', 'hydrogen', 'water', 'carbon dioxide', 'carbon monoxide'], [
		0.15,
		0.5,
		0.2,
		0.05,
		0.1,
	], 1000.0, 3e6, 'vapour'},
	TransportCase{'methane-propane dense gas', ['methane', 'propane'], [0.8, 0.2], 320.0, 1e7, 'vapour'},
]

pub const water_states = [[300.0, 3e6], [500.0, 1e6], [700.0, 3e7]] // T K, P Pa

struct Timed {
	name    string
	kind    string
	us      []f64 // each round's µs per call
	answers []f64
}

fn us_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e3
}

fn transport_run(db components.Database, c TransportCase, rounds int) !Timed {
	fl := db.setup(.pr, c.names)!.fluid
	mix := components.mixture(c.names.map(db.find(it)!))!
	label := if c.label == 'vapour' { thermo.PhaseLabel.vapor } else { thermo.PhaseLabel.liquid }
	root := if c.label == 'vapour' { thermo.Root.vapor } else { thermo.Root.liquid }
	st := fl.model.phase(c.t, c.p, c.z, root)!
	v := st.v
	one := mix.phase(label, c.t, c.p, c.z, v)
	mut answers := []f64{}
	for r in [one.viscosity, one.conductivity, one.surface_tension] {
		answers << if r.ok() { r.value } else { -1.0 } // -1 for none: a vapour's surface tension
	}
	mut us := []f64{}
	for _ in 0 .. rounds {
		t0 := time.sys_mono_now()
		_ := mix.phase(label, c.t, c.p, c.z, v)
		us << us_since(t0)
	}
	return Timed{c.name, 'phase', us, answers}
}

// Water is what the water cases need: each formulation's phase model and
// fluid, and each state's h by it, J/mol on reflux's reference.
struct Water {
	m95  thermo.PhaseModel
	f95  thermo.Fluid
	h95  []f64
	m97  thermo.PhaseModel
	f97  thermo.Fluid
	h97  []f64
	hk97 []f64 // IF97's h at each state, J/kg on IAPWS's reference
}

fn water_setup(db components.Database) !Water {
	s95 := db.iapws95('water')!
	s97 := db.if97('water')!
	mut h95, mut h97, mut hk97 := []f64{}, []f64{}, []f64{}
	for s in water_states {
		r95 := flash.pt(s95.fluid.model, s[0], s[1], [1.0])!
		h95 << s95.fluid.props(r95.phases[0].state).h
		r97 := flash.pt(s97.fluid.model, s[0], s[1], [1.0])!
		h97 << s97.fluid.props(r97.phases[0].state).h
		hk97 << if97.state_pt(s[1], s[0])!.h
	}
	return Water{s95.fluid.model, s95.fluid, h95, s97.fluid.model, s97.fluid, h97, hk97}
}

// water_call is one call of a water case: a PT call's molar volume (m³/mol)
// or h (J/kg for IF97's functions), a PH call's T.
fn water_call(w Water, name string, kind string, k int) !f64 {
	t, p := water_states[k][0], water_states[k][1]
	return match name {
		'IF97 flash' {
			if kind == 'pt' {
				flash.pt(w.m97, t, p, [1.0])!.phases[0].state.v
			} else {
				flash.ph(w.f97, p, w.h97[k], [1.0])!.t
			}
		}
		'IAPWS-95 flash' {
			if kind == 'pt' {
				flash.pt(w.m95, t, p, [1.0])!.phases[0].state.v
			} else {
				flash.ph(w.f95, p, w.h95[k], [1.0])!.t
			}
		}
		else {
			if kind == 'pt' { if97.state_pt(p, t)!.h } else { if97.state_ph(p, w.hk97[k])!.t }
		}
	}
}

fn water_run(w Water, name string, kind string, rounds int) !Timed {
	mut answers := []f64{}
	for k in 0 .. water_states.len {
		answers << water_call(w, name, kind, k)!
	}
	mut us := []f64{}
	for _ in 0 .. rounds {
		t0 := time.sys_mono_now()
		for k in 0 .. water_states.len {
			_ := water_call(w, name, kind, k)!
		}
		us << us_since(t0) / f64(water_states.len)
	}
	return Timed{name, kind, us, answers}
}

// report_run times a solved flowsheet's report, with each phase's
// transport and without: its answers are the number of phases with a
// viscosity, and with none.
fn report_run(db components.Database, file string, rounds int) ![]Timed {
	path := os.join_path(@VMODROOT, 'testdata', 'flowsheets', '${file}.toml')
	l := flowsheet.load(path, db: db)!
	sol := l.solve()!
	bare := flowsheet.Loaded{
		...l
		transport: none
	}
	mut out := []Timed{}
	for kind, loaded in {
		'with transport':    l
		'without transport': bare
	} {
		rep := loaded.report(sol)
		mut with, mut without := 0.0, 0.0
		for s in rep.streams {
			for ph in s.phases {
				if _ := ph.viscosity {
					with++
				} else {
					without++
				}
			}
		}
		mut us := []f64{}
		for _ in 0 .. rounds {
			t0 := time.sys_mono_now()
			_ := loaded.report(sol)
			us << us_since(t0)
		}
		out << Timed{file, kind, us, [with, without]}
	}
	return out
}

struct Output {
	transport []Timed
	water     []Timed
	report    []Timed
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
	mut out := Output{}
	mut tr := []Timed{}
	for c in transport_cases {
		tr << transport_run(db, c, rounds * 10) or { panic('${c.name}: ${err}') }
	}
	w := water_setup(db) or { panic(err) }
	mut wt := []Timed{}
	for name in ['IF97 flash', 'IAPWS-95 flash', 'IF97 functions'] {
		for kind in ['pt', 'ph'] {
			wt << water_run(w, name, kind, rounds) or { panic('${name} ${kind}: ${err}') }
		}
	}
	mut rt := []Timed{}
	for file in ['reference', 'rankine_if97'] {
		rt << report_run(db, file, rounds) or { panic('${file}: ${err}') }
	}
	out = Output{tr, wt, rt}
	println(json2.encode(out))
}
