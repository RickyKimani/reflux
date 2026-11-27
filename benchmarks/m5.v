module benchmarks

import json2
import os
import strings

// M5a's speed (step 10, exit criterion 10): water by IAPWS-95 and IF97 in
// reflux, DWSIM 9, CoolProp and the iapws package in one session.
// tools/gen_m5_bench.py times them all, reflux through tools/bench_m5.v,
// and writes testdata/m5_bench.json; tools/benchmarks.vsh adds it to
// docs/benchmarks.md.

pub const m5_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm5_bench.json')

// The cases, in the report's order: each formulation's PT and PH calls, at
// 300 K and 3 MPa, 500 K and 1 MPa, and 700 K and 30 MPa.
pub const m5_models = ['IAPWS-95', 'IF97']!

pub const m5_kinds = ['pt', 'ph']!

pub const m5_states = [[300.0, 3e6], [500.0, 1e6], [700.0, 3e7]]

// M5Case is one program's case: each round's µs per call, and each call's
// answer, h (kJ/kg) for a PT call and T (K) for a PH call; or the error its
// first call raised.
pub struct M5Case {
pub:
	model   string
	kind    string
	us      []f64
	answers []f64
	error   string
}

pub struct M5Source {
pub:
	generated string
	python    string
	dwsim     string
	coolprop  string
	iapws     string
	pythonnet string
}

pub struct M5Fixture {
pub:
	source        M5Source
	machine       Machine
	reflux_build  Build
	reflux_before []M5Case @[json: 'reflux_before']
	reflux_after  []M5Case @[json: 'reflux_after']
	dwsim         []M5Case
	coolprop      []M5Case
	iapws         []M5Case
	seconds       f64
}

// load_m5 reads M5a's benchmark's fixture.
pub fn load_m5() !M5Fixture {
	text := os.read_file(m5_fixture_path) or {
		return error('${m5_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m5_bench.py`.')
	}
	return json2.decode[M5Fixture](text)!
}

// m5_case is a program's case, or one with nothing.
pub fn m5_case(cases []M5Case, model string, kind string) M5Case {
	for c in cases {
		if c.model == model && c.kind == kind {
			return c
		}
	}
	return M5Case{}
}

// m5_reflux_us is reflux's median µs per call in a case: the faster of its
// runs before the others' and after.
pub fn m5_reflux_us(f M5Fixture, model string, kind string) f64 {
	a := median_of(m5_case(f.reflux_before, model, kind).us)
	b := median_of(m5_case(f.reflux_after, model, kind).us)
	return if a < b { a } else { b }
}

// m5_agreement is the largest difference between another program's answers
// and reflux's in one case.
pub fn m5_agreement(f M5Fixture, other []M5Case, model string, kind string) f64 {
	r := m5_case(f.reflux_before, model, kind).answers
	o := m5_case(other, model, kind).answers
	mut d := 0.0
	for k, v in o {
		e := if v > r[k] { v - r[k] } else { r[k] - v }
		d = if e > d { e } else { d }
	}
	return d
}

// m5_section is M5a's benchmark's part of docs/benchmarks.md.
pub fn m5_section(f M5Fixture) string {
	mut b := strings.new_builder(4096)
	b.writeln('')
	b.writeln('## Water and steam')
	b.writeln('')
	b.writeln("M5a step 10, exit criterion 10: water at 300 K and 3 MPa (a liquid), 500 K and 1 MPa (a vapour) and 700 K and 30 MPa (dense, above Tc and pc), by IAPWS-95 and IF97, in each program. `tools/gen_m5_bench.py` writes `testdata/m5_bench.json` on ${f.machine.cpu}, ${f.machine.os}, on ${f.source.generated}, with DWSIM ${dwsim_version(f.source.dwsim)}'s CoolProp and Steam Tables (IAPWS-IF97) packages, CoolProp ${f.source.coolprop}'s AbstractState (HEOS and IF97) from Python, the iapws package ${f.source.iapws}, and reflux built with ${f.reflux_build.compiler}. reflux's IAPWS-95 calls are its PT and PH flashes, phase stability and all; its IF97 calls `water.if97.state_pt` and `state_ph`. A PH call goes back from the h each program's own PT call gave. Each time is the median of its rounds, per call; reflux ran before the others and again after, and the faster is shown.")
	b.writeln('')
	b.writeln('| | reflux | CoolProp | iapws | DWSIM | DWSIM over reflux |')
	b.writeln('|---|---|---|---|---|---|')
	for model in m5_models {
		for kind in m5_kinds {
			r := m5_reflux_us(f, model, kind)
			mut cells := []string{}
			for other in [f.coolprop, f.iapws, f.dwsim] {
				c := m5_case(other, model, kind)
				cells << if c.error != '' { '—' } else { duration(median_of(c.us)) }
			}
			what := if kind == 'pt' {
				if model == 'IF97' { '(p, T)' } else { 'PT flash' }
			} else {
				if model == 'IF97' { '(p, h)' } else { 'PH flash' }
			}
			d := median_of(m5_case(f.dwsim, model, kind).us)
			b.writeln('| ${model}, ${what} | ${duration(r)} | ${cells.join(' | ')} | ${times(d, r)} |')
		}
	}
	b.writeln('')
	b.writeln("Every program's h agrees with reflux's to ${sig(m5_worst_h(f, 'IAPWS-95'), 2)} kJ/kg by IAPWS-95 and ${sig(m5_worst_h(f, 'IF97'), 2)} by IF97. Back from h, DWSIM's CoolProp package lands within ${sig(m5_agreement(f, f.dwsim, 'IAPWS-95', 'ph'), 2)} K and its Steam Tables within ${sig(m5_agreement(f, f.dwsim, 'IF97', 'ph'), 2)} K, its PH flash's tolerance; CoolProp's IF97 backend within ${sig(m5_agreement(f, f.coolprop, 'IF97', 'ph'), 2)} K, its backward equations' T unfinished (R7-97 allows 25 mK in region 1), where reflux finishes each on its basic equation (D9). DWSIM's PH flash on its CoolProp package goes through its own single-compound flash, not CoolProp's.")
	return b.str()
}

// m5_worst_h is the largest difference in h between another program's PT
// answers and reflux's, over the three others, kJ/kg.
pub fn m5_worst_h(f M5Fixture, model string) f64 {
	mut d := 0.0
	for other in [f.dwsim, f.coolprop, f.iapws] {
		e := m5_agreement(f, other, model, 'pt')
		d = if e > d { e } else { d }
	}
	return d
}
