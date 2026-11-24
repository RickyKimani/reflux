module benchmarks

import json2
import os
import strings

// M6's speed (step 13, D18, exit criterion 12): γ–φ flashes of methanol and
// water by reflux, DWSIM 9 and thermo in one session, and Peng–Robinson's
// beside them. tools/gen_m6_bench.py times all three, reflux through
// tools/bench_m6.v, and writes testdata/m6_bench.json; tools/benchmarks.vsh
// adds it to docs/benchmarks.md.

pub const m6_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm6_bench.json')

// The cases, in the report's order: each model's PT flashes and bubble
// points.
pub const m6_models = ['NRTL', 'UNIFAC', 'PR']!

pub const m6_kinds = ['pt', 'bubble']!

// M6Case is one program's case: each round's µs per call, and each call's
// answer, the vapour fraction or the bubble T; or the error its first call
// raised.
pub struct M6Case {
pub:
	model   string
	kind    string
	us      []f64
	answers []f64
	error   string
}

pub struct M6Fixture {
pub:
	source        Source
	machine       Machine
	reflux_build  Build
	reflux_before []M6Case @[json: 'reflux_before']
	reflux_after  []M6Case @[json: 'reflux_after']
	dwsim         []M6Case
	thermo        []M6Case
	seconds       f64
}

// load_m6 reads M6's benchmark's fixture.
pub fn load_m6() !M6Fixture {
	text := os.read_file(m6_fixture_path) or {
		return error('${m6_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m6_bench.py`.')
	}
	return json2.decode[M6Fixture](text)!
}

// m6_case is a program's case, or one with nothing.
pub fn m6_case(cases []M6Case, model string, kind string) M6Case {
	for c in cases {
		if c.model == model && c.kind == kind {
			return c
		}
	}
	return M6Case{}
}

// m6_reflux_us is reflux's median µs per call in a case: the faster of its
// runs before the others' and after.
pub fn m6_reflux_us(f M6Fixture, model string, kind string) f64 {
	a := median_of(m6_case(f.reflux_before, model, kind).us)
	b := median_of(m6_case(f.reflux_after, model, kind).us)
	return if a < b { a } else { b }
}

// m6_agreement is the largest difference between another program's
// answers and reflux's, over some models: in β, and in the bubble T, K.
pub fn m6_agreement(f M6Fixture, other []M6Case, models []string) (f64, f64) {
	mut beta, mut t := 0.0, 0.0
	for model in models {
		for kind in m6_kinds {
			r := m6_case(f.reflux_before, model, kind).answers
			o := m6_case(other, model, kind).answers
			for k, v in o {
				d := if v > r[k] { v - r[k] } else { r[k] - v }
				if kind == 'pt' {
					beta = if d > beta { d } else { beta }
				} else {
					t = if d > t { d } else { t }
				}
			}
		}
	}
	return beta, t
}

// m6_section is M6's benchmark's part of docs/benchmarks.md.
pub fn m6_section(f M6Fixture) string {
	mut b := strings.new_builder(4096)
	b.writeln('')
	b.writeln('## γ–φ flashes')
	b.writeln('')
	b.writeln("M6 step 13, D18, exit criterion 12: methanol and water at 1 atm, as γ–φ fluids by NRTL and UNIFAC on ChemSep's parameters (the vapour the ideal gas, Poynting on) and by Peng–Robinson on ChemSep's kᵢⱼ, in each program. `tools/gen_m6_bench.py` writes `testdata/m6_bench.json` on ${f.machine.cpu}, ${f.machine.os}, on ${f.source.generated}, with DWSIM ${dwsim_version(f.source.dwsim)} at its defaults, thermo ${f.source.thermo}'s FlashVL, and reflux built with ${f.reflux_build.compiler}. A PT flash is at 350 K (z = 0.5), 357 K (0.3) and 346 K (0.7), each two-phase; a bubble point is the T of a liquid of 0.1, 0.5 and 0.9 methanol. Each time is the median of its rounds, per call; reflux ran before the others and again after, and the faster is shown.")
	b.writeln('')
	b.writeln('| | reflux | thermo | DWSIM | thermo over reflux | DWSIM over reflux |')
	b.writeln('|---|---|---|---|---|---|')
	for model in m6_models {
		for kind in m6_kinds {
			r := m6_reflux_us(f, model, kind)
			t := m6_case(f.thermo, model, kind)
			d := m6_case(f.dwsim, model, kind)
			tu, du := median_of(t.us), median_of(d.us)
			what := if kind == 'pt' { 'PT flash' } else { 'bubble T' }
			name := if model == 'PR' { 'Peng–Robinson' } else { model }
			b.writeln('| ${name}, ${what} | ${duration(r)} | ${if t.error != '' {
				'—'
			} else {
				duration(tu)
			}} | ${if d.error != '' { '—' } else { duration(du) }} | ${times(tu, r)} | ${times(du, r)} |')
		}
	}
	nrtl := m6_reflux_us(f, 'NRTL', 'pt')
	pr := m6_reflux_us(f, 'PR', 'pt')
	b.writeln('')
	tb, tt := m6_agreement(f, f.thermo, ['NRTL', 'UNIFAC', 'PR'])
	db, dt := m6_agreement(f, f.dwsim, ['NRTL', 'UNIFAC'])
	b.writeln("reflux's γ–φ PT flash takes ${times(nrtl, pr)} of its Peng–Robinson flash's time on the same states: the liquid is xᵢγᵢPᵢˢᵃᵗφᵢˢᵃᵗ and the Poynting term, with no volume roots to find, as D18 expected. DWSIM's activity packages run ${times(median_of(m6_case(f.dwsim, 'NRTL', 'pt').us), median_of(m6_case(f.dwsim, 'PR', 'pt').us))} times slower than its own Peng–Robinson. thermo's answers are reflux's to ${sig(tb, 2)} in β and ${sig(tt, 2)} K; DWSIM's γ–φ answers differ by up to ${sig(db, 2)} in β and ${sig(dt, 2)} K, its Poynting term's sign and its vapour pressures' extra figures (M6 step 12).")
	return b.str()
}
