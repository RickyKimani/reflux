module benchmarks

import json2
import os
import strings

// M9's speed (step 11, exit criterion 9): the owner's reforming PFR, solved
// by DWSIM 9 and by reflux, one solve and the owner's 50-point study, the two
// programs taking turns. tools/gen_m9_bench.py times DWSIM and reflux,
// through tools/bench_reactor.v, and writes testdata/m9_bench.json;
// tools/benchmarks.vsh adds it to docs/benchmarks.md.

pub const reactor_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm9_bench.json')

// The reflux files timed, in the report's order: the owner's case on DWSIM's
// rates as typed, at reflux's default rtol and at 1e-4, and the published
// model.
pub const reactor_labels = ['reflux', 'reflux at rtol 1e-4', 'reflux, the published model']!

// ReactorSingle is the single solve at the case's 1000 K: DWSIM's in
// seconds, by its default integrator and its explicit Runge–Kutta; each
// reflux file's in milliseconds; and what each came to, methane's
// conversion, with reflux's integrator's steps.
pub struct ReactorSingle {
pub:
	dwsim_s      []f64
	dwsim_rk_s   []f64
	dwsim_x      f64
	dwsim_rk_x   f64
	reflux_ms    map[string][]f64
	reflux_x     map[string]f64
	reflux_steps map[string]int
}

// ReactorRun is one reflux file at one point of the study.
pub struct ReactorRun {
pub:
	ms    []f64
	x     f64
	steps int
}

// ReactorPoint is one point of the study: the feed's T, DWSIM's solve and
// conversion, and each reflux file's.
pub struct ReactorPoint {
pub:
	t       f64
	dwsim_s f64
	dwsim_x f64
	reflux  map[string]ReactorRun
}

pub struct ReactorFixture {
pub:
	source       Source
	machine      Machine
	reflux_build Build
	load_ms      map[string]f64
	single       ReactorSingle
	study        []ReactorPoint
}

// load_reactor reads the reactor benchmark's fixture.
pub fn load_reactor() !ReactorFixture {
	text := os.read_file(reactor_fixture_path) or {
		return error('${reactor_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m9_bench.py`.')
	}
	return json2.decode[ReactorFixture](text)!
}

// seconds prints seconds to three significant figures, in minutes from two
// minutes up.
fn seconds(s f64) string {
	return if round_sig(s, 3) < 120.0 { '${sig(s, 3)} s' } else { '${sig(s / 60, 3)} min' }
}

fn total(xs []f64) f64 {
	mut t := 0.0
	for x in xs {
		t += x
	}
	return t
}

// reflux_point_ms is a reflux file's median at each point of the study.
fn reflux_point_ms(rf ReactorFixture, label string) []f64 {
	return rf.study.map(median_of(it.reflux[label] or { ReactorRun{} }.ms))
}

// times is how many times longer a is than b, to two significant figures.
fn times(a f64, b f64) string {
	return if b > 0 { sig(a / b, 2) } else { '—' }
}

// reactor_section is the reactor benchmark's part of docs/benchmarks.md.
pub fn reactor_section(rf ReactorFixture) string {
	mut b := strings.new_builder(4096)
	s := rf.single
	b.writeln('')
	b.writeln("## The owner's reforming PFR")
	b.writeln('')
	b.writeln("M9 step 11, exit criterion 9: the owner's coursework case (M9 step 0), one isothermal packed bed of 1 m³ on three LHHW rates, with Ergun's drop, its feed at 1000 K and 1 atm. `tools/gen_m9_bench.py` writes `testdata/m9_bench.json`, on ${rf.machine.cpu}, ${rf.machine.os}, on ${rf.source.generated}, the two programs taking turns: the single solve in five turns, and each point of the study solved by DWSIM and then by reflux before the next, since a laptop's clock moved DWSIM's solve between 9.4 and 29 s in step 0.")
	b.writeln('')
	b.writeln("- **DWSIM ${dwsim_version(rf.source.dwsim)}** as the owner's file sets it up: sections of dV = 0.05, Ergun's drop, its default integrator (Radau); and for the single solve its explicit Runge–Kutta 4(5), which step 0 found seven times faster.")
	b.writeln("- **reflux,** built with ${rf.reflux_build.compiler}, on `testdata/flowsheets/reformer_dwsim.toml`, the same case on DWSIM's rates as typed, with reflux's own viscosity: at its default rtol of 1e-8, and at 1e-4. And on `testdata/flowsheets/reformer.toml`, the published model. A solve includes the feed's flash; reading the file, the data and the fluid took ${duration((rf.load_ms['reflux'] or {
		0.0
	}) * 1e3)}, once.")
	b.writeln('')
	b.writeln('### A single solve, at 1000 K')
	b.writeln('')
	b.writeln("| | A solve | Fastest | DWSIM's default over it | Methane's conversion | Steps |")
	b.writeln('|---|---|---|---|---|---|')
	dw := median_of(s.dwsim_s)
	b.writeln('| DWSIM, Radau, dV = 0.05 | ${seconds(dw)} | ${seconds(minimum(s.dwsim_s))} | 1 | ${decimals(s.dwsim_x, 6)} | — |')
	b.writeln('| DWSIM, explicit Runge–Kutta | ${seconds(median_of(s.dwsim_rk_s))} | ${seconds(minimum(s.dwsim_rk_s))} | ${times(dw, median_of(s.dwsim_rk_s))} | ${decimals(s.dwsim_rk_x, 6)} | — |')
	for label in reactor_labels {
		xs := s.reflux_ms[label] or { []f64{} }
		m := median_of(xs)
		b.writeln('| ${label} | ${ms(xs, median_of)} | ${ms(xs, minimum)} | ${times(dw * 1e3, m)} | ${decimals(s.reflux_x[label] or { 0.0 }, 6)} | ${s.reflux_steps[label] or {
			0
		}} |')
	}
	b.writeln('')
	b.writeln("Medians, but for the fastest. DWSIM's conversion is 1.3e-3 below the exact answer of its own equations with Ergun's drop, its sections' error (step 0). reflux's at 1e-8 is that exact answer to 2.3e-10 when it's given DWSIM's viscosity (step 10); on its own, Wilke's rule, it's 9.4e-4 above. Its steps are its integrator's accepted ones.")
	b.writeln('')
	b.writeln('### The study: 50 feed temperatures, 700 to 1000 °C')
	b.writeln('')
	b.writeln('| | All 50 | A point | Slowest point | Over reflux at 1e-8 |')
	b.writeln('|---|---|---|---|---|')
	dws := rf.study.map(it.dwsim_s)
	main_ms := reflux_point_ms(rf, 'reflux')
	b.writeln('| DWSIM, Radau, dV = 0.05 | ${seconds(total(dws))} | ${seconds(median_of(dws))} | ${seconds(maximum(dws))} | ${times(total(dws) * 1e3, total(main_ms))} |')
	for label in reactor_labels {
		xs := reflux_point_ms(rf, label)
		b.writeln('| ${label} | ${ms([total(xs)], median_of)} | ${ms(xs, median_of)} | ${ms(xs, maximum)} | ${times(total(xs), total(main_ms))} |')
	}
	b.writeln('')
	b.writeln("Each reflux point is its median of three solves. Every seventh point, and the last, with each program's conversion:")
	b.writeln('')
	b.writeln("| Feed | DWSIM | reflux | DWSIM's conversion | reflux's |")
	b.writeln('|---|---|---|---|---|')
	for i, p in rf.study {
		if i % 7 != 0 && i != rf.study.len - 1 {
			continue
		}
		r := p.reflux['reflux'] or { ReactorRun{} }
		b.writeln('| ${sig(p.t - 273.15, 4)} °C | ${seconds(p.dwsim_s)} | ${ms(r.ms, median_of)} | ${decimals(p.dwsim_x, 6)} | ${decimals(r.x, 6)} |')
	}
	return b.str()
}

fn maximum(xs []f64) f64 {
	mut m := xs[0] or { return 0.0 }
	for x in xs {
		if x > m {
			m = x
		}
	}
	return m
}
