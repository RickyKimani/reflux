module benchmarks

import json2
import os
import strings

// M8's speed (step 14, exit criterion 11): the owner's 50-point study, an
// optimisation and a spec, by DWSIM 9's own forms and by reflux, in one
// session. tools/gen_m8_bench.py times both, reflux through
// tools/bench_m8.v, and writes testdata/m8_bench.json; tools/benchmarks.vsh
// adds it to docs/benchmarks.md.

pub const m8_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm8_bench.json')

// The reflux studies timed, in the report's order.
pub const m8_studies = ['cold', 'warm', 'warm, 2 workers', 'warm, 4 workers']!

pub struct M8Study {
pub:
	ms      []f64
	methane []f64 // the product's methane at each point, mol/s
}

pub struct M8Times {
pub:
	ms []f64
}

pub struct M8Optimize {
pub:
	ms          []f64
	p           f64 @[json: 'P']
	work        f64
	iterations  int
	evaluations int
	solves      int
}

pub struct M8Spec {
pub:
	ms     []f64
	t      f64 @[json: 'T']
	passes int
}

// M8Round is everything reflux ran, once: before DWSIM's runs, or after.
pub struct M8Round {
pub:
	study     map[string]M8Study
	published M8Times
	optimize  M8Optimize
	spec      M8Spec
}

pub struct M8Row {
pub:
	t       f64 @[json: 'T']
	methane f64
}

pub struct M8DwsimStudy {
pub:
	s    f64
	rows []M8Row
}

pub struct M8DwsimOptimize {
pub:
	s           f64
	p           f64 @[json: 'P']
	evaluations int
	why         string
}

pub struct M8DwsimSpec {
pub:
	s      f64
	t      f64 @[json: 'T']
	solves int
}

pub struct M8Dwsim {
pub:
	study    M8DwsimStudy
	optimize M8DwsimOptimize
	spec     M8DwsimSpec
}

pub struct M8Fixture {
pub:
	source        Source
	machine       Machine
	reflux_build  Build
	reflux_before M8Round @[json: 'reflux_before']
	reflux_after  M8Round @[json: 'reflux_after']
	dwsim         M8Dwsim
	seconds       f64
}

// load_m8 reads M8's benchmark's fixture.
pub fn load_m8() !M8Fixture {
	text := os.read_file(m8_fixture_path) or {
		return error('${m8_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m8_bench.py`.')
	}
	return json2.decode[M8Fixture](text)!
}

// best is the fastest of a reflux run's times, before DWSIM's runs or
// after, whichever was faster: each is the best of three.
fn best(a []f64, b []f64) f64 {
	return minimum([minimum(a), minimum(b)])
}

// m8_section is M8's benchmark's part of docs/benchmarks.md.
pub fn m8_section(f M8Fixture) string {
	mut b := strings.new_builder(4096)
	before, after := f.reflux_before, f.reflux_after
	b.writeln('')
	b.writeln('## Studies, an optimisation and a spec')
	b.writeln('')
	b.writeln("M8 step 14, exit criterion 11: what M8 added, each against DWSIM's own form for it, driven unshown (M8 step 0). `tools/gen_m8_bench.py` writes `testdata/m8_bench.json`, on ${f.machine.cpu}, ${f.machine.os}, on ${f.source.generated}, with DWSIM ${dwsim_version(f.source.dwsim)} and reflux built with ${f.reflux_build.compiler}. reflux ran everything before DWSIM's runs and again after, each the best of three after an untimed run, since a laptop's clock moves DWSIM's solve between 9.4 and 29 s (M9 step 0); the faster of the two is shown.")
	b.writeln('')
	b.writeln("### The owner's study: 50 feed temperatures, 700 to 1000 °C")
	b.writeln('')
	dw := f.dwsim.study.s
	b.writeln('| | All 50 | DWSIM over it |')
	b.writeln('|---|---|---|')
	b.writeln("| DWSIM's sensitivity analysis, the owner's setup (dV = 0.05, Radau) | ${seconds(dw)} | 1 |")
	for label in m8_studies {
		x := best(before.study[label] or { M8Study{} }.ms, after.study[label] or { M8Study{} }.ms)
		b.writeln('| reflux, ${label} | ${ms([x], median_of)} | ${times(dw * 1e3, x)} |')
	}
	b.writeln('')
	b.writeln("reflux solves `testdata/flowsheets/reformer_dwsim.toml`, the same case on DWSIM's rates as typed, in chunks of 8: cold, every point from the feed; warm, each chunk's points from their neighbours (D6, D7); and warm on 2 and 4 threads (D8). The answers are the same bits on any number of workers. One solve of the published model, `testdata/flowsheets/reformer.toml`, whose K(T) is kept until T moves (step 1), takes ${ms([
		best(before.published.ms, after.published.ms),
	], median_of)}: it took 393 ms in M9.")
	b.writeln('')
	b.writeln('### An optimisation and a spec')
	b.writeln('')
	o, s := f.dwsim.optimize, f.dwsim.spec
	ro := best(before.optimize.ms, after.optimize.ms)
	rs := best(before.spec.ms, after.spec.ms)
	b.writeln('| | DWSIM | Its solves | reflux | Its solves | DWSIM over reflux |')
	b.writeln('|---|---|---|---|---|---|')
	b.writeln('| Two-stage compression, the least work (`compression_dwsim.toml`) | ${seconds(o.s)}, bounded Brent | ${o.evaluations} | ${ms([
		ro,
	], median_of)}, SQP | ${before.optimize.solves} | ${times(o.s * 1e3, ro)} |')
	b.writeln("| The reference flowsheet's spec (`reference_dwsim.toml`) | ${seconds(s.s)}, an Adjust's secant at 1e-4 | ${s.solves} | ${ms([
		rs,
	], median_of)}, with its recycle | ${before.spec.passes} passes | ${times(s.s * 1e3, rs)} |")
	b.writeln('')
	b.writeln("DWSIM's optimiser and Adjust solve the whole flowsheet at each evaluation. reflux's SQP counts every solve, its gradients' included, five a gradient; its spec is met with the recycle's tears, by Broyden on both, so its count is the block's passes, each one pass of its units.")
	return b.str()
}
