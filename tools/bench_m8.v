// bench_m8 times reflux for M8 step 14 (criterion 11), as
// tools/gen_m8_bench.py runs it between DWSIM's own runs of the same:
//
//     v -prod -cc gcc -o bench_m8.exe tools/bench_m8.v
//     bench_m8.exe --build                              # how it was built
//     bench_m8.exe study FILE WORKERS WARM REPEATS      # the owner's 50-point study
//     bench_m8.exe solve FILE REPEATS                   # one solve of a file
//     bench_m8.exe optimize FILE NAME REPEATS           # one of a file's optimisations
//     bench_m8.exe spec FILE REPEATS                    # D18's spec on the reference flowsheet
//
// Each mode loads the file once, runs once untimed, then times REPEATS runs,
// and prints one JSON object. The study is the owner's: the feed of FILE, a
// format-2 flowsheet, from 700 to 1000 °C in 50 points, in chunks of 8,
// methane's outflow at each; WARM is 1 for each chunk's points from their
// neighbours, as a study runs, and 0 for every point cold. It's a .v program,
// not a .vsh script, because a script runs whenever V compiles it
// (docs/v-for-go.md, §1).
module main

import benchmarks
import flowsheet
import json2
import os
import reactors
import study
import time

struct Timed {
	ms          []f64
	values      []f64 // what each mode measures: the study's methane, an optimum's inputs, a spec's values
	objective   f64
	iterations  int
	evaluations int
	solves      int
	passes      int
	error       string
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

fn load_file(path string, extra string) !flowsheet.Loaded {
	if extra == '' {
		return flowsheet.load(path, makers: reactors.makers())
	}
	text := os.read_file(path)!.replace('format = 1', 'format = 3').replace('format = 2', 'format = 3') +
		extra
	return flowsheet.load_text(text, makers: reactors.makers(), dir: os.dir(path))
}

// D18's spec, on the reference flowsheet: the cooler's outlet T until the
// first drum's vapour is 60 mol/s.
const vapour_spec = '
[spec.vapour]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"
'

const owners_study = '\n[study.owner]\nvary = [{ key = "feed.feed.T", from = "700 degC", to = "1000 degC", points = 50 }]\nreport = { methane = "product.n[CH4]" }\n'

fn run_study(path string, workers int, warm bool, repeats int) !Timed {
	l := load_file(path, owners_study)!
	plan := l.studies[0]
	mut ms := []f64{}
	mut last := study.run(l, plan, workers: workers, warm: warm)
	for _ in 0 .. repeats {
		t0 := time.sys_mono_now()
		last = study.run(l, plan, workers: workers, warm: warm)
		ms << ms_since(t0)
	}
	if !last.points.all(it.status == .solved) {
		return error('a point did not solve')
	}
	return Timed{
		ms:     ms
		values: last.points.map(it.outputs[0])
	}
}

fn run_solve(path string, repeats int) !Timed {
	l := load_file(path, '')!
	_ := l.solve()!
	mut ms := []f64{}
	for _ in 0 .. repeats {
		t0 := time.sys_mono_now()
		_ := l.solve()!
		ms << ms_since(t0)
	}
	return Timed{
		ms: ms
	}
}

fn run_optimize(path string, name string, repeats int) !Timed {
	l := load_file(path, '')!
	plans := l.optimizations.filter(it.name == name)
	if plans.len == 0 {
		return error('no optimisation ${name}')
	}
	mut o := study.optimize(l, plans[0])
	mut ms := []f64{}
	for _ in 0 .. repeats {
		t0 := time.sys_mono_now()
		o = study.optimize(l, plans[0])
		ms << ms_since(t0)
	}
	if !o.ok {
		return error(o.error)
	}
	return Timed{
		ms:          ms
		values:      o.inputs
		objective:   o.objective
		iterations:  o.iterations.len
		evaluations: o.evaluations
		solves:      o.solves
	}
}

fn run_spec(path string, repeats int) !Timed {
	l := load_file(path, vapour_spec)!
	mut sol := l.solve()!
	mut ms := []f64{}
	for _ in 0 .. repeats {
		t0 := time.sys_mono_now()
		sol = l.solve()!
		ms << ms_since(t0)
	}
	mut passes := 0
	for h in sol.history {
		passes += h.residuals.len
	}
	return Timed{
		ms:     ms
		values: sol.specs.map(it.value)
		passes: passes
	}
}

// dispatch runs the mode the arguments name.
fn dispatch(args []string) !Timed {
	mode := args[0] or { '' }
	if mode == 'study' {
		return run_study(args[1] or { '' }, (args[2] or { '1' }).int(), (args[3] or { '1' }) == '1',
			(args[4] or { '1' }).int())
	}
	if mode == 'solve' {
		return run_solve(args[1] or { '' }, (args[2] or { '1' }).int())
	}
	if mode == 'optimize' {
		return run_optimize(args[1] or { '' }, args[2] or { '' }, (args[3] or { '1' }).int())
	}
	if mode == 'spec' {
		return run_spec(args[1] or { '' }, (args[2] or { '1' }).int())
	}
	return error('usage: bench_m8 --build | study FILE WORKERS WARM REPEATS | solve FILE REPEATS | optimize FILE NAME REPEATS | spec FILE REPEATS')
}

fn main() {
	$if !prod {
		eprintln('build with -prod: v -prod -cc gcc -o bench_m8.exe tools/bench_m8.v')
		exit(1)
	}
	args := os.args[1..]
	if args == ['--build'] {
		println(json2.encode(benchmarks.Build{
			compiler: benchmarks.compiler()
			v:        @VCURRENTHASH
		}))
		return
	}
	t := dispatch(args) or {
		Timed{
			error: err.msg()
		}
	}
	println(json2.encode(t))
	if t.error != '' {
		exit(1)
	}
}
