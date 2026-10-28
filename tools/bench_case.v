// bench_case times reflux on one benchmark case and prints the timing as
// JSON. tools/gen_m3_bench.py builds it and runs it between thermo's and
// DWSIM's timings of each case, so that all three programs time a case
// within seconds of each other:
//
//     v -prod -cc gcc -o bench_case.exe tools/bench_case.v
//     bench_case.exe case.json    # {"system": …, "flash": …}, as in testdata/m3_bench.json
//     bench_case.exe --build      # how it was built
//
// It's a .v program, not a .vsh script, because a script runs whenever V
// compiles it (see docs/v-for-go.md, §1).
module main

import benchmarks
import json2
import os

struct Input {
	system benchmarks.System
	flash  benchmarks.Case
}

fn main() {
	$if !prod {
		eprintln('build with -prod: v -prod -cc gcc -o bench_case.exe tools/bench_case.v')
		exit(1)
	}
	if os.args.len != 2 {
		eprintln('usage: bench_case case.json, or bench_case --build')
		exit(2)
	}
	if os.args[1] == '--build' {
		println(json2.encode(benchmarks.Build{
			compiler: benchmarks.compiler()
			v:        @VCURRENTHASH
		}))
		return
	}
	timing := time_case(os.args[1]) or {
		// A failed flash is a result for the fixture, not a crash.
		benchmarks.Timing{
			error: err.msg()
		}
	}
	println(json2.encode(timing))
}

fn time_case(path string) !benchmarks.Timing {
	input := json2.decode[Input](os.read_file(path)!)!
	fl := benchmarks.fluid(input.system)!
	return benchmarks.time_calls(input.flash.kind, benchmarks.prepare(fl, input.system, input.flash)!)
}
