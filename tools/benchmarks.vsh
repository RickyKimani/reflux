#!/usr/bin/env -S v

// benchmarks.vsh writes docs/benchmarks.md from testdata/m3_bench.json,
// M3's benchmarks (decision D10, exit criterion 8), which
// tools/gen_m3_bench.py writes with all three programs' timings:
//
//     v run tools/benchmarks.vsh          # write the report
//     v -prod -cc gcc run tools/benchmarks.vsh --time
//                                         # also time reflux now, for regressions
//
// --time prints today's reflux next to the fixture's, without changing the
// report: the fixture's were timed case by case with thermo's and DWSIM's,
// and a reflux timed later is timed on a machine in another state.
import benchmarks

fx := benchmarks.load()!
write_file(benchmarks.report_path, benchmarks.report(fx))!
println('wrote ${benchmarks.report_path}')
if '--time' in arguments() {
	$if !prod {
		eprintln('time with -prod: v -prod -cc gcc run tools/benchmarks.vsh --time')
		exit(1)
	}
	println('reflux now against the fixture (${fx.source.generated}), median µs:')
	for c in fx.cases {
		s := fx.systems[c.system]
		t := benchmarks.time_calls(c.kind, benchmarks.prepare(benchmarks.fluid(s)!, s, c)!)!
		println('  ${s.name}, ${c.kind}: ${t.median_us:.1f} now, ${c.reflux.median_us:.1f} then')
	}
}
