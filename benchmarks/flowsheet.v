module benchmarks

import flowsheet
import json2
import os
import strings
import time

// M4's flowsheet benchmark (D14, exit criterion 9): the reference flowsheet on
// DWSIM's constants, solved by reflux and by DWSIM 9. tools/gen_m4_bench.py
// times DWSIM and runs time_flowsheet through tools/bench_case.v, and writes
// testdata/m4_bench.json; tools/benchmarks.vsh adds it to docs/benchmarks.md.

pub const flowsheet_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm4_bench.json')

// FlowsheetTiming is one program's solves of the flowsheet, in milliseconds,
// and what the first came to.
pub struct FlowsheetTiming {
pub:
	label string // which program, and how it solved
	// load_ms is reflux's reading of the file, the component data and the fluid,
	// each time; DWSIM builds its flowsheet in the script, untimed.
	load_ms []f64
	// cold_ms are first solves: each of a flowsheet just loaded or built.
	cold_ms []f64
	// warm_ms are solves again, each with the feed 0.01 K warmer than the last.
	// DWSIM starts each from its last answer; reflux, which has no warm start
	// until M8, from scratch.
	warm_ms     []f64
	passes      int // the loop's, in the first solve
	warm_passes int // in the last solve again
	// products are the products' flows after the first solve, mol/s.
	products map[string]f64
	error    string
}

pub struct FlowsheetFixture {
pub:
	source       Source
	machine      Machine
	reflux_build Build
	flowsheet    string
	timings      []FlowsheetTiming
}

// load_flowsheets reads the flowsheet benchmark's fixture.
pub fn load_flowsheets() !FlowsheetFixture {
	text := os.read_file(flowsheet_fixture_path) or {
		return error('${flowsheet_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m4_bench.py`.')
	}
	return json2.decode[FlowsheetFixture](text)!
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

// passes is how many passes a solution's loops took, all together.
fn passes(sol flowsheet.Solution) int {
	mut n := 0
	for h in sol.history {
		n += h.residuals.len
	}
	return n
}

// time_flowsheet times a flowsheet file solved by a method: loaded and solved
// cold times over, then solved again warm times, each with every feed 0.01 K
// warmer than the time before.
pub fn time_flowsheet(path string, method flowsheet.Method, cold int, warm int) !FlowsheetTiming {
	mut load_ms, mut cold_ms, mut warm_ms := []f64{}, []f64{}, []f64{}
	t0 := time.sys_mono_now()
	loaded := flowsheet.load(path)!
	load_ms << ms_since(t0)
	params := flowsheet.SolveParams{
		...loaded.params
		method: method
	}
	t1 := time.sys_mono_now()
	first := loaded.sheet.solve(loaded.fluid, params)!
	cold_ms << ms_since(t1)
	for _ in 1 .. cold {
		t2 := time.sys_mono_now()
		again := flowsheet.load(path)!
		load_ms << ms_since(t2)
		t3 := time.sys_mono_now()
		again.sheet.solve(again.fluid, params)!
		cold_ms << ms_since(t3)
	}
	mut last := first
	for k in 1 .. warm + 1 {
		mut feeds := []flowsheet.Feed{}
		for f in loaded.sheet.feeds {
			s := f.stream
			feeds << flowsheet.Feed{f.name, flowsheet.feed_tp(loaded.fluid, s.flows, s.t + 0.01 * k,
				s.p)!}
		}
		sheet := flowsheet.Flowsheet{
			units: loaded.sheet.units
			feeds: feeds
		}
		t4 := time.sys_mono_now()
		last = sheet.solve(loaded.fluid, params)!
		warm_ms << ms_since(t4)
	}
	mut products := map[string]f64{}
	for name in first.products {
		products[name] = first.streams[name].total()
	}
	return FlowsheetTiming{
		label:       'reflux, ${method_name(method)}'
		load_ms:     load_ms
		cold_ms:     cold_ms
		warm_ms:     warm_ms
		passes:      passes(first)
		warm_passes: passes(last)
		products:    products
	}
}

fn method_name(m flowsheet.Method) string {
	return match m {
		.direct { 'direct substitution' }
		.wegstein { 'Wegstein' }
		.broyden { 'Broyden' }
	}
}

// median_of is the middle of xs, or the mean of the middle two.
fn median_of(xs []f64) f64 {
	mut s := xs.clone()
	s.sort()
	n := s.len
	if n == 0 {
		return 0.0
	}
	return if n % 2 == 1 { s[n / 2] } else { 0.5 * (s[n / 2 - 1] + s[n / 2]) }
}

fn minimum(xs []f64) f64 {
	mut m := xs[0] or { return 0.0 }
	for x in xs {
		if x < m {
			m = x
		}
	}
	return m
}

// ms prints a statistic of some milliseconds to three significant figures, in
// seconds from a second up, or a dash for none.
fn ms(xs []f64, f fn ([]f64) f64) string {
	if xs.len == 0 {
		return '—'
	}
	x := f(xs)
	return if round_sig(x, 3) < 1000.0 { duration(x * 1e3) } else { '${sig(x / 1e3, 3)} s' }
}

// flowsheet_section is the flowsheet benchmark's part of docs/benchmarks.md.
pub fn flowsheet_section(ff FlowsheetFixture) string {
	mut b := strings.new_builder(2048)
	b.writeln('')
	b.writeln('## The reference flowsheet')
	b.writeln('')
	b.writeln("M4 step 9, decision D14, exit criterion 9: `${ff.flowsheet}`, the reference flowsheet on DWSIM's constants and kᵢⱼ, solved in each program. `tools/gen_m4_bench.py` writes `testdata/m4_bench.json`, timing DWSIM and then reflux within the same minute, on ${ff.machine.cpu}, ${ff.machine.os}, on ${ff.source.generated}.")
	b.writeln('')
	b.writeln('- **A first solve** is of a flowsheet just built: loaded from its file in reflux, built through the automation in DWSIM, after one untimed first solve so that .NET has compiled it.')
	b.writeln('- **A solve again** follows the feed 0.01 K warmer each time. DWSIM starts it from its last answer. reflux has no warm start until M8 (F5), so its solve again is a first solve in a warm process.')
	b.writeln("- **reflux** at its default tolerance, 1e-9 of each tear variable's scale, by Wegstein (its default) and by direct substitution, built with ${ff.reflux_build.compiler}. Its load is the file, the component data and the fluid.")
	b.writeln("- **DWSIM ${dwsim_version(ff.source.dwsim)}** at its default recycle tolerances, and at them scaled by 1e-8, which closes the balances to its floor (step 0). Its acceleration changes nothing here, so it's off: in effect DWSIM's recycle is direct substitution.")
	b.writeln('')
	b.writeln('| | Passes | First solve | Solve again | Fastest again | Passes again | Load |')
	b.writeln('|---|---|---|---|---|---|---|')
	for t in ff.timings {
		if t.error != '' {
			b.writeln('| ${t.label} | failed: ${t.error} | | | | | |')
			continue
		}
		b.writeln('| ${t.label} | ${t.passes} | ${ms(t.cold_ms, median_of)} | ${ms(t.warm_ms, median_of)} | ${ms(t.warm_ms, minimum)} | ${t.warm_passes} | ${ms(t.load_ms, median_of)} |')
	}
	b.writeln('')
	b.writeln('Medians, but for the fastest. What each solved, the products in mol/s after the first solve:')
	b.writeln('')
	mut names := []string{}
	for t in ff.timings {
		for name, _ in t.products {
			if name !in names {
				names << name
			}
		}
	}
	b.writeln('| | ${names.join(' | ')} |')
	b.writeln('|---${'|---'.repeat(names.len)}|')
	for t in ff.timings {
		mut cells := []string{}
		for name in names {
			cells << if v := t.products[name] { decimals(v, 6) } else { '—' }
		}
		b.writeln('| ${t.label} | ${cells.join(' | ')} |')
	}
	return b.str()
}
