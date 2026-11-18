// M8's benchmark (step 14, exit criterion 11): what each program solved,
// that reflux is faster, and that the committed report shows it. The
// timings depend on the machine, so they aren't asserted beyond that;
// tools/benchmarks.vsh writes them to docs/benchmarks.md.
import benchmarks
import flowsheet
import math
import os
import reactors
import study

// The owner's study as tools/bench_m8.v runs it.
fn owners_study() !flowsheet.Loaded {
	path := os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reformer_dwsim.toml')
	text := os.read_file(path)!.replace('format = 2', 'format = 3') +
		'\n[study.owner]\nvary = [{ key = "feed.feed.T", from = "700 degC", to = "1000 degC", points = 50 }]\nreport = { methane = "product.n[CH4]" }\n'
	return flowsheet.load_text(text, makers: reactors.makers(), dir: os.dir(path))
}

fn test_reflux_solved_what_it_says() {
	// The study's methane now, under this build, against the fixture's,
	// gcc -prod's: within 1e-9, the solves' tolerance; on 1, 2 and 4
	// workers the fixture's are the same bits (D7); cold and warm within
	// 1e-8 of each other, each a converged solve (criterion 4).
	f := benchmarks.load_m8()!
	l := owners_study()!
	r := study.run(l, l.studies[0])
	for round in [f.reflux_before, f.reflux_after] {
		warm := round.study['warm'] or { panic('no warm study') }.methane
		assert warm.len == 50
		for k, p in r.points {
			assert math.abs(p.outputs[0] - warm[k]) <= 1e-9 * warm[k], '${k}'
		}
		for label in ['warm, 2 workers', 'warm, 4 workers'] {
			assert (round.study[label] or { panic(label) }).methane == warm, label
		}
		cold := round.study['cold'] or { panic('no cold study') }.methane
		for k in 0 .. 50 {
			assert math.abs(cold[k] - warm[k]) <= 1e-8 * warm[k], '${k}'
		}
		assert round.study.keys().len == benchmarks.m8_studies.len
	}
	assert f.reflux_before.optimize.solves == f.reflux_after.optimize.solves
	assert math.abs(f.reflux_before.optimize.p - 405503.96) <= 0.5
	assert math.abs(f.reflux_before.spec.t - 219.7574877222) <= 1e-8
}

fn test_dwsim_solved_the_same_case() {
	// DWSIM's study by its own form: the owner's 50 feed temperatures, and
	// methane's outflow within 2.4e-3 of the feed's methane of reflux's, as
	// M9 found its conversions (its sections' error and reflux's own
	// viscosity, M9 steps 0 and 10). DWSIM's optimisation and spec ended at
	// step 13's answers.
	f := benchmarks.load_m8()!
	rows := f.dwsim.study.rows
	assert rows.len == 50
	assert math.abs(rows[0].t - 973.15) <= 1e-6 && math.abs(rows.last().t - 1273.15) <= 1e-6
	warm := (f.reflux_before.study['warm'] or { panic('no warm study') }).methane
	feed := 10.0 / 3.6
	mut worst := 0.0
	for k, row in rows {
		worst = math.max(worst, math.abs(row.methane - warm[k]) / feed)
	}
	println("reflux within ${worst:.2e} of DWSIM across the study, of the feed's methane")
	assert worst <= 2.4e-3
	assert math.abs(f.dwsim.optimize.p - 405502.5) <= 1.0
	assert math.abs(f.dwsim.spec.t - 219.7546705) <= 1e-6
}

fn test_reflux_is_faster() {
	// Criterion 11: faster than DWSIM's own forms, each of them, before
	// DWSIM's runs and after.
	f := benchmarks.load_m8()!
	for round in [f.reflux_before, f.reflux_after] {
		for label in benchmarks.m8_studies {
			ms := (round.study[label] or { panic(label) }).ms
			assert ms.len == 3 && f.dwsim.study.s * 1e3 > least(ms), label
		}
		assert f.dwsim.optimize.s * 1e3 > least(round.optimize.ms)
		assert f.dwsim.spec.s * 1e3 > least(round.spec.ms)
	}
}

fn test_the_committed_report_has_the_m8_section() {
	f := benchmarks.load_m8()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains(benchmarks.m8_section(f))
}

fn least(xs []f64) f64 {
	mut m := math.inf(1)
	for x in xs {
		m = math.min(m, x)
	}
	return m
}
