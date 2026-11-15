// M9's reactor benchmark (step 11, exit criterion 9): what each program
// solved, and that the committed report shows it. The timings depend on the
// machine, so they aren't asserted; tools/benchmarks.vsh writes them to
// docs/benchmarks.md.
import benchmarks
import flowsheet
import math
import os
import reactors

// load is a reflux file the benchmark timed, by its label.
fn load(label string) !flowsheet.Loaded {
	file := if label.ends_with('published model') { 'reformer.toml' } else { 'reformer_dwsim.toml' }
	mut text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', file))!
	if label.ends_with('rtol 1e-4') {
		text = text.replace('pressure_drop = "ergun"', 'pressure_drop = "ergun"\nrtol = 1e-4')
	}
	return flowsheet.load_text(text,
		makers: reactors.makers()
		dir:    os.join_path(@VMODROOT,
			'testdata', 'flowsheets')
	)
}

fn conversion(sol flowsheet.Solution) f64 {
	feed := sol.streams['feed'].flows[0]
	return (feed - sol.streams['product'].flows[0]) / feed
}

fn test_reflux_solved_what_it_says() {
	// reflux's timings were of these solves: the same conversion and steps
	// now, at 1000 K and at the study's ends.
	rf := benchmarks.load_reactor()!
	for label in benchmarks.reactor_labels {
		l := load(label)!
		sol := l.solve()!
		assert math.abs(conversion(sol) - rf.single.reflux_x[label]) <= 1e-12, label
		mut steps := 0
		for u in sol.units {
			if p := u.profile {
				steps += p.points.len - 1
			}
		}
		assert steps == rf.single.reflux_steps[label], label
		assert rf.single.reflux_ms[label].len == 50, label
		assert rf.study.all((it.reflux[label] or { benchmarks.ReactorRun{} }).ms.len == 3), label
	}
}

fn test_dwsim_solved_the_same_case() {
	// DWSIM's single solve is step 0's answer, dV = 0.05 with Ergun's drop:
	// X 0.76034558. Across the study, reflux on DWSIM's rates is within
	// 2.4e-3 of DWSIM's conversion (measured 2.29e-3): DWSIM's sections put
	// it 1.3e-3 low at 1000 K, and reflux's own viscosity 9.4e-4 high. Both
	// rise with the feed's T.
	rf := benchmarks.load_reactor()!
	assert math.abs(rf.single.dwsim_x - 0.76034558) <= 1e-8
	assert rf.study.len == 50
	assert math.abs(rf.study[0].t - 973.15) <= 1e-9 && math.abs(rf.study.last().t - 1273.15) <= 1e-9
	mut worst := 0.0
	for i, p in rf.study {
		x := (p.reflux['reflux'] or { benchmarks.ReactorRun{} }).x
		worst = math.max(worst, math.abs(x - p.dwsim_x))
		if i > 0 {
			assert p.dwsim_x > rf.study[i - 1].dwsim_x
			assert x > (rf.study[i - 1].reflux['reflux'] or { benchmarks.ReactorRun{} }).x
		}
	}
	println('reflux within ${worst:.2e} of DWSIM across the study')
	assert worst <= 2.4e-3
}

fn test_the_committed_report_has_the_reactor_section() {
	// The report shows this fixture, as tools/benchmarks.vsh writes it.
	rf := benchmarks.load_reactor()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains(benchmarks.reactor_section(rf))
}
