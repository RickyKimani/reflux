// M4's flowsheet benchmark (exit criterion 9): what each program solved, and
// that the committed report shows it. The timings depend on the machine, so
// they aren't asserted; tools/benchmarks.vsh writes them to docs/benchmarks.md.
import benchmarks
import flowsheet
import math
import os

fn test_reflux_solved_what_it_says() {
	// reflux's timings were of this solve: the same passes and products now.
	ff := benchmarks.load_flowsheets()!
	l := flowsheet.load(os.join_path(@VMODROOT, ff.flowsheet))!
	for t in ff.timings.filter(it.label.starts_with('reflux')) {
		assert t.error == '' && t.cold_ms.len > 0 && t.warm_ms.len > 0 && t.load_ms.len == t.cold_ms.len, t.label
		method := if t.label.ends_with('direct substitution') {
			flowsheet.Method.direct
		} else {
			flowsheet.Method.wegstein
		}
		sol := l.sheet.solve(l.fluid, flowsheet.SolveParams{ ...l.params, method: method })!
		assert sol.history[0].residuals.len == t.passes, t.label
		for name, f in t.products {
			assert math.abs(sol.streams[name].total() - f) <= 1e-9 * f, '${t.label} ${name}'
		}
	}
}

fn test_dwsim_solved_the_same_flowsheet() {
	// At 1e-8, DWSIM's products are reflux's within criterion 2's 1e-4. At its
	// defaults it stops short: its balances 0.2% open (step 0), its flash gas
	// 0.5% low.
	ff := benchmarks.load_flowsheets()!
	reflux := ff.timings.filter(it.label == 'reflux, Wegstein')[0]
	for t in ff.timings.filter(it.label.starts_with('DWSIM')) {
		assert t.error == '' && t.cold_ms.len > 0 && t.warm_ms.len > 0, t.label
		tol := if t.label == 'DWSIM at 1e-8' { 1e-4 } else { 6e-3 }
		for name, f in reflux.products {
			assert math.abs(t.products[name] - f) <= tol * f, '${t.label} ${name}: ${t.products[name]} against ${f}'
		}
	}
}

fn test_the_committed_report_has_every_timing() {
	ff := benchmarks.load_flowsheets()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains('## The reference flowsheet')
	// A row of each table for each, timings and products.
	lines := text.all_after('## The reference flowsheet').split_into_lines()
	for t in ff.timings {
		assert lines.filter(it.starts_with('| ${t.label} | ')).len == 2, t.label
	}
}
