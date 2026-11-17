// Studies (M8 steps 5 and 6, criteria 3 and 4): every point's outputs, or
// why one has none; a failed point among the rest; the order points are
// solved in and the chunks they start cold at; and every warm point against
// its cold solve, with what each took, which can't grow.
module study

import flowsheet
import json2
import math
import os
import reactors

fn flowsheet_text(name string) string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', name)) or { panic(err) }
}

fn load(text string) !flowsheet.Loaded {
	return flowsheet.load_text(text, makers: reactors.makers())
}

// Tutorial 2's reactor over its feed's T: CO's conversion at each point is
// the correlation's, and methanol's has no value, never −∞, saying why.
fn test_tutorial_2s_conversions() {
	l := load(flowsheet_text('methanol_study.toml'))!
	r := run(l, l.studies[0])
	assert r.points.len == 4
	for k, p in r.points {
		t_c := 200.0 + 40.0 * f64(k)
		assert p.status == .no_value && p.error == ''
		assert math.abs(p.outputs[0] - (0.85 - 0.0018 * t_c)) <= 1e-12, '${p.outputs[0]}'
		assert math.is_nan(p.outputs[1])
		assert p.reasons[1] == 'it divides by syngas.n[CH4O], which is zero here'
		assert p.reasons[0] == '' && p.reasons[2] == ''
		assert math.abs(p.outputs[2] - t_c) <= 1e-12
		assert p.warm == (k > 0)
	}
	rep := r.report('methanol_study.toml', [])
	j := json2.decode[json2.Any](rep.json())!.as_map()
	assert j['format']!.int() == 1 && j['status']!.str() == 'ran'
	point := j['points']!.as_array()[0].as_map()
	values := point['outputs']!.as_array()
	assert values[1].as_map()['value']! is json2.Null
	assert values[1].as_map()['reason']!.str().starts_with('it divides by')
	lines := rep.csv().split_into_lines()
	assert lines[0] == 'point,feed.syngas.T [K],CO_conversion,methanol_conversion,T_out,status'
	assert lines[1].contains(',NaN,') && lines[1].ends_with(',no_value')
	assert rep.table().contains('methanol_conversion, at every point: it divides by syngas.n[CH4O], which is zero here')
}

// documented is every key the study report's section of
// docs/report-format.md lists, by the object it's in: '' for the study,
// 'points[]' for a point.
fn documented() map[string][]string {
	text := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md')) or { panic(err) }
	section := text.all_after('## The study report').all_before('
## ')
	mut out := map[string][]string{}
	mut path := ''
	for line in section.split_into_lines() {
		if line == '#### The study' {
			path = ''
		} else if line.starts_with('#### `') {
			path = line.all_after('`').all_before('`')
		} else if line.starts_with('| `') {
			out[path] << line.all_after('| `').all_before('`')
		}
	}
	return out
}

// The JSON's keys are the format's, in order, and nothing more.
fn test_the_json_is_as_documented() {
	doc := documented()
	l := load(flowsheet_text('methanol_study.toml'))!
	j := json2.decode[json2.Any](run(l, l.studies[0]).report('x', []).json())!.as_map()
	assert j.keys() == doc[''], '${j.keys()}'
	assert j['points']!.as_array()[0].as_map().keys() == doc['points[]']
	assert doc.keys().len == 2
}

// The order (D7): along one input in order; along two, row by row,
// alternating direction, so that each point's neighbour is the last solved.
fn test_the_order_points_are_solved_in() {
	one := flowsheet.StudyPlan{
		axes: [flowsheet.Axis{
			values: [1.0, 2, 3]
		}]
	}
	assert order(one) == [0, 1, 2]
	two := flowsheet.StudyPlan{
		axes: [flowsheet.Axis{
			values: [1.0, 2, 3]
		}, flowsheet.Axis{
			values: [10.0, 20]
		}]
	}
	assert order(two) == [0, 1, 3, 2, 4, 5]
	assert values(two, 3) == [2.0, 20]
	// The predictor's ratio: along one input, one way, or none.
	assert ratio([1.0, 10], [2.0, 10], [4.0, 10])? == 2.0
	assert ratio([1.0, 10], [2.0, 10], [2.0, 20]) == none // a row's turn
	assert ratio([2.0, 10], [2.0, 20], [3.0, 20]) == none
	assert ratio([1.0], [2.0], [1.0]) == none // back the way it came
	assert ratio([1.0], [1.0], [2.0]) == none
}

// apart is how far two points' outputs are, each as a fraction of the
// larger, or of 1 for a pure number near zero.
fn apart(a Point, b Point) f64 {
	mut worst := 0.0
	for k, x in a.outputs {
		y := b.outputs[k]
		worst = math.max(worst, math.abs(x - y) / math.max(1.0, math.max(math.abs(x), math.abs(y))))
	}
	return worst
}

struct Case {
	name   string
	text   string
	key    string
	values string
	report string
	// How far a warm point may be from its cold one: the solves' own
	// tolerances, measured by spike/warm.
	within f64
	// The warm study's passes, spec evaluations and integrator steps,
	// measured, which can't grow.
	total []int
}

// Criterion 4: every warm point against its cold solve, point by point, and
// none taking more passes, spec evaluations or integrator steps than it.
fn test_warm_points_against_cold() {
	reference := flowsheet_text('reference.toml')
	reformer := flowsheet_text('reformer_dwsim.toml')
	pfr := reformer[reformer.index('[unit.reformer]') or { 0 }..reformer.index('[site.Ni]') or {
		0
	}]
	cases := [
		Case{'the reference flowsheet, by direct substitution', reference.replace('method = "wegstein"',
			'method = "direct"'), 'unit.cooler.outlet_T', 'from = "230 K", to = "210 K", points = 11', 'gas = "gas.F", product = "product.F", cold = "cold.T"', 1e-8, [
			256,
			0,
			0,
		]},
		Case{'the reference flowsheet, by direct substitution, unevenly spaced', reference.replace('method = "wegstein"',
			'method = "direct"'), 'unit.cooler.outlet_T', 'values = ["230 K", "229 K", "227 K", "223 K", "215 K", "214 K", "212 K", "208 K"]', 'gas = "gas.F", cold = "cold.T"', 1e-8, [
			185,
			0,
			0,
		]},
		Case{'the reference flowsheet, by Wegstein', reference, 'unit.cooler.outlet_T', 'from = "230 K", to = "210 K", points = 11', 'gas = "gas.F", duty = "cooler.duty"', 0, [
			44,
			0,
			0,
		]},
		Case{"the owner's PFR", reformer, 'feed.feed.T', 'from = "700 degC", to = "1000 degC", points = 9', 'methane = "product.n[CH4]", hydrogen = "product.y[H2]", T = "product.T"', 1e-9, [
			0,
			0,
			5101,
		]},
		Case{'a CSTR', reformer.replace(pfr, '[unit.reformer]\ntype = "cstr"\nin = "feed"\nout = "product"\nvolume = "1 m^3"\ncatalyst_loading = "0.386 kg/m^3"\n\n'), 'feed.feed.T', 'from = "700 degC", to = "1000 degC", points = 9', 'methane = "product.n[CH4]", duty = "reformer.duty"', 1e-11, [
			0,
			0,
			0,
		]},
		Case{'the reference flowsheet with a spec in its loop', reference +
			'\n[spec.vapour]\nvary = "unit.cooler.outlet_T"\nbetween = ["200 K", "260 K"]\nuntil = "gas.F = 60 mol/s"\n', 'feed.feed.P', 'from = "30 bar", to = "31 bar", points = 9', 'cold = "cold.T", product = "product.F"', 1e-8, [
			134,
			134,
			0,
		]},
		Case{'an adiabatic Gibbs reactor', reformer.replace(pfr, '[unit.reformer]\ntype = "gibbs_reactor"\nin = "feed"\nout = "product"\nadiabatic = true\n\n'), 'feed.feed.T', 'from = "700 degC", to = "1000 degC", points = 9', 'methane = "product.n[CH4]", T = "product.T"', 1e-10, [
			0,
			0,
			0,
		]},
	]
	for c in cases {
		text := c.text.replace('format = 1', 'format = 3').replace('format = 2', 'format = 3') +
			'\n[study.s]\nvary = [{ key = "${c.key}", ${c.values} }]\nreport = { ${c.report} }\n'
		l := load(text)!
		warm := run(l, l.studies[0])
		cold := run(l, l.studies[0], warm: false)
		mut total := [0, 0, 0]
		for p in warm.points {
			total[0] += p.passes
			total[1] += p.evaluations
			total[2] += p.steps
		}
		assert total[0] <= c.total[0] && total[1] <= c.total[1] && total[2] <= c.total[2], '${c.name}: ${total}'
		for k, p in warm.points {
			q := cold.points[k]
			assert p.status == .solved && q.status == .solved, '${c.name} ${k}: ${p.error}${q.error}'
			assert p.warm == (k % 8 != 0) && !q.warm, c.name
			assert apart(p, q) <= c.within, '${c.name} ${k}: ${apart(p, q)}'
			assert p.passes <= q.passes, '${c.name} ${k}: ${p.passes} passes, cold ${q.passes}'
			assert p.evaluations <= q.evaluations, '${c.name} ${k}: ${p.evaluations} evaluations, cold ${q.evaluations}'
			assert p.steps <= q.steps, '${c.name} ${k}: ${p.steps} steps, cold ${q.steps}'
		}
	}
}

// A heater's spec over its feed's T: each warm point starts from its
// neighbour's duty, and meets its target in fewer evaluations than cold.
fn test_a_spec_over_a_study() {
	text := 'format = 3\n\n[components]\nnames = ["methane", "propane", "n-pentane"]\n\n[model]\nequation = "PR"\n\n[feed.f]\nT = "300 K"\nP = "10 bar"\nflow = "100 mol/s"\nmole_fractions = { methane = 0.6, propane = 0.3, n-pentane = 0.1 }\n\n[unit.heat]\ntype = "heater"\nin = "f"\nout = "out"\nduty = "0 W"\n\n[spec.cool]\nvary = "unit.heat.duty"\nbetween = ["-2 MW", "1 MW"]\nuntil = "out.T = 250 K"\n\n[study.feed]\nvary = [{ key = "feed.f.T", from = "280 K", to = "320 K", points = 21 }]\nreport = { duty = "heat.duty", T = "out.T" }\n'
	l := load(text)!
	warm := run(l, l.studies[0], warm: true)
	cold := run(l, l.studies[0], warm: false)
	mut w, mut c := 0, 0
	for k, p in warm.points {
		q := cold.points[k]
		assert math.abs(p.outputs[1] - 250.0) <= 1e-9 * 250.0
		assert math.abs(p.specs[0] - q.specs[0]) <= 1e-8 * math.abs(q.specs[0]), '${k}'
		assert p.evaluations <= q.evaluations, '${k}: ${p.evaluations} evaluations, cold ${q.evaluations}'
		w += p.evaluations
		c += q.evaluations
	}
	// Measured: 108 evaluations warm against 168 cold, which can't grow.
	assert w <= 108 && c <= 168, '${w} ${c}'
	// The same value again: met where the last point left it, in one
	// evaluation.
	again := load(text.replace('from = "280 K", to = "320 K", points = 21', 'values = ["300 K", "300 K"]'))!
	twice := run(again, again.studies[0])
	assert twice.points[1].warm && twice.points[1].evaluations == 1
}

// A point that fails doesn't stop the study (D9): the cooler's T too low
// for its spec's target fails, its neighbours stand, and the next point
// after it starts cold.
fn test_a_failed_point_among_the_rest() {
	text := flowsheet_text('reference.toml').replace('format = 1', 'format = 3') +
		'\n[study.feed]\nvary = [{ key = "feed.feed.flow", values = ["100 mol/s", "100 mol/s", "-1 mol/s", "100 mol/s", "100 mol/s"] }]\nreport = { gas = "gas.F" }\n'
	l := load(text)!
	r := run(l, l.studies[0])
	assert r.points.map(it.status) == [Status.solved, .solved, .failed, .solved, .solved]
	assert r.points[2].error.contains('flow'), r.points[2].error
	assert math.is_nan(r.points[2].outputs[0])
	assert r.points.map(it.warm) == [false, true, false, false, true]
	assert r.points[3].outputs[0] == r.points[0].outputs[0]
	rep := r.report('x.toml', [])
	assert rep.table().contains('Failed points:\n  point 2 (feed.feed.flow = -1.0000 mol/s): ')
}

fn refused(extra string, line int, want string) {
	text := flowsheet_text('reference.toml').replace('format = 1', 'format = 3') + '\n' + extra
	load(text) or {
		assert err is flowsheet.FlowsheetError, err.msg()
		e := err as flowsheet.FlowsheetError
		assert e.reason.contains(want), err.msg()
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		return
	}
	assert false, 'it loaded: ${want}'
}

fn test_studies_refused() {
	base := '[study.s]\nvary = [{ key = "unit.cooler.outlet_T", from = "220 K", to = "230 K", points = 3 }]\nreport = { gas = "gas.F" }\n'
	// Evenly spaced values end at to exactly, which 0.2 + 0.7·1 doesn't.
	even := load(flowsheet_text('reference.toml').replace('format = 1', 'format = 3') +
		'\n[study.e]\nvary = [{ key = "unit.pump.efficiency", from = 0.2, to = 0.9, points = 2 }]\nreport = { gas = "gas.F" }\n')!
	assert even.studies[0].axes[0].values == [0.2, 0.9]
	assert 0.2 + (0.9 - 0.2) * 1.0 != 0.9
	assert load(flowsheet_text('reference.toml').replace('format = 1', 'format = 3') + '\n' +
		base)!.studies.len == 1
	refused(base.replace('outlet_T', 'duty'), 65, 'unit.cooler has no duty in the file')
	refused(base.replace('from = "220 K", to = "230 K", points = 3', 'from = "220 K", to = "230 K"'),
		65, 'give from, to and points, or values')
	refused(base.replace('points = 3', 'points = 1'), 65, 'from and to need 2 points or more')
	refused(base.replace('"220 K"', '220'), 65, 'a value of unit.cooler.outlet_T needs its unit')
	refused(base.replace('"220 K"', '"220 bar"'), 65, 'bar is a unit of pressure, not temperature')
	refused(base.replace('from = "220 K", to = "230 K", points = 3', 'values = []'), 65, 'it needs a list of values')
	refused(base.replace('vary = [', 'vary = [{ key = "unit.cooler.outlet_T", values = ["220 K"] }, '),
		65, 'unit.cooler.outlet_T is varied twice')
	refused(base.replace('gas.F', 'gs.F'), 66, 'column 1: no stream or unit gs here')
	refused(base.replace('{ gas = "gas.F" }', '{}'), 64, 'a study reports at least one expression')
	refused(base + 'chunk = 0\n', 67, 'a chunk of 0 points')
	refused(base + 'workers = 2\n', 67, 'no key workers here')
	refused(base + '[spec.t]\nvary = "unit.cooler.outlet_T"\nbetween = ["200 K", "260 K"]\nuntil = "gas.F = 60 mol/s"\n',
		65, 'spec t varies unit.cooler.outlet_T, to meet its target, so a study')
	load(flowsheet_text('reference.toml') + '\n' + base) or {
		assert err.msg().contains('studies need format 3'), err.msg()
		return
	}
	assert false
}
