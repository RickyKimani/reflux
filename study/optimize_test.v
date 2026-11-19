// Optimising a flowsheet (M8 step 10, D13, criterion 8): series reactions
// in a PFR, whose optima are known in closed form, inside the bounds, on a
// bound, on an inequality and on an equality; each verified from the
// closed form's own derivatives, as criterion 7 verifies.
module study

import flowsheet
import json2
import math
import os
import reactors

const series_file = 'format = 3

[components]
names = ["n-pentane", "isopentane", "neopentane"]

[model]
equation = "PR"

[feed.f]
T = "400 K"
P = "1 bar"
flow = "1 mol/s"
mole_fractions = { n-pentane = 1, isopentane = 0, neopentane = 0 }

[reaction.first]
equation = "n-pentane = isopentane"
law = "power_law"
k = { A = "70 mol/(m^3*s*Pa)", Ea = "50 kJ/mol" }
orders = { n-pentane = 1 }
K = "none"
per = "volume"

[reaction.second]
equation = "isopentane = neopentane"
law = "power_law"
k = { A = "0.09 mol/(m^3*s*Pa)", Ea = "30 kJ/mol" }
orders = { isopentane = 1 }
K = "none"
per = "volume"

[unit.tube]
type = "pfr"
in = "f"
out = "product"
volume = "1 m^3"
rtol = 1e-12

[optimize.interior]
vary = [{ key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] }]
maximize = "product.n[isopentane]"

[optimize.bound]
vary = [{ key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] }, { key = "feed.f.T", between = ["380 K", "420 K"] }]
maximize = "product.n[isopentane]"

[optimize.constraint]
vary = [{ key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] }]
maximize = "product.n[isopentane]"
subject_to = ["product.n[neopentane] <= 0.05 mol/s"]

[optimize.equality]
vary = [{ key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] }]
minimize = "product.n[neopentane]"
subject_to = ["product.n[\'n-pentane\'] = 0.5 mol/s"]

[optimize.offset]
vary = [{ key = \"unit.tube.volume\", between = [\"0.1 m^3\", \"2 m^3\"] }]
maximize = \"product.n[isopentane] + 100 mol/s\"

[optimize.scaled]
vary = [{ key = \"unit.tube.volume\", between = [\"0.1 m^3\", \"2 m^3\"] }]
maximize = \"product.n[isopentane]\"
subject_to = [\"1e9 * product.n[neopentane] <= 5e7 mol/s\"]

[optimize.none]
vary = [{ key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] }]
maximize = "product.n[isopentane]"
subject_to = ["product.n[neopentane] >= 2 mol/s"]
'

fn load_series() !flowsheet.Loaded {
	return flowsheet.load_text(series_file, makers: reactors.makers())
}

// Series is the closed form: K₁ and K₂ at T, per m³ of volume, and the
// outflows of A, B and C at a volume, with B's derivative by V.
struct Series {
	k1 f64
	k2 f64
}

fn series(t f64) Series {
	r := 1.380649e-23 * 6.02214076e23
	return Series{
		k1: 70.0 * math.exp(-50000.0 / (r * t)) * 1e5
		k2: 0.09 * math.exp(-30000.0 / (r * t)) * 1e5
	}
}

fn (s Series) b(v f64) f64 {
	return s.k1 / (s.k2 - s.k1) * (math.exp(-s.k1 * v) - math.exp(-s.k2 * v))
}

fn (s Series) db(v f64) f64 {
	return s.k1 / (s.k2 - s.k1) * (-s.k1 * math.exp(-s.k1 * v) + s.k2 * math.exp(-s.k2 * v))
}

fn (s Series) dc(v f64) f64 {
	// C = 1 − A − B, A = exp(−K₁V).
	return s.k1 * math.exp(-s.k1 * v) - s.db(v)
}

fn series_fixture(name string) map[string]json2.Any {
	fixture := json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata',
		'm8_series.json')) or { panic(err) }) or { panic(err) }.as_map()
	return fixture[name] or { panic(name) }.as_map()
}

fn plan_of(l flowsheet.Loaded, name string) flowsheet.OptimizePlan {
	for p in l.optimizations {
		if p.name == name {
			return p
		}
	}
	panic('no optimisation ${name}')
}

// The solves each optimisation took, measured, which can't grow: one an
// evaluation, five an input for each gradient, the start's gradient taken
// once for the scales and the SQP's first both.
const optimize_solves = {
	'interior':   45
	'bound':      112
	'constraint': 52
	'equality':   37
}

// The iterations each optimisation took, measured, which can't grow.
const optimize_iterations = {
	'interior':   7
	'bound':      10
	'constraint': 8
	'equality':   6
}

fn ran(l flowsheet.Loaded, name string) Optimum {
	o := optimize(l, plan_of(l, name))
	assert o.ok, '${name}: ${o.error}'
	println('${name}: ${o.iterations.len} iterations, ${o.evaluations} evaluations, ${o.solves} solves')
	most := optimize_iterations[name] or { 0 }
	assert o.iterations.len <= most, '${name}: ${o.iterations.len} iterations, recorded ${most}'
	solves := optimize_solves[name] or { 0 }
	assert o.solves <= solves, '${name}: ${o.solves} solves, recorded ${solves}'
	assert o.iterations.last().step == 0.0 && o.iterations.last().inputs == o.inputs
	return o
}

// B over V at 400 K, inside V's bounds: V* = ln(K₂/K₁)/(K₂ − K₁). dB/dV is
// zero there to the SQP's tolerance, 1e-6 of the gradient's scale at the
// start, which leaves V good to about 1e-6 of its range over B's curvature.
fn test_an_optimum_inside_the_bounds() {
	l := load_series()!
	o := ran(l, 'interior')
	want := series_fixture('interior')
	s := series(400.0)
	v := o.inputs[0]
	assert math.abs(v - want['V']!.f64()) <= 1e-5 * want['V']!.f64(), '${v}'
	assert math.abs(o.objective - want['B']!.f64()) <= 1e-11, '${o.objective}'
	// From scratch: B's own derivative at the answer, beside dB/dV's scale,
	// K₁·B.
	assert math.abs(s.db(v)) <= 1e-5 * s.k1 * s.b(v), '${s.db(v)}'
	assert o.lower == [0.0] && o.upper == [0.0]
	assert math.abs(o.objective - s.b(v)) <= 1e-10
}

// B over V and T: B's peak rises with T, so T sits on its upper bound,
// exactly, and V at V*(420 K). The bound's multiplier is how fast B's
// optimum rises with that bound, ∂B/∂T at (V*, 420 K) by the envelope
// theorem, and B's own ∂B/∂V is zero there.
fn test_an_optimum_on_a_bound() {
	l := load_series()!
	o := ran(l, 'bound')
	want := series_fixture('bound')
	assert o.inputs[1] == 420.0
	v := o.inputs[0]
	assert math.abs(v - want['V']!.f64()) <= 1e-5 * want['V']!.f64(), '${v}'
	assert math.abs(o.objective - want['B']!.f64()) <= 1e-11, '${o.objective}'
	s := series(420.0)
	assert math.abs(s.db(v)) <= 1e-5 * s.k1 * s.b(v)
	// ∂B/∂T by the closed form, central differences over 1e-3 K, good to
	// about 1e-9 of itself.
	dbdt := (series(420.001).b(v) - series(419.999).b(v)) / 0.002
	assert math.abs(dbdt - want['T_multiplier']!.f64()) <= 1e-6 * dbdt
	assert math.abs(o.upper[1] - dbdt) <= 1e-5 * dbdt, '${o.upper[1]} against ${dbdt}'
	assert o.lower == [0.0, 0.0] && o.upper[0] == 0.0
}

// B over V with C's outflow at most 0.05 mol/s, which C reaches before B
// peaks: V where C = 0.05, and the constraint's multiplier, how fast B's
// optimum rises as C's limit is relaxed, dB/dV over dC/dV there.
fn test_an_optimum_on_a_constraint() {
	l := load_series()!
	o := ran(l, 'constraint')
	want := series_fixture('constraint')
	v := o.inputs[0]
	c := o.constraints[0]
	assert c.active && math.abs(c.lhs - 0.05) <= 1e-6 * 0.05 && c.rhs == 0.05, '${c}'
	assert math.abs(v - want['V']!.f64()) <= 1e-5 * want['V']!.f64(), '${v}'
	assert math.abs(o.objective - want['B']!.f64()) <= 1e-6 * want['B']!.f64()
	s := series(400.0)
	// From scratch: dB/dV = λ·dC/dV, λ ≥ 0.
	assert c.multiplier > 0.0
	assert math.abs(s.db(v) - c.multiplier * s.dc(v)) <= 1e-5 * s.db(v)
	assert math.abs(c.multiplier - want['multiplier']!.f64()) <= 1e-5 * want['multiplier']!.f64(), '${c.multiplier}'
	assert o.lower == [0.0] && o.upper == [0.0]
}

// An equality: A's outflow 0.5 mol/s fixes V = ln 2/K₁, whatever the
// objective, C's outflow minimised. The multiplier is how fast C's optimum
// improves, falls, as the equality's right side grows: V = −ln(Fₐ/F)/K₁, so
// dC/dFₐ = −C′(V)/(K₁·Fₐ), and the improvement is its negative.
fn test_an_optimum_on_an_equality() {
	l := load_series()!
	o := ran(l, 'equality')
	s := series(400.0)
	v := math.ln2 / s.k1
	assert math.abs(o.inputs[0] - v) <= 1e-9 * v, '${o.inputs[0]} against ${v}'
	c := o.constraints[0]
	// Met to the SQP's tolerance, 1e-6 of its scale, 0.5 mol/s.
	assert c.active && math.abs(c.lhs - 0.5) <= 5e-7, '${c.lhs}'
	improvement := s.dc(v) / (s.k1 * 0.5)
	assert math.abs(c.multiplier - improvement) <= 1e-5 * math.abs(improvement), '${c.multiplier} against ${improvement}'
}

// Neither an objective's size nor a constraint's moves the optimum: the
// interior one with 100 mol/s added to B, and the constraint's with both
// sides times 1e9, whose multiplier is so 1e9 times smaller. Each is scaled
// at the start, the objective by its derivative, which the 100 mol/s
// doesn't change, and the constraint by its sides.
fn test_scales_dont_move_an_optimum() {
	l := load_series()!
	offset := optimize(l, plan_of(l, 'offset'))
	inside := optimize(l, plan_of(l, 'interior'))
	assert offset.ok && inside.ok
	v := series_fixture('interior')['V']!.f64()
	assert math.abs(offset.inputs[0] - v) <= 1e-5 * v, '${offset.inputs[0]}'
	assert math.abs(offset.objective - 100.0 - inside.objective) <= 1e-12
	scaled := optimize(l, plan_of(l, 'scaled'))
	limited := optimize(l, plan_of(l, 'constraint'))
	assert scaled.ok && limited.ok, scaled.error
	assert math.abs(scaled.inputs[0] - limited.inputs[0]) <= 1e-9 * limited.inputs[0]
	m := limited.constraints[0].multiplier
	assert math.abs(scaled.constraints[0].multiplier * 1e9 - m) <= 1e-6 * m
}

// The step and its side: about x inside the bounds, above at the lower, below
// at the upper, and at most a quarter of a narrow range.
fn test_the_step_and_its_side() {
	h, inside := step_of(0.75, 0.5, 1.0)
	assert h == 1e-4 && inside == .both
	_, low := step_of(0.5, 0.5, 1.0)
	_, high := step_of(1.0, 0.5, 1.0)
	_, near := step_of(0.50005, 0.5, 1.0)
	assert low == .above && high == .below && near == .above
	narrow, _ := step_of(1000.0, 1000.0, 1000.0001)
	assert narrow == 0.25 * (1000.0001 - 1000.0)
	zero, _ := step_of(0.0, 0.0, 2.0)
	assert zero == 2e-4
}

// No point has more neopentane than the feed has of anything.
fn test_no_feasible_point() {
	l := load_series()!
	o := optimize(l, plan_of(l, 'none'))
	assert !o.ok && o.error.contains('no feasible point'), o.error
	assert math.is_nan(o.objective)
}

// Two inputs' gradients on two workers give the same bits as on one.
fn test_workers_dont_change_an_optimum() {
	l := load_series()!
	one := optimize(l, plan_of(l, 'bound'))
	two := optimize(l, plan_of(l, 'bound'), workers: 2)
	assert one.ok && two.ok
	assert one.inputs == two.inputs && one.objective == two.objective
	assert one.upper == two.upper && one.iterations.len == two.iterations.len
}

fn refused_optimisation(extra string, line int, want string) {
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml')) or {
		panic(err)
	}.replace('format = 1', 'format = 3') + '
' + extra
	flowsheet.load_text(text, makers: reactors.makers()) or {
		assert err is flowsheet.FlowsheetError, err.msg()
		e := err as flowsheet.FlowsheetError
		assert e.reason.contains(want), err.msg()
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		return
	}
	assert false, 'it loaded: ${want}'
}

// [optimize.NAME] (D13, criterion 12): every key checked, every refusal at
// its line. An empty list has no place in V's TOML tree, so its refusal is
// at its table's header.
fn test_optimisations_refused() {
	base := '[optimize.o]
vary = [{ key = "unit.cooler.outlet_T", between = ["200 K", "260 K"] }]
minimize = "gas.F"
subject_to = ["gas.T <= 300 K"]
'
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml'))!.replace('format = 1',
		'format = 3') + '
' + base
	ok := flowsheet.load_text(text, makers: reactors.makers())!
	p := ok.optimizations[0]
	assert p.name == 'o' && !p.maximize && p.vary[0].lo == 200.0 && p.vary[0].hi == 260.0
	assert p.vary[0].start == 230.0 && p.constraints[0].relation == .le && p.line == 64
	started := flowsheet.load_text(text.replace('"260 K"] }', '"260 K"], start = "250 K" }'),
		makers: reactors.makers()
	)!
	assert started.optimizations[0].vary[0].start == 250.0
	refused_optimisation(base.replace('[{ key = "unit.cooler.outlet_T", between = ["200 K", "260 K"] }]',
		'[]'), 64, 'an optimisation varies at least one input')
	refused_optimisation(base.replace('outlet_T', 'duty'), 65, 'unit.cooler has no duty in the file')
	refused_optimisation(base.replace('["200 K", "260 K"]', '["200 K"]'), 65, 'it needs two values')
	refused_optimisation(base.replace('["200 K", "260 K"]', '["260 K", "200 K"]'), 65,
		'its least value must be below its most')
	refused_optimisation(base.replace('"260 K"] }', '"260 K"], start = "300 K" }'), 65,
		'its start is outside it')
	refused_optimisation(base.replace('["200 K", "260 K"]', '["100 K", "150 K"]'), 65,
		"the file's unit.cooler.outlet_T, 230 K, is outside it")
	refused_optimisation(base.replace('"260 K"] }', '"260 K"], step = 1 }'), 65, 'no key step here')
	refused_optimisation(base.replace('"260 K"] }]', '"260 K"] }, { key = "unit.cooler.outlet_T", between = ["200 K", "260 K"] }]'),
		65, 'unit.cooler.outlet_T is varied twice')
	refused_optimisation(base + 'maximize = "gas.F"
', 68, 'not both')
	refused_optimisation(base.replace('minimize = "gas.F"
', ''), 64, 'an optimisation needs an objective')
	refused_optimisation(base.replace('"gas.F"', '"gas.Q"'), 66, 'column 1: stream gas has no Q')
	refused_optimisation(base.replace('gas.T <= 300 K', 'gas.F <= 300 K'), 67, 'must have the same units')
	refused_optimisation(base.replace('gas.T <= 300 K', 'gas.T'), 67, 'an equation needs =, <= or >=')
	refused_optimisation(base + 'chunk = 8
', 68, 'no key chunk here')
	spec := '[spec.v]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"

'
	refused_optimisation(spec + base, 70, "spec v varies unit.cooler.outlet_T, to meet its target, so an optimisation can't")
	two := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml'))! + '
' + base
	flowsheet.load_text(two, makers: reactors.makers()) or {
		assert err.msg().contains('optimisations need format 3'), err.msg()
		return
	}
	assert false, 'format 1 took an optimisation'
}

// optimize_documented is each object's keys in the optimisation report's
// section of docs/report-format.md, by its path.
fn optimize_documented() map[string][]string {
	text := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md')) or { panic(err) }
	section := text.all_after('## The optimisation report').all_before('
## ')
	mut out := map[string][]string{}
	mut path := ''
	for line in section.split_into_lines() {
		if line == '#### The optimisation' {
			path = ''
		} else if line.starts_with('#### `') {
			path = line.all_after('`').all_before('`')
		} else if line.starts_with('| `') {
			out[path] << line.all_after('| `').all_before('`')
		}
	}
	return out
}

// The JSON's keys are the format's, in order, and nothing more, with an
// optimum and without.
fn test_the_json_is_as_documented() {
	doc := optimize_documented()
	l := load_series()!
	for name in ['constraint', 'none'] {
		j := json2.decode[json2.Any](optimize(l, plan_of(l, name)).report('x', []).json())!.as_map()
		assert j.keys() == doc[''], '${j.keys()}'
		assert j['inputs']!.as_array()[0].as_map().keys() == doc['inputs[]']
		assert j['constraints']!.as_array()[0].as_map().keys() == doc['constraints[]']
		assert j['iterations']!.as_array()[0].as_map().keys() == doc['iterations[]']
		assert j['objective']!.as_map().keys() == ['expression', 'unit', 'value']
	}
	assert doc.keys().len == 4
}

// A pump's work is W_s/η exactly, so the least is at η = 1, a bound the
// pump's efficiency can't step past: its maker refuses more than 1, so the
// gradient there must be one-sided. The bound's multiplier is how fast the
// work falls as it's relaxed, −dW/dη = W_s/η² = W at η = 1.
fn test_a_bound_the_input_cant_pass() {
	text := 'format = 3

[components]
names = ["n-pentane"]

[model]
equation = "PR"

[feed.f]
T = "300 K"
P = "1 bar"
flow = "1 mol/s"
mole_fractions = { n-pentane = 1 }

[unit.pump]
type = "pump"
in = "f"
out = "out"
outlet_P = "10 bar"
efficiency = 0.7

[optimize.work]
vary = [{ key = "unit.pump.efficiency", between = [0.5, 1] }]
minimize = "pump.work"
'
	l := flowsheet.load_text(text, makers: reactors.makers())!
	o := optimize(l, l.optimizations[0])
	assert o.ok, o.error
	assert o.inputs == [1.0]
	w := o.objective
	assert math.abs(o.upper[0] - w) <= 1e-6 * w, '${o.upper[0]} against ${w}'
	assert o.lower == [0.0]
	println('the pump: ${o.iterations.len} iterations, ${o.solves} solves')
	// Measured: 2 iterations and 13 solves, which can't grow.
	assert o.iterations.len <= 2 && o.solves <= 13
}

// x_of puts an input on its bound exactly at 0 and 1, which lo + u·(hi − lo)
// doesn't: 0.2 + 1·(0.9 − 0.2) isn't 0.9.
fn test_a_bound_exactly() {
	s := Sheet{
		plan: flowsheet.OptimizePlan{
			vary: [flowsheet.Varied{
				lo: 0.2
				hi: 0.9
			}]
		}
	}
	assert 0.2 + 1.0 * (0.9 - 0.2) != 0.9
	assert s.x_of([1.0]) == [0.9] && s.x_of([0.0]) == [0.2]
	assert s.x_of([0.5])[0] == 0.2 + 0.5 * (0.9 - 0.2)
}
