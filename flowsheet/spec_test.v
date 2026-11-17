// M8 steps 2 to 4 on M4's units (criteria 1, 2 and 12): a flowsheet's
// variables in expressions and their values in a solution; every number a
// file gives remade by its maker; and specs met, checked from the solution,
// against closed forms and against solves driven from outside.
import components
import flash
import flowsheet { FlowsheetError, Loaded, Solution }
import math
import numeric
import os
import thermo

const db = components.load_default() or { panic(err) }

fn flowsheet_text(name string) string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', name)) or { panic(err) }
}

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

// with_spec is the reference flowsheet in format 3, with these tables after.
fn with_spec(extra string) string {
	return flowsheet_text('reference.toml').replace('format = 1', 'format = 3') + '\n' + extra
}

// heater is a feed, a heater at a duty and its product, in format 3.
fn heater(extra string) string {
	return 'format = 3

[components]
names = ["methane", "propane", "n-pentane"]

[model]
equation = "PR"

[feed.f]
T = "300 K"
P = "10 bar"
flow = "100 mol/s"
mole_fractions = { methane = 0.6, propane = 0.3, n-pentane = 0.1 }

[unit.heat]
type = "heater"
in = "f"
out = "out"
duty = "0 W"
' +
		extra
}

fn refused(text string, line int, want string) FlowsheetError {
	load(text) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.reason.contains(want), err.msg()
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		return e
	}
	assert false, 'it loaded: ${want}'
	return FlowsheetError{}
}

fn test_variables_bind_and_take_a_solutions_values() {
	l := load(flowsheet_text('reference.toml'))!
	sol := l.solve()!
	for name, st in sol.streams {
		quoted := if name.contains(' ') { "'${name}'" } else { name }
		for v, want in {
			'T':  st.t
			'P':  st.p
			'F':  st.total()
			'h':  st.h
			'vf': st.beta
		} {
			e := l.expression('${quoted}.${v}')!
			assert l.value(e, sol)! == want, '${name}.${v}'
		}
		mut sum := 0.0
		for i, c in l.names {
			n := l.value(l.expression("${quoted}.n['${c}']")!, sol)!
			y := l.value(l.expression("${quoted}.y['${c}']")!, sol)!
			assert n == st.flows[i], '${name} ${c}'
			assert y == st.flows[i] / st.total(), '${name} ${c}'
			sum += y
		}
		assert math.abs(sum - 1.0) <= 4e-16, name
	}
	// Each unit's energy stream by its name, and none for those with none.
	for i, p in l.sheet.units {
		for e in sol.units[i].energy {
			assert l.value(l.expression('${p.unit.name}.${e.name}')!, sol)! == e.rate
		}
	}
	for name in ['mixer', 'split', 'valve'] {
		l.expression('${name}.duty') or {
			assert err.msg().contains('unit ${name} has no duty: it has no variables'), err.msg()
			continue
		}
		assert false, name
	}
	l.expression('pump.duty') or {
		assert err.msg().contains('unit pump has no duty: it has work'), err.msg()
	}
	// A member by formula or by the file's name, and a quoted stream.
	assert l.expression('gas.y[CH4]')!.plain() == 'gas.y[CH4]'
	assert l.expression("product.n['n-butane']")!.plain() == 'product.n[C4H10]'
	assert l.value(l.expression("'let down'.P")!, sol)! == 5e5
	// The vapour's methane by the liquid's, and the cooler's duty per mole.
	assert l.expression('gas.y[CH4] / liquid.y[methane]')!.dim.is_dimensionless()
	assert l.expression('cooler.duty / mixed.F')!.unicode() == 'cooler.duty/mixed.F'
}

fn test_a_value_that_isnt_there_says_why() {
	l := load(flowsheet_text('reference.toml'))!
	sol := l.solve()!
	// A stream emptied by hand: its fractions have no value, and a ratio of
	// its flows divides by zero.
	mut streams := sol.streams.clone()
	gas := streams['gas'] or { panic('no gas') }
	streams['gas'] = flowsheet.Stream{
		...gas
		flows: []f64{len: gas.flows.len}
	}
	empty := Solution{
		...sol
		streams: streams
	}
	l.value(l.expression('1 - gas.y[CH4]')!, empty) or {
		assert err.msg().contains('column 5: gas.y[CH4] has no value here: gas has no flow'), err.msg()
	}
	l.value(l.expression('(gas.n[CH4] - liquid.n[CH4]) / gas.n[CH4]')!, empty) or {
		assert err.msg().contains('it divides by gas.n[CH4], which is zero here'), err.msg()
		return
	}
	assert false
}

fn test_every_number_can_be_varied() {
	for name in ['reference.toml', 'every_unit.toml', 'no_recycle.toml'] {
		l := load(flowsheet_text(name))!
		sol := l.solve()!
		mut count := 0
		for i, p in l.sheet.units {
			for r in l.readings_of('unit', p.unit.name) {
				input := l.input(r.key)!
				assert input.value == r.value && input.index == i, r.key
				// The same value again: the same unit, to the bit on its inlets.
				u := l.with([input], [input.value])!.sheet.units[i].unit
				o := u.calculate(l.fluid, sol.units[i].inlets)!
				for k, x in o.outlets {
					y := sol.units[i].outlets[k]
					assert x.t == y.t && x.p == y.p && x.h == y.h && x.flows == y.flows, '${name} ${r.key}'
				}
				assert o.energy == sol.units[i].energy, '${name} ${r.key}'
				count++
			}
		}
		for feed in l.sheet.feeds {
			for r in l.readings_of('feed', feed.name) {
				input := l.input(r.key)!
				assert input.value == r.value
				again := l.with([input], [input.value])!
				assert again.sheet.feeds.map(it.stream.t) == l.sheet.feeds.map(it.stream.t), r.key
				count++
			}
		}
		assert count > 10, '${name}: ${count}'
	}
}

fn test_a_remade_value_is_what_the_maker_reads() {
	l := load(flowsheet_text('reference.toml'))!
	t := l.input('unit.cooler.outlet_T')!
	assert t.value == 230.0 && t.text == '230 K' && t.line == 29
	sol := l.with([t], [220.0])!.solve()!
	assert sol.streams['cold'].t == 220.0
	// The file it was made from is as it was.
	assert l.input('unit.cooler.outlet_T')!.value == 230.0 && l.write().contains('outlet_T = "230 K"')
	// A feed's T and P, and a pump's efficiency, a plain number.
	feed_t := l.input('feed.feed.T')!
	feed_p := l.input('feed.feed.P')!
	remade := l.with([feed_t, feed_p], [310.0, 31e5])!
	assert remade.sheet.feeds[0].stream.t == 310.0 && remade.sheet.feeds[0].stream.p == 31e5
	eff := l.input('unit.pump.efficiency')!
	assert eff.dim.is_dimensionless() && eff.value == 0.75
	l.with([eff], [0.7])!
	// The splitter's fractions, both at once; one alone is refused by the
	// splitter, at its line, as the file would be.
	a := l.input('unit.split.fractions[0]')!
	b := l.input('unit.split.fractions[1]')!
	split := l.with([a, b], [0.4, 0.6])!.solve()!
	assert math.abs(split.streams['back'].total() / split.streams['liquid'].total() - 0.4) < 1e-12
	l.with([a], [0.4]) or {
		assert err is FlowsheetError && (err as FlowsheetError).line == 38, err.msg()
	}
	// A temperature given in °C: varied in K, the maker's reading.
	e := load(flowsheet_text('every_unit.toml'))!
	oil := e.input('feed.oil.T')!
	assert oil.value == 300.15 && oil.text == '27 degC'
	assert e.with([oil], [301.0])!.sheet.feeds[1].stream.t == 301.0
	// The writing keeps the new value.
	assert l.with([t], [219.5])!.write().contains('outlet_T = "219.5 K"')
}

fn test_keys_that_cant_be_varied() {
	l := load(flowsheet_text('reference.toml'))!
	for key, want in {
		'unit.coolerx.outlet_T':   'no unit coolerx'
		'stream.gas.T':            'a key to vary starts with feed. or unit.'
		'unit.coolr.outlet_T':     'no unit coolr'
		'unit.cooler.duty':        'unit.cooler has no duty in the file: it gives outlet_T'
		'unit.cooler.type':        "it isn't a number the unit reads as one"
		'unit.mixer.in':           "it isn't a number the unit reads as one"
		'feed.feed.flow.x':        'feed.feed has no flow.x in the file'
		'feed.fed.T':              'no feed fed'
		'unit.split.fractions[2]': 'unit.split has no fractions[2] in the file'
	} {
		l.input(key) or {
			assert err.msg().contains(want), '${key}: ${err.msg()}'
			continue
		}
		assert false, key
	}
	l.input('unit.coolr.outlet_T') or {
		assert (err as FlowsheetError).suggestions == ['unit.cooler.outlet_T']
	}
}

// The heater's duty for an outlet T (criterion 1): a block of one unit and
// one spec, met by Brent, against the duty a PT flash at exactly that T
// gives.
fn test_a_heaters_duty_for_its_outlet_t() {
	l := load(heater('
[spec.cool]
vary = "unit.heat.duty"
between = ["-2 MW", "1 MW"]
until = "out.T = 250 K"
'))!
	assert l.plan.blocks.len == 1 && l.plan.blocks[0].specs == [0]
	sol := l.solve()!
	r := sol.specs[0]
	assert r.name == 'cool' && r.key == 'unit.heat.duty' && r.rhs == 250.0
	assert math.abs(r.lhs - 250.0) <= 1e-9 * 250.0 && math.abs(r.residual) <= 1e-9
	// The exact duty, from the feed and a flash at 250 K.
	f := l.sheet.feeds[0].stream
	at := flash.pt(l.fluid.model, 250.0, f.p, f.flows)!
	mut h := 0.0
	for ph in at.phases {
		h += ph.fraction * l.fluid.props(ph.state).h
	}
	exact := f.total() * (h - f.h)
	// T is met to tol·T, which moves the duty by F·Cp·ΔT, Cp under 200 J/(mol·K).
	assert math.abs(r.value - exact) <= f.total() * 200.0 * math.abs(r.lhs - 250.0) + 1e-9 * math.abs(exact), '${r.value} ${exact}'
	assert sol.streams['out'].t == r.lhs
	// The heater reports F·Δh, which rounds the duty it was given.
	assert math.abs(sol.units[0].energy[0].rate - r.value) <= 1e-12 * math.abs(r.value)
	// Brent's evaluations, measured: 8, which can't grow.
	assert sol.history[0].spec_residuals.len <= 8, '${sol.history[0].spec_residuals.len}'
	// The solve remade the heater in its own lists: the file's is at 0 W.
	plain := l.sheet.units[0].unit.calculate(l.fluid, [f])!
	assert math.abs(plain.energy[0].rate) <= 1e-9 * math.abs(exact)
	// Its report.
	rep := l.report(sol)
	assert rep.specs.len == 1 && rep.specs[0].until == 'out.T = 250 K' && rep.specs[0].value == r.value
	assert rep.blocks[0].specs == ['cool'] && rep.blocks[0].spec_residuals.len > 2
}

// A feed's flow for a heater's duty: linear, so exactly F = Q/(h_out − h_in).
fn test_a_feeds_flow_for_a_duty() {
	text := heater('
[spec.flow]
vary = "feed.f.flow"
between = ["10 mol/s", "1000 mol/s"]
until = "heat.duty = 3 MW"
').replace('duty = "0 W"', 'outlet_T = "400 K"')
	l := load(text)!
	sol := l.solve()!
	r := sol.specs[0]
	out := sol.streams['out']
	f := sol.streams['f']
	assert math.abs(r.value - 3e6 / (out.h - f.h)) <= 1e-9 * r.value, '${r.value}'
	assert math.abs(f.total() - r.value) <= 4e-16 * r.value // the flows sum to it
	// Measured: 10 evaluations, which can't grow.
	assert sol.history[0].spec_residuals.len <= 10, '${sol.history[0].spec_residuals.len}'
}

// Two specs in one block with no tears, each varying what the other
// measures: Broyden on their values.
fn test_two_specs_in_one_block() {
	text := heater('
[unit.trim]
type = "heater"
in = "out"
out = "trimmed"
duty = "0 W"

[spec.a]
vary = "unit.heat.duty"
between = ["-1 MW", "3 MW"]
until = "trimmed.T = 330 K"

[spec.b]
vary = "unit.trim.duty"
between = ["-1 MW", "1 MW"]
until = "out.T = 350 K"
')
	l := load(text)!
	assert l.plan.blocks.len == 1 && l.plan.blocks[0].specs == [0, 1]
	sol := l.solve()!
	assert math.abs(sol.streams['out'].t - 350.0) <= 1e-9 * 350.0
	assert math.abs(sol.streams['trimmed'].t - 330.0) <= 1e-9 * 330.0
	assert sol.specs.map(it.name) == ['a', 'b']
	// Broyden's evaluations, measured: 9, which can't grow.
	assert sol.history[0].spec_residuals.len <= 9, '${sol.history[0].spec_residuals.len}'
}

// A target the bounds' ends don't bracket, though values between them do:
// (T − 300 K)² is 400 K² at 280 and 320 K, and larger at both ends.
fn test_a_bracket_inside_the_bounds() {
	l := load(heater('
[spec.square]
vary = "unit.heat.duty"
between = ["-2 MW", "2 MW"]
until = "(out.T - 300 K)^2 = 400 K^2"
'))!
	sol := l.solve()!
	t := sol.streams['out'].t
	assert math.abs(math.abs(t - 300.0) - 20.0) <= 1e-6, '${t}'
	// Measured: 13 evaluations, nine of them the bracket's, which can't grow.
	assert sol.history[0].spec_residuals.len <= 13, '${sol.history[0].spec_residuals.len}'
}

// A target of zero: the residual is a fraction of F·RT, the block's own
// scale for a power, since the duty itself goes to zero. The heater's duty
// is zero where its outlet is the feed's T.
fn test_a_target_of_zero() {
	text := heater('
[spec.none]
vary = "unit.heat.outlet_T"
between = ["250 K", "350 K"]
until = "heat.duty = 0 W"
').replace('duty = "0 W"', 'outlet_T = "320 K"')
	sol := load(text)!.solve()!
	r := sol.specs[0]
	f := sol.streams['f']
	assert math.abs(r.lhs) <= 1e-9 * f.total() * 8.31446261815324 * 300.0, '${r.lhs}'
	assert math.abs(r.value - 300.0) <= 1e-6, '${r.value}'
	// Measured: 3 evaluations, which can't grow.
	assert sol.history[0].spec_residuals.len <= 3, '${sol.history[0].spec_residuals.len}'
}

// A tolerance below rounding: Brent's bracket shrinks to the last bit, and
// the check from the solution refuses what's left rather than passing it.
fn test_a_target_tighter_than_rounding_is_refused() {
	text := heater('
[solver]
tol = 1e-18

[spec.cool]
vary = "unit.heat.duty"
between = ["-2 MW", "1 MW"]
until = "out.T = 250 K"
')
	load(text)!.solve() or {
		assert err.msg().contains('spec.cool: at unit.heat.duty = '), err.msg()
		assert err.msg().contains('out.T = 250 K is off by'), err.msg()
		return
	}
	assert false
}

// An equation with no value at a point says why, never ±∞: a heater makes
// no methane, so its conversion divides by zero.
fn test_an_equation_with_no_value() {
	l := load(heater('
[spec.x]
vary = "unit.heat.duty"
between = ["-1 MW", "1 MW"]
until = "(f.n[methane] - out.n[methane]) / (out.n[methane] - f.n[methane]) = 1"
'))!
	l.solve() or {
		assert err.msg().contains('spec.x: at unit.heat.duty = -1000000 W'), err.msg()
		assert err.msg().contains('its equation has no value: it divides by out.n[CH4] - f.n[CH4], which is zero here'), err.msg()
		return
	}
	assert false
}

// The reference flowsheet with a spec in its recycle's block: D18's, the
// cooler's outlet T until the first drum's vapour is 60 mol/s. Against
// Brent driving whole solves from outside, each remade by with.
fn test_a_spec_in_a_recycle() {
	l := load(with_spec('
[spec.vapour]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"
'))!
	b := l.plan.blocks.filter(it.specs.len > 0)
	assert b.len == 1 && b[0].tears == ['recycle'] && b[0].specs == [0]
	sol := l.solve()!
	r := sol.specs[0]
	assert math.abs(sol.streams['gas'].total() - 60.0) <= 1e-9 * 60.0
	// From outside: Brent on whole solves.
	plain := load(flowsheet_text('reference.toml'))!
	input := plain.input('unit.cooler.outlet_T')!
	outside := numeric.brent(Outside{plain, input}, 200.0, 260.0, xtol: 1e-12)!
	println('spec ${r.value} K, outside ${outside.x} K; ${sol.history.map(it.residuals.len)} passes, ${sol.history.map(it.spec_residuals.len)} evaluations')
	assert math.abs(r.value - outside.x) <= 1e-6, '${r.value} ${outside.x}'
	// Broyden on the tear and the spec together, measured: 27 passes, which
	// can't grow.
	assert sol.history[0].residuals.len <= 27
}

struct Outside {
	l     Loaded
	input flowsheet.Input
}

fn (o Outside) at(t f64) !f64 {
	sol := o.l.with([o.input], [t])!.solve()!
	return sol.streams['gas'].total() - 60.0
}

fn test_a_target_out_of_reach() {
	l := load(with_spec('
[spec.vapour]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 200 mol/s"
'))!
	l.solve() or {
		assert err is FlowsheetError
		e := err as FlowsheetError
		assert e.kind == .spec_not_met && e.key == 'spec.vapour', err.msg()
		assert e.reason.contains('its target is out of reach: with unit.cooler.outlet_T from 200 K to 260 K, gas.F goes from'), err.msg()
		assert e.reason.contains('never 200 mol/s'), err.msg()
		// Broyden held at the bound stops after three tries, and Brent's
		// nested bracket takes ten values: measured, 15 evaluations in all.
		assert e.history[0].spec_residuals.len <= 15, '${e.history[0].spec_residuals.len}'
		return
	}
	assert false
}

fn test_specs_refused() {
	spec := '
[spec.s]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"
'
	refused(flowsheet_text('reference.toml') + spec, 65, 'specs need format 3')
	refused(with_spec(spec.replace('outlet_T"', 'duty"')), 66, 'unit.cooler has no duty in the file')
	refused(with_spec(spec.replace('["200 K", "260 K"]', '["200 K"]')), 67, 'it needs two values')
	refused(with_spec(spec.replace('["200 K", "260 K"]', '[200, 260]')), 67, 'a bound of unit.cooler.outlet_T needs its unit')
	refused(with_spec(spec.replace('["200 K", "260 K"]', '["260 K", "200 K"]')), 67, 'its least value must be below its most')
	refused(with_spec(spec.replace('["200 K", "260 K"]', '["200 K", "220 K"]')), 67, "the file's unit.cooler.outlet_T, 230 K, is outside it")
	refused(with_spec(spec.replace('["200 K", "260 K"]', '["200 bar", "260 K"]')), 67, 'bar is a unit of pressure, not temperature')
	e := refused(with_spec(spec.replace('gas.F', 'gs.F')), 68, 'column 1: no stream or unit gs here')
	assert e.suggestions[0] == 'gas.F'
	refused(with_spec(spec.replace('gas.F = 60 mol/s', 'gas.F <= 60 mol/s')), 68, "a spec's until is an equation, with =")
	refused(with_spec(spec.replace('gas.F = 60 mol/s', 'gas.F = 60 K')), 68, 'the two sides of = must have the same units')
	refused(with_spec(spec + 'tol = 1\n'), 69, 'no key tol here')
	// Two specs on one key, one reading only a feed, and one whose target
	// nothing it varies reaches: the second drum is after the recycle.
	refused(with_spec(spec + spec.replace('spec.s', 'spec.t')), 71, 'spec s varies unit.cooler.outlet_T already')
	refused(with_spec(spec.replace('gas.F', 'feed.F')), 65, 'it reads only feeds')
	refused(with_spec(spec.replace('unit.cooler.outlet_T', 'unit.second.duty').replace('["200 K", "260 K"]',
		'["-1 MW", "1 MW"]')), 65, 'nothing it measures depends on what it varies')
}

// Hot is a unit that can't start cold: it passes its inlet through only
// from a neighbour's answer, to show that a warm solve asks its units for
// calculate_from (M8 step 5, D6).
struct Hot {
	name string
}

fn (h Hot) inlet_ports() []string {
	return ['in']
}

fn (h Hot) outlet_ports() []string {
	return ['out']
}

fn (h Hot) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	return error('Hot starts only from a neighbour')
}

fn (h Hot) calculate_from(fl thermo.Fluid, inlets []flowsheet.Stream, last flowsheet.UnitResult) !flowsheet.Outcome {
	assert last.outlets.len == 1
	return flowsheet.Outcome{
		outlets: inlets
	}
}

fn test_a_warm_solve_starts_its_units_from_their_last_answers() {
	l := load(flowsheet_text('no_recycle.toml'))!
	sol := l.solve()!
	mut sheet := flowsheet.Flowsheet{
		feeds: l.sheet.feeds
	}
	feed := l.sheet.feeds[0].name
	sheet.add_unit(Hot{'hot'}, {
		'in':  feed
		'out': 'through'
	})!
	sheet.solve(l.fluid) or { assert err.msg().contains('Hot starts only from a neighbour') }
	last := flowsheet.Solution{
		...sol
		units: [flowsheet.UnitResult{
			name:    'hot'
			inlets:  [sol.streams[feed]]
			outlets: [sol.streams[feed]]
		}]
	}
	loaded := Loaded{
		...l
		sheet: sheet
		plan:  sheet.plan()!
		specs: []
	}
	warm := loaded.solve_from(flowsheet.Start{
		last: last
	})!
	assert warm.streams['through'].flows == sol.streams[feed].flows
	// Started without its units, it's cold again.
	loaded.solve_from(flowsheet.Start{ last: last, units: false }) or {
		assert err.msg().contains('Hot starts only from a neighbour')
		return
	}
	assert false
}
