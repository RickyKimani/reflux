// The tape (M9 step 3, D5, criterion 3): the tree's own bits, exact
// derivatives against SymPy and central differences, several outputs on
// one tape, and nothing allocated by an evaluation.
import expr
import json2
import math
import os
import quantity

struct Case {
	text      string
	points    [][]f64
	values    []f64
	scales    []f64
	gradients [][]f64
}

struct Fixture {
	cases []Case
}

fn load() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_expressions.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

// The fixture's context, as tools/expr_render.v has it.
fn fixture_context() expr.Context {
	pure := quantity.Dim{}
	return expr.Context{
		variables: ['a', 'b', 'x', 'y'].map(expr.Variable{it, pure})
		families:  [
			expr.Family{'p', pure, ['CH4', 'H2O'], {
				'CH4': 0
				'H2O': 1
			}},
		]
	}
}

// rate_context is a rate law's: T, P and each species' partial pressure.
fn rate_context() expr.Context {
	species := ['CH4', 'H2O', 'H2', 'CO', 'CO2']
	mut keys := map[string]int{}
	for i, s in species {
		keys[s] = i
	}
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, species, keys},
		]
	}
}

// Xu and Froment's three rates (AIChE J. 35, 88, 1989) as a user writes
// them: their constants and their one adsorption term as definitions.
fn xu_froment() !(expr.Definitions, []expr.Expr) {
	defs := expr.define(rate_context(), [
		expr.Source{'k1', '4.22e15 kmol*bar^0.5/(kg*h) * exp(-240.1 kJ/mol / (R*T))'},
		expr.Source{'k2', '1.96e6 kmol/(kg*h*bar) * exp(-67.13 kJ/mol / (R*T))'},
		expr.Source{'k3', '1.02e15 kmol*bar^0.5/(kg*h) * exp(-243.9 kJ/mol / (R*T))'},
		expr.Source{'K_CO', '8.23e-5 bar^-1 * exp(70.65 kJ/mol / (R*T))'},
		expr.Source{'K_H2', '6.12e-9 bar^-1 * exp(82.9 kJ/mol / (R*T))'},
		expr.Source{'K_CH4', '6.65e-4 bar^-1 * exp(38.28 kJ/mol / (R*T))'},
		expr.Source{'K_H2O', '1.77e5 * exp(-88.68 kJ/mol / (R*T))'},
		expr.Source{'K1', 'exp(30.42 - 27106 K / T) * 1 bar^2'},
		expr.Source{'K2', 'exp(-3.798 + 4160 K / T)'},
		expr.Source{'K3', 'exp(34.218 - 31266 K / T) * 1 bar^2'},
		expr.Source{'DEN', '1 + K_CO*p[CO] + K_H2*p[H2] + K_CH4*p[CH4] + K_H2O*p[H2O]/p[H2]'},
	])!
	ctx := rate_context()
	mut rates := []expr.Expr{}
	for text in [
		'k1 / p[H2]^2.5 * (p[CH4]*p[H2O] - p[H2]^3*p[CO]/K1) / DEN^2',
		'k2 / p[H2] * (p[CO]*p[H2O] - p[H2]*p[CO2]/K2) / DEN^2',
		'k3 / p[H2]^3.5 * (p[CH4]*p[H2O]^2 - p[H2]^4*p[CO2]/K3) / DEN^2',
	] {
		rates << expr.parse(ctx, &defs, text)!
	}
	return defs, rates
}

// A reforming state: T, P, and each partial pressure, in Pa.
const reforming = [900.0, 1.5e6, 2.5e5, 5.0e5, 3.0e5, 1.0e5, 0.8e5]

fn test_the_tape_gives_the_trees_bits() {
	// Every fixture case at every point: the tape's value is eval's, to the
	// bit, since every fusion keeps the tree's arithmetic.
	f := load()!
	ctx := fixture_context()
	defs := expr.Definitions{}
	for c in f.cases {
		e := expr.parse(ctx, &defs, c.text)!
		t := expr.compile(ctx, &defs, [e])
		mut buf := t.buffer()
		for x in c.points {
			t.eval(x, mut buf)
			assert math.f64_bits(t.output(&buf, 0)) == math.f64_bits(e.eval(&defs, x)), '${c.text} at ${x}: ${t.output(&buf, 0)} vs ${e.eval(&defs, x)}'
		}
	}
	// And Xu and Froment's rates, through their definitions.
	d, rates := xu_froment()!
	t := expr.compile(rate_context(), &d, rates)
	mut buf := t.buffer()
	t.eval(reforming, mut buf)
	for k, r in rates {
		assert math.f64_bits(t.output(&buf, k)) == math.f64_bits(r.eval(&d, reforming)), 'rate ${k + 1}'
	}
}

fn test_the_gradient_is_sympys() {
	// Every fixture case's derivative with respect to every input, at every
	// point, against SymPy's at 50 digits: to 1e-10 of the larger of the
	// derivative, the value's scale and 1. A wrong derivative's error is of
	// the derivative's own size.
	f := load()!
	ctx := fixture_context()
	defs := expr.Definitions{}
	mut worst := 0.0
	for c in f.cases {
		e := expr.parse(ctx, &defs, c.text)!
		t := expr.compile(ctx, &defs, [e])
		mut buf := t.buffer()
		mut g := []f64{len: t.inputs}
		for k, x in c.points {
			t.eval(x, mut buf)
			t.gradient(mut buf, mut g)
			for j, want in c.gradients[k] {
				err := math.abs(g[j] - want) / math.max(math.max(1.0, math.abs(want)), c.scales[k])
				worst = math.max(worst, err)
				assert err <= 1e-10, "${c.text} at ${x}, d/d input ${j}: ${g[j]} vs SymPy's ${want}"
			}
		}
	}
	println('worst gradient error against SymPy: ${worst:.2e}')
}

// central is the central difference of output k in input j, h relative.
fn central(t &expr.Tape, x []f64, k int, j int, rel f64) f64 {
	mut xp := x.clone()
	mut xm := x.clone()
	h := rel * math.max(1.0, math.abs(x[j]))
	xp[j] += h
	xm[j] -= h
	mut buf := t.buffer()
	t.eval(xp, mut buf)
	up := t.output(&buf, k)
	t.eval(xm, mut buf)
	return (up - t.output(&buf, k)) / (xp[j] - xm[j])
}

fn test_the_gradient_matches_central_differences() {
	// Xu and Froment's rates with respect to T, P and every partial
	// pressure: h = 1e-5 relative, whose truncation, h²·f‴/6, is of order
	// 1e-9 of a derivative here, and rounding, ε·f/h, of order 1e-11.
	d, rates := xu_froment()!
	t := expr.compile(rate_context(), &d, rates)
	mut buf := t.buffer()
	t.eval(reforming, mut buf)
	mut g := []f64{len: rates.len * t.inputs}
	t.gradient(mut buf, mut g)
	for k in 0 .. rates.len {
		// Measured against the largest |∂r/∂x_j · x_j| in the row, the terms
		// the rate's change is made of.
		mut row := 0.0
		for j in 0 .. t.inputs {
			row = math.max(row, math.abs(g[k * t.inputs + j] * reforming[j]))
		}
		for j in 0 .. t.inputs {
			fd := central(&t, reforming, k, j, 1e-5)
			got := g[k * t.inputs + j]
			assert math.abs(got - fd) * math.abs(reforming[j]) <= 1e-7 * row, 'rate ${k + 1}, input ${j}: ${got} vs ${fd}'
		}
	}
	// P is no rate's input, so its derivatives are 0.
	for k in 0 .. rates.len {
		assert g[k * t.inputs + 1] == 0.0
	}
}

fn test_several_outputs_share_one_tape() {
	// The three rates on one tape share the adsorption term, R·T and the
	// constants: fewer operations than three tapes, the same values to the
	// bit, and the same gradients to rounding (an adjoint's sum may come in
	// another order).
	d, rates := xu_froment()!
	ctx := rate_context()
	shared := expr.compile(ctx, &d, rates)
	mut alone := 0
	mut buf := shared.buffer()
	shared.eval(reforming, mut buf)
	mut g := []f64{len: rates.len * shared.inputs}
	shared.gradient(mut buf, mut g)
	for k, r in rates {
		t := expr.compile(ctx, &d, [r])
		alone += t.instrs.len
		mut b := t.buffer()
		t.eval(reforming, mut b)
		assert math.f64_bits(t.output(&b, 0)) == math.f64_bits(shared.output(&buf, k)), 'rate ${k + 1}'
		mut gk := []f64{len: t.inputs}
		t.gradient(mut b, mut gk)
		for j in 0 .. t.inputs {
			assert math.abs(gk[j] - g[k * t.inputs + j]) <= 1e-14 * math.max(math.abs(gk[j]), 1e-300), 'rate ${k + 1}, input ${j}'
		}
	}
	assert shared.instrs.len < alone
	println("Xu and Froment's three rates: ${shared.instrs.len} operations on one tape, ${alone} on three; ${shared.merged} merged, ${shared.folded} folded, ${shared.fused} fused")
	// Every Arrhenius and van 't Hoff term is one operation, and the
	// constants are all fused away.
	assert shared.instrs.filter(it.op == .aexp).len == 7
	assert shared.instrs.filter(it.op == .expl).len == 3
	assert shared.n_const == 0
}

fn test_a_fixed_expression_has_no_gradient() {
	ctx := rate_context()
	defs := expr.Definitions{}
	e := expr.parse(ctx, &defs, '2 * R * 300 K')!
	t := expr.compile(ctx, &defs, [e])
	assert t.first == t.n_const && t.n_const == 1 // one constant, no inputs, no operations
	mut buf := t.buffer()
	t.eval(reforming, mut buf)
	assert t.output(&buf, 0) == 2.0 * 8.31446261815324 * 300.0
	mut g := []f64{len: t.inputs, init: 1.0}
	t.gradient(mut buf, mut g)
	assert g.all(it == 0.0)
}

fn test_an_evaluation_allocates_nothing() {
	// 100,000 evaluations and gradients of the three rates, into one buffer.
	d, rates := xu_froment()!
	t := expr.compile(rate_context(), &d, rates)
	mut buf := t.buffer()
	mut g := []f64{len: rates.len * t.inputs}
	mut x := reforming.clone()
	before := gc_heap_usage().total_bytes
	mut sum := 0.0
	for i in 0 .. 100000 {
		x[0] = 800.0 + f64(i % 500)
		t.eval(x, mut buf)
		t.gradient(mut buf, mut g)
		sum += t.output(&buf, 0) + g[0]
	}
	allocated := gc_heap_usage().total_bytes - before
	assert allocated == 0, '${allocated} bytes'
	assert math.is_finite(sum)
}

fn test_a_definition_is_compiled_once() {
	// 60 definitions, each using the next twice: compiled once each, not
	// 2⁶⁰ times, into an input, a quotient and 59 sums.
	mut chain := []expr.Source{}
	for i in 0 .. 59 {
		chain << expr.Source{'d${i}', 'd${i + 1} + d${i + 1}'}
	}
	chain << expr.Source{'d59', 'T / (1 K)'}
	ctx := rate_context()
	defs := expr.define(ctx, chain)!
	e := expr.parse(ctx, &defs, 'd0')!
	t := expr.compile(ctx, &defs, [e])
	assert t.instrs.len == 1 + 1 + 59
	mut buf := t.buffer()
	t.eval(reforming, mut buf)
	assert t.output(&buf, 0) == 900.0 * math.pow(2.0, 59)
}

fn test_whole_powers_are_products() {
	// A whole power is a repeated product, the same bits on every platform,
	// in the tree and on the tape: 0.01⁴, where this machine's math.pow
	// gives 1e-8, and a negative power as 1 over that.
	pure := quantity.Dim{}
	ctx := expr.Context{
		variables: [expr.Variable{'x', pure}]
	}
	defs := expr.Definitions{}
	for text, want in {
		'x^4':  0.01 * 0.01 * 0.01 * 0.01
		'x^-2': 1.0 / (0.01 * 0.01)
	} {
		e := expr.parse(ctx, &defs, text)!
		t := expr.compile(ctx, &defs, [e])
		mut buf := t.buffer()
		t.eval([0.01], mut buf)
		assert math.f64_bits(e.eval(&defs, [0.01])) == math.f64_bits(want), text
		assert math.f64_bits(t.output(&buf, 0)) == math.f64_bits(want), text
	}
}

fn test_a_scaled_exponential_isnt_scaled_again() {
	// 4·exp(d/T) is one operation; times 2 it stays (4·exp)·2, the tree's
	// arithmetic, rather than becoming 8·exp, which rounds otherwise.
	ctx := rate_context()
	defs := expr.Definitions{}
	for text in ['4 * exp(1000 K / T) * 2', '3 * (5 * exp(-1000 K / T))', '0.1 * exp(27 K / T) * 0.3'] {
		e := expr.parse(ctx, &defs, text)!
		t := expr.compile(ctx, &defs, [e])
		mut buf := t.buffer()
		for temp in [300.0, 777.7, 1234.5] {
			mut x := reforming.clone()
			x[0] = temp
			t.eval(x, mut buf)
			assert math.f64_bits(t.output(&buf, 0)) == math.f64_bits(e.eval(&defs, x)), '${text} at ${temp} K'
		}
	}
}
