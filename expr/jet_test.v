// The jet pass (M6 step 2): an expression's value with its first two
// derivatives along one input. Its value is eval's, to the bit; its first
// derivatives are SymPy's (testdata/m9_expressions.json), as the tape's
// are; its second derivatives are central differences of the first, which
// are exact, and closed forms where the expression has one.
import expr
import json2
import math
import numeric
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

// seeded is x as jets along input j.
fn seeded(x []f64, j int) []numeric.Jet {
	return []numeric.Jet{len: x.len, init: if index == j {
		numeric.along(x[index])
	} else {
		numeric.constant(x[index])
	}}
}

fn test_the_jet_is_evals_and_sympys() {
	f := load()!
	ctx := fixture_context()
	defs := expr.Definitions{}
	mut worst_d, mut worst_dd := 0.0, 0.0
	mut checked := 0
	for c in f.cases {
		e := expr.parse(ctx, &defs, c.text)!
		for k, x in c.points {
			for j, want in c.gradients[k] {
				jet := e.jet(&defs, seeded(x, j))
				assert math.f64_bits(jet.v) == math.f64_bits(e.eval(&defs, x)), '${c.text} at ${x}'
				err := math.abs(jet.d - want) / math.max(math.max(1.0, math.abs(want)), c.scales[k])
				worst_d = math.max(worst_d, err)
				// The second derivative against a central difference of the
				// first, step h relative to the input: truncation h²·f⁗/6,
				// rounding ε·f′/h.
				h := 1e-5 * math.max(math.abs(x[j]), 1.0)
				mut up := x.clone()
				mut down := x.clone()
				up[j] += h
				down[j] -= h
				d2 := (e.jet(&defs, seeded(up, j)).d - e.jet(&defs, seeded(down, j)).d) / (2.0 * h)
				if math.is_finite(d2) && math.is_finite(jet.dd) {
					scale := math.max(math.max(1.0, math.abs(jet.dd)), math.max(math.abs(jet.d), c.scales[k]))
					worst_dd = math.max(worst_dd, math.abs(jet.dd - d2) / scale)
				}
				checked++
			}
		}
	}
	println('${checked} derivatives: first against SymPy ${worst_d:.1e}, second against differences ${worst_dd:.1e}')
	// Measured: 6.5e-16 over 6,240 derivatives, and 3.9e-7. The first held as
	// the tape's is, to 1e-10 of its scale (tape_test.v); the second is the
	// differences' own error, their truncation h²·f⁗/6 where an expression's
	// fourth derivative is large.
	assert worst_d <= 1e-10
	assert worst_dd <= 2e-6
}

fn test_closed_forms_in_t() {
	ctx := expr.Context{
		variables: [expr.Variable{'T', quantity.Dimension.temperature.dim()}]
	}
	defs := expr.Definitions{}
	t := 333.15
	at := [numeric.along(t)]
	// DECHEMA's NRTL τ, A/(R·T) with its own R: τ′ = −τ/T, τ″ = 2τ/T².
	tau := expr.parse(ctx, &defs, '-189.0469 K / (1.98721*T)')!.jet(&defs, at)
	v := -189.0469 / (1.98721 * t)
	assert math.abs(tau.v - v) <= 1e-16 * math.abs(v)
	assert math.abs(tau.d + v / t) <= 1e-15 * math.abs(v / t)
	assert math.abs(tau.dd - 2.0 * v / (t * t)) <= 1e-15 * math.abs(2.0 * v / (t * t))
	// A reading on a unit's scale moves as T does: (T in degC)² has 2(T −
	// 273.15) and 2.
	r := expr.parse(ctx, &defs, '(T in degC)^2')!.jet(&defs, at)
	assert math.abs(r.d - 2.0 * (t - 273.15)) <= 1e-12 && math.abs(r.dd - 2.0) <= 1e-15
	// A whole power below the line, and ln: (T/(1 K))^-2 and ln(T/(1 K)).
	p := expr.parse(ctx, &defs, '(T / (1 K))^-2')!.jet(&defs, at)
	assert math.abs(p.dd - 6.0 / (t * t * t * t)) <= 1e-15 * 6.0 / (t * t * t * t)
	l := expr.parse(ctx, &defs, 'ln(T / (1 K))')!.jet(&defs, at)
	assert math.abs(l.d - 1.0 / t) <= 1e-18 && math.abs(l.dd + 1.0 / (t * t)) <= 1e-20
	// And through a definition.
	d := expr.define(ctx, [expr.Source{'k', 'exp(-1000 K / T)'}])!
	k := expr.parse(ctx, &d, '2 * k')!.jet(&d, at)
	e := math.exp(-1000.0 / t)
	assert math.abs(k.d - 2.0 * e * 1000.0 / (t * t)) <= 1e-15 * 2.0 * e * 1000.0 / (t * t)
}
