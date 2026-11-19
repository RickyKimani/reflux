// A flowsheet's names in expressions, the values they take and the
// equations they make (M8 steps 2 and 4, criteria 2 and 3): dotted names
// bind with their dimensions or are refused with their place and a
// suggestion; each renders back; an expression with no value at a point
// says why, never ±∞; and an equation's two sides share a dimension.
import expr { Kind }
import math
import quantity
import rand

// sheet is a flowsheet's context as M8's D3 has it: streams' T, P, F and
// flows, a unit's duty, and a stream whose name has a space in it.
fn sheet() expr.Context {
	members := ['CH4', 'H2O', 'CH3OH']
	mut keys := map[string]int{}
	for i, m in members {
		keys[m] = i
	}
	keys['methane'] = 0
	keys['n-butane'] = 1 // a key that must be quoted
	t := quantity.Dimension.temperature.dim()
	p := quantity.Dimension.pressure.dim()
	f := quantity.Dimension.molar_flow.dim()
	return expr.Context{
		variables: [
			expr.Variable{'feed.T', t},
			expr.Variable{'feed.P', p},
			expr.Variable{'feed.F', f},
			expr.Variable{'product.T', t},
			expr.Variable{'product.F', f},
			expr.Variable{'let down.T', t},
			expr.Variable{'s1.2.T', t}, // a scope with a dot in its name
			expr.Variable{'reactor.duty', quantity.Dimension.power.dim()},
		]
		families:  [
			expr.Family{'feed.n', f, members, keys},
			expr.Family{'product.n', f, members, keys},
			expr.Family{'product.y', quantity.Dim{}, members, keys},
			expr.Family{'let down.y', quantity.Dim{}, members, keys},
		]
		scopes:    [
			expr.Scope{'feed', 'stream'},
			expr.Scope{'product', 'stream'},
			expr.Scope{'let down', 'stream'},
			expr.Scope{'s1.2', 'stream'},
			expr.Scope{'reactor', 'unit'},
			expr.Scope{'mixer', 'unit'},
		]
	}
}

const none_defined = expr.Definitions{}

fn on_sheet(text string) !expr.Expr {
	return expr.parse(sheet(), &none_defined, text)
}

fn refused_on_sheet(text string, column int, want string) expr.ExprError {
	on_sheet(text) or {
		assert err is expr.ExprError, err.msg()
		e := err as expr.ExprError
		assert e.reason.contains(want), '${text}: ${err.msg()}'
		assert e.line == 1 && e.column == column, '${text}: column ${e.column}, not ${column}: ${err.msg()}'
		return e
	}
	assert false, '"${text}" was read'
	return expr.ExprError{}
}

fn test_dotted_names_bind_with_their_dimensions() {
	ctx := sheet()
	// Each variable is its slot, each member after the variables.
	for i, v in ctx.variables {
		text := if v.name.contains(' ') || v.name.count('.') > 1 {
			"'${v.name.all_before_last('.')}'.${v.name.all_after_last('.')}"
		} else {
			v.name
		}
		e := on_sheet(text)!
		assert e.nodes.len == 1 && e.nodes[0].kind == .var && e.nodes[0].index == i, text
		assert e.dim.eq(v.dim), text
	}
	m := on_sheet('product.n[methane]')!
	assert m.nodes[0].kind == .member && m.nodes[0].index == ctx.slot(1, 0)
	assert m.nodes[0].text == 'CH4' // as the context names the member
	assert m.dim.eq(quantity.Dimension.molar_flow.dim())
	q := on_sheet('\'let down\'.y["n-butane"]')!
	assert q.nodes[0].index == ctx.slot(3, 1) && q.dim.is_dimensionless()
	// A product's conversion: a pure number.
	x := on_sheet('(feed.n[CH3OH] - product.n[CH3OH]) / feed.n[CH3OH]')!
	assert x.dim.is_dimensionless()
	assert x.nodes.map(it.kind) == [Kind.member, .member, .sub, .member, .div]
	// A duty over a flow is a molar enthalpy; a temperature difference is
	// one of temperature.
	assert on_sheet('reactor.duty / feed.F')!.dim.eq(quantity.Dimension.molar_enthalpy.dim())
	assert on_sheet('product.T - feed.T')!.dim.eq(quantity.Dimension.temperature.dim())
	// A number with a unit beside a dotted name, and a reading of one.
	assert on_sheet('feed.T - 300 K')!.nodes.len == 3
	assert on_sheet('(feed.T in degC)')!.dim.is_dimensionless()
}

fn test_dotted_names_refused() {
	// A sum of two dimensions, as M9's.
	refused_on_sheet('product.T - 300', 11, 'the two sides of − must have the same units')
	// A scope that isn't there: the scopes like it, with the variable.
	e := refused_on_sheet('prodct.T', 1, 'no stream or unit prodct here')
	assert e.suggestions == ['product.T'], '${e.suggestions}'
	// A variable the scope hasn't: what it has, and the closest.
	f := refused_on_sheet('product.TT', 1, 'stream product has no TT: it has T, F, n[…], y[…]')
	assert f.suggestions[0] == 'product.T', '${f.suggestions}'
	refused_on_sheet("'let down'.P", 1, "stream 'let down' has no P: it has T, y[…]")
	refused_on_sheet("'s1.2'.P", 1, "stream 's1.2' has no P: it has T")
	refused_on_sheet('reactor.work', 1, 'unit reactor has no work: it has duty')
	refused_on_sheet('mixer.duty', 1, 'unit mixer has no duty: it has no variables')
	// A scope alone, a family alone, and a family's name without its scope.
	refused_on_sheet('2*product', 3, 'product is a stream: name one of its variables, as product.T')
	refused_on_sheet('mixer', 1, 'mixer is a unit, with no variables')
	refused_on_sheet('product.n', 1, 'product.n is a family of variables: say which, as product.n[CH4]')
	refused_on_sheet('n[CH4]', 1, "n is a stream's: say whose, as feed.n[CH4]")
	refused_on_sheet('duty', 1, "duty is a unit's: say whose, as reactor.duty")
	refused_on_sheet('product.n[CO]', 11, 'no CO in product.n[ ]')
	refused_on_sheet("'let down'.y[]", 14, "'let down'.y[ ] needs a name inside")
	refused_on_sheet('product.T[CH4]', 1, 'stream product has no T')
	// A quoted name needs its variable, and its closing quote.
	refused_on_sheet("'let down'", 11, "a name in quotes is a stream's or a unit's")
	refused_on_sheet("'let down.T", 1, 'a quote with no quote to close it')
	refused_on_sheet("''.T", 1, 'empty quotes')
	// No space at the dot: product .T is product, then something else.
	refused_on_sheet('product .T', 1, 'product is a stream')
	// A number beside a dotted name.
	refused_on_sheet('2 product.T', 3, 'write 2*product.T to multiply')
	refused_on_sheet('product.T(2)', 10, "product.T isn't a function")
	// Where nothing has a dot, a dotted name says so.
	expr.parse(expr.Context{}, &none_defined, 'a.b') or {
		assert err.msg().contains('nothing here has a dot in its name'), err.msg()
		return
	}
	assert false
}

fn test_dotted_names_render_back() {
	cases := {
		'product.T - feed.T':                                 ['product.T − feed.T',
			'\\mathrm{product}.T - \\mathrm{feed}.T']
		"'let down'.y[CH4]":                                  ["'let down'.y_CH₄",
			'\\mathrm{let\\ down}.y_{\\mathrm{CH_{4}}}']
		'(feed.n[CH3OH] - product.n[CH3OH]) / feed.n[CH3OH]': [
			'(feed.n_CH₃OH − product.n_CH₃OH)/feed.n_CH₃OH',
			'',
		]
		"product.y['n-butane']":                              ['product.y_H₂O', '']
		'reactor.duty / feed.F':                              ['reactor.duty/feed.F',
			'\\frac{\\mathrm{reactor}.duty}{\\mathrm{feed}.F}']
	}
	for text, want in cases {
		e := on_sheet(text)!
		// Plain text reads back to the same tree.
		again := on_sheet(e.plain())!
		assert again.nodes.map(it.kind) == e.nodes.map(it.kind), text
		assert again.nodes.map(it.index) == e.nodes.map(it.index), text
		assert again.plain() == e.plain(), text
		assert e.unicode() == want[0], '${text}: ${e.unicode()}'
		if want[1] != '' {
			assert e.latex() == want[1], '${text}: ${e.latex()}'
		}
	}
	// A member keyed by a name that isn't a word renders by the context's.
	assert on_sheet("product.y['n-butane']")!.plain() == 'product.y[H2O]'
	assert on_sheet("'let down'.T")!.plain() == "'let down'.T"
}

// value is checked against eval: the same bits wherever there's a value.
fn test_a_value_is_evals_bits() {
	ctx := sheet()
	texts := [
		'(feed.n[CH3OH] - product.n[CH3OH]) / feed.n[CH3OH]',
		'ln(product.y[CH4]) + exp(-feed.T/(500 K)) * sqrt(feed.P/(1 bar))',
		'reactor.duty / feed.F / (R * feed.T) + (feed.T in degC)^1.5',
		'abs(product.T - feed.T) / (2 K)',
		'(product.F / (1 mol/s)) [%]',
	]
	rand.seed([u32(8), 2])
	for text in texts {
		e := on_sheet(text)!
		for _ in 0 .. 200 {
			mut x := []f64{len: ctx.slots()}
			for k in 0 .. x.len {
				x[k] = 300.0 + 1000.0 * rand.f64() // every T above 0 °C
			}
			got := e.value(&none_defined, x, [])!
			want := e.eval(&none_defined, x)
			assert math.f64_bits(got) == math.f64_bits(want), text
		}
	}
}

fn no_value(text string, x []f64, why []string, column int, want string) {
	e := on_sheet(text) or { panic(err) }
	e.value(&none_defined, x, why) or {
		assert err is expr.ExprError, err.msg()
		f := err as expr.ExprError
		assert f.reason.contains(want), '${text}: ${f.reason}'
		assert f.line == 1 && f.column == column, '${text}: column ${f.column}, not ${column}: ${f.reason}'
		// eval itself gives what value refuses: ±∞ or NaN.
		assert !math.is_finite(e.eval(&none_defined, x)), text
		return
	}
	assert false, '${text} had a value'
}

fn test_no_value_says_why() {
	ctx := sheet()
	mut x := []f64{len: ctx.slots(), init: 2.0}
	// Tutorial 2's methanol conversion: the feed has none.
	x[ctx.slot(0, 2)] = 0.0
	x[ctx.slot(1, 2)] = 0.0
	x[ctx.slot(2, 2)] = 0.0
	no_value('(feed.n[CH3OH] - product.n[CH3OH]) / feed.n[CH3OH]', x, [], 36, 'it divides by feed.n[CH3OH], which is zero here')
	// x/0 with x not zero, which eval gives as ∞.
	no_value('product.F / feed.n[CH3OH]', x, [], 11, 'it divides by feed.n[CH3OH], which is zero here')
	// A logarithm, a root and a power.
	no_value('ln(product.y[CH3OH])', x, [], 1, 'it takes ln of product.y[CH3OH], which is 0 here: a logarithm needs a number above zero')
	no_value('log10(product.y[CH4] - 2)', x, [], 1, 'it takes log10 of product.y[CH4] - 2, which is 0 here')
	no_value('sqrt(product.y[CH4] - 3)', x, [], 1, 'it takes the square root of product.y[CH4] - 3, which is -1 here')
	no_value('product.y[CH3OH]^-1', x, [], 17, 'product.y[CH3OH] is zero here, and its power -1 is negative')
	no_value('(product.y[CH4] - 3)^0.5', x, [], 21, 'is -1 here, and its power 0.5 is 0.5, not a whole number')
	// Overflow, in exp and in a product.
	no_value('exp(product.y[CH4] * 1000)', x, [], 1, 'exp(product.y[CH4] * 1000) is too large for a number here: its argument is 2000')
	no_value('product.y[CH4] * 1e308 * 1e308', x, [], 16, 'product.y[CH4] * 1e308 is too large for a number here')
	// An input with no value, with its reason where the caller has one.
	mut y := x.clone()
	y[ctx.slot(2, 0)] = math.nan()
	mut why := []string{len: ctx.slots()}
	why[ctx.slot(2, 0)] = 'product has no flow'
	no_value('1 + product.y[CH4]', y, why, 5, 'product.y[CH4] has no value here: product has no flow')
	no_value('1 + product.y[CH4]', y, [], 5, 'product.y[CH4] has no value here')
	// A definition with no value says so, and why.
	defs := expr.define(ctx, [expr.Source{'conversion', '(feed.n[CH3OH] - product.n[CH3OH]) / feed.n[CH3OH]'}])!
	c := expr.parse(ctx, &defs, '100 * conversion')!
	c.value(&defs, x, []) or {
		assert err.msg().contains('column 7: conversion has no value here: it divides by feed.n[CH3OH]'), err.msg()
		return
	}
	assert false
}

fn equation_of(text string) !expr.Equation {
	return expr.equation(sheet(), &none_defined, text)
}

fn refused_equation(text string, column int, want string) {
	equation_of(text) or {
		assert err is expr.ExprError, err.msg()
		e := err as expr.ExprError
		assert e.reason.contains(want), '${text}: ${e.reason}'
		assert e.text == text, text
		assert e.line == 1 && e.column == column, '${text}: column ${e.column}, not ${column}: ${e.reason}'
		return
	}
	assert false, '"${text}" was read'
}

fn test_equations() {
	q := equation_of('product.F = 60 mol/s')!
	assert q.relation == .eq && q.column == 11
	assert q.lhs.nodes[q.lhs.root].kind == .var && q.rhs.nodes[q.rhs.root].kind == .num
	assert q.rhs.eval(&none_defined, []) == 60.0
	assert q.lhs.text == q.text && q.rhs.text == q.text
	assert q.plain() == 'product.F = 60 mol/s'
	for text, rel in {
		'product.T <= 1100 K':     expr.Relation.le
		'product.T >= 300 K':      .ge
		'product.T ≤ 1100 K':      .le
		'1100 K ≥ product.T':      .ge
		"'let down'.y[CH4] = 5 %": .eq
	} {
		e := equation_of(text)!
		assert e.relation == rel, text
		again := equation_of(e.plain())!
		assert again.relation == rel && again.plain() == e.plain(), text
	}
	assert equation_of('product.T <= 1100 K')!.unicode() == 'product.T ≤ 1100 K'
	assert equation_of('product.T >= 300 K')!.latex() == '\\mathrm{product}.T \\ge 300\\,\\mathrm{K}'
	// A quoted key with = in it is a key, not a relation.
	k := expr.Context{
		families: [
			expr.Family{'y', quantity.Dim{}, ['a=b'], {
				'a=b': 0
			}},
		]
	}
	assert expr.equation(k, &none_defined, "y['a=b'] = 0.5")!.column == 10
	assert expr.equation(k, &none_defined, 'y["a=b"] = 0.5')!.column == 10
}

fn test_equations_refused() {
	refused_equation('product.F', 1, 'an equation needs =, <= or >= between its two sides')
	refused_equation('product.F = 60 mol/s = 2', 22, 'a second relation')
	refused_equation('product.T < 300 K', 11, 'a strict <: write <=')
	refused_equation('product.T > 300 K', 11, 'a strict >: write >=')
	refused_equation('product.T == 300 K', 11, 'an == here: write = once')
	refused_equation('product.T =< 300 K', 11, 'write <= with the < first')
	refused_equation('exp(product.T = 3)', 15, 'a relation inside brackets')
	refused_equation(' = 300 K', 2, 'the left side of the equation is empty')
	refused_equation('product.T = ', 11, 'the right side of the equation is empty')
	// Each side's own refusal, at its place in the whole text.
	refused_equation('product.F = 60 mol/s + prodct.F', 24, 'no stream or unit prodct here')
	refused_equation('product.F = 60 mol/ss', 20, 'in the unit mol/ss')
	// The two sides must share a dimension.
	refused_equation('product.F = 5 %', 11, 'the two sides of = must have the same units, but the left is a molar flow')
	refused_equation('product.T <= 1100', 11, 'to use a temperature')
}
