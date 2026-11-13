// Readings (D2, as amended before M9 step 10): a value read on a unit's
// scale, (T in degC), and a group with a unit after it, (…) [%], so that a
// correlation goes in as its source prints it, the offset applied by reflux
// and never typed.
import expr
import math
import quantity

fn context() expr.Context {
	species := ['CH4', 'H2O']
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, species, {
				'CH4': 0
				'H2O': 1
			}},
		]
	}
}

const no_defs = expr.Definitions{}

fn read(text string) !expr.Expr {
	return expr.parse(context(), &no_defs, text)
}

// value is an expression's value at T and P, both partial pressures 1 bar.
fn value(text string, t f64, p f64) f64 {
	e := read(text) or { panic('${text}: ${err.msg()}') }
	return e.eval(&no_defs, [t, p, 1e5, 1e5])
}

fn refusal(text string) string {
	read(text) or { return err.msg() }
	return 'accepted'
}

fn test_a_reading_is_the_number_a_scale_shows() {
	// Points where the scales are known exactly, independently of how
	// quantity converts: 0 °C is 273.15 K, −40 °C is −40 °F, 1 barg is a
	// bar above a standard atmosphere, which is 760 torr, and a conventional
	// millimetre of mercury is 133.322387415 Pa (NIST SP 811). Each within
	// 1e-12 absolute, the rounding of numbers near 500.
	cases := [
		['T in degC', '273.15', '0'],
		['T in degC', '593.15', '320'],
		['T in °C', '233.15', '-40'],
		['T in degF', '233.15', '-40'],
		['T in degF', '273.15', '32'],
		['T in K', '593.15', '593.15'],
		['T in R', '100', '180'],
		['P in barg', '201325', '1'],
		['P in barg', '101325', '0'],
		['P in torr', '101325', '760'],
		['P in mmHg', '133.322387415', '1'],
		['P in bar', '7e6', '70'],
		['P in psig', '108219.757293168', '1'],
	]
	for c in cases {
		x := c[1].f64()
		got := value(c[0], x, x)
		assert math.abs(got - c[2].f64()) <= 1e-12, '${c[0]} at ${c[1]}: ${got}, not ${c[2]}'
	}
	// Brackets say which unit, as after a number; the reading of a member
	// and of a definition needs no offset.
	assert value('T in [degC]', 273.15, 0) == 0.0
	assert value('p[CH4] in bar', 300, 0) == 1.0
}

fn test_in_binds_to_the_value_before_it() {
	// As a unit binds to its number: 0.0018*T in degC is 0.0018·(T in °C),
	// T in degC^2 has the unit degC^2, and a power applies to the reading.
	a := read('0.85 - 0.0018*T in degC')!
	b := read('0.85 - 0.0018*(T in degC)')!
	assert a.plain() == b.plain()
	assert a.plain() == '0.85 - 0.0018 * (T in degC)'
	assert value('0.85 - 0.0018*T in degC', 593.15, 0) == value('0.85 - 0.0018*(T in degC)',
		593.15, 0)
	assert math.abs(value('T in degC ^ 2', 283.15, 0) - 100.0) <= 1e-12
	assert refusal('T in degC^2').contains("can't be read in degC^2")
	// A name that starts with in is a name: Tin is one, and "T inlet" two
	// values side by side.
	assert refusal('T inlet').contains('no name inlet') || refusal('T inlet').contains('side by side')
	// After a number, in alone is still the inch.
	assert math.abs(value('P * 2 in', 0, 1.0) - 0.0508) <= 1e-17
}

fn test_the_tutorials_correlation_reads_as_printed() {
	// The owner's Tutorial 2: X = 0.85 − 0.0018 T, T in °C, and DWSIM's form,
	// X% = 85 − 0.18 (T − 273.15) with T in K, both 27.4% at 320 °C.
	assert math.abs(value('0.85 - 0.0018*(T in degC)', 593.15, 0) - 0.274) <= 1e-15
	pct := '(85 - 0.18*(T in degC)) [%]'
	assert math.abs(value(pct, 593.15, 0) - 0.274) <= 1e-15
	assert read(pct)!.dim.is_dimensionless()
	// % on each number, as before, gives the same.
	assert math.abs(value('85 % - 0.18 % * (T in degC)', 593.15, 0) - 0.274) <= 1e-15
	// A unit after a group carries its dimension.
	g := read('(T - 300 K) [kJ/(mol*K)]')!
	assert g.dim.eq(quantity.unit('kJ/mol')!.dim)
	assert value('(T - 300 K) [kJ/(mol*K)]', 301, 0) == 1000.0
	// A definition names a reading once, for every use.
	defs := expr.define(context(), [expr.Source{'t', 'T in degC'}])!
	e := expr.parse(context(), &defs, '0.9 - 2.1e-3*t + 1.5e-6*t^2')!
	x := e.eval(&defs, [593.15, 0, 0, 0])
	t := 593.15 - 273.15
	assert math.abs(x - (0.9 - 2.1e-3 * t + 1.5e-6 * t * t)) <= 1e-15
}

fn test_what_a_reading_refuses_and_why() {
	// An offset reads only a context variable; each other case says why.
	assert refusal('(T - 300 K) in degC').contains('not an expression, which may be a difference')
	assert refusal('p[CH4] in barg').contains('not a member of a family')
	assert refusal('(300 K) in degC').contains('not a number')
	assert refusal('2 in degC').contains('a number before in')
	why := refusal('T in kPa')
	assert why.contains("can't be read in kPa"), why
	assert refusal('T in').contains('in needs a unit after it')
	assert refusal('T in * 2').contains('in needs a unit after it')
	assert refusal('T in degK').contains('in the unit degK')
	assert refusal('(T) [degC]').contains('has an offset')
	assert refusal('T + 2 degC').contains('read a variable on its scale, as (T in degC)')
	defs := expr.define(context(), [expr.Source{'Tw', '350 K'}])!
	expr.parse(context(), &defs, 'Tw in degC') or {
		assert err.msg().contains('not a definition, which may be a difference')
	}
	assert expr.parse(context(), &defs, 'Tw in K')!.eval(&defs, [0.0, 0, 0, 0]) == 350.0
	expr.define(context(), [expr.Source{'in', '2'}]) or {
		assert err.msg().contains('in is a keyword')
		return
	}
	assert false, 'a definition named in'
}

fn test_a_forgotten_reading_is_refused_with_the_fix() {
	// The source's equation typed as printed: T is a temperature, so the
	// sum is refused, and the message says how to read it.
	why := refusal('0.85 - 0.0018*T')
	assert why.contains('the two sides of −'), why
	assert why.contains('(T in degC)'), why
	assert refusal('1 + P / 1e5').contains('(P in bar)')
	// And where a pure number is asked for.
	e := read('0.0018*T')!
	e.must_be(quantity.Dim{}, 'a conversion') or {
		assert err.msg().contains('(T in degC)'), err.msg()
		return
	}
	assert false
}

fn test_readings_on_the_tape() {
	// The tape gives the tree's bits, and the derivatives are the scales':
	// 1 per kelvin for °C, 1.8 for °F, 1e-5 per pascal for barg, and the
	// tutorial's −0.0018 per kelvin.
	cases := [
		['T in degC', '1'],
		['T in degF', '1.8'],
		['P in barg', '1e-5'],
		['(85 - 0.18*(T in degC)) [%]', '-0.0018'],
		['0.85 - 0.0018*T in degC', '-0.0018'],
	]
	ctx := context()
	for c in cases {
		e := read(c[0])!
		tape := expr.compile(ctx, &no_defs, [e])
		mut buf := tape.buffer()
		x := [593.15, 201325.0, 1e5, 1e5]
		tape.eval(x, mut buf)
		assert math.f64_bits(tape.output(&buf, 0)) == math.f64_bits(e.eval(&no_defs, x)), c[0]
		mut g := []f64{len: ctx.slots()}
		tape.gradient(mut buf, mut g)
		d := g[0] + g[1]
		assert math.abs(d - c[1].f64()) <= 1e-15 * math.max(1.0, math.abs(c[1].f64())), '${c[0]}: ${d}'
	}
}

fn test_readings_render_and_read_back() {
	// Always in parentheses, °C for degC, and the plain form reads back to
	// the same tree.
	cases := [
		['0.85 - 0.0018*T in degC', '0.85 - 0.0018 * (T in degC)', '0.85 − 0.0018·(T in °C)',
			'0.85 - 0.0018 \\cdot \\left(T\\ \\mathrm{in}\\ \\mathrm{{}^{\\circ}C}\\right)'],
		['(85 - 0.18*(T in degC)) [%]', '(85 - 0.18 * (T in degC)) [%]', '(85 − 0.18·(T in °C)) %',
			'\\left(85 - 0.18 \\cdot \\left(T\\ \\mathrm{in}\\ \\mathrm{{}^{\\circ}C}\\right)\\right)\\,\\mathrm{\\%}'],
		['log10(P in mmHg)', 'log10((P in mmHg))', 'log₁₀((P in mmHg))',
			'\\log_{10}\\left(\\left(P\\ \\mathrm{in}\\ \\mathrm{mmHg}\\right)\\right)'],
		['T in [degF]', '(T in [degF])', '(T in °F)',
			'\\left(T\\ \\mathrm{in}\\ \\mathrm{{}^{\\circ}F}\\right)'],
		['P / (P - 1 bar) [K/Pa]', 'P / (P - 1 bar) [K/Pa]', 'P/((P − 1 bar) K/Pa)',
			'\\frac{P}{\\left(P - 1\\,\\mathrm{bar}\\right)\\,\\mathrm{K/Pa}}'],
	]
	for c in cases {
		e := read(c[0])!
		assert e.plain() == c[1], '${c[0]}: plain ${e.plain()}'
		assert e.unicode() == c[2], '${c[0]}: unicode ${e.unicode()}'
		assert e.latex() == c[3], '${c[0]}: latex ${e.latex()}'
		back := read(e.plain())!
		assert back.nodes.len == e.nodes.len && back.root == e.root, c[0]
		for i, n in e.nodes {
			m := back.nodes[i]
			assert n.kind == m.kind && n.a == m.a && n.b == m.b && n.unit == m.unit
				&& math.f64_bits(n.value) == math.f64_bits(m.value) && n.pre == m.pre
				&& n.post == m.post, '${c[0]}: node ${i}'
		}
	}
}
