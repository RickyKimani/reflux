// numeric.Jet (M6 step 2): each operation's value is the f64 operation's,
// to the bit, and its first two derivatives are its closed forms', to
// rounding.
module numeric

import math

fn close(a f64, b f64) bool {
	return math.abs(a - b) <= 4e-16 * math.max(math.abs(b), 1.0)
}

fn check(j Jet, v f64, d f64, dd f64) {
	assert math.f64_bits(j.v) == math.f64_bits(v), '${j.v} against ${v}'
	assert close(j.d, d), 'first: ${j.d} against ${d}'
	assert close(j.dd, dd), 'second: ${j.dd} against ${dd}'
}

fn test_each_operation_against_its_closed_form() {
	x := 1.7
	j := along(x)
	check(j.exp(), math.exp(x), math.exp(x), math.exp(x))
	check(j.log(), math.log(x), 1.0 / x, -1.0 / (x * x))
	check(j.log10(), math.log10(x), 1.0 / (x * math.ln10), -1.0 / (x * x * math.ln10))
	check(j.sqrt(), math.sqrt(x), 0.5 / math.sqrt(x), -0.25 / (x * math.sqrt(x)))
	check(j.pow(2.5), math.pow(x, 2.5), 2.5 * math.pow(x, 1.5), 3.75 * math.pow(x, 0.5))
	check(j * j, x * x, 2.0 * x, 2.0)
	check(constant(1.0) / j, 1.0 / x, -1.0 / (x * x), 2.0 / (x * x * x))
	check(j + constant(3.0), x + 3.0, 1.0, 0.0)
	check(constant(3.0) - j, 3.0 - x, -1.0, 0.0)
	check(j.neg(), -x, -1.0, 0.0)
	check(along(-x).abs(), x, -1.0, 0.0)
	// xˣ = e^(x ln x): its first derivative xˣ(ln x + 1), its second
	// xˣ((ln x + 1)² + 1/x).
	xx := math.pow(x, x)
	l := math.log(x) + 1.0
	check(j.pow_jet(j), xx, xx * l, xx * (l * l + 1.0 / x))
	// A constant doesn't move.
	check(constant(2.0) * constant(5.0), 10.0, 0.0, 0.0)
}

fn test_a_composition_by_the_chain_rule() {
	// f(x) = exp(−a/x), a τ-like factor: f′ = f·a/x², f″ = f·(a²/x⁴ − 2a/x³).
	a, x := 900.0, 350.0
	f := (constant(-a) / along(x)).exp()
	e := math.exp(-a / x)
	assert math.f64_bits(f.v) == math.f64_bits(e)
	assert math.abs(f.d - e * a / (x * x)) <= 1e-15 * math.abs(e * a / (x * x))
	want := e * (a * a / (x * x * x * x) - 2.0 * a / (x * x * x))
	assert math.abs(f.dd - want) <= 1e-14 * math.abs(want), '${f.dd} against ${want}'
}
