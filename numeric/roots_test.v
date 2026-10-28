module numeric

import math

fn test_bisect_converges_to_machine_precision() {
	r := bisect(fn (x f64) f64 {
		return x * x - 2.0
	}, 0.0, 2.0)!
	assert math.abs(r - math.sqrt2) <= 2.5e-16
}

fn test_bisect_with_a_capturing_closure() {
	// The closure captures `target`; V requires the capture list [target].
	target := 0.3
	r := bisect(fn [target] (x f64) f64 {
		return math.cos(x) - target
	}, 0.0, 3.0)!
	assert math.tolerance(r, math.acos(target), 1e-15)
}

fn test_bisect_accepts_a_root_on_the_bracket() {
	r := bisect(fn (x f64) f64 {
		return x - 1.0
	}, 1.0, 5.0)!
	assert r == 1.0
}

fn test_bisect_rejects_a_bracket_without_a_sign_change() {
	if _ := bisect(fn (x f64) f64 {
		return x * x + 1.0
	}, -1.0, 1.0)
	{
		assert false, 'expected an error'
	}
}

fn test_brent_finds_smooth_roots_fast() {
	// Wallis's cubic x³ − 2x − 5 and cos x = x, to round-off in a few
	// evaluations (measured: 8 and 9; bisection needs about 50).
	for c in [[2.0, 3.0, 2.0945514815423265], [0.0, 1.0, 0.7390851332151607]] {
		r := brent(fn [c] (x f64) !f64 {
			return if c[2] > 1.0 { x * x * x - 2.0 * x - 5.0 } else { math.cos(x) - x }
		}, c[0], c[1])!
		assert math.abs(r.x - c[2]) <= 4e-16 * c[2], '${r.x} vs ${c[2]}'
		assert r.evaluations <= 10, '${r.evaluations} evaluations'
	}
}

fn test_brent_closes_in_on_a_jump() {
	// No root, only a sign change at 0.3: the bracket shrinks onto it.
	r := brent(fn (x f64) !f64 {
		return if x < 0.3 { -1.0 } else { 2.0 }
	}, 0.0, 1.0)!
	assert r.a < 0.3 && r.b >= 0.3 && r.b - r.a <= 1e-15, '${r.a} ${r.b}'
}

fn test_brent_stops_at_ftol_and_xtol() {
	// x·e⁻ˣ = 1e-3 has its root near 0.001001: with ftol the answer meets
	// it, and with xtol the last bracket is no wider than xtol and holds the
	// root (measured 5e-7 wide at xtol 1e-6).
	f := fn (x f64) !f64 {
		return x * math.exp(-x) - 1e-3
	}
	root := brent(f, -1.0, 4.0)!.x
	r := brent(f, -1.0, 4.0, ftol: 1e-6)!
	assert math.abs(r.f) <= 1e-6
	s := brent(f, -1.0, 4.0, xtol: 1e-6)!
	assert s.b - s.a <= 1e-6 && s.a <= root && root <= s.b, '${s.a} ${s.b}'
}

fn test_brent_is_never_much_slower_than_bisection() {
	// e^(−1/x) − 1 is flat near 0 and jumps to 1 at 0.7: interpolation alone
	// creeps there, and Brent's guard falls back to bisection. Halving the
	// bracket from 4 to round-off takes about 52 steps; Brent takes 54.
	r := brent(fn (x f64) !f64 {
		return if x < 0.7 { math.exp(-1.0 / x) - 1.0 } else { 1.0 }
	}, 0.01, 4.0)!
	assert r.b - r.a <= 1e-15 && r.evaluations <= 60, '${r.evaluations} evaluations'
}

fn test_brent_refuses_and_passes_errors_on() {
	if _ := brent(fn (x f64) !f64 {
		return x * x + 1.0
	}, -1.0, 1.0)
	{
		assert false, 'no sign change, yet a root'
	}
	if _ := brent(fn (x f64) !f64 {
		if x > 0.5 {
			return error('out of range at ${x}')
		}
		return x - 0.9
	}, 0.0, 0.4)
	{
		assert false, 'no sign change, yet a root'
	}
	brent(fn (x f64) !f64 {
		if x > 0.5 && x < 1.0 {
			return error('out of range')
		}
		return x - 0.7
	}, 0.0, 1.0) or {
		assert err.msg().contains('out of range'), err.msg()
		return
	}
	assert false, 'the error from f was lost'
}
