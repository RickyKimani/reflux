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
