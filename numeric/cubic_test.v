// Declaring `module numeric` makes this an internal test, like a Go test
// file in `package numeric`: it can call private functions such as eval.
module numeric

import math

fn residuals_are_small(roots []f64, c2 f64, c1 f64, c0 f64) bool {
	for r in roots {
		if math.abs(eval(r, c2, c1, c0)) > 1e-12 * (1.0 + math.abs(r * r * r)) {
			return false
		}
	}
	return true
}

fn test_three_distinct_real_roots() {
	// (x − 1)(x − 2)(x − 3) = x³ − 6x² + 11x − 6
	roots := cubic_real_roots(-6.0, 11.0, -6.0)
	assert roots.len == 3
	assert math.tolerance(roots[0], 1.0, 1e-14)
	assert math.tolerance(roots[1], 2.0, 1e-14)
	assert math.tolerance(roots[2], 3.0, 1e-14)
}

fn test_one_real_root() {
	// (x − 2)(x² + 1) = x³ − 2x² + x − 2
	roots := cubic_real_roots(-2.0, 1.0, -2.0)
	assert roots.len == 1
	assert math.tolerance(roots[0], 2.0, 1e-14)
}

fn test_exact_triple_root() {
	// (x − 3/8)³: van der Waals at its own critical point.
	z := 3.0 / 8.0
	roots := cubic_real_roots(-3.0 * z, 3.0 * z * z, -z * z * z)
	assert roots.len >= 1
	for r in roots {
		assert math.tolerance(r, z, 1e-6), 'root ${r} is not ${z}'
	}
}

fn test_widely_separated_roots() {
	// Liquid-like and vapor-like Z at low pressure differ by ~10⁴.
	// (x − 1e-4)(x − 0.5)(x − 1)
	a, b, c := 1e-4, 0.5, 1.0
	c2 := -(a + b + c)
	c1 := a * b + a * c + b * c
	c0 := -a * b * c
	roots := cubic_real_roots(c2, c1, c0)
	assert roots.len == 3
	assert math.tolerance(roots[0], a, 1e-15)
	assert residuals_are_small(roots, c2, c1, c0)
}

fn test_roots_come_back_sorted() {
	for c0 in [-6.0, -1.0, 0.0, 0.5, 2.0] {
		roots := cubic_real_roots(-1.5, -2.0, c0)
		for i in 1 .. roots.len {
			assert roots[i - 1] <= roots[i]
		}
		assert residuals_are_small(roots, -1.5, -2.0, c0)
	}
}
