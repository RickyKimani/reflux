module numeric

import math

fn test_cholesky_solves_a_spd_system() {
	// A = MᵀM + I for a fixed M is symmetric positive definite.
	n := 5
	mut m := []f64{len: n * n, init: math.sin(f64(index + 1))}
	mut a := []f64{len: n * n}
	for i in 0 .. n {
		for j in 0 .. n {
			mut s := if i == j { 1.0 } else { 0.0 }
			for k in 0 .. n {
				s += m[k * n + i] * m[k * n + j]
			}
			a[i * n + j] = s
		}
	}
	b := []f64{len: n, init: f64(index) - 2.0}
	x := cholesky(a, n)!.solve(b)!
	for i in 0 .. n {
		mut r := -b[i]
		for j in 0 .. n {
			r += a[i * n + j] * x[j]
		}
		assert math.abs(r) <= 1e-13
	}
	// It agrees with LU.
	y := lu_factor(a, n)!.solve(b)!
	for i in 0 .. n {
		assert math.abs(x[i] - y[i]) <= 1e-13 * math.max(1.0, math.abs(y[i]))
	}
}

fn test_cholesky_refuses_what_is_not_positive_definite() {
	for a in [[1.0, 2.0, 2.0, 1.0], [1.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 1.0]] {
		if _ := cholesky(a, 2) {
			assert false, '${a} is not positive definite'
		}
	}
	if _ := cholesky([1.0, 2.0, 3.0], 2) {
		assert false, 'wrong size'
	}
	// [[2, 1], [1, 2]] is: L = [[√2, 0], [1/√2, √1.5]].
	c := cholesky([2.0, 1.0, 1.0, 2.0], 2)!
	assert math.abs(c.l[0] - math.sqrt2) <= 1e-15 && c.l[1] == 0.0
	assert math.abs(c.l[2] - 1.0 / math.sqrt2) <= 1e-15 && math.abs(c.l[3] - math.sqrt(1.5)) <= 1e-15
}
