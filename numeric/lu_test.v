module numeric

import math

// A small deterministic generator (xorshift64*), so that the "random"
// matrices are the same on every machine and every V version.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) next() f64 {
	r.state ^= r.state >> 12
	r.state ^= r.state << 25
	r.state ^= r.state >> 27
	return f64((r.state * u64(2685821657736338717)) >> 11) / 9007199254740992.0 * 2.0 - 1.0
}

fn matvec(a []f64, x []f64, n int) []f64 {
	mut y := []f64{len: n}
	for i in 0 .. n {
		for j in 0 .. n {
			y[i] += a[i * n + j] * x[j]
		}
	}
	return y
}

fn norm_inf(x []f64) f64 {
	mut m := 0.0
	for v in x {
		m = math.max(m, math.abs(v))
	}
	return m
}

// mat_norm_inf is the largest absolute row sum.
fn mat_norm_inf(a []f64, n int) f64 {
	mut m := 0.0
	for i in 0 .. n {
		mut s := 0.0
		for j in 0 .. n {
			s += math.abs(a[i * n + j])
		}
		m = math.max(m, s)
	}
	return m
}

// residual_ok checks the backward-stability bound the plan asks for:
// ‖A·x − b‖ ≤ 1e-12·‖A‖·‖x‖.
fn residual_ok(a []f64, x []f64, b []f64, n int) bool {
	ax := matvec(a, x, n)
	mut r := []f64{len: n}
	for i in 0 .. n {
		r[i] = ax[i] - b[i]
	}
	return norm_inf(r) <= 1e-12 * mat_norm_inf(a, n) * norm_inf(x)
}

fn test_random_systems_solve_accurately() {
	mut rng := Rng{
		state: 0x9E3779B97F4A7C15
	}
	for n in [1, 2, 3, 5, 10, 50] {
		for _ in 0 .. 5 {
			// Random entries plus n on the diagonal: well conditioned, but
			// pivoting still reorders rows.
			mut a := []f64{len: n * n, init: rng.next()}
			for i in 0 .. n {
				a[i * n + i] += f64(n) * rng.next()
			}
			x_true := []f64{len: n, init: rng.next()}
			b := matvec(a, x_true, n)
			x := lu_factor(a, n)!.solve(b)!
			assert residual_ok(a, x, b, n), 'n = ${n}: residual too large'
			for i in 0 .. n {
				assert math.abs(x[i] - x_true[i]) <= 1e-9, 'n = ${n}: x[${i}] = ${x[i]}, expected ${x_true[i]}'
			}
		}
	}
}

fn test_pivoting_handles_zero_and_tiny_leading_entries() {
	// A zero in the (0, 0) position: without a row swap the first pivot is 0.
	x1 := lu_factor([0.0, 1.0, 1.0, 0.0], 2)!.solve([3.0, 7.0])!
	assert x1 == [7.0, 3.0]
	// A tiny (0, 0) entry: without pivoting, 1 − 1e20 swamps the other row
	// and the answer is badly wrong. With pivoting it's right to rounding.
	x2 := lu_factor([1e-20, 1.0, 1.0, 1.0], 2)!.solve([1.0, 2.0])!
	assert math.tolerance(x2[0], 1.0, 1e-15) && math.tolerance(x2[1], 1.0, 1e-15), '${x2}'
}

fn test_singular_matrices_are_errors_not_nans() {
	for name, m in {
		'rank 1 2×2': [1.0, 2.0, 2.0, 4.0]
		'zero 2×2':   [0.0, 0.0, 0.0, 0.0]
		'rank 2 3×3': [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
	} {
		n := int(math.sqrt(f64(m.len)))
		if f := lu_factor(m, n) {
			assert false, '${name}: expected an error, got det ${f.det()}'
		}
	}
}

fn hilbert(n int) []f64 {
	return []f64{len: n * n, init: 1.0 / f64(index / n + index % n + 1)}
}

fn test_an_ill_conditioned_system_still_has_a_small_residual() {
	// Hilbert(8) has an ∞-norm condition number around 3e10. LU with partial
	// pivoting is backward stable, so the residual stays tiny. The forward
	// error, though, can grow to about cond·ε ≈ 1e-5, and that's the
	// matrix's fault, not the solver's.
	n := 8
	a := hilbert(n)
	x_true := []f64{len: n, init: 1.0}
	b := matvec(a, x_true, n)
	x := lu_factor(a, n)!.solve(b)!
	assert residual_ok(a, x, b, n)
	mut err := 0.0
	for i in 0 .. n {
		err = math.max(err, math.abs(x[i] - 1.0))
	}
	assert err < 1e-4, 'forward error ${err}'
}

fn test_determinant() {
	assert lu_factor([2.0, 0.0, 0.0, 3.0], 2)!.det() == 6.0
	assert lu_factor([0.0, 1.0, 1.0, 0.0], 2)!.det() == -1.0 // one row swap
	// det(Hilbert(3)) = 1/2160.
	assert math.tolerance(lu_factor(hilbert(3), 3)!.det(), 1.0 / 2160.0, 1e-13)
}

fn test_solve_into_matches_solve_and_reuses_the_factorisation() {
	a := [4.0, 1.0, 2.0, 1.0, 5.0, 3.0, 2.0, 3.0, 6.0]
	f := lu_factor(a, 3)!
	mut x := []f64{len: 3}
	for b in [[1.0, 2.0, 3.0], [-1.0, 0.5, 4.0]] {
		f.solve_into(b, mut x)!
		assert x == f.solve(b)!
		assert residual_ok(a, x, b, 3)
	}
}

fn test_input_is_left_alone_and_sizes_are_checked() {
	a := [0.0, 1.0, 1.0, 0.0]
	_ := lu_factor(a, 2)!
	assert a == [0.0, 1.0, 1.0, 0.0], 'lu_factor modified its input'
	if _ := lu_factor([1.0, 2.0, 3.0], 2) {
		assert false, 'wrong matrix size accepted'
	}
	f := lu_factor(a, 2)!
	if _ := f.solve([1.0]) {
		assert false, 'wrong vector size accepted'
	}
}
