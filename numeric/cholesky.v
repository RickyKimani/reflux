module numeric

import math

// Cholesky is the factorisation A = L·Lᵀ of a symmetric positive definite
// matrix, with L lower triangular. It exists only when A is positive
// definite, which makes cholesky the cheap test for that, about half the
// work of an LU.
pub struct Cholesky {
pub:
	n int
mut:
	l []f64 // L, row-major, zero above the diagonal
}

// cholesky factorises the n×n row-major symmetric matrix a, reading only
// its lower triangle. It fails when a isn't positive definite to working
// precision: some pivot is no larger than n·ε times the largest diagonal
// entry.
pub fn cholesky(a []f64, n int) !Cholesky {
	if n <= 0 || a.len != n * n {
		return error('cholesky: expected ${n}×${n} = ${n * n} entries, got ${a.len}')
	}
	mut diag := 0.0
	for i in 0 .. n {
		diag = math.max(diag, math.abs(a[i * n + i]))
	}
	tol := f64(n) * f64_epsilon * diag
	mut l := []f64{len: n * n}
	for j in 0 .. n {
		mut d := a[j * n + j]
		for k in 0 .. j {
			d -= l[j * n + k] * l[j * n + k]
		}
		if !(d > tol) {
			return error('cholesky: not positive definite (pivot ${d} in column ${j})')
		}
		ljj := math.sqrt(d)
		l[j * n + j] = ljj
		for i in j + 1 .. n {
			mut s := a[i * n + j]
			for k in 0 .. j {
				s -= l[i * n + k] * l[j * n + k]
			}
			l[i * n + j] = s / ljj
		}
	}
	return Cholesky{
		n: n
		l: l
	}
}

// solve returns x with A·x = b, by L·y = b and then Lᵀ·x = y.
pub fn (c Cholesky) solve(b []f64) ![]f64 {
	n := c.n
	if b.len != n {
		return error('cholesky solve: expected a vector of length ${n}, got ${b.len}')
	}
	mut x := []f64{len: n}
	for i in 0 .. n {
		mut s := b[i]
		for k in 0 .. i {
			s -= c.l[i * n + k] * x[k]
		}
		x[i] = s / c.l[i * n + i]
	}
	for i := n - 1; i >= 0; i-- {
		mut s := x[i]
		for k in i + 1 .. n {
			s -= c.l[k * n + i] * x[k]
		}
		x[i] = s / c.l[i * n + i]
	}
	return x
}
