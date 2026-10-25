module numeric

import math

// LU is the factorisation P·A = L·U of a square matrix, with partial
// pivoting: at each column the row with the largest pivot is swapped up,
// which keeps the multipliers in L at most 1 in size and the solve
// backward-stable.
//
// This is sized for the systems thermodynamics produces: nc×nc Newton steps
// in a flash, where nc is a few to a few dozen. Large or sparse systems
// (columns, equation-oriented flowsheets) will use a C backend, as planned
// in docs/roadmap.md.
pub struct LU {
pub:
	n int
mut:
	// lu holds L strictly below the diagonal (its unit diagonal is implied)
	// and U on and above it, row-major.
	lu []f64
	// perm[i] is the row of A that ended up in row i of P·A.
	perm []int
	// swaps counts row exchanges, for the sign of the determinant.
	swaps int
}

// lu_factor factorises the n×n row-major matrix a. a itself isn't changed.
// It fails if a is singular to working precision: if some pivot is no
// larger than n·ε·max|aᵢⱼ|, a solve would amplify rounding without bound.
pub fn lu_factor(a []f64, n int) !LU {
	if n <= 0 || a.len != n * n {
		return error('lu_factor: expected ${n}×${n} = ${n * n} entries, got ${a.len}')
	}
	mut lu := a.clone()
	mut perm := []int{len: n, init: index}
	mut swaps := 0
	mut scale := 0.0
	for x in a {
		scale = math.max(scale, math.abs(x))
	}
	tol := f64(n) * f64_epsilon * scale
	for k in 0 .. n {
		// Partial pivoting: the largest |entry| in column k, rows k..n−1.
		mut p := k
		for i in k + 1 .. n {
			if math.abs(lu[i * n + k]) > math.abs(lu[p * n + k]) {
				p = i
			}
		}
		pivot := lu[p * n + k]
		if math.abs(pivot) <= tol {
			return error('lu_factor: matrix is singular to working precision (pivot ${pivot} in column ${k})')
		}
		if p != k {
			for j in 0 .. n {
				lu[k * n + j], lu[p * n + j] = lu[p * n + j], lu[k * n + j]
			}
			perm[k], perm[p] = perm[p], perm[k]
			swaps++
		}
		for i in k + 1 .. n {
			l := lu[i * n + k] / pivot
			lu[i * n + k] = l
			if l != 0.0 {
				for j in k + 1 .. n {
					lu[i * n + j] -= l * lu[k * n + j]
				}
			}
		}
	}
	return LU{
		n:     n
		lu:    lu
		perm:  perm
		swaps: swaps
	}
}

// solve returns x with A·x = b.
pub fn (f LU) solve(b []f64) ![]f64 {
	mut x := []f64{len: f.n}
	f.solve_into(b, mut x)!
	return x
}

// solve_into writes the solution of A·x = b into x, without allocating:
// the form to use inside an iteration. Both slices must have length n, and
// they must not be the same slice.
pub fn (f LU) solve_into(b []f64, mut x []f64) ! {
	n := f.n
	if b.len != n || x.len != n {
		return error('solve_into: expected vectors of length ${n}, got b ${b.len} and x ${x.len}')
	}
	// Forward substitution, L·y = P·b (L has a unit diagonal).
	for i in 0 .. n {
		mut s := b[f.perm[i]]
		for j in 0 .. i {
			s -= f.lu[i * n + j] * x[j]
		}
		x[i] = s
	}
	// Back substitution, U·x = y.
	for i := n - 1; i >= 0; i-- {
		mut s := x[i]
		for j in i + 1 .. n {
			s -= f.lu[i * n + j] * x[j]
		}
		x[i] = s / f.lu[i * n + i]
	}
}

// det returns the determinant of A: the product of U's diagonal, with the
// sign of the row permutation.
pub fn (f LU) det() f64 {
	mut d := if f.swaps % 2 == 0 { 1.0 } else { -1.0 }
	for i in 0 .. f.n {
		d *= f.lu[i * f.n + i]
	}
	return d
}

// f64_epsilon is the gap between 1.0 and the next f64, 2⁻⁵².
const f64_epsilon = 2.220446049250313e-16
