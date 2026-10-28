module numeric

import math

// bisect finds a root of f in [lo, hi], where f(lo) and f(hi) have opposite
// signs. It halves the bracket until floating point can't split it any
// further, so the result is as accurate as f can be evaluated. It's slow
// (one bit per step) but it can't fail once the bracket is valid, which
// makes it the right tool for one-off constants. Iterative solvers inside
// flash and flowsheet loops will use Newton or Brent instead.
//
// `f fn (f64) f64` is a function type, as in Go. Callers usually pass a
// closure, and V closures must list what they capture: fn [a, b] (x f64) f64.
pub fn bisect(f fn (f64) f64, lo f64, hi f64) !f64 {
	mut a, mut b := lo, hi
	mut fa := f(a)
	fb := f(b)
	if fa == 0.0 {
		return a
	}
	if fb == 0.0 {
		return b
	}
	if (fa > 0.0) == (fb > 0.0) {
		return error('bisect: f(${lo}) = ${fa} and f(${hi}) = ${fb} have the same sign')
	}
	for _ in 0 .. 200 {
		mid := 0.5 * (a + b)
		if mid <= a || mid >= b {
			break // a and b are adjacent floats
		}
		fm := f(mid)
		if fm == 0.0 {
			return mid
		}
		if (fm > 0.0) == (fa > 0.0) {
			a, fa = mid, fm
		} else {
			b = mid
		}
	}
	return 0.5 * (a + b)
}

// ScalarFn is a function of one variable that can fail, as a flash can.
pub type ScalarFn = fn (x f64) !f64

// BrentParams are brent's options.
@[params]
pub struct BrentParams {
pub:
	// Converged when the bracket is narrower than xtol + 4ε|x|, or when
	// |f| ≤ ftol. Both 0 by default: the bracket then shrinks to round-off.
	xtol     f64
	ftol     f64
	max_iter int = 100
}

// BrentResult is a root found by brent.
pub struct BrentResult {
pub:
	x f64 // the root
	f f64 // f(x)
	// a and b are the last bracket, with f(a) and f(b) of opposite signs
	// (or one of them zero); x is one of them.
	a           f64
	b           f64
	evaluations int
}

// brent finds a root of f in [a, b], where f(a) and f(b) have opposite
// signs, by Brent's (1973) method: inverse quadratic interpolation or the
// secant where they make progress, bisection where they don't, so it
// converges superlinearly on a smooth f and never much more slowly than
// bisection on any f. On a function with a jump it closes in on the jump.
// An error from f stops it and is passed on. Its tolerance is relative,
// 4ε|x|, plus xtol/2, so a root at x = 0 needs an xtol. It's slow on a
// multiple root, as Brent's method is (x³ = 0 takes more than 100 steps
// to 1e-12), which a flash's residual never has.
pub fn brent(f ScalarFn, a0 f64, b0 f64, params BrentParams) !BrentResult {
	if !(params.xtol >= 0.0) || !(params.ftol >= 0.0) || params.max_iter < 1 {
		return error('brent: xtol and ftol must be non-negative and max_iter positive')
	}
	mut a, mut b := a0, b0
	mut fa := f(a)!
	mut fb := f(b)!
	mut evaluations := 2
	if fa == 0.0 || math.abs(fa) <= params.ftol {
		return BrentResult{a, fa, a, b, evaluations}
	}
	if fb == 0.0 || math.abs(fb) <= params.ftol {
		return BrentResult{b, fb, a, b, evaluations}
	}
	if (fa > 0.0) == (fb > 0.0) {
		return error('brent: f(${a0}) = ${fa} and f(${b0}) = ${fb} have the same sign')
	}
	// c is the point that, with b, brackets the root; b is the best guess.
	mut c, mut fc := a, fa
	mut d := b - a
	mut e := d
	for _ in 0 .. params.max_iter {
		if (fb > 0.0) == (fc > 0.0) {
			c, fc = a, fa
			d = b - a
			e = d
		}
		if math.abs(fc) < math.abs(fb) {
			a, b, c = b, c, b
			fa, fb, fc = fb, fc, fb
		}
		tol := 4.0 * f64_epsilon * math.abs(b) + 0.5 * params.xtol
		m := 0.5 * (c - b)
		if math.abs(m) <= tol || fb == 0.0 || math.abs(fb) <= params.ftol {
			return BrentResult{b, fb, math.min(b, c), math.max(b, c), evaluations}
		}
		if math.abs(e) >= tol && math.abs(fa) > math.abs(fb) {
			s := fb / fa
			mut p, mut q := 0.0, 0.0
			if a == c {
				p = 2.0 * m * s // the secant
				q = 1.0 - s
			} else {
				// Inverse quadratic interpolation through a, b and c.
				qa := fa / fc
				r := fb / fc
				p = s * (2.0 * m * qa * (qa - r) - (b - a) * (r - 1.0))
				q = (qa - 1.0) * (r - 1.0) * (s - 1.0)
			}
			if p > 0.0 {
				q = -q
			} else {
				p = -p
			}
			if 2.0 * p < math.min(3.0 * m * q - math.abs(tol * q), math.abs(e * q)) {
				e = d
				d = p / q
			} else {
				d = m // interpolation would leave the bracket or crawl: bisect
				e = d
			}
		} else {
			d = m
			e = d
		}
		a, fa = b, fb
		b += if math.abs(d) > tol {
			d
		} else if m > 0.0 {
			tol
		} else {
			-tol
		}
		fb = f(b)!
		evaluations++
	}
	return error('brent: no convergence in ${params.max_iter} steps (bracket ${math.min(b, c)} to ${math.max(b, c)})')
}
