// Module numeric holds plain numerics with no chemistry in it.
module numeric

import math

// cubic_real_roots returns the real roots of the monic cubic
//
//     x³ + c2·x² + c1·x + c0 = 0
//
// in ascending order. It returns one root when the other two are complex
// and three when all are real. Near a double or triple root, rounding can
// tip the discriminant either way, so callers must not rely on the count
// there.
pub fn cubic_real_roots(c2 f64, c1 f64, c0 f64) []f64 {
	// Substitute x = t − c2/3 to get the depressed cubic t³ + p·t + q = 0.
	shift := c2 / 3.0
	p := c1 - c2 * shift
	q := 2.0 * shift * shift * shift - shift * c1 + c0

	half_q := q / 2.0
	third_p := p / 3.0
	disc := half_q * half_q + third_p * third_p * third_p

	mut roots := []f64{cap: 3}
	if disc > 0.0 {
		// One real root. Take the cube root of the larger-magnitude term
		// so the sum u + v does not cancel.
		mag := math.cbrt(math.abs(half_q) + math.sqrt(disc))
		u := if q >= 0.0 { -mag } else { mag }
		t := if u == 0.0 { 0.0 } else { u - third_p / u }
		roots << t - shift
	} else if p == 0.0 {
		// disc <= 0 and p == 0 force q == 0: a triple root.
		roots << -shift
	} else {
		// Three real roots: the trigonometric (Viète) form.
		m := 2.0 * math.sqrt(-third_p)
		arg := math.clamp(3.0 * q / (2.0 * p) * math.sqrt(-3.0 / p), -1.0, 1.0)
		theta := math.acos(arg) / 3.0
		for k in 0 .. 3 {
			roots << m * math.cos(theta - 2.0 * math.pi * f64(k) / 3.0) - shift
		}
	}

	for mut r in roots {
		r = polish(r, c2, c1, c0)
	}
	roots.sort()
	return roots
}

// polish takes a few Newton steps on the original (undepressed) cubic.
// The closed forms lose digits when the roots differ by orders of
// magnitude, as liquid and vapor Z do at low pressure; Newton gets them
// back. A step is only kept if it reduces the residual, so a root sitting
// on a near-zero derivative cannot be thrown somewhere worse.
fn polish(x0 f64, c2 f64, c1 f64, c0 f64) f64 {
	mut x := x0
	mut f := eval(x, c2, c1, c0)
	for _ in 0 .. 3 {
		df := (3.0 * x + 2.0 * c2) * x + c1
		if f == 0.0 || df == 0.0 {
			break
		}
		next := x - f / df
		f_next := eval(next, c2, c1, c0)
		if !math.is_finite(next) || math.abs(f_next) >= math.abs(f) {
			break
		}
		x, f = next, f_next
	}
	return x
}

@[inline]
fn eval(x f64, c2 f64, c1 f64, c0 f64) f64 {
	return ((x + c2) * x + c1) * x + c0
}
