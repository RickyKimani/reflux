module numeric

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
