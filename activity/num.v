module activity

import math
import numeric

// The generic formulas' arithmetic, for N = f64 or numeric.Jet. Constants
// are lifted to N, so that an f64 constant times a jet is a jet.

@[inline]
fn lift[N](c f64) N {
	$if N is f64 {
		return c
	} $else {
		return numeric.constant(c)
	}
}

@[inline]
fn nlog[N](a N) N {
	$if N is f64 {
		return math.log(a)
	} $else {
		return a.log()
	}
}
