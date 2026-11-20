module numeric

import math

// Jet is a number with its first two derivatives along one direction:
// v, dv/ds and d²v/ds². Arithmetic on jets is forward-mode differentiation,
// exact to rounding, so a formula written once, generic over f64 and Jet,
// gives its value on f64 and its derivatives on Jet (M6's activity models,
// whose T-derivatives and composition derivatives come this way, and
// expr's jet pass). The value is computed as the f64 arithmetic computes
// it, so a jet's v is the plain formula's bits.
pub struct Jet {
pub:
	v  f64
	d  f64
	dd f64
}

// constant is c, which doesn't move: its derivatives are zero.
@[inline]
pub fn constant(c f64) Jet {
	return Jet{
		v: c
	}
}

// along is x, the direction itself: dx/ds = 1.
@[inline]
pub fn along(x f64) Jet {
	return Jet{
		v: x
		d: 1.0
	}
}

pub fn (a Jet) + (b Jet) Jet {
	return Jet{a.v + b.v, a.d + b.d, a.dd + b.dd}
}

pub fn (a Jet) - (b Jet) Jet {
	return Jet{a.v - b.v, a.d - b.d, a.dd - b.dd}
}

pub fn (a Jet) * (b Jet) Jet {
	return Jet{a.v * b.v, a.d * b.v + a.v * b.d, a.dd * b.v + 2.0 * a.d * b.d + a.v * b.dd}
}

// a/b is q = a/b with q' = (a' − q·b')/b and q'' = (a'' − 2q'·b' − q·b'')/b,
// from a = q·b differentiated twice.
pub fn (a Jet) / (b Jet) Jet {
	q := a.v / b.v
	d := (a.d - q * b.d) / b.v
	return Jet{q, d, (a.dd - 2.0 * d * b.d - q * b.dd) / b.v}
}

// neg is −a.
pub fn (a Jet) neg() Jet {
	return Jet{-a.v, -a.d, -a.dd}
}

// chain is f(a) given f, f' and f'' at a.v: (f∘a)' = f'·a' and
// (f∘a)'' = f''·a'² + f'·a''.
@[inline]
fn (a Jet) chain(f f64, f1 f64, f2 f64) Jet {
	return Jet{f, f1 * a.d, f2 * a.d * a.d + f1 * a.dd}
}

// exp is eᵃ.
pub fn (a Jet) exp() Jet {
	e := math.exp(a.v)
	return a.chain(e, e, e)
}

// log is ln a.
pub fn (a Jet) log() Jet {
	return a.chain(math.log(a.v), 1.0 / a.v, -1.0 / (a.v * a.v))
}

// log10 is log₁₀ a, as math.log10 gives its value.
pub fn (a Jet) log10() Jet {
	k := 1.0 / math.ln10
	return a.chain(math.log10(a.v), k / a.v, -k / (a.v * a.v))
}

// sqrt is √a.
pub fn (a Jet) sqrt() Jet {
	s := math.sqrt(a.v)
	return a.chain(s, 0.5 / s, -0.25 / (s * a.v))
}

// abs is |a|, differentiated as ±a away from 0.
pub fn (a Jet) abs() Jet {
	return if a.v < 0.0 { a.neg() } else { a }
}

// pow is aᵇ for a b that doesn't move, as math.pow gives its value:
// b·aᵇ⁻¹ and b(b − 1)·aᵇ⁻² its derivatives.
pub fn (a Jet) pow(b f64) Jet {
	p := math.pow(a.v, b)
	p1 := b * math.pow(a.v, b - 1.0)
	return a.chain(p, p1, (b - 1.0) * p1 / a.v)
}

// pow_jet is aᵇ where both move, for a > 0: e^(b·ln a), its value as
// math.pow gives it.
pub fn (a Jet) pow_jet(b Jet) Jet {
	e := b * a.log()
	return Jet{
		...e.exp()
		v: math.pow(a.v, b.v)
	}
}
