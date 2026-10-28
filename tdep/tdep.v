// Module tdep evaluates temperature-dependent property correlations: the
// numbered equations that ChemSep and DIPPR use for vapour pressure, heat
// capacities, densities, latent heats and transport properties.
//
// It is pure mathematics in the correlation's own units. Converting the
// result to SI per mol is the caller's job (the components module), so the
// coefficients here are always exactly what the data source gives.
module tdep

import math
import numeric

// Equation numbers this module can evaluate. ChemSep numbers its own forms
// 1–45, and the 100-series are the DIPPR equations.
const supported = [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 45, 100, 101, 102, 104, 105,
	106, 107, 114, 115, 116]!

// Forms written in reduced temperature, which need a critical temperature.
const needs_tc = [106, 114, 116]!

// How many coefficients each form actually uses. Coefficients beyond this
// are ignored, as in ChemSep.
const used_coeffs = {
	1:   1
	2:   2
	3:   3
	4:   4
	5:   5
	10:  3
	11:  1
	12:  2
	13:  3
	14:  4
	15:  5
	16:  5
	17:  5
	45:  5
	100: 5
	101: 5
	102: 4
	104: 5
	105: 4
	106: 5
	107: 5
	114: 4
	115: 5
	116: 5
}

const integral_rtol = 1e-13

// Correlation is one temperature-dependent correlation, Y(T). T is in K;
// Y is in whatever units the source used.
pub struct Correlation {
pub:
	eq   int // equation number
	a    f64
	b    f64
	c    f64
	d    f64
	e    f64
	tmin f64 // K; the range the fit is valid over
	tmax f64 // K
	tc   f64 // K; only for the reduced-temperature forms (106, 114, 116)
}

// is_supported says whether equation `eq` can be evaluated.
pub fn is_supported(eq int) bool {
	return eq in supported
}

// new builds a correlation from equation number, coefficients A, B, … (up
// to five; missing ones are zero), and its temperature range. tc is needed
// only by forms 106, 114 and 116; pass 0 otherwise.
//
// It refuses forms it can't evaluate (such as ChemSep's eq. 200), so a
// Correlation can always be evaluated, and ranges that run past Tc for 114
// and 116, which aren't defined there. Form 106 may run past Tc (ChemSep
// has 101 such records); see eval.
pub fn new(eq int, coeffs []f64, tmin f64, tmax f64, tc f64) !Correlation {
	if !is_supported(eq) {
		return error('equation ${eq} is not supported')
	}
	if coeffs.len > 5 {
		return error('equation ${eq}: at most 5 coefficients, got ${coeffs.len}')
	}
	if eq in needs_tc && tc <= 0.0 {
		return error('equation ${eq} needs a critical temperature')
	}
	// Some ChemSep solid-phase records start at 0 K; that's a fine bound.
	if !(tmin >= 0.0 && tmin <= tmax) {
		return error('equation ${eq}: bad temperature range [${tmin}, ${tmax}]')
	}
	if eq in [114, 116] && tmax > tc {
		return error('equation ${eq}: range ends at ${tmax} K, above Tc = ${tc} K')
	}
	k := []f64{len: 5, init: if index < coeffs.len { coeffs[index] } else { 0.0 }}
	return Correlation{
		eq:   eq
		a:    k[0]
		b:    k[1]
		c:    k[2]
		d:    k[3]
		e:    k[4]
		tmin: tmin
		tmax: tmax
		tc:   tc
	}
}

// unused_nonzero says whether coeffs has non-zero values beyond what
// equation eq uses. That usually means the equation number in the source is
// wrong, so the importer reports it.
pub fn unused_nonzero(eq int, coeffs []f64) bool {
	n := used_coeffs[eq] or { return false }
	return coeffs.len > n && coeffs[n..].any(it != 0.0)
}

// in_range says whether t lies inside the fit's temperature range.
// Evaluating outside it extrapolates; the caller decides whether that's
// acceptable.
pub fn (c Correlation) in_range(t f64) bool {
	return t >= c.tmin && t <= c.tmax
}

// poly evaluates the polynomial with coefficients a, b, c, d, e (lowest
// power first), truncated to n terms, by Horner's rule.
fn (c Correlation) poly(t f64, n int) f64 {
	k := [c.a, c.b, c.c, c.d, c.e]
	mut s := 0.0
	for i := n - 1; i >= 0; i-- {
		s = s * t + k[i]
	}
	return s
}

// dpoly is the derivative of poly.
fn (c Correlation) dpoly(t f64, n int) f64 {
	k := [c.a, c.b, c.c, c.d, c.e]
	mut s := 0.0
	for i := n - 1; i >= 1; i-- {
		s = s * t + f64(i) * k[i]
	}
	return s
}

// eval returns Y(T).
//
// Form 106 describes properties that vanish at Tc (heat of vaporization,
// surface tension), so it returns 0 at and above Tc. Forms 114 and 116
// aren't defined above Tc and return NaN there.
pub fn (c Correlation) eval(t f64) f64 {
	return match c.eq {
		1, 2, 3, 4, 5 {
			c.poly(t, c.eq)
		}
		10 {
			math.exp(c.a - c.b / (t + c.c))
		}
		11, 12, 13, 14, 15 {
			math.exp(c.poly(t, c.eq - 10))
		}
		16 {
			c.a + math.exp(c.b / t + c.c + c.d * t + c.e * t * t)
		}
		17 {
			c.a + math.exp(c.b + c.c * t + c.d * t * t + c.e * t * t * t)
		}
		45 {
			// ∫ of the eq. 5 polynomial: AT + BT²/2 + CT³/3 + DT⁴/4 + ET⁵/5
			t * (c.a + t * (c.b / 2.0 + t * (c.c / 3.0 + t * (c.d / 4.0 + t * c.e / 5.0))))
		}
		100 {
			c.poly(t, 5)
		}
		101 {
			math.exp(c.a + c.b / t + c.c * math.log(t) + c.d * math.pow(t, c.e))
		}
		102 {
			c.a * math.pow(t, c.b) / (1.0 + c.c / t + c.d / (t * t))
		}
		104 {
			c.a + c.b / t + c.c / math.pow(t, 3) + c.d / math.pow(t, 8) + c.e / math.pow(t, 9)
		}
		105 {
			c.a / math.pow(c.b, 1.0 + math.pow(1.0 - t / c.c, c.d))
		}
		106 {
			tr := t / c.tc
			if tr >= 1.0 {
				0.0
			} else {
				c.a * math.pow(1.0 - tr, c.b + tr * (c.c + tr * (c.d + tr * c.e)))
			}
		}
		107 {
			c.a + c.b * x_over_sinh_sq(c.c / t) + c.d * x_over_cosh_sq(c.e / t)
		}
		114 {
			tau := 1.0 - t / c.tc
			if tau < 0.0 {
				math.nan()
			} else {
				c.a * c.a / tau + c.b - 2.0 * c.a * c.c * tau - c.a * c.d * tau * tau - c.c * c.c * math.pow(tau,
					3) / 3.0 - c.c * c.d * math.pow(tau, 4) / 2.0 - c.d * c.d * math.pow(tau, 5) / 5.0
			}
		}
		115 {
			math.exp(c.a + c.b / t + c.c * math.log(t) + c.d * t * t + c.e / (t * t))
		}
		116 {
			tau := 1.0 - t / c.tc
			if tau < 0.0 {
				math.nan()
			} else {
				c.a + c.b * math.pow(tau, 0.35) + c.c * math.pow(tau, 2.0 / 3.0) + c.d * tau +
					c.e * math.pow(tau, 4.0 / 3.0)
			}
		}
		else {
			math.nan() // unreachable: new() only accepts supported forms
		}
	}
}

// deriv returns dY/dT, in closed form for every supported equation. Above
// Tc it follows eval: 0 for form 106, NaN for 114 and 116.
pub fn (c Correlation) deriv(t f64) f64 {
	return match c.eq {
		1, 2, 3, 4, 5 {
			c.dpoly(t, c.eq)
		}
		10 {
			c.eval(t) * c.b / ((t + c.c) * (t + c.c))
		}
		11, 12, 13, 14, 15 {
			c.eval(t) * c.dpoly(t, c.eq - 10)
		}
		16 {
			math.exp(c.b / t + c.c + c.d * t + c.e * t * t) * (-c.b / (t * t) + c.d + 2.0 * c.e * t)
		}
		17 {
			math.exp(c.b + c.c * t + c.d * t * t + c.e * t * t * t) * (c.c + 2.0 * c.d * t +
				3.0 * c.e * t * t)
		}
		45 {
			c.poly(t, 5)
		}
		100 {
			c.dpoly(t, 5)
		}
		101 {
			c.eval(t) * (-c.b / (t * t) + c.c / t + c.d * c.e * math.pow(t, c.e - 1.0))
		}
		102 {
			den := 1.0 + c.c / t + c.d / (t * t)
			c.eval(t) * (c.b / t + (c.c / (t * t) + 2.0 * c.d / (t * t * t)) / den)
		}
		104 {
			-c.b / (t * t) - 3.0 * c.c / math.pow(t, 4) - 8.0 * c.d / math.pow(t, 9) - 9.0 * c.e / math.pow(t,
				10)
		}
		105 {
			// Y = A·B^−(1+τᴰ), τ = 1 − T/C
			tau := 1.0 - t / c.c
			c.eval(t) * math.log(c.b) * c.d * math.pow(tau, c.d - 1.0) / c.c
		}
		106 {
			// ln Y = ln A + h·ln(1 − Tr), h = B + C·Tr + D·Tr² + E·Tr³
			tr := t / c.tc
			if tr >= 1.0 {
				0.0
			} else {
				h := c.b + tr * (c.c + tr * (c.d + tr * c.e))
				dh := (c.c + tr * (2.0 * c.d + 3.0 * c.e * tr)) / c.tc
				c.eval(t) * (dh * math.log(1.0 - tr) - h / (c.tc * (1.0 - tr)))
			}
		}
		107 {
			// d/dT[x²/sinh²x] with x = C/T is −2x²/sinh²x·(1 − x·coth x)/T,
			// and likewise for the cosh term with tanh.
			x := c.c / t
			y := c.e / t
			-2.0 * c.b * x_over_sinh_sq(x) * (1.0 - x_coth_x(x)) / t - 2.0 * c.d * x_over_cosh_sq(y) * (1.0 - y * math.tanh(y)) / t
		}
		114 {
			tau := 1.0 - t / c.tc
			if tau < 0.0 {
				math.nan()
			} else {
				dy_dtau := -c.a * c.a / (tau * tau) - 2.0 * c.a * c.c - 2.0 * c.a * c.d * tau - c.c * c.c * tau * tau - 2.0 * c.c * c.d * math.pow(tau,
					3) - c.d * c.d * math.pow(tau, 4)
				-dy_dtau / c.tc
			}
		}
		115 {
			c.eval(t) * (-c.b / (t * t) + c.c / t + 2.0 * c.d * t - 2.0 * c.e / (t * t * t))
		}
		116 {
			tau := 1.0 - t / c.tc
			if tau < 0.0 {
				math.nan()
			} else {
				dy_dtau := 0.35 * c.b * math.pow(tau, -0.65) + 2.0 / 3.0 * c.c * math.pow(tau, -1.0 / 3.0) +
					c.d + 4.0 / 3.0 * c.e * math.pow(tau, 1.0 / 3.0)
				-dy_dtau / c.tc
			}
		}
		else {
			math.nan()
		}
	}
}

// integral returns ∫ Y dT from t1 to t2, in closed form where one exists
// and by adaptive quadrature otherwise. Enthalpies come from this. It fails
// where Y isn't finite on the interval (114 at Tc, or 114 and 116 above it).
pub fn (c Correlation) integral(t1 f64, t2 f64) !f64 {
	return match c.eq {
		1, 2, 3, 4, 5, 45, 100, 104, 107, 114, 116 {
			finite(c.antiderivative(t2) - c.antiderivative(t1), c.eq, t1, t2)!
		}
		else {
			numeric.integrate(c, t1, t2, integral_rtol)!
		}
	}
}

// integral_over_t returns ∫ Y/T dT from t1 to t2. Entropies come from this.
pub fn (c Correlation) integral_over_t(t1 f64, t2 f64) !f64 {
	return match c.eq {
		1, 2, 3, 4, 5, 45, 100, 104, 107 {
			finite(c.antiderivative_over_t(t2) - c.antiderivative_over_t(t1), c.eq, t1, t2)!
		}
		else {
			numeric.integrate(OverT{c}, t1, t2, integral_rtol)!
		}
	}
}

// at is Y(T), which makes a Correlation a numeric.Function to integrate.
pub fn (c Correlation) at(t f64) f64 {
	return c.eval(t)
}

// OverT is Y(T)/T, for the integrals entropies need.
struct OverT {
	c Correlation
}

fn (f OverT) at(t f64) f64 {
	return f.c.eval(t) / t
}

// finite passes x through, or fails if a closed form came out non-finite,
// matching numeric.integrate, which fails the same way.
fn finite(x f64, eq int, t1 f64, t2 f64) !f64 {
	if !math.is_finite(x) {
		return error('equation ${eq}: integral is not finite on [${t1}, ${t2}]')
	}
	return x
}

// antiderivative is a closed-form ∫ Y dT, for the forms that have one.
fn (c Correlation) antiderivative(t f64) f64 {
	return match c.eq {
		1, 2, 3, 4, 5, 100 {
			n := if c.eq == 100 { 5 } else { c.eq }
			k := [c.a, c.b, c.c, c.d, c.e]
			mut s := 0.0
			for i := n - 1; i >= 0; i-- {
				s = s * t + k[i] / f64(i + 1)
			}
			s * t
		}
		45 {
			t * t * (c.a / 2.0 + t * (c.b / 6.0 + t * (c.c / 12.0 + t * (c.d / 20.0 + t * c.e / 30.0))))
		}
		104 {
			c.a * t + c.b * math.log(t) - c.c / (2.0 * t * t) - c.d / (7.0 * math.pow(t, 7)) - c.e / (8.0 * math.pow(t,
				8))
		}
		107 {
			// d/dT[C·coth(C/T)] = (C/T)²/sinh²(C/T); d/dT[−E·tanh(E/T)] =
			// (E/T)²/cosh²(E/T). x·coth x is written as C·coth(C/T)/T·T.
			mut s := c.a * t
			if c.b != 0.0 {
				s += c.b * t * x_coth_x(c.c / t)
			}
			if c.d != 0.0 {
				s -= c.d * c.e * math.tanh(c.e / t)
			}
			s
		}
		114 {
			// In τ = 1 − T/Tc, dT = −Tc·dτ.
			tau := 1.0 - t / c.tc
			-c.tc * (c.a * c.a * math.log(tau) + c.b * tau - c.a * c.c * tau * tau - c.a * c.d * math.pow(tau,
				3) / 3.0 - c.c * c.c * math.pow(tau, 4) / 12.0 - c.c * c.d * math.pow(tau, 5) / 10.0 - c.d * c.d * math.pow(tau,
				6) / 30.0)
		}
		116 {
			tau := 1.0 - t / c.tc
			-c.tc * (c.a * tau + c.b * math.pow(tau, 1.35) / 1.35 + c.c * math.pow(tau, 5.0 / 3.0) * 0.6 +
				c.d * tau * tau / 2.0 + c.e * math.pow(tau, 7.0 / 3.0) * 3.0 / 7.0)
		}
		else {
			math.nan()
		}
	}
}

// antiderivative_over_t is a closed-form ∫ Y/T dT, for the forms that have one.
fn (c Correlation) antiderivative_over_t(t f64) f64 {
	return match c.eq {
		1, 2, 3, 4, 5, 100 {
			n := if c.eq == 100 { 5 } else { c.eq }
			k := [c.a, c.b, c.c, c.d, c.e]
			mut s := 0.0
			for i := n - 1; i >= 1; i-- {
				s = s * t + k[i] / f64(i)
			}
			k[0] * math.log(t) + s * t
		}
		45 {
			t * (c.a + t * (c.b / 4.0 + t * (c.c / 9.0 + t * (c.d / 16.0 + t * c.e / 25.0))))
		}
		104 {
			c.a * math.log(t) - c.b / t - c.c / (3.0 * math.pow(t, 3)) - c.d / (8.0 * math.pow(t, 8)) - c.e / (9.0 * math.pow(t,
				9))
		}
		107 {
			// With x = C/T: d/dT[x·coth x − ln sinh x] = x²/(T·sinh²x); with
			// y = E/T: d/dT[−(y·tanh y − ln cosh y)] = y²/(T·cosh²y).
			mut s := c.a * math.log(t)
			if c.b != 0.0 && c.c == 0.0 {
				s += c.b * math.log(t) // (x/sinh x)² → 1 as C → 0
			} else if c.b != 0.0 {
				x := c.c / t
				s += c.b * (x_coth_x(x) - ln_sinh(x))
			}
			if c.d != 0.0 {
				y := c.e / t
				s -= c.d * (y * math.tanh(y) - ln_cosh(y))
			}
			s
		}
		else {
			math.nan()
		}
	}
}

// Overflow-safe pieces of the Aly–Lee form (DIPPR 107). sinh and cosh
// overflow above about 710, which C/T reaches at low T, so these are
// written in terms of e^(−|x|). All of them are even in x (ln_sinh takes
// |sinh x|), so a negative C or E is harmless too.

// x_over_sinh_sq returns (x/sinh x)².
fn x_over_sinh_sq(x f64) f64 {
	ax := math.abs(x)
	if ax == 0.0 {
		return 1.0
	}
	if ax < 1.0 {
		r := ax / math.sinh(ax)
		return r * r
	}
	e := math.exp(-ax)
	r := 2.0 * ax * e / (1.0 - e * e)
	return r * r
}

// x_over_cosh_sq returns (x/cosh x)².
fn x_over_cosh_sq(x f64) f64 {
	ax := math.abs(x)
	e := math.exp(-ax)
	r := 2.0 * ax * e / (1.0 + e * e)
	return r * r
}

// x_coth_x returns x·coth x, which tends to 1 as x → 0.
fn x_coth_x(x f64) f64 {
	if x == 0.0 {
		return 1.0
	}
	return x / math.tanh(x)
}

// ln_sinh returns ln|sinh x| for x ≠ 0.
fn ln_sinh(x f64) f64 {
	ax := math.abs(x)
	if ax < 1.0 {
		return math.log(math.sinh(ax))
	}
	return ax - math.ln2 + math.log1p(-math.exp(-2.0 * ax))
}

// ln_cosh returns ln(cosh x). Near 0, cosh x − 1 = 2·sinh²(x/2) keeps its
// relative precision.
fn ln_cosh(x f64) f64 {
	ax := math.abs(x)
	if ax < 1.0 {
		s := math.sinh(0.5 * ax)
		return math.log1p(2.0 * s * s)
	}
	return ax - math.ln2 + math.log1p(math.exp(-2.0 * ax))
}
