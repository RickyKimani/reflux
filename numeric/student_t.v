module numeric

import math

// student_t_cdf is P(T ≤ t) for Student's t with ν degrees of freedom:
// 1 − ½·I_x(ν/2, ½) for t ≥ 0, x = ν/(ν + t²), and by symmetry below zero,
// I the regularized incomplete beta function. 1 − x is taken as t²/(ν + t²),
// not by subtraction, which would lose its digits where x is near 1.
pub fn student_t_cdf(t f64, nu f64) f64 {
	x := nu / (nu + t * t)
	y := t * t / (nu + t * t)
	tail := 0.5 * incomplete_beta(0.5 * nu, 0.5, x, y)
	return if t >= 0.0 { 1.0 - tail } else { tail }
}

// student_t_quantile is the t with P(T ≤ t) = p, for p in (0, 1) and ν > 0:
// the 95% interval's half-width over a standard error is the quantile at
// 0.975 with ν = n − p (M8's D14). Bisection on the CDF, to the doubles'
// resolution of t.
pub fn student_t_quantile(p f64, nu f64) !f64 {
	if !(p > 0.0 && p < 1.0) || !(nu > 0.0) || !math.is_finite(nu) {
		return error('student_t_quantile: p must be in (0, 1) and ν positive, not ${p} and ${nu}')
	}
	if p == 0.5 {
		return 0.0
	}
	if p < 0.5 {
		return -student_t_quantile(1.0 - p, nu)!
	}
	mut lo := 0.0
	mut hi := 1.0
	for student_t_cdf(hi, nu) < p {
		lo = hi
		hi *= 2.0
		if hi > 1e300 {
			return error('student_t_quantile: no t below 1e300 has P = ${p} at ν = ${nu}')
		}
	}
	for _ in 0 .. 200 {
		mid := 0.5 * (lo + hi)
		if mid == lo || mid == hi {
			break
		}
		if student_t_cdf(mid, nu) < p {
			lo = mid
		} else {
			hi = mid
		}
	}
	return 0.5 * (lo + hi)
}

// incomplete_beta is I_x(a, b), the regularized incomplete beta function,
// given x and y = 1 − x, by its continued fraction (Numerical Recipes, §6.4,
// by Lentz's method) on the side of x where it converges fast, and the
// symmetry I_x(a, b) = 1 − I_y(b, a) on the other. Its front factor,
// x^a·y^b/(a·B(a, b)), takes each logarithm from the smaller of x and y.
fn incomplete_beta(a f64, b f64, x f64, y f64) f64 {
	if x <= 0.0 {
		return 0.0
	}
	if y <= 0.0 {
		return 1.0
	}
	ln_x := if y < 0.5 { math.log1p(-y) } else { math.log(x) }
	ln_y := if x < 0.5 { math.log1p(-x) } else { math.log(y) }
	front := math.exp(ln_beta_inverse(a, b) + a * ln_x + b * ln_y)
	if x < (a + 1.0) / (a + b + 2.0) {
		return front * beta_fraction(a, b, x) / a
	}
	return 1.0 - front * beta_fraction(b, a, y) / b
}

// ln_beta_inverse is ln(Γ(a + b)/(Γ(a)·Γ(b))), the larger of a and b's part
// by Stirling's series: ln Γ(z + c) − ln Γ(z) is a difference of numbers far
// larger than itself for a large z, and the series takes it with log1p.
fn ln_beta_inverse(a f64, b f64) f64 {
	big, small := if a >= b { a, b } else { b, a }
	if big < 10.0 {
		return math.log_gamma(a + b) - math.log_gamma(a) - math.log_gamma(b)
	}
	// ln Γ(z) = (z − ½)·ln z − z + ½·ln 2π + 1/(12z) − 1/(360z³) + 1/(1260z⁵) − …
	s := big + small
	diff := (big - 0.5) * math.log1p(small / big) + small * math.log(s) - small +
		stirling_tail(s) - stirling_tail(big)
	return diff - math.log_gamma(small)
}

// stirling_tail is Stirling's series' terms in 1/z, for z ≥ 10.
fn stirling_tail(z f64) f64 {
	z2 := z * z
	return (1.0 / 12.0 - (1.0 / 360.0 - (1.0 / 1260.0 - 1.0 / (1680.0 * z2)) / z2) / z2) / z
}

// beta_fraction is the incomplete beta's continued fraction at x, to the
// doubles' resolution.
fn beta_fraction(a f64, b f64, x f64) f64 {
	tiny := 1e-300
	qab := a + b
	qap := a + 1.0
	qam := a - 1.0
	mut c := 1.0
	mut d := 1.0 - qab * x / qap
	if math.abs(d) < tiny {
		d = tiny
	}
	d = 1.0 / d
	mut h := d
	for m in 1 .. 1000 {
		fm := f64(m)
		m2 := 2.0 * fm
		mut aa := fm * (b - fm) * x / ((qam + m2) * (a + m2))
		d = 1.0 + aa * d
		if math.abs(d) < tiny {
			d = tiny
		}
		c = 1.0 + aa / c
		if math.abs(c) < tiny {
			c = tiny
		}
		d = 1.0 / d
		h *= d * c
		aa = -(a + fm) * (qab + fm) * x / ((a + m2) * (qap + m2))
		d = 1.0 + aa * d
		if math.abs(d) < tiny {
			d = tiny
		}
		c = 1.0 + aa / c
		if math.abs(c) < tiny {
			c = tiny
		}
		d = 1.0 / d
		del := d * c
		h *= del
		if math.abs(del - 1.0) < 1e-16 {
			break
		}
	}
	return h
}
