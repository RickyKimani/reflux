// Internal test (module thermo), so the overflow-safe helpers can be checked
// directly as well as through the public CpForm methods.
module thermo

import math

// Cp/R = A + B·T + C·T² + D·T⁻², the four-constant form, as Polynomial terms.
fn four_constant(a f64, b f64, c f64, d f64) CpForm {
	return Polynomial{
		terms: [
			PowerTerm{r_gas * a, 0},
			PowerTerm{r_gas * b, 1},
			PowerTerm{r_gas * c, 2},
			PowerTerm{r_gas * d, -2},
		]
	}
}

// Forms to test. The four-constant values are Smith, Van Ness & Abbott
// Table C.1 (n-butane and CO₂). The Aly–Lee coefficients are illustrative,
// of the size DIPPR tabulates for a light hydrocarbon. The last polynomial
// has a T⁻¹ term and a constant, the two powers whose integrals are
// logarithms.
fn forms() map[string]CpForm {
	return {
		'n-butane poly': four_constant(1.935, 36.915e-3, -11.402e-6, 0.0)
		'CO2 poly':      four_constant(5.457, 1.045e-3, 0.0, -1.157e5)
		'Aly-Lee':       CpForm(AlyLee{
			a: 80.0
			b: 240.0
			c: 1630.0
			d: 150.0
			e: 730.0
		})
		'log powers':    CpForm(Polynomial{
			terms: [PowerTerm{30.0, 0}, PowerTerm{2000.0, -1}, PowerTerm{0.01, 1}]
		})
	}
}

// simpson integrates f over [a, b] with n (even) intervals.
fn simpson(f fn (f64) f64, a f64, b f64, n int) f64 {
	h := (b - a) / f64(n)
	mut s := f(a) + f(b)
	for i in 1 .. n {
		s += f(a + f64(i) * h) * if i % 2 == 1 { 4.0 } else { 2.0 }
	}
	return s * h / 3.0
}

fn test_cp_integrals_match_quadrature() {
	for name, form in forms() {
		for r in [[50.0, 300.0], [298.15, 1500.0], [200.0, 5000.0]] {
			t1, t2 := r[0], r[1]
			q_h := simpson(fn [form] (t f64) f64 {
				return form.cp(t)
			}, t1, t2, 20000)
			q_s := simpson(fn [form] (t f64) f64 {
				return form.cp(t) / t
			}, t1, t2, 20000)
			assert math.tolerance(form.int_cp(t1, t2), q_h, 1e-10), '${name} ∫Cp dT on [${t1}, ${t2}]: ${form.int_cp(t1, t2)} vs ${q_h}'
			assert math.tolerance(form.int_cp_over_t(t1, t2), q_s, 1e-10), '${name} ∫Cp/T dT on [${t1}, ${t2}]'
		}
	}
}

fn test_aly_lee_limits_and_extreme_temperatures() {
	f := CpForm(AlyLee{
		a: 80.0
		b: 240.0
		c: 1630.0
		d: 150.0
		e: 730.0
	})
	// T → 0: both bracketed terms vanish, Cp → A. T → ∞: the sinh term → 1
	// and the cosh term → 0, so Cp → A + B.
	assert math.tolerance(f.cp(1.0), 80.0, 1e-12)
	assert math.tolerance(f.cp(1e6), 80.0 + 240.0, 1e-5)
	// At 1 K, C/T = 1630: naive sinh and cosh would overflow to Inf.
	for t in [1.0, 2.0, 10.0, 5000.0, 1e5] {
		for x in [f.cp(t), f.int_cp(t, 300.0), f.int_cp_over_t(t, 300.0)] {
			assert math.is_finite(x), 'not finite at T = ${t}'
		}
	}
}

// ref_ln_cosh is an accurate reference: the Taylor series where cosh x is
// too close to 1 for log(cosh x) to keep its digits, the naive form elsewhere.
fn ref_ln_cosh(x f64) f64 {
	if x < 0.1 {
		x2 := x * x
		return x2 / 2.0 - x2 * x2 / 12.0 + x2 * x2 * x2 / 45.0 - 17.0 * x2 * x2 * x2 * x2 / 2520.0
	}
	return math.log(math.cosh(x))
}

fn test_overflow_safe_helpers_agree_with_the_naive_forms() {
	for x in [0.001, 0.01, 0.5, 0.999, 1.0, 1.001, 3.0, 20.0, 300.0] {
		assert math.tolerance(x_over_sinh_sq(x), math.pow(x / math.sinh(x), 2), 1e-13), 'x/sinh at ${x}'
		assert math.tolerance(x_over_cosh_sq(x), math.pow(x / math.cosh(x), 2), 1e-13), 'x/cosh at ${x}'
		assert math.tolerance(ln_sinh(x), math.log(math.sinh(x)), 1e-13), 'ln sinh at ${x}'
		assert math.tolerance(ln_cosh(x), ref_ln_cosh(x), 1e-13), 'ln cosh at ${x}: ${ln_cosh(x)} vs ${ref_ln_cosh(x)}'
	}
	assert x_over_sinh_sq(0.0) == 1.0
	assert x_over_sinh_sq(2000.0) == 0.0 // underflows cleanly instead of Inf/Inf
}
