// Internal test (module thermo). The overflow-safe hyperbolic helpers behind
// AlyLee live in tdep now, and are tested there.
module thermo

import math
import tdep

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
		// ChemSep's n-butane ideal-gas Cp, eq. 16 in J/(kmol·K): no closed-form
		// integral, so this exercises the quadrature.
		'eq. 16':        CpForm(Correlated{
			c:     tdep.new(16, [44749.9, -338.141, 11.8145, 0.00097744, -3.35913e-07], 80.0,
				1500.0, 0.0) or { panic(err) }
			scale: 1e-3
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

fn test_correlated_scales_its_correlation() {
	// The same Aly–Lee Cp, once in J/(mol·K) and once as an eq. 107
	// correlation in J/(kmol·K) scaled back down.
	direct := CpForm(AlyLee{
		a: 80.0
		b: 240.0
		c: 1630.0
		d: 150.0
		e: 730.0
	})
	scaled := CpForm(Correlated{
		c:     tdep.new(107, [80e3, 240e3, 1630.0, 150e3, 730.0], 50.0, 1500.0, 0.0)!
		scale: 1e-3
	})
	for t in [100.0, 298.15, 1000.0] {
		assert math.tolerance(scaled.cp(t), direct.cp(t), 1e-14)
		assert math.tolerance(scaled.int_cp(298.15, t), direct.int_cp(298.15, t), 1e-13)
		assert math.tolerance(scaled.int_cp_over_t(298.15, t), direct.int_cp_over_t(298.15, t),
			1e-13)
	}
}

fn test_an_integral_that_cannot_be_done_is_nan() {
	// exp(1000/T) overflows below about 1.4 K, so ∫Cp dT down to 1 K has no
	// finite value. NaN says so; a number here would be wrong.
	f := CpForm(Correlated{
		c: tdep.new(16, [0.0, 1000.0], 100.0, 1000.0, 0.0)!
	})
	assert math.is_nan(f.int_cp(1.0, 300.0))
	assert math.is_nan(f.int_cp_over_t(1.0, 300.0))
	assert math.is_finite(f.int_cp(200.0, 300.0))
}
