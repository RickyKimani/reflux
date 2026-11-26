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
		'n-butane poly':   four_constant(1.935, 36.915e-3, -11.402e-6, 0.0)
		'CO2 poly':        four_constant(5.457, 1.045e-3, 0.0, -1.157e5)
		'Aly-Lee':         CpForm(AlyLee{
			a: 80.0
			b: 240.0
			c: 1630.0
			d: 150.0
			e: 730.0
		})
		'log powers':      CpForm(Polynomial{
			terms: [PowerTerm{30.0, 0}, PowerTerm{2000.0, -1}, PowerTerm{0.01, 1}]
		})
		// ChemSep's n-butane ideal-gas Cp, eq. 16 in J/(kmol·K): no closed-form
		// integral, so this exercises the quadrature.
		'eq. 16':          CpForm(Correlated{
			c:     tdep.new(16, [44749.9, -338.141, 11.8145, 0.00097744, -3.35913e-07], 80.0,
				1500.0, 0.0) or { panic(err) }
			scale: 1e-3
		})
		'Planck-Einstein': iapws95_like()
	}
}

// iapws95_like is IAPWS-95's ideal-gas Cp (its Table 1, R = 8.314371 J/(mol·K)
// and Tc = 647.096 K), as M5a's water module builds it, but written here
// from the release so that thermo's test needs no other module.
fn iapws95_like() CpForm {
	r := 8.314371357587
	tc := 647.096
	n := [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
	g := [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]
	return PlanckEinstein{
		a:     r * (1.0 + 3.00632)
		terms: []EinsteinTerm{len: 5, init: EinsteinTerm{r * n[index], g[index] * tc}}
	}
}

fn test_planck_einstein_is_aly_lee_s_sinh_term() {
	// (x/sinh x)² = E(2x): a one-term PlanckEinstein with θ = 2C is Aly–Lee
	// with no cosh term, which tdep evaluates and integrates by its own
	// closed forms. Measured: Cp's bits, and 2.5e-15 and 1.6e-15 in the
	// integrals, each a difference of two antiderivatives.
	al := CpForm(AlyLee{
		a: 33.4
		b: 26.8
		c: 2610.5
		d: 0.0
		e: 1.0
	})
	pe := CpForm(PlanckEinstein{
		a:     33.4
		terms: [EinsteinTerm{26.8, 2.0 * 2610.5}]
	})
	for t in [20.0, 100.0, 298.15, 1000.0, 5000.0, 1e5] {
		assert math.tolerance(pe.cp(t), al.cp(t), 1e-14), 'Cp at ${t} K: ${pe.cp(t)} vs ${al.cp(t)}'
		assert math.tolerance(pe.int_cp(298.15, t), al.int_cp(298.15, t), 1e-12)
		assert math.tolerance(pe.int_cp_over_t(298.15, t), al.int_cp_over_t(298.15, t), 1e-12)
	}
}

fn test_planck_einstein_limits_and_extreme_temperatures() {
	f := iapws95_like()
	pe := f as PlanckEinstein
	mut total := pe.a
	for term in pe.terms {
		total += term.c
	}
	// T → 0: every term freezes out, Cp → a. T → ∞: each tends to cᵢ, as
	// E(u) = 1 − u²/12 + …, so at 1e12 K the gap is below 1e-16.
	assert math.tolerance(f.cp(1.0), pe.a, 1e-14)
	assert math.tolerance(f.cp(1e12), total, 1e-14)
	// At 1 K the largest θ/T is 1.8e4, where eᵘ overflows a float.
	for t in [1.0, 2.0, 10.0, 5000.0, 1e5, 1e9] {
		for x in [f.cp(t), f.int_cp(t, 300.0), f.int_cp_over_t(t, 300.0)] {
			assert math.is_finite(x), 'not finite at T = ${t}'
		}
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
