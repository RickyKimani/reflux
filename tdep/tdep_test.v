module tdep

import math
import numeric

fn test_new_refuses_what_it_cant_evaluate() {
	// ChemSep's eq. 200 (5 vapour pressures) and DIPPR 124 (1 liquid Cp)
	// are the forms in the data that reflux doesn't evaluate.
	for eq in [0, 6, 118, 120, 121, 124, 200] {
		if _ := new(eq, [1.0], 200.0, 300.0, 400.0) {
			assert false, 'equation ${eq} should be refused'
		} else {
			assert err.msg() == 'equation ${eq} is not supported'
		}
	}
	if _ := new(100, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 200.0, 300.0, 0.0) {
		assert false, 'six coefficients should be refused'
	}
	for eq in needs_tc {
		if _ := new(eq, [1.0, 0.5], 200.0, 300.0, 0.0) {
			assert false, 'equation ${eq} without Tc should be refused'
		}
	}
	for r in [[-1.0, 300.0], [300.0, 200.0], [math.nan(), 300.0], [200.0, math.nan()]] {
		if _ := new(100, [1.0], r[0], r[1], 0.0) {
			assert false, 'range ${r} should be refused'
		}
	}
	// 114 and 116 aren't defined above Tc, so a range past it is refused;
	// 106 is extended by zero there and allowed (ChemSep has 101 such records).
	for eq in [114, 116] {
		if _ := new(eq, [1.0, 1.0, 1.0, 1.0], 200.0, 401.0, 400.0) {
			assert false, 'equation ${eq} past Tc should be refused'
		}
		new(eq, [1.0, 1.0, 1.0, 1.0], 200.0, 400.0, 400.0)!
	}
	new(106, [1.0, 0.4], 200.0, 401.0, 400.0)!
	// A range starting at 0 K is fine, and so is a single temperature.
	new(1, [1.0], 0.0, 5.0, 0.0)!
	new(1, [1.0], 298.15, 298.15, 0.0)!
}

fn test_missing_coefficients_are_zero() {
	c := new(100, [1.0, 2.0], 200.0, 300.0, 0.0)!
	assert c.c == 0.0 && c.d == 0.0 && c.e == 0.0
	assert c.eval(250.0) == 501.0
	assert c.deriv(250.0) == 2.0
}

fn test_unused_nonzero_flags_suspect_records() {
	// The two kinds ChemSep has: eq. 2 solid densities with a C, and eq. 4
	// liquid heat capacities (amines) with an E.
	assert unused_nonzero(2, [44.0, -0.1, 3.0, 0.0, 0.0])
	assert unused_nonzero(4, [78653.4, 311.4, -0.06, -2.7e-6, -1.8e-5])
	assert unused_nonzero(10, [20.0, 700.0, 4.0, 1.0])
	assert unused_nonzero(105, [2.6, 0.25, 132.5, 0.27, 1.0])
	assert !unused_nonzero(4, [78653.4, 311.4, -0.06, -2.7e-6, 0.0])
	assert !unused_nonzero(100, [1.0, 2.0, 3.0, 4.0, 5.0])
	assert !unused_nonzero(2, [44.0])
	assert !unused_nonzero(200, [1.0, 2.0, 3.0, 4.0, 5.0])
}

fn test_in_range_includes_both_ends() {
	c := new(100, [1.0], 200.0, 300.0, 0.0)!
	assert c.in_range(200.0) && c.in_range(250.0) && c.in_range(300.0)
	assert !c.in_range(199.999) && !c.in_range(300.001)
}

fn test_integrals_change_sign_when_reversed() {
	// A closed form (100) and a quadrature (101).
	for c in [new(100, [29562.3, -7.16, 0.0216, -1.4e-5, 2.9e-9], 50.0, 1500.0, 0.0)!,
		new(101, [66.343, -4363.2, -7.046, 9.4509e-6, 2.0], 134.86, 425.12, 0.0)!] {
		assert c.integral(300.0, 400.0)! == -c.integral(400.0, 300.0)!
		assert math.tolerance(c.integral_over_t(300.0, 400.0)!, -c.integral_over_t(400.0, 300.0)!,
			1e-15)
		assert c.integral(350.0, 350.0)! == 0.0
	}
}

fn test_closed_forms_agree_with_quadrature() {
	// Independent of the fixtures: each closed-form antiderivative against
	// numeric.integrate of eval, on made-up coefficients that use every term.
	cases := [
		new(5, [3.0e4, -10.0, 0.05, -3e-5, 6e-9], 200.0, 1000.0, 0.0)!,
		new(45, [3.0e4, 50.0, 0.1, -1e-4, 3e-8], 250.0, 800.0, 0.0)!,
		new(100, [-8018.8, 1284.8, -11.519, 0.046329, -6.5817e-05], 15.0, 265.9, 0.0)!,
		// In real second virial coefficients the D and E terms are ~1e-14 of
		// the total, invisible to a test; here every term counts.
		new(104, [0.5, -1.0, -0.3, 0.05, -0.01], 1.0, 3.0, 0.0)!,
		new(107, [33298.3, 26618.3, 2572.35, 8630.72, 1138.84], 50.0, 2000.0, 0.0)!,
		new(114, [67.155, 105580.0, 310.21, -490.54], 154.56, 359.98, 386.6)!,
		new(116, [0.3, 0.8, -0.4, 0.2, 0.1], 250.0, 700.0, 726.0)!,
	]
	for c in cases {
		a, b := c.tmin, c.tmax
		want := numeric.integrate(fn [c] (t f64) f64 {
			return c.eval(t)
		}, a, b, 1e-14)!
		assert math.tolerance(c.integral(a, b)!, want, 1e-12), 'eq ${c.eq}: ${c.integral(a, b)!} vs ${want}'
		if c.eq !in [114, 116] {
			want_t := numeric.integrate(fn [c] (t f64) f64 {
				return c.eval(t) / t
			}, a, b, 1e-14)!
			assert math.tolerance(c.integral_over_t(a, b)!, want_t, 1e-12), 'eq ${c.eq}: ${c.integral_over_t(a, b)!} vs ${want_t}'
		}
	}
}

fn test_form_106_vanishes_at_and_above_tc() {
	// n-butane's ChemSep heat of vaporization, whose range ends at 425.16 K,
	// past its Tc of 425.12 K.
	tc := 425.12
	c := new(106, [3.6258e7, 0.83741, -0.83676, 0.41526, -0.007606], 134.86, 425.16, tc)!
	assert c.eval(tc) == 0.0 && c.eval(425.16) == 0.0 && c.eval(500.0) == 0.0
	assert c.deriv(tc) == 0.0 && c.deriv(500.0) == 0.0
	assert c.eval(tc - 1e-6) > 0.0
	// Integrating across Tc adds nothing.
	below := c.integral(400.0, tc)!
	assert math.tolerance(c.integral(400.0, 425.16)!, below, 1e-12)
	assert math.tolerance(c.integral(400.0, 450.0)!, below, 1e-12)
}

fn test_forms_114_and_116_are_undefined_above_tc() {
	tc := 386.6
	c114 := new(114, [67.155, 105580.0, 310.21, -490.54], 154.56, 359.98, tc)!
	c116 := new(116, [0.3, 0.8, -0.4, 0.2, 0.1], 250.0, 386.6, tc)!
	for c in [c114, c116] {
		assert math.is_nan(c.eval(tc + 1.0))
		assert math.is_nan(c.deriv(tc + 1.0))
		if _ := c.integral(300.0, tc + 1.0) {
			assert false, 'eq ${c.eq}: integrating past Tc should fail'
		}
	}
	// 114 has A²/τ, so it's infinite at Tc and its integral to Tc diverges.
	assert math.is_inf(c114.eval(tc), 1)
	if _ := c114.integral(300.0, tc) {
		assert false, 'the 114 integral to Tc diverges'
	}
	// 116 is finite at Tc: Y(Tc) = A.
	assert c116.eval(tc) == 0.3
	c116.integral(300.0, tc)!
}

fn test_aly_lee_survives_extreme_temperatures() {
	// C/T and E/T reach 10⁶ at 1 mK, where sinh and cosh overflow a float.
	// Both hyperbolic terms vanish there, leaving Y = A.
	c := new(107, [33298.3, 26618.3, 2572.35, 8630.72, 1138.84], 50.0, 2000.0, 0.0)!
	for t in [1e-3, 0.1, 1.0, 1e6] {
		for v in [c.eval(t), c.deriv(t), c.integral(t, 300.0)!, c.integral_over_t(t, 300.0)!] {
			assert math.is_finite(v), 'not finite at ${t} K'
		}
	}
	assert c.eval(1e-3) == 33298.3
	assert c.deriv(1e-3) == 0.0
	// At high T both terms tend to their limits: (x/sinh x)² → 1 and
	// (x/cosh x)² → 0, so Y → A + B.
	assert math.tolerance(c.eval(1e9), 33298.3 + 26618.3, 1e-9)
}

fn test_aly_lee_edge_coefficients() {
	// C = 0 makes the sinh term a constant B, and E = 0 drops the cosh
	// term: Y = A + B everywhere.
	c := new(107, [100.0, 20.0, 0.0, 5.0, 0.0], 100.0, 1000.0, 0.0)!
	assert c.eval(300.0) == 120.0
	assert c.deriv(300.0) == 0.0
	assert math.tolerance(c.integral(200.0, 400.0)!, 120.0 * 200.0, 1e-14)
	assert math.tolerance(c.integral_over_t(200.0, 400.0)!, 120.0 * math.log(2.0), 1e-14)
	// Every term is even in C and in E, so negative ones change nothing.
	pos := new(107, [33298.3, 26618.3, 2572.35, 8630.72, 1138.84], 50.0, 2000.0, 0.0)!
	neg := new(107, [33298.3, 26618.3, -2572.35, 8630.72, -1138.84], 50.0, 2000.0, 0.0)!
	for t in [60.0, 300.0, 1500.0] {
		assert pos.eval(t) == neg.eval(t)
		assert pos.deriv(t) == neg.deriv(t)
	}
	assert math.tolerance(pos.integral(100.0, 1000.0)!, neg.integral(100.0, 1000.0)!, 1e-14)
	assert math.tolerance(pos.integral_over_t(100.0, 1000.0)!, neg.integral_over_t(100.0, 1000.0)!,
		1e-14)
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

fn test_hyperbolic_helpers_agree_with_the_naive_forms() {
	// Moved here from thermo in M2 step 5, when thermo.AlyLee began using
	// these helpers through tdep instead of keeping copies.
	for x in [0.001, 0.01, 0.5, 0.999, 1.0, 1.001, 3.0, 20.0, 300.0] {
		assert math.tolerance(x_over_sinh_sq(x), math.pow(x / math.sinh(x), 2), 1e-13), 'x/sinh at ${x}'
		assert math.tolerance(x_over_cosh_sq(x), math.pow(x / math.cosh(x), 2), 1e-13), 'x/cosh at ${x}'
		assert math.tolerance(ln_sinh(x), math.log(math.sinh(x)), 1e-13), 'ln sinh at ${x}'
		assert math.tolerance(ln_cosh(x), ref_ln_cosh(x), 1e-13), 'ln cosh at ${x}: ${ln_cosh(x)} vs ${ref_ln_cosh(x)}'
	}
	assert x_over_sinh_sq(2000.0) == 0.0 // underflows cleanly instead of Inf/Inf
}

fn test_hyperbolic_helpers_near_zero_and_switchover() {
	// Taylor series near 0, and continuity where the helpers switch method
	// at x = 1.
	for x in [1e-8, 1e-4] {
		assert math.tolerance(x_over_sinh_sq(x), 1.0 - x * x / 3.0, 1e-15)
		assert math.tolerance(x_over_cosh_sq(x), x * x * (1.0 - x * x), 1e-15)
		assert math.tolerance(x_coth_x(x), 1.0 + x * x / 3.0, 1e-15)
		assert math.tolerance(ln_sinh(x), math.log(x) + x * x / 6.0, 1e-15)
		assert math.tolerance(ln_cosh(x), x * x / 2.0 - math.pow(x, 4) / 12.0, 1e-15)
	}
	for f in [x_over_sinh_sq, ln_sinh, ln_cosh] {
		assert math.tolerance(f(1.0 - 1e-12), f(1.0), 1e-11)
	}
	// Large x: sinh would overflow, the helpers don't.
	assert x_over_sinh_sq(1000.0) == 0.0
	assert x_over_cosh_sq(1000.0) == 0.0
	assert math.tolerance(ln_sinh(1000.0), 1000.0 - math.ln2, 1e-15)
	assert math.tolerance(ln_cosh(1000.0), 1000.0 - math.ln2, 1e-15)
}
