// Checks for the α functions: derivatives, normalisation, and the special
// cases where one form reduces to another. An external test.
import math
import thermo { Alpha, MathiasCopeman, Twu }

// The Twu and Mathias–Copeman parameters are illustrative values of
// typical size; the tests check mathematics, not data.
fn alphas() map[string]Alpha {
	return {
		'unity':          Alpha(thermo.Unity{})
		'1/sqrt(Tr)':     thermo.InvSqrtTr{}
		'Soave':          thermo.soave_pr(0.2)
		'PRSV':           thermo.prsv(0.2, 0.05)
		'Twu':            Twu{
			l: 0.3
			m: 0.9
			n: 2.0
		}
		'MC truncated':   MathiasCopeman{
			c1: 0.7
			c2: -0.3
			c3: 0.5
		}
		'MC untruncated': MathiasCopeman{
			c1:        0.7
			c2:        -0.3
			c3:        0.5
			truncated: false
		}
	}
}

fn close(analytic f64, fd f64) bool {
	return math.abs(analytic - fd) <= 1e-6 * math.max(1e-3, math.max(math.abs(analytic),
		math.abs(fd)))
}

fn test_alpha_derivatives_match_finite_differences() {
	for name, a in alphas() {
		// Both sides of Tr = 0.7 (PRSV's κ₁ factor changes sign there) and
		// of Tr = 1 (truncated Mathias–Copeman changes form there). The step
		// never crosses either point.
		for tr in [0.3, 0.6999, 0.7, 0.7001, 0.999, 1.001, 1.5, 4.0] {
			_, d1, d2 := a.eval(tr)
			h := 1e-6 * tr
			up, d1_up, _ := a.eval(tr + h)
			dn, d1_dn, _ := a.eval(tr - h)
			fd1 := (up - dn) / (2.0 * h)
			fd2 := (d1_up - d1_dn) / (2.0 * h)
			assert close(d1, fd1), '${name} at Tr=${tr}: α′ ${d1} vs ${fd1}'
			assert close(d2, fd2), '${name} at Tr=${tr}: α″ ${d2} vs ${fd2}'
		}
	}
}

fn test_every_alpha_is_one_at_the_critical_temperature() {
	// a(Tc) = a_c is what Ωa was fitted to, so α(1) = 1 for every form.
	for name, a in alphas() {
		alpha, _, _ := a.eval(1.0)
		assert math.tolerance(alpha, 1.0, 1e-15), '${name}: α(1) = ${alpha}'
	}
}

fn test_mathias_copeman_with_only_c1_is_soave() {
	soave := Alpha(thermo.Soave{
		m: 0.65
	})
	for truncated in [true, false] {
		mc := Alpha(MathiasCopeman{
			c1:        0.65
			truncated: truncated
		})
		for tr in [0.3, 0.9, 1.0, 1.7] {
			a1, d1, e1 := mc.eval(tr)
			a2, d2, e2 := soave.eval(tr)
			assert math.tolerance(a1, a2, 1e-15) && math.tolerance(d1, d2, 1e-14)
				&& math.tolerance(e1, e2, 1e-14), 'truncated=${truncated} at Tr=${tr}'
		}
	}
}

fn test_twu_with_n_zero_is_unity() {
	twu := Alpha(Twu{
		l: 0.4
		m: 0.8
		n: 0.0
	})
	for tr in [0.3, 1.0, 3.0] {
		a, d1, d2 := twu.eval(tr)
		assert a == 1.0 && d1 == 0.0 && d2 == 0.0, 'Tr=${tr}: ${a}, ${d1}, ${d2}'
	}
}

fn test_twu_stays_positive_and_falls_at_high_tr() {
	twu := Alpha(Twu{
		l: 0.3
		m: 0.9
		n: 2.0
	})
	mut prev := 2.0
	for tr in [1.0, 2.0, 5.0, 10.0, 50.0] {
		a, _, _ := twu.eval(tr)
		assert a > 0.0 && a < prev, 'Tr=${tr}: α = ${a}'
		prev = a
	}
}

fn test_truncated_mathias_copeman_is_smooth_to_first_order_at_tc() {
	// At Tr = 1, u = 0, so dropping c₂ and c₃ leaves α and α′ unchanged but
	// removes 2c₂·u′² from k″. With k = 1 and u′ = −1/2 that makes α″ fall
	// by exactly c₂ across Tc.
	c2 := -0.3
	mc := Alpha(MathiasCopeman{
		c1: 0.7
		c2: c2
		c3: 0.5
	})
	eps := 1e-9
	a_below, d_below, e_below := mc.eval(1.0 - eps)
	a_above, d_above, e_above := mc.eval(1.0 + eps)
	assert math.abs(a_below - a_above) < 1e-8
	assert math.abs(d_below - d_above) < 1e-8
	assert math.abs((e_below - e_above) - c2) < 1e-7, 'α″ jump ${e_below - e_above}, expected ${c2}'
}

fn test_root_is_the_square_root_of_alpha() {
	// root gives √α and its derivatives, and eval gives α: r² = α,
	// 2r·r′ = α′ and 2(r′² + r·r″) = α″, for every form, at every Tr.
	// Tr = 20 is past the Soave and PRSV zeros, where r is negative.
	for name, a in alphas() {
		for tr in [0.3, 0.7, 0.999, 1.001, 1.5, 4.0, 20.0] {
			alpha, d1, d2 := a.eval(tr)
			r, r1, r2 := a.root(tr)
			assert close(r * r, alpha), '${name} at Tr=${tr}: r² ${r * r} vs α ${alpha}'
			assert close(2.0 * r * r1, d1), '${name} at Tr=${tr}: 2rr′ vs α′'
			assert close(2.0 * (r1 * r1 + r * r2), d2), '${name} at Tr=${tr}: α″'
			// And root's own derivatives, by central differences.
			h := 1e-6 * tr
			up, r1_up, _ := a.root(tr + h)
			dn, r1_dn, _ := a.root(tr - h)
			assert close(r1, (up - dn) / (2.0 * h)), '${name} at Tr=${tr}: r′'
			assert close(r2, (r1_up - r1_dn) / (2.0 * h)), '${name} at Tr=${tr}: r″'
		}
	}
}

fn test_soave_root_passes_through_zero() {
	// With m = 0.5, k = 1 + m(1 − √Tr) is exactly zero at Tr = 9. √α = k
	// crosses zero there with a finite slope, and α = k² touches zero and
	// rises again, as the published form does.
	a := Alpha(thermo.Soave{
		m: 0.5
	})
	r, r1, r2 := a.root(9.0)
	assert r == 0.0
	assert r1 == -0.5 / 6.0 && r2 == 0.5 / 108.0
	alpha, d1, d2 := a.eval(9.0)
	assert alpha == 0.0 && d1 == 0.0 && d2 == 2.0 * r1 * r1
	below, _, _ := a.root(8.99)
	above, _, _ := a.root(9.01)
	assert below > 0.0 && above < 0.0
	alpha_above, _, _ := a.eval(9.01)
	assert alpha_above > 0.0
}
