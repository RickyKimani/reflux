module numeric

import math

// A closure wrapped here is never freed, which a test can afford (see
// functions.v).
// F wraps a function value as a Function.
struct F {
	f fn (f64) f64 @[required]
}

fn (w F) at(x f64) f64 {
	return w.f(x)
}

fn test_gk15_is_exact_for_polynomials_up_to_degree_22() {
	// One panel, no refinement: ∫₀¹ xᵏ dx = 1/(k+1) for k ≤ 22.
	for k in 0 .. 23 {
		got := gk15(F{fn [k] (x f64) f64 {
			return math.pow(x, f64(k))
		}}, 0.0, 1.0).value
		assert math.abs(got - 1.0 / f64(k + 1)) <= 2e-16, 'degree ${k}: ${got}'
	}
}

fn test_the_weights_are_consistent() {
	// Both rules integrate 1 exactly: the weights sum to 2 on [−1, 1].
	mut sk := wgk[7]
	mut sg := wg[3]
	for j in 0 .. 7 {
		sk += 2.0 * wgk[j]
		if j % 2 == 1 {
			sg += 2.0 * wg[j / 2]
		}
	}
	assert math.abs(sk - 2.0) <= 4e-16
	assert math.abs(sg - 2.0) <= 4e-16
}

fn test_smooth_integrands_reach_the_requested_accuracy() {
	cases := [
		// ∫₀³ eˣ dx
		[0.0, 3.0, math.exp(3.0) - 1.0],
		// ∫₁¹⁰ dx/x
		[1.0, 10.0, math.log(10.0)],
		// ∫₀^π sin x dx
		[0.0, math.pi, 2.0],
	]
	fns := [fn (x f64) f64 {
		return math.exp(x)
	}, fn (x f64) f64 {
		return 1.0 / x
	}, fn (x f64) f64 {
		return math.sin(x)
	}]
	for i, c in cases {
		got := integrate(F{fns[i]}, c[0], c[1], 1e-13)!
		assert math.tolerance(got, c[2], 1e-13), 'case ${i}: ${got} vs ${c[2]}'
	}
}

fn test_a_sharp_peak_forces_refinement() {
	// A narrow Lorentzian: one panel can't resolve it, refinement must.
	w := 1e-3
	got := integrate(F{fn [w] (x f64) f64 {
		return w / (x * x + w * w)
	}}, -1.0, 1.0, 1e-12)!
	want := 2.0 * math.atan(1.0 / w)
	assert math.tolerance(got, want, 1e-11), '${got} vs ${want}'
}

fn test_reversed_and_empty_ranges() {
	f := fn (x f64) f64 {
		return x * x
	}
	assert integrate(F{f}, 2.0, 2.0, 1e-12)! == 0.0
	assert math.tolerance(integrate(F{f}, 1.0, 0.0, 1e-12)!, -1.0 / 3.0, 1e-15)
}

struct Counter {
mut:
	n int
}

fn test_an_unreachable_rtol_gets_roundoff_accuracy_quickly() {
	// 1e-17 is below what f64 can deliver. integrate should settle for
	// roundoff-limited accuracy in one panel, not split forever.
	//
	// V closures capture by value, so the count goes through a reference:
	// `fn [mut calls]` on an int would bump the closure's own copy.
	mut calls := &Counter{}
	got := integrate(F{fn [mut calls] (x f64) f64 {
		calls.n++
		return math.exp(x)
	}}, 0.0, 3.0, 1e-17)!
	assert math.tolerance(got, math.exp(3.0) - 1.0, 1e-14)
	assert calls.n == 15, '${calls.n} evaluations'
}

fn test_an_endpoint_singularity_converges() {
	// ∫₀¹ x^0.4 dx = 1/1.4: the shape of DIPPR 106 near Tc, where the
	// derivative is infinite. Only the panel at the singularity is split.
	got := integrate(F{fn (x f64) f64 {
		return math.pow(1.0 - x, 0.4)
	}}, 0.0, 1.0, 1e-13)!
	assert math.tolerance(got, 1.0 / 1.4, 1e-13), '${got}'
}

fn test_a_divergent_integral_is_an_error() {
	// ∫₀¹ dx/x diverges. The panel at 0 keeps being halved; once it's down
	// to denormal widths, a node rounds onto 0 itself and f(0) = ∞ ends it.
	if v := integrate(F{fn (x f64) f64 {
		return 1.0 / x
	}}, 0.0, 1.0, 1e-10)
	{
		assert false, 'a divergent integral returned ${v}'
	}
}

fn test_noise_exhausts_the_panel_budget() {
	// A hash of x looks like noise at every scale, so no panel's estimate
	// ever settles and the budget runs out.
	if v := integrate(F{fn (x f64) f64 {
		s := math.sin(x * 12.9898e3) * 43758.5453
		return s - math.floor(s)
	}}, 0.0, 1.0, 1e-10)
	{
		assert false, 'noise returned ${v}'
	} else {
		assert err.msg().contains('with ${max_panels} panels')
	}
}

fn test_non_finite_integrands_are_errors() {
	if _ := integrate(F{fn (x f64) f64 {
		return 1.0 / (x - 0.5)
	}}, 0.0, 1.0, 1e-12)
	{
		assert false, 'a pole should not integrate quietly'
	}
	if _ := integrate(F{fn (x f64) f64 {
		return math.log(x)
	}}, -1.0, 1.0, 1e-12)
	{
		assert false, 'NaN in the integrand should be an error'
	}
}
