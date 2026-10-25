// No `module thermo` line, so this is an external (black-box) test, like a
// Go file in `package thermo_test`: only `pub` names are reachable.
import math
import thermo { Component, CubicKind, Model }

// Smith, Van Ness & Abbott, Table B.1.
const butane = Component{
	name:  'n-butane'
	tc:    425.1
	pc:    37.96e5
	omega: 0.200
}

const methane = Component{
	name:  'methane'
	tc:    190.6
	pc:    45.99e5
	omega: 0.012
}

const kinds = [CubicKind.vdw, .rk, .srk, .pr]

// mixture is methane + n-butane with a non-zero kᵢⱼ and a PRSV α on one
// component, so every code path in the residual is exercised.
fn mixture() thermo.Cubic {
	return thermo.new_cubic(.pr, [methane, butane],
		kij:    [0.0, 0.02, 0.02, 0.0]
		alphas: [thermo.soave_pr(methane.omega), thermo.prsv(butane.omega, 0.04)]
	) or { panic(err) }
}

// ---- Agreement with the book ------------------------------------------------

fn test_matches_sva_example_3_9() {
	// Saturated n-butane at 350 K, 9.4573 bar. The book tabulates V in
	// cm³/mol to four figures (experiment: 2482 vapor, 115.0 liquid).
	t, p := 350.0, 9.4573e5
	book := {
		CubicKind.vdw: [191.0, 2667.0]
		.rk:           [133.3, 2555.0]
		.srk:          [127.8, 2520.0]
		.pr:           [112.6, 2486.0]
	}
	for kind in kinds {
		eos := thermo.new_cubic(kind, [butane])!
		want := book[kind] or { panic('no book value for ${kind}') }
		liq := thermo.solve_tp(eos, t, p, [1.0], .liquid)!.molar_volume() * 1e6
		vap := thermo.solve_tp(eos, t, p, [1.0], .vapor)!.molar_volume() * 1e6
		assert math.abs(liq - want[0]) <= 0.1, '${eos.name()} liquid: ${liq}'
		assert math.abs(vap - want[1]) <= 1.0, '${eos.name()} vapor: ${vap}'
	}
}

fn test_pr_saturation_pressure_matches_experiment() {
	// SVA Example 3.9 gives Psat(350 K) = 9.4573 bar for n-butane.
	eos := thermo.new_cubic(.pr, [butane])!
	psat := thermo.saturation_pressure(eos, 350.0)!
	assert math.abs(psat / 9.4573e5 - 1.0) < 0.02, 'PR Psat = ${psat / 1e5} bar'
	liq := thermo.solve_tp(eos, 350.0, psat, [1.0], .liquid)!
	vap := thermo.solve_tp(eos, 350.0, psat, [1.0], .vapor)!
	assert math.tolerance(liq.ln_phi[0], vap.ln_phi[0], 1e-8)
}

// ---- The generic solver against the cubic's own polynomial ------------------

fn test_newton_volume_solver_matches_polynomial_roots() {
	mix := mixture()
	mut cases := []TestCase{}
	for kind in kinds {
		cases << TestCase{thermo.new_cubic(kind, [butane])!, [1.0]}
	}
	cases << TestCase{mix, [0.3, 0.7]}
	for c in cases {
		for t in [250.0, 350.0, 420.0, 600.0] {
			for p in [1e3, 1e5, 9e5, 4e6, 2e7] {
				roots := c.eos.z_roots(t, p, c.n)
				liq := thermo.solve_tp(c.eos, t, p, c.n, .liquid)!
				vap := thermo.solve_tp(c.eos, t, p, c.n, .vapor)!
				assert math.tolerance(liq.z, roots.first(), 1e-9), '${c.eos.name()} T=${t} P=${p} liquid'
				assert math.tolerance(vap.z, roots.last(), 1e-9), '${c.eos.name()} T=${t} P=${p} vapor'
			}
		}
	}
}

struct TestCase {
	eos thermo.Cubic
	n   []f64
}

// ---- Thermodynamic identities -----------------------------------------------

fn test_derivatives_match_finite_differences() {
	eos := mixture()
	n := [0.3, 0.7]
	for t in [250.0, 400.0] {
		for v in [1.2e-4, 3e-3] {
			r := eos.residual(t, v, n)
			hv := 1e-6 * v
			assert close(r.f_v, (f_at(eos, t, v + hv, n) - f_at(eos, t, v - hv, n)) / (2.0 * hv))
			assert close(r.f_vv, (eos.volume_derivs(t, v + hv, n).f_v - eos.volume_derivs(t, v - hv,
				n).f_v) / (2.0 * hv))
			ht := 1e-4
			assert close(r.f_t, (f_at(eos, t + ht, v, n) - f_at(eos, t - ht, v, n)) / (2.0 * ht))
			for i in 0 .. n.len {
				mut up := n.clone()
				mut dn := n.clone()
				up[i] += 1e-6
				dn[i] -= 1e-6
				assert close(r.f_n[i], (f_at(eos, t, v, up) - f_at(eos, t, v, dn)) / 2e-6), 'F_n[${i}]'
			}
		}
	}
}

fn f_at(m Model, t f64, v f64, n []f64) f64 {
	return m.volume_derivs(t, v, n).f
}

fn close(a f64, b f64) bool {
	return math.abs(a - b) <= 1e-6 * math.max(math.abs(a), math.abs(b)) + 1e-12
}

fn test_f_is_homogeneous_in_volume_and_moles() {
	// F(T, λV, λn) = λF, so by Euler's theorem F = V·F_V + Σ nᵢ·F_nᵢ.
	eos := mixture()
	n := [0.3, 0.7]
	for v in [1.2e-4, 3e-3] {
		r := eos.residual(350.0, v, n)
		euler := v * r.f_v + n[0] * r.f_n[0] + n[1] * r.f_n[1]
		assert math.tolerance(r.f, euler, 1e-10)
	}
}

fn test_fugacities_add_up_to_residual_gibbs_energy() {
	// Σ nᵢ ln φᵢ = Gʳ/RT = F + n(Z − 1) − n ln Z
	eos := mixture()
	n := [0.3, 0.7]
	for root in [thermo.Root.liquid, .vapor] {
		s := thermo.solve_tp(eos, 350.0, 2e6, n, root)!
		r := eos.residual(350.0, s.v, n)
		gr := r.f + (s.z - 1.0) - math.log(s.z)
		assert math.tolerance(n[0] * s.ln_phi[0] + n[1] * s.ln_phi[1], gr, 1e-9)
	}
}

fn test_results_are_intensive() {
	// Doubling every mole number doubles V and leaves Z and ln φ alone.
	eos := mixture()
	a := thermo.solve_tp(eos, 350.0, 2e6, [0.3, 0.7], .vapor)!
	b := thermo.solve_tp(eos, 350.0, 2e6, [0.6, 1.4], .vapor)!
	assert math.tolerance(b.v, 2.0 * a.v, 1e-10)
	assert math.tolerance(b.z, a.z, 1e-10)
	assert math.tolerance(b.ln_phi[1], a.ln_phi[1], 1e-10)
}

fn test_a_component_mixed_with_itself_is_the_pure_component() {
	pure := thermo.new_cubic(.pr, [butane])!
	twice := thermo.new_cubic(.pr, [butane, butane])!
	sp := thermo.solve_tp(pure, 350.0, 5e5, [1.0], .vapor)!
	sm := thermo.solve_tp(twice, 350.0, 5e5, [0.4, 0.6], .vapor)!
	assert math.tolerance(sm.z, sp.z, 1e-12)
	assert math.tolerance(sm.ln_phi[0], sp.ln_phi[0], 1e-10)
	assert math.tolerance(sm.ln_phi[1], sp.ln_phi[0], 1e-10)
}

fn test_ideal_gas_limit() {
	for kind in kinds {
		eos := thermo.new_cubic(kind, [butane])!
		s := thermo.solve_tp(eos, 500.0, 1.0, [1.0], .vapor)!
		assert math.tolerance(s.z, 1.0, 1e-6), eos.name()
		assert math.tolerance(s.ln_phi[0], 0.0, 1e-6), eos.name()
	}
}

// ---- Saturation and phase choice --------------------------------------------

fn test_saturation_pressure_converges_over_the_range() {
	for kind in kinds {
		eos := thermo.new_cubic(kind, [butane])!
		for tr in [0.5, 0.7, 0.9, 0.98] {
			thermo.saturation_pressure(eos, tr * butane.tc) or {
				assert false, '${eos.name()} at Tr = ${tr}: ${err}'
			}
		}
	}
}

fn test_stable_root_switches_across_psat() {
	eos := thermo.new_cubic(.pr, [butane])!
	psat := thermo.saturation_pressure(eos, 350.0)!
	below := thermo.solve_tp(eos, 350.0, 0.8 * psat, [1.0], .stable)!
	above := thermo.solve_tp(eos, 350.0, 1.2 * psat, [1.0], .stable)!
	assert below.z > 0.5 // vapor
	assert above.z < 0.1 // liquid
}

// ---- PRSV --------------------------------------------------------------------

fn test_prsv_with_zero_kappa1_is_close_to_pr() {
	// κ₀(ω) is a refit of PR's m(ω), so with κ₁ = 0 the two should agree
	// closely but not exactly.
	pr := thermo.new_cubic(.pr, [butane])!
	prsv := thermo.new_cubic(.pr, [butane], alphas: [thermo.prsv(butane.omega, 0.0)])!
	p1 := thermo.saturation_pressure(pr, 350.0)!
	p2 := thermo.saturation_pressure(prsv, 350.0)!
	assert p1 != p2
	assert math.abs(p2 / p1 - 1.0) < 0.01, 'PR ${p1}, PRSV ${p2}'
}

fn test_prsv_alpha_derivative_is_smooth_through_tr_0_7() {
	a := thermo.prsv(0.2, 0.05)
	for tr in [0.3, 0.6999, 0.7, 0.7001, 1.5] {
		_, d := a.eval(tr)
		h := 1e-6
		up, _ := a.eval(tr + h)
		dn, _ := a.eval(tr - h)
		assert close(d, (up - dn) / (2.0 * h)), 'Tr = ${tr}'
	}
}

// ---- Errors ------------------------------------------------------------------

fn test_bad_input_is_a_typed_error() {
	broken := Component{
		name:  'broken'
		tc:    0.0
		pc:    1e6
		omega: 0.0
	}
	thermo.new_cubic(.pr, [broken]) or {
		// `is` checks the concrete error type: Go's errors.As.
		assert err is thermo.InputError
		return
	}
	assert false, 'expected an InputError'
}

fn test_bad_options_are_rejected() {
	if _ := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.1]) {
		assert false, 'kij of the wrong size'
	}
	if _ := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.1, 0.2, 0.0]) {
		assert false, 'asymmetric kij'
	}
	if _ := thermo.new_cubic(.pr, [methane, butane], alphas: [thermo.soave_pr(0.0)]) {
		assert false, 'too few alphas'
	}
	if _ := thermo.solve_tp(mixture(), 300.0, 1e5, [1.0], .vapor) {
		assert false, 'composition of the wrong length'
	}
}

fn test_no_saturation_pressure_above_tc_or_for_mixtures() {
	pure := thermo.new_cubic(.pr, [butane])!
	if _ := thermo.saturation_pressure(pure, 450.0) {
		assert false, 'expected an error above Tc'
	}
	if _ := thermo.saturation_pressure(mixture(), 300.0) {
		assert false, 'expected an error for a mixture'
	}
}
