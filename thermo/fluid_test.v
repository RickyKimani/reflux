// Checks for Fluid.props: total properties = ideal gas + residual. An
// external test, through the pub API only.
import math
import thermo { Component, Root }

// Critical constants: Smith, Van Ness & Abbott, Table B.1. Ideal-gas Cp:
// the same book's Table C.1 (Cp/R = A + BT + CT² + DT⁻²). ΔH°f: Table C.4.
// S° values are illustrative: the tests check identities, not data, and
// real values arrive with M2.
const butane = Component{
	name:  'n-butane'
	tc:    425.1
	pc:    37.96e5
	omega: 0.200
	mw:    0.058123
}

const methane = Component{
	name:  'methane'
	tc:    190.6
	pc:    45.99e5
	omega: 0.012
	mw:    0.016043
}

const co2 = Component{
	name:  'carbon dioxide'
	tc:    304.2
	pc:    73.83e5
	omega: 0.224
	mw:    0.044010
}

fn four_constant(a f64, b f64, c f64, d f64) thermo.CpForm {
	r := thermo.r_gas
	return thermo.Polynomial{
		terms: [
			thermo.PowerTerm{r * a, 0},
			thermo.PowerTerm{r * b, 1},
			thermo.PowerTerm{r * c, 2},
			thermo.PowerTerm{r * d, -2},
		]
	}
}

fn butane_ig() thermo.IdealGas {
	return thermo.IdealGas{
		cp:  four_constant(1.935, 36.915e-3, -11.402e-6, 0.0)
		dhf: -125.79e3
		s0:  310.0
	}
}

fn methane_ig() thermo.IdealGas {
	return thermo.IdealGas{
		cp:  four_constant(1.702, 9.081e-3, -2.164e-6, 0.0)
		dhf: -74.52e3
		s0:  186.3
	}
}

fn co2_ig() thermo.IdealGas {
	return thermo.IdealGas{
		cp:  four_constant(5.457, 1.045e-3, 0.0, -1.157e5)
		dhf: -393.51e3
		s0:  213.8
	}
}

struct Case {
	label string
	fluid thermo.Fluid
	n     []f64
	t     f64
	p     f64
	root  Root
}

fn pure_butane() thermo.Fluid {
	eos := thermo.new_cubic(.pr, [butane]) or { panic(err) }
	return thermo.new_fluid(eos, [butane_ig()]) or { panic(err) }
}

fn ternary() thermo.Fluid {
	eos := thermo.new_cubic(.pr, [methane, butane, co2],
		kij: [0.0, 0.02, 0.10, 0.02, 0.0, 0.13, 0.10, 0.13, 0.0]
	) or { panic(err) }
	return thermo.new_fluid(eos, [methane_ig(), butane_ig(), co2_ig()]) or { panic(err) }
}

fn cases() []Case {
	return [
		Case{'PR n-butane vapor', pure_butane(), [1.0], 400.0, 5e5, .vapor},
		Case{'PR n-butane liquid', pure_butane(), [1.0], 300.0, 30e5, .liquid},
		Case{'PR ternary vapor', ternary(), [0.2, 0.5, 0.3], 400.0, 20e5, .vapor},
	]
}

fn solve(c Case, t f64, p f64) thermo.State {
	return thermo.solve_tp(c.fluid.model, t, p, c.n, c.root) or { panic('${c.label}: ${err}') }
}

fn props_at(c Case, t f64, p f64) thermo.Props {
	return c.fluid.props(solve(c, t, p))
}

// state_tv solves at the temperature t and the volume v, by finding the
// pressure the model gives there and solving back from (t, P).
fn state_tv(c Case, t f64, v f64) thermo.State {
	p := thermo.pressure(c.fluid.model, t, v, c.n)
	st := solve(c, t, p)
	assert math.tolerance(st.v, v, 1e-9), '${c.label}: state_tv landed on another root'
	return st
}

fn total(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

fn test_reference_state() {
	// At 298.15 K and 1 bar the ideal-gas part is exactly the formation
	// data, plus −R Σ nᵢ ln xᵢ in S for a mixture.
	for c in [cases()[0], cases()[2]] {
		st := solve(c, thermo.t_ref, thermo.p_ref)
		pr := c.fluid.props(st)
		rp := thermo.residual_props(c.fluid.model, st)
		mut h0 := 0.0
		mut s0 := 0.0
		for i, ig in c.fluid.ideal {
			x := c.n[i] / total(c.n)
			h0 += c.n[i] * (ig.dhf or { 0.0 })
			s0 += c.n[i] * ((ig.s0 or { 0.0 }) - thermo.r_gas * math.log(x))
		}
		assert math.tolerance(pr.h - rp.h, h0, 1e-12), c.label
		assert math.tolerance(pr.s - rp.s, s0, 1e-12), c.label
		assert pr.absolute
	}
}

fn test_missing_formation_data_shifts_the_zero_only() {
	c := cases()[0]
	no_data := thermo.new_fluid(c.fluid.model, [
		thermo.IdealGas{
			cp: butane_ig().cp
		},
	])!
	st1 := solve(c, 400.0, 5e5)
	st2 := solve(c, 450.0, 8e5)
	with1, with2 := c.fluid.props(st1), c.fluid.props(st2)
	wo1, wo2 := no_data.props(st1), no_data.props(st2)
	assert !wo1.absolute
	assert math.tolerance(with2.h - with1.h, wo2.h - wo1.h, 1e-10)
	assert math.tolerance(with2.s - with1.s, wo2.s - wo1.s, 1e-10)
	assert math.tolerance(with1.cp, wo1.cp, 1e-14)
	assert math.tolerance(with1.w, wo1.w, 1e-14)
}

fn test_h_and_s_integrate_cp_along_an_isobar() {
	// ΔH = ∫Cp dT and ΔS = ∫Cp/T dT at constant P, for the real fluid.
	for c in [cases()[0], cases()[2]] {
		t1, t2 := 400.0, 500.0
		steps := 200
		h := (t2 - t1) / f64(steps)
		mut q_h := 0.0
		mut q_s := 0.0
		for i in 0 .. steps + 1 {
			t := t1 + f64(i) * h
			w := if i == 0 || i == steps {
				1.0
			} else if i % 2 == 1 {
				4.0
			} else {
				2.0
			}
			cp := props_at(c, t, c.p).cp
			q_h += w * cp
			q_s += w * cp / t
		}
		q_h *= h / 3.0
		q_s *= h / 3.0
		p1, p2 := props_at(c, t1, c.p), props_at(c, t2, c.p)
		assert math.tolerance(p2.h - p1.h, q_h, 1e-9), '${c.label}: ΔH ${p2.h - p1.h} vs ∫Cp ${q_h}'
		assert math.tolerance(p2.s - p1.s, q_s, 1e-9), '${c.label}: ΔS ${p2.s - p1.s} vs ∫Cp/T ${q_s}'
	}
}

fn test_cp_and_cv_match_finite_differences() {
	for c in cases() {
		st := solve(c, c.t, c.p)
		pr := c.fluid.props(st)
		ht := 1e-5 * c.t
		fd_cp := (props_at(c, c.t + ht, c.p).h - props_at(c, c.t - ht, c.p).h) / (2.0 * ht)
		assert math.tolerance(pr.cp, fd_cp, 1e-6), '${c.label}: Cp ${pr.cp} vs ${fd_cp}'
		up := c.fluid.props(state_tv(c, c.t + ht, st.v))
		dn := c.fluid.props(state_tv(c, c.t - ht, st.v))
		fd_cv := (up.u - dn.u) / (2.0 * ht)
		assert math.tolerance(pr.cv, fd_cv, 1e-6), '${c.label}: Cv ${pr.cv} vs ${fd_cv}'
	}
}

fn test_speed_of_sound_matches_an_isentropic_finite_difference() {
	// w² = (∂P/∂ρ)_S. Change V by ±δ, find the T that keeps S fixed (Newton
	// with (∂S/∂T)_V = Cv/T), and difference P against the mass density.
	for c in cases() {
		st := solve(c, c.t, c.p)
		pr := c.fluid.props(st)
		mut mass := 0.0
		for i, comp in c.fluid.model.components() {
			mass += c.n[i] * comp.mw
		}
		mut p_at := [0.0, 0.0]
		mut rho := [0.0, 0.0]
		for k, sign in [1.0, -1.0] {
			v := st.v * (1.0 + sign * 1e-5)
			mut t := c.t
			for _ in 0 .. 20 {
				q := c.fluid.props(state_tv(c, t, v))
				dt := (q.s - pr.s) * t / q.cv
				t -= dt
				if math.abs(dt) < 1e-12 * t {
					break
				}
			}
			p_at[k] = thermo.pressure(c.fluid.model, t, v, c.n)
			rho[k] = mass / v
		}
		w2_fd := (p_at[0] - p_at[1]) / (rho[0] - rho[1])
		assert math.tolerance(pr.w * pr.w, w2_fd, 1e-6), '${c.label}: w² ${pr.w * pr.w} vs ${w2_fd}'
	}
}

fn test_joule_thomson_matches_an_isenthalpic_finite_difference() {
	// μ_JT = (∂T/∂P)_H. Change P by ±δ and find the T that keeps H fixed
	// (Newton with (∂H/∂T)_P = Cp).
	for c in cases() {
		pr := props_at(c, c.t, c.p)
		mut t_at := [0.0, 0.0]
		for k, sign in [1.0, -1.0] {
			p := c.p * (1.0 + sign * 1e-5)
			mut t := c.t
			for _ in 0 .. 20 {
				q := props_at(c, t, p)
				dt := (q.h - pr.h) / q.cp
				t -= dt
				if math.abs(dt) < 1e-12 * t {
					break
				}
			}
			t_at[k] = t
		}
		mu_fd := (t_at[0] - t_at[1]) / (2e-5 * c.p)
		assert math.tolerance(pr.mu_jt, mu_fd, 1e-6), '${c.label}: μ_JT ${pr.mu_jt} vs ${mu_fd}'
	}
}

fn test_ideal_gas_limits() {
	// At 1 Pa: Cp − Cv = nR and w² = (Cp/Cv)·RT/M̄.
	for c in [cases()[0], cases()[2]] {
		pr := props_at(c, 400.0, 1.0)
		ntot := total(c.n)
		mut mass := 0.0
		for i, comp in c.fluid.model.components() {
			mass += c.n[i] * comp.mw
		}
		assert math.tolerance(pr.cp - pr.cv, ntot * thermo.r_gas, 1e-6), c.label
		w2 := pr.cp / pr.cv * ntot * thermo.r_gas * 400.0 / mass
		assert math.tolerance(pr.w * pr.w, w2, 1e-6), c.label
	}
}

fn test_a_component_mixed_with_itself_only_adds_mixing_entropy() {
	pure := pure_butane()
	eos2 := thermo.new_cubic(.pr, [butane, butane])!
	twice := thermo.new_fluid(eos2, [butane_ig(), butane_ig()])!
	n := [0.4, 0.6]
	p1 := pure.props(thermo.solve_tp(pure.model, 400.0, 5e5, [1.0], .vapor)!)
	p2 := twice.props(thermo.solve_tp(twice.model, 400.0, 5e5, n, .vapor)!)
	s_mix := -thermo.r_gas * (0.4 * math.log(0.4) + 0.6 * math.log(0.6))
	assert math.tolerance(p2.s - p1.s, s_mix, 1e-10)
	assert math.tolerance(p2.h, p1.h, 1e-12)
	assert math.tolerance(p2.cp, p1.cp, 1e-12)
	assert math.tolerance(p2.w, p1.w, 1e-12)
}

fn test_new_fluid_checks_the_data_length() {
	eos := thermo.new_cubic(.pr, [methane, butane])!
	if _ := thermo.new_fluid(eos, [butane_ig()]) {
		assert false, 'expected an error'
	}
}
