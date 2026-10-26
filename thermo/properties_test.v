// Checks for the model-independent properties in properties.v. Every
// derivative is compared with finite differences of re-solved states, and
// the identities are checked on top. An external test: pub API only.
import math
import thermo { Component, Root }

// Smith, Van Ness & Abbott, Table B.1.
const methane = Component{
	name:  'methane'
	tc:    190.6
	pc:    45.99e5
	omega: 0.012
	mw:    0.016043
}

const butane = Component{
	name:  'n-butane'
	tc:    425.1
	pc:    37.96e5
	omega: 0.200
	mw:    0.058123
}

const co2 = Component{
	name:  'carbon dioxide'
	tc:    304.2
	pc:    73.83e5
	omega: 0.224
	mw:    0.044010
}

struct Point {
	t    f64
	p    f64
	root Root
}

struct Case {
	label  string
	eos    thermo.Cubic
	n      []f64
	points []Point
}

// Each case has a vapor point and a dense point (compressed liquid, or a
// single dense root near the mixture critical point), both away
// from phase boundaries so that small T and P steps stay on the same root.
fn cases() []Case {
	kij := [0.0, 0.02, 0.10, 0.02, 0.0, 0.13, 0.10, 0.13, 0.0]
	pure := [Point{400.0, 5e5, .vapor}, Point{300.0, 30e5, .liquid}]
	mixed := [Point{400.0, 20e5, .vapor}, Point{250.0, 60e5, .liquid}]
	return [
		Case{
			label:  'PR n-butane'
			eos:    thermo.new_cubic(.pr, [butane]) or { panic(err) }
			n:      [1.0]
			points: pure
		},
		Case{
			label:  'vdW n-butane'
			eos:    thermo.new_cubic(.vdw, [butane]) or { panic(err) }
			n:      [1.0]
			points: pure
		},
		Case{
			label:  'PR ternary, PRSV butane'
			eos:    thermo.new_cubic(.pr, [methane, butane, co2],
				kij:    kij
				alphas: [thermo.soave_pr(methane.omega), thermo.prsv(butane.omega, 0.04),
					thermo.soave_pr(co2.omega)]
			) or { panic(err) }
			n:      [0.2, 0.5, 0.3]
			points: mixed
		},
		Case{
			label:  'SRK ternary, 2 mol'
			eos:    thermo.new_cubic(.srk, [methane, butane, co2], kij: kij) or { panic(err) }
			n:      [1.2, 0.2, 0.6]
			points: mixed
		},
	]
}

fn solve(c Case, t f64, p f64, n []f64, root Root) thermo.State {
	return thermo.solve_tp(c.eos, t, p, n, root) or { panic('${c.label}: ${err}') }
}

fn total(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

fn max_abs(x []f64) f64 {
	mut m := 0.0
	for v in x {
		m = math.max(m, math.abs(v))
	}
	return m
}

// close judges an analytic value against a finite-difference estimate,
// relative to `scale` (the size of the quantity's family) or to the
// values themselves, whichever is larger.
fn close(analytic f64, fd f64, scale f64) bool {
	return math.abs(analytic - fd) <= 1e-6 * math.max(scale, math.max(math.abs(analytic),
		math.abs(fd)))
}

fn test_pressure_derivs_match_finite_differences() {
	for c in cases() {
		for pt in c.points {
			st := solve(c, pt.t, pt.p, c.n, pt.root)
			t, v := st.t, st.v
			d := thermo.pressure_derivs(c.eos, t, v, c.n)
			at := '${c.label} ${pt.root} at ${t} K'
			assert math.tolerance(d.p, pt.p, 1e-9), '${at}: P'
			hv := 1e-6 * v
			fd_v := (thermo.pressure(c.eos, t, v + hv, c.n) - thermo.pressure(c.eos, t, v - hv,
				c.n)) / (2.0 * hv)
			assert close(d.dp_dv, fd_v, 0.0), '${at}: dP/dV ${d.dp_dv} vs ${fd_v}'
			ht := 1e-4 * t
			fd_t := (thermo.pressure(c.eos, t + ht, v, c.n) - thermo.pressure(c.eos, t - ht, v,
				c.n)) / (2.0 * ht)
			assert close(d.dp_dt, fd_t, 0.0), '${at}: dP/dT ${d.dp_dt} vs ${fd_t}'
			hn := 1e-6 * total(c.n)
			for i in 0 .. c.n.len {
				mut up := c.n.clone()
				mut dn := c.n.clone()
				up[i] += hn
				dn[i] -= hn
				fd_n := (thermo.pressure(c.eos, t, v, up) - thermo.pressure(c.eos, t, v, dn)) / (2.0 * hn)
				assert close(d.dp_dn[i], fd_n, max_abs(d.dp_dn)), '${at}: dP/dn[${i}]'
			}
			// V is homogeneous of degree 1 in n at fixed T and P: Σ nᵢV̄ᵢ = V.
			mut sum_v := 0.0
			for i in 0 .. c.n.len {
				sum_v += c.n[i] * d.v_bar[i]
			}
			assert math.tolerance(sum_v, v, 1e-10), '${at}: Σ nᵢV̄ᵢ = ${sum_v}, V = ${v}'
		}
	}
}

fn test_ln_phi_derivs_match_finite_differences() {
	for c in cases() {
		nc := c.n.len
		for pt in c.points {
			st := solve(c, pt.t, pt.p, c.n, pt.root)
			d := thermo.ln_phi_derivs(c.eos, st)
			at := '${c.label} ${pt.root} at ${pt.t} K'

			// 1e-5·T, not 1e-4·T: near a mixture critical point (the SRK ternary
			// at 250 K, 60 bar) properties curve so sharply in T that a larger step
			// leaves a truncation error above 1e-6.
			ht := 1e-5 * pt.t
			tp := solve(c, pt.t + ht, pt.p, c.n, pt.root)
			tm := solve(c, pt.t - ht, pt.p, c.n, pt.root)
			hp := 1e-5 * pt.p
			pp := solve(c, pt.t, pt.p + hp, c.n, pt.root)
			pm := solve(c, pt.t, pt.p - hp, c.n, pt.root)
			for i in 0 .. nc {
				fd_t := (tp.ln_phi[i] - tm.ln_phi[i]) / (2.0 * ht)
				assert close(d.dt[i], fd_t, max_abs(d.dt)), '${at}: dlnφ/dT[${i}] ${d.dt[i]} vs ${fd_t}'
				fd_p := (pp.ln_phi[i] - pm.ln_phi[i]) / (2.0 * hp)
				assert close(d.dp[i], fd_p, max_abs(d.dp)), '${at}: dlnφ/dP[${i}] ${d.dp[i]} vs ${fd_p}'
			}

			// n·∂ln φᵢ/∂nⱼ at constant T and P: perturb nⱼ and re-solve.
			ntot := total(c.n)
			hn := 1e-5 * ntot
			scale := math.max(max_abs(d.dn), 1.0)
			for j in 0 .. nc {
				mut up := c.n.clone()
				mut dn := c.n.clone()
				up[j] += hn
				dn[j] -= hn
				su := solve(c, pt.t, pt.p, up, pt.root)
				sd := solve(c, pt.t, pt.p, dn, pt.root)
				for i in 0 .. nc {
					fd := ntot * (su.ln_phi[i] - sd.ln_phi[i]) / (2.0 * hn)
					assert close(d.dn[i * nc + j], fd, scale), '${at}: n·dlnφ[${i}]/dn[${j}] ${d.dn[i * nc + j]} vs ${fd}'
				}
			}
		}
	}
}

fn test_composition_derivatives_obey_gibbs_duhem_and_symmetry() {
	for c in cases() {
		nc := c.n.len
		for pt in c.points {
			d := thermo.ln_phi_derivs(c.eos, solve(c, pt.t, pt.p, c.n, pt.root))
			for j in 0 .. nc {
				mut s := 0.0
				mut size := 0.0
				for i in 0 .. nc {
					s += c.n[i] * d.dn[i * nc + j]
					size += math.abs(c.n[i] * d.dn[i * nc + j])
					assert math.abs(d.dn[i * nc + j] - d.dn[j * nc + i]) <= 1e-12 * math.max(size,
						1.0), '${c.label}: dn not symmetric'
				}
				assert math.abs(s) <= 1e-12 * math.max(size, 1.0), '${c.label}: Gibbs–Duhem column ${j} sums to ${s}'
			}
		}
	}
}

fn test_summed_ln_phi_derivatives_give_the_residual_properties() {
	// Σ nᵢ·∂ln φᵢ/∂T = −Hʳ/RT² and Σ nᵢ·∂ln φᵢ/∂P = (V − nRT/P)/RT.
	for c in cases() {
		for pt in c.points {
			st := solve(c, pt.t, pt.p, c.n, pt.root)
			d := thermo.ln_phi_derivs(c.eos, st)
			rp := thermo.residual_props(c.eos, st)
			rt := thermo.r_gas * pt.t
			mut sum_t := 0.0
			mut sum_p := 0.0
			for i in 0 .. c.n.len {
				sum_t += c.n[i] * d.dt[i]
				sum_p += c.n[i] * d.dp[i]
			}
			assert math.tolerance(sum_t, -rp.h / (rt * pt.t), 1e-10), c.label
			v_res := st.v - total(c.n) * rt / pt.p
			assert math.tolerance(sum_p, v_res / rt, 1e-10), c.label
		}
	}
}

fn test_residual_props_match_finite_differences() {
	for c in cases() {
		for pt in c.points {
			st := solve(c, pt.t, pt.p, c.n, pt.root)
			rp := thermo.residual_props(c.eos, st)
			at := '${c.label} ${pt.root} at ${pt.t} K'
			// 1e-5·T, not 1e-4·T: near a mixture critical point (the SRK ternary
			// at 250 K, 60 bar) properties curve so sharply in T that a larger step
			// leaves a truncation error above 1e-6.
			ht := 1e-5 * pt.t
			up := thermo.residual_props(c.eos, solve(c, pt.t + ht, pt.p, c.n, pt.root))
			dn := thermo.residual_props(c.eos, solve(c, pt.t - ht, pt.p, c.n, pt.root))
			// At constant P: Cpʳ = ∂Hʳ/∂T and Sʳ = −∂Gʳ/∂T.
			fd_cp := (up.h - dn.h) / (2.0 * ht)
			assert close(rp.cp, fd_cp, 0.0), '${at}: Cpʳ ${rp.cp} vs ${fd_cp}'
			fd_s := -(up.g - dn.g) / (2.0 * ht)
			assert close(rp.s, fd_s, 0.0), '${at}: Sʳ ${rp.s} vs ${fd_s}'
			// At constant T: ∂Gʳ/∂P = V − nRT/P.
			hp := 1e-5 * pt.p
			gp := thermo.residual_props(c.eos, solve(c, pt.t, pt.p + hp, c.n, pt.root)).g
			gm := thermo.residual_props(c.eos, solve(c, pt.t, pt.p - hp, c.n, pt.root)).g
			v_res := st.v - total(c.n) * thermo.r_gas * pt.t / pt.p
			assert close(v_res, (gp - gm) / (2.0 * hp), 0.0), '${at}: ∂Gʳ/∂P'
			// At constant V: Cvʳ = ∂Uʳ/∂T, with Uʳ = −RT²·F_T.
			ur := fn [c, st] (t f64) f64 {
				return -thermo.r_gas * t * t * c.eos.residual(t, st.v, c.n).f_t
			}
			fd_cv := (ur(pt.t + ht) - ur(pt.t - ht)) / (2.0 * ht)
			assert close(rp.cv, fd_cv, 0.0), '${at}: Cvʳ ${rp.cv} vs ${fd_cv}'
		}
	}
}

fn test_residual_gibbs_energy_is_the_fugacity_sum() {
	for c in cases() {
		for pt in c.points {
			st := solve(c, pt.t, pt.p, c.n, pt.root)
			rp := thermo.residual_props(c.eos, st)
			mut sum := 0.0
			for i in 0 .. c.n.len {
				sum += c.n[i] * st.ln_phi[i]
			}
			assert math.tolerance(rp.g, thermo.r_gas * pt.t * sum, 1e-10), c.label
		}
	}
}

fn test_residual_props_vanish_at_low_pressure() {
	for c in cases() {
		st := solve(c, 400.0, 1.0, c.n, .vapor)
		rp := thermo.residual_props(c.eos, st)
		nr := total(c.n) * thermo.r_gas
		for name, x in {
			'Hʳ/RT': rp.h / 400.0
			'Gʳ/RT': rp.g / 400.0
			'Uʳ/RT': rp.u / 400.0
			'Aʳ/RT': rp.a / 400.0
			'Sʳ':    rp.s
			'Cvʳ':   rp.cv
			'Cpʳ':   rp.cp
		} {
			assert math.abs(x) <= 1e-5 * nr, '${c.label}: ${name} = ${x / nr} nR at 1 Pa'
		}
	}
}
