// Checks for the volume-shift wrapper. The shift must move volumes and
// ln φ by exactly the predicted amounts and leave phase equilibrium alone,
// and its derivatives must be as consistent as any other model's.
import math
import thermo { Component, Model, Root }

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

struct Case {
	label   string
	inner   thermo.Cubic
	shifted thermo.Translated
	n       []f64
	c       []f64
}

// The shifts are illustrative, a few cm³/mol, which is the size Péneloux
// corrections usually have. One is negative, to check the sign handling.
fn cases() []Case {
	pr := thermo.new_cubic(.pr, [butane]) or { panic(err) }
	mix := thermo.new_cubic(.pr, [methane, butane, co2],
		kij: [0.0, 0.02, 0.10, 0.02, 0.0, 0.13, 0.10, 0.13, 0.0]
	) or { panic(err) }
	c_mix := [1e-6, 5e-6, -2e-6]
	return [
		Case{'PR n-butane', pr, thermo.translate(pr, [5e-6]) or { panic(err) }, [1.0], [
			5e-6,
		]},
		Case{'PR ternary', mix, thermo.translate(mix, c_mix) or { panic(err) }, [0.2, 0.5, 0.3], c_mix},
	]
}

struct Point {
	t    f64
	p    f64
	root Root
}

const points = [Point{400.0, 20e5, .vapor}, Point{250.0, 60e5, .liquid}, Point{300.0, 30e5, .liquid}]

fn solve(m Model, pt Point, n []f64) thermo.State {
	return thermo.solve_tp(m, pt.t, pt.p, n, pt.root) or { panic(err) }
}

fn dot(a []f64, b []f64) f64 {
	mut s := 0.0
	for i, x in a {
		s += x * b[i]
	}
	return s
}

fn test_volume_moves_by_exactly_c() {
	for c in cases() {
		big_c := dot(c.n, c.c)
		for pt in points {
			a := solve(c.inner, pt, c.n)
			b := solve(c.shifted, pt, c.n)
			assert math.tolerance(b.v, a.v - big_c, 1e-10), '${c.label} ${pt.root} at ${pt.t} K'
		}
	}
}

fn test_ln_phi_moves_by_exactly_c_p_over_rt() {
	for c in cases() {
		for pt in points {
			a := solve(c.inner, pt, c.n)
			b := solve(c.shifted, pt, c.n)
			rt := thermo.r_gas * pt.t
			for i in 0 .. c.n.len {
				want := a.ln_phi[i] - c.c[i] * pt.p / rt
				assert math.abs(b.ln_phi[i] - want) <= 1e-10 * math.max(1.0, math.abs(want)), '${c.label}: ln φ[${i}]'
			}
		}
	}
}

fn test_saturation_pressure_is_unchanged() {
	c := cases()[0]
	for t in [250.0, 300.0, 350.0, 400.0] {
		p1 := thermo.saturation_pressure(c.inner, t)!
		p2 := thermo.saturation_pressure(c.shifted, t)!
		assert math.tolerance(p1, p2, 1e-9), 'Psat at ${t} K: ${p1} vs ${p2}'
	}
}

fn test_residual_enthalpy_moves_by_p_c_and_nothing_else_does() {
	// With constant c, ln φ′ᵢ = ln φᵢ − cᵢP/RT gives G′ʳ = Gʳ − PC. Since PC
	// doesn't depend on T, Sʳ is unchanged and H′ʳ = Hʳ − PC, so Cpʳ is
	// unchanged too.
	for c in cases() {
		big_c := dot(c.n, c.c)
		for pt in points {
			a := thermo.residual_props(c.inner, solve(c.inner, pt, c.n))
			b := thermo.residual_props(c.shifted, solve(c.shifted, pt, c.n))
			scale := math.max(math.abs(a.h), thermo.r_gas * pt.t)
			assert math.abs(b.g - (a.g - pt.p * big_c)) <= 1e-9 * scale, '${c.label}: Gʳ'
			assert math.abs(b.h - (a.h - pt.p * big_c)) <= 1e-9 * scale, '${c.label}: Hʳ'
			assert math.abs(b.s - a.s) <= 1e-9 * math.max(math.abs(a.s), thermo.r_gas), '${c.label}: Sʳ'
			assert math.tolerance(b.cp, a.cp, 1e-8), '${c.label}: Cpʳ ${b.cp} vs ${a.cp}'
		}
	}
}

fn test_two_shifts_add() {
	c := cases()[1]
	once := thermo.translate(c.inner, [2e-6, 6e-6, 1e-6])!
	twice := thermo.translate(thermo.translate(c.inner, [1e-6, 4e-6, 3e-6])!, [1e-6, 2e-6, -2e-6])!
	for pt in points {
		a := solve(once, pt, c.n)
		b := solve(twice, pt, c.n)
		assert math.tolerance(a.v, b.v, 1e-10)
		for i in 0 .. c.n.len {
			assert math.abs(a.ln_phi[i] - b.ln_phi[i]) <= 1e-10 * math.max(1.0, math.abs(a.ln_phi[i]))
		}
		r1 := once.residual2(pt.t, a.v, c.n)
		r2 := twice.residual2(pt.t, a.v, c.n)
		for k in 0 .. r1.f_nn.len {
			assert math.abs(r1.f_nn[k] - r2.f_nn[k]) <= 1e-9 * math.max(1.0, math.abs(r1.f_nn[k]))
		}
	}
}

fn close(analytic f64, fd f64, scale f64) bool {
	return math.abs(analytic - fd) <= 1e-6 * math.max(scale, math.max(math.abs(analytic),
		math.abs(fd)))
}

fn max_abs(x []f64) f64 {
	mut m := 0.0
	for v in x {
		m = math.max(m, math.abs(v))
	}
	return m
}

fn test_derivatives_match_finite_differences() {
	for c in cases() {
		m := c.shifted
		nc := c.n.len
		for pt in points {
			v := solve(m, pt, c.n).v
			t := pt.t
			r := m.residual2(t, v, c.n)
			hv := 1e-6 * v
			ht := 1e-5 * t
			at := '${c.label} ${pt.root} at ${t} K'
			vp, vm := m.residual(t, v + hv, c.n), m.residual(t, v - hv, c.n)
			tp, tm := m.residual(t + ht, v, c.n), m.residual(t - ht, v, c.n)
			// First derivatives against F itself.
			fv := (m.volume_derivs(t, v + hv, c.n).f - m.volume_derivs(t, v - hv, c.n).f) / (2.0 * hv)
			assert close(r.f_v, fv, 0.0), '${at}: F_V'
			assert close(r.f_vv, (vp.f_v - vm.f_v) / (2.0 * hv), 0.0), '${at}: F_VV'
			assert close(r.f_tv, (tp.f_v - tm.f_v) / (2.0 * ht), 0.0), '${at}: F_TV'
			assert close(r.f_tt, (tp.f_t - tm.f_t) / (2.0 * ht), 0.0), '${at}: F_TT'
			hn := 1e-6
			for i in 0 .. nc {
				assert close(r.f_vn[i], (vp.f_n[i] - vm.f_n[i]) / (2.0 * hv), max_abs(r.f_vn)), '${at}: F_Vn[${i}]'
				assert close(r.f_tn[i], (tp.f_n[i] - tm.f_n[i]) / (2.0 * ht), max_abs(r.f_tn)), '${at}: F_Tn[${i}]'
				mut up := c.n.clone()
				mut dn := c.n.clone()
				up[i] += hn
				dn[i] -= hn
				fp := m.volume_derivs(t, v, up).f
				fm := m.volume_derivs(t, v, dn).f
				assert close(r.f_n[i], (fp - fm) / (2.0 * hn), max_abs(r.f_n)), '${at}: F_n[${i}]'
				np, nm := m.residual(t, v, up), m.residual(t, v, dn)
				for j in 0 .. nc {
					fd := (np.f_n[j] - nm.f_n[j]) / (2.0 * hn)
					assert close(r.f_nn[j * nc + i], fd, max_abs(r.f_nn)), '${at}: F_nn[${j}][${i}]'
				}
			}
		}
	}
}

fn test_identities_still_hold() {
	// The same identities as for any Model: F's homogeneity, Gibbs–Duhem,
	// and the summed ln φ derivatives.
	for c in cases() {
		m := c.shifted
		nc := c.n.len
		for pt in points {
			st := solve(m, pt, c.n)
			r := m.residual2(pt.t, st.v, c.n)
			mut euler := st.v * r.f_v
			for i in 0 .. nc {
				euler += c.n[i] * r.f_n[i]
			}
			assert math.abs(euler - r.f) <= 1e-10 * math.max(1.0, math.abs(r.f)), '${c.label}: F = V·F_V + Σ n·F_n'
			d := thermo.ln_phi_derivs(m, st)
			rp := thermo.residual_props(m, st)
			mut sum_t := 0.0
			for i in 0 .. nc {
				sum_t += c.n[i] * d.dt[i]
				mut gd := 0.0
				for k in 0 .. nc {
					gd += c.n[k] * d.dn[k * nc + i]
				}
				assert math.abs(gd) <= 1e-10, '${c.label}: Gibbs–Duhem column ${i}'
			}
			rt := thermo.r_gas * pt.t
			assert math.tolerance(sum_t, -rp.h / (rt * pt.t), 1e-9), '${c.label}: Σ n·∂lnφ/∂T'
		}
	}
}

fn test_packing_limit_moves_by_c() {
	for c in cases() {
		assert math.tolerance(c.shifted.min_volume(c.n), c.inner.min_volume(c.n) - dot(c.n,
			c.c), 1e-14), c.label
	}
}

fn test_dense_liquid_below_the_inner_covolume_is_found() {
	// With a shift of half the covolume, a compressed liquid's shifted volume
	// lies *below* the inner model's covolume B. The solver only finds it if
	// min_volume really moved down by C.
	pr := thermo.new_cubic(.pr, [butane])!
	b := pr.min_volume([1.0])
	big := thermo.translate(pr, [0.5 * b])!
	pt := Point{300.0, 500e5, .liquid}
	a := solve(pr, pt, [1.0])
	s := solve(big, pt, [1.0])
	assert s.v < b, 'test set-up: the shifted liquid volume should be below B'
	assert math.tolerance(s.v, a.v - 0.5 * b, 1e-10), 'shifted ${s.v}, expected ${a.v - 0.5 * b}'
}

fn test_translate_rejects_bad_shifts() {
	pr := thermo.new_cubic(.pr, [butane])!
	b := pr.min_volume([1.0])
	if _ := thermo.translate(pr, [b]) {
		assert false, 'a shift as large as the covolume should be rejected'
	}
	if _ := thermo.translate(pr, [1e-6, 1e-6]) {
		assert false, 'wrong number of shifts should be rejected'
	}
}
