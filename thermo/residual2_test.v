// Checks for residual2: every second derivative of F against finite
// differences of the first derivatives, plus the identities that follow from
// homogeneity. An external test, so it sees thermo only through its pub API.
import math
import thermo { Component, CubicKind }

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
	label string
	eos   thermo.Cubic
	n     []f64
}

// cases covers every α variant and both pure and mixed models. The ternaries
// have different kᵢⱼ off the diagonal, so a transposed index shows up.
fn cases() []Case {
	mut out := []Case{}
	for kind in [CubicKind.vdw, .rk, .srk, .pr] {
		out << Case{
			label: '${kind} pure'
			eos:   thermo.new_cubic(kind, [butane]) or { panic(err) }
			n:     [1.0]
		}
	}
	kij := [0.0, 0.02, 0.10, 0.02, 0.0, 0.13, 0.10, 0.13, 0.0]
	out << Case{
		label: 'PR ternary, PRSV butane'
		eos:   thermo.new_cubic(.pr, [methane, butane, co2],
			kij:    kij
			alphas: [thermo.soave_pr(methane.omega), thermo.prsv(butane.omega, 0.04),
				thermo.soave_pr(co2.omega)]
		) or { panic(err) }
		n:     [0.2, 0.5, 0.3]
	}
	// Twu on methane, Mathias–Copeman on CO₂. At 250 K CO₂ is below its Tc
	// (304 K) and at 400 K above it, so both MC branches are exercised.
	// The parameter values are illustrative.
	out << Case{
		label: 'PR ternary, Twu + Mathias-Copeman'
		eos:   thermo.new_cubic(.pr, [methane, butane, co2],
			kij:    kij
			alphas: [thermo.Alpha(thermo.Twu{
				l: 0.15
				m: 0.90
				n: 1.9
			}), thermo.soave_pr(butane.omega), thermo.MathiasCopeman{
				c1: 0.70
				c2: -0.3
				c3: 0.5
			}]
		) or { panic(err) }
		n:     [0.2, 0.5, 0.3]
	}
	out << Case{
		label: 'SRK ternary'
		eos:   thermo.new_cubic(.srk, [methane, butane, co2], kij: kij) or { panic(err) }
		n:     [0.6, 0.1, 0.3]
	}
	return out
}

// close compares an analytic value with a finite-difference estimate. scale
// is the size of the quantity's family (e.g. the largest |F_nᵢnⱼ|), so that a
// near-zero entry is judged against its neighbours, not against itself.
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

// Each case is checked at a liquid-like and a gas-like volume, relative to
// the covolume B, and at two temperatures.
fn states(c Case) [][]f64 {
	b := c.eos.min_volume(c.n)
	mut out := [][]f64{}
	for t in [250.0, 400.0] {
		for ratio in [1.6, 40.0] {
			out << [t, ratio * b]
		}
	}
	return out
}

fn test_first_order_fields_agree_with_residual() {
	for c in cases() {
		for st in states(c) {
			r1 := c.eos.residual(st[0], st[1], c.n)
			r2 := c.eos.residual2(st[0], st[1], c.n)
			assert r2.f == r1.f && r2.f_v == r1.f_v && r2.f_vv == r1.f_vv, c.label
			assert math.tolerance(r2.f_t, r1.f_t, 1e-14), c.label
			for i in 0 .. c.n.len {
				assert math.tolerance(r2.f_n[i], r1.f_n[i], 1e-14), c.label
			}
		}
	}
}

fn test_second_derivatives_match_finite_differences() {
	for c in cases() {
		nc := c.n.len
		for st in states(c) {
			t, v := st[0], st[1]
			r := c.eos.residual2(t, v, c.n)
			ht := 1e-4 * t
			hv := 1e-6 * v
			tp := c.eos.residual(t + ht, v, c.n)
			tm := c.eos.residual(t - ht, v, c.n)
			vp := c.eos.residual(t, v + hv, c.n)
			vm := c.eos.residual(t, v - hv, c.n)
			at := '${c.label} at T=${t}, V=${v}'

			fd_tt := (tp.f_t - tm.f_t) / (2.0 * ht)
			assert close(r.f_tt, fd_tt, 0.0), '${at}: F_TT ${r.f_tt} vs ${fd_tt}'
			fd_tv := (tp.f_v - tm.f_v) / (2.0 * ht)
			assert close(r.f_tv, fd_tv, 0.0), '${at}: F_TV ${r.f_tv} vs ${fd_tv}'

			vn_scale := max_abs(r.f_vn)
			tn_scale := max_abs(r.f_tn)
			nn_scale := max_abs(r.f_nn)
			for i in 0 .. nc {
				// F_Vnᵢ two ways: d/dV of F_nᵢ, and d/dnᵢ of F_V (Schwarz).
				fd_vn := (vp.f_n[i] - vm.f_n[i]) / (2.0 * hv)
				assert close(r.f_vn[i], fd_vn, vn_scale), '${at}: F_Vn[${i}] ${r.f_vn[i]} vs ${fd_vn}'
				fd_tn := (tp.f_n[i] - tm.f_n[i]) / (2.0 * ht)
				assert close(r.f_tn[i], fd_tn, tn_scale), '${at}: F_Tn[${i}] ${r.f_tn[i]} vs ${fd_tn}'

				hn := 1e-6 * numeric_sum(c.n)
				mut up := c.n.clone()
				mut dn := c.n.clone()
				up[i] += hn
				dn[i] -= hn
				np := c.eos.residual(t, v, up)
				nm := c.eos.residual(t, v, dn)
				fd_vn2 := (np.f_v - nm.f_v) / (2.0 * hn)
				assert close(r.f_vn[i], fd_vn2, vn_scale), '${at}: F_Vn[${i}] via F_V'
				for j in 0 .. nc {
					// Column i of the Hessian: d/dnᵢ of F_nⱼ.
					fd_nn := (np.f_n[j] - nm.f_n[j]) / (2.0 * hn)
					assert close(r.f_nn[j * nc + i], fd_nn, nn_scale), '${at}: F_nn[${j}][${i}] ${r.f_nn[j * nc + i]} vs ${fd_nn}'
				}
			}
		}
	}
}

fn numeric_sum(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

fn test_hessian_is_symmetric() {
	for c in cases() {
		nc := c.n.len
		for st in states(c) {
			r := c.eos.residual2(st[0], st[1], c.n)
			for i in 0 .. nc {
				for j in 0 .. i {
					assert math.tolerance(r.f_nn[i * nc + j], r.f_nn[j * nc + i], 1e-14), c.label
				}
			}
		}
	}
}

fn test_euler_identities_for_the_second_derivatives() {
	// F_V and F_nᵢ are homogeneous of degree 0 in (V, n), and F_T of
	// degree 1, so:
	//   V·F_VV  + Σᵢ nᵢ·F_Vnᵢ   = 0
	//   V·F_Vnᵢ + Σⱼ nⱼ·F_nᵢnⱼ = 0
	//   V·F_TV  + Σᵢ nᵢ·F_Tnᵢ   = F_T
	for c in cases() {
		nc := c.n.len
		for st in states(c) {
			v := st[1]
			r := c.eos.residual2(st[0], v, c.n)
			mut s1 := v * r.f_vv
			mut size1 := math.abs(v * r.f_vv)
			mut s3 := v * r.f_tv
			mut size3 := math.abs(v * r.f_tv)
			for i in 0 .. nc {
				s1 += c.n[i] * r.f_vn[i]
				size1 += math.abs(c.n[i] * r.f_vn[i])
				s3 += c.n[i] * r.f_tn[i]
				size3 += math.abs(c.n[i] * r.f_tn[i])
				mut s2 := v * r.f_vn[i]
				mut size2 := math.abs(s2)
				for j in 0 .. nc {
					s2 += c.n[j] * r.f_nn[i * nc + j]
					size2 += math.abs(c.n[j] * r.f_nn[i * nc + j])
				}
				assert math.abs(s2) <= 1e-12 * size2, '${c.label}: V·F_Vn + Σ n·F_nn = ${s2}'
			}
			assert math.abs(s1) <= 1e-12 * size1, '${c.label}: V·F_VV + Σ n·F_Vn = ${s1}'
			assert math.abs(s3 - r.f_t) <= 1e-12 * size3, '${c.label}: V·F_TV + Σ n·F_Tn − F_T = ${s3 - r.f_t}'
		}
	}
}

// The α derivative checks are in alpha_test.v.
