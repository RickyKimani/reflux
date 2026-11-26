// Internal test (module water): IAPWS-95's φʳ and F against independent
// references (criterion 2): the coefficients against the iapws package's;
// φʳ's six derivatives against mpmath.diff at 50 digits
// (testdata/m5_iapws95.json, tools/gen_m5_iapws95.py), within what rounding
// allows; and F's derivatives against Euler's theorem and differences of F.
module water

import json2
import math
import os

struct Residual95 {
	nr1      []f64
	d1       []f64
	t1       []f64
	nr2      []f64
	c2       []f64
	d2       []f64
	t2       []f64
	nr3      []f64
	d3       []f64
	t3       []f64
	alfa3    []f64
	beta3    []f64
	gamma3   []f64
	epsilon3 []f64
	nr4      []f64
	a4       []f64
	b4       []f64
	aa       []f64 @[json: 'A']
	bb       []f64 @[json: 'B']
	cc       []f64 @[json: 'C']
	dd       []f64 @[json: 'D']
	beta4    []f64
}

struct Ideal95 {
	ao_log []f64
	ao_pow []f64
	ao_exp []f64
	titao  []f64
}

struct Coefficients {
	residual Residual95
	ideal    Ideal95
	r_molar  f64 @[json: 'R_molar']
}

struct DerivRow {
	tau    f64
	delta  f64
	values []f64
	scales []f64
}

struct DerivFixture {
	coefficients Coefficients
	derivatives  []DerivRow
}

fn deriv_fixture() !DerivFixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm5_iapws95.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m5_iapws95.py`.')
	}
	return json2.decode[DerivFixture](text)!
}

fn same(a []f64, b []f64) bool {
	return a.len == b.len && a == b
}

fn ints(a []int) []f64 {
	return a.map(f64(it))
}

fn test_the_coefficients_are_the_iapws_package_s() {
	// D10: water/coefficients.v is transcribed from the release; the iapws
	// package is an independent transcription. Every number equal.
	c := deriv_fixture()!.coefficients
	r := c.residual
	assert same(poly_n, r.nr1) && same(ints(poly_d), r.d1) && same(poly_t, r.t1)
	assert same(expo_n, r.nr2) && same(ints(expo_c), r.c2) && same(ints(expo_d), r.d2)
	assert same(ints(expo_t), r.t2)
	assert same(gauss_n, r.nr3) && same(ints(gauss_d), r.d3) && same(ints(gauss_t), r.t3)
	assert same(gauss_alpha, r.alfa3) && same(gauss_beta, r.beta3)
	assert same(gauss_gamma, r.gamma3) && same(gauss_epsilon, r.epsilon3)
	assert same(nonan_n, r.nr4) && same(nonan_b, r.b4) && same(nonan_c, r.cc)
	assert same(nonan_dd, r.dd)
	for i in 0 .. 2 {
		assert r.a4[i] == nonan_a && r.bb[i] == nonan_bb && r.aa[i] == nonan_aa
		assert r.beta4[i] == nonan_beta
	}
	// Eq. (5): the package writes n₁° and n₂° as ao_pow, ln δ's and ln τ's
	// coefficients as ao_log, and γᵢ° as titao.
	i := c.ideal
	assert same(i.ao_pow, [ideal_n1, ideal_n2]) && same(i.ao_log, [1.0, ideal_n3])
	assert same(i.ao_exp, ideal_n) && same(i.titao, ideal_gamma)
	// And R per mole, which the package rounds to 12 digits.
	assert math.abs(r95 / c.r_molar - 1.0) < 1e-12
}

fn test_phi_r_against_mpmath() {
	// Each derivative within c·ε·scale of mpmath's at 50 digits, the scale
	// being the sum of its terms' magnitudes: what rounding the terms can
	// cost. That sum is at least the value's magnitude, but the fixture takes
	// it at 20 digits, which near the critical point (δ = 1 ± 1e-10) comes
	// out below it, so the scale is the larger of the two. Measured, over
	// 291 points: at most 18.3 ε·scale, in φʳ_δδ at 251.165 K and δ = 1.8,
	// inside the dome, where the τ⁵⁰ terms cancel to 1e-4 of their size;
	// φʳ_ττ at δ = 1 ± 1e-10 within 1 ulp of its value.
	rows := deriv_fixture()!.derivatives
	assert rows.len > 200, 'the fixture has ${rows.len} derivative rows'
	mut worst := []f64{len: 6}
	for row in rows {
		p := phi_residual(row.tau, row.delta)
		got := [p.f, p.d, p.dd, p.t, p.tt, p.dt]
		for k in 0 .. 6 {
			scale := math.max(row.scales[k], math.abs(row.values[k]))
			err := math.abs(got[k] - row.values[k]) / (f64_eps * scale)
			worst[k] = math.max(worst[k], err)
			assert err < 32.0, 'derivative ${k} at τ = ${row.tau}, δ = ${row.delta}: ${got[k]} against ${row.values[k]} (${err} ε·scale)'
		}
	}
	println('φʳ against mpmath, in ε·scale: ${worst.map('${it:.1f}')}')
}

const f64_eps = 2.220446049250313e-16

fn test_f_is_homogeneous() {
	// Euler's theorem for F, of degree 1 in (V, n), and so for its
	// derivatives: F = nF_n + VF_V, nF_nn + VF_Vn = 0, nF_Vn + VF_VV = 0,
	// F_T = nF_Tn + VF_TV. Each holds to rounding of its terms.
	m := iapws95()
	for t in [280.0, 450.0, 647.0, 900.0, 2000.0] {
		for rho in [1e-3, 50.0, 1e4, 3e4, 5.5e4] {
			n := 1.3
			v := n / rho
			r := m.residual2(t, v, [n])
			for pair in [
				[r.f, n * r.f_n[0] + v * r.f_v],
				[0.0, n * r.f_nn[0] + v * r.f_vn[0]],
				[0.0, n * r.f_vn[0] + v * r.f_vv],
				[r.f_t, n * r.f_tn[0] + v * r.f_tv],
			] {
				scale := math.abs(n * r.f_n[0]) + math.abs(v * r.f_v) + math.abs(r.f) +
					math.abs(n * r.f_nn[0]) + math.abs(v * r.f_vv) + math.abs(r.f_t)
				assert math.abs(pair[0] - pair[1]) <= 1e-13 * scale, 'Euler at ${t} K, ${rho} mol/m³: ${pair}'
			}
		}
	}
}

fn test_f_against_its_differences() {
	// Each first and second derivative of F against central differences of
	// the one below it, in V, T and n, by Richardson's extrapolation from h
	// and h/2: the D3 mapping from φʳ's derivatives, checked by F itself.
	m := iapws95()
	for t in [280.0, 450.0, 647.0, 900.0, 2000.0] {
		for rho in [1e-3, 50.0, 1e4, 3e4, 5.5e4] {
			n := 1.0
			v := n / rho
			r := m.residual2(t, v, [n])
			fv := fn [m, t, n] (x f64) []f64 {
				q := m.residual2(t, x, [n])
				return [q.f, q.f_v, q.f_t, q.f_n[0]]
			}
			ft := fn [m, v, n] (x f64) []f64 {
				q := m.residual2(x, v, [n])
				return [q.f, q.f_t]
			}
			fnn := fn [m, t, v] (x f64) []f64 {
				q := m.residual2(t, v, [x])
				return [q.f, q.f_n[0]]
			}
			dv := richardson(fv, v, 1e-4 * v)
			dt := richardson(ft, t, 1e-4 * t)
			dn := richardson(fnn, n, 1e-4 * n)
			checks := [
				[dv[0], r.f_v],
				[dv[1], r.f_vv],
				[dv[2], r.f_tv],
				[dv[3], r.f_vn[0]],
				[dt[0], r.f_t],
				[dt[1], r.f_tt],
				[dn[0], r.f_n[0]],
				[dn[1], r.f_nn[0]],
			]
			for k, c in checks {
				assert math.abs(c[0] - c[1]) <= 1e-7 * math.max(math.abs(c[1]), 1e-6), 'check ${k} at ${t} K, ${rho} mol/m³: ${c[0]} by differences, ${c[1]}'
			}
		}
	}
}

// richardson is the derivative of each of f's outputs at x, by central
// differences at h and h/2 combined to cancel the h² term.
fn richardson(f fn (f64) []f64, x f64, h f64) []f64 {
	a := f(x + h)
	b := f(x - h)
	c := f(x + 0.5 * h)
	d := f(x - 0.5 * h)
	return []f64{len: a.len, init: (4.0 * (c[index] - d[index]) / h - (a[index] - b[index]) / (2.0 * h)) / 3.0}
}
