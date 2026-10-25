// An internal test (module thermo), so it can reach the private
// critical_constants and CubicKind.consts.
module thermo

import math

fn test_vdw_constants_are_the_exact_fractions() {
	k := CubicKind.vdw.consts()
	assert math.tolerance(k.omega_a, 27.0 / 64.0, 1e-15)
	assert math.tolerance(k.omega_b, 1.0 / 8.0, 1e-15)
	assert math.tolerance(k.zc, 3.0 / 8.0, 1e-15)
}

fn test_srk_constants_match_their_closed_form() {
	// For δ₁ = 1, δ₂ = 0 the equation in B is (3B + 1)³ = 2.
	c := math.cbrt(2.0)
	k := CubicKind.srk.consts()
	assert math.tolerance(k.omega_b, (c - 1.0) / 3.0, 1e-14)
	assert math.tolerance(k.omega_a, 1.0 / (9.0 * (c - 1.0)), 1e-14)
	assert math.tolerance(k.zc, 1.0 / 3.0, 1e-14)
	assert CubicKind.rk.consts() == k
}

fn test_pr_constants_match_the_literature() {
	k := CubicKind.pr.consts()
	// Peng & Robinson's η = b/Vc has the closed form
	// (−1 + ∛(6√2 + 8) − ∛(6√2 − 8)) / 3, and Ωb = η·Zc.
	eta := (-1.0 + math.cbrt(6.0 * math.sqrt2 + 8.0) - math.cbrt(6.0 * math.sqrt2 - 8.0)) / 3.0
	assert math.tolerance(k.omega_b, eta * k.zc, 1e-13)
	// The ten-figure values quoted in the literature (and used by `thermo`).
	assert math.abs(k.omega_a - 0.4572355289) < 1e-10
	assert math.abs(k.omega_b - 0.0777960739) < 1e-10
	assert math.abs(k.zc - 0.3074013087) < 1e-10
}

fn test_every_cubic_has_a_triple_root_at_its_critical_point() {
	// At (Tc, Pc) the polynomial is (Z − Zc)³ up to rounding in the last bit.
	// A triple root moves with the cube root of any coefficient error, so
	// 1e-16 of rounding shows up as about 5e-6 in Z.
	c := Component{
		name:  'test'
		tc:    400.0
		pc:    40e5
		omega: 0.1
		mw:    0.1
	}
	for kind in [CubicKind.vdw, .rk, .srk, .pr] {
		eos := new_cubic(kind, [c])!
		zc := kind.consts().zc
		roots := eos.z_roots(c.tc, c.pc, [1.0])
		assert roots.len >= 1, '${kind}: no root'
		for z in roots {
			assert math.abs(z - zc) < 1e-5, '${kind}: root ${z}, Zc ${zc}'
		}
	}
}

fn test_equal_deltas_are_van_der_waals_shifted_in_volume() {
	// With δ₁ = δ₂ = δ the attraction is a/(V + δb)², i.e. van der Waals in
	// the shifted volume V + δb. So Ωa is always 27/64, and only Ωb moves.
	for d in [1.0, 3.0, 10.0] {
		k := critical_constants(d, d)!
		assert math.tolerance(k.omega_a, 27.0 / 64.0, 1e-14)
		assert math.tolerance(k.omega_b, 1.0 / (8.0 * (1.0 + d)), 1e-14)
	}
}

fn test_critical_constants_reject_impossible_deltas() {
	// With δ₁ = δ₂ = −1 the attraction a/(V − b)² diverges at V = b, so there
	// is no critical point. That must be an error, not garbage.
	if _ := critical_constants(-1.0, -1.0) {
		assert false, 'expected an error'
	}
}
