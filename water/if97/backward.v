module if97

import math

// The backward equations, as published: R7-97's T(p, h) and T(p, s) for
// regions 1 and 2, and SR3-03's T and v of (p, h) and (p, s) for region 3,
// with their boundaries. Each is within its release's tolerance of the
// basic equation (25 mK in region 1, 10 mK in region 2, 25 mK and 0.01%
// in v in region 3), not to rounding: the steam-table functions finish
// each answer on the basic equation (D9).

// sc is the critical entropy, SR3-03's boundary between 3a and 3b in (p, s).
pub const sc = 4.41202148223476e3 // J/(kg·K)

// s2bc is R7-97's boundary between subregions 2b and 2c in (p, s).
pub const s2bc = 5.85e3 // J/(kg·K)

// t1_ph is Eq. (11): θ = Σ nᵢ π^Iᵢ (η + 1)^Jᵢ, p* = 1 MPa, h* = 2500 kJ/kg.
pub fn t1_ph(p f64, h f64) f64 {
	return sum_xy(t1ph_i, t1ph_j, t1ph_n, p / 1e6, h / 2500e3 + 1.0)
}

// t1_ps is Eq. (13): θ = Σ nᵢ π^Iᵢ (σ + 2)^Jᵢ, s* = 1 kJ/(kg·K).
pub fn t1_ps(p f64, s f64) f64 {
	return sum_xy(t1ps_i, t1ps_j, t1ps_n, p / 1e6, s / 1e3 + 2.0)
}

// t2_ph is Eqs. (22) to (24): subregion 2a up to 4 MPa, then 2b and 2c
// either side of B2bc; h* = 2000 kJ/kg.
pub fn t2_ph(p f64, h f64) f64 {
	pi := p / 1e6
	eta := h / 2000e3
	if p <= 4e6 {
		return sum_xy(t2aph_i, t2aph_j, t2aph_n, pi, eta - 2.1)
	}
	if h >= h2bc(p) {
		return sum_xy(t2bph_i, t2bph_j, t2bph_n, pi - 2.0, eta - 2.6)
	}
	return sum_xy(t2cph_i, t2cph_j, t2cph_n, pi + 25.0, eta - 1.8)
}

// t2_ps is Eqs. (25) to (27): subregion 2a up to 4 MPa, then 2b and 2c
// either side of s = 5.85 kJ/(kg·K). 2a's exponents Iᵢ are quarters, taken
// as integer powers of π^0.25.
pub fn t2_ps(p f64, s f64) f64 {
	pi := p / 1e6
	if p <= 4e6 {
		q := math.sqrt(math.sqrt(pi))
		y := s / 2e3 - 2.0
		mut sum := 0.0
		for k in 0 .. t2aps_n.len {
			sum += t2aps_n[k] * ipow(q, int(4.0 * t2aps_i[k])) * ipow(y, t2aps_j[k])
		}
		return sum
	}
	if s >= s2bc {
		return sum_xy(t2bps_i, t2bps_j, t2bps_n, pi, 10.0 - s / 0.7853e3)
	}
	return sum_xy(t2cps_i, t2cps_j, t2cps_n, pi, 2.0 - s / 2.9251e3)
}

// h3ab is SR3-03's Eq. (1), the boundary between 3a and 3b in (p, h):
// η = n₁ + n₂π + n₃π² + n₄π³, p* = 1 MPa, h* = 1 kJ/kg.
pub fn h3ab(p f64) f64 {
	pi := p / 1e6
	return (h3ab_n[0] + pi * (h3ab_n[1] + pi * (h3ab_n[2] + pi * h3ab_n[3]))) * 1e3
}

// t3_ph is SR3-03's Eqs. (2) and (3): p* = 100 MPa; 3a's T* = 760 K and
// h* = 2300 kJ/kg, 3b's 860 K and 2800 kJ/kg.
pub fn t3_ph(p f64, h f64) f64 {
	pi := p / 100e6
	if h <= h3ab(p) {
		return 760.0 * sum_xy(t3aph_i, t3aph_j, t3aph_n, pi + 0.240, h / 2300e3 - 0.615)
	}
	return 860.0 * sum_xy(t3bph_i, t3bph_j, t3bph_n, pi + 0.298, h / 2800e3 - 0.720)
}

// v3_ph is SR3-03's Eqs. (4) and (5): 3a's v* = 0.0028 m³/kg and h* =
// 2100 kJ/kg, 3b's 0.0088 m³/kg and 2800 kJ/kg.
pub fn v3_ph(p f64, h f64) f64 {
	pi := p / 100e6
	if h <= h3ab(p) {
		return 0.0028 * sum_xy(v3aph_i, v3aph_j, v3aph_n, pi + 0.128, h / 2100e3 - 0.727)
	}
	return 0.0088 * sum_xy(v3bph_i, v3bph_j, v3bph_n, pi + 0.0661, h / 2800e3 - 0.720)
}

// t3_ps is SR3-03's Eqs. (6) and (7): 3a's T* = 760 K and s* = 4.4
// kJ/(kg·K), 3b's 860 K and 5.3 kJ/(kg·K), split at sc.
pub fn t3_ps(p f64, s f64) f64 {
	pi := p / 100e6
	if s <= sc {
		return 760.0 * sum_xy(t3aps_i, t3aps_j, t3aps_n, pi + 0.240, s / 4.4e3 - 0.703)
	}
	return 860.0 * sum_xy(t3bps_i, t3bps_j, t3bps_n, pi + 0.760, s / 5.3e3 - 0.818)
}

// v3_ps is SR3-03's Eqs. (8) and (9).
pub fn v3_ps(p f64, s f64) f64 {
	pi := p / 100e6
	if s <= sc {
		return 0.0028 * sum_xy(v3aps_i, v3aps_j, v3aps_n, pi + 0.187, s / 4.4e3 - 0.755)
	}
	return 0.0088 * sum_xy(v3bps_i, v3bps_j, v3bps_n, pi + 0.298, s / 5.3e3 - 0.816)
}

// p3sat_h is SR3-03's Eq. (10), the saturation pressure from h along
// region 3's part of the saturation line: p* = 22 MPa, h* = 2600 kJ/kg.
pub fn p3sat_h(h f64) f64 {
	eta := h / 2600e3
	return 22e6 * sum_xy(p3sath_i, p3sath_j, p3sath_n, eta - 1.02, eta - 0.608)
}

// p3sat_s is SR3-03's Eq. (11): p* = 22 MPa, s* = 5.2 kJ/(kg·K).
pub fn p3sat_s(s f64) f64 {
	sigma := s / 5.2e3
	return 22e6 * sum_xy(p3sats_i, p3sats_j, p3sats_n, sigma - 1.03, sigma - 0.699)
}
