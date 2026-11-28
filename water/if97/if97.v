// IAPWS-IF97, the industrial formulation for water and steam (M5a's D9):
// the basic equations of its five regions, its backward equations and the
// region-3 supplements, and the steam-table functions over them. Everything
// is per kg in SI, as the release is (J/kg, J/(kg·K), m³/kg), with IF97's
// own gas constant, on IAPWS's reference (u = s = 0 for the triple point's
// liquid). docs/design.md, "IAPWS-IF97", has the equations.
module if97

import math

// R7-97's constants, Eqs. (1) to (4).
pub const r = 461.526 // J/(kg·K), IF97's gas constant
pub const tc = 647.096 // K
pub const pc = 22.064e6 // Pa
pub const rhoc = 322.0 // kg/m³

// The formulation's range (R7-97, Section 3): 273.15 K to 1073.15 K up to
// 100 MPa, and to 2273.15 K up to 50 MPa.
pub const t_min = 273.15
pub const t_max = 2273.15
pub const p_max = 100e6
pub const t5 = 1073.15 // region 5 above it
pub const p5_max = 50e6
pub const t13 = 623.15 // regions 1 and 3 meet on this isotherm
pub const t23_max = 863.15 // B23 reaches 100 MPa here

// A state: the region, T and p, and the properties per kg. In region 4, the
// two-phase mixture, x is the vapour's mass fraction and cp, cv, w and
// drhodp are NaN unless x is 0 or 1, a saturated phase; elsewhere x is NaN.
pub struct State {
pub:
	region int
	t      f64 // K
	p      f64 // Pa
	v      f64 // m³/kg
	h      f64 // J/kg
	s      f64 // J/(kg·K)
	u      f64 // J/kg
	cp     f64 // J/(kg·K)
	cv     f64 // J/(kg·K)
	w      f64 // m/s
	drhodp f64 // (∂ρ/∂p)_T, kg/(m³·Pa): R15-11's industrial λ needs it
	dvdt   f64 // (∂v/∂T)_p, m³/(kg·K): a phase's, for M5b's phase model
	x      f64 = math.nan()
}

pub fn (st State) rho() f64 {
	return 1.0 / st.v
}

pub enum Failure {
	out_of_range // the state isn't in the formulation's range
	no_answer    // a search that didn't meet its tolerance
}

// If97Error is a refusal: a state outside the range (with the range), or a
// search that failed (with how far it got).
pub struct If97Error {
	Error
pub:
	kind Failure
	what string
}

pub fn (e If97Error) msg() string {
	return e.what
}

fn out_of_range(what string) IError {
	return If97Error{
		kind: .out_of_range
		what: what
	}
}

fn no_answer(what string) IError {
	return If97Error{
		kind: .no_answer
		what: what
	}
}

// ipow is x^n by repeated squaring, for the equations' integer exponents.
@[inline]
fn ipow(x f64, n int) f64 {
	mut k := if n < 0 { -n } else { n }
	mut b := x
	mut acc := 1.0
	for k > 0 {
		if k & 1 == 1 {
			acc *= b
		}
		b *= b
		k >>= 1
	}
	return if n < 0 { 1.0 / acc } else { acc }
}

// sum_xy is Σ nᵢ x^Iᵢ y^Jᵢ, the form of every backward equation.
fn sum_xy(ii []int, jj []int, nn []f64, x f64, y f64) f64 {
	mut s := 0.0
	for k in 0 .. nn.len {
		s += nn[k] * ipow(x, ii[k]) * ipow(y, jj[k])
	}
	return s
}

// psat97 is region 4's saturation pressure, Eq. (30), with no range check:
// p in Pa from T in K.
pub fn psat97(t f64) f64 {
	theta := t + r4_n[8] / (t - r4_n[9])
	a := theta * theta + r4_n[0] * theta + r4_n[1]
	b := r4_n[2] * theta * theta + r4_n[3] * theta + r4_n[4]
	c := r4_n[5] * theta * theta + r4_n[6] * theta + r4_n[7]
	q := 2.0 * c / (-b + math.sqrt(b * b - 4.0 * a * c))
	return q * q * q * q * 1e6
}

// tsat97 is region 4's saturation temperature, Eq. (31), with no range
// check: T in K from p in Pa.
pub fn tsat97(p f64) f64 {
	beta := math.sqrt(math.sqrt(p / 1e6))
	e := beta * beta + r4_n[2] * beta + r4_n[5]
	f := r4_n[0] * beta * beta + r4_n[3] * beta + r4_n[6]
	g := r4_n[1] * beta * beta + r4_n[4] * beta + r4_n[7]
	d := 2.0 * g / (-f - math.sqrt(f * f - 4.0 * e * g))
	nd := r4_n[9] + d
	return 0.5 * (nd - math.sqrt(nd * nd - 4.0 * (r4_n[8] + r4_n[9] * d)))
}

// psat is the saturation pressure from 273.15 K to Tc, refused outside.
pub fn psat(t f64) !f64 {
	if !(t >= t_min && t <= tc) {
		return out_of_range('IF97 has a saturation pressure from ${t_min} K to ${tc} K, not at ${t} K')
	}
	return psat97(t)
}

// tsat is the saturation temperature from psat(273.15 K) to pc.
pub fn tsat(p f64) !f64 {
	if !(p >= psat97(t_min) && p <= pc) {
		return out_of_range('IF97 has a saturation temperature from ${psat97(t_min)} Pa to ${pc} Pa, not at ${p} Pa')
	}
	return tsat97(p)
}

// p23 is B23, the boundary between regions 2 and 3, Eq. (5): p in Pa
// from T in K, 623.15 K to 863.15 K.
pub fn p23(t f64) f64 {
	return (b23_n[0] + b23_n[1] * t + b23_n[2] * t * t) * 1e6
}

// t23 is B23's inverse, Eq. (6).
pub fn t23(p f64) f64 {
	return b23_n[3] + math.sqrt((p / 1e6 - b23_n[4]) / b23_n[2])
}

// p2bc is B2bc, between subregions 2b and 2c, Eq. (20): p in Pa from h in
// J/kg.
pub fn p2bc(h f64) f64 {
	eta := h / 1e3
	return (b2bc_n[0] + b2bc_n[1] * eta + b2bc_n[2] * eta * eta) * 1e6
}

// h2bc is B2bc's inverse, Eq. (21).
pub fn h2bc(p f64) f64 {
	return (b2bc_n[3] + math.sqrt((p / 1e6 - b2bc_n[4]) / b2bc_n[2])) * 1e3
}
