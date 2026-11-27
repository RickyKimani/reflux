module water

import math
import thermo

// IAPWS-95's reference constants (its Eqs. (1) to (3)), and the molar mass
// IAPWS G5-01 gives for VSMOW, which turns its per-kg formulation into
// reflux's per-mole one.
pub const critical_t = 647.096 // K

pub const critical_p = 22.064e6 // Pa, IAPWS's critical pressure (R2-83)

pub const critical_rho_mass = 322.0 // kg/m³

pub const molar_mass = 0.018015268 // kg/mol

// critical_rho is ρc in mol/m³.
pub const critical_rho = critical_rho_mass / molar_mass

// r95 is IAPWS-95's own gas constant, R = 0.46151805 kJ/(kg·K), per mole:
// 8.314371357587 J/(mol·K). The release requires it: its coefficients were
// fitted with it. It's 1.1e-5 below reflux's, and D4 carries the difference
// (docs/m5-plan.md).
pub const r95 = 461.51805 * molar_mass

// k95 is R₉₅/R, the factor F takes φʳ by (D4).
const k95 = r95 / thermo.r_gas

// acentric is ω from IAPWS-95's own vapour pressure at 0.7·Tc,
// −1 − log₁₀(Psat/Pc), rounded: 0.3443 (water_test.v checks it). Only
// estimates use it, Wilson's K and Edmister's start.
pub const acentric = 0.3443

// delta_max is the reduced density of IAPWS-95's min_volume (D7, as
// amended in step 2): there P is 1.16 GPa or more at every T from 251 K,
// above the formulation's 1 GPa, and four such volumes are 1/ρc, the volume
// thermo.phase_label takes for the critical one.
const delta_max = 4.0

// Iapws95 is the IAPWS Formulation 1995 for water, as a thermo.Model (D3):
//
//     F(T, V, n) = n·(R₉₅/R)·φʳ(τ, δ) + n·(R₉₅/R − 1)·ln δ
//
// with δ = n/(Vρc) and τ = Tc/T. The second term carries the difference
// between IAPWS-95's gas constant and reflux's (D4), so that P, and every
// property at (T, ρ), are the release's. It has one component, water.
pub struct Iapws95 {
	comps []thermo.Component
}

// iapws95 makes the model.
pub fn iapws95() Iapws95 {
	return Iapws95{
		comps: [
			thermo.Component{
				name:  'water'
				tc:    critical_t
				pc:    critical_p
				omega: acentric
				mw:    molar_mass
			},
		]
	}
}

// components implements thermo.Model.
pub fn (m Iapws95) components() []thermo.Component {
	return m.comps
}

// min_volume implements thermo.Model. IAPWS-95 has no covolume: P rises
// without bound as δ grows, like δ⁵. It's put at δ = 4 (delta_max).
pub fn (m Iapws95) min_volume(n []f64) f64 {
	return n[0] / (delta_max * critical_rho)
}

// volume_function implements thermo.Model: the parts of φʳ in τ alone are
// computed once, here.
pub fn (m Iapws95) volume_function(t f64, n []f64) thermo.VolumeFunction {
	return AtT{
		tau:  tau_parts(critical_t / t)
		t:    t
		ntot: n[0]
	}
}

// branches implements thermo.Branches (D7): below Tc the vapour's branch
// lies below ρc and the liquid's above. Above Tc an isotherm has one
// branch, and no spurious loops (M5a step 0). Between the vapour's spinodal
// and ρc, P rises only where it's hugely negative (−1e27 to −5e8 Pa up to
// 600 K) or, at 625 to 640 K, at pressures whose ideal-gas start lies far
// below that stretch, on the vapour's own branch (M5a step 2), so the split
// is enough: the vapour solver needs no start of its own.
//
// Below 251.165 K, the lowest point of the melting curve (IAPWS R14-08),
// reflux's IAPWS-95 has no state (M5a step 2). The release's range starts at
// the melting curve. Below it the liquid branch is an extrapolation whose
// enthalpy falls as T rises (at 562 MPa, H at 215.7 K is above H at
// 273.16 K), which a search in T can't use; and the gas, which the release
// endorses to 130 K below the sublimation curve, gave a UV flash a vapour at
// 233.8 K and 522 Pa, far above that curve, with the U and V of a two-phase
// state at 273.16 K. Without R14-08's curves (not in M5a), no state at all.
pub fn (m Iapws95) branches(t f64, n []f64) thermo.BranchLimits {
	if t >= critical_t {
		return thermo.BranchLimits{}
	}
	return thermo.BranchLimits{
		split:    n[0] / critical_rho
		no_state: t < liquid_t_min
	}
}

// liquid_t_min is the lowest T of reflux's IAPWS-95, K: the melting curve's
// lowest point, where ice Ih, ice III and the liquid meet at 208.566 MPa.
pub const liquid_t_min = 251.165

// residual implements thermo.Model.
pub fn (m Iapws95) residual(t f64, v f64, n []f64) thermo.Residual {
	return m.volume_function(t, n).residual(v)
}

// residual2 implements thermo.Model, from φʳ's six derivatives (D3). With
// ψ = k·φʳ + (k − 1)·ln δ, k = R₉₅/R, F = nψ and
//
//     F_V  = −(nδ/V)·ψ_δ           F_VV = (n/V²)(2δψ_δ + δ²ψ_δδ)
//     F_T  = −(nτ/T)·ψ_τ           F_TT = (n/T²)(τ²ψ_ττ + 2τψ_τ)
//     F_TV = (nδτ/(VT))·ψ_δτ       F_n  = ψ + δψ_δ
//     F_Vn = −(δ/V)(2ψ_δ + δψ_δδ)  F_Tn = −(τ/T)(ψ_τ + δψ_δτ)
//     F_nn = (δ/n)(2ψ_δ + δψ_δδ)
pub fn (m Iapws95) residual2(t f64, v f64, n []f64) thermo.Residual2 {
	ntot := n[0]
	tau := critical_t / t
	delta := ntot / (v * critical_rho)
	p := phi_r(tau_parts(tau), delta, true)
	s := psi_of(p, delta)
	two := 2.0 * s.d + delta * s.dd
	return thermo.Residual2{
		Residual: thermo.Residual{
			f:    ntot * s.f
			f_v:  -ntot * delta / v * s.d
			f_vv: ntot / (v * v) * delta * two
			f_t:  -ntot * tau / t * s.t
			f_n:  [s.f + delta * s.d]
		}
		f_tt:     ntot / (t * t) * (tau * tau * s.tt + 2.0 * tau * s.t)
		f_tv:     ntot * delta * tau / (v * t) * s.dt
		f_vn:     [-delta / v * two]
		f_tn:     [-tau / t * (s.t + delta * s.dt)]
		f_nn:     [delta / ntot * two]
	}
}

// AtT is the model at a fixed T and n: its thermo.VolumeFunction.
struct AtT {
	tau  TauParts
	t    f64
	ntot f64
}

fn (s AtT) at(v f64) thermo.VolumeDerivs {
	delta := s.ntot / (v * critical_rho)
	p := phi_r(s.tau, delta, false)
	ps := psi_of(p, delta)
	return thermo.VolumeDerivs{
		f:    s.ntot * ps.f
		f_v:  -s.ntot * delta / v * ps.d
		f_vv: s.ntot / (v * v) * delta * (2.0 * ps.d + delta * ps.dd)
	}
}

fn (s AtT) residual(v f64) thermo.Residual {
	delta := s.ntot / (v * critical_rho)
	p := phi_r(s.tau, delta, true)
	ps := psi_of(p, delta)
	tau := s.tau.tau
	return thermo.Residual{
		f:    s.ntot * ps.f
		f_v:  -s.ntot * delta / v * ps.d
		f_vv: s.ntot / (v * v) * delta * (2.0 * ps.d + delta * ps.dd)
		f_t:  -s.ntot * tau / s.t * ps.t
		f_n:  [ps.f + delta * ps.d]
	}
}

// psi_of is ψ = k·φʳ + (k − 1)·ln δ with its derivatives, k = R₉₅/R (D4).
fn psi_of(p Phi, delta f64) Phi {
	return Phi{
		f:  k95 * p.f + (k95 - 1.0) * math.log(delta)
		d:  k95 * p.d + (k95 - 1.0) / delta
		dd: k95 * p.dd - (k95 - 1.0) / (delta * delta)
		t:  k95 * p.t
		tt: k95 * p.tt
		dt: k95 * p.dt
	}
}

// Phi is a dimensionless Helmholtz energy, φ° or φʳ, with its derivatives:
// f, ∂/∂δ, ∂²/∂δ², ∂/∂τ, ∂²/∂τ², ∂²/∂δ∂τ.
pub struct Phi {
pub:
	f  f64
	d  f64
	dd f64
	t  f64
	tt f64
	dt f64
}

// TauParts is what φʳ needs of τ alone, at one T: each term's coefficient
// times its power of τ, the Gaussian terms' factors in τ, and the
// non-analytic terms' 1 − τ and τ − 1.
struct TauParts {
	tau f64
mut:
	poly  [7]f64 // nᵢτ^tᵢ
	expo  [44]f64
	gauss [3]f64 // nᵢτ^tᵢ·exp(−βᵢ(τ − γᵢ)²)
	gx    [3]f64 // ∂ln(…)/∂τ = tᵢ/τ − 2βᵢ(τ − γᵢ)
	gxx   [3]f64 // gxᵢ² − tᵢ/τ² − 2βᵢ
}

fn tau_parts(tau f64) TauParts {
	mut p := TauParts{
		tau: tau
	}
	for i in 0 .. 7 {
		p.poly[i] = poly_n[i] * math.pow(tau, poly_t[i])
	}
	for i in 0 .. 44 {
		p.expo[i] = expo_n[i] * math.pow(tau, f64(expo_t[i]))
	}
	for i in 0 .. 3 {
		dt := tau - gauss_gamma[i]
		t := f64(gauss_t[i])
		p.gauss[i] = gauss_n[i] * math.pow(tau, t) * math.exp(-gauss_beta[i] * dt * dt)
		p.gx[i] = t / tau - 2.0 * gauss_beta[i] * dt
		p.gxx[i] = p.gx[i] * p.gx[i] - t / (tau * tau) - 2.0 * gauss_beta[i]
	}
	return p
}

// phi_r is IAPWS-95's residual part at (τ, δ), from its Table 5: φʳ, φʳ_δ and
// φʳ_δδ, and with all set the τ-derivatives too.
//
// Two places differ from the table as printed, for its numbers' sake:
// - the exponential terms' φ_δδ has the factor d − 1 − cδᶜ written as it
//   stands, not as (d − cδᶜ) − 1, which cancels at low density (M5a step 0:
//   1.2e-7 at δ = 1.7e-8);
// - the non-analytic terms' ∂²Δ/∂δ² divides ∂Δ/∂δ by δ − 1, which is 0/0
//   at δ = 1 for every τ; written out, the quotient is a sum of positive
//   powers of (δ − 1)², finite and 0 there (D8).
//
// At the critical point itself, where Δ = 0, the non-analytic terms and
// their first derivatives take their limits, 0, as does their part of
// φʳ_δδ and φʳ_δτ, and φʳ_ττ is −∞, IAPWS-95's own divergence (D8).
fn phi_r(tp TauParts, delta f64, all bool) Phi {
	tau := tp.tau
	// δ⁰ to δ¹⁵, the integer powers every term but the Gaussians' uses.
	mut dp := [16]f64{}
	dp[0] = 1.0
	for i in 1 .. 16 {
		dp[i] = dp[i - 1] * delta
	}
	mut f := 0.0
	mut fd := 0.0
	mut fdd := 0.0
	mut ft := 0.0
	mut ftt := 0.0
	mut fdt := 0.0
	for i in 0 .. 7 {
		d := poly_d[i]
		t := poly_t[i]
		base := tp.poly[i] * dp[d]
		f += base
		fd += tp.poly[i] * f64(d) * dp[d - 1]
		if d >= 2 {
			fdd += tp.poly[i] * f64(d * (d - 1)) * dp[d - 2]
		}
		if all {
			ft += base * t / tau
			ftt += base * t * (t - 1.0) / (tau * tau)
			fdt += tp.poly[i] * f64(d) * dp[d - 1] * t / tau
		}
	}
	// e^(−δᶜ) for the five c the exponential terms have.
	mut ec := [7]f64{}
	ec[1] = math.exp(-dp[1])
	ec[2] = math.exp(-dp[2])
	ec[3] = math.exp(-dp[3])
	ec[4] = math.exp(-dp[4])
	ec[6] = math.exp(-dp[6])
	for i in 0 .. 44 {
		c := expo_c[i]
		d := expo_d[i]
		t := f64(expo_t[i])
		dc := dp[c]
		a := f64(d) - f64(c) * dc
		head := tp.expo[i] * ec[c] // nτ^t·e^(−δᶜ)
		f += head * dp[d]
		fd += head * dp[d - 1] * a
		// δ^(d−2)·(a(d − 1 − cδᶜ) − c²δᶜ); for d = 1 the bracket carries a
		// factor δᶜ, so dividing δ^d by δ² stays exact enough and finite.
		fdd += head * dp[d] / (delta * delta) * (a * (f64(d) - 1.0 - f64(c) * dc) - f64(c * c) * dc)
		if all {
			ft += head * dp[d] * t / tau
			ftt += head * dp[d] * t * (t - 1.0) / (tau * tau)
			fdt += head * dp[d - 1] * a * t / tau
		}
	}
	for i in 0 .. 3 {
		d := gauss_d[i]
		de := delta - gauss_epsilon[i]
		base := tp.gauss[i] * dp[d] * math.exp(-gauss_alpha[i] * de * de)
		xd := f64(d) / delta - 2.0 * gauss_alpha[i] * de
		f += base
		fd += base * xd
		fdd += base * (xd * xd - f64(d) / (delta * delta) - 2.0 * gauss_alpha[i])
		if all {
			ft += base * tp.gx[i]
			ftt += base * tp.gxx[i]
			fdt += base * xd * tp.gx[i]
		}
	}
	// The non-analytic terms. With d2 = (δ − 1)² and q = d2^(1/(2β) − 1),
	// d2^(1/2β) = d2·q and d2^(1/β − 1) = d2·q², so one power gives all three.
	dm := delta - 1.0
	d2 := dm * dm
	q := if d2 > 0.0 { math.pow(d2, 0.5 / nonan_beta - 1.0) } else { 0.0 }
	r := if d2 > 0.0 { math.pow(d2, nonan_a - 1.0) } else { 0.0 }
	theta := (1.0 - tau) + nonan_aa * d2 * q
	big_delta := theta * theta + nonan_bb * d2 * r
	// ∂Δ/∂δ over (δ − 1), and ∂²Δ/∂δ² with that quotient in it.
	g := nonan_aa * theta * 2.0 / nonan_beta * q + 2.0 * nonan_bb * nonan_a * r
	delta_d := dm * g
	delta_dd := g + 4.0 * nonan_bb * nonan_a * (nonan_a - 1.0) * r +
		2.0 * nonan_aa * nonan_aa / (nonan_beta * nonan_beta) * d2 * q * q +
		nonan_aa * theta * 4.0 / nonan_beta * (0.5 / nonan_beta - 1.0) * q
	mut tt_at_critical := 0.0 // φʳ_ττ's sign at Δ = 0, from the smallest b
	mut b_least := math.inf(1)
	for i in 0 .. 2 {
		n := nonan_n[i]
		b := nonan_b[i]
		cc := nonan_c[i]
		dd := nonan_dd[i]
		tm := tau - 1.0
		psi := math.exp(-cc * d2 - dd * tm * tm)
		psi_d := -2.0 * cc * dm * psi
		psi_dd := (2.0 * cc * d2 - 1.0) * 2.0 * cc * psi
		if big_delta == 0.0 {
			// The critical point: Δ^b and its first derivatives are 0 in the
			// limit, and so are its parts of φʳ_δδ and φʳ_δτ; φʳ_ττ diverges
			// as Δ^(b−1), so the term of smallest b sets its sign.
			if all && b < b_least {
				b_least = b
				tt_at_critical = if n < 0.0 { -1.0 } else { 1.0 }
			}
			continue
		}
		db := math.pow(big_delta, b)
		db1 := db / big_delta
		db2 := db1 / big_delta
		dbd := b * db1 * delta_d
		dbdd := b * (db1 * delta_dd + (b - 1.0) * db2 * delta_d * delta_d)
		f += n * db * delta * psi
		fd += n * (db * (psi + delta * psi_d) + dbd * delta * psi)
		fdd += n * (db * (2.0 * psi_d + delta * psi_dd) + 2.0 * dbd * (psi + delta * psi_d) +
			dbdd * delta * psi)
		if all {
			psi_t := -2.0 * dd * tm * psi
			psi_tt := (2.0 * dd * tm * tm - 1.0) * 2.0 * dd * psi
			psi_dt := 4.0 * cc * dd * dm * tm * psi
			dbt := -2.0 * theta * b * db1
			dbtt := 2.0 * b * db1 + 4.0 * theta * theta * b * (b - 1.0) * db2
			dbdt := -nonan_aa * b * 2.0 / nonan_beta * db1 * dm * q -
				2.0 * theta * b * (b - 1.0) * db2 * delta_d
			ft += n * delta * (dbt * psi + db * psi_t)
			ftt += n * delta * (dbtt * psi + 2.0 * dbt * psi_t + db * psi_tt)
			fdt += n * (db * (psi_t + delta * psi_dt) + delta * dbd * psi_t +
				dbt * (psi + delta * psi_d) + dbdt * delta * psi)
		}
	}
	if tt_at_critical != 0.0 {
		ftt = tt_at_critical * math.inf(1)
	}
	return Phi{
		f:  f
		d:  fd
		dd: fdd
		t:  ft
		tt: ftt
		dt: fdt
	}
}

// phi_ideal is IAPWS-95's ideal-gas part φ° at (τ, δ), Eq. (5), with its
// derivatives (Table 4). φ°_δτ is 0.
pub fn phi_ideal(tau f64, delta f64) Phi {
	mut f := math.log(delta) + ideal_n1 + ideal_n2 * tau + ideal_n3 * math.log(tau)
	mut t := ideal_n2 + ideal_n3 / tau
	mut tt := -ideal_n3 / (tau * tau)
	for i in 0 .. 5 {
		gt := ideal_gamma[i] * tau
		e := math.exp(-gt)
		om := -math.expm1(-gt) // 1 − e^(−γτ)
		f += ideal_n[i] * math.log(om)
		t += ideal_n[i] * ideal_gamma[i] * (1.0 / om - 1.0)
		tt -= ideal_n[i] * ideal_gamma[i] * ideal_gamma[i] * e / (om * om)
	}
	return Phi{
		f:  f
		d:  1.0 / delta
		dd: -1.0 / (delta * delta)
		t:  t
		tt: tt
	}
}

// phi_residual is φʳ at (τ, δ) with all six derivatives, for tests and the
// steam-table functions: the model's own evaluation.
pub fn phi_residual(tau f64, delta f64) Phi {
	return phi_r(tau_parts(tau), delta, true)
}

// ideal_cp is IAPWS-95's ideal-gas Cp as reflux's ideal part (D5): a
// thermo.PlanckEinstein with a = R + R₉₅·n₃°, cᵢ = R₉₅·nᵢ° and θᵢ = γᵢ°·Tc.
// R − R₉₅ in a is D4's: with it the ideal gas's Cv is IAPWS-95's, and the
// Helmholtz energy is the release's but for a term linear in T.
pub fn ideal_cp() thermo.CpForm {
	return thermo.PlanckEinstein{
		a:     thermo.r_gas + r95 * ideal_n3
		terms: []thermo.EinsteinTerm{len: 5, init: thermo.EinsteinTerm{r95 * ideal_n[index], ideal_gamma[index] * critical_t}}
	}
}

// fluid is IAPWS-95 water as a thermo.Fluid, its H and S on reflux's
// reference through water's ΔH°f (J/mol) and S° (J/(mol·K)) at 298.15 K
// and 1 bar, which the caller takes from its data (D5).
pub fn fluid(dhf ?f64, s0 ?f64) !thermo.Fluid {
	return thermo.new_fluid(iapws95(), [
		thermo.IdealGas{
			cp:  ideal_cp()
			dhf: dhf
			s0:  s0
		},
	])
}
