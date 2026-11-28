module water

import math
import transport
import water.if97

// Water's transport properties (M5a's D11): the viscosity by IAPWS R12-08,
// the thermal conductivity by R15-11 and the surface tension by R1-76(2014),
// each transcribed from its release with its tables cited. The scientific
// forms take IAPWS-95's derivatives for their critical enhancements; the
// industrial forms take IF97's states. docs/design.md, "Water's transport",
// has the equations.

// TransportError is a state outside a correlation's range, saying which.
pub struct TransportError {
	Error
pub:
	what string
}

pub fn (e TransportError) msg() string {
	return e.what
}

// The reducing constants (R12-08, Eqs. 1 to 4; R15-11, Eqs. 1 to 6).
const t_star = 647.096 // K
const rho_star = 322.0 // kg/m³
const p_star = 22.064e6 // Pa
const r_spec = 461.51805 // J/(kg·K), IAPWS-95's, R15-11's Eq. (6)

// R12-08's Table 1: μ₀'s Hᵢ.
const visc_h = [1.67752, 2.20462, 0.6366564, -0.241605]

// R12-08's Table 2: μ₁'s Hᵢⱼ, those not zero, as (i, j, Hᵢⱼ).
const visc_hij_i = [0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5]
const visc_hij_j = [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6]
const visc_hij = [5.20094e-1, 8.50895e-2, -1.08374, -2.89555e-1, 2.22531e-1, 9.99115e-1, 1.88797,
	1.26613, 1.20573e-1, -2.81378e-1, -9.06851e-1, -7.72479e-1, -4.89837e-1, -2.57040e-1, 1.61913e-1,
	2.57399e-1, -3.25372e-2, 6.98452e-2, 8.72102e-3, -4.35673e-3, -5.93264e-4]

// R12-08's Table 3 and R15-11's Table 3: the critical region's constants.
const x_mu = 0.068
const visc_qc = 1.0 / 1.9 // nm⁻¹
const visc_qd = 1.0 / 1.1 // nm⁻¹
const crit_nu = 0.630
const crit_gamma = 1.239
const crit_xi0 = 0.13 // nm
const crit_gamma0 = 0.06
const crit_tr = 1.5
const cond_lambda = 177.8514
const cond_qd = 1.0 / 0.40 // nm⁻¹

// R15-11's Table 1: λ₀'s Lₖ.
const cond_l = [2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4]

// R15-11's Table 2: λ₁'s Lᵢⱼ, row i, column j.
const cond_lij = [
	[1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258],
	[2.33771842, -2.78843778, 1.53616167, -0.463045512, 0.0832827019, -0.00719201245],
	[2.19650529, -4.54580785, 3.55777244, -1.40944978, 0.275418278, -0.0205938816],
	[-1.21051378, 1.60812989, -0.621178141, 0.0716373224, 0.0, 0.0],
	[-2.7203370, 4.57586331, -3.18369245, 1.1168348, -0.19268305, 0.012913842],
]

// R15-11's Table 6 and Eq. (26): the industrial ζ(T_R, ρ)'s Aᵢⱼ, a column
// per range of ρ̄, each column's six Aᵢ from i = 0.
const cond_a_rho = [0.310559006, 0.776397516, 1.242236025, 1.863354037]
const cond_a = [
	[6.53786807199516, -5.61149954923348, 3.39624167361325, -2.27492629730878, 10.2631854662709,
		1.97815050331519],
	[6.52717759281799, -6.30816983387575, 8.08379285492595, -9.82240510197603, 12.1358413791395,
		-5.54349664571295],
	[5.35500529896124, -3.96415689925446, 8.91990208918795, -12.0338729505790, 9.19494865194302,
		-2.16866274479712],
	[1.55225959906681, 0.464621290821181, 8.93237374861479, -11.0321960061126, 6.16780999933360,
		-0.965458722086812],
	[1.11999926419994, 0.595748562571649, 9.88952565078920, -10.3255051147040, 4.66861294457414,
		-0.503243546373828],
]

// R1-76(2014): σ = B τ^μ (1 + bτ), τ = 1 − T/Tc.
const sigma_b = 235.8e-3 // N/m
const sigma_bb = -0.625
const sigma_mu = 1.256

// The ranges (D11, as read in step 6): Eqs. (9) of R12-08 and (14) of
// R15-11 run from the melting curve, which M5a doesn't have, to 1173.15 K,
// and both releases endorse their extrapolation to the rest of IAPWS-95's
// range (to 1273.15 K), and from 250 K for the vapour and the subcooled
// liquid. Above 1273.15 K they endorse only the dilute gas, the viscosity's
// to 2500 K and the conductivity's to about 1500 K; a gas is read as dilute
// to ρ = 0.01ρc, where μ₁ and λ₁ are within about 1% of 1.
const transport_t_min = 250.0
const transport_t_dense = 1273.15
const visc_t_dilute = 2500.0
const cond_t_dilute = 1500.0
const dilute_rho = 0.01 * rho_star

fn check_range(t f64, rho f64, t_dilute f64, what string) ! {
	if !(rho >= 0.0) || !(t >= transport_t_min) {
		return TransportError{
			what: 'IAPWS gives ${what} from 250 K, not ${transport.shown(t)} K and ${transport.shown(rho)} kg/m³'
		}
	}
	if t > transport_t_dense && !(t <= t_dilute && rho <= dilute_rho) {
		return TransportError{
			what: 'IAPWS gives ${what} above 1273.15 K only for a dilute gas (to ${dilute_rho} kg/m³) up to ${t_dilute} K, not ${transport.shown(t)} K and ${transport.shown(rho)} kg/m³'
		}
	}
}

// mu0 and mu1 are R12-08's Eqs. (11) and (12), dimensionless (μ* = 1 µPa·s).
fn mu0(tr f64) f64 {
	mut s := 0.0
	for i, h in visc_h {
		s += h / math.pow(tr, f64(i))
	}
	return 100.0 * math.sqrt(tr) / s
}

fn mu1(tr f64, dr f64) f64 {
	a := 1.0 / tr - 1.0
	b := dr - 1.0
	mut s := 0.0
	for k, h in visc_hij {
		s += h * math.pow(a, f64(visc_hij_i[k])) * math.pow(b, f64(visc_hij_j[k]))
	}
	return math.exp(dr * s)
}

// zeta95 is ζ = (∂ρ̄/∂p̄)_T by IAPWS-95 (R12-08's Eq. 21a).
fn zeta95(t f64, rho f64) f64 {
	delta := rho / critical_rho_mass
	r := phi_residual(critical_t / t, delta)
	dpdrho := r_spec * t * (1.0 + 2.0 * delta * r.d + delta * delta * r.dd)
	return p_star / (rho_star * dpdrho)
}

// xi is the correlation length in nm (R12-08's Eqs. 20 and 21, R15-11's
// 22 and 23) from ζ at T and at T_R·T*, Δχ̄ set to 0 where it's negative.
fn xi(tr f64, dr f64, zeta f64, zeta_r f64) f64 {
	dchi := dr * (zeta - zeta_r * crit_tr / tr)
	if !(dchi > 0.0) {
		return 0.0
	}
	return crit_xi0 * math.pow(dchi / crit_gamma0, crit_nu / crit_gamma)
}

// mu2 is R12-08's critical enhancement, Eqs. (14) to (19), from ξ (nm).
fn mu2(x f64) f64 {
	qc := visc_qc * x
	qd := visc_qd * x
	if x <= 0.3817016416 {
		y := 0.2 * qc * math.pow(qd, 5.0) * (1.0 - qc + qc * qc - 765.0 / 504.0 * qd * qd)
		return math.exp(x_mu * y)
	}
	psi := math.acos(1.0 / math.sqrt(1.0 + qd * qd))
	w := math.sqrt(math.abs((qc - 1.0) / (qc + 1.0))) * math.tan(0.5 * psi)
	l := if qc > 1.0 { math.log((1.0 + w) / (1.0 - w)) } else { 2.0 * math.atan(math.abs(w)) }
	y := math.sin(3.0 * psi) / 12.0 - math.sin(2.0 * psi) / (4.0 * qc) +
		(1.0 - 1.25 * qc * qc) * math.sin(psi) / (qc * qc) - ((1.0 - 1.5 * qc * qc) * psi - math.pow(math.abs(qc * qc - 1.0), 1.5) * l) / (qc * qc * qc)
	return math.exp(x_mu * y)
}

// viscosity is R12-08's μ at (T, ρ), Pa·s, its critical enhancement from
// IAPWS-95 (Eq. 10).
pub fn viscosity(t f64, rho f64) !f64 {
	check_range(t, rho, visc_t_dilute, 'the viscosity')!
	return viscosity_parts(t, rho).mu
}

// Viscosity is R12-08's μ in its factors, and the correlation length.
pub struct Viscosity {
pub:
	mu  f64 // Pa·s
	mu0 f64 // dimensionless
	mu1 f64
	mu2 f64
	xi  f64 // nm
}

// viscosity_parts is R12-08's μ at (T, ρ) and its parts, with no range
// check.
pub fn viscosity_parts(t f64, rho f64) Viscosity {
	tr := t / t_star
	dr := rho / rho_star
	// At ρ = 0, the dilute gas, IAPWS-95's ζ diverges and Δχ̄ is 0.
	x := if rho == 0.0 { 0.0 } else { xi(tr, dr, zeta95(t, rho), zeta95(crit_tr * t_star, rho)) }
	m0, m1, m2 := mu0(tr), mu1(tr, dr), mu2(x)
	return Viscosity{
		mu:  m0 * m1 * m2 * 1e-6
		mu0: m0
		mu1: m1
		mu2: m2
		xi:  x
	}
}

// viscosity_industrial is R12-08's industrial μ, Pa·s (its Section 3): μ₂
// = 1, at a density IF97 gives.
pub fn viscosity_industrial(t f64, rho f64) !f64 {
	check_range(t, rho, visc_t_dilute, 'the viscosity')!
	tr := t / t_star
	return mu0(tr) * mu1(tr, rho / rho_star) * 1e-6
}

// lambda0 and lambda1 are R15-11's Eqs. (16) and (17), dimensionless (λ* =
// 1 mW/(m·K)).
fn lambda0(tr f64) f64 {
	mut s := 0.0
	for k, l in cond_l {
		s += l / math.pow(tr, f64(k))
	}
	return math.sqrt(tr) / s
}

fn lambda1(tr f64, dr f64) f64 {
	a := 1.0 / tr - 1.0
	b := dr - 1.0
	mut s := 0.0
	for i, row in cond_lij {
		mut inner := 0.0
		for j, l in row {
			inner += l * math.pow(b, f64(j))
		}
		s += math.pow(a, f64(i)) * inner
	}
	return math.exp(dr * s)
}

// lambda2 is R15-11's critical enhancement, Eqs. (18) to (21),
// dimensionless: cp in J/(kg·K), κ = cp/cv, μ̄ the dimensionless viscosity
// and ξ in nm.
fn lambda2(tr f64, dr f64, cp f64, kappa f64, mubar f64, x f64) f64 {
	return cond_lambda * dr * (cp / r_spec) * tr / mubar * z_of(cond_qd * x, kappa, dr)
}

// z_of is R15-11's Z(y), Eqs. (19) and (21).
fn z_of(y f64, kappa f64, dr f64) f64 {
	if y < 1.2e-7 {
		return 0.0
	}
	return 2.0 / (math.pi * y) * (((1.0 - 1.0 / kappa) * math.atan(y) + y / kappa) - (1.0 - math.exp(-1.0 / (1.0 / y +
		y * y / (3.0 * dr * dr)))))
}

// Conductivity is R15-11's λ in its parts.
pub struct Conductivity {
pub:
	lambda  f64 // W/(m·K)
	lambda0 f64 // dimensionless
	lambda1 f64
	lambda2 f64
	xi      f64 // nm
}

// thermal_conductivity is R15-11's λ at (T, ρ), W/(m·K), Eq. (15), its
// critical enhancement from IAPWS-95's cp, cv and ζ and R12-08's μ.
pub fn thermal_conductivity(t f64, rho f64) !f64 {
	check_range(t, rho, cond_t_dilute, 'the thermal conductivity')!
	return conductivity_parts(t, rho).lambda
}

// conductivity_parts is R15-11's λ at (T, ρ) and its parts, with no range
// check.
pub fn conductivity_parts(t f64, rho f64) Conductivity {
	tr := t / t_star
	dr := rho / rho_star
	l0, l1 := lambda0(tr), lambda1(tr, dr)
	v := viscosity_parts(t, rho)
	// At ρ = 0 IAPWS-95's derivatives diverge, and λ₂ is 0 (R15-11, the
	// note to its Table 4).
	l2 := if rho == 0.0 {
		0.0
	} else {
		pr := props95(t, rho / molar_mass)
		lambda2(tr, dr, pr.cp / molar_mass, pr.cp / pr.cv, v.mu * 1e6, v.xi)
	}
	return Conductivity{
		lambda:  (l0 * l1 + l2) * 1e-3
		lambda0: l0
		lambda1: l1
		lambda2: l2
		xi:      v.xi
	}
}

// zeta_tr_industrial is R15-11's Eq. (25), the industrial ζ(T_R, ρ).
fn zeta_tr_industrial(dr f64) f64 {
	mut j := cond_a_rho.len
	for k, edge in cond_a_rho {
		if dr <= edge {
			j = k
			break
		}
	}
	mut s := 0.0
	for i, a in cond_a[j] {
		s += a * math.pow(dr, f64(i))
	}
	return 1.0 / s
}

// conductivity_industrial is R15-11's industrial λ (its Section 3) at an
// IF97 state, and its parts: cp, cv and (∂ρ/∂p)_T from IF97, ζ(T_R) by Eq.
// (25), μ by R12-08's industrial form; ζ and cp taken as 1e13 where IF97
// gives them negative or above (its footnote 2), and λ₂ = 0 in region 5
// (footnote 3).
pub fn conductivity_industrial(st if97.State) !Conductivity {
	if st.region == 4 && !(st.x == 0.0 || st.x == 1.0) {
		return TransportError{
			what: 'a two-phase mixture has no thermal conductivity'
		}
	}
	rho := 1.0 / st.v
	check_range(st.t, rho, cond_t_dilute, 'the thermal conductivity')!
	tr := st.t / t_star
	dr := rho / rho_star
	l0, l1 := lambda0(tr), lambda1(tr, dr)
	mut l2 := 0.0
	mut x := 0.0
	if st.region != 5 {
		mut zeta := p_star / rho_star * st.drhodp
		if !(zeta >= 0.0 && zeta <= 1e13) {
			zeta = 1e13
		}
		mut cp := st.cp
		if !(cp >= 0.0 && cp <= 1e13) {
			cp = 1e13
		}
		x = xi(tr, dr, zeta, zeta_tr_industrial(dr))
		mubar := mu0(tr) * mu1(tr, dr)
		l2 = lambda2(tr, dr, cp, cp / st.cv, mubar, x)
	}
	return Conductivity{
		lambda:  (l0 * l1 + l2) * 1e-3
		lambda0: l0
		lambda1: l1
		lambda2: l2
		xi:      x
	}
}

// surface_tension is R1-76's σ, N/m, from 248.15 K (−25 °C, the release's
// endorsed extrapolation into the supercooled liquid) to Tc.
pub fn surface_tension(t f64) !f64 {
	if !(t >= 248.15 && t <= critical_t) {
		return TransportError{
			what: 'IAPWS gives the surface tension from 248.15 K to ${critical_t} K, not ${transport.shown(t)} K'
		}
	}
	tau := 1.0 - t / critical_t
	return sigma_b * math.pow(tau, sigma_mu) * (1.0 + sigma_bb * tau)
}
