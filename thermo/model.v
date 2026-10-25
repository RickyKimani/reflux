module thermo

import math
import numeric

// Model is anything that can supply the reduced residual Helmholtz energy
//
//     F(T, V, n) = Aʳ(T, V, n) / RT
//
// and its derivatives, where V is the total volume (m³) and n the mole
// numbers (mol). This follows Michelsen & Mollerup, "Thermodynamic Models:
// Fundamentals and Computational Aspects". Pressure, Z, ln φᵢ, Hʳ and Sʳ
// all follow from F, so a model only has to provide F. Cubics, CPA (a cubic
// plus an association term), PC-SAFT and multiparameter reference
// equations all fit, and none of the code below knows which one it has.
//
// This is an interface rather than a sum type because the set of models is
// open: adding one should never require editing this module.
pub interface Model {
	// components lists the components, in the order used by n.
	components() []Component
	// min_volume is the smallest physically allowed V for composition n:
	// the covolume Σnᵢbᵢ for a cubic, close packing for SAFT. Every volume
	// root lies above it.
	min_volume(n []f64) f64
	// volume_derivs returns F and its first two V-derivatives. The volume
	// solver calls it on every iteration, so it should be cheap.
	volume_derivs(t f64, v f64, n []f64) VolumeDerivs
	// residual returns F with its derivatives in V, T and n.
	residual(t f64, v f64, n []f64) Residual
}

// VolumeDerivs is the part of the residual the volume solver needs. Since
// F is in mol, F_V is in mol/m³ and F_VV in mol/m⁶.
pub struct VolumeDerivs {
pub:
	f    f64
	f_v  f64
	f_vv f64
}

// Residual is F with the derivatives needed for P, ln φᵢ, Hʳ and Sʳ.
// Second derivatives in T and n (for Cp and for Newton flash) come later.
pub struct Residual {
pub:
	f    f64   // F, mol
	f_v  f64   // ∂F/∂V, mol/m³
	f_vv f64   // ∂²F/∂V², mol/m⁶
	f_t  f64   // ∂F/∂T, mol/K
	f_n  []f64 // ∂F/∂nᵢ, dimensionless
}

// Root says which volume root to return.
pub enum Root {
	liquid // smallest root
	vapor  // largest root
	stable // the root with the lower Gibbs energy
}

// State is one phase at a given T, P and n.
pub struct State {
pub:
	t      f64   // K
	p      f64   // Pa
	n      []f64 // mol
	v      f64   // total volume, m³
	z      f64   // compressibility factor
	ln_phi []f64 // ln φᵢ
}

// molar_volume returns V/n in m³/mol.
pub fn (s State) molar_volume() f64 {
	return s.v / numeric.sum(s.n)
}

// pressure evaluates P(T, V, n) = RT(n/V − ∂F/∂V), in Pa.
pub fn pressure(m Model, t f64, v f64, n []f64) f64 {
	d := m.volume_derivs(t, v, n)
	return r_gas * t * (numeric.sum(n) / v - d.f_v)
}

// solve_tp finds the requested volume root at (T, P, n) and the fugacity
// coefficients there. If only one root exists it is returned whatever
// `root` asks for, so asking for `.liquid` where only vapor exists gives
// you the vapor root.
pub fn solve_tp(m Model, t f64, p f64, n []f64, root Root) !State {
	check_tpn(m, t, p, n)!
	v := match root {
		.liquid {
			liquid_volume(m, t, p, n) or { vapor_volume(m, t, p, n) or { return no_root(t, p) } }
		}
		.vapor {
			vapor_volume(m, t, p, n) or { liquid_volume(m, t, p, n) or { return no_root(t, p) } }
		}
		.stable {
			stable_volume(m, t, p, n)!
		}
	}
	return state_at(m, t, p, v, n)
}

fn state_at(m Model, t f64, p f64, v f64, n []f64) State {
	r := m.residual(t, v, n)
	z := p * v / (numeric.sum(n) * r_gas * t)
	ln_z := math.log(z)
	return State{
		t:      t
		p:      p
		n:      n.clone()
		v:      v
		z:      z
		ln_phi: r.f_n.map(it - ln_z) // ln φᵢ = ∂F/∂nᵢ − ln Z
	}
}

fn stable_volume(m Model, t f64, p f64, n []f64) !f64 {
	vl := liquid_volume(m, t, p, n) or {
		return vapor_volume(m, t, p, n) or { return no_root(t, p) }
	}
	vv := vapor_volume(m, t, p, n) or { return vl }
	return if gibbs_residual(m, t, p, vl, n) < gibbs_residual(m, t, p, vv, n) { vl } else { vv }
}

// gibbs_residual returns Gʳ/RT = F + n(Z − 1) − n·ln Z at a volume root.
// For a pure component this is ln φ, which is why the root with the lower
// value is the stable one.
fn gibbs_residual(m Model, t f64, p f64, v f64, n []f64) f64 {
	ntot := numeric.sum(n)
	z := p * v / (ntot * r_gas * t)
	return m.volume_derivs(t, v, n).f + ntot * (z - 1.0) - ntot * math.log(z)
}

fn no_root(t f64, p f64) IError {
	return error('no volume root at T = ${t} K, P = ${p} Pa')
}

// The volume solvers below use plain Newton steps, started on the side of
// the root where the iteration cannot overshoot, and limited in size so
// that they cannot jump across the unstable (van der Waals loop) region
// onto the other branch. Each returns `none` when its branch does not
// reach the requested pressure. `?f64` is an Option: a value or `none`,
// with no error message, because "no liquid root here" is an answer,
// not a failure.

const volume_max_iter = 100
const volume_rtol = 1e-12

// liquid_volume finds the smallest root of P(T, V, n) = p. It starts just
// above the minimum volume, where P is huge, and walks right. P(V) is convex
// on the liquid branch, so Newton approaches the root from the left without
// overshooting.
fn liquid_volume(m Model, t f64, p f64, n []f64) ?f64 {
	ntot := numeric.sum(n)
	vmin := m.min_volume(n)
	mut v := 1.05 * vmin
	for _ in 0 .. 60 {
		if p_dpdv(m, t, v, n, ntot).p > p {
			break
		}
		v = vmin + 0.5 * (v - vmin)
	}
	for _ in 0 .. volume_max_iter {
		s := p_dpdv(m, t, v, n, ntot)
		if s.dpdv >= 0.0 {
			return none // past the liquid spinodal: the liquid branch never reaches p
		}
		mut dv := -(s.p - p) / s.dpdv
		dv = math.clamp(dv, -0.5 * (v - vmin), 0.5 * (v - vmin))
		v += dv
		if math.abs(dv) <= volume_rtol * v {
			return v
		}
	}
	return none
}

// vapor_volume finds the largest root. It works in density ρ = n/V and
// starts from the ideal gas. P(ρ) is concave on the vapor branch, so Newton
// approaches the root from the dilute side without overshooting. At high
// pressure the ideal-gas density can exceed the packing limit n/B, so the
// start is capped at half of it.
fn vapor_volume(m Model, t f64, p f64, n []f64) ?f64 {
	ntot := numeric.sum(n)
	rho_max := ntot / m.min_volume(n)
	mut rho := math.min(p / (r_gas * t), 0.5 * rho_max)
	for _ in 0 .. volume_max_iter {
		v := ntot / rho
		s := p_dpdv(m, t, v, n, ntot)
		dp_drho := -s.dpdv * v * v / ntot
		if dp_drho <= 0.0 {
			return none // past the vapor spinodal
		}
		mut drho := -(s.p - p) / dp_drho
		drho = math.clamp(drho, -0.5 * rho, math.min(0.5 * rho, 0.5 * (rho_max - rho)))
		rho += drho
		if math.abs(drho) <= volume_rtol * rho {
			if rho_max - rho <= 1e-9 * rho_max {
				return none // crept up to the packing limit: not a root
			}
			return ntot / rho
		}
	}
	return none
}

struct PDpdv {
	p    f64
	dpdv f64
}

fn p_dpdv(m Model, t f64, v f64, n []f64, ntot f64) PDpdv {
	d := m.volume_derivs(t, v, n)
	rt := r_gas * t
	return PDpdv{
		p:    rt * (ntot / v - d.f_v)
		dpdv: -rt * (d.f_vv + ntot / (v * v))
	}
}

fn check_tpn(m Model, t f64, p f64, n []f64) ! {
	if t <= 0.0 {
		return InputError{
			what:  'temperature (K)'
			value: t
		}
	}
	if p <= 0.0 {
		return InputError{
			what:  'pressure (Pa)'
			value: p
		}
	}
	nc := m.components().len
	if n.len != nc {
		return error('composition has ${n.len} entries but the model has ${nc} components')
	}
	if n.any(it < 0.0) || numeric.sum(n) <= 0.0 {
		return error('mole numbers must be non-negative with a positive total, got ${n}')
	}
}
