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
	// volume_function fixes T and n and returns F as a function of V alone,
	// with its first two V-derivatives. The volume solver calls that function
	// on every iteration, so whatever depends only on T and n belongs here,
	// computed once: a cubic's a(T, n) and b(n), CPA's association
	// strengths, SAFT's segment diameters.
	volume_function(t f64, n []f64) VolumeFunction
	// residual returns F with its derivatives in V, T and n.
	residual(t f64, v f64, n []f64) Residual
	// residual2 returns everything residual does plus every second
	// derivative. The composition Hessian makes it O(nc²), so only callers
	// that need it (Cp, Newton flash, Jacobians) should use it.
	residual2(t f64, v f64, n []f64) Residual2
}

// Branches is what a model may know of its own isotherms' two branches, the
// vapour's from ρ → 0 and the liquid's from high density, which the volume
// solvers then keep to (M5a's D7). It's optional, as M8's flowsheet.Warm
// is: a model without it, a cubic, takes the solvers' paths unchanged.
//
// IAPWS-95 needs it. Inside its dome its isotherms have spurious loops whose
// roots are mechanically stable with Gibbs energies up to 1e19 RT below the
// real phase's, and near Tc its metastable windows are narrow enough for
// Newton to step across them onto the other branch.
pub interface Branches {
	// branches gives the limits at T and n.
	branches(t f64, n []f64) BranchLimits
}

// BranchLimits are a model's limits on its branches at one T and n.
pub struct BranchLimits {
pub:
	// split is a volume between the branches, where the isotherm has two:
	// every vapour root lies above it, and every liquid root below. Below
	// Tc a pure fluid's spinodal has its vapour side below ρc and its liquid
	// side above, so V = n/ρc is one.
	split ?f64 // m³
	// no_state says the model has no state at this T, liquid or vapour: an
	// isotherm outside its range. IAPWS-95's starts at the melting curve,
	// 251.165 K at its lowest; below it its liquid branch is an extrapolation
	// whose enthalpy falls as T rises, and its vapour is endorsed only below
	// the sublimation curve, which reflux doesn't have.
	no_state bool
}

// VolumeDerivs is the part of the residual the volume solver needs. Since
// F is in mol, F_V is in mol/m³ and F_VV in mol/m⁶.
pub struct VolumeDerivs {
pub:
	f    f64
	f_v  f64
	f_vv f64
}

// VolumeFunction is F at a fixed T and n, as a function of V. It's an
// interface, implemented by a small struct per model, rather than a struct
// of closures: V 0.5.2 never frees a closure, and the volume solver makes
// one of these for every root it finds (see docs/v-for-go.md).
pub interface VolumeFunction {
	// at is F, F_V and F_VV at V: what the volume solver needs.
	at(v f64) VolumeDerivs
	// residual is everything Model.residual gives, at V, without computing
	// the parts in T and n again: for the state at the root found.
	residual(v f64) Residual
}

// volume_derivs is F and its first two V-derivatives at one (T, V, n). To
// evaluate several volumes at the same T and n, call m.volume_function once.
pub fn volume_derivs(m Model, t f64, v f64, n []f64) VolumeDerivs {
	return m.volume_function(t, n).at(v)
}

// Residual is F with the derivatives needed for P, ln φᵢ, Hʳ and Sʳ.
pub struct Residual {
pub:
	f    f64   // F, mol
	f_v  f64   // ∂F/∂V, mol/m³
	f_vv f64   // ∂²F/∂V², mol/m⁶
	f_t  f64   // ∂F/∂T, mol/K
	f_n  []f64 // ∂F/∂nᵢ, dimensionless
}

// Residual2 is a Residual plus every second derivative of F. Embedding
// Residual works like Go's struct embedding: r.f_v, r.f_n and the rest are
// reachable directly on a Residual2.
pub struct Residual2 {
	Residual
pub:
	f_tt f64   // ∂²F/∂T², mol/K²
	f_tv f64   // ∂²F/∂T∂V, mol/(K·m³)
	f_vn []f64 // ∂²F/∂V∂nᵢ, 1/m³
	f_tn []f64 // ∂²F/∂T∂nᵢ, 1/K
	f_nn []f64 // ∂²F/∂nᵢ∂nⱼ, nc×nc row-major, 1/mol
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
	d := volume_derivs(m, t, v, n)
	return r_gas * t * (numeric.sum(n) / v - d.f_v)
}

// solve_tp finds the requested volume root at (T, P, n) and the fugacity
// coefficients there. If only one root exists it is returned whatever
// `root` asks for, so asking for `.liquid` where only vapor exists gives
// you the vapor root.
pub fn solve_tp(m Model, t f64, p f64, n []f64, root Root) !State {
	check_tpn(m, t, p, n)!
	fx := fix(m, t, n)
	v := match root {
		.liquid {
			liquid_volume(fx, p) or { vapor_volume(fx, p) or { return no_root(t, p) } }
		}
		.vapor {
			vapor_volume(fx, p) or { liquid_volume(fx, p) or { return no_root(t, p) } }
		}
		.stable {
			stable_volume(fx, p)!
		}
	}
	return state_at(fx, p, v, n)
}

fn state_at(fx Fixed, p f64, v f64, n []f64) State {
	t := fx.t
	r := fx.f.residual(v)
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

// Fixed is a model at a fixed T and n, which is what the volume solvers
// iterate on: its volume function, n's total and the smallest volume, and
// its Branches limits, 0 where it has none.
struct Fixed {
	f        VolumeFunction
	t        f64
	ntot     f64
	vmin     f64
	split    f64 // BranchLimits.split, or 0
	no_state bool
}

fn fix(m Model, t f64, n []f64) Fixed {
	f := m.volume_function(t, n)
	vmin := m.min_volume(n)
	if m is Branches {
		b := m.branches(t, n)
		return Fixed{
			f:        f
			t:        t
			ntot:     numeric.sum(n)
			vmin:     vmin
			split:    b.split or { 0.0 }
			no_state: b.no_state
		}
	}
	return Fixed{
		f:    f
		t:    t
		ntot: numeric.sum(n)
		vmin: vmin
	}
}

fn stable_volume(fx Fixed, p f64) !f64 {
	vl := liquid_volume(fx, p) or { return vapor_volume(fx, p) or { return no_root(fx.t, p) } }
	vv := vapor_volume(fx, p) or { return vl }
	return if gibbs_residual(fx, p, vl) < gibbs_residual(fx, p, vv) { vl } else { vv }
}

// gibbs_residual returns Gʳ/RT = F + n(Z − 1) − n·ln Z at a volume root.
// For a pure component this is ln φ, which is why the root with the lower
// value is the stable one.
fn gibbs_residual(fx Fixed, p f64, v f64) f64 {
	z := p * v / (fx.ntot * r_gas * fx.t)
	return fx.f.at(v).f + fx.ntot * (z - 1.0) - fx.ntot * math.log(z)
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

// f64_epsilon is the gap between 1.0 and the next f64, 2⁻⁵².
const f64_epsilon = 2.220446049250313e-16

// liquid_volume finds the smallest root of P(T, V, n) = p. It starts just
// above the minimum volume, where P is huge, and walks right. P(V) is convex
// on the liquid branch, so Newton approaches the root from the left without
// overshooting.
//
// A model whose P is finite at its minimum volume (IAPWS-95's, M5a's D7)
// may not reach p there at all: then there's no liquid root, rather than
// one at vmin, where the clamped steps would otherwise converge. And with a
// Branches split, an iterate past it has left the liquid's branch: there's
// no root on it either. A cubic reaches neither: its P → ∞ at its covolume,
// and it has no split.
fn liquid_volume(fx Fixed, p f64) ?f64 {
	if fx.no_state {
		return none
	}
	vmin := fx.vmin
	mut v := 1.05 * vmin
	mut reached := false
	for _ in 0 .. 60 {
		if p_dpdv(fx, v).p > p {
			reached = true
			break
		}
		v = vmin + 0.5 * (v - vmin)
	}
	if !reached {
		return none // the liquid's branch doesn't reach p above vmin
	}
	for _ in 0 .. volume_max_iter {
		if fx.split > 0.0 && v > fx.split {
			return none // past the split: off the liquid's branch
		}
		s := p_dpdv(fx, v)
		if s.dpdv >= 0.0 {
			return none // past the liquid spinodal: the liquid branch never reaches p
		}
		mut dv := -(s.p - p) / s.dpdv
		if math.abs(s.p - p) <= s.noise && math.abs(dv) > volume_rtol * v {
			return v // rounding has the last word (see PDpdv.noise)
		}
		dv = math.clamp(dv, -0.5 * (v - vmin), 0.5 * (v - vmin))
		v += dv
		if math.abs(dv) <= volume_rtol * v {
			if fx.split > 0.0 && v > fx.split {
				return none
			}
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
//
// With a Branches split (M5a's D7), an iterate past it has left the
// vapour's branch, and there's no root on it: IAPWS-95's ideal-gas density
// at 42 MPa and 300 K lies inside its dome, where a spurious root sits just
// above ρc. A cubic has no split.
fn vapor_volume(fx Fixed, p f64) ?f64 {
	if fx.no_state {
		return none
	}
	ntot := fx.ntot
	rho_max := ntot / fx.vmin
	mut rho := math.min(p / (r_gas * fx.t), 0.5 * rho_max)
	for _ in 0 .. volume_max_iter {
		if fx.split > 0.0 && rho > ntot / fx.split {
			return none // past the split: off the vapour's branch
		}
		v := ntot / rho
		s := p_dpdv(fx, v)
		dp_drho := -s.dpdv * v * v / ntot
		if dp_drho <= 0.0 {
			return none // past the vapor spinodal
		}
		mut drho := -(s.p - p) / dp_drho
		if math.abs(s.p - p) <= s.noise && math.abs(drho) > volume_rtol * rho {
			return v // rounding has the last word (see PDpdv.noise)
		}
		drho = math.clamp(drho, -0.5 * rho, math.min(0.5 * rho, 0.5 * (rho_max - rho)))
		rho += drho
		if math.abs(drho) <= volume_rtol * rho {
			if rho_max - rho <= 1e-9 * rho_max {
				return none // crept up to the packing limit: not a root
			}
			if fx.split > 0.0 && rho > ntot / fx.split {
				return none
			}
			return ntot / rho
		}
	}
	return none
}

struct PDpdv {
	p    f64
	dpdv f64
	// noise is how closely P can be computed here: a few rounding errors of
	// its two terms, RT·n/V and RT·F_V. Normally Newton's step test is met
	// long before. But at a critical point P(V) − p has a triple root, where
	// dP/dV and d²P/dV² vanish too, so Newton slows to linear convergence and
	// V is only determined to about ε^(1/3) relative. Close to one, dP/dV is
	// tiny and V is poorly determined too. Either way rounding takes over
	// before the step test can be met, and a residual at this level is then
	// as converged as V can get. Without this, the solver failed on half of
	// the states within 1e-6 of a pure component's Tc and 1e-8 of its Pc.
	noise f64
}

fn p_dpdv(fx Fixed, v f64) PDpdv {
	d := fx.f.at(v)
	rt := r_gas * fx.t
	ntot := fx.ntot
	return PDpdv{
		p:     rt * (ntot / v - d.f_v)
		dpdv:  -rt * (d.f_vv + ntot / (v * v))
		noise: 16.0 * f64_epsilon * rt * (ntot / v + math.abs(d.f_v))
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
