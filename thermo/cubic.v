module thermo

import math
import numeric

// CubicKind picks the two constants δ₁, δ₂ of the generalized cubic. The
// Ωa, Ωb that place its critical point at (Tc, Pc) follow from them.
pub enum CubicKind {
	vdw
	rk
	srk
	pr
}

struct CubicConsts {
	delta1  f64
	delta2  f64
	omega_a f64
	omega_b f64
	zc      f64 // the equation's own critical compressibility
}

fn (k CubicKind) consts() CubicConsts {
	d1, d2 := match k {
		.vdw { 0.0, 0.0 }
		.rk, .srk { 1.0, 0.0 }
		.pr { 1.0 + math.sqrt2, 1.0 - math.sqrt2 }
	}
	// These four δ pairs always have a solution, so a failure here is a bug
	// in reflux, not bad input. panic is for bugs; `!` is for things callers
	// can get wrong.
	return critical_constants(d1, d2) or { panic(err) }
}

// critical_constants finds Ωa, Ωb and Zc for a cubic with constants δ₁, δ₂.
//
// At Tr = Pr = 1 (so α = 1, A = Ωa, B = Ωb), the critical point is an
// inflection with zero slope, which means the cubic in Z must be exactly
// (Z − Zc)³. Matching the coefficients of
//
//     Z³ + [(δ₁+δ₂−1)B − 1]·Z² + [δ₁δ₂B² − (δ₁+δ₂)B(B+1) + A]·Z
//        − [δ₁δ₂B²(B+1) + AB]
//
// with Z³ − 3Zc·Z² + 3Zc²·Z − Zc³ gives Zc and A in terms of B, and leaves
// one equation in B. Computing the constants this way, instead of using
// rounded textbook values, puts each equation's critical point exactly at
// (Tc, Pc).
fn critical_constants(delta1 f64, delta2 f64) !CubicConsts {
	s := delta1 + delta2
	m := delta1 * delta2
	resid := fn [s, m] (b f64) f64 {
		zc := (1.0 - (s - 1.0) * b) / 3.0
		a := 3.0 * zc * zc - m * b * b + s * b * (b + 1.0)
		return m * b * b * (b + 1.0) + a * b - zc * zc * zc
	}
	b := numeric.bisect(resid, 1e-6, 0.3)!
	zc := (1.0 - (s - 1.0) * b) / 3.0
	return CubicConsts{
		delta1:  delta1
		delta2:  delta2
		omega_a: 3.0 * zc * zc - m * b * b + s * b * (b + 1.0)
		omega_b: b
		zc:      zc
	}
}

// default_alpha is the α function each equation was published with.
fn (k CubicKind) default_alpha(omega f64) Alpha {
	return match k {
		.vdw { Unity{} }
		.rk { InvSqrtTr{} }
		.srk { soave_srk(omega) }
		.pr { soave_pr(omega) }
	}
}

// Cubic is a cubic equation of state for a mixture,
//
//     P = nRT/(V − B) − A(T) / ((V + δ₁B)(V + δ₂B))
//
// with van der Waals one-fluid mixing, B = Σ nᵢbᵢ and
// A = ΣΣ nᵢnⱼ √(aᵢaⱼ)(1 − kᵢⱼ). Following Michelsen & Mollerup,
//
//     F = −n·g(V, B) − A(T)/(RT) · f(V, B)
//     g = ln(1 − B/V)
//     f = ln((V + δ₁B)/(V + δ₂B)) / (B(δ₁ − δ₂))
//
// Cubic is one Model among many. Nothing outside this file knows that it
// is cubic, apart from the z_roots shortcut.
pub struct Cubic {
pub:
	kind   CubicKind
	comps  []Component
	alphas []Alpha // one per component
	kij    []f64   // nc×nc, row-major, symmetric
	ac     []f64   // aᵢ at Tr = 1, Pa·m⁶/mol²
	b      []f64   // bᵢ, m³/mol
	delta1 f64
	delta2 f64
}

// CubicOptions holds optional arguments to new_cubic. @[params] lets the
// caller pass them by name, or leave them all out:
//
//     thermo.new_cubic(.pr, comps)!
//     thermo.new_cubic(.pr, comps, kij: [0.0, 0.02, 0.02, 0.0])!
@[params]
pub struct CubicOptions {
pub:
	// kij is the nc×nc binary interaction matrix, row-major. Empty means zero.
	kij []f64
	// alphas sets each component's α. Empty means the equation's own.
	alphas []Alpha
}

// new_cubic builds a cubic model for the given components.
pub fn new_cubic(kind CubicKind, comps []Component, opts CubicOptions) !Cubic {
	nc := comps.len
	if nc == 0 {
		return error('a cubic needs at least one component')
	}
	for c in comps {
		c.check()!
	}
	kij := if opts.kij.len == 0 { []f64{len: nc * nc} } else { opts.kij.clone() }
	if kij.len != nc * nc {
		return error('kij has ${kij.len} entries, expected ${nc}×${nc} = ${nc * nc}')
	}
	for i in 0 .. nc {
		for j in 0 .. i {
			if kij[i * nc + j] != kij[j * nc + i] {
				return error('kij must be symmetric: k[${i}][${j}] = ${kij[i * nc + j]} but k[${j}][${i}] = ${kij[j * nc + i]}')
			}
		}
	}
	alphas := if opts.alphas.len == 0 {
		comps.map(kind.default_alpha(it.omega))
	} else {
		opts.alphas.clone()
	}
	if alphas.len != nc {
		return error('got ${alphas.len} alpha functions for ${nc} components')
	}
	k := kind.consts()
	return Cubic{
		kind:   kind
		comps:  comps.clone()
		alphas: alphas
		kij:    kij
		ac:     comps.map(k.omega_a * r_gas * r_gas * it.tc * it.tc / it.pc)
		b:      comps.map(k.omega_b * r_gas * it.tc / it.pc)
		delta1: k.delta1
		delta2: k.delta2
	}
}

// name returns the equation's usual name.
pub fn (c Cubic) name() string {
	return match c.kind {
		.vdw { 'van der Waals' }
		.rk { 'Redlich-Kwong' }
		.srk { 'Soave-Redlich-Kwong' }
		.pr { 'Peng-Robinson' }
	}
}

// components implements Model.
pub fn (c Cubic) components() []Component {
	return c.comps
}

// min_volume implements Model: the covolume B = Σ nᵢbᵢ.
pub fn (c Cubic) min_volume(n []f64) f64 {
	return numeric.dot(n, c.b)
}

// volume_derivs implements Model.
pub fn (c Cubic) volume_derivs(t f64, v f64, n []f64) VolumeDerivs {
	ntot := numeric.sum(n)
	q := c.vb_terms(v, c.min_volume(n))
	a_rt := c.attraction(t, n).a / (r_gas * t)
	return VolumeDerivs{
		f:    -ntot * q.g - a_rt * q.f
		f_v:  -ntot * q.g_v - a_rt * q.f_v
		f_vv: -ntot * q.g_vv - a_rt * q.f_vv
	}
}

// residual implements Model.
pub fn (c Cubic) residual(t f64, v f64, n []f64) Residual {
	ntot := numeric.sum(n)
	q := c.vb_terms(v, c.min_volume(n))
	att := c.attraction(t, n)
	rt := r_gas * t
	a_rt := att.a / rt
	// ∂F/∂nᵢ = −g − n·g_B·bᵢ − f·Aᵢ/RT − (A/RT)·f_B·bᵢ
	mut f_n := []f64{len: n.len}
	for i in 0 .. n.len {
		f_n[i] = -q.g - ntot * q.g_b * c.b[i] - q.f * att.a_n[i] / rt - a_rt * q.f_b * c.b[i]
	}
	return Residual{
		f:    -ntot * q.g - a_rt * q.f
		f_v:  -ntot * q.g_v - a_rt * q.f_v
		f_vv: -ntot * q.g_vv - a_rt * q.f_vv
		f_t:  -q.f / rt * (att.a_t - att.a / t)
		f_n:  f_n
	}
}

// z_roots solves the cubic directly as a polynomial in Z and returns the
// roots with Z > B·P/(nRT), ascending. It is a shortcut only a cubic has.
// solve_tp finds the same roots by Newton iteration, which is what the
// tests check.
pub fn (c Cubic) z_roots(t f64, p f64, n []f64) []f64 {
	ntot := numeric.sum(n)
	x := n.map(it / ntot)
	rt := r_gas * t
	big_a := c.attraction(t, x).a * p / (rt * rt)
	big_b := c.min_volume(x) * p / rt
	s := c.delta1 + c.delta2
	m := c.delta1 * c.delta2
	c2 := (s - 1.0) * big_b - 1.0
	c1 := m * big_b * big_b - s * big_b * (big_b + 1.0) + big_a
	c0 := -(m * big_b * big_b * (big_b + 1.0) + big_a * big_b)
	return numeric.cubic_real_roots(c2, c1, c0).filter(it > big_b)
}

// VbTerms holds g and f with their derivatives at one (V, B).
struct VbTerms {
	g    f64
	g_v  f64
	g_vv f64
	g_b  f64
	f    f64
	f_v  f64
	f_vv f64
	f_b  f64
}

fn (c Cubic) vb_terms(v f64, bb f64) VbTerms {
	d1 := v + c.delta1 * bb
	d2 := v + c.delta2 * bb
	// For van der Waals δ₁ = δ₂ and f takes its limit, 1/(V + δB).
	f := if c.delta1 == c.delta2 {
		1.0 / d1
	} else {
		math.log(d1 / d2) / (bb * (c.delta1 - c.delta2))
	}
	f_v := -1.0 / (d1 * d2)
	return VbTerms{
		g:    math.log(1.0 - bb / v)
		g_v:  bb / (v * (v - bb))
		g_vv: 1.0 / (v * v) - 1.0 / ((v - bb) * (v - bb))
		g_b:  -1.0 / (v - bb)
		f:    f
		f_v:  f_v
		f_vv: (d1 + d2) / (d1 * d1 * d2 * d2)
		// f is homogeneous of degree −1 in (V, B), so V·f_V + B·f_B = −f.
		f_b:  -(f + v * f_v) / bb
	}
}

// Attraction is A(T) for the mixture with ∂A/∂T and ∂A/∂nᵢ.
struct Attraction {
	a   f64
	a_t f64
	a_n []f64
}

fn (c Cubic) attraction(t f64, n []f64) Attraction {
	nc := c.comps.len
	// sᵢ = √aᵢ(T), so that aᵢⱼ = (1 − kᵢⱼ)·sᵢ·sⱼ.
	mut s := []f64{len: nc}
	mut s_t := []f64{len: nc}
	for i in 0 .. nc {
		tc := c.comps[i].tc
		alpha, dalpha_dtr := c.alphas[i].eval(t / tc)
		s[i] = math.sqrt(c.ac[i] * alpha)
		s_t[i] = c.ac[i] * dalpha_dtr / tc / (2.0 * s[i])
	}
	mut a := 0.0
	mut a_t := 0.0
	mut a_n := []f64{len: nc}
	for i in 0 .. nc {
		// qᵢ = Σⱼ (1 − kᵢⱼ)·sⱼ·nⱼ. Symmetry of kᵢⱼ gives Aᵢ = 2sᵢqᵢ and
		// ∂A/∂T = 2 Σᵢ nᵢ·∂sᵢ/∂T·qᵢ.
		mut q := 0.0
		for j in 0 .. nc {
			q += (1.0 - c.kij[i * nc + j]) * s[j] * n[j]
		}
		a += n[i] * s[i] * q
		a_t += 2.0 * n[i] * s_t[i] * q
		a_n[i] = 2.0 * s[i] * q
	}
	return Attraction{
		a:   a
		a_t: a_t
		a_n: a_n
	}
}
