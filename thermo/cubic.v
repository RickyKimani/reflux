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

// residual2 implements Model. Each line matches a row of the F table in
// docs/design.md, "Second derivatives of the cubic".
pub fn (c Cubic) residual2(t f64, v f64, n []f64) Residual2 {
	nc := n.len
	ntot := numeric.sum(n)
	q := c.vb_terms(v, c.min_volume(n))
	att := c.attraction2(t, n)
	rt := r_gas * t
	d := att.a
	d_rt := d / rt
	dd := att.a_t - d / t // D_T − D/T, which every T-derivative contains
	mut f_n := []f64{len: nc}
	mut f_vn := []f64{len: nc}
	mut f_tn := []f64{len: nc}
	mut f_nn := []f64{len: nc * nc}
	for i in 0 .. nc {
		bi := c.b[i]
		di := att.a_n[i]
		f_n[i] = -q.g - ntot * q.g_b * bi - q.f * di / rt - d_rt * q.f_b * bi
		f_vn[i] = -q.g_v - ntot * q.g_bv * bi - q.f_v * di / rt - d_rt * q.f_bv * bi
		f_tn[i] = -q.f / rt * (att.a_nt[i] - di / t) - q.f_b * bi / rt * dd
		for j in 0 .. nc {
			bj := c.b[j]
			dij := 2.0 * (1.0 - c.kij[i * nc + j]) * att.s[i] * att.s[j]
			f_nn[i * nc + j] = -q.g_b * (bi + bj) - ntot * q.g_bb * bi * bj - (dij * q.f +
				q.f_b * (di * bj + att.a_n[j] * bi) + d * q.f_bb * bi * bj) / rt
		}
	}
	return Residual2{
		Residual: Residual{
			f:    -ntot * q.g - d_rt * q.f
			f_v:  -ntot * q.g_v - d_rt * q.f_v
			f_vv: -ntot * q.g_vv - d_rt * q.f_vv
			f_t:  -q.f / rt * dd
			f_n:  f_n
		}
		f_tt:     -q.f / rt * (att.a_tt - 2.0 * att.a_t / t + 2.0 * d / (t * t))
		f_tv:     -q.f_v / rt * dd
		f_vn:     f_vn
		f_tn:     f_tn
		f_nn:     f_nn
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
	g_bv f64
	g_bb f64
	f    f64
	f_v  f64
	f_vv f64
	f_b  f64
	f_bv f64
	f_bb f64
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
	// f is homogeneous of degree −1 in (V, B) and f_B of degree −2, which
	// gives f_B and f_BB without dividing by δ₁ − δ₂.
	f_b := -(f + v * f_v) / bb
	f_bv := (c.delta1 * d2 + c.delta2 * d1) / (d1 * d1 * d2 * d2)
	vb := v - bb
	return VbTerms{
		g:    math.log(1.0 - bb / v)
		g_v:  bb / (v * vb)
		g_vv: 1.0 / (v * v) - 1.0 / (vb * vb)
		g_b:  -1.0 / vb
		g_bv: 1.0 / (vb * vb)
		g_bb: -1.0 / (vb * vb)
		f:    f
		f_v:  f_v
		f_vv: (d1 + d2) / (d1 * d1 * d2 * d2)
		f_b:  f_b
		f_bv: f_bv
		f_bb: -(2.0 * f_b + v * f_bv) / bb
	}
}

// sqrt_a returns sᵢ = √aᵢ(T) for every component, with ds/dT and d²s/dT².
// Working with s instead of a makes aᵢⱼ = (1 − kᵢⱼ)·sᵢ·sⱼ.
fn (c Cubic) sqrt_a(t f64) ([]f64, []f64, []f64) {
	nc := c.comps.len
	mut s := []f64{len: nc}
	mut s_t := []f64{len: nc}
	mut s_tt := []f64{len: nc}
	for i in 0 .. nc {
		tc := c.comps[i].tc
		alpha, d1, d2 := c.alphas[i].eval(t / tc)
		s[i] = math.sqrt(c.ac[i] * alpha)
		s_t[i] = c.ac[i] * d1 / tc / (2.0 * s[i])
		s_tt[i] = (0.5 * c.ac[i] * d2 / (tc * tc) - s_t[i] * s_t[i]) / s[i]
	}
	return s, s_t, s_tt
}

// Attraction is A(T) for the mixture with ∂A/∂T and ∂A/∂nᵢ.
struct Attraction {
	a   f64
	a_t f64
	a_n []f64
}

fn (c Cubic) attraction(t f64, n []f64) Attraction {
	nc := c.comps.len
	s, s_t, _ := c.sqrt_a(t)
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

// Attraction2 adds the second derivatives of A. It keeps s so that
// Aᵢⱼ = 2(1 − kᵢⱼ)sᵢsⱼ can be formed where it's used, instead of storing
// an nc×nc matrix.
struct Attraction2 {
	a    f64   // D
	a_t  f64   // D_T
	a_tt f64   // D_TT
	a_n  []f64 // Dᵢ
	a_nt []f64 // D_iT
	s    []f64 // √aᵢ
}

fn (c Cubic) attraction2(t f64, n []f64) Attraction2 {
	nc := c.comps.len
	s, s_t, s_tt := c.sqrt_a(t)
	mut a := 0.0
	mut a_t := 0.0
	mut a_tt := 0.0
	mut a_n := []f64{len: nc}
	mut a_nt := []f64{len: nc}
	for i in 0 .. nc {
		// qᵢ = Σⱼ (1 − kᵢⱼ)·sⱼ·nⱼ and rᵢ = Σⱼ (1 − kᵢⱼ)·s′ⱼ·nⱼ
		mut q := 0.0
		mut r := 0.0
		for j in 0 .. nc {
			w := 1.0 - c.kij[i * nc + j]
			q += w * s[j] * n[j]
			r += w * s_t[j] * n[j]
		}
		a += n[i] * s[i] * q
		a_t += 2.0 * n[i] * s_t[i] * q
		a_tt += 2.0 * n[i] * (s_tt[i] * q + s_t[i] * r)
		a_n[i] = 2.0 * s[i] * q
		a_nt[i] = 2.0 * (s_t[i] * q + s[i] * r)
	}
	return Attraction2{
		a:    a
		a_t:  a_t
		a_tt: a_tt
		a_n:  a_n
		a_nt: a_nt
		s:    s
	}
}
