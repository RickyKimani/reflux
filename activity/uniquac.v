module activity

import math
import numeric

// Uniquac is Abrams and Prausnitz's UNIQUAC (1975), in DECHEMA's form, q in
// both terms:
//
//     ln γᵢ = ln(rᵢ/R̄) + 5qᵢ·ln((qᵢ/Q̄)/(rᵢ/R̄)) + lᵢ − (rᵢ/R̄)·Σⱼ xⱼlⱼ
//             + qᵢ·(1 − ln Sᵢ − Σⱼ θⱼτᵢⱼ/Sⱼ),
//     R̄ = Σ xⱼrⱼ,   Q̄ = Σ xⱼqⱼ,   θⱼ = xⱼqⱼ/Q̄,   Sᵢ = Σⱼ θⱼτⱼᵢ,
//     lᵢ = 5(rᵢ − qᵢ) − (rᵢ − 1),
//
// the coordination number being 10, and τᵢᵢ = 1. ln(rᵢ/R̄) is ln(Φᵢ/xᵢ),
// written so that an absent component's xᵢ = 0 divides nothing. Each τᵢⱼ is
// an expression in T (D7); DECHEMA's is exp(−Aᵢⱼ/RT).
pub struct Uniquac {
pub:
	names []string
	r     []f64
	q     []f64
	tau   []Param // nc×nc, row-major; the diagonal unused
mut:
	l []f64
}

// UniquacPair is one pair's parameters, i before j: τᵢⱼ and τⱼᵢ.
pub struct UniquacPair {
pub:
	i      int
	j      int
	tau_ij Param
	tau_ji Param
}

// new_uniquac is UNIQUAC for components named names, with their volume and
// area parameters r and q, every pair given once.
pub fn new_uniquac(names []string, r []f64, q []f64, pairs []UniquacPair) !Uniquac {
	nc := names.len
	if r.len != nc || q.len != nc {
		return pair_error('UNIQUAC', '${r.len} r and ${q.len} q for ${nc} components')
	}
	for i in 0 .. nc {
		if !(r[i] > 0.0) || !(q[i] > 0.0) || !math.is_finite(r[i]) || !math.is_finite(q[i]) {
			return pair_error('UNIQUAC', '${names[i]} has r = ${r[i]} and q = ${q[i]}')
		}
	}
	one := param('1', expr_definitions(), 'τ')!
	mut tau := []Param{len: nc * nc, init: one}
	mut given := []bool{len: nc * nc}
	for p in pairs {
		check_pair('UNIQUAC', names, p.i, p.j, mut given)!
		tau[p.i * nc + p.j] = p.tau_ij
		tau[p.j * nc + p.i] = p.tau_ji
	}
	every_pair('UNIQUAC', names, given)!
	return Uniquac{
		names: names.clone()
		r:     r.clone()
		q:     q.clone()
		tau:   tau
		l:     []f64{len: nc, init: 5.0 * (r[index] - q[index]) - (r[index] - 1.0)}
	}
}

pub fn (m &Uniquac) size() int {
	return m.names.len
}

pub fn (m &Uniquac) name() string {
	return 'UNIQUAC'
}

// work holds τ, nc², and θ and S.
pub fn (m &Uniquac) work() Work {
	nc := m.names.len
	return new_work(nc, nc * nc, 2 * nc)
}

pub fn (m &Uniquac) params(t f64, mut w Work) ! {
	nc := m.names.len
	for i in 0 .. nc {
		for j in 0 .. nc {
			k := i * nc + j
			w.pj[k] = if i == j {
				numeric.constant(1.0)
			} else {
				finite_param(m, &m.tau[k], t, 'τ of ${m.names[i]} with ${m.names[j]}')!
			}
			if !usable(w.pj[k]) {
				return unusable(m, 'τ of ${m.names[i]} with ${m.names[j]}', w.pj[k], t)
			}
			w.pv[k] = w.pj[k].v
		}
	}
}

pub fn (m &Uniquac) lng(p []f64, x []f64, mut s []f64, mut out []f64) {
	uniquac[f64](m.r, m.q, m.l, p, x, mut s, mut out)
}

pub fn (m &Uniquac) lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet) {
	uniquac[numeric.Jet](m.r, m.q, m.l, p, x, mut s, mut out)
}

// uniquac is UNIQUAC's ln γ: p holds τ, s takes θ then S.
fn uniquac[N](r []f64, q []f64, l []f64, p []N, x []N, mut s []N, mut out []N) {
	nc := r.len
	mut rbar := lift[N](0.0)
	mut qbar := lift[N](0.0)
	mut xl := lift[N](0.0)
	for j in 0 .. nc {
		rbar = rbar + x[j] * lift[N](r[j])
		qbar = qbar + x[j] * lift[N](q[j])
		xl = xl + x[j] * lift[N](l[j])
	}
	for j in 0 .. nc {
		s[j] = x[j] * lift[N](q[j]) / qbar
	}
	for i in 0 .. nc {
		mut si := lift[N](0.0)
		for j in 0 .. nc {
			si = si + s[j] * p[j * nc + i]
		}
		s[nc + i] = si
	}
	for i in 0 .. nc {
		ri := lift[N](r[i])
		qi := lift[N](q[i])
		v := ri / rbar
		f := qi / qbar
		mut res := lift[N](0.0)
		for j in 0 .. nc {
			res = res + s[j] * p[i * nc + j] / s[nc + j]
		}
		out[i] = nlog[N](v) + lift[N](5.0) * qi * nlog[N](f / v) + lift[N](l[i]) - v * xl +
			qi * (lift[N](1.0) - nlog[N](s[nc + i]) - res)
	}
}
