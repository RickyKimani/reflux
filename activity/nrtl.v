module activity

import numeric

// Nrtl is Renon and Prausnitz's NRTL (1968):
//
//     ln γᵢ = Aᵢ/Bᵢ + Σⱼ (xⱼGᵢⱼ/Bⱼ)·(τᵢⱼ − Aⱼ/Bⱼ),
//     Aⱼ = Σₖ xₖτₖⱼGₖⱼ,   Bⱼ = Σₖ xₖGₖⱼ,   Gᵢⱼ = exp(−αᵢⱼτᵢⱼ),
//
// with τᵢᵢ = 0. Each τᵢⱼ and αᵢⱼ is an expression in T (D7); αᵢⱼ = αⱼᵢ.
// DECHEMA's form, which ChemSep's bank gives, is τᵢⱼ = Aᵢⱼ/RT.
pub struct Nrtl {
pub:
	names []string
	tau   []Param // nc×nc, row-major; the diagonal unused
	alpha []Param // nc×nc, symmetric; the diagonal unused
}

// NrtlPair is one pair's parameters, i before j: τᵢⱼ, τⱼᵢ and α.
pub struct NrtlPair {
pub:
	i      int
	j      int
	tau_ij Param
	tau_ji Param
	alpha  Param
}

// new_nrtl is NRTL for components named names, every pair given once.
pub fn new_nrtl(names []string, pairs []NrtlPair) !Nrtl {
	nc := names.len
	zero := param('0', expr_definitions(), 'τ')!
	mut tau := []Param{len: nc * nc, init: zero}
	mut alpha := []Param{len: nc * nc, init: zero}
	mut given := []bool{len: nc * nc}
	for p in pairs {
		check_pair('NRTL', names, p.i, p.j, mut given)!
		tau[p.i * nc + p.j] = p.tau_ij
		tau[p.j * nc + p.i] = p.tau_ji
		alpha[p.i * nc + p.j] = p.alpha
		alpha[p.j * nc + p.i] = p.alpha
	}
	every_pair('NRTL', names, given)!
	return Nrtl{
		names: names.clone()
		tau:   tau
		alpha: alpha
	}
}

pub fn (m &Nrtl) size() int {
	return m.names.len
}

pub fn (m &Nrtl) name() string {
	return 'NRTL'
}

// work holds τ and G, nc² each, and A and B.
pub fn (m &Nrtl) work() Work {
	nc := m.names.len
	return new_work(nc, 2 * nc * nc, 2 * nc)
}

pub fn (m &Nrtl) params(t f64, mut w Work) ! {
	nc := m.names.len
	nn := nc * nc
	for i in 0 .. nc {
		for j in 0 .. nc {
			k := i * nc + j
			if i == j {
				w.pj[k] = numeric.constant(0.0)
				w.pj[nn + k] = numeric.constant(1.0)
			} else {
				tau := finite_param(m, &m.tau[k], t, 'τ of ${m.names[i]} with ${m.names[j]}')!
				alpha := finite_param(m, &m.alpha[k], t, 'α of ${m.names[i]} with ${m.names[j]}')!
				w.pj[k] = tau
				g := (alpha * tau).neg().exp()
				if !usable(g) {
					return unusable(m, 'G of ${m.names[i]} with ${m.names[j]}, exp(−ατ),',
						g, t)
				}
				w.pj[nn + k] = g
			}
			w.pv[k] = w.pj[k].v
			w.pv[nn + k] = w.pj[nn + k].v
		}
	}
}

pub fn (m &Nrtl) lng(p []f64, x []f64, mut s []f64, mut out []f64) {
	nrtl[f64](m.names.len, p, x, mut s, mut out)
}

pub fn (m &Nrtl) lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet) {
	nrtl[numeric.Jet](m.names.len, p, x, mut s, mut out)
}

// nrtl is NRTL's ln γ: p holds τ then G, s takes A then B.
fn nrtl[N](nc int, p []N, x []N, mut s []N, mut out []N) {
	nn := nc * nc
	for j in 0 .. nc {
		mut a := lift[N](0.0)
		mut b := lift[N](0.0)
		for k in 0 .. nc {
			g := p[nn + k * nc + j]
			b = b + x[k] * g
			a = a + x[k] * p[k * nc + j] * g
		}
		s[j] = a
		s[nc + j] = b
	}
	for i in 0 .. nc {
		mut v := s[i] / s[nc + i]
		for j in 0 .. nc {
			v = v + x[j] * p[nn + i * nc + j] / s[nc + j] * (p[i * nc + j] - s[j] / s[nc + j])
		}
		out[i] = v
	}
}
