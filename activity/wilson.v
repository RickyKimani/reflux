module activity

import numeric

// Wilson is Wilson's equation (1964):
//
//     ln γᵢ = 1 − ln Sᵢ − Σₖ xₖΛₖᵢ/Sₖ,   Sᵢ = Σⱼ xⱼΛᵢⱼ,
//
// with Λᵢᵢ = 1. Each Λᵢⱼ is an expression in T (D7); DECHEMA's is
// (Vⱼ/Vᵢ)·exp(−Aᵢⱼ/RT), V each component's liquid volume. It can't predict
// two liquids.
pub struct Wilson {
pub:
	names  []string
	lambda []Param // nc×nc, row-major; the diagonal unused
}

// WilsonPair is one pair's parameters, i before j: Λᵢⱼ and Λⱼᵢ.
pub struct WilsonPair {
pub:
	i         int
	j         int
	lambda_ij Param
	lambda_ji Param
}

// new_wilson is Wilson's equation for components named names, every pair
// given once.
pub fn new_wilson(names []string, pairs []WilsonPair) !Wilson {
	nc := names.len
	one := param('1', expr_definitions(), 'Λ')!
	mut lambda := []Param{len: nc * nc, init: one}
	mut given := []bool{len: nc * nc}
	for p in pairs {
		check_pair('Wilson', names, p.i, p.j, mut given)!
		lambda[p.i * nc + p.j] = p.lambda_ij
		lambda[p.j * nc + p.i] = p.lambda_ji
	}
	every_pair('Wilson', names, given)!
	return Wilson{
		names:  names.clone()
		lambda: lambda
	}
}

pub fn (m &Wilson) size() int {
	return m.names.len
}

pub fn (m &Wilson) name() string {
	return 'Wilson'
}

// work holds Λ, nc², and S.
pub fn (m &Wilson) work() Work {
	nc := m.names.len
	return new_work(nc, nc * nc, nc)
}

pub fn (m &Wilson) params(t f64, mut w Work) ! {
	nc := m.names.len
	for i in 0 .. nc {
		for j in 0 .. nc {
			k := i * nc + j
			w.pj[k] = if i == j {
				numeric.constant(1.0)
			} else {
				finite_param(m, &m.lambda[k], t, 'Λ of ${m.names[i]} with ${m.names[j]}')!
			}
			if !(w.pj[k].v > 0.0) {
				return refuse(m, 'Λ of ${m.names[i]} with ${m.names[j]} is ${w.pj[k].v} at ${t} K')
			}
			w.pv[k] = w.pj[k].v
		}
	}
}

pub fn (m &Wilson) lng(p []f64, x []f64, mut s []f64, mut out []f64) {
	wilson[f64](m.names.len, p, x, mut s, mut out)
}

pub fn (m &Wilson) lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet) {
	wilson[numeric.Jet](m.names.len, p, x, mut s, mut out)
}

// wilson is Wilson's ln γ: p holds Λ, s takes S.
fn wilson[N](nc int, p []N, x []N, mut s []N, mut out []N) {
	for i in 0 .. nc {
		mut si := lift[N](0.0)
		for j in 0 .. nc {
			si = si + x[j] * p[i * nc + j]
		}
		s[i] = si
	}
	for i in 0 .. nc {
		mut v := lift[N](1.0) - nlog[N](s[i])
		for k in 0 .. nc {
			v = v - x[k] * p[k * nc + i] / s[k]
		}
		out[i] = v
	}
}
