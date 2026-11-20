module activity

import math
import numeric

// Variant is a UNIFAC: its combinatorial term and its table.
pub enum Variant {
	original // Fredenslund et al., for VLE, with Hansen et al.'s table (1991)
	lle      // the original's equations with Magnussen et al.'s LLE table (1981)
	dortmund // Weidlich and Gmehling's modification, with Gmehling et al.'s published table
}

pub fn (v Variant) name() string {
	return match v {
		.original { 'UNIFAC' }
		.lle { 'UNIFAC-LLE' }
		.dortmund { 'UNIFAC-Dortmund' }
	}
}

// Subgroup is one of a table's subgroups: its main group, and its volume
// and area parameters R and Q.
pub struct Subgroup {
pub:
	id   int
	name string
	main int
	r    f64
	q    f64
	mw   f64 // g/mol, as the table gives it
}

// Interaction is the parameter of main group m with main group n:
// ψₘₙ = exp(−(a + b·T + c·T²)/T), a in K, b dimensionless, c in 1/K. The
// original and LLE tables have b = c = 0.
pub struct Interaction {
pub:
	a f64
	b f64
	c f64
}

// Table is a variant's subgroups, main groups and interaction parameters,
// by the numbers its source gives them.
pub struct Table {
pub:
	variant   Variant
	subgroups map[int]Subgroup
	mains     map[int]string
	// interactions by m·1000 + n (interaction_key), aₘₙ and not aₙₘ.
	interactions map[int]Interaction
}

// interaction_key is the key of main group m's parameter with n.
pub fn interaction_key(m int, n int) int {
	return m * 1000 + n
}

// Unifac is a UNIFAC model of some components, each a set of subgroups:
//
//     ln γᵢ = ln γᵢᶜ + Σₖ νₖ⁽ⁱ⁾·(ln Γₖ − ln Γₖ⁽ⁱ⁾),
//     ln Γₖ = Qₖ·(1 − ln Σₘ Θₘψₘₖ − Σₘ Θₘψₖₘ/Σₙ Θₙψₙₘ),
//     Θₘ = QₘXₘ/Σₙ QₙXₙ,   Xₘ = Σᵢ νₘ⁽ⁱ⁾xᵢ / Σᵢ Σₙ νₙ⁽ⁱ⁾xᵢ,
//
// ln Γₖ⁽ⁱ⁾ being ln Γₖ in pure i. The combinatorial term is
//
//     ln γᵢᶜ = 1 − Vᵢ′ + ln Vᵢ′ − 5qᵢ·(1 − Vᵢ/Fᵢ + ln(Vᵢ/Fᵢ)),
//     Vᵢ = rᵢ/Σⱼ xⱼrⱼ,   Fᵢ = qᵢ/Σⱼ xⱼqⱼ,
//
// with Vᵢ′ = Vᵢ in the original and LLE forms, and rᵢ^(3/4)/Σⱼ xⱼrⱼ^(3/4)
// in Dortmund's; rᵢ = Σₖ νₖ⁽ⁱ⁾Rₖ and qᵢ = Σₖ νₖ⁽ⁱ⁾Qₖ.
pub struct Unifac {
pub:
	variant Variant
	names   []string
	groups  []int // the subgroups present, by the table's numbers
	gq      []f64 // each group's Q
	nu      []f64 // νₖ⁽ⁱ⁾, nc×ng: component i's count of group k
	r       []f64 // each component's r
	q       []f64 // and q
	rv      []f64 // the r in V′: r, or r^(3/4)
	ia      []f64 // each pair of groups' a, ng×ng
	ib      []f64
	ic      []f64
}

// new_unifac is UNIFAC for components named names, each given as its
// subgroups' counts, by the table's numbers. A subgroup the table doesn't
// have, or two groups whose main groups have no parameters, is refused.
pub fn new_unifac(table Table, names []string, groups []map[int]int) !Unifac {
	model := table.variant.name()
	nc := names.len
	if groups.len != nc {
		return pair_error(model, '${groups.len} group lists for ${nc} components')
	}
	mut present := []int{}
	for i, g in groups {
		if g.len == 0 {
			return pair_error(model, '${names[i]} has no ${model} groups')
		}
		for id, count in g {
			sg := table.subgroups[id] or {
				return pair_error(model, "${names[i]}'s subgroup ${id} isn't in the table")
			}
			if count <= 0 {
				return pair_error(model, "${names[i]}'s subgroup ${id} (${sg.name}) counts ${count}")
			}
			if id !in present {
				present << id
			}
		}
	}
	present.sort()
	ng := present.len
	mut nu := []f64{len: nc * ng}
	mut r := []f64{len: nc}
	mut q := []f64{len: nc}
	for i, g in groups {
		for k, id in present {
			count := g[id] or { 0 }
			nu[i * ng + k] = f64(count)
		}
		// Each sum in the groups' order, so that it's the same whichever
		// way the map was built.
		for k, id in present {
			sg := table.subgroups[id] or { continue }
			r[i] += nu[i * ng + k] * sg.r
			q[i] += nu[i * ng + k] * sg.q
		}
	}
	mut ia := []f64{len: ng * ng}
	mut ib := []f64{len: ng * ng}
	mut ic := []f64{len: ng * ng}
	for k, idk in present {
		mk := table.subgroups[idk] or { continue }.main
		for l, idl in present {
			ml := table.subgroups[idl] or { continue }.main
			if mk == ml {
				continue
			}
			p := table.interactions[interaction_key(mk, ml)] or {
				return pair_error(model, 'no parameters for main group ${mk} (${table.mains[mk] or {
					'?'
				}}) with ${ml} (${table.mains[ml] or { '?' }})')
			}
			ia[k * ng + l] = p.a
			ib[k * ng + l] = p.b
			ic[k * ng + l] = p.c
		}
	}
	return Unifac{
		variant: table.variant
		names:   names.clone()
		groups:  present
		gq:      present.map(table.subgroups[it] or { Subgroup{} }.q)
		nu:      nu
		r:       r
		q:       q
		rv:      if table.variant == .dortmund { r.map(math.pow(it, 0.75)) } else { r.clone() }
		ia:      ia
		ib:      ib
		ic:      ic
	}
}

pub fn (m &Unifac) size() int {
	return m.names.len
}

pub fn (m &Unifac) name() string {
	return m.variant.name()
}

// work holds ψ, ng², then ln Γₖ⁽ⁱ⁾, nc×ng; its scratch is a composition,
// X, Θ, the sums Σₘ Θₘψₘₖ and ln Γₖ.
pub fn (m &Unifac) work() Work {
	nc := m.names.len
	ng := m.groups.len
	return new_work(nc, ng * ng + nc * ng, nc + 4 * ng)
}

// params fills ψ at t, then each pure component's ln Γₖ⁽ⁱ⁾ by the same
// arithmetic as ln Γₖ in a mixture of i alone, so that a pure component's γ
// is exactly 1.
pub fn (m &Unifac) params(t f64, mut w Work) ! {
	nc := m.names.len
	ng := m.groups.len
	tj := numeric.along(t)
	for k in 0 .. ng * ng {
		e := (numeric.constant(m.ia[k]) + numeric.constant(m.ib[k]) * tj +
			numeric.constant(m.ic[k]) * tj * tj) / tj
		w.pj[k] = e.neg().exp()
		w.pv[k] = w.pj[k].v
		if !(w.pv[k] > 0.0) || !math.is_finite(w.pv[k]) {
			return refuse(m, 'ψ is ${w.pv[k]} at ${t} K for groups ${m.groups[k / ng]} and ${m.groups[k % ng]}')
		}
	}
	for i in 0 .. nc {
		for j in 0 .. nc {
			w.sf[j] = if j == i { 1.0 } else { 0.0 }
			w.sj[j] = numeric.constant(w.sf[j])
		}
		group_lngamma[f64](m.nu, m.gq, nc, w.pv, mut w.sf)
		group_lngamma[numeric.Jet](m.nu, m.gq, nc, w.pj, mut w.sj)
		off := nc + 3 * ng
		for k in 0 .. ng {
			w.pv[ng * ng + i * ng + k] = w.sf[off + k]
			w.pj[ng * ng + i * ng + k] = w.sj[off + k]
		}
	}
}

pub fn (m &Unifac) lng(p []f64, x []f64, mut s []f64, mut out []f64) {
	unifac[f64](m, p, x, mut s, mut out)
}

pub fn (m &Unifac) lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet) {
	unifac[numeric.Jet](m, p, x, mut s, mut out)
}

// group_lngamma is each group's ln Γₖ for the composition in s[0..nc), with
// ψ in p, into s[nc + 3ng..nc + 4ng), using the rest of s as it goes: X,
// Θ, then the sums Σₘ Θₘψₘₖ.
fn group_lngamma[N](nu []f64, gq []f64, nc int, p []N, mut s []N) {
	ng := gq.len
	ox, ot, os, og := nc, nc + ng, nc + 2 * ng, nc + 3 * ng
	mut total := lift[N](0.0)
	for k in 0 .. ng {
		mut xk := lift[N](0.0)
		for i in 0 .. nc {
			xk = xk + lift[N](nu[i * ng + k]) * s[i]
		}
		s[ox + k] = xk
		total = total + xk
	}
	mut qx := lift[N](0.0)
	for k in 0 .. ng {
		s[ox + k] = s[ox + k] / total
		qx = qx + lift[N](gq[k]) * s[ox + k]
	}
	for k in 0 .. ng {
		s[ot + k] = lift[N](gq[k]) * s[ox + k] / qx
	}
	for k in 0 .. ng {
		mut sk := lift[N](0.0)
		for mm in 0 .. ng {
			sk = sk + s[ot + mm] * p[mm * ng + k]
		}
		s[os + k] = sk
	}
	for k in 0 .. ng {
		mut v := lift[N](0.0)
		for mm in 0 .. ng {
			v = v + s[ot + mm] * p[k * ng + mm] / s[os + mm]
		}
		s[og + k] = lift[N](gq[k]) * (lift[N](1.0) - nlog[N](s[os + k]) - v)
	}
}

// unifac is UNIFAC's ln γ.
fn unifac[N](m &Unifac, p []N, x []N, mut s []N, mut out []N) {
	nc := m.names.len
	ng := m.groups.len
	mut rbar := lift[N](0.0)
	mut rvbar := lift[N](0.0)
	mut qbar := lift[N](0.0)
	for j in 0 .. nc {
		rbar = rbar + x[j] * lift[N](m.r[j])
		rvbar = rvbar + x[j] * lift[N](m.rv[j])
		qbar = qbar + x[j] * lift[N](m.q[j])
		s[j] = x[j]
	}
	group_lngamma[N](m.nu, m.gq, nc, p, mut s)
	og := nc + 3 * ng
	for i in 0 .. nc {
		v := lift[N](m.r[i]) / rbar
		vv := lift[N](m.rv[i]) / rvbar
		vf := v / (lift[N](m.q[i]) / qbar)
		mut res := lift[N](0.0)
		for k in 0 .. ng {
			nk := m.nu[i * ng + k]
			if nk != 0.0 {
				res = res + lift[N](nk) * (s[og + k] - p[ng * ng + i * ng + k])
			}
		}
		out[i] = lift[N](1.0) - vv + nlog[N](vv) - lift[N](5.0 * m.q[i]) * (lift[N](1.0) - vf +
			nlog[N](vf)) + res
	}
}
