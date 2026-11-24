module thermo

import activity
import math
import numeric
import tdep

// GammaPhi is a γ–φ fluid (M6's D4 to D6): a liquid by its activity
// coefficients, each component's vapour pressure and liquid volume, and a
// vapour that's the ideal gas or an equation of state's. It's a PhaseModel,
// with no F(T, V, n) behind its liquid, so the flash takes it as it takes a
// cubic, and what works in F(T, V, n) refuses it (Fluid.helmholtz).
//
// A solvent's liquid fugacity is
//
//     fᵢ = xᵢ γᵢ Pᵢˢᵃᵗ φᵢˢᵃᵗ exp(Vᵢᴸ(P − Pᵢˢᵃᵗ)/RT),
//
// γ from the activity model over the solvents alone, on their own mole
// fractions; φˢᵃᵗ the vapour model's pure φ at Pˢᵃᵗ, 1 for the ideal gas;
// the exponential, Poynting's, if asked for. A Henry's-law solute's is
// fᵢ = xᵢ Hᵢ,mix, its H from the solvents that have parameters for it.
// Everything the phase gives (ln φ and its derivatives, H, S, Cp) comes
// from those fugacities and their exact T-derivatives (D5).
pub struct GammaPhi {
pub:
	comps []Component
	// liquid is the activity model over the solvents, in solvents' order.
	liquid   activity.Model
	solvents []int
	// Each component's vapour pressure, Pa, and liquid molar density,
	// mol/m³, from their correlations; a solute's are unused.
	psat    []Scaled
	density []Scaled
	// vapour is the vapour's equation of state, or none for the ideal gas;
	// pures are each component alone in it, for φˢᵃᵗ.
	vapour   ?Model
	pures    []Model
	poynting bool
	solutes  []Solute
mut:
	work ?&GammaWork // a workspace, for one thread at a time (scoped)
}

// Scaled is a correlation in SI: scale times its value.
pub struct Scaled {
pub:
	c     tdep.Correlation
	scale f64
}

// Solute is a Henry's-law component, with its constants in each solvent
// that has them (D6).
pub struct Solute {
pub:
	i      int
	henrys []Henry
}

// Henry is a solute's constant in one solvent, between tmin and tmax, in
// one of two forms (HenryForm).
pub struct Henry {
pub:
	solvent int
	form    HenryForm
	a       f64
	b       f64
	c       f64
	d       f64 // equation 121's alone
	e       f64 // equation 121's alone
	tmin    f64
	tmax    f64
}

// HenryForm is how a Henry's-law constant depends on T.
pub enum HenryForm {
	// ChemSep's equation 121, ln(H/Pa) = A + B/T + C ln T + D·Tᴱ.
	chemsep_121
	// IAPWS's for a gas in water (G7-04, Fernández-Prini, Alvarez and
	// Harvey, 2003): ln(kH/p*) = A/Tr + B·τ^0.355/Tr + C·Tr^−0.41·e^τ, with
	// Tr = T/Tc and τ = 1 − Tr on water's critical point, and p* water's
	// vapour pressure by Wagner and Pruss.
	iapws
}

// Water's critical point and Wagner and Pruss's vapour pressure, IAPWS's
// supplementary release of 1992, which IAPWS G7-04 takes p* from:
// ln(p*/pc) = (Tc/T)·Σ aᵢτ^bᵢ.
const water_tc = 647.096
const water_pc = 22.064e6
const water_psat_a = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]
const water_psat_b = [1.0, 1.5, 3.0, 3.5, 4.0, 7.5]

// jet is ln(H/Pa) at T with its first two T-derivatives: t is T along
// itself (numeric.along(T)). Within tmin and tmax, which are below water's
// critical point for IAPWS's form, it is finite.
pub fn (h Henry) jet(t numeric.Jet) numeric.Jet {
	c := numeric.constant
	return match h.form {
		.chemsep_121 {
			c(h.a) + c(h.b) / t + c(h.c) * t.log() + c(h.d) * t.pow(h.e)
		}
		.iapws {
			tr := t / c(water_tc)
			tau := c(1.0) - tr
			mut s := c(0.0)
			for k, ak in water_psat_a {
				s = s + c(ak) * tau.pow(water_psat_b[k])
			}
			c(math.log(water_pc)) + s / tr + c(h.a) / tr + c(h.b) * tau.pow(0.355) / tr +
				c(h.c) * tr.pow(-0.41) * tau.exp()
		}
	}
}

// GammaWork keeps what depends on T alone until T moves (D4): each
// component's ln Pˢᵃᵗ, Vᴸ and ln φˢᵃᵗ, and each solute's ln H in each of
// its solvents, as jets in T; and the activity model's own workspace.
@[heap]
pub struct GammaWork {
mut:
	t        f64 // 0 before any, and after a refusal
	lnpsat   []numeric.Jet
	vl       []numeric.Jet
	lnphisat []numeric.Jet
	lnh      [][]numeric.Jet
	act      activity.Work
	x        []f64 // the solvents' own mole fractions
}

// new_gamma_phi checks a γ–φ fluid's parts and makes it.
pub fn new_gamma_phi(g GammaPhi) !GammaPhi {
	nc := g.comps.len
	if g.psat.len != nc || g.density.len != nc {
		return error('γ–φ: ${g.psat.len} vapour pressures and ${g.density.len} densities for ${nc} components')
	}
	if g.liquid.size() != g.solvents.len {
		return error('γ–φ: an activity model of ${g.liquid.size()} components for ${g.solvents.len} solvents')
	}
	mut role := []int{len: nc} // 1 a solvent, 2 a solute
	for i in g.solvents {
		if i < 0 || i >= nc || role[i] != 0 {
			return error('γ–φ: solvent ${i} of ${nc} components, or listed twice')
		}
		role[i] = 1
	}
	for s in g.solutes {
		if s.i < 0 || s.i >= nc || role[s.i] != 0 {
			return error('γ–φ: solute ${s.i} of ${nc} components, or listed twice')
		}
		role[s.i] = 2
		if s.henrys.len == 0 {
			return error("γ–φ: ${g.comps[s.i].name} is a Henry's-law solute with no constants")
		}
		for h in s.henrys {
			if h.solvent < 0 || h.solvent >= nc || role[h.solvent] == 2 {
				return error("γ–φ: ${g.comps[s.i].name}'s Henry's-law constant is in component ${h.solvent}, which isn't a solvent")
			}
		}
	}
	if role.any(it == 0) {
		return error('γ–φ: every component must be a solvent or a solute')
	}
	for i in g.solvents {
		for c in [g.psat[i].c, g.density[i].c] {
			if c.eq !in tdep.jet_forms {
				return error('γ–φ: ${g.comps[i].name} has a correlation of form ${c.eq}, which has no exact derivatives here (only ${tdep.jet_forms})')
			}
		}
	}
	if m := g.vapour {
		if m.components().len != nc || g.pures.len != nc {
			return error('γ–φ: a vapour model of ${m.components().len} components, and ${g.pures.len} pure ones, for ${nc}')
		}
	}
	return GammaPhi{
		...g
		work: none
	}
}

// scoped is this model with a workspace of its own, which keeps what
// depends on T alone between calls: for one thread, one flash at a time.
// Without one, each call makes its own.
pub fn (g GammaPhi) scoped() PhaseModel {
	return GammaPhi{
		...g
		work: g.new_work()
	}
}

fn (g GammaPhi) new_work() &GammaWork {
	nc := g.comps.len
	return &GammaWork{
		lnpsat:   []numeric.Jet{len: nc}
		vl:       []numeric.Jet{len: nc}
		lnphisat: []numeric.Jet{len: nc}
		lnh:      g.solutes.map([]numeric.Jet{len: it.henrys.len})
		act:      g.liquid.work()
		x:        []f64{len: g.solvents.len}
	}
}

fn (g GammaPhi) workspace() &GammaWork {
	if w := g.work {
		return w
	}
	return g.new_work()
}

// components implements PhaseModel.
pub fn (g GammaPhi) components() []Component {
	return g.comps
}

fn (g GammaPhi) refuse(reason string) IError {
	return error('γ–φ: ${reason}')
}

// at_t fills w's T-dependent parts at t, or refuses a T outside a
// correlation's range (D4): never extrapolated.
fn (g GammaPhi) at_t(t f64, mut w GammaWork) ! {
	if w.t == t {
		return
	}
	w.t = 0.0
	if !(t > 0.0) || !math.is_finite(t) {
		return g.refuse('a temperature of ${t} K')
	}
	tj := numeric.along(t)
	for i in g.solvents {
		c := g.comps[i]
		for which, cor in {
			'vapour pressure': g.psat[i]
			'liquid density':  g.density[i]
		} {
			if !cor.c.in_range(t) {
				return g.refuse("${c.name}'s ${which} is fitted from ${cor.c.tmin} K to ${cor.c.tmax} K, not at ${t} K")
			}
		}
		p := g.psat[i].c.jet(tj)!
		w.lnpsat[i] = (numeric.constant(g.psat[i].scale) * p).log()
		rho := numeric.constant(g.density[i].scale) * g.density[i].c.jet(tj)!
		w.vl[i] = numeric.constant(1.0) / rho
		if !(w.lnpsat[i].v > -700.0 && w.vl[i].v > 0.0) || !finite_jet(w.lnpsat[i]) || !finite_jet(w.vl[i]) {
			return g.refuse('${c.name} has Pˢᵃᵗ = ${math.exp(w.lnpsat[i].v)} Pa and Vᴸ = ${w.vl[i].v} m³/mol at ${t} K')
		}
		w.lnphisat[i] = if g.vapour != none {
			g.ln_phi_sat(i, t, w.lnpsat[i])!
		} else {
			numeric.constant(0.0)
		}
	}
	for k, s in g.solutes {
		for q, h in s.henrys {
			if !(t >= h.tmin && t <= h.tmax) {
				return g.refuse("${g.comps[s.i].name}'s Henry's-law constant in ${g.comps[h.solvent].name} is fitted from ${h.tmin} K to ${h.tmax} K, not at ${t} K")
			}
			w.lnh[k][q] = h.jet(tj)
		}
	}
	w.t = t
}

fn finite_jet(j numeric.Jet) bool {
	return math.is_finite(j.v) && math.is_finite(j.d) && math.is_finite(j.dd)
}

// ln_phi_sat is component i's pure φ in the vapour model at its own Pˢᵃᵗ,
// on the vapour root, as a jet in T: its first derivative exact,
// ∂ln φ/∂T + ∂ln φ/∂P·dPˢᵃᵗ/dT, and its second by Richardson's
// extrapolation of central differences of the first, since a cubic has no
// third derivatives of F (design.md, "The γ–φ liquid").
fn (g GammaPhi) ln_phi_sat(i int, t f64, lnpsat numeric.Jet) !numeric.Jet {
	v, d := g.phi_sat_slope(i, t, lnpsat)!
	h := 1e-4 * t
	mut s := []f64{len: 4}
	for k, dt in [h, -h, h / 2.0, -h / 2.0] {
		lp := (numeric.constant(g.psat[i].scale) * g.psat[i].c.jet(numeric.along(t + dt))!).log()
		_, s[k] = g.phi_sat_slope(i, t + dt, lp)!
	}
	dd := (4.0 * (s[2] - s[3]) / h - (s[0] - s[1]) / (2.0 * h)) / 3.0
	return numeric.Jet{v, d, dd}
}

fn (g GammaPhi) phi_sat_slope(i int, t f64, lnpsat numeric.Jet) !(f64, f64) {
	ps := math.exp(lnpsat.v)
	m := g.pures[i]
	fx := fix(m, t, [1.0])
	v := vapor_volume(fx, ps) or {
		return g.refuse('${g.comps[i].name} has no vapour root at its vapour pressure, ${ps} Pa, at ${t} K, so no φˢᵃᵗ')
	}
	st := state_at(fx, ps, v, [1.0])
	d := ln_phi_derivs(m, st)
	return st.ln_phi[0], d.dt[0] + d.dp[0] * ps * lnpsat.d
}

// liquid_volume is the liquid's volume, Σ nᵢVᵢᴸ over the solvents (a
// solute adds none), with its T-derivative.
fn (g GammaPhi) liquid_volume(n []f64, w &GammaWork) (f64, f64) {
	mut v, mut dv := 0.0, 0.0
	for i in g.solvents {
		v += n[i] * w.vl[i].v
		dv += n[i] * w.vl[i].d
	}
	return v, dv
}

// Liquid is the liquid's ln φ at (T, P, x) as jets in T, and, if asked,
// n·∂ln φᵢ/∂nⱼ.
struct Liquid {
	lnphi []numeric.Jet
	dn    []f64
}

fn (g GammaPhi) liquid_at(t f64, p f64, x []f64, want activity.Want, mut w GammaWork) !Liquid {
	g.at_t(t, mut w)!
	nc := g.comps.len
	mut s := 0.0
	for i in g.solvents {
		s += x[i]
	}
	if !(s > 0.0) {
		return g.refuse("a liquid of Henry's-law solutes alone")
	}
	for k, i in g.solvents {
		w.x[k] = x[i] / s
	}
	activity.eval(g.liquid, t, w.x, want, mut w.act)!
	c := numeric.constant
	mut lnphi := []numeric.Jet{len: nc}
	lnp := math.log(p)
	for k, i in g.solvents {
		gam := numeric.Jet{w.act.ln_gamma[k], w.act.dt[k], w.act.dtt[k]}
		mut f := gam + w.lnpsat[i] + w.lnphisat[i] - c(lnp)
		if g.poynting {
			// V(P − Pˢᵃᵗ)/RT, V and Pˢᵃᵗ moving with T.
			ps := w.lnpsat[i].exp()
			f = f + w.vl[i] * (c(p) - ps) / (c(r_gas) * numeric.along(t))
		}
		lnphi[i] = f
	}
	mut dn := []f64{}
	if want.has(.n) {
		dn = []f64{len: nc * nc}
		for k, i in g.solvents {
			for q, j in g.solvents {
				dn[i * nc + j] = w.act.dn[k * g.solvents.len + q] / s
			}
		}
	}
	for k, sol in g.solutes {
		// ln H mixed over the solvents that have constants, Σ xⱼ ln Hᵢⱼ / Σ xⱼ,
		// an absent one weighing nothing.
		mut num := c(0.0)
		mut den := 0.0
		for q, h in sol.henrys {
			num = num + c(x[h.solvent]) * w.lnh[k][q]
			den += x[h.solvent]
		}
		if !(den > 0.0) {
			return g.refuse("${g.comps[sol.i].name} has a Henry's-law constant in none of the solvents present")
		}
		lnhmix := num / c(den)
		lnphi[sol.i] = lnhmix - c(lnp)
		if want.has(.n) {
			// n·∂ln Hᵢ,mix/∂nⱼ, an absent solvent's too: a little of it moves
			// the mix towards its own H.
			for q, h in sol.henrys {
				dn[sol.i * nc + h.solvent] = (w.lnh[k][q].v - lnhmix.v) / den
			}
		}
	}
	return Liquid{lnphi, dn}
}

// phase implements PhaseModel: the liquid, the vapour, or the one of lower
// G. Where the liquid can't be had (a T outside a correlation's range), the
// refusal is the answer for .liquid and .stable alike, never the vapour in
// its place (D4).
pub fn (g GammaPhi) phase(t f64, p f64, n []f64, root Root) !State {
	if !(p > 0.0) || !math.is_finite(p) {
		return g.refuse('a pressure of ${p} Pa')
	}
	nc := g.comps.len
	if n.len != nc {
		return g.refuse('${n.len} amounts for ${nc} components')
	}
	mut total := 0.0
	for v in n {
		if !(v >= 0.0) || !math.is_finite(v) {
			return g.refuse('an amount of ${v} mol')
		}
		total += v
	}
	if !(total > 0.0) {
		return g.refuse('no moles')
	}
	if root == .vapor {
		return g.vapour_phase(t, p, n, total)
	}
	mut w := g.workspace()
	x := n.map(it / total)
	liq := g.liquid_at(t, p, x, activity.ln_gamma_only, mut w)!
	vl, _ := g.liquid_volume(n, w)
	liquid := State{
		t:      t
		p:      p
		n:      n.clone()
		v:      vl
		z:      p * vl / (total * r_gas * t)
		ln_phi: liq.lnphi.map(it.v)
	}
	if root == .liquid {
		return liquid
	}
	vapour := g.vapour_phase(t, p, n, total)!
	mut gl, mut gv := 0.0, 0.0
	for i in 0 .. nc {
		gl += n[i] * liquid.ln_phi[i]
		gv += n[i] * vapour.ln_phi[i]
	}
	return if gl < gv { liquid } else { vapour }
}

fn (g GammaPhi) vapour_phase(t f64, p f64, n []f64, total f64) !State {
	if m := g.vapour {
		st := solve_tp(m, t, p, n, .vapor)!
		if st.v == g.maybe_liquid_volume(t, n) {
			return g.refuse("the vapour model gave a phase with the liquid's volume")
		}
		return st
	}
	if !(t > 0.0) || !math.is_finite(t) {
		return g.refuse('a temperature of ${t} K')
	}
	v := total * r_gas * t / p
	return State{
		t:      t
		p:      p
		n:      n.clone()
		v:      v
		z:      1.0
		ln_phi: []f64{len: n.len}
	}
}

// maybe_liquid_volume is the liquid's volume at t, or NaN where there's
// none: what tells the model's two phases apart.
fn (g GammaPhi) maybe_liquid_volume(t f64, n []f64) f64 {
	mut w := g.workspace()
	g.at_t(t, mut w) or { return math.nan() }
	v, _ := g.liquid_volume(n, w)
	return v
}

// is_liquid says which of the model's phases a state is: the liquid's
// volume is Σ nᵢVᵢᴸ to the bit, and the vapour's never is.
fn (g GammaPhi) is_liquid(st State) bool {
	return st.v == g.maybe_liquid_volume(st.t, st.n)
}

// ln_phi_derivs implements PhaseModel.
pub fn (g GammaPhi) ln_phi_derivs(st State) !LnPhiDerivs {
	nc := g.comps.len
	if !g.is_liquid(st) {
		if m := g.vapour {
			return ln_phi_derivs(m, st)
		}
		return LnPhiDerivs{
			dt: []f64{len: nc}
			dp: []f64{len: nc}
			dn: []f64{len: nc * nc}
		}
	}
	mut w := g.workspace()
	total := numeric.sum(st.n)
	liq := g.liquid_at(st.t, st.p, st.n.map(it / total), .t | .n, mut w)!
	mut dp := []f64{len: nc}
	rt := r_gas * st.t
	for i in 0 .. nc {
		dp[i] = -1.0 / st.p
	}
	if g.poynting {
		for i in g.solvents {
			dp[i] += w.vl[i].v / rt
		}
	}
	return LnPhiDerivs{
		dt: liq.lnphi.map(it.d)
		dp: dp
		dn: liq.dn
	}
}

// phase_label implements PhaseModel: by the model that made the phase (D9).
pub fn (g GammaPhi) phase_label(st State) PhaseLabel {
	return if g.is_liquid(st) { .liquid } else { .vapor }
}

// phase_residual implements PhaseModel. The liquid's, from its fugacities
// (D5): Hʳ = −RT²Σnᵢ∂ln φᵢ/∂T, Gʳ = RTΣnᵢ ln φᵢ, Sʳ = (Hʳ − Gʳ)/T,
// Cpʳ = ∂Hʳ/∂T. An incompressible liquid has no Cv, and so no speed of
// sound: those are NaN (PhaseModel).
pub fn (g GammaPhi) phase_residual(st State) PhaseResidual {
	nan := math.nan()
	if !g.is_liquid(st) {
		if m := g.vapour {
			return VolumeRoots{m}.phase_residual(st)
		}
		total := numeric.sum(st.n)
		return PhaseResidual{
			dv_dt: total * r_gas / st.p
			dp_dv: -st.p / st.v
		}
	}
	mut w := g.workspace()
	total := numeric.sum(st.n)
	liq := g.liquid_at(st.t, st.p, st.n.map(it / total), .t, mut w) or {
		return PhaseResidual{nan, nan, nan, nan, nan, nan}
	}
	t := st.t
	mut d1, mut d2, mut gr := 0.0, 0.0, 0.0
	for i, ni in st.n {
		// An absent component's ln φ is finite (eval answers finitely or
		// refuses), so it adds nothing.
		d1 += ni * liq.lnphi[i].d
		d2 += ni * liq.lnphi[i].dd
		gr += ni * liq.lnphi[i].v
	}
	rt := r_gas * t
	h := -rt * t * d1
	_, dv := g.liquid_volume(st.n, w)
	return PhaseResidual{
		h:     h
		s:     (h - rt * gr) / t
		cp:    -r_gas * t * (2.0 * d1 + t * d2)
		cv:    nan
		dv_dt: dv
		dp_dv: math.inf(-1)
	}
}
