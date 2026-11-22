// Package activity is M6's liquid activity models: a mixture's excess Gibbs
// energy, as each component's activity coefficient γᵢ at a temperature and
// composition (D2). It knows nothing of phases, pressures or volumes: the
// γ–φ liquid in thermo is built on it (D4).
//
// Each model's formula is written once, generic over f64 and numeric.Jet,
// so its value comes from f64 arithmetic and its derivatives from the same
// formula on jets, by forward-mode differentiation, exactly to rounding:
// the T-derivatives with each parameter's own T-derivatives, and the
// composition derivatives one mole fraction at a time.
module activity

import expr
import math
import numeric
import quantity

// r_gas is the molar gas constant, J/(mol·K): N_A·k_B, exact since the 2019
// SI, as thermo.r_gas and expr's R.
pub const r_gas = 8.314462618

// least_normal is the smallest normal double, 2⁻¹⁰²². Below it a number has
// lost digits, and at 0 all of them, so a parameter there (a G, τ, Λ or ψ,
// which the formulas divide by or take the log of) is refused: its ln γ
// could be finite and wrong.
const least_normal = 2.2250738585072014e-308

// sum_slack is how far from 1 a composition's mole fractions may sum. It's
// far above normalising's rounding, a few ε for each component, and moves
// a ln γ by about as much as itself.
const sum_slack = 1e-12

// Model is a liquid's excess Gibbs energy. A model holds nothing about an
// evaluation, so threads may share one; each thread brings its own Work,
// which keeps the parameters that depend on T alone until T moves (M8's
// D11).
pub interface Model {
	// size is the number of components.
	size() int
	// name is the model's, for messages: "NRTL".
	name() string
	// work makes a workspace sized for the model.
	work() Work
	// params fills w's T-dependent parameters at t: their values, and
	// their first two T-derivatives as jets. It refuses a T where one isn't
	// what the formula can take.
	params(t f64, mut w Work) !
	// lng is ln γ at x, given the parameters' values p, into out, with s
	// scratch room: the formula on f64.
	lng(p []f64, x []f64, mut s []f64, mut out []f64)
	// lng_jet is the same formula on jets.
	lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet)
}

// Want asks eval for derivatives beyond ln γ.
@[flag]
pub enum Want {
	t // ∂ln γᵢ/∂T and ∂²ln γᵢ/∂T², at constant composition
	n // n·∂ln γᵢ/∂nⱼ, at constant T
}

// ln_gamma_only asks eval for ln γ alone, as a flash's substitution steps
// need it. A flag enum has no zero member, and a number cast to an enum
// needs unsafe, so it's named once here.
pub const ln_gamma_only = unsafe { Want(0) }

// Work is a model's workspace: what eval gives, and what it keeps.
pub struct Work {
pub mut:
	ln_gamma []f64 // ln γᵢ
	dt       []f64 // ∂ln γᵢ/∂T, 1/K
	dtt      []f64 // ∂²ln γᵢ/∂T², 1/K²
	// dn is n·∂ln γᵢ/∂nⱼ, nc×nc row-major, dimensionless. Gibbs–Duhem makes
	// Σᵢ xᵢ·dn[i][j] = 0, and it's symmetric, being n times the Hessian of
	// nGᴱ/RT in the mole numbers.
	dn []f64
mut:
	t  f64 // the T the parameters are at; 0 before any, and after a refusal
	pv []f64
	pj []numeric.Jet // the parameters with their T-derivatives
	pc []numeric.Jet // their values as jets that don't move
	sf []f64
	sj []numeric.Jet
	xj []numeric.Jet
	oj []numeric.Jet
	dx []f64 // ∂ln γᵢ/∂xⱼ of the formula, the x's taken as independent
}

// new_work is a workspace for nc components, np parameters and ns of
// scratch.
fn new_work(nc int, np int, ns int) Work {
	return Work{
		ln_gamma: []f64{len: nc}
		dt:       []f64{len: nc}
		dtt:      []f64{len: nc}
		dn:       []f64{len: nc * nc}
		pv:       []f64{len: np}
		pj:       []numeric.Jet{len: np}
		pc:       []numeric.Jet{len: np}
		sf:       []f64{len: ns}
		sj:       []numeric.Jet{len: ns}
		xj:       []numeric.Jet{len: nc}
		oj:       []numeric.Jet{len: nc}
		dx:       []f64{len: nc * nc}
	}
}

// ActivityError is a model refusing what it was given.
pub struct ActivityError {
	Error
pub:
	model  string
	reason string
}

pub fn (e ActivityError) msg() string {
	return '${e.model}: ${e.reason}'
}

fn refuse(m Model, reason string) IError {
	return ActivityError{
		model:  m.name()
		reason: reason
	}
}

// eval is m's ln γ at (t, x) into w.ln_gamma, and with want its
// T-derivatives (w.dt, w.dtt) and composition derivatives (w.dn). x is mole
// fractions, normalised by the caller, and refused where they don't sum to
// 1; a component may be absent, at 0.
//
// It answers or refuses: what it was asked for is finite, or it returns an
// ActivityError, and then w holds no answer.
//
// The composition derivatives come from the formula's derivatives with the
// x's taken as independent, ∂ᵢⱼ: since xₖ = nₖ/n, n·∂ln γᵢ/∂nⱼ = ∂ᵢⱼ −
// Σₖ xₖ·∂ᵢₖ, whatever the formula does off the simplex.
pub fn eval(m Model, t f64, x []f64, want Want, mut w Work) ! {
	nc := m.size()
	if !(t > 0.0) || !math.is_finite(t) {
		return refuse(m, 'a temperature of ${t} K')
	}
	if x.len != nc {
		return refuse(m, '${x.len} mole fractions for ${nc} components')
	}
	mut total := 0.0
	for v in x {
		if !(v >= 0.0) || !math.is_finite(v) {
			return refuse(m, 'a mole fraction of ${v}')
		}
		total += v
	}
	// All zero among them. UNIQUAC, UNIFAC and Wilson would answer an
	// unnormalised x wrongly, and finitely.
	if !(math.abs(total - 1.0) <= sum_slack) {
		return refuse(m, 'mole fractions that sum to ${total}, not to 1')
	}
	if w.t != t {
		// Not kept until every parameter is in: params writes them as it
		// goes, so a refusal part of the way through leaves none of a T's.
		w.t = 0.0
		m.params(t, mut w)!
		for k, v in w.pv {
			w.pc[k] = numeric.constant(v)
		}
		w.t = t
	}
	m.lng(w.pv, x, mut w.sf, mut w.ln_gamma)
	finite_answer(m, w.ln_gamma, 'ln γ', t, x)!
	if want.has(.t) {
		for i in 0 .. nc {
			w.xj[i] = numeric.constant(x[i])
		}
		m.lng_jet(w.pj, w.xj, mut w.sj, mut w.oj)
		for i in 0 .. nc {
			w.dt[i] = w.oj[i].d
			w.dtt[i] = w.oj[i].dd
		}
		finite_answer(m, w.dt, '∂ln γ/∂T', t, x)!
		finite_answer(m, w.dtt, '∂²ln γ/∂T²', t, x)!
	}
	if want.has(.n) {
		for j in 0 .. nc {
			for i in 0 .. nc {
				w.xj[i] = if i == j { numeric.along(x[i]) } else { numeric.constant(x[i]) }
			}
			m.lng_jet(w.pc, w.xj, mut w.sj, mut w.oj)
			for i in 0 .. nc {
				w.dx[i * nc + j] = w.oj[i].d
			}
		}
		for i in 0 .. nc {
			mut s := 0.0
			for k in 0 .. nc {
				s += x[k] * w.dx[i * nc + k]
			}
			for j in 0 .. nc {
				w.dn[i * nc + j] = w.dx[i * nc + j] - s
			}
		}
		finite_answer(m, w.dn, 'n·∂ln γ/∂n', t, x)!
	}
}

// finite_answer refuses an answer that isn't finite, naming its component:
// values is one for each, or a row of x.len for each.
fn finite_answer(m Model, values []f64, what string, t f64, x []f64) ! {
	for k, v in values {
		if !math.is_finite(v) {
			return refuse(m, '${what} of component ${k * x.len / values.len} is ${v} at ${t} K and x = ${x}')
		}
	}
}

// Excess is a mixture's excess properties, per mole, from its ln γ.
pub struct Excess {
pub:
	g  f64 // Gᴱ, J/mol: RT·Σ xᵢ ln γᵢ
	h  f64 // Hᴱ, J/mol: −RT²·Σ xᵢ ∂ln γᵢ/∂T
	s  f64 // Sᴱ, J/(mol·K): (Hᴱ − Gᴱ)/T
	cp f64 // Cpᴱ, J/(mol·K): ∂Hᴱ/∂T
}

// excess is the excess properties at (t, x) from a Work that eval filled
// with .t there.
pub fn excess(t f64, x []f64, w &Work) Excess {
	mut g, mut d1, mut d2 := 0.0, 0.0, 0.0
	for i, v in x {
		if v == 0.0 {
			continue
		}
		g += v * w.ln_gamma[i]
		d1 += v * w.dt[i]
		d2 += v * w.dtt[i]
	}
	rt := r_gas * t
	h := -rt * t * d1
	return Excess{
		g:  rt * g
		h:  h
		s:  (h - rt * g) / t
		cp: -r_gas * t * (2.0 * d1 + t * d2)
	}
}

// Param is a parameter that depends on T: an expression in T, dimensionless,
// as a file or a bank gives it (D7). A number is the expression of one.
pub struct Param {
pub:
	e    expr.Expr
	defs expr.Definitions
}

// temperature is the one variable a parameter's expression may use, T.
pub fn temperature() expr.Context {
	return expr.Context{
		variables: [expr.Variable{'T', quantity.Dimension.temperature.dim()}]
	}
}

// param reads a parameter's text, with a file's definitions, refusing one
// that isn't dimensionless: what it is says what it's for, in a message.
pub fn param(text string, defs expr.Definitions, what string) !Param {
	e := expr.parse(temperature(), &defs, text)!
	e.must_be(quantity.Dim{}, what)!
	return Param{
		e:    e
		defs: defs
	}
}

// plain is param with no definitions: a bank's parameter.
pub fn plain(text string, what string) !Param {
	return param(text, expr.Definitions{}, what)
}

// at is the parameter at t with its first two T-derivatives.
pub fn (p &Param) at(t f64) numeric.Jet {
	return p.e.jet(&p.defs, [numeric.along(t)])
}

// text is the parameter as written.
pub fn (p &Param) text() string {
	return p.e.text
}

// finite_param is p at t, refused where it isn't finite.
fn finite_param(m Model, p &Param, t f64, what string) !numeric.Jet {
	j := p.at(t)
	if !math.is_finite(j.v) || !math.is_finite(j.d) || !math.is_finite(j.dd) {
		return refuse(m, '${what}, ${p.text()}, is ${j.v} at ${t} K')
	}
	return j
}

// usable is whether j can be kept as a parameter that a formula divides by
// or takes the log of: a positive normal double, with finite T-derivatives.
// The second derivative's check is the value's and the first's too: j is
// either finite_param's, all three finite, or an exp, eᵃ, whose second,
// (eᵃa′)a′ + eᵃa″, isn't finite where eᵃ or eᵃa′ isn't (∞·0 being NaN).
fn usable(j numeric.Jet) bool {
	return j.v >= least_normal && math.is_finite(j.dd)
}

// unusable is the refusal of such a parameter, named by what. (It's only
// called to refuse, so that a fill that succeeds builds no message.)
fn unusable(m Model, what string, j numeric.Jet, t f64) IError {
	return refuse(m, '${what} is ${j.v} at ${t} K, with T-derivatives ${j.d} and ${j.dd}: it must be a positive normal double, and its derivatives finite')
}
