module flash

import math
import numeric
import thermo

// azeotrope_t finds the azeotrope near x0 at pressure P: the T and the
// composition x at which the vapour in equilibrium with the liquid is x
// itself (M6 step 6, criterion 5). azeotrope_p is the same at a given T.
//
// A vapour-fraction flash can't find one, since its liquid and vapour are
// the same there, which it takes for the trivial solution. Here the
// unknowns are x over the components present but the last, which is
// 1 − Σ of the rest, and τ = ln T (or ln P), and the equations are equal
// fugacities with y = x:
//
//     gᵢ = ln φᵢᴸ(x) − ln φᵢⱽ(x) = 0   for each component present,
//
// as many as the unknowns. Newton runs on them with the analytic Jacobian:
// ∂gᵢ/∂xₖ = (Φᴸᵢₖ − Φᴸᵢₗ) − (Φⱽᵢₖ − Φⱽᵢₗ), Φ = n·∂ln φᵢ/∂nⱼ and l the last
// component, and ∂gᵢ/∂τ = T(∂ln φᵢᴸ/∂T − ∂ln φᵢⱽ/∂T) (or P·∂/∂P). It starts
// from x0's bubble point, or one 2% from it where x0 is an azeotrope
// already. The answer is shown: each |gᵢ| ≤ tol, every xᵢ in (0, 1), the
// liquid and the vapour two phases, not one root twice, and the liquid
// stable with the vapour as a trial. A FlashResult gives the vapour and the
// liquid, at β = 0, each of composition x.
pub fn azeotrope_t(model thermo.PhaseModel, p f64, x0 []f64, params VfParams) !FlashResult {
	if !(p > 0.0) || !math.is_finite(p) {
		return fail(.bad_input, 'P = ${p} Pa is not a pressure')
	}
	return find_azeotrope(model, .p, p, x0, params)
}

// azeotrope_p is azeotrope_t at a given temperature, finding P.
pub fn azeotrope_p(model thermo.PhaseModel, t f64, x0 []f64, params VfParams) !FlashResult {
	if !(t > 0.0) || !math.is_finite(t) {
		return fail(.bad_input, 'T = ${t} K is not a temperature')
	}
	return find_azeotrope(model, .t, t, x0, params)
}

// Azeotrope is one search, a numeric.System in (x, τ).
struct Azeotrope {
	m       thermo.PhaseModel
	given   Given
	value   f64
	nc      int
	present []int
	params  VfParams
}

fn find_azeotrope(model thermo.PhaseModel, given Given, value f64, x0 []f64, params VfParams) !FlashResult {
	if !(params.tol > 0.0) || params.max_steps < 1 {
		return fail(.bad_input, 'tol and max_steps must be positive')
	}
	m := thermo.for_one_thread(model)
	check_feed(m, x0)!
	x := normalised(x0)
	present := []int{len: x.len, init: index}.filter(x[it] > 0.0)
	if present.len < 2 {
		return fail(.bad_input, 'an azeotrope needs two components or more, and x0 has ${present.len}')
	}
	a := Azeotrope{
		m:       m
		given:   given
		value:   value
		nc:      x.len
		present: present
		params:  params
	}
	start := saturate(m, given, value, 0.0, x, params) or { a.nearby(x, err)! }
	mut v := present[..present.len - 1].map(x[it])
	v << math.log(if given == .t { start.p } else { start.t })
	sol := numeric.newton(a, v, ftol: params.tol, max_iter: params.max_steps) or {
		return fail(.no_convergence, '${a.what()}: ${err.msg()}')
	}
	return a.result(sol.x, sol.iterations)!
}

// nearby is the start where x0's own bubble point fails. At an azeotrope,
// or within about 1e-4 of one, the bubble point collapses onto the feed,
// since K = 1 is the trivial solution it refuses; 2% from x0 towards each
// component in turn, ln K is far from that, and T or P there is close.
fn (a Azeotrope) nearby(x []f64, first IError) !FlashResult {
	for i in a.present {
		mut near := x.map(0.98 * it)
		near[i] += 0.02
		if r := saturate(a.m, a.given, a.value, 0.0, near, a.params) {
			return r
		}
	}
	return fail(.no_convergence, '${a.what()}: no bubble point at x0 or 2% from it to start from: ${first.msg()}')
}

fn (a Azeotrope) what() string {
	return if a.given == .t {
		'the azeotrope at ${a.value} K'
	} else {
		'the azeotrope at ${a.value} Pa'
	}
}

// state is (x, T, P) at v.
fn (a Azeotrope) state(v []f64) !([]f64, f64, f64) {
	k := a.present.len
	mut x := []f64{len: a.nc}
	mut last := 1.0
	for j in 0 .. k - 1 {
		xi := v[j]
		if !(xi > 0.0 && xi < 1.0) {
			return error('a mole fraction left (0, 1)')
		}
		x[a.present[j]] = xi
		last -= xi
	}
	if !(last > 0.0) {
		return error('a mole fraction left (0, 1)')
	}
	x[a.present[k - 1]] = last
	tau := math.exp(v[k - 1])
	t, p := if a.given == .t { a.value, tau } else { tau, a.value }
	if !(t >= 1.0 && t <= 1e5) || !(p >= 1e-6 && p <= 1e12) {
		return error('T = ${t} K, P = ${p} Pa is out of range')
	}
	return x, t, p
}

// eval makes Azeotrope a numeric.System.
fn (a Azeotrope) eval(v []f64, mut f []f64, mut jac []f64) ! {
	x, t, p := a.state(v)!
	liq := a.m.phase(t, p, x, .liquid)!
	vap := a.m.phase(t, p, x, .vapor)!
	dl := a.m.ln_phi_derivs(liq)!
	dv := a.m.ln_phi_derivs(vap)!
	k := a.present.len
	nc := a.nc
	l := a.present[k - 1]
	for r, i in a.present {
		f[r] = liq.ln_phi[i] - vap.ln_phi[i]
		for c in 0 .. k - 1 {
			j := a.present[c]
			jac[r * k + c] = (dl.dn[i * nc + j] - dl.dn[i * nc + l]) - (dv.dn[i * nc + j] - dv.dn[i * nc + l])
		}
		jac[r * k + k - 1] = if a.given == .t {
			p * (dl.dp[i] - dv.dp[i])
		} else {
			t * (dl.dt[i] - dv.dt[i])
		}
	}
}

fn (a Azeotrope) result(v []f64, iterations int) !FlashResult {
	x, t, p := a.state(v) or { return fail(.no_convergence, '${a.what()}: ${err.msg()}') }
	liq := a.m.phase(t, p, x, .liquid)!
	vap := a.m.phase(t, p, x, .vapor)!
	mut worst := 0.0
	for i in a.present {
		worst = math.max(worst, math.abs(liq.ln_phi[i] - vap.ln_phi[i]))
	}
	if !(worst <= a.params.tol) {
		return fail(.no_convergence, '${a.what()}: Newton stopped at a fugacity residual of ${worst}, above tol = ${a.params.tol}')
	}
	// Where an equation of state has one root, its liquid and vapour are one
	// phase, and every gᵢ is 0 at any x: the trivial solution, not an
	// azeotrope. (The two phases' G are equal at an azeotrope, so which is
	// lower says nothing.)
	if liq.v == vap.v {
		return fail(.no_convergence, '${a.what()}: at ${t} K and ${p} Pa the model has one phase, its own vapour trivially, so this is no azeotrope')
	}
	// The liquid's own stability, the vapour's composition (its own) aside:
	// a liquid that would split is a heterogeneous azeotrope, which needs
	// three phases (step 7).
	check := stability(a.m, t, p, x, tol: verify_tol(a.params.tol))!
	if !check.stable {
		tm := if inc := check.incipient() { inc.tm } else { 0.0 }
		return fail(.more_phases, '${a.what()}: the liquid at ${t} K and ${p} Pa is unstable (tm = ${tm}): a heterogeneous azeotrope, with two liquids')
	}
	return FlashResult{
		t:          t
		p:          p
		z:          x
		beta:       0.0
		phases:     [Phase{.vapor, vap, 0.0}, Phase{.liquid, liq, 1.0}]
		path:       .newton
		iterations: iterations
		stability:  check
	}
}
