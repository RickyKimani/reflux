module flash

import math
import numeric
import thermo

// VfParams are the vapour-fraction flashes' options.
@[params]
pub struct VfParams {
pub:
	// Converged when every |ln fᵢⱽ − ln fᵢᴸ| ≤ tol, and Σy − Σx is within
	// tol of 0, so that each component's mass balance holds to tol.
	tol f64 = 1e-12
	// verify re-tests each phase for stability, with the other as a trial,
	// as the PT flash does.
	verify    bool = true
	max_steps int  = 50
	// estimate is where to start the unknown, P in Pa for tvf or T in K for
	// pvf, instead of Wilson's estimate: with the K-values of a PT split
	// there if there is one, else Wilson's. Where there are two answers,
	// such as the two dew pressures of a gas condensate at one temperature,
	// it picks which one is found.
	estimate f64
}

// tvf finds the pressure at which a feed z at temperature T has vapour
// fraction β. β = 0 is the bubble point and β = 1 the dew point.
//
// The unknowns are ln Kᵢ for each component present and ln P. The equations
// are equal fugacities and Rachford–Rice at the given β:
//
//     ln Kᵢ + ln φᵢ(y) − ln φᵢ(x) = 0,   Σᵢ (yᵢ − xᵢ) = 0,
//     xᵢ = zᵢ/(1 + β(Kᵢ − 1)),   yᵢ = Kᵢxᵢ.
//
// It runs Newton's method with the analytic Jacobian (Michelsen & Mollerup,
// chapter 12) from Wilson's K-values, at the pressure where they satisfy
// Rachford–Rice. If that fails, it starts again from the K-values of the
// nearest two-phase PT flash it finds at that temperature.
// (Successive substitution first, as for a PT flash, slid to the trivial
// solution K = 1 from Wilson's poor starts near a critical point.) Each phase must end on its volume root of lower Gibbs
// energy, and phase y (fraction β) must come out the less dense by mass
// (D6).
//
// Where two pressures fit, as at a temperature between a mixture's critical
// point and its cricondentherm, it returns the one Newton reaches from its
// start: Wilson's estimate, which usually leads to the lower, or
// params.estimate. A pure component's answer is its saturation pressure,
// whatever β.
//
// The answer is a FlashResult with both phases, the vapour first, even when
// one has no moles (β = 0 or 1). Failure is a FlashError: the iteration
// collapsing onto the feed (K = 1), as it does where no such point exists,
// above the cricondenbar or cricondentherm, say; the phases coming out the
// wrong way round (.phases_swapped); or no convergence.
pub fn tvf(m thermo.Model, t f64, beta f64, z []f64, params VfParams) !FlashResult {
	if !(t > 0.0) || !math.is_finite(t) {
		return fail(.bad_input, 'T = ${t} K is not a temperature')
	}
	return saturate(m, .t, t, beta, z, params)
}

// pvf finds the temperature at which a feed z at pressure P has vapour
// fraction β: as tvf, with ln T as the unknown in place of ln P. Where two
// temperatures fit, as at a pressure between a mixture's critical point and
// its cricondenbar, Wilson's estimate usually leads to the higher.
pub fn pvf(m thermo.Model, p f64, beta f64, z []f64, params VfParams) !FlashResult {
	if !(p > 0.0) || !math.is_finite(p) {
		return fail(.bad_input, 'P = ${p} Pa is not a pressure')
	}
	return saturate(m, .p, p, beta, z, params)
}

// bubble_p is the bubble-point pressure at T: tvf with β = 0.
pub fn bubble_p(m thermo.Model, t f64, z []f64, params VfParams) !FlashResult {
	return tvf(m, t, 0.0, z, params)
}

// dew_p is the dew-point pressure at T: tvf with β = 1.
pub fn dew_p(m thermo.Model, t f64, z []f64, params VfParams) !FlashResult {
	return tvf(m, t, 1.0, z, params)
}

// bubble_t is the bubble-point temperature at P: pvf with β = 0.
pub fn bubble_t(m thermo.Model, p f64, z []f64, params VfParams) !FlashResult {
	return pvf(m, p, 0.0, z, params)
}

// dew_t is the dew-point temperature at P: pvf with β = 1.
pub fn dew_t(m thermo.Model, p f64, z []f64, params VfParams) !FlashResult {
	return pvf(m, p, 1.0, z, params)
}

// Given says which of T and P a vapour-fraction flash was given; the other
// is the unknown.
enum Given {
	t
	p
}

// Saturation is one vapour-fraction flash.
struct Saturation {
	m       thermo.Model
	given   Given
	value   f64 // the given T in K, or P in Pa
	beta    f64
	z       []f64
	present []int
	params  VfParams
}

fn saturate(m thermo.Model, given Given, value f64, beta f64, z []f64, params VfParams) !FlashResult {
	if !(beta >= 0.0 && beta <= 1.0) {
		return fail(.bad_input, 'the vapour fraction β = ${beta} is outside [0, 1]')
	}
	if !(params.tol > 0.0) || params.max_steps < 1 {
		return fail(.bad_input, 'tol and max_steps must be positive')
	}
	if !(params.estimate >= 0.0) || !math.is_finite(params.estimate) {
		return fail(.bad_input, 'estimate = ${params.estimate} must be a positive T or P, or 0 to start from Wilson')
	}
	check_feed(m, z)!
	feed := normalised(z)
	s := Saturation{
		m:       m
		given:   given
		value:   value
		beta:    beta
		z:       feed
		present: []int{len: feed.len, init: index}.filter(feed[it] > 0.0)
		params:  params
	}
	if s.present.len == 1 {
		return s.pure()!
	}
	return s.mixture()!
}

// tp turns the unknown τ = ln P or ln T into (T, P), refusing values no
// model is meant for.
fn (s Saturation) tp(tau f64) !(f64, f64) {
	x := math.exp(tau)
	t, p := if s.given == .t { s.value, x } else { x, s.value }
	if !(t >= 1.0 && t <= 1e5) || !(p >= 1e-6 && p <= 1e12) {
		return error('T = ${t} K, P = ${p} Pa is out of range')
	}
	return t, p
}

// at is the state of the iteration at v = (ln K over the components present,
// τ): both phases' compositions, unnormalised as Rachford–Rice at the given
// β makes them, and solved, x on its liquid root and y on its vapour root.
// (With each on its lower-Gibbs root instead, an incipient phase can start
// on the feed's kind of root, and the iteration slides to the trivial
// solution K = 1; mixture checks the roots at the end.)
struct SatPoint {
	t     f64
	p     f64
	d     []f64 // 1 + β(Kᵢ − 1)
	x     []f64
	y     []f64
	sum_x f64
	sum_y f64
	sx    thermo.State // at x/Σx
	sy    thermo.State // at y/Σy
}

fn (s Saturation) at(v []f64) !SatPoint {
	k := s.present.len
	nc := s.z.len
	t, p := s.tp(v[k])!
	mut d := []f64{len: nc, init: 1.0}
	mut x := []f64{len: nc}
	mut y := []f64{len: nc}
	for j, i in s.present {
		kk := math.exp(v[j])
		if !(kk > 0.0) || !math.is_finite(kk) {
			return error('a K-value left (0, ∞)')
		}
		d[i] = 1.0 + s.beta * (kk - 1.0)
		x[i] = s.z[i] / d[i]
		y[i] = kk * x[i]
	}
	sum_x, sum_y := numeric.sum(x), numeric.sum(y)
	sx := thermo.solve_tp(s.m, t, p, x.map(it / sum_x), .liquid)!
	sy := thermo.solve_tp(s.m, t, p, y.map(it / sum_y), .vapor)!
	return SatPoint{t, p, d, x, y, sum_x, sum_y, sx, sy}
}

// tau_slope is ∂ln φᵢ/∂τ at fixed composition: P·∂ln φᵢ/∂P, or T·∂ln φᵢ/∂T.
fn (s Saturation) tau_slope(d thermo.LnPhiDerivs, i int, cur SatPoint) f64 {
	return if s.given == .t { cur.p * d.dp[i] } else { cur.t * d.dt[i] }
}

// eval is the system newton solves: F and its Jacobian at v.
fn (s Saturation) eval(v []f64, mut f []f64, mut jac []f64) ! {
	k := s.present.len
	n := k + 1
	nc := s.z.len
	cur := s.at(v)!
	dx := thermo.ln_phi_derivs(s.m, cur.sx)
	dy := thermo.ln_phi_derivs(s.m, cur.sy)
	mut h := 0.0
	for j, i in s.present {
		f[j] = v[j] + cur.sy.ln_phi[i] - cur.sx.ln_phi[i]
		h += cur.y[i] - cur.x[i]
		for l, q in s.present {
			// ∂yq/∂ln Kq = yq(1 − β)/dq and ∂xq/∂ln Kq = −βKqxq/dq; ln φ
			// depends on mole numbers through the mole fractions, so its
			// derivative there is Φᵢq/Σ, with Φ = n·∂ln φᵢ/∂nq at one mole.
			kq := math.exp(v[l])
			dyq := cur.y[q] * (1.0 - s.beta) / cur.d[q]
			dxq := -s.beta * kq * cur.x[q] / cur.d[q]
			mut a := dy.dn[i * nc + q] / cur.sum_y * dyq - dx.dn[i * nc + q] / cur.sum_x * dxq
			if j == l {
				a += 1.0
			}
			jac[j * n + l] = a
		}
		jac[j * n + k] = s.tau_slope(dy, i, cur) - s.tau_slope(dx, i, cur)
	}
	f[k] = h
	for l, q in s.present {
		kq := math.exp(v[l])
		jac[k * n + l] = s.z[q] * kq / (cur.d[q] * cur.d[q])
	}
	jac[k * n + k] = 0.0
}

// WilsonGap is Rachford–Rice's sum at β with Wilson's K-values: in ln P at
// a given T, or in 1/T at a given P, for start's bisection. The solvers take
// structs, not closures: see numeric/functions.v.
struct WilsonGap {
	s     Saturation
	comps []thermo.Component
	t     f64 // the given T, or 0 to search in 1/T at p
	p     f64
}

fn (g WilsonGap) at(x f64) f64 {
	t, p := if g.t > 0.0 { g.t, math.exp(x) } else { 1.0 / x, g.p }
	k := wilson_k(g.comps, t, p)
	mut sum := 0.0
	for i in g.s.present {
		sum += g.s.z[i] * (k[i] - 1.0) / (1.0 + g.s.beta * (k[i] - 1.0))
	}
	return sum
}

// start is Wilson's K-values at the T or P where they satisfy
// Rachford–Rice at β, or at params.estimate.
fn (s Saturation) start() ![]f64 {
	comps := s.m.components()
	mut t, mut p := 0.0, 0.0
	if s.params.estimate > 0.0 {
		t, p = if s.given == .t { s.value, s.params.estimate } else { s.params.estimate, s.value }
		// A split there has far better K-values than Wilson's.
		if r := pt(s.m, t, p, s.z, verify: false) {
			if r.phases.len == 2 {
				y, x := r.phases[0].state.n, r.phases[1].state.n
				mut v := s.present.map(math.log(y[it] / x[it]))
				v << math.log(if s.given == .t { p } else { t })
				return v
			}
		}
	} else if s.given == .t {
		// Every Kᵢ = aᵢ/P, so the sum falls as P rises: bisect in ln P
		// between where all K > 1 and where all K < 1.
		t = s.value
		a := wilson_k(comps, t, 1.0)
		mut amin, mut amax := math.inf(1), 0.0
		for i in s.present {
			amin = math.min(amin, a[i])
			amax = math.max(amax, a[i])
		}
		lo, hi := math.log(amin / 2.0), math.log(amax * 2.0)
		lnp := numeric.bisect(WilsonGap{
			s:     s
			comps: comps
			t:     t
		}, lo, hi) or { return fail(.bad_input, 'Wilson finds no pressure at ${t} K: ${err}') }
		p = math.exp(lnp)
	} else {
		// Every ln Kᵢ is linear in 1/T and rises with T: bisect in 1/T.
		p = s.value
		inv := numeric.bisect(WilsonGap{
			s:     s
			comps: comps
			p:     p
		}, 1e-5, 1.0) or { return fail(.bad_input, 'Wilson finds no temperature at ${p} Pa: ${err}') }
		t = 1.0 / inv
	}
	k := wilson_k(comps, t, p)
	mut v := s.present.map(math.log(k[it]))
	v << math.log(if s.given == .t { p } else { t })
	return v
}

// pt_step is pt_start's step in τ: 10% in P, 2% in T.
fn (s Saturation) pt_step() f64 {
	return if s.given == .t { 0.1 } else { 0.02 }
}

// pt_start is a second start, for when Newton fails from Wilson's: the
// K-values of a two-phase PT flash at the given T or P. It searches out from
// Wilson's estimate of the unknown (or params.estimate), a step at a time each way, nearest
// first, up to 20 steps (P by 10%, T by 2%), for one that splits. (Which
// way the two-phase region lies can't be read off a single phase's label:
// a dense fluid past the critical point can be called a liquid either side.)
fn (s Saturation) pt_start(v0 []f64) ?[]f64 {
	k := s.present.len
	for i in 1 .. 21 {
		for dir in [1.0, -1.0] {
			tau := v0[k] + dir * f64(i) * s.pt_step()
			t, p := s.tp(tau) or { continue }
			r := pt(s.m, t, p, s.z, verify: false) or { continue }
			if r.phases.len == 2 {
				y, x := r.phases[0].state.n, r.phases[1].state.n
				mut v := s.present.map(math.log(y[it] / x[it]))
				v << tau
				return v
			}
		}
	}
	return none
}

// mixture runs Newton from Wilson's start and, if that fails, from
// pt_start's. A failure from both reports both. Where Newton failed and no
// PT flash in pt_start's range splits, the point is very likely outside the
// two-phase envelope.
fn (s Saturation) mixture() !FlashResult {
	v0 := s.start()!
	return s.attempt(v0) or {
		wilson := err
		v1 := s.pt_start(v0) or {
			if wilson is FlashError && wilson.kind == .no_convergence {
				k := s.present.len
				lo, hi := math.exp(v0[k] - 20.0 * s.pt_step()), math.exp(v0[k] + 20.0 * s.pt_step())
				range := if s.given == .t {
					'P from ${lo:.4g} to ${hi:.4g} Pa'
				} else {
					'T from ${lo:.4g} to ${hi:.4g} K'
				}
				return fail(.outside_envelope, '${s.what()}: Newton failed from Wilson K-values (${wilson.reason.all_after('${s.what()}: ')}), and no PT flash splits for ${range}. There is very likely no such point: above the cricondenbar or past the cricondentherm, say')
			}
			return wilson
		}
		return s.attempt(v1) or {
			if err is FlashError {
				at := math.exp(v1[s.present.len])
				unit := if s.given == .t { 'Pa' } else { 'K' }
				return fail(err.kind, '${err.reason}, starting from a PT split at ${at} ${unit} (and from Wilson K-values: ${wilson.msg()}). The two-phase region reaches here, but not, it seems, a point of this kind: near a critical point both ends of an isobar can be bubble points')
			}
			return err
		}
	}
}

// attempt is Newton from v0, and every check on where it ends.
fn (s Saturation) attempt(v0 []f64) !FlashResult {
	k := s.present.len
	// eval makes Saturation a numeric.System.
	sol := numeric.newton(s, v0, ftol: s.params.tol, max_iter: s.params.max_steps) or {
		return fail(.no_convergence, '${s.what()}: ${err.msg()}')
	}
	v := sol.x
	mut collapsed := true
	for j in 0 .. k {
		if math.abs(v[j]) >= split_trivial_tol {
			collapsed = false
		}
	}
	if collapsed {
		return fail(.no_convergence, '${s.what()}: the iteration collapsed onto the feed (every |ln K| < ${split_trivial_tol}). There may be no such point, as above the cricondenbar or cricondentherm, or it may be close to a critical point')
	}
	cur := s.at(v) or { return fail(.no_volume_root, '${s.what()}: ${err}') }
	x := cur.x.map(it / cur.sum_x)
	y := cur.y.map(it / cur.sum_y)
	s.on_stable_roots(cur)!
	mut worst := 0.0
	for i in s.present {
		worst = math.max(worst, math.abs(math.log(y[i]) + cur.sy.ln_phi[i] - math.log(x[i]) - cur.sx.ln_phi[i]))
	}
	if worst > s.params.tol || math.abs(cur.sum_y - cur.sum_x) > s.params.tol {
		return fail(.no_convergence, '${s.what()}: Newton stopped at a fugacity residual of ${worst}, above tol = ${s.params.tol}')
	}
	return s.result(cur.t, cur.p, x, cur.sx, y, cur.sy, sol.iterations)!
}

// on_stable_roots checks that each phase is on its root of lower Gibbs
// energy. Otherwise one of them is metastable and this isn't an
// equilibrium. (Where there's a single root, solving from the liquid and
// the vapour side can land an ulp or so apart, while distinct roots differ
// by much more.)
fn (s Saturation) on_stable_roots(cur SatPoint) ! {
	for st in [cur.sx, cur.sy] {
		stable := thermo.solve_tp(s.m, cur.t, cur.p, st.n, .stable) or {
			return fail(.no_volume_root, '${s.what()}: ${err}')
		}
		if math.abs(stable.v - st.v) > 1e-6 * st.v {
			return fail(.no_convergence, '${s.what()}: at ${cur.t} K and ${cur.p} Pa a phase sits on its volume root of higher Gibbs energy, so this boundary is between metastable phases')
		}
	}
}

// what describes the flash for an error message.
fn (s Saturation) what() string {
	name := if s.beta == 0.0 {
		'the bubble point'
	} else if s.beta == 1.0 {
		'the dew point'
	} else {
		'vapour fraction ${s.beta}'
	}
	return if s.given == .t { '${name} at ${s.value} K' } else { '${name} at ${s.value} Pa' }
}

// pure finds a pure component's saturation point: ln φ equal on the liquid
// and vapour roots, by Newton in τ. It first moves τ until both roots
// exist: with the vapour root alone P is too low (or T too high), with the
// liquid root alone the reverse.
// roots is the pure component's liquid and vapour roots at τ = ln P or ln T.
fn (s Saturation) roots(tau f64) !(thermo.State, thermo.State) {
	t, p := s.tp(tau)!
	liq := thermo.solve_tp(s.m, t, p, s.z, .liquid)!
	vap := thermo.solve_tp(s.m, t, p, s.z, .vapor)!
	return liq, vap
}

// PureSat is ln φᵢᴸ = ln φᵢⱽ in τ, for component i alone, as a
// numeric.System.
struct PureSat {
	s Saturation
	i int
}

fn (q PureSat) eval(x []f64, mut f []f64, mut jac []f64) ! {
	s, i := q.s, q.i
	liq, vap := s.roots(x[0])!
	if !(vap.v - liq.v > 1e-9 * vap.v) {
		return error('only one volume root')
	}
	dl := thermo.ln_phi_derivs(s.m, liq)
	dv := thermo.ln_phi_derivs(s.m, vap)
	f[0] = liq.ln_phi[i] - vap.ln_phi[i]
	jac[0] = if s.given == .t {
		liq.p * (dl.dp[i] - dv.dp[i])
	} else {
		liq.t * (dl.dt[i] - dv.dt[i])
	}
}

fn (s Saturation) pure() !FlashResult {
	i := s.present[0]
	comps := s.m.components()
	c := comps[i]
	mut tau := 0.0
	if s.params.estimate > 0.0 {
		tau = math.log(s.params.estimate)
	} else if s.given == .t {
		tau = math.log(wilson_k([c], s.value, 1.0)[0])
	} else {
		// Wilson's K = 1 at P: 1/T = (1 − ln(P/Pc)/(5.373(1 + ω)))/Tc.
		u := (1.0 - math.log(s.value / c.pc) / (5.373 * (1.0 + c.omega))) / c.tc
		if !(u > 0.0) {
			return fail(.bad_input, '${s.what()}: P is far above the critical pressure')
		}
		tau = math.log(1.0 / u)
	}
	mut step, mut last_dir := 0.05, 0.0
	mut found := false
	for _ in 0 .. 200 {
		liq, vap := s.roots(tau) or { return fail(.no_volume_root, '${s.what()}: ${err}') }
		if vap.v - liq.v > 1e-9 * vap.v {
			found = true
			break
		}
		// One root: raise P (lower T) if it's vapour-like, else the reverse.
		up := thermo.phase_label(s.m, liq) == .vapor
		dir := if up == (s.given == .t) { 1.0 } else { -1.0 }
		if last_dir != 0.0 && dir != last_dir {
			step /= 2.0
		}
		last_dir = dir
		tau += dir * step
	}
	if !found {
		return fail(.no_convergence, '${s.what()}: found no state with both a liquid and a vapour root; past the critical point?')
	}
	sol := numeric.newton(PureSat{s, i}, [tau], ftol: s.params.tol, max_iter: s.params.max_steps) or {
		return fail(.no_convergence, '${s.what()}: ${err.msg()}')
	}
	liq, vap := s.roots(sol.x[0]) or { return fail(.no_volume_root, '${s.what()}: ${err}') }
	stab := stability(s.m, liq.t, liq.p, s.z)!
	return FlashResult{
		t:          liq.t
		p:          liq.p
		z:          s.z
		beta:       s.beta
		phases:     [Phase{.vapor, vap, s.beta}, Phase{.liquid, liq, 1.0 - s.beta}]
		path:       .newton
		iterations: sol.iterations
		stability:  stab
	}
}

// result checks the converged point and builds the FlashResult: phase y,
// of fraction β, must be the less dense, and with params.verify neither
// phase may be unstable. stability is the feed's test at the answer: for a
// bubble or dew point that's the feed phase's, with the incipient phase as
// a trial.
fn (s Saturation) result(t f64, p f64, x []f64, sx thermo.State, y []f64, sy thermo.State, iterations int) !FlashResult {
	if !(mass_density(s.m, sy) < mass_density(s.m, sx)) {
		return fail(.phases_swapped, '${s.what()}: the iteration reached a boundary at ${t} K and ${p} Pa where the phase of fraction β is the denser, so its vapour fraction is ${1.0 - s.beta}')
	}
	mut feed_test := Stability{}
	mut have_feed_test := false
	if s.params.verify {
		for k, pair in [[x, y], [y, x]] {
			check := stability(s.m, t, p, pair[0], trials: [pair[1]], tol: verify_tol(s.params.tol))!
			if !check.stable {
				tm := if inc := check.incipient() { inc.tm } else { 0.0 }
				return fail(.more_phases, '${s.what()}: a phase at ${t} K and ${p} Pa is unstable (tm = ${tm}): three or more phases, or another split')
			}
			// The feed is x at a bubble point and y at a dew point.
			if (s.beta == 0.0 && k == 0) || (s.beta == 1.0 && k == 1) {
				feed_test = check
				have_feed_test = true
			}
		}
	}
	if !have_feed_test {
		other := if s.beta == 0.0 { y } else { x }
		feed_test = stability(s.m, t, p, s.z, trials: [other], tol: verify_tol(s.params.tol))!
	}
	return FlashResult{
		t:          t
		p:          p
		z:          s.z
		beta:       s.beta
		phases:     [Phase{.vapor, sy, s.beta}, Phase{.liquid, sx, 1.0 - s.beta}]
		path:       .newton
		iterations: iterations
		stability:  feed_test
	}
}
