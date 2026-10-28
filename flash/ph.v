module flash

import math
import numeric
import thermo

// StateParams are the PH and PS flashes' options.
@[params]
pub struct StateParams {
pub:
	// The PT flashes' fugacity tolerance.
	tol f64 = 1e-12
	// The energy balance's: converged when |H − H_spec| ≤ etol·RT, or
	// |S − S_spec| ≤ etol·R.
	etol f64 = 1e-10
	// verify re-tests each phase of a split for stability, as the PT flash
	// does.
	verify    bool = true
	max_steps int  = 100
	// estimate is a temperature to start from, in K, instead of the ideal
	// gas's estimate.
	estimate f64
}

// Spec says which state function a PH or PS flash was given.
enum Spec {
	h
	s
}

// ph flashes a feed z at pressure P to its molar enthalpy H, in J per mole
// of feed, on the Fluid's reference state (thermo.IdealGas). It finds T and
// the split.
//
// It's Michelsen's (1999) state-function flash (D5). At fixed P and H the
// equilibrium is the saddle point of
//
//     Q(T, v) = (G(T, P, v) − H_spec)/RT,
//
// a minimum in the vapour mole numbers v and a maximum in T. At fixed T,
// minimising over v is the PT flash; what's left, Q*(T), has
// dQ*/dT = −(H(T) − H_spec)/RT² with H(T) the equilibrium enthalpy, which
// never falls as T rises. So:
//
// - from the ideal gas's T for H_spec (or params.estimate), T is stepped out
//   by 50% at a time until H(T) − H_spec changes sign, with a PT flash at
//   each T;
// - Brent's method narrows that bracket to 1e-3 in T;
// - if the PT flash there splits, Newton on Q's stationary conditions in
//   (v, ln T) together, equal fugacities and H = H_spec, with Q's Hessian,
//   takes it to round-off. (Otherwise, and if that Newton fails, Brent goes
//   on to the end.)
// - The answer is a verified PT flash at the final T, held to
//   |H − H_spec| ≤ etol·RT.
//
// A pure component's H jumps at its saturation temperature, by the heat of
// vaporisation. An H_spec inside the jump is the two-phase state at Tsat
// with the vapour fraction β that the enthalpy balance gives. One outside
// it is a single phase, found on its own side of Tsat.
//
// iterations counts PT flashes and Newton steps. Failure is a FlashError:
// bad_input for an H that no temperature from 10 K to 10⁴ K reaches.
pub fn ph(fl thermo.Fluid, p f64, h f64, z []f64, params StateParams) !FlashResult {
	return state_flash(fl, .h, h, p, z, params)
}

// ps flashes a feed z at pressure P to its molar entropy S, in J/(mol·K),
// on the Fluid's reference state. It's ph's method with
// Q(T, v) = (G(T, P, v) + T·S_spec)/R, whose dQ*/dT = S_spec − S(T): the
// same saddle, with the energy balance S = S_spec held to etol·R.
pub fn ps(fl thermo.Fluid, p f64, s f64, z []f64, params StateParams) !FlashResult {
	return state_flash(fl, .s, s, p, z, params)
}

// StateFlash is one PH or PS flash.
struct StateFlash {
	fl      thermo.Fluid
	spec    Spec
	value   f64 // H_spec in J/mol or S_spec in J/(mol·K)
	p       f64
	z       []f64
	present []int
	params  StateParams
}

// The temperatures a PH or PS flash may search.
const state_t_min = 10.0
const state_t_max = 1e4

fn state_flash(fl thermo.Fluid, spec Spec, value f64, p f64, z []f64, params StateParams) !FlashResult {
	if !(p > 0.0) || !math.is_finite(p) {
		return fail(.bad_input, 'P = ${p} Pa is not a pressure')
	}
	if !math.is_finite(value) {
		return fail(.bad_input, 'the specified ${spec} = ${value} is not a number')
	}
	if !(params.tol > 0.0) || !(params.etol > 0.0) || params.max_steps < 1 {
		return fail(.bad_input, 'tol, etol and max_steps must be positive')
	}
	if !(params.estimate >= 0.0) || !math.is_finite(params.estimate) {
		return fail(.bad_input, 'estimate = ${params.estimate} must be a positive T, or 0 to start from the ideal gas')
	}
	check_feed(fl.model, z)!
	feed := normalised(z)
	s := StateFlash{
		fl:      fl
		spec:    spec
		value:   value
		p:       p
		z:       feed
		present: []int{len: feed.len, init: index}.filter(feed[it] > 0.0)
		params:  params
	}
	if s.present.len == 1 {
		return s.pure()!
	}
	return s.mixture()!
}

// what describes the flash for an error message.
fn (s StateFlash) what() string {
	return match s.spec {
		.h { 'the PH flash to H = ${s.value} J/mol at ${s.p} Pa' }
		.s { 'the PS flash to S = ${s.value} J/(mol K) at ${s.p} Pa' }
	}
}

// of is a phase's molar H or S.
fn (s StateFlash) of(st thermo.State) f64 {
	pr := s.fl.props(st)
	return match s.spec {
		.h { pr.h }
		.s { pr.s }
	}
}

// residual is how far a flash result misses the specification, scaled to
// be dimensionless: (H − H_spec)/RT or (S − S_spec)/R.
fn (s StateFlash) residual(r FlashResult) f64 {
	mut x := 0.0
	for ph in r.phases {
		x += ph.fraction * s.of(ph.state)
	}
	return (x - s.value) / s.scale(r.t)
}

fn (s StateFlash) scale(t f64) f64 {
	return match s.spec {
		.h { thermo.r_gas * t }
		.s { thermo.r_gas }
	}
}

fn (s StateFlash) flash_at(t f64) !FlashResult {
	if !(t >= state_t_min && t <= state_t_max) {
		return fail(.bad_input, "${s.what()}: no temperature from ${state_t_min} K to ${state_t_max} K reaches it, so the specification is outside the model's range")
	}
	return pt(s.fl.model, t, s.p, s.z, tol: s.params.tol, verify: s.params.verify)
}

// ideal_t is the temperature at which the feed as an ideal gas has the
// specified H or S: where the search starts.
fn (s StateFlash) ideal_t() f64 {
	if s.params.estimate > 0.0 {
		return s.params.estimate
	}
	g := fn [s] (t f64) f64 {
		mut x := 0.0
		for i in s.present {
			ig := s.fl.ideal[i]
			x += s.z[i] * match s.spec {
				.h {
					(ig.dhf or { 0.0 }) + ig.cp.int_cp(thermo.t_ref, t)
				}
				.s {
					(ig.s0 or { 0.0 }) + ig.cp.int_cp_over_t(thermo.t_ref, t) - thermo.r_gas * math.log(s.p / thermo.p_ref) - thermo.r_gas * math.log(s.z[i])
				}
			}
		}
		return x - s.value
	}
	return numeric.bisect(g, state_t_min, state_t_max) or { 300.0 }
}

// bracket steps T out from t0 by 50% at a time until the residual changes
// sign, and returns the two temperatures either side.
fn (s StateFlash) bracket(t0 f64) !(f64, f64, int) {
	mut t := t0
	f0 := s.residual(s.flash_at(t)!)
	if f0 == 0.0 {
		return t, t, 1
	}
	up := f0 < 0.0 // H or S too low: raise T
	for k in 0 .. 60 {
		// 1.5 and 2 take the fewest steps on the reference set; 1.25 takes 3%
		// more.
		next := if up { t * 1.5 } else { t / 1.5 }
		fn_ := s.residual(s.flash_at(next)!)
		if fn_ == 0.0 || (fn_ > 0.0) != (f0 > 0.0) {
			return if up { t, next, k + 2 } else { next, t, k + 2 }
		}
		t = next
	}
	return fail(.no_convergence, '${s.what()}: no change of sign in 60 steps from ${t0} K')
}

fn (s StateFlash) mixture() !FlashResult {
	lo, hi, mut steps := s.bracket(s.ideal_t())!
	f := fn [s] (t f64) !f64 {
		return s.residual(s.flash_at(t)!)
	}
	coarse := numeric.brent(f, lo, hi, xtol: 1e-3 * lo, max_iter: s.params.max_steps)!
	steps += coarse.evaluations
	mut r := s.flash_at(coarse.x)!
	steps++
	if r.phases.len == 2 {
		if t, n := s.newton(r) {
			steps += n
			if q := s.flash_at(t) {
				r = q
				steps++
			}
		}
	}
	if math.abs(s.residual(r)) > s.params.etol {
		fine := numeric.brent(f, coarse.a, coarse.b,
			ftol:     s.params.etol
			max_iter: s.params.max_steps
		)!
		steps += fine.evaluations
		r = s.flash_at(fine.x)!
		steps++
	}
	res := s.residual(r)
	if math.abs(res) > s.params.etol {
		return fail(.no_convergence, '${s.what()}: stopped at ${r.t} K with the energy balance off by ${res} (etol ${s.params.etol})')
	}
	return FlashResult{
		...r
		iterations: steps
	}
}

// newton solves Q's stationary conditions in (v, ln T) together, from the
// PT split r: gᵢ = ln fᵢ(y) − ln fᵢ(x) = 0 and the energy balance
// r_E = 0, with r_E = (H − H_spec)/RT or (S − S_spec)/R. The Jacobian is
// Q's Hessian: in v the PT flash's (gibbs), and, with cᵢ = T·∂(ln φᵢⱽ −
// ln φᵢᴸ)/∂T,
//
//     ∂gᵢ/∂ln T = cᵢ,
//     ∂r_E/∂vᵢ = −cᵢ (PH) or −gᵢ − cᵢ (PS),
//     ∂r_E/∂ln T = Cp/R − r_E (PH) or Cp/R (PS),
//
// with Cp the two phases' heat capacity at fixed composition. It returns
// the temperature and the steps taken.
fn (s StateFlash) newton(r FlashResult) !(f64, int) {
	y := r.phases[0].state.n
	mut v0 := s.present.map(r.beta * y[it])
	v0 << math.log(r.t)
	system := fn [s] (x []f64, mut f []f64, mut jac []f64) ! {
		s.eval(x, mut f, mut jac)!
	}
	sol := numeric.newton(system, v0, ftol: s.params.tol, max_iter: s.params.max_steps)!
	return math.exp(sol.x[s.present.len]), sol.iterations
}

// eval is newton's system at (v, ln T).
fn (s StateFlash) eval(v []f64, mut f []f64, mut jac []f64) ! {
	nc := s.z.len
	k := s.present.len
	n := k + 1
	t := math.exp(v[k])
	if !(t >= state_t_min && t <= state_t_max) {
		return error('T left the range searched')
	}
	mut vt, mut lt := 0.0, 0.0
	for j, i in s.present {
		if !(v[j] > 0.0) || !(v[j] < s.z[i]) {
			return error('a mole number left (0, zᵢ)')
		}
		vt += v[j]
		lt += s.z[i] - v[j]
	}
	mut x := []f64{len: nc}
	mut y := []f64{len: nc}
	for j, i in s.present {
		y[i] = v[j] / vt
		x[i] = (s.z[i] - v[j]) / lt
	}
	m := s.fl.model
	sy := thermo.solve_tp(m, t, s.p, y, .stable)!
	sx := thermo.solve_tp(m, t, s.p, x, .stable)!
	dy := thermo.ln_phi_derivs(m, sy)
	dx := thermo.ln_phi_derivs(m, sx)
	py := s.fl.props(sy)
	px := s.fl.props(sx)
	for j, i in s.present {
		f[j] = math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i]
		for l, q in s.present {
			mut hv := dy.dn[i * nc + q] - 1.0
			mut hl := dx.dn[i * nc + q] - 1.0
			if i == q {
				hv += 1.0 / y[i]
				hl += 1.0 / x[i]
			}
			jac[j * n + l] = hv / vt + hl / lt
		}
		c := t * (dy.dt[i] - dx.dt[i])
		jac[j * n + k] = c
		jac[k * n + j] = match s.spec {
			.h { -c }
			.s { -f[j] - c }
		}
	}
	cp := vt * py.cp + lt * px.cp
	match s.spec {
		.h {
			f[k] = (vt * py.h + lt * px.h - s.value) / (thermo.r_gas * t)
			jac[k * n + k] = cp / thermo.r_gas - f[k]
		}
		.s {
			f[k] = (vt * py.s + lt * px.s - s.value) / thermo.r_gas
			jac[k * n + k] = cp / thermo.r_gas
		}
	}
}

// pure flashes a single component. Below its critical pressure it first
// finds the saturation temperature (pvf): an H or S between the saturated
// liquid's and vapour's is the two-phase state there, with β from the
// balance, and one outside is a single phase on its own side of Tsat.
fn (s StateFlash) pure() !FlashResult {
	m := s.fl.model
	if sat := pvf(m, s.p, 0.5, s.z, tol: s.params.tol) {
		vap, liq := sat.phases[0].state, sat.phases[1].state
		xv, xl := s.of(vap), s.of(liq)
		// Within etol of a saturated phase's H or S is that phase: a single
		// phase found there would sit at Tsat, where its two roots tie.
		tol := s.params.etol * s.scale(sat.t)
		if s.value >= xl - tol && s.value <= xv + tol {
			beta := math.max(0.0, math.min(1.0, (s.value - xl) / (xv - xl)))
			return FlashResult{
				...sat
				beta:       beta
				phases:     [Phase{.vapor, vap, beta}, Phase{.liquid, liq, 1.0 - beta}]
				iterations: sat.iterations
			}
		}
		// A liquid below Tsat, or a vapour above it.
		if s.value < xl {
			return s.single(.liquid, sat.t)!
		}
		return s.single(.vapor, sat.t)!
	}
	// No saturation: past the critical pressure, one phase at every T.
	return s.single(.stable, s.ideal_t())!
}

// single is a pure component's single phase on one root, by Brent's method
// in T. A single phase's H and S rise with T, so the search for a change of
// sign goes up from t0 where they're too low, and down where too high.
fn (s StateFlash) single(root thermo.Root, t0 f64) !FlashResult {
	m := s.fl.model
	f := fn [s, m, root] (t f64) !f64 {
		if !(t >= state_t_min && t <= state_t_max) {
			return fail(.bad_input, "${s.what()}: no temperature from ${state_t_min} K to ${state_t_max} K reaches it, so the specification is outside the model's range")
		}
		st := thermo.solve_tp(m, t, s.p, s.z, root)!
		return (s.of(st) - s.value) / s.scale(t)
	}
	mut a, mut b := t0, t0
	fa := f(a)!
	mut found := fa == 0.0
	for _ in 0 .. 60 {
		if found {
			break
		}
		b = if fa < 0.0 { b * 1.5 } else { b / 1.5 }
		fb := f(b)!
		found = fb == 0.0 || (fb > 0.0) != (fa > 0.0)
	}
	if !found {
		return fail(.no_convergence, '${s.what()}: no change of sign in 60 steps from ${t0} K')
	}
	lo, hi := math.min(a, b), math.max(a, b)
	sol := numeric.brent(f, lo, hi, ftol: s.params.etol, max_iter: s.params.max_steps)!
	// Brent also stops when the bracket can't shrink, short of ftol.
	if math.abs(sol.f) > s.params.etol {
		return fail(.no_convergence, '${s.what()}: stopped at ${sol.x} K with the energy balance off by ${sol.f} (etol ${s.params.etol})')
	}
	// The root is the stable one by construction: the liquid below Tsat,
	// the vapour above it, and with no saturation the lower-Gibbs root.
	st := thermo.solve_tp(m, sol.x, s.p, s.z, root)!
	label := thermo.phase_label(m, st)
	return FlashResult{
		t:          sol.x
		p:          s.p
		z:          s.z
		beta:       if label == .vapor { 1.0 } else { 0.0 }
		phases:     [Phase{label, st, 1.0}]
		path:       .single_phase
		iterations: sol.evaluations
		stability:  stability(m, sol.x, s.p, s.z)!
	}
}
