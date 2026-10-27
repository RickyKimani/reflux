module flash

import math
import numeric
import thermo

// TrialStart says where a stability trial began.
pub enum TrialStart {
	vapor_wilson  // Wᵢ = zᵢKᵢ with Wilson's K: a vapour-like trial phase
	liquid_wilson // Wᵢ = zᵢ/Kᵢ: a liquid-like one
	given         // a composition the caller passed in, such as a phase of a split
}

// Trial is one search for a stationary point of the tangent-plane distance.
pub struct Trial {
pub:
	start TrialStart
	// w is where the search ended, as mole fractions.
	w []f64
	// tm is the tangent-plane distance at w: the Gibbs energy, divided by RT,
	// of moving a little of the feed into a new phase of composition w.
	// Negative means the feed is unstable.
	tm f64
	// trivial is true when the search ended at the feed itself, w = z, where
	// tm = 0 and nothing is learned.
	trivial bool
	// converged is true when w is a stationary point. A trial that didn't
	// converge still proves instability if its tm is negative.
	converged  bool
	newton     bool // whether the second-order finish was needed
	iterations int
}

// Stability is the verdict on a feed at (T, P), with every trial behind it.
pub struct Stability {
pub:
	stable bool
	trials []Trial
}

// incipient returns the non-trivial trial with the lowest tm: for an
// unstable feed, the first estimate of the phase that would appear, and a
// good start for the split. It returns none for a stable feed.
pub fn (s Stability) incipient() ?Trial {
	if s.stable {
		return none
	}
	mut best := -1
	for i, tr in s.trials {
		if !tr.trivial && (best < 0 || tr.tm < s.trials[best].tm) {
			best = i
		}
	}
	return if best < 0 { none } else { s.trials[best] }
}

// StabilityParams are stability's options.
@[params]
pub struct StabilityParams {
pub:
	// trials are extra starting compositions, as mole fractions: the phases
	// of a split being re-checked, for instance (D3).
	trials [][]f64
	// The feed is unstable when some trial reaches tm < −tol.
	tol f64 = 1e-10
	// max_steps caps each trial's successive-substitution steps; the Newton
	// finish has 80 of its own. A trial that runs out still proves the feed
	// unstable if its tm is below −tol. Otherwise stability returns a
	// FlashError rather than call the feed stable.
	max_steps int = 500
}

// Successive substitution steps before the Newton finish takes over.
const ss_before_newton = 20
// A stationary point has every |√Wᵢ·rᵢ| below this: the gradient of tm*
// in α = 2√W.
const stationary_tol = 1e-12
// A trial within this of the feed in every ln(Wᵢ/zᵢ) is the feed itself.
const trivial_tol = 1e-5

// stability tests whether a feed z at (T, P) is stable as a single phase,
// by Michelsen's (1982) tangent-plane analysis. A trial phase of mole numbers
// W, compared with the feed's tangent plane, has
//
//     tm*(W) = 1 + Σᵢ Wᵢ(ln Wᵢ + ln φᵢ(w) − dᵢ − 1),   dᵢ = ln zᵢ + ln φᵢ(z),
//
// where w = W/ΣW. At a stationary point ln Wᵢ + ln φᵢ(w) = dᵢ, and there
// tm* = 1 − ΣW, while the normalised distance is tm = −ln ΣW. The feed is
// unstable if any stationary point has tm < 0. Each phase takes its volume
// root of lower Gibbs energy.
//
// Two trials start from Wilson's K, one vapour-like and one liquid-like,
// plus any in params.trials. Each runs successive substitution,
// ln Wᵢ ← dᵢ − ln φᵢ(w), which lowers tm* at every step, accelerated every
// fifth step by Michelsen's dominant-eigenvalue extrapolation when that
// lowers tm* further. A trial still going after 20 steps, typically near a
// critical point or a phase boundary where substitution crawls, switches to
// Newton's method in αᵢ = 2√Wᵢ, whose Hessian stays well conditioned.
//
// A trial that ends at the feed itself is trivial and is never offered as a
// new phase. The feed is stable only if every trial converged without
// reaching tm < −tol; a trial that runs out of steps without proving
// instability is a FlashError, not a guess.
//
// A pure component is always stable here: it can't split by composition,
// and solve_tp has already taken its root of lower Gibbs energy.
pub fn stability(m thermo.Model, t f64, p f64, z []f64, params StabilityParams) !Stability {
	comps := m.components()
	nc := comps.len
	if z.len != nc {
		return fail(.bad_input, 'z has ${z.len} entries for ${nc} components')
	}
	if !(t > 0.0) || !(p > 0.0) || !math.is_finite(t) || !math.is_finite(p) {
		return fail(.bad_input, 'T = ${t} K and P = ${p} Pa must be positive')
	}
	if params.max_steps < 1 {
		return fail(.bad_input, 'max_steps must be at least 1')
	}
	mut total := 0.0
	for v in z {
		if !(v >= 0.0) || !math.is_finite(v) {
			return fail(.bad_input, 'z = ${z} is not a composition')
		}
		total += v
	}
	if total == 0.0 {
		return fail(.bad_input, 'z is all zeros')
	}
	feed_x := z.map(it / total)
	present := []int{len: nc, init: index}.filter(feed_x[it] > 0.0)
	if present.len < 2 {
		return Stability{
			stable: true
		}
	}
	feed_state := thermo.solve_tp(m, t, p, feed_x, .stable) or {
		return fail(.no_volume_root, 'the feed at ${t} K and ${p} Pa: ${err}')
	}
	mut d := []f64{len: nc}
	for i in present {
		d[i] = math.log(feed_x[i]) + feed_state.ln_phi[i]
	}
	f := Feed{
		m:       m
		t:       t
		p:       p
		z:       feed_x
		d:       d
		present: present
	}
	k := wilson_k(comps, t, p)
	mut starts := [
		Start{.vapor_wilson, []f64{len: nc, init: feed_x[index] * k[index]}},
		Start{.liquid_wilson, []f64{len: nc, init: feed_x[index] / k[index]}},
	]
	for w in params.trials {
		if w.len != nc {
			return fail(.bad_input, 'a given trial has ${w.len} entries for ${nc} components')
		}
		starts << Start{.given, w}
	}
	mut trials := []Trial{}
	mut stable := true
	mut unconverged := 0
	for s in starts {
		tr := f.trial(s, params.max_steps)!
		if !tr.trivial && tr.tm < -params.tol {
			stable = false
		}
		if !tr.converged {
			unconverged++
		}
		trials << tr
	}
	if stable && unconverged > 0 {
		return fail(.no_convergence, 'stability at ${t} K and ${p} Pa: ${unconverged} trial(s) did not converge, and none showed instability')
	}
	return Stability{
		stable: stable
		trials: trials
	}
}

struct Start {
	kind TrialStart
	w    []f64
}

// Feed is one stability test: the feed and its tangent plane.
struct Feed {
	m       thermo.Model
	t       f64
	p       f64
	z       []f64
	d       []f64 // dᵢ = ln zᵢ + ln φᵢ(z), for the components present
	present []int // the components with zᵢ > 0; the others stay out of every trial
}

// Point is a trial W and everything the iterations need there.
struct Point {
	u   []f64 // ln Wᵢ
	w   []f64 // Wᵢ
	sum f64   // ΣW
	r   []f64 // rᵢ = ln Wᵢ + ln φᵢ(w) − dᵢ, zero at a stationary point
	st  thermo.State
	// tm_star is the unnormalised tm*(W), which substitution lowers at every
	// step; tm is the normalised tangent-plane distance of the composition w.
	tm_star f64
	tm      f64
}

// at evaluates the trial with ln W = u.
fn (f Feed) at(u []f64) !Point {
	nc := f.z.len
	mut w := []f64{len: nc}
	for i in f.present {
		// exp(709) is the largest finite f64. A trial that big has long since
		// shown the feed to be unstable.
		w[i] = math.exp(math.min(u[i], 700.0))
	}
	sum := numeric.sum(w)
	if !(sum > 0.0) || !math.is_finite(sum) {
		return fail(.no_convergence, 'a stability trial left the composition space (ΣW = ${sum})')
	}
	x := w.map(it / sum)
	st := thermo.solve_tp(f.m, f.t, f.p, x, .stable) or {
		return fail(.no_volume_root, 'a trial phase at ${f.t} K and ${f.p} Pa: ${err}')
	}
	mut r := []f64{len: nc}
	mut tm_star := 1.0
	mut tm := -math.log(sum)
	for i in f.present {
		r[i] = u[i] + st.ln_phi[i] - f.d[i]
		tm_star += w[i] * (r[i] - 1.0)
		tm += x[i] * r[i]
	}
	return Point{
		u:       u.clone()
		w:       w
		sum:     sum
		r:       r
		st:      st
		tm_star: tm_star
		tm:      tm
	}
}

// stationary reports whether the gradient of tm* in α = 2√W, √Wᵢ·rᵢ, has
// vanished.
fn (f Feed) stationary(pt Point) bool {
	for i in f.present {
		if math.abs(math.sqrt(pt.w[i]) * pt.r[i]) > stationary_tol {
			return false
		}
	}
	return true
}

// trivial reports whether W is the feed itself.
fn (f Feed) trivial(pt Point) bool {
	for i in f.present {
		if math.abs(pt.u[i] - math.log(f.z[i])) > trivial_tol {
			return false
		}
	}
	return true
}

fn (f Feed) result(start TrialStart, pt Point, converged bool, newton bool, iterations int) Trial {
	return Trial{
		start:      start
		w:          pt.w.map(it / pt.sum)
		tm:         pt.tm
		trivial:    f.trivial(pt)
		converged:  converged
		newton:     newton
		iterations: iterations
	}
}

// trial runs one search from s: successive substitution with DEM, then the
// Newton finish if substitution is slow.
fn (f Feed) trial(s Start, max_steps int) !Trial {
	nc := f.z.len
	mut u := []f64{len: nc}
	for i in f.present {
		u[i] = math.log(math.max(s.w[i], 1e-300))
	}
	mut pt := f.at(u)!
	mut best := pt // the lowest tm*, the proof of instability if nothing converges
	mut last_step := []f64{}
	mut since_dem := 0
	mut tried_newton := false
	for iter in 1 .. max_steps + 1 {
		if f.stationary(pt) || f.trivial(pt) {
			return f.result(s.kind, pt, true, tried_newton, iter - 1)
		}
		if iter > ss_before_newton && !tried_newton {
			tried_newton = true
			if done := f.newton(pt) {
				return f.result(s.kind, done.point, true, true, iter - 1 + done.iterations)
			}
			// Newton failed from here: carry on substituting.
		}
		// Successive substitution: ln Wᵢ ← dᵢ − ln φᵢ(w), a step of −rᵢ.
		step := []f64{len: nc, init: -pt.r[index]}
		next_u := []f64{len: nc, init: pt.u[index] + step[index]}
		mut next := ?Point(none)
		since_dem++
		if since_dem >= 5 && last_step.len > 0 {
			since_dem = 0
			// Michelsen's dominant-eigenvalue method. Substitution converges
			// like a geometric series with ratio λ, the dominant eigenvalue
			// of its iteration matrix; summing the rest of the series in one
			// go is a step of λ/(1 − λ) times the last one. Michelsen and
			// Mollerup keep the extrapolation only if it lowers tm*; on 6,000
			// random states that test changed no verdict, so it isn't made.
			lambda := numeric.dot(step, step) / numeric.dot(last_step, step)
			if lambda > 0.0 && lambda < 1.0 {
				jump := lambda / (1.0 - lambda)
				acc_u := []f64{len: nc, init: next_u[index] + jump * step[index]}
				next = f.at(acc_u) or { none }
			}
		}
		last_step = step.clone()
		pt = next or { f.at(next_u)! }
		if pt.tm_star < best.tm_star {
			best = pt
		}
	}
	// Out of steps. A negative tm still proves instability.
	return f.result(s.kind, best, false, tried_newton, max_steps)
}

struct NewtonFinish {
	point      Point
	iterations int
}

// newton finishes a trial by minimising tm* in αᵢ = 2√Wᵢ with Newton's
// method (Michelsen & Mollerup, chapter 10), with gradient and Hessian
//
//     gᵢ = √Wᵢ·rᵢ,   Hᵢⱼ = δᵢⱼ(1 + rᵢ/2) + √(WᵢWⱼ)·(n∂ln φᵢ/∂nⱼ)/ΣW.
//
// The δᵢⱼ term dominates, so H stays well conditioned even with traces,
// which is why α beats ln W here. It minimises tm* itself (numeric.minimize)
// rather than solving g = 0: near a critical point the way to the trivial
// solution is a long, nearly flat valley, along which ‖g‖ has dips that
// aren't stationary points, and a root-finder stops in them.
fn (f Feed) newton(start Point) !NewtonFinish {
	objective := fn [f] (a []f64, mut g []f64, mut h []f64) !f64 {
		return f.objective(a, mut g, mut h)
	}
	a0 := f.present.map(2.0 * math.sqrt(start.w[it]))
	// tm* = 1 + Σ Wᵢ(rᵢ − 1) is a small difference of terms of size about 1.
	sol := numeric.minimize(objective, a0, gtol: stationary_tol, fscale: 2.0, max_iter: 80)!
	return NewtonFinish{
		point:      f.at(f.u_of(sol.x))!
		iterations: sol.iterations
	}
}

// objective is tm* as a function of α, over the components present, with
// its gradient and Hessian, for newton.
fn (f Feed) objective(a []f64, mut g []f64, mut h []f64) !f64 {
	for v in a {
		if !(v > 0.0) {
			return error('α must stay positive')
		}
	}
	pt := f.at(f.u_of(a))!
	der := thermo.ln_phi_derivs(f.m, pt.st)
	k := f.present.len
	nc := f.z.len
	for j, i in f.present {
		g[j] = math.sqrt(pt.w[i]) * pt.r[i]
		for l, q in f.present {
			h[j * k + l] = math.sqrt(pt.w[i] * pt.w[q]) * der.dn[i * nc + q] / pt.sum
		}
		h[j * k + j] += 1.0 + 0.5 * pt.r[i]
	}
	return pt.tm_star
}

// u_of maps α, over the components present, to ln W over all of them:
// Wᵢ = αᵢ²/4.
fn (f Feed) u_of(a []f64) []f64 {
	mut u := []f64{len: f.z.len}
	for j, i in f.present {
		u[i] = math.log(0.25 * a[j] * a[j])
	}
	return u
}
