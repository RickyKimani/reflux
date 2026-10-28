module flash

import math
import numeric
import thermo

// Phase is one phase of a flash result.
pub struct Phase {
pub:
	label thermo.PhaseLabel
	// state is the phase at its own composition, one mole in all.
	state thermo.State
	// fraction is the share of the feed's moles in this phase.
	fraction f64
}

// Path says how a flash reached its answer.
pub enum Path {
	single_phase // the feed was stable
	substitution // successive substitution converged
	newton       // Newton on the Gibbs energy finished the split
}

// FlashResult is a converged flash (D7).
pub struct FlashResult {
pub:
	t f64
	p f64
	z []f64 // the feed, normalised
	// beta is the vapour fraction: in a split, the fraction of the less
	// dense phase (D6); for a single phase, 1 if it's labelled vapour and 0
	// if liquid.
	beta f64
	// phases are in order of increasing mass density, so a split's vapour
	// comes first. M3 finds one or two; the type allows more, for M6.
	phases []Phase
	path   Path
	// iterations counts the split's successive-substitution and Newton
	// steps, not the stability test's.
	iterations int
	stability  Stability // the feed's stability test
}

// FlashParams are the PT flash's options.
@[params]
pub struct FlashParams {
pub:
	// A split is converged when every |ln fᵢⱽ − ln fᵢᴸ| ≤ tol.
	tol f64 = 1e-12
	// verify re-tests each phase of a split for stability, with the other
	// phase as a trial (D3). An unstable phase means a third phase or a
	// better split, which is a FlashError (.more_phases) since M3 handles
	// two phases. It roughly doubles the cost of a split.
	verify    bool = true
	max_steps int  = 200
}

// Successive substitution steps before Newton on G takes over.
const split_ss_before_newton = 20
// A split whose every |ln Kᵢ| is below this has collapsed onto the feed.
const split_trivial_tol = 1e-4

// pt flashes a feed z at (T, P), for any thermo.Model (D4).
//
// It first tests the feed's stability. A stable feed is one phase, labelled
// by its phase identification parameter. An unstable one is split into two
// phases by successive substitution on the K-values,
//
//     ln Kᵢ ← ln φᵢᴸ(x) − ln φᵢⱽ(y),
//
// with x, y and β from Rachford–Rice at each step and Michelsen's
// dominant-eigenvalue extrapolation every fifth step. Rachford–Rice may
// return β outside [0, 1] along the way (the negative flash), which keeps
// the iteration going where a clipped β would stall. A split still going
// after 20 steps, as near a critical point, switches to Newton's method on
// the Gibbs energy in the vapour mole numbers, minimising G itself.
//
// The first K-values come from the stability test's incipient phase, and
// if that split fails, from its other unstable trials in turn.
// Each phase takes its volume root of lower Gibbs energy, and the phase of
// lower mass density is the vapour (D6).
//
// A split is accepted only with every |ln fᵢⱽ − ln fᵢᴸ| ≤ tol and
// 0 < β < 1, and, with params.verify, only if neither phase is unstable.
// Anything else is a FlashError.
pub fn pt(m thermo.Model, t f64, p f64, z []f64, params FlashParams) !FlashResult {
	if !(params.tol > 0.0) || params.max_steps < 1 {
		return fail(.bad_input, 'tol and max_steps must be positive')
	}
	stab := stability(m, t, p, z)!
	feed := normalised(z)
	if stab.stable {
		st := thermo.solve_tp(m, t, p, feed, .stable) or {
			return fail(.no_volume_root, 'the feed at ${t} K and ${p} Pa: ${err}')
		}
		label := thermo.phase_label(m, st)
		return FlashResult{
			t:         t
			p:         p
			z:         feed
			beta:      if label == .vapor { 1.0 } else { 0.0 }
			phases:    [Phase{label, st, 1.0}]
			path:      .single_phase
			stability: stab
		}
	}
	sp := Splitter{
		m:       m
		t:       t
		p:       p
		z:       feed
		present: []int{len: feed.len, init: index}.filter(feed[it] > 0.0)
		params:  params
	}
	mut last := IError(fail(.no_convergence, 'no starting K-values'))
	for k0 in sp.starts(stab) {
		split := sp.split(k0) or {
			last = err
			continue
		}
		return sp.result(split, stab)!
	}
	return last
}

fn normalised(z []f64) []f64 {
	total := numeric.sum(z)
	return z.map(it / total)
}

// Splitter is one two-phase split at (T, P).
struct Splitter {
	m       thermo.Model
	t       f64
	p       f64
	z       []f64
	present []int
	params  FlashParams
}

// Split is a converged two-phase split, before labelling.
struct Split {
	beta       f64 // the fraction in phase y
	liquid     f64 // 1 − β, kept separately for a trace
	x          []f64
	y          []f64
	sx         thermo.State
	sy         thermo.State
	path       Path
	iterations int
}

// starts lists the K-values to start from: one from each unstable trial,
// lowest tm first. From an unstable trial phase W (tm = −ln ΣW) the
// stationarity condition Wᵢφᵢ(w) = zᵢφᵢ(z) gives Kᵢ = Wᵢ/zᵢ, taking W as
// phase y. It doesn't matter whether W is the lighter phase: K and 1/K
// describe the same split with the roles swapped (Rachford–Rice gives
// β′ = 1 − β, x′ = y, y′ = x), every step after mirrors the other, and
// result labels the phases by density at the end.
fn (sp Splitter) starts(stab Stability) [][]f64 {
	mut trials := stab.trials.filter(!it.trivial && it.tm < 0.0)
	trials.sort(a.tm < b.tm)
	mut out := [][]f64{}
	for tr in trials {
		scale := math.exp(-tr.tm) // ΣW
		out << []f64{len: sp.z.len, init: if sp.z[index] == 0.0 {
			1.0
		} else {
			scale * tr.w[index] / sp.z[index]
		}}
	}
	return out
}

// state is a phase of composition x at the root of lower Gibbs energy.
fn (sp Splitter) state(x []f64) !thermo.State {
	return thermo.solve_tp(sp.m, sp.t, sp.p, x, .stable) or {
		return fail(.no_volume_root, 'a phase at ${sp.t} K and ${sp.p} Pa: ${err}')
	}
}

// residual is max |ln fᵢ(y) − ln fᵢ(x)|.
fn (sp Splitter) residual(x []f64, sx thermo.State, y []f64, sy thermo.State) f64 {
	mut worst := 0.0
	for i in sp.present {
		worst = math.max(worst, math.abs(math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i]))
	}
	return worst
}

// split runs successive substitution with DEM from K-values k0, and Newton
// on G if substitution is slow.
fn (sp Splitter) split(k0 []f64) !Split {
	nc := sp.z.len
	mut lnk := k0.map(math.log(it))
	mut last_step := []f64{}
	mut since_dem := 0
	mut tried_newton := false
	for iter in 1 .. sp.params.max_steps + 1 {
		rr := numeric.rachford_rice(sp.z, lnk.map(math.exp(it))) or {
			return fail(.no_convergence, 'the K-values left the two-phase range: ${err.msg()}')
		}
		sx := sp.state(rr.x)!
		sy := sp.state(rr.y)!
		if sp.residual(rr.x, sx, rr.y, sy) <= sp.params.tol {
			if !(rr.beta > 0.0 && rr.liquid > 0.0) {
				return fail(.no_convergence, 'the split converged to β = ${rr.beta}, outside (0, 1)')
			}
			return Split{
				beta:       rr.beta
				liquid:     rr.liquid
				x:          rr.x
				y:          rr.y
				sx:         sx
				sy:         sy
				path:       .substitution
				iterations: iter - 1
			}
		}
		mut collapsed := true
		for i in sp.present {
			if math.abs(lnk[i]) >= split_trivial_tol {
				collapsed = false
			}
		}
		if collapsed {
			return fail(.no_convergence, 'the split collapsed onto the feed (every |ln K| < ${split_trivial_tol})')
		}
		if iter > split_ss_before_newton && !tried_newton && rr.beta > 0.0 && rr.liquid > 0.0 {
			tried_newton = true
			v0 := sp.present.map(rr.beta * rr.y[it])
			if done := sp.newton(v0) {
				// minimize can also stop where G stops showing progress, so
				// the split is held to the same test as substitution's.
				if sp.residual(done.x, done.sx, done.y, done.sy) <= sp.params.tol {
					return Split{
						...done
						iterations: iter - 1 + done.iterations
					}
				}
			}
			// Newton failed from here: carry on substituting.
		}
		mut next := []f64{len: nc}
		for i in sp.present {
			next[i] = sx.ln_phi[i] - sy.ln_phi[i]
		}
		step := []f64{len: nc, init: next[index] - lnk[index]}
		since_dem++
		if since_dem >= 5 && last_step.len > 0 {
			since_dem = 0
			// Michelsen's dominant-eigenvalue extrapolation, as in stability.
			lambda := numeric.dot(step, step) / numeric.dot(last_step, step)
			if lambda > 0.0 && lambda < 1.0 {
				jump := lambda / (1.0 - lambda)
				acc := []f64{len: nc, init: next[index] + jump * step[index]}
				if _ := numeric.rachford_rice(sp.z, acc.map(math.exp(it))) {
					next = acc.clone()
				}
			}
		}
		last_step = step.clone()
		lnk = next.clone()
	}
	return fail(.no_convergence, 'the split did not converge in ${sp.params.max_steps} steps')
}

// newton minimises the Gibbs energy over the vapour mole numbers v (one mole
// of feed; the liquid has l = z − v), from v0 over the components present.
// With V = Σv, y = v/V, L = 1 − V and x = l/L (Michelsen & Mollerup,
// chapter 10):
//
//     G/RT = Σᵢ vᵢ(ln yᵢ + ln φᵢ(y)) + Σᵢ lᵢ(ln xᵢ + ln φᵢ(x)),
//     gᵢ = ln fᵢ(y) − ln fᵢ(x),
//     Hᵢⱼ = (δᵢⱼ/yᵢ − 1 + Φⱽᵢⱼ)/V + (δᵢⱼ/xᵢ − 1 + Φᴸᵢⱼ)/L,
//
// where Φ is each phase's n·∂ln φᵢ/∂nⱼ. At a two-phase minimum H is
// positive definite.
fn (sp Splitter) newton(v0 []f64) !Split {
	objective := fn [sp] (v []f64, mut g []f64, mut h []f64) !f64 {
		return sp.gibbs(v, mut g, mut h)
	}
	// G is a sum of terms vᵢ ln fᵢ that don't vanish at the minimum, and its
	// rounding is a few ε times their size.
	mut g := []f64{len: v0.len}
	mut h := []f64{len: v0.len * v0.len}
	sp.gibbs(v0, mut g, mut h)!
	scale := sp.gibbs_scale(v0)!
	sol := numeric.minimize(objective, v0,
		gtol:     sp.params.tol
		fscale:   scale
		max_iter: 50
	)!
	nc := sp.z.len
	mut x := []f64{len: nc}
	mut y := []f64{len: nc}
	mut vt, mut lt := 0.0, 0.0
	for j, i in sp.present {
		vt += sol.x[j]
		lt += sp.z[i] - sol.x[j]
	}
	for j, i in sp.present {
		y[i] = sol.x[j] / vt
		x[i] = (sp.z[i] - sol.x[j]) / lt
	}
	sx := sp.state(x)!
	sy := sp.state(y)!
	return Split{
		beta:       vt
		liquid:     lt
		x:          x
		y:          y
		sx:         sx
		sy:         sy
		path:       .newton
		iterations: sol.iterations
	}
}

// gibbs is G/RT of the split with vapour mole numbers v, with its gradient
// and Hessian, for newton.
fn (sp Splitter) gibbs(v []f64, mut g []f64, mut h []f64) !f64 {
	nc := sp.z.len
	k := sp.present.len
	mut vt, mut lt := 0.0, 0.0
	for j, i in sp.present {
		if !(v[j] > 0.0) || !(v[j] < sp.z[i]) {
			return error('a mole number left (0, zᵢ)')
		}
		vt += v[j]
		lt += sp.z[i] - v[j]
	}
	mut x := []f64{len: nc}
	mut y := []f64{len: nc}
	for j, i in sp.present {
		y[i] = v[j] / vt
		x[i] = (sp.z[i] - v[j]) / lt
	}
	sy := sp.state(y)!
	sx := sp.state(x)!
	dy := thermo.ln_phi_derivs(sp.m, sy).dn
	dx := thermo.ln_phi_derivs(sp.m, sx).dn
	mut gibbs := 0.0
	for j, i in sp.present {
		fy := math.log(y[i]) + sy.ln_phi[i]
		fx := math.log(x[i]) + sx.ln_phi[i]
		gibbs += v[j] * fy + (sp.z[i] - v[j]) * fx
		g[j] = fy - fx
		for l, q in sp.present {
			mut hv := dy[i * nc + q] - 1.0
			mut hl := dx[i * nc + q] - 1.0
			if i == q {
				hv += 1.0 / y[i]
				hl += 1.0 / x[i]
			}
			h[j * k + l] = hv / vt + hl / lt
		}
	}
	return gibbs
}

// gibbs_scale is the size of the terms G is summed from, Σ|vᵢ ln fᵢⱽ| +
// |lᵢ ln fᵢᴸ|, for minimize's fscale.
fn (sp Splitter) gibbs_scale(v []f64) !f64 {
	nc := sp.z.len
	mut vt, mut lt := 0.0, 0.0
	for j, i in sp.present {
		vt += v[j]
		lt += sp.z[i] - v[j]
	}
	mut x := []f64{len: nc}
	mut y := []f64{len: nc}
	for j, i in sp.present {
		y[i] = v[j] / vt
		x[i] = (sp.z[i] - v[j]) / lt
	}
	sy := sp.state(y)!
	sx := sp.state(x)!
	mut s := 0.0
	for j, i in sp.present {
		s += math.abs(v[j] * (math.log(y[i]) + sy.ln_phi[i])) +
			math.abs((sp.z[i] - v[j]) * (math.log(x[i]) + sx.ln_phi[i]))
	}
	return s
}

// mass_density is a phase's mass per volume, kg/m³. It's what D6's "less
// dense" means. Molar volume isn't a substitute: in a mixture of very
// different sizes, the methane-rich phase of methane and n-decane at 300 K
// and 200 bar has the smaller molar volume, yet weighs 168 kg/m³ against
// the decane-rich phase's 585.
fn mass_density(m thermo.Model, st thermo.State) f64 {
	comps := m.components()
	mut mass := 0.0
	for i, n in st.n {
		mass += n * comps[i].mw
	}
	return mass / st.v
}

// verify_tol is the tm a converged phase may show against its own split
// before it counts as unstable: 1e-8, or ten times a looser flash
// tolerance, since an answer good only to tol has tm of about that size.
fn verify_tol(tol f64) f64 {
	return math.max(1e-8, 10.0 * tol)
}

// result labels a split's phases, the less dense by mass being the vapour,
// and with params.verify re-tests each for stability.
fn (sp Splitter) result(s Split, stab Stability) !FlashResult {
	if sp.params.verify {
		for pair in [[s.x, s.y], [s.y, s.x]] {
			check := stability(sp.m, sp.t, sp.p, pair[0],
				trials: [pair[1]]
				tol:    verify_tol(sp.params.tol)
			)!
			if !check.stable {
				tm := if inc := check.incipient() { inc.tm } else { 0.0 }
				return fail(.more_phases, 'a phase of the split at ${sp.t} K and ${sp.p} Pa is unstable (tm = ${tm}): three or more phases, or a better split')
			}
		}
	}
	// y is the vapour when it's the less dense; otherwise swap.
	y_lighter := mass_density(sp.m, s.sy) < mass_density(sp.m, s.sx)
	vap, liq := if y_lighter { s.sy, s.sx } else { s.sx, s.sy }
	beta, liquid := if y_lighter { s.beta, s.liquid } else { s.liquid, s.beta }
	return FlashResult{
		t:          sp.t
		p:          sp.p
		z:          sp.z
		beta:       beta
		phases:     [Phase{.vapor, vap, beta}, Phase{.liquid, liq, liquid}]
		path:       s.path
		iterations: s.iterations
		stability:  stab
	}
}
