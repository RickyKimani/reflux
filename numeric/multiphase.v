module numeric

import math

// MultiphaseSplit is the phase fractions multiphase_rr finds for fixed
// K-values, and each phase's mole fractions.
pub struct MultiphaseSplit {
pub:
	// beta is each phase's fraction of the feed, β ≥ 0. The phases present
	// (β > 0) sum to 1.
	beta []f64
	// x[j] is phase j's mole fractions: zᵢKᵢⱼ/Eᵢ, normalised.
	x [][]f64
	// sum_x[j] is Σᵢ zᵢKᵢⱼ/Eᵢ before normalising, 1 − ∂Q/∂βⱼ: 1 for a phase
	// present, at most 1 for one held at β = 0, and above 1 only for a phase
	// that should appear but whose K-values put it nowhere.
	sum_x      []f64
	iterations int
}

// The most Newton steps multiphase_rr takes. Q is convex and its Newton
// step exact on a quadratic, so it needs a handful; the bound on β can cost
// one step for each phase that leaves or joins.
const multiphase_max_steps = 100

// multiphase_rr is Michelsen's (1994) multiphase Rachford–Rice, as a
// minimisation: for phases j with Kᵢⱼ = 1/φᵢⱼ (given as ln Kᵢⱼ, so that
// their ratios are a split's K-values),
//
//     Q(β) = Σⱼ βⱼ − Σᵢ zᵢ ln Eᵢ,   Eᵢ = Σⱼ βⱼKᵢⱼ,   β ≥ 0,
//
// and xᵢⱼ = zᵢKᵢⱼ/Eᵢ. Its gradient is ∂Q/∂βⱼ = 1 − Σᵢ xᵢⱼ and its Hessian
// Σᵢ zᵢKᵢⱼKᵢₖ/Eᵢ², positive semi-definite, so Q is convex and its minimum
// over β ≥ 0 unique where the Kⱼ are independent. There, each phase
// present has Σᵢ xᵢⱼ = 1, so the phases' fractions sum to 1 and
// Σⱼ βⱼxᵢⱼ = zᵢ; a phase held at β = 0 has Σᵢ xᵢⱼ ≤ 1. It's what a
// two-phase Rachford–Rice is with the negative flash cut off at 0 and 1, and
// the phases' roles symmetric, so a phase leaves by its fraction reaching 0
// and nothing else.
//
// Each component's K-values are scaled by its largest, so that none
// overflows; Q moves by a constant, and x not at all. Newton's method runs
// on the free phases (β > 0, or β = 0 with ∂Q/∂β < 0), the step cut short
// to keep β ≥ 0 and halved until Q falls, and stops when every free
// gradient is within tol (a phase held at 0 has ∂Q/∂β ≥ 0 already). z is the
// components present, each above 0; beta0 is the start, β ≥ 0, not all 0.
pub fn multiphase_rr(z []f64, ln_k [][]f64, beta0 []f64, tol f64) !MultiphaseSplit {
	np := ln_k.len
	nc := z.len
	if np < 1 || beta0.len != np || !(tol > 0.0) {
		return error('multiphase_rr: ${np} phases, ${beta0.len} starting fractions, tol ${tol}')
	}
	for i, zi in z {
		if !(zi > 0.0) || !math.is_finite(zi) {
			return error('multiphase_rr: z[${i}] = ${zi}')
		}
	}
	mut k := [][]f64{len: np, init: []f64{len: nc}}
	for i in 0 .. nc {
		mut top := math.inf(-1)
		for j in 0 .. np {
			if ln_k[j].len != nc || !math.is_finite(ln_k[j][i]) {
				return error('multiphase_rr: ln K of component ${i} in phase ${j} is not finite')
			}
			top = math.max(top, ln_k[j][i])
		}
		for j in 0 .. np {
			k[j][i] = math.exp(ln_k[j][i] - top)
		}
	}
	mut beta := beta0.clone()
	if beta.any(!(it >= 0.0) || !math.is_finite(it)) || !(sum(beta) > 0.0) {
		return error('multiphase_rr: a start of ${beta0}')
	}
	mut e := []f64{len: nc}
	mut g := []f64{len: np}
	for iter in 0 .. multiphase_max_steps {
		q := multiphase_q(z, k, beta, mut e) or { return error('multiphase_rr: ${err.msg()}') }
		for j in 0 .. np {
			mut s := 0.0
			for i in 0 .. nc {
				s += z[i] * k[j][i] / e[i]
			}
			g[j] = 1.0 - s
		}
		mut loose := []int{len: np, init: index}.filter(beta[it] > 0.0 || g[it] < 0.0)
		// A phase outside the free set has β = 0 and ∂Q/∂β ≥ 0 already, so the
		// free ones' gradients are the whole test.
		mut done := true
		for j in loose {
			if math.abs(g[j]) > tol {
				done = false
			}
		}
		if done {
			return multiphase_split(z, k, beta, e, iter)
		}
		free_set := loose.clone()
		// Newton on the free phases, H_FF d = −g_F. A phase at β = 0 that the
		// step would take below 0 is held there and the step taken again
		// without it, so that every phase moves or stays inside β ≥ 0.
		mut moving := loose.clone()
		mut d := []f64{}
		for {
			d = multiphase_step(z, k, e, g, moving) or {
				return error('multiphase_rr: ${err.msg()} at β = ${beta}')
			}
			held := []int{len: moving.len, init: index}.filter(beta[moving[it]] == 0.0 && d[it] < 0.0)
			if held.len == 0 {
				break
			}
			moving = []int{len: moving.len, init: index}.filter(it !in held).map(moving[it])
			if moving.len == 0 {
				return error('multiphase_rr: no phase can move at β = ${beta}')
			}
		}
		loose = moving.clone()
		// The longest step that keeps β ≥ 0, then halve it until Q falls, or
		// rises by no more than its rounding: near the minimum a Newton step
		// moves Q by about the gradient squared, far below a few ulps of the
		// terms it's summed from, and is the step to take.
		//
		// The phase that sets the longest step reaches its bound, so it goes
		// to 0 exactly: β + (−β/d)·d rounds to a residue of about β·ε, not 0,
		// and a phase leaving by residues shrinks by ε a step until −β/d
		// underflows to 0 and no step moves at all (two components in three
		// phases, whose Hessian is singular, did so in 20 steps).
		mut alpha := 1.0
		mut block := -1
		for a, ja in loose {
			if d[a] < 0.0 && -beta[ja] / d[a] <= alpha {
				alpha = -beta[ja] / d[a]
				block = a
			}
		}
		mut scale := sum(beta)
		for i in 0 .. nc {
			scale += z[i] * math.abs(math.log(e[i]))
		}
		mut moved := false
		for halving in 0 .. 60 {
			mut trial := beta.clone()
			for a, ja in loose {
				trial[ja] = if a == block && halving == 0 {
					0.0
				} else {
					math.max(0.0, beta[ja] + alpha * d[a])
				}
			}
			if qt := multiphase_q(z, k, trial, mut e) {
				if trial != beta && qt <= q + 8.0 * f64_epsilon * scale {
					beta = trial.clone()
					moved = true
					break
				}
			}
			alpha *= 0.5
		}
		if !moved {
			worst := []f64{len: np, init: if index in free_set {
				math.abs(g[index])
			} else {
				-g[index]
			}}
			return error('multiphase_rr: Q falls no further at β = ${beta}, with a gradient of ${worst} against tol = ${tol}')
		}
	}
	return error('multiphase_rr: no minimum in ${multiphase_max_steps} steps from β = ${beta0}')
}

// multiphase_q is Q(β), less a constant, with each Eᵢ written into e.
fn multiphase_q(z []f64, k [][]f64, beta []f64, mut e []f64) !f64 {
	mut q := sum(beta)
	for i in 0 .. z.len {
		mut ei := 0.0
		for j in 0 .. beta.len {
			ei += beta[j] * k[j][i]
		}
		if !(ei > 0.0) {
			return error('component ${i} is in no phase at β = ${beta}')
		}
		e[i] = ei
		q -= z[i] * math.log(ei)
	}
	return q
}

fn multiphase_split(z []f64, k [][]f64, beta []f64, e []f64, iter int) MultiphaseSplit {
	np := beta.len
	mut x := [][]f64{len: np, init: []f64{len: z.len}}
	mut sums := []f64{len: np}
	for j in 0 .. np {
		mut s := 0.0
		for i in 0 .. z.len {
			x[j][i] = z[i] * k[j][i] / e[i]
			s += x[j][i]
		}
		sums[j] = s
		for i in 0 .. z.len {
			x[j][i] /= s
		}
	}
	// The phases present sum to 1 at the minimum; dividing by the sum takes
	// out what's left of the gradient's tolerance.
	total := sum(beta)
	return MultiphaseSplit{
		beta:       beta.map(it / total)
		x:          x
		sum_x:      sums
		iterations: iter
	}
}

// multiphase_step is Newton's step on the phases in moving: H d = −g there,
// H positive definite unless two phases have the same K-values, when a
// small shift on its diagonal still descends.
fn multiphase_step(z []f64, k [][]f64, e []f64, g []f64, moving []int) ![]f64 {
	nf := moving.len
	mut h := []f64{len: nf * nf}
	mut rhs := []f64{len: nf}
	for a, ja in moving {
		rhs[a] = -g[ja]
		for b, jb in moving {
			mut s := 0.0
			for i in 0 .. z.len {
				s += z[i] * k[ja][i] * k[jb][i] / (e[i] * e[i])
			}
			h[a * nf + b] = s
		}
	}
	for shift in [0.0, 1e-12, 1e-8, 1e-4] {
		mut hs := h.clone()
		for a in 0 .. nf {
			hs[a * nf + a] += shift * (1.0 + hs[a * nf + a])
		}
		c := cholesky(hs, nf) or { continue }
		return c.solve(rhs) or { continue }
	}
	return error("the Hessian can't be factored")
}
