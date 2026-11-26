module flash

import math
import numeric
import thermo

// The most phases a PT flash gives (M6's D8). A phase found unstable in a
// three-phase answer is a fourth, refused as .more_phases.
pub const max_phases = 3

// Successive substitution steps on a split of three before Newton on G.
const multi_ss_before_newton = 30

// A phase held at no moles this many steps running leaves the split.
const multi_absent_steps = 3

// Multi is a split of the feed into up to max_phases phases at (T, P), the
// stage after a two-phase split one of whose phases is unstable (D8):
// Michelsen's stage-wise flash.
struct Multi {
	m       thermo.PhaseModel
	t       f64
	p       f64
	z       []f64
	present []int
	params  FlashParams
	// newton_after is the substitution steps before Newton on G: a test
	// can bring it forward, to finish from a real iterate.
	newton_after int = multi_ss_before_newton
}

// MultiSplit is a converged split into two or more phases, before labelling.
struct MultiSplit {
	beta       []f64
	x          [][]f64
	st         []thermo.State
	path       Path
	iterations int
}

fn (mu Multi) state(x []f64) !thermo.State {
	return mu.m.phase(mu.t, mu.p, x, .stable) or {
		return fail(.no_volume_root, 'a phase at ${mu.t} K and ${mu.p} Pa: ${err}')
	}
}

// residual is max |ln fᵢⱼ − ln fᵢ₀| over the phases present, against the
// first, and the components above trace_floor in both: equal fugacities
// are equal pairwise, so any phase present serves.
fn (mu Multi) residual(beta []f64, x [][]f64, st []thermo.State) f64 {
	mut worst := 0.0
	for j in 1 .. beta.len {
		if !(beta[j] > 0.0) || !(beta[0] > 0.0) {
			return math.inf(1)
		}
		for i in mu.present {
			if x[j][i] < trace_floor || x[0][i] < trace_floor {
				continue
			}
			worst = math.max(worst, math.abs(math.log(x[j][i]) + st[j].ln_phi[i] - math.log(x[0][i]) - st[0].ln_phi[i]))
		}
	}
	return worst
}

// split runs successive substitution from the phases x0 with fractions
// beta0: at each step, each phase's φ at its composition, then the phases'
// fractions and compositions from Michelsen's multiphase Rachford–Rice on
// Kᵢⱼ = 1/φᵢⱼ (numeric.multiphase_rr), which holds every fraction ≥ 0. A
// phase held at no moles multi_absent_steps steps running leaves, and of two
// phases that come together one does. A split still going after
// multi_ss_before_newton steps goes on by Newton on G. The answer has every
// |ln fᵢⱼ − ln fᵢᵣ| ≤ tol and every fraction above 0.
fn (mu Multi) split(x0 [][]f64, beta0 []f64) !MultiSplit {
	nc := mu.z.len
	zp := mu.present.map(mu.z[it])
	mut x := [][]f64{}
	for row in x0 {
		x << row.clone()
	}
	mut beta := beta0.clone()
	mut absent := []int{len: x.len}
	mut tried_newton := false
	mut balanced := false
	rr_tol := math.max(1e-15, 0.01 * mu.params.tol)
	for iter in 1 .. mu.params.max_steps + 1 {
		mut st := []thermo.State{}
		for xj in x {
			st << mu.state(xj)!
		}
		// Only fractions and compositions from one Rachford–Rice solve hold
		// the balance: after a phase leaves, not before the next.
		if balanced && mu.residual(beta, x, st) <= mu.params.tol && beta.all(it > 0.0) {
			mu.balance(beta, x)!
			return MultiSplit{beta, x, st, .substitution, iter - 1}
		}
		if iter > mu.newton_after && !tried_newton && beta.all(it > 0.0) {
			tried_newton = true
			if done := mu.newton(beta, x) {
				if mu.residual(done.beta, done.x, done.st) <= mu.params.tol {
					mu.balance(done.beta, done.x)!
					return MultiSplit{
						...done
						iterations: iter - 1 + done.iterations
					}
				}
			}
		}
		mut ln_k := [][]f64{}
		for s in st {
			ln_k << mu.present.map(-s.ln_phi[it])
		}
		start := beta.map(if it > 0.0 { it } else { 0.0 })
		rr := numeric.multiphase_rr(zp, ln_k, start, rr_tol) or {
			return fail(.no_convergence, 'the phases fractions: ${err.msg()}')
		}
		beta = rr.beta.clone()
		balanced = true
		for j in 0 .. x.len {
			mut xj := []f64{len: nc}
			for a, i in mu.present {
				xj[i] = rr.x[j][a]
			}
			x[j] = xj.clone()
			absent[j] = if beta[j] > 0.0 { 0 } else { absent[j] + 1 }
		}
		// A phase gone for good, or two phases become one.
		mut leaving := []int{len: x.len, init: index}.filter(absent[it] >= multi_absent_steps)
		for j in 0 .. x.len {
			for l in j + 1 .. x.len {
				if j in leaving || l in leaving {
					continue
				}
				mut same := true
				for i in mu.present {
					if math.abs(math.log(x[j][i] / x[l][i])) >= split_trivial_tol {
						same = false
						break
					}
				}
				if same {
					leaving << if beta[j] < beta[l] { j } else { l }
				}
			}
		}
		if leaving.len > 0 {
			keep := []int{len: x.len, init: index}.filter(it !in leaving)
			if keep.len < 2 {
				return fail(.no_convergence, 'the split collapsed to one phase')
			}
			x = keep.map(x[it].clone())
			beta = keep.map(beta[it])
			balanced = false
			absent = keep.map(absent[it])
			total := numeric.sum(beta)
			beta = beta.map(it / total)
		}
	}
	return fail(.no_convergence, 'the split into ${x.len} phases did not converge in ${mu.params.max_steps} steps')
}

// balance refuses a split whose phases don't add up to the feed, each
// component to 1e-12: what a split must show besides its fugacities.
fn (mu Multi) balance(beta []f64, x [][]f64) ! {
	for i in mu.present {
		mut b := 0.0
		for j in 0 .. beta.len {
			b += beta[j] * x[j][i]
		}
		if !(math.abs(b - mu.z[i]) <= 1e-12) {
			return fail(.no_convergence, 'the split into ${beta.len} phases is off the feed by ${b - mu.z[i]} in component ${i}')
		}
	}
}

// MultiGibbs is a split's G/RT as a numeric.Objective, for newton, in
// variables scaled by the square root of each starting mole number, uₐ =
// nₐ/sₐ, sₐ = √nₐ⁰: a trace's Hessian entry, 1/n, is 1e29 for water's decane
// beside entries of order 1, and Cholesky can't factor that; scaled, each
// diagonal entry is of order 1 (Michelsen and Mollerup's scaling).
struct MultiGibbs {
	mu  Multi
	np  int
	ref []int // each component present's dependent phase
	s   []f64 // each variable's scale
}

fn (o MultiGibbs) eval(u []f64, mut g []f64, mut h []f64) !f64 {
	n := u.len
	v := []f64{len: n, init: u[index] * o.s[index]}
	f := o.mu.gibbs(o.np, o.ref, v, mut g, mut h)!
	for a in 0 .. n {
		g[a] *= o.s[a]
		for b in 0 .. n {
			h[a * n + b] *= o.s[a] * o.s[b]
		}
	}
	return f
}

// slot is where component a's moles in phase j (not its dependent phase,
// ref[a]) are among Newton's variables.
fn slot(np int, ref []int, a int, j int) int {
	return a * (np - 1) + if j < ref[a] { j } else { j - 1 }
}

// amounts are each phase's mole numbers from Newton's variables: each
// component present's moles in every phase but its dependent one, which
// takes the rest of the feed. That phase is the one holding the most of
// it, so the rest is a difference of nearly equal numbers nowhere: water's
// trace of decane, 1e-29, is a variable of its own, where z − Σ would be
// rounding.
fn (mu Multi) amounts(np int, ref []int, v []f64) ![][]f64 {
	mut n := [][]f64{len: np, init: []f64{len: mu.z.len}}
	for a, i in mu.present {
		mut rest := mu.z[i]
		for j in 0 .. np {
			if j == ref[a] {
				continue
			}
			nij := v[slot(np, ref, a, j)]
			if !(nij > 0.0) {
				return error('a mole number left (0, zᵢ)')
			}
			n[j][i] = nij
			rest -= nij
		}
		if !(rest > 0.0) {
			return error('a mole number left (0, zᵢ)')
		}
		n[ref[a]][i] = rest
	}
	return n
}

// gibbs is G/RT = Σⱼ Σᵢ nᵢⱼ ln fᵢⱼ over the phases' mole numbers, with its
// gradient and Hessian in Newton's variables (Michelsen & Mollerup, chapter
// 10, for more than two phases, each component with a dependent phase rᵢ
// of its own):
//
//     gᵢⱼ = ln fᵢⱼ − ln fᵢᵣᵢ,
//     Hᵢⱼ,ₖₘ = (δⱼₘ − δⱼᵣₖ)Aʲᵢₖ − (δᵣᵢₘ − δᵣᵢᵣₖ)Aʳⁱᵢₖ,
//     Aʲᵢₖ = (δᵢₖ/xᵢⱼ − 1 + Φʲᵢₖ)/Nⱼ,
//
// Nⱼ a phase's moles and Φʲ its n·∂ln φᵢ/∂nₖ. With one dependent phase l
// for all, H = δⱼₘAʲ + Aˡ, the two-phase flash's form.
fn (mu Multi) gibbs(np int, ref []int, v []f64, mut g []f64, mut h []f64) !f64 {
	nc := mu.z.len
	k := mu.present.len
	n := mu.amounts(np, ref, v)!
	mut lnf := [][]f64{len: np, init: []f64{len: nc}}
	mut a := [][]f64{len: np, init: []f64{len: k * k}}
	mut gibbs := 0.0
	for j in 0 .. np {
		total := numeric.sum(n[j])
		x := n[j].map(it / total)
		st := mu.state(x)!
		dn := mu.m.ln_phi_derivs(st)!.dn
		for p, i in mu.present {
			lnf[j][i] = math.log(x[i]) + st.ln_phi[i]
			gibbs += n[j][i] * lnf[j][i]
			for q, l in mu.present {
				mut aa := dn[i * nc + l] - 1.0
				if i == l {
					aa += 1.0 / x[i]
				}
				a[j][p * k + q] = aa / total
			}
		}
	}
	size := (np - 1) * k
	for p, i in mu.present {
		for j in 0 .. np {
			if j == ref[p] {
				continue
			}
			row := slot(np, ref, p, j)
			g[row] = lnf[j][i] - lnf[ref[p]][i]
			for q in 0 .. k {
				for m in 0 .. np {
					if m == ref[q] {
						continue
					}
					mut hv := 0.0
					if j == m {
						hv += a[j][p * k + q]
					}
					if j == ref[q] {
						hv -= a[j][p * k + q]
					}
					if ref[p] == m {
						hv -= a[ref[p]][p * k + q]
					}
					if ref[p] == ref[q] {
						hv += a[ref[p]][p * k + q]
					}
					h[row * size + slot(np, ref, q, m)] = hv
				}
			}
		}
	}
	return gibbs
}

// newton minimises G over the phases' mole numbers from a split's, as the
// two-phase flash's Newton does with two.
fn (mu Multi) newton(beta []f64, x [][]f64) !MultiSplit {
	np := beta.len
	k := mu.present.len
	mut ref := []int{len: k}
	for a, i in mu.present {
		for j in 1 .. np {
			if beta[j] * x[j][i] > beta[ref[a]] * x[ref[a]][i] {
				ref[a] = j
			}
		}
	}
	mut v := []f64{len: (np - 1) * k}
	for a, i in mu.present {
		for j in 0 .. np {
			if j != ref[a] {
				v[slot(np, ref, a, j)] = beta[j] * x[j][i]
			}
		}
	}
	n0 := mu.amounts(np, ref, v)!
	mut scale := 0.0
	for j in 0 .. np {
		total := numeric.sum(n0[j])
		xj := n0[j].map(it / total)
		st := mu.state(xj)!
		for i in mu.present {
			scale += math.abs(n0[j][i] * (math.log(xj[i]) + st.ln_phi[i]))
		}
	}
	sv := v.map(math.sqrt(it))
	mut smallest := sv[0]
	for sa in sv {
		smallest = math.min(smallest, sa)
	}
	// The scaled gradient is s·g, so |g| ≤ tol everywhere is |s·g| ≤ s·tol:
	// held at the smallest s, which leaves the stop to an unshifted step too
	// small to move anything, or a fall G's rounding can't show; split then
	// holds the answer to tol in every ln f itself.
	sol := numeric.minimize(MultiGibbs{mu, np, ref, sv},
		[]f64{len: v.len, init: v[index] / sv[index]},
		gtol:     mu.params.tol * smallest
		fscale:   scale
		max_iter: 50
	)!
	n := mu.amounts(np, ref, []f64{len: v.len, init: sol.x[index] * sv[index]})!
	mut fractions := []f64{}
	mut xs := [][]f64{}
	mut st := []thermo.State{}
	for j in 0 .. np {
		total := numeric.sum(n[j])
		fractions << total
		xs << n[j].map(it / total)
		st << mu.state(xs[j])!
	}
	return MultiSplit{fractions, xs, st, .newton, sol.iterations}
}

// grown is the split of three from a split of two one of whose phases is
// unstable: its two phases and the unstable trial's composition, the trial
// with 1% of the feed to begin with. At no moles, a component that lives
// almost only in the trial would be in no phase at all, its Eᵢ in the
// multiphase Rachford–Rice underflowing and its Hessian overflowing.
fn (mu Multi) grown(x [][]f64, beta []f64, w []f64, stab Stability) !FlashResult {
	if x.len >= max_phases {
		return fail(.more_phases, 'a phase of the ${x.len}-phase split at ${mu.t} K and ${mu.p} Pa is unstable: a ${x.len + 1}th phase, and the flash finds ${max_phases} at most')
	}
	mut xs := x.map(it.clone())
	xs << w.clone()
	mut b := beta.map(0.99 * it)
	b << 0.01
	s := mu.split(xs, b)!
	return mu.result(s, stab)
}

// result verifies each phase of a split, with the others as trials, and
// grows it by an unstable one's trial up to max_phases; then labels it.
fn (mu Multi) result(s MultiSplit, stab Stability) !FlashResult {
	if mu.params.verify {
		for j in 0 .. s.x.len {
			others := [][]f64{len: s.x.len, init: s.x[index]}.filter(it != s.x[j])
			check := stability(mu.m, mu.t, mu.p, s.x[j],
				trials: others
				tol:    verify_tol(mu.params.tol)
			)!
			if !check.stable {
				inc := check.incipient() or {
					return fail(.more_phases, 'a phase of the split at ${mu.t} K and ${mu.p} Pa is unstable')
				}
				return mu.grown(s.x, s.beta, inc.w, stab)
			}
		}
	}
	return labelled(mu.m, mu.t, mu.p, mu.z, s.st, s.beta, s.path, s.iterations, stab)
}

// labelled orders a flash's phases by increasing mass density and names
// them (M6's D9): each by its model's phase_label (an equation of state's
// by Π and the four-covolume rule, a γ–φ fluid's by the model that made
// it); of the phases called vapour, only the least dense stays so, the rest
// being liquids; and of the liquids, the less dense is the liquid and the
// denser liquid2. A third liquid has no name in D9, and is refused. β is
// the vapour's fraction, or 0 with none.
fn labelled(m thermo.PhaseModel, t f64, p f64, z []f64, states []thermo.State, fractions []f64, path Path, iterations int, stab Stability) !FlashResult {
	// By insertion, not sort_with_compare: a comparison closure would be one
	// more allocation that V never frees, on every flash.
	dens := states.map(mass_density(m, it))
	mut order := []int{len: states.len, init: index}
	for a in 1 .. order.len {
		mut b := a
		for b > 0 && dens[order[b]] < dens[order[b - 1]] {
			order[b], order[b - 1] = order[b - 1], order[b]
			b--
		}
	}
	mut phases := []Phase{}
	mut vapour := false
	mut liquids := 0
	mut beta := 0.0
	for j in order {
		mut label := m.phase_label(states[j])
		if label == .vapor && vapour {
			label = .liquid
		}
		if label == .vapor {
			vapour = true
			beta = fractions[j]
		} else {
			liquids++
			if liquids > 2 {
				return fail(.more_phases, 'three liquids at ${t} K and ${p} Pa, and D9 names two')
			}
			label = if liquids == 1 { thermo.PhaseLabel.liquid } else { .liquid2 }
		}
		phases << Phase{label, states[j], fractions[j]}
	}
	return FlashResult{
		t:          t
		p:          p
		z:          z
		beta:       beta
		phases:     phases
		path:       path
		iterations: iterations
		stability:  stab
	}
}
