// Module equilibrium is M9's chemical equilibrium (D9, D10): Gibbs
// minimisation over every species of one phase, ideal gas or a cubic, and
// the equilibrium over chosen reactions, each by the RAND method in log
// space (rand.v). An answer is returned only once it's shown right: its
// constraints met, every species' chemical potential its constraints'
// (or every reaction's Q its K), and, for a cubic, the phase stable.
module equilibrium

import flash
import kinetics
import math
import numeric
import thermo

// Failure is what kind of failure an EquilibriumError is.
pub enum Failure {
	bad_input     // a state, an amount or a species that can't be
	no_data       // a species without the formula or ideal-gas data it needs
	dependent     // reactions that aren't independent
	not_converged // the minimisation didn't converge, or its answer failed a check
	two_phases    // the answer would split into phases, which M6 brings
}

// EquilibriumError is every failure here.
pub struct EquilibriumError {
	Error
pub:
	kind   Failure
	reason string
}

// msg is the reason, after the module's name.
pub fn (e EquilibriumError) msg() string {
	return 'equilibrium: ${e.reason}'
}

fn fail(kind Failure, reason string) IError {
	return EquilibriumError{
		kind:   kind
		reason: reason
	}
}

// Mixture is what's at equilibrium: the species, with their formulas and
// ideal-gas data, and the residual model of their one phase, in the
// species' order, or none for the ideal gas.
pub struct Mixture {
pub:
	species []kinetics.Species
	model   ?thermo.Model
}

// Result is an equilibrium: each species' amount, mole fraction and ln φ,
// and the constraints' potentials π (the elements', for a Gibbs
// minimisation), with the RAND steps and φ updates it took.
pub struct Result {
pub:
	n      []f64 // mol
	y      []f64
	ln_phi []f64
	pi     []f64
	steps  int
	outer  int
}

// phi_tol is how far ln φ may move between the last two updates: an answer
// to rounding.
const phi_tol = 1e-12

// max_outer bounds the φ updates (D9).
const max_outer = 100

// check_state refuses what no equilibrium is at.
fn check_state(m Mixture, t f64, p f64, n0 []f64) ! {
	if !(t > 0) || math.is_inf(t, 1) || !(p > 0) || math.is_inf(p, 1) {
		return fail(.bad_input, 'no equilibrium at T = ${kinetics.sig(t)} K and P = ${kinetics.sig(p)} Pa')
	}
	if n0.len != m.species.len {
		return fail(.bad_input, '${n0.len} amounts for ${m.species.len} species')
	}
	mut total := 0.0
	for i, x in n0 {
		if !(x >= 0) || math.is_inf(x, 1) {
			return fail(.bad_input, "${m.species[i].name}'s amount, ${kinetics.sig(x)} mol, isn't one")
		}
		total += x
	}
	if total <= 0 {
		return fail(.bad_input, 'nothing to bring to equilibrium')
	}
	if model := m.model {
		if model.components().len != m.species.len {
			return fail(.bad_input, 'a model of ${model.components().len} components for ${m.species.len} species')
		}
	}
}

// standard is each species' μ°/RT at T: (H° − TS°)/RT from its ideal-gas
// data, as kinetics' standard properties take them (D7).
fn standard(m Mixture, t f64, active []int) ![]f64 {
	mut g := []f64{len: m.species.len}
	for j in active {
		s := m.species[j]
		gas := s.gas or { return fail(.no_data, '${s.name} has no ideal-gas data') }
		dhf := gas.dhf or { return fail(.no_data, '${s.name} has no ΔH°f') }
		s0 := gas.s0 or { return fail(.no_data, '${s.name} has no S°') }
		h := dhf + gas.cp.int_cp(thermo.t_ref, t)
		e := s0 + gas.cp.int_cp_over_t(thermo.t_ref, t)
		if !math.is_finite(h) || !math.is_finite(e) {
			return fail(.no_data, "${s.name}'s Cp doesn't integrate from ${thermo.t_ref} to ${t} K")
		}
		g[j] = (h - t * e) / (kinetics.gas_constant * t)
	}
	return g
}

// independent keeps a largest set of linearly independent rows, in order,
// by elimination with partial pivoting: the rows are small whole numbers or
// simple fractions, so 1e-9 tells a dependent row from rounding.
fn independent(rows [][]f64) [][]f64 {
	mut kept := [][]f64{}
	mut reduced := [][]f64{}
	mut pivots := []int{}
	for row in rows {
		mut r := row.clone()
		for q, p in pivots {
			f := r[p] / reduced[q][p]
			if f != 0 {
				for j in 0 .. r.len {
					r[j] -= f * reduced[q][j]
				}
			}
		}
		mut best := -1
		for j, x in r {
			if math.abs(x) > 1e-9 && (best < 0 || math.abs(x) > math.abs(r[best])) {
				best = j
			}
		}
		if best >= 0 {
			kept << row
			reduced << r
			pivots << best
		}
	}
	return kept
}

// Constraints are a minimisation's rows over the species, what they hold,
// and which species can be present at all.
struct Constraints {
	rows   [][]f64
	b      []f64
	active []int
}

// active_rows drops what the feed can't make: a conserved row with every
// coefficient of one sign and nothing in it fixes its species at zero (an
// element the feed hasn't, or an inert species it hasn't), so those species
// go. extra are more conserved rows, in the span of rows, looked at for
// that alone: the elements, where rows are a reaction set's null space,
// whose basis mixes them. What's left of rows is reduced to independent
// rows over the species that can be present.
fn active_rows(rows [][]f64, extra [][]f64, n0 []f64) Constraints {
	s := n0.len
	mut alive := []bool{len: s, init: true}
	mut changed := true
	mut looked := rows.clone()
	looked << extra
	for changed {
		changed = false
		for row in looked {
			mut b := 0.0
			mut pos, mut neg := false, false
			for j in 0 .. s {
				if !alive[j] || row[j] == 0 {
					continue
				}
				b += row[j] * n0[j]
				if row[j] > 0 {
					pos = true
				} else {
					neg = true
				}
			}
			if (pos != neg) && math.abs(b) == 0 {
				for j in 0 .. s {
					if alive[j] && row[j] != 0 {
						alive[j] = false
						changed = true
					}
				}
			}
		}
	}
	active := []int{len: s, init: index}.filter(alive[it])
	mut over := [][]f64{}
	for row in rows {
		r := active.map(row[it])
		if r.any(it != 0) {
			over << r
		}
	}
	kept := independent(over)
	mut b := []f64{len: kept.len}
	for k, r in kept {
		for q, j in active {
			b[k] += r[q] * n0[j]
		}
	}
	return Constraints{kept, b, active}
}

// element_rows is the element matrix over the species, an element per row,
// and a row of its own for each species without a formula or held inert.
fn element_rows(m Mixture, inert []int) ![][]f64 {
	s := m.species.len
	mut zs := []int{}
	for i, sp in m.species {
		if i in inert {
			continue
		}
		f := sp.formula or { continue }
		for c in f.counts {
			if c.z !in zs {
				zs << c.z
			}
		}
	}
	zs.sort()
	mut rows := [][]f64{}
	for z in zs {
		mut row := []f64{len: s}
		for i, sp in m.species {
			if i in inert {
				continue
			}
			if f := sp.formula {
				row[i] = f64(f.count(z))
			}
		}
		rows << row
	}
	for i, sp in m.species {
		if i in inert || sp.formula == none {
			mut row := []f64{len: s}
			row[i] = 1
			rows << row
		}
	}
	return rows
}

// solve minimises over the active species, updating φ outside RAND for a
// cubic (D9) until it stops moving, and returns the amounts over every
// species with the last π.
fn solve(m Mixture, cons Constraints, g []f64, t f64, p f64, n0 []f64) !Result {
	s := m.species.len
	a := cons.active.len
	if a == 0 {
		return fail(.bad_input, 'no species can be present')
	}
	mut total0 := 0.0
	for j in cons.active {
		total0 += n0[j]
	}
	pr := Problem{
		c:    cons.rows
		b:    cons.b
		g:    cons.active.map(g[it])
		ln_p: math.log(p / thermo.p_ref)
	}
	mut pt := Point{
		ln_n:     []f64{len: a, init: math.log(total0 / f64(a))}
		ln_total: math.log(total0)
	}
	mut ln_phi := []f64{len: a}
	mut full_phi := []f64{len: s}
	mut outer := 0
	for {
		outer++
		pr.minimise(ln_phi, mut pt)!
		model := m.model or { break }
		n := amounts(cons, pt.ln_n, s)
		st := thermo.solve_tp(model, t, p, n, .stable) or {
			return fail(.not_converged, 'no phase at the answer: ${err.msg()}')
		}
		mut moved := 0.0
		for q, j in cons.active {
			moved = math.max(moved, math.abs(st.ln_phi[j] - ln_phi[q]))
			ln_phi[q] = st.ln_phi[j]
		}
		full_phi = st.ln_phi.clone()
		if moved <= phi_tol {
			break
		}
		if outer >= max_outer {
			return fail(.not_converged, 'φ moved by ${moved} after ${max_outer} updates')
		}
	}
	n := amounts(cons, pt.ln_n, s)
	total := numeric.sum(n)
	return Result{
		n:      n
		y:      n.map(it / total)
		ln_phi: full_phi
		pi:     pt.pi
		steps:  pt.steps
		outer:  outer
	}
}

// amounts puts the active species' amounts in their places, and zero for
// the species nothing in the feed can make. A species that a constraint
// holds alone, an inert one, is set to the feed's amount exactly, rather
// than to the minimisation's to rounding.
fn amounts(cons Constraints, ln_n []f64, s int) []f64 {
	mut n := []f64{len: s}
	for q, j in cons.active {
		n[j] = math.exp(ln_n[q])
	}
	for k, row in cons.rows {
		mut only := -1
		for q, c in row {
			if c != 0 {
				if only >= 0 {
					only = -2
					break
				}
				only = q
			}
		}
		if only >= 0 {
			n[cons.active[only]] = cons.b[k] / row[only]
		}
	}
	return n
}

// gibbs is the equilibrium over every species of the phase at T and P from
// the amounts n0: G at its minimum over the amounts that keep the feed's
// elements (D9). inert lists species held at their feed amounts; a species
// without a formula is inert anyway. For a cubic, the answer's phase must
// be stable: one that would split is refused, a second phase being M6's.
pub fn gibbs(m Mixture, t f64, p f64, n0 []f64, inert []int) !Result {
	check_state(m, t, p, n0)!
	rows := element_rows(m, inert)!
	cons := active_rows(rows, [], n0)
	reacting := cons.active.filter(it !in inert && m.species[it].formula != none)
	g := standard(m, t, reacting)!
	r := solve(m, cons, g, t, p, n0)!
	verify_balance(rows, n0, r)!
	verify_stationary(m, cons, g, p, r, inert)!
	verify_stable(m, t, p, r)!
	return r
}

// verify_balance holds the answer to every row of the feed, to rounding:
// 1e-12 of the row's largest term.
fn verify_balance(rows [][]f64, n0 []f64, r Result) ! {
	for row in rows {
		mut want, mut got, mut scale := 0.0, 0.0, 0.0
		for j, c in row {
			want += c * n0[j]
			got += c * r.n[j]
			scale = math.max(scale, math.abs(c * n0[j]))
		}
		if math.abs(got - want) > 1e-12 * scale {
			return fail(.not_converged, 'the answer holds ${got} of a conserved quantity where the feed has ${want}')
		}
	}
}

// stationary_tol is how far a species' μ/RT may be from its constraints'
// potentials at the answer: rounding in μ/RT values of order 100.
const stationary_tol = 1e-9

// verify_stationary holds every species that can be present to μⱼ/RT =
// Σᵢ cᵢⱼπᵢ, with ln φ at the answer: G at a minimum over what the
// constraints allow.
fn verify_stationary(m Mixture, cons Constraints, g []f64, p f64, r Result, inert []int) ! {
	ln_p := math.log(p / thermo.p_ref)
	for q, j in cons.active {
		if j in inert || m.species[j].formula == none {
			continue
		}
		mu := g[j] + ln_p + math.log(r.y[j]) + if r.ln_phi.len > 0 { r.ln_phi[j] } else { 0.0 }
		mut cpi := 0.0
		for i, row in cons.rows {
			cpi += row[q] * r.pi[i]
		}
		if math.abs(mu - cpi) > stationary_tol {
			return fail(.not_converged, "${m.species[j].name}'s μ/RT is ${mu}, its elements' potentials give ${cpi}")
		}
	}
}

// verify_stable refuses an answer whose one phase would split: one phase
// for now (D9), a second being M6's.
fn verify_stable(m Mixture, t f64, p f64, r Result) ! {
	model := m.model or { return }
	st := flash.stability(model, t, p, r.y) or {
		return fail(.not_converged, "the answer's phase can't be tested: ${err.msg()}")
	}
	if !st.stable {
		return fail(.two_phases, 'the equilibrium mixture at ${kinetics.sig(t)} K and ${kinetics.sig(p)} Pa would split into two phases, which needs M6')
	}
}

// over_reactions is the equilibrium over chosen reactions alone (D10, as
// amended in step 7): each reaction's Q its K, ln_k at 1 bar on fugacities,
// with every species no reaction touches at its feed amount. The reactions
// must be independent. It's the RAND minimisation with what the reactions
// conserve as its constraints, the null space of their stoichiometry, and a
// stand-in potential for each species, g = −νᵀ(ννᵀ)⁻¹ ln K, so that each
// reaction's Σⱼ νⱼgⱼ is −ln K: the same Newton on the same log variables as
// gibbs, so a product the feed hasn't is no trouble.
pub fn over_reactions(m Mixture, reactions []kinetics.Reaction, ln_k []f64, t f64, p f64, n0 []f64) !Result {
	check_state(m, t, p, n0)!
	if reactions.len == 0 || ln_k.len != reactions.len {
		return fail(.bad_input, '${ln_k.len} equilibrium constants for ${reactions.len} reactions')
	}
	d := kinetics.dependence(reactions, m.species.len) or { return fail(.bad_input, err.msg()) }
	if d.combinations.len > 0 {
		return fail(.dependent, 'the reactions must be independent, but ${d.combinations.map(it.describe(reactions)).join('; ')}')
	}
	s := m.species.len
	mut nu := [][]f64{}
	for r in reactions {
		mut row := []f64{len: s}
		for j, k in r.species {
			row[k] = r.nu[j]
		}
		nu << row
	}
	rows := conserved(nu, s)
	g := stand_in(nu, ln_k)
	// Every element is conserved too: kinetics.parse_equation refuses an
	// equation that doesn't balance, or has a species without a formula.
	cons := active_rows(rows, element_rows(m, [])!, n0)
	r := solve(m, cons, g, t, p, n0)!
	verify_balance(rows, n0, r)!
	verify_k(reactions, nu, ln_k, p, r)!
	verify_stable(m, t, p, r)!
	return r
}

// conserved is a basis of what the reactions conserve: the vectors c over
// the species with c·νᵣ = 0 for every reaction, the null space of ν, by
// reduction to row echelon form.
fn conserved(nu [][]f64, s int) [][]f64 {
	mut rows := [][]f64{len: nu.len, init: nu[index].clone()}
	mut pivots := []int{}
	mut r := 0
	for c in 0 .. s {
		if r == rows.len {
			break
		}
		mut best := -1
		for i in r .. rows.len {
			if math.abs(rows[i][c]) > 1e-12 && (best < 0 || math.abs(rows[i][c]) > math.abs(rows[best][c])) {
				best = i
			}
		}
		if best < 0 {
			continue
		}
		if best != r {
			top := rows[r].clone()
			rows[r] = rows[best].clone()
			rows[best] = top
		}
		lead := rows[r][c]
		for j in 0 .. s {
			rows[r][j] /= lead
		}
		for i in 0 .. rows.len {
			if i != r && rows[i][c] != 0 {
				f := rows[i][c]
				for j in 0 .. s {
					rows[i][j] -= f * rows[r][j]
				}
			}
		}
		pivots << c
		r++
	}
	mut basis := [][]f64{}
	for f in 0 .. s {
		if f in pivots {
			continue
		}
		mut v := []f64{len: s}
		v[f] = 1
		for i, p in pivots {
			v[p] = -rows[i][f]
		}
		basis << v
	}
	return basis
}

// stand_in is the least-norm potentials g with Σⱼ νᵣⱼgⱼ = −ln Kᵣ for every
// reaction: g = −νᵀw with (ννᵀ)w = ln K.
fn stand_in(nu [][]f64, ln_k []f64) []f64 {
	r := nu.len
	s := nu[0].len
	mut a := []f64{len: r * r}
	for i in 0 .. r {
		for k in 0 .. r {
			mut x := 0.0
			for j in 0 .. s {
				x += nu[i][j] * nu[k][j]
			}
			a[i * r + k] = x
		}
	}
	lu := numeric_lu(a, r)
	w := lu.solve(ln_k) or { []f64{len: r} }
	mut g := []f64{len: s}
	for j in 0 .. s {
		for i in 0 .. r {
			g[j] -= nu[i][j] * w[i]
		}
	}
	return g
}

// k_tol is how far a reaction's ln Q may be from its ln K at the answer.
const k_tol = 1e-9

// verify_k holds every reaction to ln Q = ln K at the answer, Q on
// fugacities at 1 bar.
fn verify_k(reactions []kinetics.Reaction, nu [][]f64, ln_k []f64, p f64, r Result) ! {
	ln_p := math.log(p / thermo.p_ref)
	for i, row in nu {
		mut ln_q := 0.0
		for j, v in row {
			if v == 0 {
				continue
			}
			ln_q += v * (math.log(r.y[j]) + ln_p + if r.ln_phi.len > 0 { r.ln_phi[j] } else { 0.0 })
		}
		if math.abs(ln_q - ln_k[i]) > k_tol {
			return fail(.not_converged, "${reactions[i].name}'s ln Q is ${ln_q} at the answer, its ln K ${ln_k[i]}")
		}
	}
}

// ln_k_gibbs is each reaction's ln K at 1 bar from the Gibbs energies (D7).
pub fn ln_k_gibbs(reactions []kinetics.Reaction, species []kinetics.Species, t f64) ![]f64 {
	mut out := []f64{}
	for r in reactions {
		out << r.standard(species, t)!.ln_k
	}
	return out
}
