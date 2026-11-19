module equilibrium

import expr
import kinetics
import math
import numeric
import thermo

// The consistency check (kinetics.md's check 4, criterion 4): from a state,
// the equilibrium the kinetics' reactions reach, and each reversible rate
// held to zero there and to running towards it at the state.
//
// A rate is judged by where it stops: along its own reaction from that
// equilibrium, the extent at which it's zero, and its reaction's ΔG/RT =
// ln(Q/K) there, Q on fugacities and K from the Gibbs energies. That's zero
// for a rate consistent with the thermodynamics. A K fitted apart from the
// Gibbs energies stops where its own Q/K is 1, ln(K_fitted/K) from
// equilibrium; a rate on partial pressures in a gas that isn't ideal stops
// where Q on pressures is K, Σν ln φ from it; and a rate typed with a wrong
// species or a missing (P/p°)^Δn stops wherever that puts it.

// stop_tol is how far from zero a rate's stop may be, in ΔG/RT, and still be
// at equilibrium. A consistent rate is zero where its own Q is K, whatever
// the equilibrium's accuracy, so what's left is rounding in its Q/K near 1
// and the bisection's last step: 2.7e-15 at worst on Xu and Froment's. A
// constant typed to five figures is 1e-5 off, so 1e-10 tells the two apart.
const stop_tol = 1e-10

// stop_limit is how far a rate may stop from equilibrium before the check
// warns rather than notes: a factor of 2 in Q/K, as a fitted K is judged
// (kinetics' check_fitted).
const stop_limit = math.ln2

// Judged is one rate's consistency with equilibrium.
pub struct Judged {
pub:
	rate     int    // its index among the kinetics' rates
	reaction int    // its reaction's index
	found    bool   // whether it's zero anywhere along its reaction from the equilibrium
	stops    f64    // its reaction's ΔG/RT where it's zero, if found
	at       []f64  // the amounts there, mol, if found
	affinity f64    // its reaction's ΔG/RT at the state
	at_state f64    // the rate at the state
	skipped  string // why it isn't judged, where it isn't
}

// Consistency is the check's whole result: the equilibrium from the state,
// over the independent reactions, and each reversible rate judged by it.
pub struct Consistency {
pub:
	equilibrium Result
	over        []int // the reactions brought to equilibrium
	judged      []Judged
}

// Line evaluates rates and ΔG/RT along a reaction at fixed T and P.
struct Line {
	k     kinetics.Kinetics
	model ?thermo.Model
	t     f64
	p     f64
	ln_k  []f64   // each reaction's K at 1 bar, from the Gibbs energies
	nu    [][]f64 // each reaction's coefficients over the species
mut:
	pt  kinetics.Point
	buf expr.Buffer
	out []f64
}

// state is the rates' state at amounts n: their mole fractions, with φ and
// Z from the model, or the ideal gas's.
fn (l &Line) state(n []f64) !kinetics.State {
	total := numeric.sum(n)
	y := n.map(it / total)
	model := l.model or {
		return kinetics.State{
			t: l.t
			p: l.p
			y: y
		}
	}
	st := thermo.solve_tp(model, l.t, l.p, n, .stable)!
	return kinetics.State{
		t:   l.t
		p:   l.p
		y:   y
		phi: st.ln_phi.map(math.exp(it))
		z:   st.z
	}
}

// rate is rate j at a state, NaN where it can't be evaluated.
fn (mut l Line) rate(j int, s kinetics.State) f64 {
	l.k.fill(s, mut l.pt) or { return math.nan() }
	l.k.rates(&l.pt, mut l.buf, mut l.out)
	return l.out[j]
}

// rate_at is rate j at amounts n, NaN where it can't be evaluated.
fn (mut l Line) rate_at(j int, n []f64) f64 {
	s := l.state(n) or { return math.nan() }
	return l.rate(j, s)
}

// affinity is reaction i's ΔG/RT = ln(Q/K) at a state, Q on fugacities at
// 1 bar.
fn (l &Line) affinity(i int, s kinetics.State) f64 {
	ln_p := math.log(l.p / thermo.p_ref)
	mut a := -l.ln_k[i]
	for j, v in l.nu[i] {
		if v == 0 {
			continue
		}
		ln_phi := if s.phi.len > 0 { math.log(s.phi[j]) } else { 0.0 }
		a += v * (math.log(s.y[j]) + ln_p + ln_phi)
	}
	return a
}

// along is the amounts n moved by an extent x of reaction i.
fn (l &Line) along(n []f64, i int, x f64) []f64 {
	mut out := n.clone()
	for j, v in l.nu[i] {
		out[j] += v * x
	}
	return out
}

// reach is how far reaction i can run from n before a species runs out:
// the largest extent forwards, or the most negative backwards.
fn (l &Line) reach(n []f64, i int, forwards bool) f64 {
	mut most := math.inf(1)
	for j, v in l.nu[i] {
		if (forwards && v < 0) || (!forwards && v > 0) {
			most = math.min(most, n[j] / math.abs(v))
		}
	}
	return if forwards { most } else { -most }
}

// fractions are where stop probes, as fractions of the way to where a
// species runs out: from 2⁻⁵⁰, doubling to a half, then halving what's left
// to 2⁻⁵⁰ of it.
const fractions = probe_fractions()

fn probe_fractions() []f64 {
	mut out := []f64{}
	for e := 50; e >= 1; e-- {
		out << math.pow(2, -e)
	}
	for e in 2 .. 51 {
		out << 1 - math.pow(2, -e)
	}
	return out
}

// stop finds where rate j, of reaction i, is zero along its reaction from n,
// where it's r0: probing outwards both ways at once, the way it runs first,
// until its sign changes, then bisecting to rounding. Both ways, because a
// rate with its sign reversed is zero at equilibrium but runs away from it.
// It returns the amounts there, or none if the rate keeps its sign both ways
// until a species runs out or it can't be evaluated. A rate of exactly zero
// at n is taken as running backwards, and stops at once forwards or back.
fn (mut l Line) stop(j int, i int, n []f64, r0 f64) ?[]f64 {
	sign := if r0 > 0 { 1.0 } else { -1.0 }
	ends := [l.reach(n, i, sign > 0), l.reach(n, i, sign < 0)]
	mut open := ends.map(math.is_finite(it))
	mut last := [0.0, 0.0]
	for f in fractions {
		for q in 0 .. 2 {
			if !open[q] {
				continue
			}
			x := f * ends[q]
			r := l.rate_at(j, l.along(n, i, x))
			if !math.is_finite(r) {
				open[q] = false
			} else if r * sign > 0 {
				last[q] = x
			} else {
				return l.bisect(j, i, n, last[q], x, sign)
			}
		}
		if !open[0] && !open[1] {
			break
		}
	}
	return none
}

// bisect narrows a sign change of rate j between extents lo, where its sign
// is sign, and hi, where it isn't, to rounding, and returns the amounts at
// its middle.
fn (mut l Line) bisect(j int, i int, n []f64, from f64, to f64, sign f64) ?[]f64 {
	mut lo, mut hi := from, to
	for _ in 0 .. 200 {
		mid := 0.5 * (lo + hi)
		if mid == lo || mid == hi {
			break
		}
		r := l.rate_at(j, l.along(n, i, mid))
		if !math.is_finite(r) {
			return none
		}
		if r * sign > 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	return l.along(n, i, 0.5 * (lo + hi))
}

// reversible says whether a rate law can be zero short of a species running
// out: a structured law with a K, or an expression, which may be.
fn reversible(law kinetics.RateLaw) bool {
	return match law {
		kinetics.PowerLaw, kinetics.Lhhw, kinetics.EleyRideal {
			law.equilibrium.kind != .irreversible
		}
		kinetics.Expression { true }
	}
}

// fitted says whether a rate law's K is a fitted correlation.
fn fitted(law kinetics.RateLaw) bool {
	return match law {
		kinetics.PowerLaw, kinetics.Lhhw, kinetics.EleyRideal {
			law.equilibrium.kind == .correlation
		}
		kinetics.Expression { false }
	}
}

// consistency brings the kinetics' independent reactions to equilibrium from
// a state, at its T and P, from its mole fractions as amounts, in the
// model's gas or the ideal gas, and judges each reversible rate there and at
// the state.
pub fn consistency(k &kinetics.Kinetics, model ?thermo.Model, s kinetics.State) !Consistency {
	n := k.species.len
	if s.y.len != n {
		return fail(.bad_input, '${s.y.len} mole fractions for ${n} species')
	}
	d := kinetics.dependence(k.reactions, n) or { return fail(.bad_input, err.msg()) }
	dependent := d.combinations.map(it.reaction)
	over := []int{len: k.reactions.len, init: index}.filter(it !in dependent)
	chosen := over.map(k.reactions[it])
	m := Mixture{
		species: k.species
		model:   model
	}
	eq := over_reactions(m, chosen, ln_k_gibbs(chosen, k.species, s.t)!, s.t, s.p, s.y)!
	mut nu := [][]f64{}
	for r in k.reactions {
		mut row := []f64{len: n}
		for q, j in r.species {
			row[j] = r.nu[q]
		}
		nu << row
	}
	mut l := Line{
		k:     *k
		model: model
		t:     s.t
		p:     s.p
		ln_k:  ln_k_gibbs(k.reactions, k.species, s.t)!
		nu:    nu
		pt:    k.point()
		buf:   k.buffer()
		out:   []f64{len: k.rates.len}
	}
	// The state's own φ and Z, from the model, for its rates and ΔG.
	at_state := l.state(s.y.clone())!
	at_eq := l.state(eq.n)!
	mut judged := []Judged{}
	for j, law in k.laws {
		if !reversible(law) {
			continue
		}
		i := k.rated[j]
		absent := k.reactions[i].species.filter(eq.n[it] <= 0)
		if absent.len > 0 {
			judged << Judged{
				rate:     j
				reaction: i
				skipped:  'no ${k.species[absent[0]].name} at the equilibrium this state reaches'
			}
			continue
		}
		affinity := l.affinity(i, at_state)
		rs := l.rate(j, at_state)
		r0 := l.rate(j, at_eq)
		if !math.is_finite(r0) {
			judged << Judged{
				rate:     j
				reaction: i
				skipped:  'its rate is ${r0} at the equilibrium this state reaches'
			}
			continue
		}
		if at := l.stop(j, i, eq.n, r0) {
			judged << Judged{
				rate:     j
				reaction: i
				found:    true
				stops:    l.affinity(i, l.state(at)!)
				at:       at
				affinity: affinity
				at_state: rs
			}
		} else {
			judged << Judged{
				rate:     j
				reaction: i
				affinity: affinity
				at_state: rs
			}
		}
	}
	return Consistency{
		equilibrium: eq
		over:        over
		judged:      judged
	}
}

// consistent says whether a judged rate stops at equilibrium and runs
// towards it at the state.
fn (j Judged) consistent() bool {
	return j.skipped == '' && j.found && math.abs(j.stops) <= stop_tol && !j.runs_away()
}

// runs_away says whether a rate runs away from equilibrium at the state:
// forwards where its reaction's ΔG is above zero, or backwards where it's
// below, and not because it stops on the far side of equilibrium, as a
// fitted K can.
fn (j Judged) runs_away() bool {
	a := j.affinity
	if !math.is_finite(a) || !math.is_finite(j.at_state) || math.abs(a) <= stop_tol
		|| j.at_state * a <= 0 {
		return false
	}
	between := j.found && ((a > 0 && a < j.stops) || (a < 0 && a > j.stops))
	return !between
}

fn direction(r f64) string {
	return if r > 0 { 'forwards' } else { 'backwards' }
}

// stopped writes where a rate stops: as Q/K with ΔG/RT beside it, or as
// ΔG/RT alone where Q/K is 1 to five figures.
fn stopped(a f64) string {
	q := kinetics.sig(math.exp(a))
	if q == '1' {
		return 'ΔG/RT is ${plus(kinetics.sig(a))}'
	}
	return 'Q/K is ${q} (ΔG/RT ${plus(kinetics.sig(a))})'
}

fn plus(s string) string {
	return if s.starts_with('-') { s } else { '+${s}' }
}

// listed writes names as a list: "a", "a and b", "a, b and c".
fn listed(names []string) string {
	if names.len <= 1 {
		return names.join('')
	}
	return '${names[..names.len - 1].join(', ')} and ${names.last()}'
}

// findings are the check's findings: a note naming the rates that are
// consistent, and a finding for each that isn't, or can't be judged.
pub fn (c Consistency) findings(k &kinetics.Kinetics) []kinetics.Finding {
	mut out := []kinetics.Finding{}
	mut fine := []string{}
	for j in c.judged {
		name := k.reactions[j.reaction].name
		if j.skipped != '' {
			out << kinetics.Finding{
				severity: .note
				reaction: j.reaction
				message:  "${name}'s rate can't be held to equilibrium: ${j.skipped}"
			}
			continue
		}
		if j.consistent() {
			fine << "${name}'s"
			continue
		}
		if !j.found {
			expression := k.laws[j.rate] is kinetics.Expression
			out << kinetics.Finding{
				severity: if expression { kinetics.Severity.note } else { .warning }
				reaction: j.reaction
				message:  "${name}'s rate is never zero along its reaction from the equilibrium this state reaches: it runs one way until a species runs out${if expression {
					', irreversible as written'
				} else {
					''
				}}"
			}
		} else if math.abs(j.stops) > stop_tol {
			out << kinetics.Finding{
				severity: if math.abs(j.stops) > stop_limit {
					kinetics.Severity.warning
				} else {
					.note
				}
				reaction: j.reaction
				message:  "${name}'s rate isn't zero at the equilibrium this state reaches: ${if fitted(k.laws[j.rate]) {
					'on its fitted K '
				} else {
					''
				}}it stops where ${stopped(j.stops)}, so a reactor on it heads for a different equilibrium"
			}
		}
		if j.runs_away() {
			out << kinetics.Finding{
				severity: .warning
				reaction: j.reaction
				message:  "${name}'s rate runs ${direction(j.at_state)} at the state, where its reaction's ΔG/RT is ${plus(kinetics.sig(j.affinity))}: away from equilibrium"
			}
		}
	}
	if fine.len > 0 {
		out.prepend(kinetics.Finding{
			severity: .note
			message:  '${listed(fine)} ${if fine.len == 1 { 'rate is' } else { 'rates are' }} zero at the equilibrium this state reaches, over ${listed(c.over.map(k.reactions[it].name))}, and run${if fine.len == 1 {
				's'
			} else {
				''
			}} towards it at the state'
		})
	}
	return out
}

// check is the consistency check's findings at a state, or a note saying
// why it can't run there.
pub fn check(k &kinetics.Kinetics, model ?thermo.Model, s kinetics.State) []kinetics.Finding {
	c := consistency(k, model, s) or {
		return [
			kinetics.Finding{
				severity: .note
				message:  "the rates can't be held to equilibrium from this state: ${err.msg()}"
			},
		]
	}
	return c.findings(k)
}
