module sweep

// The sweep's own checks, each shown to refuse a wrong answer. A check that
// passes every flash proves nothing until it's seen to fail.
import components
import flash
import math
import thermo

// first is the first of the seeded points with the given class.
fn first(db components.Database, class Class) !Point {
	for pt in points(candidates(db)!, 100) {
		if run(db, pt).class == class {
			return pt
		}
	}
	return error('no ${class} point in the first 100')
}

// refused asserts that check_split refuses r, for a reason that mentions
// want.
fn refused(m thermo.Model, pt Point, r flash.FlashResult, want string) ! {
	check_split(m, pt, r) or {
		assert err.msg().contains(want), err.msg()
		return
	}
	assert false, 'a wrong split passed: ${want}'
}

// with_phases is r with its phases' compositions replaced by y and x.
fn with_phases(m thermo.Model, r flash.FlashResult, y []f64, x []f64) !flash.FlashResult {
	return flash.FlashResult{
		...r
		phases: [
			flash.Phase{
				...r.phases[0]
				state: thermo.solve_tp(m, r.t, r.p, y, .stable)!
			},
			flash.Phase{
				...r.phases[1]
				state: thermo.solve_tp(m, r.t, r.p, x, .stable)!
			},
		]
	}
}

// overlapping_split is the first seeded split with a component above 0.05
// in both phases, where a small change can break one check alone.
fn overlapping_split(db components.Database) !(Point, thermo.Model, flash.FlashResult) {
	for pt in points(candidates(db)!, 200) {
		m := db.fluid(.pr, pt.cas)!.model
		r := flash.pt(m, pt.t, pt.p, pt.z) or { continue }
		if r.phases.len == 2 {
			y, x := r.phases[0].state.n, r.phases[1].state.n
			for i in 0 .. x.len {
				if math.min(x[i], y[i]) > 0.05 {
					return pt, m, r
				}
			}
		}
	}
	return error('no split with overlapping phases in the first 200 points')
}

fn test_check_split_refuses_a_wrong_split() {
	db := components.load_default()!
	pt, m, r := overlapping_split(db)!
	check_split(m, pt, r)!
	y, x := r.phases[0].state.n, r.phases[1].state.n
	// Moving δ of a component from the vapour into the liquid, as moles
	// rather than fractions, keeps the mass balance exact: x gains δ and y
	// loses c·δ, where c = (1 − β)/β. Do it to the two components best
	// present in both phases.
	c := (1.0 - r.beta) / r.beta
	mut order := []int{len: x.len, init: index}
	order.sort_with_compare(fn [x, y] (a &int, b &int) int {
		ma, mb := math.min(x[*a], y[*a]), math.min(x[*b], y[*b])
		return if ma > mb {
			-1
		} else if ma < mb {
			1
		} else {
			0
		}
	})
	k, j := order[0], order[1]
	// The mass balance, with β off by 1e-9: each component's is then off by
	// 1e-9·|yᵢ − xᵢ|, far above 1e-12.
	refused(m, pt, flash.FlashResult{ ...r, beta: r.beta + 1e-9 }, 'mass balance')!
	// The sums, with δ = 1e-12/max(1, c) moved: Σx or Σy is then 1e-12
	// from 1, while ln fₖ changes by δ/xₖ + c·δ/yₖ ≤ 1e-12·(1/xₖ + 1/yₖ),
	// below 1e-10 since both are above 0.05.
	assert math.min(x[k], y[k]) > 0.05
	s := 1e-12 / math.max(1.0, c)
	mut x1, mut y1 := x.clone(), y.clone()
	x1[k] += s
	y1[k] -= c * s
	refused(m, pt, with_phases(m, r, y1, x1)!, 'Σx')!
	// The fugacities, with δ of k moved for δ of j: the sums and the mass
	// balance hold, while ln fₖ changes by δ/xₖ + c·δ/yₖ ≥ 1e-3.
	d := 1e-3 * math.min(math.min(x[k], y[k] / c), math.min(x[j], y[j] / c))
	mut x2, mut y2 := x.clone(), y.clone()
	x2[k] += d
	x2[j] -= d
	y2[k] -= c * d
	y2[j] += c * d
	refused(m, pt, with_phases(m, r, y2, x2)!, 'ln fV')!
}

fn test_check_split_refuses_a_split_that_is_not_one() {
	db := components.load_default()!
	// Two copies of a stable feed pass every check but these two.
	pt := first(db, .single)!
	m := db.fluid(.pr, pt.cas)!.model
	r := flash.pt(m, pt.t, pt.p, pt.z)!
	st := r.phases[0].state
	two := flash.FlashResult{
		...r
		beta:   0.5
		phases: [
			flash.Phase{
				label:    .vapor
				state:    st
				fraction: 0.5
			},
			flash.Phase{
				label:    .liquid
				state:    st
				fraction: 0.5
			},
		]
	}
	refused(m, pt, two, 'the same')!
	refused(m, pt, flash.FlashResult{ ...two, beta: 1.0 }, 'β')!
}

fn test_check_split_refuses_an_unstable_phase() {
	// Unverified, the flash returns a split that has an unstable phase at a
	// point the verified flash refuses. Only stability catches it.
	db := components.load_default()!
	pt := first(db, .three_phase)!
	m := db.fluid(.pr, pt.cas)!.model
	r := flash.pt(m, pt.t, pt.p, pt.z, verify: false)!
	refused(m, pt, r, 'unstable')!
}

fn test_check_single_refuses_an_unstable_feed() {
	db := components.load_default()!
	pt := first(db, .split)!
	m := db.fluid(.pr, pt.cas)!.model
	r := flash.pt(m, pt.t, pt.p, pt.z)!
	check_single(m, pt, r) or {
		assert err.msg().contains('below its tangent plane'), err.msg()
		return
	}
	assert false, 'a two-phase feed passed as one phase'
}

fn test_prove_more_phases_needs_an_unstable_phase() {
	db := components.load_default()!
	three := first(db, .three_phase)!
	proof := prove_more_phases(db.fluid(.pr, three.cas)!.model, three)!
	assert proof.contains('tm = -'), proof
	// A verified split and a single phase prove nothing.
	split := first(db, .split)!
	if p := prove_more_phases(db.fluid(.pr, split.cas)!.model, split) {
		assert false, 'a stable split proved a third phase: ${p}'
	} else {
		assert err.msg().contains('neither phase'), err.msg()
	}
	single := first(db, .single)!
	if p := prove_more_phases(db.fluid(.pr, single.cas)!.model, single) {
		assert false, 'one phase proved a third: ${p}'
	} else {
		assert err.msg().contains('one phase'), err.msg()
	}
}

fn test_the_report_counts_every_class() {
	// The CI subset's 1,000 one-phase points, then one of each other class:
	// the CI column must count only the first 1,000.
	pt := Point{
		names: ['A', 'B']
		z:     [0.5, 0.5]
		t:     300.0
		p:     1e5
	}
	mut outcomes := []Outcome{len: ci_count, init: Outcome{pt, .single, ''}}
	for class in [Class.split, .three_phase, .near_critical, .unexplained] {
		outcomes << Outcome{pt, class, 'why'}
	}
	text := report(640, outcomes)
	assert text.contains('| Outcome | All 1004 | First 1000 (CI) |')
	assert text.contains('| One phase, checked stable | 1000 | 1000 |')
	assert text.contains('| Two phases, every invariant checked | 1 | 0 |')
	assert text.contains('| Explained: a third phase (M6) | 1 | 0 |')
	assert text.contains('| Explained: at a critical point | 1 | 0 |')
	assert text.contains('| **Unexplained** | **1** | **0** |')
	assert text.contains('## Unexplained: 1\n\n| Point | T (K) | P (bar) | Feed | Reason |\n|---|---|---|---|---|\n| 0 | 300.0 | 1.0 | A 0.5, B 0.5 | why |')
}

fn test_round_sig() {
	assert round_sig(-0.044454836393, 4) == -0.04445
	assert round_sig(-12.605847891038922, 4) == -12.61
	assert round_sig(-0.0013403364909060096, 4) == -0.00134
	assert round_sig(0.0, 4) == 0.0
}
