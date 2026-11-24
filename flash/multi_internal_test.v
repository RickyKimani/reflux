module flash

import math
import thermo

// pt_test's water, n-decane and methane by PR, a vapour and two liquids at
// 350 K and 5 bar.
fn vll() !(thermo.PhaseModel, []f64) {
	comps := [
		thermo.Component{
			name:  'water'
			tc:    647.1
			pc:    220.55e5
			omega: 0.345
			mw:    0.018015
		},
		thermo.Component{
			name:  'n-decane'
			tc:    617.7
			pc:    21.10e5
			omega: 0.492
			mw:    0.142285
		},
		thermo.Component{
			name:  'methane'
			tc:    190.6
			pc:    45.99e5
			omega: 0.012
			mw:    0.016043
		},
	]
	m := thermo.new_cubic(.pr, comps, kij: [0.0, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.04, 0.0])!
	return thermo.by_roots(m), [0.5, 0.3, 0.2]
}

// Newton on G over three phases runs only where substitution is slow,
// which no test state is, so its claims are checked here from scratch:
// its gradient and Hessian against Richardson's differences of G and of
// the gradient, at a point off the answer; and from a split nudged off the
// answer, it comes back to it.
fn test_the_gibbs_derivatives_match_differences() {
	m, z := vll()!
	mu := Multi{
		m:       m
		t:       350.0
		p:       5e5
		z:       z
		present: [0, 1, 2]
		params:  FlashParams{}
	}
	r := pt(m, 350.0, 5e5, z)!
	assert r.phases.len == 3
	// Each component's dependent phase, the one holding the most of it.
	mut ref := []int{len: 3}
	for i in 0 .. 3 {
		for k in 1 .. 3 {
			if r.phases[k].fraction * r.phases[k].state.n[i] > r.phases[ref[i]].fraction * r.phases[ref[i]].state.n[i] {
				ref[i] = k
			}
		}
	}
	mut v := []f64{}
	for i in 0 .. 3 {
		for k in 0 .. 3 {
			if k != ref[i] {
				v << r.phases[k].fraction * r.phases[k].state.n[i] * (1.0 - 0.01 * f64(i + k + 1))
			}
		}
	}
	n := v.len
	mut g := []f64{len: n}
	mut h := []f64{len: n * n}
	mu.gibbs(3, ref, v, mut g, mut h)!
	mut worst_g, mut worst_h := 0.0, 0.0
	for a in 0 .. n {
		step := 1e-4 * v[a]
		mut gs := []f64{}
		mut grads := [][]f64{}
		for d in [step, -step, step / 2.0, -step / 2.0] {
			mut w := v.clone()
			w[a] += d
			mut gw := []f64{len: n}
			mut hw := []f64{len: n * n}
			gs << mu.gibbs(3, ref, w, mut gw, mut hw)!
			grads << gw
		}
		// A trace amount, water's 1e-29 of decane, moves G far below its
		// rounding, so G's differences see only the others; the gradient's,
		// moving as 1/n, see every one.
		if v[a] > 1e-6 {
			fd := (4.0 * (gs[2] - gs[3]) / step - (gs[0] - gs[1]) / (2.0 * step)) / 3.0
			worst_g = math.max(worst_g, math.abs(g[a] - fd) / math.max(math.abs(g[a]), 1e-3))
		}
		for b in 0 .. n {
			// In a trace amount's column the other entries move the gradient
			// below its rounding; its diagonal, as 1/n, doesn't, and the
			// rest of its row is checked in the others' columns.
			if v[a] <= 1e-6 && b != a {
				continue
			}
			fh := (4.0 * (grads[2][b] - grads[3][b]) / step - (grads[0][b] - grads[1][b]) / (2.0 * step)) / 3.0
			worst_h = math.max(worst_h, math.abs(h[b * n + a] - fh) / math.max(math.abs(h[b * n + a]),
				1.0))
		}
	}
	println('Gibbs: gradient ${worst_g:.1e}, Hessian ${worst_h:.1e}')
	// Measured: 2.1e-6 and 1.6e-7, the differences' rounding: G's ulps, 4e-16
	// of its 2.8, over steps down to 1e-7, against a floor of 1e-3.
	assert worst_g <= 1e-5 && worst_h <= 1e-6, '${worst_g} ${worst_h}'
	// Newton from a substitution iterate, the fifth, as it would finish a
	// slow split, comes to the answer substitution reaches by itself. The
	// start is pt's: its split of two and the unstable phase's trial.
	two := pt(m, 350.0, 5e5, z, verify: false)!
	check := stability(m, 350.0, 5e5, two.phases[1].state.n, trials: [two.phases[0].state.n])!
	mut w := check.incipient()?.w
	if check.stable {
		check0 := stability(m, 350.0, 5e5, two.phases[0].state.n, trials: [two.phases[1].state.n])!
		w = check0.incipient()?.w
	}
	// The trial, nearly pure water, first: with one dependent phase for
	// every component it would hold the others' traces as differences.
	start := [w, two.phases[0].state.n, two.phases[1].state.n]
	fractions := [0.01, 0.99 * two.phases[0].fraction, 0.99 * two.phases[1].fraction]
	early := Multi{
		...mu
		newton_after: 5
	}.split(start, fractions)!
	assert early.path == .newton, '${early.path}'
	slow := mu.split(start, fractions)!
	assert slow.path == .substitution
	for k in 0 .. 3 {
		assert math.abs(early.beta[k] - slow.beta[k]) <= 1e-9, '${early.beta} ${slow.beta}'
	}
	// A third phase that isn't there leaves: methane and decane alone are a
	// vapour and a liquid, and a phase between them, at 1%, is held at no
	// moles until it goes, leaving their split, as PT finds it.
	vl := pt(m, 350.0, 5e5, [0.0, 0.6, 0.4])!
	two_only := Multi{
		...mu
		z:       [0.0, 0.6, 0.4]
		present: [1, 2]
	}
	extra := two_only.split([vl.phases[0].state.n, vl.phases[1].state.n, [0.0, 0.5, 0.5]],
		[0.99 * vl.phases[0].fraction, 0.99 * vl.phases[1].fraction, 0.01])!
	assert extra.beta.len == 2, '${extra.beta}'
	assert math.abs(extra.beta[0] - vl.phases[0].fraction) <= 1e-9, '${extra.beta}'
	// A phase that stays away leaves too, after multi_absent_steps steps at
	// no moles: water at 1% of the feed all dissolves in the vapour, and a
	// trial of nearly pure water, begun at 1%, goes to no moles and to a
	// water-rich liquid of its own, never the same as another phase.
	zw := [0.01, 0.59, 0.4]
	vlw := pt(m, 350.0, 5e5, zw)!
	assert vlw.phases.len == 2
	dry := Multi{
		...mu
		z: zw
	}
	gone := dry.split([vlw.phases[0].state.n, vlw.phases[1].state.n, [0.98, 0.01, 0.01]],
		[0.99 * vlw.phases[0].fraction, 0.99 * vlw.phases[1].fraction, 0.01])!
	assert gone.beta.len == 2, '${gone.beta}'
	assert math.abs(gone.beta[0] - vlw.phases[0].fraction) <= 1e-9, '${gone.beta}'
}

// D9's labels from states alone: of two vapour-like phases only the less
// dense stays vapour; liquids lighter first, the denser liquid2; a third
// liquid refused; β the vapour's fraction, or 0.
fn test_the_labels() {
	m, z := vll()!
	r := pt(m, 350.0, 5e5, z)!
	sts := r.phases.map(it.state)
	fr := r.phases.map(it.fraction)
	// The same three in another order come out in density's order.
	back := labelled(m, 350.0, 5e5, z, [sts[2], sts[0], sts[1]], [fr[2], fr[0], fr[1]], .substitution,
		0, r.stability)!
	assert back.phases.map(it.label) == [thermo.PhaseLabel.vapor, .liquid, .liquid2]
	assert back.phases.map(it.state.v) == sts.map(it.v) && back.beta == fr[0]
	// Two copies of the vapour: the denser copy can't be a second vapour.
	two := labelled(m, 350.0, 5e5, z, [sts[0], sts[0]], [0.5, 0.5], .substitution, 0, r.stability)!
	assert two.phases.map(it.label) == [thermo.PhaseLabel.vapor, .liquid] && two.beta == 0.5
	// Three liquids have no name.
	if _ := labelled(m, 350.0, 5e5, z, [sts[1], sts[2], sts[2]], [0.3, 0.3, 0.4], .substitution, 0,
		r.stability)
	{
		assert false, 'three liquids'
	} else {
		assert err is FlashError && (err as FlashError).kind == .more_phases, err.msg()
	}
	// No vapour: β is 0.
	ll := labelled(m, 350.0, 5e5, z, [sts[1], sts[2]], [0.4, 0.6], .substitution, 0, r.stability)!
	assert ll.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2] && ll.beta == 0.0
}
