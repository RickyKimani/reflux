// Checks the stability test against testdata/m3_stability.json, written by
// tools/gen_m3_fixtures.py (run by hand with uv). For the binaries the
// reference is a brute-force scan of the tangent-plane distance over the
// whole composition range; for every state it's also Caleb Bell's `thermo`.
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct Source {
	generated string
	generator string
	python    string
	thermo    string
}

struct FixtureComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

struct ThermoResult {
	phases       int
	tol          f64
	betas        []f64
	compositions [][]f64
	error        string
}

struct Brute {
	tm_min f64
	w      []f64
}

struct StabilityCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	thermo     ThermoResult
	witness_tm ?f64
	brute      ?Brute
	stable     bool
}

struct StabilityFixture {
	source     Source
	note       string
	components []FixtureComponent
	cases      []StabilityCase
}

fn read_fixture() !StabilityFixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_stability.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py`.') }
	return json2.decode[StabilityFixture](text)!
}

fn model_for(f StabilityFixture, c StabilityCase) !thermo.Cubic {
	comps := c.components.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if c.eos == 'pr' { CubicKind.pr } else { CubicKind.srk }
	return thermo.new_cubic(kind, comps, kij: c.kij)!
}

// tpd computes the normalised tangent-plane distance of w against feed z
// directly from the model, independently of flash.stability:
// Σ wᵢ(ln wᵢ + ln φᵢ(w) − ln zᵢ − ln φᵢ(z)), each phase at its root of lower
// Gibbs energy.
fn tpd(m thermo.Model, t f64, p f64, z []f64, w []f64) !f64 {
	feed := thermo.solve_tp(m, t, p, z, .stable)!
	trial := thermo.solve_tp(m, t, p, w, .stable)!
	mut sum := 0.0
	for i in 0 .. z.len {
		if w[i] > 0.0 {
			sum += w[i] * (math.log(w[i]) + trial.ln_phi[i] - math.log(z[i]) - feed.ln_phi[i])
		}
	}
	return sum
}

fn test_stability_agrees_with_the_references() {
	// 513 states: binaries on grids of T, P and z, states near the
	// methane/n-butane critical locus (some 0.1 bar apart, where substitution
	// alone needs up to 465 steps), a ternary and a 10-component gas
	// condensate. Binaries are judged against a brute-force scan of the
	// whole composition range, the rest against thermo's FlashVL (which
	// agrees with the scan on every binary). For an unstable binary,
	// reflux's lowest tm is the scan's global minimum: they agree to about
	// 2e-14. Every trial converges within 60 steps: the near-critical ones
	// only because the second-order finish takes over.
	f := read_fixture()!
	assert f.cases.len == 513
	mut unstable := 0
	for c in f.cases {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		m := model_for(f, c)!
		s := flash.stability(m, c.t, c.p, c.z)!
		assert s.stable == c.stable, '${at}: reflux says stable ${s.stable}'
		if b := c.brute {
			if !c.stable {
				inc := s.incipient() or { panic('${at}: unstable with no incipient phase') }
				assert math.abs(inc.tm - b.tm_min) <= 1e-11, '${at}: tm ${inc.tm} vs the scan ${b.tm_min}'
			}
		}
		for tr in s.trials {
			assert tr.converged, '${at}: ${tr.start} did not converge'
			assert tr.iterations <= 60, '${at}: ${tr.start} took ${tr.iterations}'
		}
		if !s.stable {
			unstable++
		}
	}
	assert unstable == 166
}

fn test_the_iteration_budget() {
	// Recorded, and can't grow by more than 5%: 25,631 steps over the
	// fixture's 2,399 trials (Wilson's two and the pure-component ones),
	// 266 of which needed the second-order finish. Dominant-eigenvalue
	// acceleration is what keeps it there: without it the same trials take
	// 35,138 steps and 895 finishes.
	f := read_fixture()!
	mut total, mut finishes := 0, 0
	for c in f.cases {
		m := model_for(f, c)!
		s := flash.stability(m, c.t, c.p, c.z)!
		for tr in s.trials {
			total += tr.iterations
			if tr.newton {
				finishes++
			}
		}
	}
	assert total <= 26900, '${total} steps'
	assert finishes <= 290, '${finishes} second-order finishes'
}

fn test_an_unstable_verdict_comes_with_its_proof() {
	// For every state reflux calls unstable, the incipient phase it returns
	// has a negative tangent-plane distance, computed here from scratch.
	f := read_fixture()!
	for c in f.cases {
		if c.stable {
			continue
		}
		m := model_for(f, c)!
		s := flash.stability(m, c.t, c.p, c.z)!
		inc := s.incipient() or { panic('unstable with no incipient phase') }
		check := tpd(m, c.t, c.p, c.z, inc.w)!
		assert check < -1e-10, '${c.components} ${c.z} at ${c.t} K, ${c.p} Pa: tm ${check}'
		assert math.abs(check - inc.tm) <= 1e-12 * math.max(1.0, math.abs(check))
		// And it's a stationary point: every ln wᵢ + ln φᵢ(w) − dᵢ equals tm
		// (weighted by √wᵢ, like the gradient the iteration converges). Near
		// a critical point the minimum is so flat that tm* stops changing
		// measurably while the gradient is still 5e-9, so the finish's last
		// steps go by the gradient (numeric.minimize's polish).
		feed := thermo.solve_tp(m, c.t, c.p, c.z, .stable)!
		trial := thermo.solve_tp(m, c.t, c.p, inc.w, .stable)!
		for i, w in inc.w {
			term := math.log(w) + trial.ln_phi[i] - math.log(c.z[i]) - feed.ln_phi[i]
			assert math.sqrt(w) * math.abs(term - check) <= 1e-11, '${c.components} ${c.z} at ${c.t} K, ${c.p} Pa: component ${i}'
		}
	}
}

fn test_a_stable_feed_has_no_negative_stationary_point() {
	// Every trial on a stable feed ends trivial, or at a stationary point
	// with tm ≥ −1e-10, and incipient() has nothing to offer.
	f := read_fixture()!
	mut trivial, mut nontrivial := 0, 0
	for c in f.cases {
		if !c.stable {
			continue
		}
		m := model_for(f, c)!
		s := flash.stability(m, c.t, c.p, c.z)!
		assert s.incipient() == none
		for tr in s.trials {
			if tr.trivial {
				trivial++
			} else {
				nontrivial++
				assert tr.tm >= -1e-10, '${c.components} ${c.z} at ${c.t} K, ${c.p} Pa: ${tr.start} tm ${tr.tm}'
			}
		}
	}
	assert trivial > 0 && nontrivial > 0
}

fn test_the_trivial_solution_is_flagged_and_never_offered() {
	// A trial flagged trivial ended at the feed; one not flagged didn't.
	// The trivial ones never become the incipient phase.
	f := read_fixture()!
	for c in f.cases {
		m := model_for(f, c)!
		s := flash.stability(m, c.t, c.p, c.z)!
		for tr in s.trials {
			mut gap := 0.0
			for i, z in c.z {
				gap = math.max(gap, math.abs(math.log(tr.w[i] / z)))
			}
			if tr.trivial {
				assert gap < 1e-4, '${c.z} at ${c.t} K, ${c.p} Pa: trivial ${tr.start} is ${gap} from z'
				assert math.abs(tr.tm) < 1e-10
			} else {
				assert gap > 1e-6
			}
		}
		if inc := s.incipient() {
			assert !inc.trivial
		}
	}
}

fn test_the_phases_of_a_split_are_stable() {
	// Re-checking a split (D3): each phase of thermo's equilibrium split,
	// taken as a feed, is stable, and the other phase given as a trial is a
	// stationary point with tm ≈ 0, the common tangent plane.
	f := read_fixture()!
	mut checked := 0
	for c in f.cases {
		if c.thermo.phases != 2 || c.components.len > 3 {
			continue
		}
		m := model_for(f, c)!
		x, y := c.thermo.compositions[0], c.thermo.compositions[1]
		for pair in [[x, y], [y, x]] {
			s := flash.stability(m, c.t, c.p, pair[0], trials: [pair[1]])!
			assert s.stable, '${c.components} at ${c.t} K, ${c.p} Pa: a phase of the split is unstable'
			given := s.trials.filter(it.start == .given)[0]
			assert !given.trivial && math.abs(given.tm) < 1e-9, 'tm ${given.tm}'
		}
		checked++
	}
	assert checked > 100
}

fn test_stability_edge_cases() {
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	butane := Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	decane := Component{
		name:  'n-decane'
		tc:    617.7
		pc:    21.10e5
		omega: 0.492
		mw:    0.142285
	}
	m := thermo.new_cubic(.pr, [methane, butane, decane])!
	// A pure component can't split by composition: stable, with no trials.
	s := flash.stability(m, 250.0, 30e5, [0.0, 1.0, 0.0])!
	assert s.stable && s.trials.len == 0
	// An absent component stays absent in every trial, and the rest is the
	// methane/n-butane binary, unstable here.
	b := flash.stability(m, 250.0, 30e5, [0.4, 0.6, 0.0])!
	assert !b.stable
	for tr in b.trials {
		assert tr.w[2] == 0.0
	}
	// z is normalised, so its scale doesn't matter.
	scaled := flash.stability(m, 250.0, 30e5, [4.0, 6.0, 0.0])!
	assert scaled.stable == b.stable
	assert math.abs(scaled.incipient()?.tm - b.incipient()?.tm) <= 1e-13
	// Bad input is a FlashError, never a panic.
	for args in [[250.0, 30e5, 0.4, 0.6], [250.0, 30e5, 0.4, -0.6, 0.0], [250.0, 30e5, 0.0, 0.0,
		0.0], [0.0, 30e5, 0.4, 0.6, 0.0], [250.0, -1.0, 0.4, 0.6, 0.0]] {
		flash.stability(m, args[0], args[1], args[2..]) or {
			assert err is flash.FlashError
			assert (err as flash.FlashError).kind == .bad_input
			continue
		}
		assert false, '${args} should be refused'
	}
	if _ := flash.stability(m, 250.0, 30e5, [0.4, 0.6, 0.0], trials: [[0.5, 0.5]]) {
		assert false, 'a given trial of the wrong length should be refused'
	}
	if _ := flash.stability(m, 250.0, 30e5, [0.4, 0.6, 0.0], max_steps: 0) {
		assert false, 'max_steps = 0 should be refused'
	}
}

fn test_a_trial_out_of_steps() {
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	butane := Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	m := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.02, 0.02, 0.0])!
	// A stable feed whose trials need more than one step: out of steps
	// without proving anything, it's a FlashError, never "stable".
	flash.stability(m, 400.0, 10e5, [0.5, 0.5], max_steps: 1) or {
		assert err is flash.FlashError
		assert (err as flash.FlashError).kind == .no_convergence
		// With room to finish, the same feed is stable.
		assert flash.stability(m, 400.0, 10e5, [0.5, 0.5])!.stable
		return
	}
	assert false, 'one step should not settle this feed'
}

fn test_an_unconverged_trial_still_proves_instability() {
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	butane := Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	m := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.02, 0.02, 0.0])!
	// Two steps from Wilson's vapour-like start are far from converged, but
	// already below the tangent plane: any W with tm < 0 is a proof.
	s := flash.stability(m, 250.0, 30e5, [0.4, 0.6], max_steps: 2)!
	assert !s.stable
	inc := s.incipient()?
	assert !inc.converged && inc.tm < -1e-10
	assert tpd(m, 250.0, 30e5, [0.4, 0.6], inc.w)! < -1e-10
}

fn test_pure_trials_find_what_wilson_misses() {
	// The methane-rich phase of a water/n-decane/methane split at 350 K and
	// 5 bar (PR, a water–hydrocarbon kᵢⱼ of 0.5). Wilson's vapour- and
	// liquid-like trials never reach a water-rich liquid and call it stable,
	// but nearly pure water lies far below its tangent plane: a third phase.
	// A trial phase of nearly pure water finds it.
	water := Component{
		name:  'water'
		tc:    647.1
		pc:    220.55e5
		omega: 0.345
		mw:    0.018015
	}
	decane := Component{
		name:  'n-decane'
		tc:    617.7
		pc:    21.10e5
		omega: 0.492
		mw:    0.142285
	}
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	m := thermo.new_cubic(.pr, [water, decane, methane],
		kij: [0.0, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.04, 0.0]
	)!
	x := [0.7066684313380365, 0.00809338270313943, 0.28523818595882405]
	assert flash.stability(m, 350.0, 5e5, x, pure_trials: false)!.stable
	s := flash.stability(m, 350.0, 5e5, x)!
	assert !s.stable
	inc := s.incipient()?
	assert inc.start == .pure && inc.w[0] > 0.99, '${inc.w}'
	check := tpd(m, 350.0, 5e5, x, inc.w)!
	assert check < -1.0 && math.abs(check - inc.tm) <= 1e-10
}

fn test_wilson_k() {
	// Kᵢ = (Pcᵢ/P)·exp(5.373(1 + ωᵢ)(1 − Tcᵢ/T)): at T = Tc and P = Pc, K = 1.
	c := Component{
		name:  'x'
		tc:    400.0
		pc:    40e5
		omega: 0.3
		mw:    0.05
	}
	assert flash.wilson_k([c], 400.0, 40e5) == [1.0]
	k := flash.wilson_k([c], 300.0, 10e5)[0]
	assert math.abs(k - 4.0 * math.exp(5.373 * 1.3 * (1.0 - 4.0 / 3.0))) <= 1e-15 * k
}
