// Three phases (M6 step 7, criterion 6) against thermo 0.6.1's FlashVLN,
// polished on its own phases (testdata/m6_three_phase.json, from
// tools/gen_m6_three_phase.py): NRTL on water and 1-butanol, and on water,
// 1-butanol and methanol, two liquids and a vapour and two liquids; and
// Peng–Robinson on water, n-decane and methane, and water and n-decane.
// Each answer is checked from scratch too, and a fourth phase is refused.
import components
import flash
import json2
import math
import os
import thermo

struct RefPhase {
	kind     string
	fraction f64
	x        []f64
}

struct RefFlashVln {
	phase_count int
}

struct Case {
	model    string
	system   []string
	t        f64 @[json: 'T']
	p        f64 @[json: 'P']
	z        []f64
	flashvln RefFlashVln
	phases   []RefPhase
}

struct PrComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

struct Fixture {
	pr_components []PrComponent
	pr_kij        [][]f64
	cases         []Case
}

fn fixture() !Fixture {
	return json2.decode[Fixture](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_three_phase.json'))!)!
}

fn pr_model(f Fixture, n int) !thermo.PhaseModel {
	comps := f.pr_components[..n].map(thermo.Component{
		name:  it.name
		tc:    it.tc
		pc:    it.pc
		omega: it.omega
		mw:    it.mw
	})
	mut kij := []f64{}
	for i in 0 .. n {
		for j in 0 .. n {
			kij << f.pr_kij[i][j]
		}
	}
	return thermo.by_roots(thermo.new_cubic(.pr, comps, kij: kij)!)
}

// checked verifies an answer from scratch: each phase on its phase of lower
// G, every fugacity within 1e-10 of the first phase's, the balance to
// 1e-12, the fractions summing to 1, and no composition below any phase's
// tangent plane, by a scan of the composition space (steps of 1/1000 for
// two components, 1/40 for three). It returns the scan's lowest TPD.
fn checked(m thermo.PhaseModel, r flash.FlashResult, at string) !f64 {
	nc := r.z.len
	mut sts := []thermo.State{}
	for ph in r.phases {
		st := m.phase(r.t, r.p, ph.state.n, .stable)!
		assert st.v == ph.state.v, '${at}: a phase not on its phase of lower G'
		sts << st
	}
	mut total := 0.0
	for ph in r.phases {
		total += ph.fraction
	}
	assert math.abs(total - 1.0) <= 1e-14, at
	for i in 0 .. nc {
		mut z := 0.0
		for k, ph in r.phases {
			z += ph.fraction * ph.state.n[i]
			a := r.phases[0].state.n
			// Below the trace floor a component's amount is zero to double
			// precision, and its fugacities can't agree (flash.trace_floor).
			if k > 0 && ph.state.n[i] >= flash.trace_floor && a[i] >= flash.trace_floor {
				f := math.log(ph.state.n[i]) + sts[k].ln_phi[i] - math.log(a[i]) - sts[0].ln_phi[i]
				assert math.abs(f) <= 1e-10, '${at}: ln f of ${i} in phase ${k} is ${f} off'
			}
		}
		assert math.abs(z - r.z[i]) <= 1e-12, '${at}: the balance of ${i}'
	}
	mut grid := [][]f64{}
	if nc == 2 {
		for k in 1 .. 1000 {
			grid << [f64(k) / 1000.0, 1.0 - f64(k) / 1000.0]
		}
	} else {
		for a in 0 .. 41 {
			for b in 0 .. 41 - a {
				grid << [f64(a) / 40.0, f64(b) / 40.0, f64(40 - a - b) / 40.0].map(math.max(it, 1e-12))
			}
		}
	}
	mut lowest := math.inf(1)
	for w in grid {
		sw := m.phase(r.t, r.p, w, .stable)!
		for k, ph in r.phases {
			x := ph.state.n
			mut tpd := 0.0
			for i in 0 .. nc {
				// A component absent from the phase gives its own tangent
				// plane no bound: the scan is of the components present.
				if x[i] < flash.trace_floor {
					continue
				}
				tpd += w[i] * (math.log(w[i]) + sw.ln_phi[i] - math.log(x[i]) - sts[k].ln_phi[i])
			}
			lowest = math.min(lowest, tpd)
		}
	}
	return lowest
}

fn test_three_phases_against_flashvln() {
	f := fixture()!
	db := components.load_default()!
	mut worst_f, mut worst_x, mut lowest := 0.0, 0.0, 0.0
	mut counted := map[string]int{}
	for c in f.cases {
		at := '${c.model} ${c.system} at ${c.t} K, ${c.p} Pa'
		m := if c.model == 'PR' {
			pr_model(f, c.z.len)!
		} else {
			db.gamma_phi(c.system)!.fluid.model
		}
		r := flash.pt(m, c.t, c.p, c.z)!
		assert c.flashvln.phase_count == c.phases.len, '${at}: FlashVLN has ${c.flashvln.phase_count} phases'
		assert r.phases.len == c.phases.len, '${at}: ${r.phases.len} phases, thermo ${c.phases.len}'
		counted['${c.model} ${r.phases.map(it.label.name()).join(', ')}']++
		for k, ph in r.phases {
			want := c.phases[k]
			worst_f = math.max(worst_f, math.abs(ph.fraction - want.fraction))
			for i in 0 .. c.z.len {
				worst_x = math.max(worst_x, math.abs(ph.state.n[i] - want.x[i]))
			}
			// A γ–φ fluid's phases are named by the model that made them, as
			// thermo's are by its kind; an equation of state's by Π (D9).
			if c.model == 'NRTL' {
				assert (ph.label == .vapor) == (want.kind == 'gas'), '${at}: phase ${k}'
			}
		}
		// D9: the lightest first, at most one vapour, and of two liquids
		// the denser liquid2.
		labels := r.phases.map(it.label)
		assert labels == [thermo.PhaseLabel.vapor, .liquid, .liquid2] || labels == [
			thermo.PhaseLabel.liquid,
			.liquid2,
		]
			|| labels == [thermo.PhaseLabel.vapor, .liquid], '${at}: ${labels}'
		assert r.beta == if labels[0] == .vapor { r.phases[0].fraction } else { 0.0 }
		lowest = math.min(lowest, checked(m, r, at)!)
	}
	println('three phases: ${counted}; fractions ${worst_f:.1e}, compositions ${worst_x:.1e}, lowest TPD ${lowest:.1e}')
	// Measured: within the reference's own convergence, 1e-14 in ln f, over
	// each phase's conditioning. Criterion 6 asks FlashVLN's answers.
	assert worst_f <= 1e-9 && worst_x <= 1e-9, '${worst_f} ${worst_x}'
	assert lowest >= -1e-10, 'a TPD of ${lowest}'
	assert counted['PR vapour, liquid, liquid2'] == 4 && counted['NRTL vapour, liquid, liquid2'] == 3
	assert counted['NRTL liquid, liquid2'] == 8 && counted['PR liquid, liquid2'] == 2
}

// A fourth phase is refused, with its reason (D8): water, n-decane, methane
// and a heavier hydrocarbon that decane won't mix with (k_ij 0.5, as with
// water), at 350 K and 5 bar: a vapour and three liquids.
fn test_a_fourth_phase_is_refused() {
	f := fixture()!
	mut comps := f.pr_components.map(thermo.Component{
		name:  it.name
		tc:    it.tc
		pc:    it.pc
		omega: it.omega
		mw:    it.mw
	})
	comps << thermo.Component{
		name:  'n-dodecane'
		tc:    658.0
		pc:    18.2e5
		omega: 0.576
		mw:    0.170338
	}
	kij := [0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.04, 0.5, 0.5, 0.04, 0.0, 0.04, 0.5, 0.5, 0.04, 0.0]
	m := thermo.by_roots(thermo.new_cubic(.pr, comps, kij: kij)!)
	if r := flash.pt(m, 350.0, 5e5, [0.4, 0.25, 0.15, 0.2]) {
		assert false, '${r.phases.map(it.label.name())}'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
		// The reason is the fourth phase, found by the split of three's own
		// verification, not three liquids named after a split of four.
		assert err.msg().contains('-phase split at 350.0 K and 500000.0 Pa is unstable: a 4th phase'), err.msg()
	}
}

// The sweep's point 4386, mercury and 1-methylnaphthalene by PR at 111.2 K
// and 1.49 bar: the split of three loses a phase that comes together with
// another, and must not return before the fractions are solved again. It
// returned with the balance 3.3e-3 off. Two liquids, each nearly pure.
fn test_a_phase_that_leaves_takes_nothing_with_it() {
	db := components.load_default()!
	m := db.fluid(.pr, ['mercury', '1-methylnaphthalene'])!.model
	z := [0.4180605441827563, 0.5819394558172438]
	r := flash.pt(m, 111.22294815822565, 149310.3429963657, z)!
	assert r.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2]
	lowest := checked(m, r, 'mercury and 1-methylnaphthalene')!
	assert lowest >= -1e-10, '${lowest}'
}

// The sweep's point 7199, six compounds at 174 K and 0.15 bar, whose trial
// phase holds a component almost alone: begun at no moles, that component
// was in no phase, and the multiphase Rachford–Rice's Hessian overflowed.
// Begun at 1%, the flash finds the split of three and its fourth phase.
fn test_a_new_phase_begins_with_moles() {
	db := components.load_default()!
	m := db.fluid(.pr, ['107-21-1', '1333-74-0o', '2532-58-3', '301-00-8', '74-90-8', '837-08-1'])!.model
	z := [0.13577159252155996, 0.4639482447800575, 0.04161989373543508, 0.010792028865167379,
		0.31988838930532576, 0.027979850792454342]
	if r := flash.pt(m, 174.20364778502181, 15219.488523282738, z) {
		assert false, '${r.phases.map(it.label.name())}'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
	}
}
