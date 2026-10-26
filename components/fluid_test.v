module components

import json2
import math
import os
import thermo

// M2 exit criterion 9, end to end: db.fluid builds a working thermo.Fluid,
// and with database constants it reproduces M1's n-butane results.

const db9 = load_default() or { panic(err) }

fn root_of(phase string) thermo.Root {
	return if phase == 'liquid' { thermo.Root.liquid } else { thermo.Root.vapor }
}

fn test_setup_wires_the_database_into_thermo() {
	s := db9.setup(.pr, ['methane', 'n-butane'])!
	assert s.records.map(it.ids.cas) == ['74-82-8', '106-97-8']
	// kᵢⱼ from pr.ipd: methane/n-butane's first listing, line 133.
	assert s.kij.k == [0.0, 0.0244, 0.0244, 0.0]
	assert (s.kij.pairs[0].used or { Kij{} }).line == 133
	assert s.alphas == ["the equation's own", "the equation's own"]
	// The same fluid built by hand, from the same numbers, gives the same
	// properties, bit for bit.
	mut comps := []thermo.Component{}
	mut ideal := []thermo.IdealGas{}
	for r in s.records {
		comps << r.component()!
		ideal << r.ideal_gas()!
	}
	eos := thermo.new_cubic(.pr, comps, kij: [0.0, 0.0244, 0.0244, 0.0])!
	by_hand := thermo.new_fluid(eos, ideal)!
	for root in [thermo.Root.liquid, .vapor] {
		st := thermo.solve_tp(s.fluid.model, 300.0, 10e5, [0.3, 0.7], root)!
		st2 := thermo.solve_tp(eos, 300.0, 10e5, [0.3, 0.7], root)!
		assert s.fluid.props(st) == by_hand.props(st2)
	}
	// db.fluid is the same thing without the record.
	fl := db9.fluid(.pr, ['methane', 'n-butane'])!
	assert fl.props(thermo.solve_tp(fl.model, 300.0, 10e5, [0.3, 0.7], .vapor)!).h == s.fluid.props(thermo.solve_tp(s.fluid.model,
		300.0, 10e5, [0.3, 0.7], .vapor)!).h
}

fn test_kij_options() {
	off := db9.setup(.pr, ['methane', 'n-butane'], kij: false)!
	assert off.kij.k == [0.0, 0.0, 0.0, 0.0] && off.kij.missing().len == 1
	hand := db9.setup(.pr, ['methane', 'n-butane'],
		overrides: [Kij{
			cas1: '106-97-8'
			cas2: '74-82-8'
			k:    0.0185
		}]
	)!
	assert hand.kij.k[1] == 0.0185 && hand.kij.pairs[0].source == .override
	// vdW has no kᵢⱼ bank, so its pairs are 0.
	vdw := db9.setup(.vdw, ['methane', 'n-butane'])!
	assert vdw.kij.k == [0.0, 0.0, 0.0, 0.0] && vdw.kij.bank == 'no kij file'
}

fn test_mathias_copeman_alpha() {
	s := db9.setup(.pr, ['n-butane', '3-chloropropene'], alpha: .mathias_copeman)!
	assert s.alphas[0] == 'Mathias–Copeman, mc-PR.lib line 99'
	// 3-chloropropene isn't in mc-PR.lib, so it keeps PR's own α, and says so.
	assert s.alphas[1] == "the equation's own: mc-PR.lib has no 3-Chloropropene"
	m := s.fluid.model
	if m is thermo.Cubic {
		assert m.alphas[0] == thermo.Alpha(thermo.MathiasCopeman{
			c1: 0.652753
			c2: 0.0998805
			c3: 0.0
		})
		assert m.alphas[1] == thermo.CubicKind.pr.default_alpha(m.comps[1].omega)
	} else {
		assert false, 'not a cubic'
	}
	// mc-PR.lib is fitted for PR.
	if _ := db9.setup(.srk, ['n-butane'], alpha: .mathias_copeman) {
		assert false, 'Mathias–Copeman for SRK'
	}
	// Mistakes in the request are errors.
	if _ := db9.setup(.pr, []) {
		assert false, 'an empty fluid'
	}
	if _ := db9.setup(.pr, ['n-butane', 'unobtanium']) {
		assert false, 'an unknown compound'
	}
	if _ := db9.setup(.pr, ['3967-55-3']) {
		assert false, 'a compound without Tc, Pc and ω'
	}
}

// M1's CoolProp fixture, as thermo/fixtures_test.v reads it.
struct M1Deviation {
	v_rel  f64
	h_dep  f64
	s_dep  f64
	cp_rel f64
	w_rel  f64
}

struct M1Point {
	t            f64
	p            f64
	phase        string
	v            f64
	h_dep        f64
	s_dep        f64
	cp0          f64
	cp           f64
	w            f64
	pr_deviation M1Deviation
}

struct M1Fixture {
	points []M1Point
}

fn test_m1_results_with_database_constants() {
	// M1 recorded PR's deviations from CoolProp's n-butane reference
	// equation at six states, with Smith, Van Ness & Abbott's constants
	// (Tc 425.1 K, ω 0.200). ChemSep's differ slightly (425.12 K, 0.199), so
	// the database-built PR must land close to those recorded deviations:
	// the same model, with the database's constants. As in M1, the ideal-gas
	// Cp is CoolProp's own at each point, so only the residual model is
	// compared.
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm1_coolprop_butane.json'))!
	fx := json2.decode[M1Fixture](text)!
	assert fx.points.len == 6
	eos := db9.fluid(.pr, ['n-butane'])!.model
	for pt in fx.points {
		at := '${pt.phase} ${pt.t} K, ${pt.p / 1e5} bar'
		st := thermo.solve_tp(eos, pt.t, pt.p, [1.0], root_of(pt.phase))!
		rp := thermo.residual_props(eos, st)
		fl := thermo.new_fluid(eos, [thermo.IdealGas{
			cp: thermo.Polynomial{
				terms: [thermo.PowerTerm{pt.cp0, 0}]
			}
		}])!
		pr := fl.props(st)
		m1 := pt.pr_deviation
		// Within what a 0.02 K and 0.001 change in the constants moves them.
		assert math.abs((st.v - pt.v) / pt.v - m1.v_rel) < 2e-3, '${at}: V'
		assert math.abs(rp.h - pt.h_dep - m1.h_dep) < 50.0, '${at}: Hʳ'
		assert math.abs(rp.s - pt.s_dep - m1.s_dep) < 0.15, '${at}: Sʳ'
		assert math.abs((pr.cp - pt.cp) / pt.cp - m1.cp_rel) < 2e-3, '${at}: Cp'
		assert math.abs((pr.w - pt.w) / pt.w - m1.w_rel) < 2e-3, '${at}: w'
	}
}

// The n-butane section of testdata/m2_coolprop.json.
struct ButanePsat {
	t f64
	p f64
}

struct ButaneState {
	label string
	phase string
	t     f64
	p     f64
	h     f64
}

struct ButaneRefs {
	psat   []ButanePsat
	states []ButaneState
}

struct ButaneFixture {
	n_butane ButaneRefs
}

fn butane_refs() !ButaneRefs {
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm2_coolprop.json'))!
	return json2.decode[ButaneFixture](text)!.n_butane
}

// recorded checks a deviation against the one recorded when this test was
// written: the same, and not larger.
fn recorded(what string, got f64, want f64, tol f64) {
	assert math.abs(got - want) <= tol, '${what}: deviation ${got}, recorded ${want}'
	assert math.abs(got) <= math.abs(want) + tol, '${what}: deviation grew from ${want} to ${got}'
}

fn test_psat_end_to_end() {
	eos := db9.fluid(.pr, ['n-butane'])!.model
	// Smith, Van Ness & Abbott's Example 3.9, measured: 9.4573 bar at 350 K.
	recorded('Psat(350 K) against SVA', thermo.saturation_pressure(eos, 350.0)! / 9.4573e5 - 1.0,
		0.0018438, 1e-6)
	// Against CoolProp's reference equation.
	want := {
		250.0: 0.0121818
		300.0: 0.0005105
		350.0: 0.0034937
		400.0: 0.0084255
	}
	refs := butane_refs()!
	assert refs.psat.len == 4
	for p in refs.psat {
		dev := thermo.saturation_pressure(eos, p.t)! / p.p - 1.0
		recorded('Psat(${p.t} K) against CoolProp', dev, want[p.t] or { panic(p.t) }, 1e-6)
	}
}

fn test_enthalpy_changes_end_to_end() {
	// Total enthalpy changes of n-butane from db.fluid (PR, ChemSep's
	// ideal-gas Cp) against CoolProp, as fractions of CoolProp's change. The
	// recorded deviations are PR's and ChemSep's, not reflux bugs: liquid
	// compression is where PR's poor liquid volumes show.
	fl := db9.fluid(.pr, ['n-butane'])!
	refs := butane_refs()!
	mut ours := map[string]f64{}
	mut theirs := map[string]f64{}
	for s in refs.states {
		st := thermo.solve_tp(fl.model, s.t, s.p, [1.0], root_of(s.phase))!
		ours[s.label] = fl.props(st).h
		theirs[s.label] = s.h
	}
	for change in [
		['vapor 300 K 1 bar', 'vapor 450 K 1 bar', '0.01123'],
		['liquid 250 K 20 bar', 'liquid 350 K 20 bar', '0.01034'],
		['saturated liquid 300 K', 'saturated vapor 300 K', '0.00551'],
		['liquid 300 K 5 bar', 'liquid 300 K 100 bar', '-0.13708'],
		['supercritical 450 K 50 bar', 'supercritical 500 K 50 bar', '-0.04451'],
		['liquid 300 K 5 bar', 'vapor 450 K 1 bar', '0.00525'],
	] {
		a, b := change[0], change[1]
		dh := theirs[b] - theirs[a]
		recorded('ΔH ${a} → ${b}', (ours[b] - ours[a] - dh) / dh, change[2].f64(), 1e-4)
	}
}
