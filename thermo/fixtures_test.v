// Checks reflux against outside references stored in testdata/, written by
// tools/gen_m1_fixtures.py (run by hand with uv). The tests only read JSON:
// no Python, CoolProp or network access at test time.
import json2
import math
import os
import thermo { Component, CubicKind, Root }

// The structs mirror the JSON. json2.decode fills fields by name, so the
// field names are the JSON keys.
struct Source {
	generated string
	generator string
	python    string
	thermo    string
	coolprop  string
}

struct FixtureComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

fn (c FixtureComponent) to_component() Component {
	return Component{
		name:  c.name
		tc:    c.tc
		pc:    c.pc
		omega: c.omega
		mw:    c.mw
	}
}

struct ThermoPhase {
	phase  string
	z      f64
	v      f64
	ln_phi []f64
	h_dep  f64
	s_dep  f64
	cp_dep f64
}

struct ThermoCase {
	eos        string
	components []int
	kij        []f64
	n          []f64
	t          f64
	p          f64
	phases     []ThermoPhase
}

struct ThermoFixture {
	source     Source
	note       string
	components []FixtureComponent
	cases      []ThermoCase
}

struct Deviation {
	v_rel  f64
	h_dep  f64
	s_dep  f64
	cp_rel f64
	w_rel  f64
}

struct CoolPropPoint {
	t            f64
	p            f64
	phase        string
	v            f64
	z            f64
	h_dep        f64
	s_dep        f64
	cp0          f64
	cp           f64
	cv           f64
	w            f64
	pr_deviation Deviation
}

struct CoolPropFixture {
	source    Source
	note      string
	component FixtureComponent
	points    []CoolPropPoint
}

fn read_fixture[T](name string) !T {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', name)
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m1_fixtures.py`.') }
	return json2.decode[T](text)!
}

fn root_of(phase string) Root {
	return if phase == 'liquid' { Root.liquid } else { Root.vapor }
}

// agree checks a reflux value against a reference to a relative tol, with
// `floor` as the smallest scale to judge by, so that values near zero
// aren't held to a relative standard they can't meet.
fn agree(got f64, want f64, tol f64, floor f64) bool {
	return math.abs(got - want) <= tol * math.max(floor, math.abs(want))
}

fn test_fixtures_load_and_have_provenance() {
	a := read_fixture[ThermoFixture]('m1_thermo.json')!
	b := read_fixture[CoolPropFixture]('m1_coolprop_butane.json')!
	assert a.source.thermo != '' && a.source.generated != ''
	assert b.source.coolprop != '' && b.points.len > 0
	assert a.cases.len > 0
}

fn test_same_models_as_thermo() {
	// reflux and `thermo` implement the same equations with the same
	// constants, so they should agree to rounding: 1e-9 is the plan's exit
	// criterion 4.
	fx := read_fixture[ThermoFixture]('m1_thermo.json')!
	mut checked := 0
	for c in fx.cases {
		comps := c.components.map(fx.components[it].to_component())
		kind := if c.eos == 'pr' { CubicKind.pr } else { CubicKind.srk }
		eos := thermo.new_cubic(kind, comps, kij: c.kij)!
		for ph in c.phases {
			at := '${c.eos} ${comps.map(it.name)} ${ph.phase} at ${c.t} K, ${c.p / 1e5} bar'
			st := thermo.solve_tp(eos, c.t, c.p, c.n, root_of(ph.phase))!
			rp := thermo.residual_props(eos, st)
			rt := thermo.r_gas * c.t
			assert agree(st.z, ph.z, 1e-9, 0.0), '${at}: Z ${st.z} vs ${ph.z}'
			assert agree(st.molar_volume(), ph.v, 1e-9, 0.0), '${at}: V'
			for i, lp in ph.ln_phi {
				assert agree(st.ln_phi[i], lp, 1e-9, 1.0), '${at}: ln φ[${i}] ${st.ln_phi[i]} vs ${lp}'
			}
			assert agree(rp.h, ph.h_dep, 1e-9, rt), '${at}: Hʳ ${rp.h} vs ${ph.h_dep}'
			assert agree(rp.s, ph.s_dep, 1e-9, thermo.r_gas), '${at}: Sʳ ${rp.s} vs ${ph.s_dep}'
			assert agree(rp.cp, ph.cp_dep, 1e-9, thermo.r_gas), '${at}: Cpʳ ${rp.cp} vs ${ph.cp_dep}'
			checked++
		}
	}
	assert checked >= 20, 'only ${checked} phases checked'
}

fn test_pr_is_no_further_from_the_reference_equation_than_recorded() {
	// The fixture records how far PR (via thermo) was from CoolProp's
	// reference equation for n-butane when it was generated. reflux's PR is
	// the same model, so its deviations must match those to rounding. In
	// particular they can't have grown: any change that makes PR worse
	// fails here. The ideal-gas Cp is CoolProp's own at each point, so only
	// the residual model is being compared.
	fx := read_fixture[CoolPropFixture]('m1_coolprop_butane.json')!
	butane := fx.component.to_component()
	eos := thermo.new_cubic(.pr, [butane])!
	for pt in fx.points {
		at := '${pt.phase} ${pt.t} K, ${pt.p / 1e5} bar'
		st := thermo.solve_tp(eos, pt.t, pt.p, [1.0], root_of(pt.phase))!
		rp := thermo.residual_props(eos, st)
		fl := thermo.new_fluid(eos, [
			thermo.IdealGas{
				cp: thermo.Polynomial{
					terms: [thermo.PowerTerm{pt.cp0, 0}]
				}
			},
		])!
		pr := fl.props(st)
		dev := Deviation{
			v_rel:  (st.v - pt.v) / pt.v
			h_dep:  rp.h - pt.h_dep
			s_dep:  rp.s - pt.s_dep
			cp_rel: (pr.cp - pt.cp) / pt.cp
			w_rel:  (pr.w - pt.w) / pt.w
		}
		rec := pt.pr_deviation
		for name, pair in {
			'V':  [dev.v_rel, rec.v_rel, 1e-8]
			'Hʳ': [dev.h_dep, rec.h_dep, 1e-6]
			'Sʳ': [dev.s_dep, rec.s_dep, 1e-8]
			'Cp': [dev.cp_rel, rec.cp_rel, 1e-8]
			'w':  [dev.w_rel, rec.w_rel, 1e-8]
		} {
			got, want, tol := pair[0], pair[1], pair[2]
			assert math.abs(got - want) <= tol, '${at}: ${name} deviation ${got}, recorded ${want}'
			assert math.abs(got) <= math.abs(want) + tol, '${at}: ${name} deviation grew from ${want} to ${got}'
		}
	}
}
