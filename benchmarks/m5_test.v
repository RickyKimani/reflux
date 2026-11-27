// M5a's benchmark (step 10, exit criterion 10): what each program solved,
// that reflux is faster than DWSIM's own packages, and that the committed
// report shows it. The timings depend on the machine, so they aren't
// asserted beyond that comparison.
import benchmarks
import flash
import math
import os
import thermo
import water
import water.if97

// Each case's answers solved again under this build, against the fixture's
// (gcc -prod): h to 1e-9 kJ/kg and T to 1e-8 K, the PH flash's tolerance
// through two compilers; and the same bits before the others' runs and
// after.
fn test_reflux_solved_what_it_says() {
	f := benchmarks.load_m5()!
	fl := water.fluid(none, none)!
	m := thermo.by_roots(water.iapws95())
	for model in benchmarks.m5_models {
		for kind in benchmarks.m5_kinds {
			want := benchmarks.m5_case(f.reflux_before, model, kind)
			assert want.answers == benchmarks.m5_case(f.reflux_after, model, kind).answers
			assert want.us.len == 200, '${model} ${kind}'
			for k, s in benchmarks.m5_states {
				t, p := s[0], s[1]
				got := if model == 'IAPWS-95' {
					r := flash.pt(m, t, p, [1.0])!
					if kind == 'pt' {
						water.props95(t, 1.0 / r.phases[0].state.v).h / water.molar_mass * 1e-3
					} else {
						flash.ph(fl, p, fl.props(r.phases[0].state).h, [1.0])!.t
					}
				} else {
					st := if97.state_pt(p, t)!
					if kind == 'pt' { st.h * 1e-3 } else { if97.state_ph(p, st.h)!.t }
				}
				bound := if kind == 'pt' { 1e-9 } else { 1e-8 }
				assert math.abs(got - want.answers[k]) <= bound, '${model} ${kind} ${k}: ${got} against ${want.answers[k]}'
			}
		}
	}
}

// The others solved the same problems: every h within 3e-10 kJ/kg of
// reflux's (1.9e-10 measured, DWSIM's CoolProp at 700 K), every T back
// within each one's own iteration: IAPWS-95's 1e-8 K (5.7e-9, DWSIM's);
// IF97's DWSIM 1e-4 K and CoolProp's IF97 backend 18 mK, which leaves the
// backward equations' T unfinished.
fn test_the_others_solved_the_same_problems() {
	f := benchmarks.load_m5()!
	for others in [f.dwsim, f.coolprop, f.iapws] {
		for c in others {
			assert c.error == '' && c.us.len > 0, '${c.model} ${c.kind}: ${c.error}'
		}
	}
	assert benchmarks.m5_worst_h(f, 'IAPWS-95') <= 3e-10
	assert benchmarks.m5_worst_h(f, 'IF97') <= 1e-11
	for others in [f.dwsim, f.coolprop, f.iapws] {
		assert benchmarks.m5_agreement(f, others, 'IAPWS-95', 'ph') <= 1e-8
	}
	assert benchmarks.m5_agreement(f, f.iapws, 'IF97', 'ph') <= 1e-9
	assert benchmarks.m5_agreement(f, f.dwsim, 'IF97', 'ph') <= 2e-4
	c97 := benchmarks.m5_agreement(f, f.coolprop, 'IF97', 'ph')
	assert c97 > 1e-3 && c97 <= 0.025, '${c97}'
}

// Criterion 10: reflux is faster than DWSIM's own packages at every case.
fn test_reflux_is_faster_than_dwsim() {
	f := benchmarks.load_m5()!
	for model in benchmarks.m5_models {
		for kind in benchmarks.m5_kinds {
			r := benchmarks.m5_reflux_us(f, model, kind)
			mut d := benchmarks.m5_case(f.dwsim, model, kind).us.clone()
			d.sort()
			assert r < d[d.len / 2], '${model} ${kind}: reflux ${r} µs'
		}
	}
}

fn test_the_committed_report_has_the_m5_section() {
	f := benchmarks.load_m5()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains(benchmarks.m5_section(f))
}

// m5_worst_h takes every program's h, not only DWSIM's, which is the
// furthest in the fixture: here CoolProp's is put 0.5 kJ/kg off.
fn test_the_worst_h_is_over_every_program() {
	f := benchmarks.load_m5()!
	r := benchmarks.m5_case(f.reflux_before, 'IAPWS-95', 'pt').answers
	g := benchmarks.M5Fixture{
		...f
		coolprop: [benchmarks.M5Case{
			model:   'IAPWS-95'
			kind:    'pt'
			answers: r.map(it + 0.5)
		}]
	}
	assert math.abs(benchmarks.m5_worst_h(g, 'IAPWS-95') - 0.5) <= 1e-12
}
