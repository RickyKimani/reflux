// M6's benchmark (step 13, D18, exit criterion 12): what each program
// solved, and that the committed report shows it. The timings depend on the
// machine, so they aren't asserted beyond being there; tools/benchmarks.vsh
// writes them to docs/benchmarks.md.
import benchmarks
import components
import flash
import math
import os
import thermo

const p_atm = 101325.0

fn fluid(db components.Database, model string) !thermo.Fluid {
	names := ['methanol', 'water']
	return match model {
		'NRTL' { db.gamma_phi(names)!.fluid }
		'UNIFAC' { db.gamma_phi(names, liquid: .unifac)!.fluid }
		else { db.fluid(.pr, names)! }
	}
}

// Each case's answers solved again under this build, against the fixture's
// (gcc -prod): to 1e-9, the flashes' tolerance through two compilers; and
// the same bits before the others' runs and after.
fn test_reflux_solved_what_it_says() {
	f := benchmarks.load_m6()!
	db := components.load_default()!
	pt := [[350.0, 0.5], [357.0, 0.3], [346.0, 0.7]]
	bubble := [0.1, 0.5, 0.9]
	for model in benchmarks.m6_models {
		fl := fluid(db, model)!
		for kind in benchmarks.m6_kinds {
			want := benchmarks.m6_case(f.reflux_before, model, kind)
			assert want.answers == benchmarks.m6_case(f.reflux_after, model, kind).answers
			assert want.us.len == 200, '${model} ${kind}'
			for k in 0 .. 3 {
				got := if kind == 'pt' {
					flash.pt(fl.model, pt[k][0], p_atm, [pt[k][1], 1.0 - pt[k][1]])!.beta
				} else {
					flash.bubble_t(fl.model, p_atm, [bubble[k], 1.0 - bubble[k]])!.t
				}
				assert math.abs(got - want.answers[k]) <= 1e-9 * math.max(want.answers[k], 1.0), '${model} ${kind} ${k}'
			}
		}
	}
}

// thermo solved the same problems, to its tolerances: 1.7e-7 in β, 7.8e-9 K
// in the bubble T, measured. DWSIM's γ–φ answers differ by its Poynting
// term's sign and its vapour pressures' extra figures (step 12): up to
// 2.3e-3 in β; its Peng–Robinson's are reflux's to 1e-4.
fn test_the_others_solved_the_same_problems() {
	f := benchmarks.load_m6()!
	tb, tt := benchmarks.m6_agreement(f, f.thermo, ['NRTL', 'UNIFAC', 'PR'])
	assert tb <= 2e-7 && tt <= 1e-8, '${tb} ${tt}'
	db, dt := benchmarks.m6_agreement(f, f.dwsim, ['NRTL', 'UNIFAC'])
	assert db > 1e-4 && db <= 5e-3 && dt <= 0.1, '${db} ${dt}'
	pb, pt := benchmarks.m6_agreement(f, f.dwsim, ['PR'])
	assert pb <= 1e-3 && pt <= 0.05, '${pb} ${pt}'
	for c in f.dwsim {
		assert c.error == '' && c.us.len == 20, '${c.model} ${c.kind}'
	}
	for c in f.thermo {
		assert c.error == '' && c.us.len == 20, '${c.model} ${c.kind}'
	}
}

fn test_the_committed_report_has_the_m6_section() {
	f := benchmarks.load_m6()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains(benchmarks.m6_section(f))
}
