// M5b's benchmark (step 11, exit criterion 10): what reflux solved, that
// it's faster than DWSIM's own functions at every transport case, and that
// the committed report shows it. The timings depend on the machine, so
// they aren't asserted beyond that comparison.
import benchmarks
import components
import flash
import flowsheet
import math
import os
import thermo
import water.if97

const db = components.load_default() or { panic(err) }

struct BenchPhase {
	names []string
	z     []f64
	t     f64
	p     f64
	vapor bool
}

const phases = {
	'benzene-toluene liquid':    BenchPhase{['benzene', 'toluene'], [0.4, 0.6], 300.0, 101325.0, false}
	'reforming gas':             BenchPhase{['methane', 'hydrogen', 'water', 'carbon dioxide',
		'carbon monoxide'], [
		0.15,
		0.5,
		0.2,
		0.05,
		0.1,
	], 1000.0, 3e6, true}
	'methane-propane dense gas': BenchPhase{['methane', 'propane'], [0.8, 0.2], 320.0, 1e7, true}
}

const water_states = [[300.0, 3e6], [500.0, 1e6], [700.0, 3e7]]

fn close(got f64, want f64, rel f64, at string) {
	assert math.abs(got - want) <= rel * math.abs(want), '${at}: ${got} against ${want}'
}

// Each case's answers solved again under this build, against the fixture's
// (gcc -prod): a phase's values to 1e-13 of themselves, each a few
// correlations and a mixing rule; a PT flash's volume to 1e-13; a PH
// flash's T to 1e-8 K, its tolerance through two compilers; IF97's h to
// 1e-12 and T to 1e-10 K; a report's phase counts exactly; and the same
// answers before the others' runs and after.
fn test_reflux_solved_what_it_says() {
	f := benchmarks.load_m5b()!
	for name in benchmarks.m5b_transport_cases {
		want := benchmarks.m5b_case(f.reflux_before.transport, name, 'phase')
		assert want.answers == benchmarks.m5b_case(f.reflux_after.transport, name, 'phase').answers
		assert want.us.len == 2000, name
		c := phases[name] or { panic(name) }
		fl := db.setup(.pr, c.names)!.fluid
		mix := components.mixture(c.names.map(db.find(it)!))!
		root := if c.vapor { thermo.Root.vapor } else { thermo.Root.liquid }
		label := if c.vapor { thermo.PhaseLabel.vapor } else { thermo.PhaseLabel.liquid }
		st := fl.model.phase(c.t, c.p, c.z, root)!
		tr := mix.phase(label, c.t, c.p, c.z, st.v)
		for k, r in [tr.viscosity, tr.conductivity, tr.surface_tension] {
			got := if r.ok() { r.value } else { -1.0 }
			close(got, want.answers[k], 1e-13, '${name} ${k}')
		}
	}
	s95 := db.iapws95('water')!
	s97 := db.if97('water')!
	for name in benchmarks.m5b_water {
		for kind in benchmarks.m5b_kinds {
			want := benchmarks.m5b_case(f.reflux_before.water, name, kind)
			assert want.answers == benchmarks.m5b_case(f.reflux_after.water, name, kind).answers
			assert want.us.len == 200
			for k, s in water_states {
				t, p := s[0], s[1]
				at := '${name} ${kind} ${k}'
				if name == 'IF97 functions' {
					st := if97.state_pt(p, t)!
					if kind == 'pt' {
						close(st.h, want.answers[k], 1e-12, at)
					} else {
						assert math.abs(if97.state_ph(p, st.h)!.t - want.answers[k]) <= 1e-10, at
					}
					continue
				}
				fl := if name == 'IF97 flash' { s97.fluid } else { s95.fluid }
				r := flash.pt(fl.model, t, p, [1.0])!
				if kind == 'pt' {
					close(r.phases[0].state.v, want.answers[k], 1e-13, at)
				} else {
					got := flash.ph(fl, p, fl.props(r.phases[0].state).h, [1.0])!.t
					assert math.abs(got - want.answers[k]) <= 1e-8, at
				}
			}
		}
	}
	for name in benchmarks.m5b_reports {
		l := flowsheet.load(os.join_path(@VMODROOT, 'testdata', 'flowsheets', '${name}.toml'),
			db: db
		)!
		rep := l.report(l.solve()!)
		mut with, mut without := 0.0, 0.0
		for s in rep.streams {
			for ph in s.phases {
				if _ := ph.viscosity {
					with++
				} else {
					without++
				}
			}
		}
		want := benchmarks.m5b_case(f.reflux_before.report, name, 'with transport').answers
		assert [with, without] == want, name
	}
}

// The others solved their cases: no error, and DWSIM's Steam Tables' h
// IF97's own (kJ/kg) to its flash's rounding.
fn test_the_others_solved_their_cases() {
	f := benchmarks.load_m5b()!
	for others in [f.dwsim, f.thermo, f.dwsim_steam_tables] {
		for c in others {
			assert c.error == '' && c.us.len > 0, '${c.name} ${c.kind}: ${c.error}'
		}
	}
	h := benchmarks.m5b_case(f.dwsim_steam_tables, 'Steam Tables', 'pt').answers
	mine := benchmarks.m5b_case(f.reflux_before.water, 'IF97 functions', 'pt').answers
	for k, v in h {
		assert math.abs(v * 1e3 - mine[k]) <= 1e-6 * math.abs(mine[k]), '${k}: ${v} kJ/kg against ${mine[k]} J/kg'
	}
}

// Criterion 10: reflux's phase is faster than DWSIM's functions at every
// case.
fn test_reflux_is_faster_than_dwsim() {
	f := benchmarks.load_m5b()!
	for name in benchmarks.m5b_transport_cases {
		r := benchmarks.m5b_reflux_us(f.reflux_before.transport, f.reflux_after.transport,
			name, 'phase')
		mut d := benchmarks.m5b_case(f.dwsim, name, 'phase').us.clone()
		d.sort()
		assert r < d[d.len / 2], '${name}: reflux ${r} µs'
	}
}

fn test_the_committed_report_has_the_m5b_section() {
	f := benchmarks.load_m5b()!
	text := os.read_file(benchmarks.report_path)!
	assert text.contains(benchmarks.m5b_section(f))
}
