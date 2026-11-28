module main

// The command line (D11, exit criterion 8): `reflux solve` writes its reports
// and exits 0, a failure exits 1 with its report, and a wrong command exits 2.
// tools/check.vsh runs the built program on the reference flowsheet too.
import flowsheet
import json2
import math
import os

const reference = os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml')
const scratch_dir = os.join_path(os.temp_dir(), 'reflux_main_test_${os.getpid()}')

fn testsuite_end() {
	os.rmdir_all(scratch_dir) or {}
}

fn scratch(name string) string {
	os.mkdir_all(scratch_dir) or { panic(err) }
	return os.join_path(scratch_dir, name)
}

// quiet runs a command, and returns its exit status and what it printed.
fn quiet(args []string) (int, string, string) {
	mut c := Console{}
	code := run(args, mut c)
	return code, c.out.str(), c.err.str()
}

fn decode(path string) !flowsheet.Report {
	return json2.decode[flowsheet.Report](os.read_file(path)!)!
}

fn test_solve_prints_the_table_and_writes_its_reports() {
	json_path, csv_path := scratch('results.json'), scratch('streams.csv')
	code, out, err := quiet(['solve', reference, '--json', json_path, '--csv=${csv_path}'])
	assert code == 0 && err == '', err
	r := decode(json_path)!
	assert r.status == 'solved' && r.file == reference && r.streams.len == 11
	assert out == r.table() + '\n'
	assert os.read_file(csv_path)! == r.csv()
}

// A flowsheet with a reactor (M9, format 2): `reflux solve` reads its
// reactions and its PFR, and the report gives the PFR's profile, from the
// feed's state to the product's. The product, hydrogen-rich at 1000 K, is
// vapour (M3's D6, as amended in M9 step 9).
fn test_solve_a_reactor() {
	file := os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reformer.toml')
	json_path := scratch('reformer.json')
	code, _, err := quiet(['solve', file, '--json', json_path])
	assert code == 0 && err == '', err
	r := decode(json_path)!
	assert r.format == 5 && r.units.len == 1 && r.units[0].kind == 'pfr'
	profile := r.units[0].profile or { panic('no profile') }
	product := r.streams.filter(it.name == 'product')[0]
	assert profile[0].volume == 0 && profile.last().volume == 1 && profile.last().length == 1
	assert profile.last().p == product.p && profile.last().molar_flows == product.molar_flows
	assert product.vapour_fraction == 1 && product.phases[0].phase == 'vapour'
}

// A spec (M8 step 4, format 3): `reflux solve` meets it, prints its value
// and writes it in the report; one out of reach exits 1, at its line.
fn test_solve_meets_a_spec() {
	spec := '
[spec.vapour]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"
'
	file := scratch('spec.toml')
	os.write_file(file, os.read_file(reference)!.replace('format = 1', 'format = 3') + spec)!
	json_path := scratch('spec.json')
	code, out, err := quiet(['solve', file, '--json', json_path, '--trace'])
	assert code == 0 && err == '', err
	assert out.contains('  Spec vapour: unit.cooler.outlet_T = 219.76 K, where gas.F = 60 mol/s, to '), out
	assert out.contains('The specs vapour, with mixer, cooler, drum, split, pump:
'), out
	r := decode(json_path)!
	assert r.specs.len == 1 && r.specs[0].vary == 'unit.cooler.outlet_T' && r.specs[0].unit == 'K'
	os.write_file(file, os.read_file(reference)!.replace('format = 1', 'format = 3') +
		spec.replace('60 mol/s', '200 mol/s'))!
	code2, _, err2 := quiet(['solve', file])
	assert code2 == 1
	assert err2.contains('spec.toml:64: spec.vapour: its target is out of reach'), err2
}

// A study (M8 step 6): `reflux study` prints its points and why an output
// has no value, writes its JSON and CSV, and exits 0 though an output had
// none; a wrong study name or option exits 2, and a file with none 1.
fn test_study_prints_its_points_and_writes_its_reports() {
	file := os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'methanol_study.toml')
	json_path, csv_path := scratch('study.json'), scratch('study.csv')
	code, out, err := quiet(['study', file, '--json', json_path, '--csv', csv_path])
	assert code == 0 && err == '', err
	assert out.contains('study temperature: 4 points, 0 solved, 4 with an output that has no value, 0 failed; in chunks of 8, 3 from a neighbour and 1 cold.'), out
	assert out.contains('methanol_conversion, at every point: it divides by syngas.n[CH4O], which is zero here'), out
	j := os.read_file(json_path)!
	assert j.contains('"format": 1') && j.contains('"value": null')
	assert os.read_file(csv_path)!.split_into_lines().len == 5
	named, _, _ := quiet(['study', file, 'temperature'])
	assert named == 0
	two, said, _ := quiet(['study', file, '--workers', '2'])
	assert two == 0 && said.contains(', on 2 workers.'), said
	for args, want in {
		'temprature':     'has no study temprature (did you mean temperature?)'
		'--trace':        "--trace is solve's"
		'a b':            'one study at a time'
		'--workerz 2':    'no option --workerz'
		'--workers 0':    '--workers needs a whole number from 1 to 1024, not "0"'
		'--workers=two':  'not "two"'
		'--workers 1025': 'not "1025"'
		'--workers 2x':   'not "2x"'
	} {
		mut full := ['study', file]
		full << args.split(' ')
		wrong, _, why := quiet(full)
		assert wrong == 2 && why.contains(want), '${args}: ${why}'
	}
	none_, _, why := quiet(['study', reference])
	assert none_ == 1 && why.contains('has no studies'), why
}

// An optimisation (M8 step 10): `reflux optimize` prints the optimum, its
// active bound and constraint and its iterations, writes its JSON and CSV,
// and exits 0; one with no feasible point exits 1, still writing its
// report; a wrong name or a second one exits 2, and a file with none 1.
fn test_optimize_prints_the_optimum_and_writes_its_reports() {
	file := os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'series.toml')
	json_path, csv_path := scratch('optimum.json'), scratch('iterations.csv')
	code, out, err := quiet(['optimize', file, 'best', '--json', json_path, '--csv', csv_path,
		'--workers', '2'])
	assert code == 0 && err == '', err
	assert out.contains('optimisation best: an optimum in '), out
	assert out.contains('maximize product.n[isopentane] = 0.54141 mol/s'), out
	assert out.contains('feed.f.T [K]   380.00  420.00  400.00   420.00  upper'), out
	assert out.contains('product.n[neopentane] <= 0.2 mol/s  0.20000  0.20000     yes'), out
	assert out.contains(', on 2 workers.'), out
	j := json2.decode[json2.Any](os.read_file(json_path)!)!.as_map()
	assert j['format']!.int() == 1 && j['status']!.str() == 'optimal' && j['goal']!.str() == 'maximize'
	iterations := j['iterations']!.as_array().len
	lines := os.read_file(csv_path)!.split_into_lines()
	assert lines.len == iterations + 1
	assert lines[0] == 'iteration,unit.tube.volume [m^3],feed.f.T [K],maximize [mol/s],violation,stationarity,step'
	none_json := scratch('none.json')
	no, _, said := quiet(['optimize', file, 'impossible', '--json', none_json])
	assert no == 1 && said == '', said
	failed := json2.decode[json2.Any](os.read_file(none_json)!)!.as_map()
	assert failed['status']!.str() == 'failed' && failed['reason']!.str().contains('no feasible point')
	for args, want in {
		'':         'has optimisations best, impossible: name one'
		'bset':     'has no optimisation bset (did you mean best?)'
		'best two': 'one optimisation at a time'
		'--trace':  "--trace is solve's"
	} {
		mut full := ['optimize', file]
		if args != '' {
			full << args.split(' ')
		}
		wrong, _, why := quiet(full)
		assert wrong == 2 && why.contains(want), '${args}: ${why}'
	}
	none_, _, why := quiet(['optimize', reference])
	assert none_ == 1 && why.contains('has no optimisations: add a [optimize.NAME]'), why
}

// A fit (M8 steps 11 and 12): `reflux fit` prints each parameter in the
// file's unit with its interval and the fitted constants, writes its JSON
// and CSV, and exits 0; a fit that fails exits 1, still writing its report;
// a wrong name or --workers exits 2, and a file with no fits 1.
fn test_fit_prints_the_constants_and_writes_its_reports() {
	file := os.join_path(@VMODROOT, 'testdata', 'm8_fit', 'ethylene.toml')
	json_path, csv_path := scratch('fit.json'), scratch('residuals.csv')
	code, out, err := quiet(['fit', file, '--json', json_path, '--csv', csv_path])
	assert code == 0 && err == '', err
	assert out.contains('fit rates: 7 parameters to 80 measurements in ethylene.csv'), out
	assert out.contains('k1.Ea      56374       2838.8          50716 to 62031            J/mol'), out
	assert out.contains('K_E.A      4.7050e-05  1.8222e-05      1.0735e-05 to 8.3366e-05  bar^-1'), out
	assert out.contains('{ A = "45480.71'), out
	j := json2.decode[json2.Any](os.read_file(json_path)!)!.as_map()
	assert j['status']!.str() == 'fitted' && j['kind']!.str() == 'rates'
	assert os.read_file(csv_path)!.split_into_lines().len == 81
	vle, said, _ := quiet(['fit', os.join_path(@VMODROOT, 'testdata', 'm8_fit', 'ethane_propane.toml')])
	assert vle == 0 && said.contains('kij[ethane, propane] = 0.0115408'), said
	// A fit whose data aren't there fails, and says why.
	broken := scratch('broken.toml')
	os.write_file(broken, os.read_file(file)!.replace('data = "ethylene.csv"', 'data = "none.csv"'))!
	failed_json := scratch('failed.json')
	no, told, _ := quiet(['fit', broken, '--json', failed_json])
	assert no == 1 && told.contains('no fit: ') && told.contains("can't be read"), told
	assert os.read_file(failed_json)!.contains('"status": "failed"')
	for args, want in {
		'rats':        'has no fit rats (did you mean rates?)'
		'--workers 2': "--workers is study's and optimize's"
		'rates other': 'one fit at a time'
	} {
		mut full := ['fit', file]
		full << args.split(' ')
		wrong, _, why := quiet(full)
		assert wrong == 2 && why.contains(want), '${args}: ${why}'
	}
	none_, _, why := quiet(['fit', os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'series.toml')])
	assert none_ == 1 && why.contains('has no fits: add a [fit.NAME]'), why
}

// The workers by default are the processors (M8's D8).
fn test_the_processors() {
	assert cpus() >= 1
	$if windows {
		assert cpus() == os.getenv('NUMBER_OF_PROCESSORS').int()
	}
}

fn test_trace_prints_every_pass_after_the_table() {
	code, out, _ := quiet(['solve', reference, '--trace'])
	assert code == 0
	lines := out.all_after('The loop of mixer, cooler, drum, split, pump, torn at recycle:\n').split_into_lines()
	assert lines[0].trim_space().starts_with('pass')
	assert lines[1..5].map(it.trim_space().all_before(' ')) == ['1', '2', '3', '4']
	assert out.index('Unit  ') or { -1 } < out.index('torn at recycle:\n') or { -1 }
}

fn test_a_file_that_does_not_load_exits_1_with_its_report() {
	// Its line on stderr, as a compiler gives one, its JSON written, and the CSV
	// of an earlier run gone.
	bad := scratch('typo.toml')
	os.write_file(bad, os.read_file(reference)!.replace('type = "heater"', 'type = "heatr"'))!
	json_path, csv_path := scratch('failed.json'), scratch('stale.csv')
	os.write_file(csv_path, 'from an earlier run')!
	code, out, err := quiet(['solve', bad, '--json', json_path, '--csv', csv_path])
	assert code == 1 && out == ''
	assert err == 'reflux: ${bad}:26: unit.cooler.type: no unit type heatr (did you mean heater?)\n'
	assert !os.exists(csv_path)
	r := decode(json_path)!
	e := r.error or { panic('no error') }
	assert r.status == 'failed' && r.file == bad && e.line == 26 && e.key == 'unit.cooler.type'

	missing := scratch('missing.toml')
	code2, _, err2 := quiet(['solve', missing])
	assert code2 == 1
	assert err2.starts_with("reflux: ${missing}: it can't be read: "), err2
}

fn test_a_flowsheet_that_does_not_solve_exits_1_with_its_passes() {
	slow := scratch('slow.toml')
	os.write_file(slow, os.read_file(reference)!.replace('method = "wegstein"', 'method = "direct"\nmax_iterations = 3'))!
	json_path := scratch('slow.json')
	code, out, err := quiet(['solve', slow, '--json', json_path])
	assert code == 1 && out == ''
	assert err.starts_with('reflux: ${slow}: tears recycle: the tears are still moving after 3 passes'), err
	assert err.ends_with('\n(--trace shows every pass)\n'), err
	r := decode(json_path)!
	assert r.status == 'failed' && r.blocks[0].residuals.len == 3
	// With --trace, the passes are on stdout and there's no hint.
	code2, out2, err2 := quiet(['solve', slow, '--trace'])
	assert code2 == 1 && !err2.contains('--trace')
	assert out2 == r.trace() + '\n'
}

fn test_a_failing_unit_is_named() {
	// A valve can't raise P.
	path := scratch('valve.toml')
	os.write_file(path, os.read_file(reference)!.replace('outlet_P = "5 bar"', 'outlet_P = "40 bar"'))!
	code, _, err := quiet(['solve', path])
	assert code == 1 && err.starts_with('reflux: ${path}: valve: '), err
}

fn test_a_report_that_cannot_be_written_exits_1() {
	// A directory where the file should go: the table is printed first.
	dir := scratch('a directory')
	os.mkdir_all(dir)!
	for flag in ['--json', '--csv'] {
		code, out, err := quiet(['solve', reference, flag, dir])
		assert code == 1 && out.contains('Stream  '), flag
		assert err.starts_with('reflux: ${dir}: '), err
	}
}

fn test_a_demo_that_fails_exits_1() {
	// No component data where it's looked for.
	old := os.getenv('REFLUX_DATA')
	os.setenv('REFLUX_DATA', scratch('no data'), true)
	code, _, err := quiet(['demo'])
	if old == '' {
		os.unsetenv('REFLUX_DATA')
	} else {
		os.setenv('REFLUX_DATA', old, true)
	}
	assert code == 1 && err.starts_with('reflux demo: '), err
}

fn test_a_wrong_command_exits_2() {
	for args, want in {
		'':                               'usage:'
		'solv':                           'reflux: no command solv (did you mean solve?)'
		'solve':                          'reflux solve: it needs a flowsheet file'
		'solve a.toml b.toml':            'reflux solve: one flowsheet at a time, not a.toml and b.toml'
		'solve a.toml --json':            'reflux solve: --json needs a file name'
		'solve a.toml --csv --trace':     'reflux solve: --csv needs a file name'
		'solve a.toml --json=':           'reflux solve: --json needs a file name'
		'solve a.toml --jsn x':           'reflux solve: no option --jsn (did you mean --json, --csv?)'
		'solve a.toml --trace=yes':       'reflux solve: no option --trace=yes (did you mean --trace?)'
		'solve a.toml --json x --json y': 'reflux solve: --json is given twice'
		'solve a.toml --csv x --csv=y':   'reflux solve: --csv is given twice'
		'solve a.toml --json x --csv x':  'reflux solve: --json and --csv name the same file, x'
		'demo abc':                       'reflux demo: argument 1 (abc): '
		'demo 300 -1':                    'reflux demo: argument 2 (-1): T and P must be above zero'
		'demo 300 nan':                   'reflux demo: argument 2 (nan): T and P must be above zero'
		'demo inf':                       'reflux demo: argument 1 (inf): T and P must be above zero'
		'demo 300 2 1':                   'reflux demo: it takes T and P, not 300 2 1'
	} {
		code, out, err := quiet(if args == '' { []string{} } else { args.split(' ') })
		assert code == 2 && out == '', args
		assert err.starts_with(want), '${args}: ${err}'
		assert err.ends_with('${usage}\n'), args
	}
	// Two names for one file that doesn't exist yet: on Windows its case is
	// the user's, and os.real_path doesn't change it.
	$if windows {
		_, _, same := quiet(['solve', 'a.toml', '--json', 'out.json', '--csv', 'OUT.json'])
		assert same.starts_with('reflux solve: --json and --csv name the same file, out.json'), same
	}
	// Nothing near enough to suggest.
	_, _, err := quiet(['solve', 'a.toml', '-q'])
	assert err == 'reflux solve: no option -q\n\n${usage}\n'
	_, _, err2 := quiet(['zzz'])
	assert err2 == 'reflux: no command zzz\n\n${usage}\n'
}

fn test_a_report_never_overwrites_the_flowsheet() {
	path := scratch('plant.toml')
	text := os.read_file(reference)!
	os.write_file(path, text)!
	mut same := [path, os.join_path(scratch_dir, '.', 'plant.toml')]
	$if windows {
		same << path.to_upper()
	}
	for flag in ['--json', '--csv'] {
		for other in same {
			code, _, err := quiet(['solve', path, flag, other])
			assert code == 2 && err.contains('would overwrite the flowsheet'), '${flag} ${other}: ${err}'
		}
	}
	assert os.read_file(path)! == text
}

fn test_the_demo_runs() {
	code, out, err := quiet(['demo', '300', '2.5'])
	assert code == 0 && err == ''
	assert out.starts_with('N-butane (106-97-8) at T = 300.0 K, P = 2.5 bar'), out
	assert out.contains("PH from the split's H")
}

fn test_help_and_version() {
	for h in ['help', '--help', '-h'] {
		code, out, err := quiet([h])
		assert code == 0 && out == '${usage}\n' && err == ''
	}
	code, out, _ := quiet(['version'])
	assert code == 0 && out == 'reflux ${flowsheet.version()}\n'
	code2, out2, _ := quiet(['--version'])
	assert code2 == 0 && out2 == out
}

const kinetics_dir = os.join_path(@VMODROOT, 'testdata', 'kinetics')

// `reflux kinetics check` prints a file's findings and exits 0; given a
// state, a fitted K is judged at its T, each typed constant says how far it
// moves the rates, and each reversible rate is held to the equilibrium the
// reactions reach from it.
fn test_kinetics_check() {
	file := os.join_path(kinetics_dir, 'xu_froment.toml')
	code, out, err := quiet(['kinetics', 'check', file])
	assert code == 0 && err == '', err
	lines := out.trim_space().split_into_lines()
	assert lines[0] == '${file}: 3 reactions, 3 rates, 2 findings'
	assert lines[1].starts_with('note: overall = reforming + shift')
	assert lines[2].starts_with('warning: K_H2O rises with T')
	state := ['900 K', '15 bar', 'CH4=0.2', 'H2O=0.45', 'H2=0.2', 'CO=0.1', 'CO2=0.05']
	mut with_state := ['kinetics', 'check', file]
	with_state << state
	code1, out1, _ := quiet(with_state)
	assert code1 == 0
	lines1 := out1.trim_space().split_into_lines()
	assert lines1[0] == '${file}: 3 reactions, 3 rates, 3 findings'
	assert lines1[3] == "note: reforming's, shift's and overall's rates are zero at the equilibrium this state reaches, over reforming and shift, and run towards it at the state"
	dwsim := os.join_path(kinetics_dir, 'dwsim_reforming.toml')
	mut dwsim_state := ['kinetics', 'check', dwsim]
	dwsim_state << state
	code2, out2, _ := quiet(dwsim_state)
	assert code2 == 0
	assert out2.contains("warning: k1 has 8.314 J/(mol*K), 5.564e-5 below R, the molar gas constant, N_A·k_B: write R unless the paper fitted with this number; at 900 K it moves reforming's rate by -0.17838% from what R gives")
	assert out2.contains("warning: overall's rate isn't zero at the equilibrium this state reaches: on its fitted K it stops where Q/K is 0.20173 (ΔG/RT -1.6008), so a reactor on it heads for a different equilibrium")
}

// `reflux kinetics explain` prints each rate and its factor table.
fn test_kinetics_explain() {
	file := os.join_path(kinetics_dir, 'vanden_bussche_froment.toml')
	code, out, err := quiet(['kinetics', 'explain', file, '513.15 K', '50 bar', 'CO2=0.08',
		'hydrogen=0.65', 'CH3OH=0.04', 'H2O=0.03', 'CO=0.08', '124-38-9=0'])
	assert code == 1, out
	assert err.contains('Carbon dioxide is given twice')
	code2, out2, err2 := quiet(['kinetics', 'explain', file, '513.15 K', '50 bar', 'CO2=0.2',
		'hydrogen=0.65', 'CH3OH=0.04', 'H2O=0.03', 'CO=0.08'])
	assert code2 == 0 && err2 == '', err2
	assert out2.starts_with("methanol: CO2 + 3 H2 = CH3OH + H2O\n  LHHW over site Cu's adsorption term cubed on partial pressures, per mass of catalyst, reversible, K fitted\n  r = k₁·p_CO₂·p_H₂·(1 − p_CH₃OH·p_H₂O/(p_CO₂·p_H₂³)/K₁)/θ_Cu³")
	assert out2.contains('\n\nRWGS: CO2 + H2 = CO + H2O\n')
}

// A file that doesn't read, or a state that's wrong, exits 1 with why; a
// wrong command exits 2.
fn test_kinetics_failures() {
	file := os.join_path(kinetics_dir, 'xu_froment.toml')
	missing := os.join_path(kinetics_dir, 'no such file.toml')
	code, _, err := quiet(['kinetics', 'check', missing])
	assert code == 1 && err.starts_with('${missing}: '), err
	bad := scratch('bad_kinetics.toml')
	os.write_file(bad, os.read_file(file)!.replace('law = "lhhw"', 'law = "lhw"'))!
	code2, _, err2 := quiet(['kinetics', 'check', bad])
	assert code2 == 1
	assert err2.trim_space() == '${bad}:23: reaction.reforming.law: no law lhw: one of power_law, lhhw, eley_rideal, expression (did you mean lhhw?)'
	for args in [
		['kinetics', 'explain', file, '900 K', '15 bar', 'CH5=1'],
		['kinetics', 'explain', file, '900 K', '15 bar', 'CH4=0.5'],
		['kinetics', 'explain', file, '900 C', '15 bar', 'CH4=1'],
		['kinetics', 'explain', file, '900 K', '15 bar', 'CH4'],
		['kinetics', 'explain', file, '900 K', '15 bar', 'CH4=x'],
	] {
		c, _, e := quiet(args)
		assert c == 1 && e != '', args.str()
	}
	_, _, e1 := quiet(['kinetics', 'explain', file, '900 K', '15 bar', 'CH4=0.5'])
	assert e1.trim_space() == 'reflux kinetics explain: the mole fractions sum to 0.5, not 1'
	for args in [
		['kinetics'],
		['kinetics', 'chek', file],
		['kinetics', 'check'],
		['kinetics', 'explain', file],
		['kinetics', 'check', file, '900 K'],
	] {
		c, _, _ := quiet(args)
		assert c == 2, args.str()
	}
	_, _, e2 := quiet(['kinetics', 'chek', file])
	assert e2.starts_with('reflux kinetics: no kinetics command chek (did you mean check?)')
	c3, out3, _ := quiet(['kinetics', 'help'])
	assert c3 == 0 && out3.starts_with('usage:\n  reflux kinetics check FILE')
}

// `reflux water` (M5a, D13, criterion 9): a state from two SPECs, by
// IAPWS-95 or IF97, as a table and as JSON, and each refusal with why.

fn spec(arg string) Spec {
	return read_spec(arg) or { panic(err) }
}

fn spec_of(g Given, value f64) Spec {
	return Spec{
		given: g
		value: value
	}
}

fn water_at(a Spec, b Spec, use_if97 bool) WaterState {
	return steam(a, b, use_if97) or { panic('${a.given}, ${b.given}: ${err.msg()}') }
}

fn fixture_rows(key string) [][]f64 {
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5_transport.json')) or { panic(err) }
	m := (json2.decode[json2.Any](text) or { panic(err) }).as_map()
	rows := m[key] or { panic('no ${key}') }
	return rows.as_array().map(it.as_array().map(it.f64()))
}

fn rel(a f64, b f64) f64 {
	return math.abs(a - b) / math.abs(b)
}

fn test_water_prints_a_state_and_writes_its_report() {
	path := scratch('water.json')
	code, out, err := quiet(['water', 'T=450degC', 'P = 40 bar', '--json', path])
	assert code == 0 && err == '', err
	st := water_at(spec('T=450degC'), spec('P = 40 bar'), false)
	assert out == st.table() + '\n'
	text := os.read_file(path)!
	assert text == st.json()
	j := json2.decode[json2.Any](text)!.as_map()
	assert j['formulation']!.str() == 'IAPWS-95' && j['phase']!.str() == 'vapour'
	assert j['T']!.f64() == 723.15 && j['P']!.f64() == 4e6 && 'x' !in j
	per_kg := j['per_kg']!.as_map()
	per_mole := j['per_mole']!.as_map()
	assert per_kg['h']!.f64() == st.whole.h && per_kg['mu']!.f64() == st.whole.mu
	assert per_mole['h']!.f64() == st.whole.h * 0.018015268 && 'mu' !in per_mole
	// Per mole: a density divided by the molar mass, everything else times it.
	assert per_mole['rho']!.f64() == st.whole.rho / 0.018015268
	assert per_mole['v']!.f64() == st.whole.v * 0.018015268
	given := j['given']!.as_array()
	assert given.len == 2 && given[1].as_map()['text']!.str() == '40 bar'
	// A wet state has both phases and the mixture, its x and σ; IF97 its
	// region, 4.
	wet := scratch('wet.json')
	code2, out2, err2 := quiet(['water', 'P=10kPa', 'x=0.9', '--if97', '--json=${wet}'])
	assert code2 == 0 && err2 == '', err2
	assert out2.contains('liquid and vapour, x = 0.9') && out2.contains('region 4')
	w := json2.decode[json2.Any](os.read_file(wet)!)!.as_map()
	assert w['x']!.f64() == 0.9 && w['region']!.int() == 4 && w['sigma']!.f64() > 0.0
	assert w['liquid']!.as_map()['cp']!.f64() > 0.0 && 'cp' !in w['per_kg']!.as_map()
	// Every key is one docs/report-format.md lists for the water report.
	doc := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md'))!.all_after('## The water report').all_before('\n## ')
	mut keys := []string{}
	for line in doc.split_into_lines().filter(it.starts_with('| `')) {
		for i, part in line.all_after('| ').all_before(' |').split('`') {
			if i % 2 == 1 {
				keys << part
			}
		}
	}
	for m in [j, w, w['per_kg']!.as_map(), w['per_mole']!.as_map(), w['liquid']!.as_map()] {
		for k in m.keys() {
			assert k in keys, 'the water report has ${k}, which report-format.md lacks'
		}
	}
	// help
	code3, out3, _ := quiet(['water', 'help'])
	assert code3 == 0 && out3.contains('T=450degC')
}

// Every single-phase state of the iapws package's grid (tools/
// gen_m5_transport.py), from its T and P: the density, viscosity and
// thermal conductivity, by IAPWS-95 and IF97. The package's T and P are
// exact, so this holds the PT flash's root and IF97's region choice.
//
// Measured: IAPWS-95's ρ within 1.9e-12, μ 1.3e-12 and λ 7.0e-12, the
// worst at 647 K and 22 MPa, where the package's density solve is loosest
// (ρ is 1/κ-sensitive there); IF97's ρ within 5.4e-13, μ 4.9e-13 and λ
// 2.9e-11 (region 5, as water/transport_test.v finds on the package's own
// ρ).
fn test_water_against_the_iapws_package() {
	for key, bounds in {
		'iapws95': [1e-11, 1e-11, 5e-11]
		'if97':    [5e-12, 5e-12, 1e-10]
	} {
		mut n := 0
		for r in fixture_rows(key) {
			st := water_at(spec_of(.t, r[0]), spec_of(.p, r[1] * 1e6), key == 'if97')
			for k, d in [rel(st.whole.rho, r[2]), rel(st.whole.mu, r[3]), rel(st.whole.lambda, r[4])] {
				assert d <= bounds[k], '${key} at ${r[0]} K and ${r[1]} MPa: ${['rho', 'mu', 'lambda'][k]} off by ${d}'
			}
			n++
		}
		assert n == if key == 'iapws95' { 791 } else { 775 }
	}
}

// Each pair that fixes a state gives back the state it came from: from a
// grid's (T, P), its (P, h), (P, s) and (T, ρ); from a wet state's (T, x),
// its (P, x), (P, h), (P, s) and (T, ρ). IF97 the same, but (T, ρ).
//
// Measured: T back within 4.3e-9 K by IAPWS-95's (P, h) and (P, s), whose
// flashes stop at 1e-10·RT of their balances, T·R/cp·1e-10 in T, under
// 2.5e-11·T for steam's cp of at least 4R; IF97's within 5.5e-12 K, its
// backward equations finished by Newton. P by (T, ρ) within 3.9% of a
// bound of ρ's rounding through the bulk modulus; a wet state's T within
// 6.4e-11 K and x within 3.9e-12, the balance's 1e-10·RT over h_vap.
fn test_water_pairs_give_back_their_state() {
	mut worst_t := [0.0, 0.0]
	mut worst_p := 0.0
	mut worst_wet := [0.0, 0.0]
	for use_if97 in [false, true] {
		for t in [300.0, 450.0, 640.0, 700.0, 1000.0] {
			for p in [1e4, 1e6, 3e7] {
				st := water_at(spec_of(.t, t), spec_of(.p, p), use_if97)
				ph := water_at(spec_of(.p, p), spec_of(.h, st.whole.h), use_if97)
				ps := water_at(spec_of(.p, p), spec_of(.s, st.whole.s), use_if97)
				for back in [ph, ps] {
					worst_t[int(use_if97)] = math.max(worst_t[int(use_if97)], math.abs(back.t - t))
					assert math.abs(back.t - t) <= 3e-11 * t, '${t} K, ${p} Pa: back at ${back.t} K'
					assert back.phase == st.phase
				}
				if !use_if97 {
					// A liquid's P moves by its bulk modulus, under 3 GPa
					// here, times ρ's rounding.
					tr := water_at(spec_of(.t, t), spec_of(.rho, st.whole.rho), false)
					worst_p = math.max(worst_p, math.abs(tr.p - p) / (1e-12 * p + 3e9 * 1e-13))
					assert math.abs(tr.p - p) <= 1e-12 * p + 3e9 * 1e-13, '${t} K, ${p} Pa: back at ${tr.p} Pa'
					assert tr.phase == st.phase
				}
			}
		}
		for t in [300.0, 450.0, 640.0] {
			for x in [0.0, 0.3, 1.0] {
				st := water_at(spec_of(.t, t), spec_of(.x, x), use_if97)
				assert st.x == x && st.sigma > 0.0
				assert st.phase == match x {
					0.0 { 'saturated liquid' }
					1.0 { 'saturated vapour' }
					else { 'liquid and vapour' }
				}
				px := water_at(spec_of(.p, st.p), spec_of(.x, x), use_if97)
				worst_t[int(use_if97)] = math.max(worst_t[int(use_if97)], math.abs(px.t - t))
				assert math.abs(px.t - t) <= 1e-9, '${t} K, x = ${x}: (P, x) at ${px.t} K'
				if x == 0.3 {
					mut backs := [
						water_at(spec_of(.p, st.p), spec_of(.h, st.whole.h), use_if97),
						water_at(spec_of(.p, st.p), spec_of(.s, st.whole.s), use_if97),
					]
					if !use_if97 {
						backs << water_at(spec_of(.t, t), spec_of(.rho, st.whole.rho), false)
					}
					for back in backs {
						worst_wet[0] = math.max(worst_wet[0], math.abs(back.t - t))
						worst_wet[1] = math.max(worst_wet[1], math.abs(back.x - x))
						assert math.abs(back.t - t) <= 1e-9 && math.abs(back.x - x) <= 1e-10, '${t} K: back at ${back.t} K, x = ${back.x}'
						l, v := back.liquid or { panic('no liquid') }, back.vapour or { panic('no vapour') }
						assert l.rho > v.rho && math.is_nan(back.whole.cp)
					}
				}
			}
		}
	}
	println('T back: IAPWS-95 ${worst_t[0]:.2e} K, IF97 ${worst_t[1]:.2e} K; P by (T, rho) ${worst_p:.2e} of its bound; wet T ${worst_wet[0]:.2e} K, x ${worst_wet[1]:.2e}')
}

// A SPEC's units: offsets, gauge pressures, molar values and percentages
// read as the file's numbers do.
fn test_water_reads_specs_in_any_unit() {
	assert spec('T=450degC').value == 723.15
	assert spec('T = 450 degC').value == 723.15
	assert spec('P=40bar').value == 4e6
	assert spec('p=1atm').value == 101325.0
	assert spec('x=50%').value == 0.5
	assert rel(spec('rho=1g/cm3').value, 1000.0) <= 1e-15
	assert math.abs(spec('h=54.0412kJ/mol').value - 54.0412e3 / 0.018015268) <= 1e-9
	assert math.abs(spec('s=6.5kJ/(kg*K)').value - 6500.0) <= 1e-12
	assert math.abs(spec('s=117.1J/(mol*K)').value - 117.1 / 0.018015268) <= 1e-9
	assert spec('rho=2e5mol/m3').value == 2e5 * 0.018015268
	assert spec('h=3e3kJ/kg').value == 3e6
}

// Each refusal says why: a wrong SPEC or pair is a wrong command (2), a
// state water or the formulation lacks is 1.
fn test_water_refusals() {
	for args, want in {
		'T=300K':                     'it needs two SPECs'
		'T=300K P=1bar h=1kJ/kg':     'it needs two SPECs, not 3'
		'T=300K T=310K':              'T is given twice'
		'T=300 P=1bar':               'it needs a unit'
		'T=300K P=1kg':               'P=1kg'
		'T=300K q=1':                 'no property q'
		'T=300K x=2':                 'x is from 0 to 1'
		'P=1bar h=3000kJ/K':          'not of kJ/kg or kJ/mol'
		'T=300K rho=-1kg/m3':         'a density is above zero'
		'T=300K P=1bar --json':       '--json needs a file name'
		'T=300K P=1bar --jsn x':      'no option --jsn'
		'T=300K h=100kJ/kg':          'T and h fix no state'
		'T=275K s=0.1kJ/(kg*K)':      'T and s can fix two states'
		'P=1bar rho=1000kg/m3':       'P and rho can fix two states'
		'h=2000kJ/kg x=0.5':          'h and x can fix two states'
		'h=2000kJ/kg s=6kJ/(kg*K)':   "doesn't solve for T and P from h and s"
		'T=500K rho=800kg/m3 --if97': 'IF97 has no state from T and rho'
	} {
		code, out, err := quiet(['water', ...args.split(' ')])
		assert code == 2 && out == '', args
		assert err.contains(want), '${args}: ${err}'
	}
	for args, want in {
		'T=240K P=1bar':         'no state below 251.165 K'
		'T=1300K P=1bar':        'IAPWS-95 holds to 1273.15 K (R6-95): not 1300 K'
		'T=500K P=1100MPa':      'IAPWS-95 holds to 1000 MPa (R6-95): not 1100 MPa'
		'P=1100MPa h=1000kJ/kg': 'IAPWS-95 holds to 1000 MPa'
		'T=650K x=0.5':          'no saturated states at or above its critical point'
		'P=23MPa x=0.5':         'no saturated states at or above its critical point'
		'T=260K P=500MPa':       'ice V at 260 K and 500 MPa'
		'T=300K P=1000MPa':      'ice VI'
		'T=500K rho=1300kg/m3':  'denser than 4 ρc'
		'T=2300K P=1bar --if97': 'no state by IF97'
		'T=270K P=1bar --if97':  'no state by IF97'
	} {
		code, out, err := quiet(['water', ...args.split(' ')])
		assert code == 1 && out == '', args
		assert err.contains(want), '${args}: ${err}'
	}
	// Ice Ih makes a liquid below 273.16 K and its melting pressure
	// metastable: the state stands, with a note.
	cold := water_at(spec('T=265K'), spec('P=1bar'), false)
	assert cold.phase == 'liquid' && cold.notes.any(it.contains('ice Ih is the stable phase'))
	assert water_at(spec('T=265K'), spec('P=200MPa'), false).notes.len == 0
	// A vapour below its sublimation pressure is stable.
	assert water_at(spec('T=265K'), spec('P=100Pa'), false).notes.len == 0
}

// R14-08's melting and sublimation pressures against its Table 3, to the
// printed digits.
fn test_ice_curves_against_r14_08() {
	assert math.abs(ih_melting_pressure(260.0) - 138.268e6) <= 500.0
	for t, want in {
		254.0: [268.685e6, 500.0]
		265.0: [479.640e6, 500.0]
		320.0: [1356.76e6, 5000.0]
	} {
		top, _ := melting_top(t) or { panic('none at ${t}') }
		assert math.abs(top - want[0]) <= want[1], '${t} K: ${top}'
	}
	assert math.abs(sublimation_pressure(230.0) - 8.94735) <= 0.000005
}

fn test_water_numbers_are_shown_to_seven_figures() {
	for x, want in {
		723.15:          '723.15'
		4e6:             '4000000'
		1.23456789e7:    '1.234568e7'
		0.00012345678:   '0.0001234568'
		1.5e-5:          '1.5e-5'
		-0.5:            '-0.5'
		9999999.6:       '1e7'
		0.99999996:      '1'
		318.95632892394: '318.9563'
	} {
		assert shown(x) == want, '${x}: ${shown(x)}'
	}
}
