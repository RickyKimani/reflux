module main

// The command line (D11, exit criterion 8): `reflux solve` writes its reports
// and exits 0, a failure exits 1 with its report, and a wrong command exits 2.
// tools/check.vsh runs the built program on the reference flowsheet too.
import flowsheet
import json2
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
// state, a fitted K is judged at its T and each typed constant says how far
// it moves the rates.
fn test_kinetics_check() {
	file := os.join_path(kinetics_dir, 'xu_froment.toml')
	code, out, err := quiet(['kinetics', 'check', file])
	assert code == 0 && err == '', err
	lines := out.trim_space().split_into_lines()
	assert lines[0] == '${file}: 3 reactions, 3 rates, 2 findings'
	assert lines[1].starts_with('note: overall = reforming + shift')
	assert lines[2].starts_with('warning: K_H2O rises with T')
	dwsim := os.join_path(kinetics_dir, 'dwsim_reforming.toml')
	code2, out2, _ := quiet(['kinetics', 'check', dwsim, '900 K', '15 bar', 'CH4=0.2', 'H2O=0.45',
		'H2=0.2', 'CO=0.1', 'CO2=0.05'])
	assert code2 == 0
	assert out2.contains("warning: k1 has 8.314 J/(mol*K), 5.564e-5 below R, the molar gas constant, N_A·k_B: write R unless the paper fitted with this number; at 900 K it moves reforming's rate by -0.17838% from what R gives")
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
