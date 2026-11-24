// Fits (M8 steps 11 and 12, D14, D15, criterion 9): rate constants and a
// kij fitted by reflux against SciPy's fits of the same data
// (tools/gen_m8_fits.py), parameters and intervals to the optimum's
// tolerance.
module fit

import flowsheet
import json2
import math
import numeric
import os
import quantity
import thermo

fn fixture() map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm8_fits.json')) or {
		panic(err)
	}) or { panic(err) }.as_map()
}

fn floats(v json2.Any) []f64 {
	return v.as_array().map(it.f64())
}

fn loaded(name string) !flowsheet.FitFile {
	return flowsheet.load_fits(os.join_path(@VMODROOT, 'testdata', 'm8_fit', name))
}

// The rate fit: seven constants of an LHHW rate and a power law against
// SciPy's, each value, standard error and interval bound. SciPy stops at
// its own tolerance, 1e-15 on the cost's fall, and reflux at the gradient
// test's 1e-10; the cost is flat near its minimum, so each value agrees to
// about 1e-6 of its standard error, which the test holds at 1e-4.
fn test_rate_constants_against_scipy() {
	want := fixture()['rates']!.as_map()
	f := loaded('ethylene.toml')!
	r := run(f, f.fits[0])
	assert r.ok, r.error
	names := want['names']!.as_array().map(it.str())
	assert f.fits[0].parameters.map(it.name) == names
	// SciPy's numbers are in the file's units, bar; reflux's in SI.
	scales := f.fits[0].parameters.map(scale_of(it))
	values := times(floats(want['value']!), scales)
	ses := times(floats(want['standard_error']!), scales)
	for q, name in names {
		assert math.abs(r.values[q] - values[q]) <= 1e-4 * ses[q], '${name}: ${r.values[q]}, SciPy ${values[q]} ± ${ses[q]}'
		assert math.abs(r.standard_errors[q] - ses[q]) <= 1e-6 * ses[q], '${name}: se ${r.standard_errors[q]}, SciPy ${ses[q]}'
	}
	lower := times(floats(want['lower']!), scales)
	upper := times(floats(want['upper']!), scales)
	for q, name in names {
		assert math.abs(r.lower[q] - lower[q]) <= 1e-4 * ses[q], '${name}: lower'
		assert math.abs(r.upper[q] - upper[q]) <= 1e-4 * ses[q], '${name}: upper'
	}
	corr := want['correlation']!.as_array()
	for p in 0 .. names.len {
		row := floats(corr[p])
		for q in 0 .. names.len {
			assert math.abs(r.correlation[p][q] - row[q]) <= 1e-6, '${names[p]}, ${names[q]}'
		}
	}
	assert math.abs(r.ssr - want['ssr']!.f64()) <= 1e-9 * want['ssr']!.f64()
	assert r.dof == want['dof']!.int() && math.abs(r.t - want['t']!.f64()) <= 1e-12
	assert r.centred.len == 3 && math.abs(r.centred[0].t_ref - want['t_ref']!.f64()) <= 1e-12
	println('rates: ${r.iterations} iterations, ${r.evaluations} evaluations, gradient test ${r.gradient}')
	// Measured: 10 iterations, which can't grow.
	assert r.iterations <= 10
}

// The kij fit: one kij against SciPy's, on thermo's PR with reflux's
// constants. Both PRs agree to about 1e-11 in a bubble point (M3), and
// SciPy's Jacobian is by differences too, so the values agree to about
// 1e-6 of the standard error.
fn test_kij_against_scipy() {
	want := fixture()['vle']!.as_map()
	f := loaded('ethane_propane.toml')!
	r := run(f, f.fits[0])
	assert r.ok, r.error
	v, se := want['value']!.f64(), want['standard_error']!.f64()
	assert math.abs(r.values[0] - v) <= 1e-3 * se, 'kij ${r.values[0]}, SciPy ${v} ± ${se}'
	assert math.abs(r.standard_errors[0] - se) <= 1e-4 * se, 'se ${r.standard_errors[0]}, SciPy ${se}'
	assert math.abs(r.lower[0] - want['lower']!.f64()) <= 1e-3 * se
	assert math.abs(r.upper[0] - want['upper']!.f64()) <= 1e-3 * se
	assert math.abs(r.ssr - want['ssr']!.f64()) <= 1e-6 * want['ssr']!.f64()
	println('kij: ${r.values[0]} ± ${r.standard_errors[0]}, ${r.iterations} iterations, ${r.evaluations} evaluations')
	// Measured: 4 iterations under tcc and gcc, which can't grow.
	assert r.iterations <= 4
}

// Student's t against mpmath's, at four probabilities and degrees of
// freedom from 1 to 10⁶. Measured: 1.2e-14 of t at ν ≤ 1000, 3.2e-12 at
// 10⁶, where the incomplete beta's continued fraction is longest; the test
// holds 5e-14 and 1e-11.
fn test_student_t() {
	for qv in fixture()['t_quantiles']!.as_array() {
		q := qv.as_map()
		p, nu, want := q['p']!.f64(), q['nu']!.f64(), q['t']!.f64()
		got := numeric_t(p, nu)
		tol := if nu <= 1000 { 5e-14 } else { 1e-11 }
		assert math.abs(got - want) <= tol * math.abs(want), 't(${p}, ${nu}) = ${got}, mpmath ${want}'
	}
	assert numeric_t(0.5, 3.0) == 0.0
	assert numeric_t(0.025, 10.0) == -numeric_t(0.975, 10.0)
}

// scale_of is a parameter's unit's SI scale, 1 for a pure number.
fn scale_of(p flowsheet.FitParameter) f64 {
	if p.unit == '' {
		return 1.0
	}
	return (quantity.unit(p.unit) or { panic(err) }).scale
}

fn numeric_t(p f64, nu f64) f64 {
	return numeric.student_t_quantile(p, nu) or { panic(err) }
}

// times is each value times its scale.
fn times(v []f64, s []f64) []f64 {
	mut out := []f64{len: v.len}
	for i in 0 .. v.len {
		out[i] = v[i] * s[i]
	}
	return out
}

// fit_documented is each object's keys in the fit report's section of
// docs/report-format.md, by its path.
fn fit_documented() map[string][]string {
	text := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md')) or { panic(err) }
	section := text.all_after('## The fit report').all_before('\n## ')
	mut out := map[string][]string{}
	mut path := ''
	for line in section.split_into_lines() {
		if line == '#### The fit' {
			path = ''
		} else if line.starts_with('#### `') {
			path = line.all_after('`').all_before('`')
		} else if line.starts_with('| `') {
			out[path] << line.all_after('| `').all_before('`')
		}
	}
	return out
}

// The JSON's keys are the format's, in order, and nothing more, for a fit
// and for one that failed.
fn test_the_json_is_as_documented() {
	doc := fit_documented()
	f := loaded('ethylene.toml')!
	good := run(f, f.fits[0])
	j := json2.decode[json2.Any](good.report('x', f).json())!.as_map()
	assert j.keys() == doc[''], '${j.keys()}'
	assert j['parameters']!.as_array()[0].as_map().keys() == doc['parameters[]']
	assert j['residuals']!.as_array()[0].as_map().keys() == doc['residuals[]']
	assert j['centred']!.as_array()[0].as_map().keys() == ['constant', 't_ref']
	bad := failed(f.fits[0], 'no data')
	k := json2.decode[json2.Any](bad.report('x', f).json())!.as_map()
	assert k.keys() == doc['']
	assert k['parameters']!.as_array()[0].as_map()['value']! is json2.Null
	assert doc.keys().len == 3
}

// The same rates in mol/(kg·h), 3600 times the numbers, give the same fit:
// each column is read in its own unit.
fn test_a_rate_in_another_unit() {
	src := os.join_path(@VMODROOT, 'testdata', 'm8_fit')
	dir := os.join_path(os.temp_dir(), 'reflux_fit_units_${os.getpid()}')
	os.mkdir_all(dir)!
	defer {
		os.rmdir_all(dir) or {}
	}
	lines := os.read_file(os.join_path(src, 'ethylene.csv'))!.split_into_lines()
	mut out := [lines[0].replace('[mol/(kg*s)]', '[mol/(kg*h)]')]
	for line in lines[1..] {
		mut cells := line.split(',')
		for k in [6, 7] {
			cells[k] = (cells[k].f64() * 3600.0).str()
		}
		out << cells.join(',')
	}
	os.write_file(os.join_path(dir, 'ethylene.csv'), out.join('\n') + '\n')!
	text := os.read_file(os.join_path(src, 'ethylene.toml'))!
	hourly := flowsheet.load_fits_text(text, dir: dir)!
	per_second := loaded('ethylene.toml')!
	a := run(hourly, hourly.fits[0])
	b := run(per_second, per_second.fits[0])
	assert a.ok && b.ok, a.error
	for q in 0 .. a.values.len {
		assert math.abs(a.values[q] - b.values[q]) <= 1e-6 * b.standard_errors[q], '${q}'
	}
}

// A kij starts from the model's own, the database's where the file takes
// it, and 0 where it doesn't.
fn test_a_kij_starts_from_the_model() {
	vle := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm8_fit', 'ethane_propane.toml'))!
	off := flowsheet.load_fits_text(vle, dir: os.join_path(@VMODROOT, 'testdata', 'm8_fit'))!
	assert off.fits[0].parameters[0].start == 0.0
	on := flowsheet.load_fits_text(vle.replace('kij = false', 'kij = true'),
		dir: os.join_path(@VMODROOT,
			'testdata', 'm8_fit')
	)!
	fl := on.fluid or { panic('no fluid') }
	c := fl.helmholtz()! as thermo.Cubic
	assert c.kij[1] != 0.0 && on.fits[0].parameters[0].start == c.kij[1]
}
