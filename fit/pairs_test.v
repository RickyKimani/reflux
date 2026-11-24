// An activity model's pairs fitted by reflux (M6 step 11, criterion 9):
// NRTL's two energies to a published set of methanol–water bubble points,
// and NRTL's with α and UNIQUAC's to ethanol–water's drawn from thermo's
// models (testdata/m6_fit/), against SciPy's fits of the same data with
// thermo's models (testdata/m6_fits.json, from tools/gen_m6_fits.py). Then
// where each starts, what it writes back, and what it refuses.
import components
import fit
import flowsheet
import json2
import math
import os

const dir = os.join_path(@VMODROOT, 'testdata', 'm6_fit')

fn fixture() !map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_fits.json'))!)!.as_map()
}

fn fits(file string) !flowsheet.FitFile {
	return flowsheet.load_fits_text(os.read_file(os.join_path(dir, file))!, dir: dir)!
}

fn floats(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

// Each value, standard error and interval against SciPy's: both meet their
// tests at the same minimum of the same residuals, on models that agree to
// 1e-12 (criterion 3), so the values agree to a sliver of their standard
// errors, and the errors, from the same central differences, closely.
fn test_pairs_against_scipy() {
	want := fixture()!['fits']!.as_map()
	mut worst_value, mut worst_se := 0.0, 0.0
	for file, names in {
		'methanol_water.toml':        ['nrtl']
		'ethanol_water.toml':         ['nrtl_alpha']
		'ethanol_water_uniquac.toml': ['uniquac']
	} {
		ff := fits(file)!
		for name in names {
			plan := ff.fits.filter(it.name == name)[0]
			r := fit.run(ff, plan)
			assert r.ok, '${name}: ${r.error}'
			w := want[name]!.as_map()
			values, se := floats(w['value']!), floats(w['standard_error']!)
			lower, upper := floats(w['lower']!), floats(w['upper']!)
			assert r.values.len == values.len
			assert r.dof == w['dof']!.int(), name
			assert math.abs(r.ssr - w['ssr']!.f64()) <= 1e-9 * w['ssr']!.f64(), '${name}: ${r.ssr}'
			assert math.abs(r.t - w['t']!.f64()) <= 1e-12 * r.t, name
			for q in 0 .. values.len {
				worst_value = math.max(worst_value, math.abs(r.values[q] - values[q]) / se[q])
				worst_se = math.max(worst_se, math.abs(r.standard_errors[q] - se[q]) / se[q])
				worst_value = math.max(worst_value, math.abs(r.lower[q] - lower[q]) / se[q])
				worst_value = math.max(worst_value, math.abs(r.upper[q] - upper[q]) / se[q])
			}
			println('${name}: ${r.values} ± ${r.standard_errors}, s² ${r.ssr / f64(r.dof):.4f}, ${r.iterations} iterations')
		}
	}
	println('against SciPy: values ${worst_value:.1e} and errors ${worst_se:.1e} of the errors')
	// Measured: 1.5e-6 and 4.0e-7, the two least squares' stopping tests and
	// the models' 1e-12, through J.
	assert worst_value <= 1e-5 && worst_se <= 1e-5, '${worst_value} ${worst_se}'
}

// A pair starts from the model's own: ChemSep's A₁₂ in cal/mol over
// DECHEMA's R, so A = τRT is A₁₂·R/1.98721 at any T, J/mol.
fn test_a_pair_starts_from_the_model() {
	text := os.read_file(os.join_path(dir, 'methanol_water.toml'))!.replace('"A[water, methanol]"]',
		'"A[water, methanol]", "alpha[water, methanol]"]')
	ff := flowsheet.load_fits_text(text, dir: dir)!
	ps := ff.fits[0].parameters
	r := flowsheet.molar_gas_constant()
	// nrtl.ipd: Methanol/Water, A₁₂ = −189.0469, A₂₁ = 792.8020, α = 0.2999.
	assert math.abs(ps[0].start - -189.0469 * r / 1.98721) <= 1e-12 * 400.0
	assert math.abs(ps[1].start - 792.8020 * r / 1.98721) <= 1e-12 * 3400.0
	assert ps[2].start == 0.2999 && ps[2].part == .alpha
	assert ps[0].part == .pair_energy && ps[0].unit == 'J/mol'
	// A start given in the file is used, in its unit, as J/mol.
	e := fits('ethanol_water_uniquac.toml')!.fits[0].parameters
	assert e.map(it.start) == [0.0, 0.0]
	u := flowsheet.load_fits_text(text.replace('liquid = "NRTL"', 'liquid = "UNIQUAC"').replace(', "alpha[water, methanol]"',
		''), dir: dir)!.fits[0].parameters
	// uniquac.ipd's pair, the same way: τ = exp(−A/RT).
	db := components.load_default()!
	b := db.uniquac(['67-56-1', '7732-18-5'])!
	l := b.pairs[0].listing or { panic('no listing') }
	assert math.abs(u[0].start - l.values[0] * r / 1.98721) <= 1e-12 * 1000.0
}

// A pair's start is its τ at 298.15 K as A = τRT, whatever its form: here
// a file's pair that isn't A/RT, τ = 0.5 + 300 K / T.
fn test_a_pairs_start_is_at_298_15_k() {
	text := os.read_file(os.join_path(dir, 'methanol_water.toml'))!.replace('liquid = "NRTL"
',
		'liquid = "NRTL"

[[model.pairs]]
a = "methanol"
b = "water"
tau_ab = "0.5 + 300 K / T"
tau_ba = "1"
alpha = "0.3"
')
	ps := flowsheet.load_fits_text(text, dir: dir)!.fits[0].parameters
	r := flowsheet.molar_gas_constant()
	t := 298.15
	assert math.abs(ps[0].start - (0.5 + 300.0 / t) * r * t) <= 1e-12 * 3000.0
	assert math.abs(ps[1].start - r * t) <= 1e-12 * 3000.0
	// A Henry's-law solute's pair isn't the activity model's.
	solute := os.read_file(os.join_path(dir, 'methanol_water.toml'))!.replace('names = ["methanol", "water"]',
		'names = ["methanol", "water", "hydrogen"]').replace('"A[water, methanol]"]', '"A[hydrogen, water]"]')
	if _ := flowsheet.load_fits_text(solute, dir: dir) {
		assert false, 'a solute'
	} else {
		assert err.msg().contains("hydrogen is a Henry's-law solute, outside the activity model"), err.msg()
	}
}

// The published set doesn't determine NRTL's α besides its two energies,
// but the fit still meets its test, and says so in α's interval: with
// central differences of 1e-6 it couldn't (the step's reason in fit/vle.v).
fn test_an_undetermined_alpha_still_fits() {
	text := os.read_file(os.join_path(dir, 'methanol_water.toml'))!.replace('"A[water, methanol]"]',
		'"A[water, methanol]", "alpha[methanol, water]"]')
	ff := flowsheet.load_fits_text(text, dir: dir)!
	r := fit.run(ff, ff.fits[0])
	// Where its residuals' rounding lies depends on the platform's libm: a
	// fit this badly determined meets its test on Windows, under tcc and
	// gcc, and stopped at its rounding on CI's Ubuntu (a gradient test of
	// 2.3e-10, on d8a1e0c), saying so. Either is honest; claiming an α
	// isn't.
	if !r.ok {
		assert r.error.contains('the steps fell below rounding'), r.error
		return
	}
	// Measured on Windows: α 0.199 ± 5.29, its interval −10.7 to 11.1.
	assert r.standard_errors[2] > 1.0 && r.lower[2] < 0.0 && r.upper[2] > 1.0, '${r.values} ${r.standard_errors}'
}

// What a fit writes back is the file's own [[model.pairs]], and the file
// with it loads to the model at the answer: its pair from the file, and
// its bubble points the fit's.
fn test_a_fitted_pair_writes_back() {
	for file in ['methanol_water.toml', 'ethanol_water.toml', 'ethanol_water_uniquac.toml'] {
		ff := fits(file)!
		r := fit.run(ff, ff.fits[0])
		assert r.ok, '${file}: ${r.error}'
		rep := r.report(file, ff)
		names := ff.names
		assert rep.section.starts_with('[[model.pairs]]\na = "${names[0]}"\nb = "${names[1]}"\ntau_ab = "'), rep.section
		if file == 'methanol_water.toml' {
			// α wasn't fitted, so it's the bank's, as it was.
			assert rep.section.ends_with('alpha = "0.2999"\n'), rep.section
		}
		text := os.read_file(os.join_path(dir, file))!
		liquid := if file.contains('uniquac') { 'UNIQUAC' } else { 'NRTL' }
		back := text.replace('[model]\nliquid = "${liquid}"\n', '[model]\nliquid = "${liquid}"\n\n' +
			rep.section).all_before('start = {') + text.all_after('start = {').all_after('}\n')
		again := flowsheet.load_fits_text(back, dir: dir)!
		g := again.gamma_phi or { panic('no γ–φ') }
		assert g.pairs[0].source == .file
		// It starts where the fit ended, to the written values' 10 figures,
		// the pair written back being the file's, and comes back there.
		for q, p in again.fits[0].parameters {
			assert math.abs(p.start - r.values[q]) <= 1e-9 * math.max(math.abs(r.values[q]), 1.0), '${file} ${p.name}'
		}
		r2 := fit.run(again, again.fits[0])
		assert r2.ok
		for q in 0 .. r.values.len {
			assert math.abs(r2.values[q] - r.values[q]) <= 1e-6 * r.standard_errors[q], '${r2.values} ${r.values}'
		}
	}
}

// What a fit of pairs refuses, at the parameters' line.
fn test_what_a_fit_of_pairs_refuses() {
	text := os.read_file(os.join_path(dir, 'methanol_water.toml'))!
	for bad, want in {
		'["A[methanol, water]", "kij[methanol, water]"]':       "kij is an equation of state's"
		'["alpha[methanol, water]", "alpha[water, methanol]"]': 'listed twice, as alpha[methanol, water]'
		'["A[methanol, methanol]"]':                            'a pair needs two components'
		'["A[methanol, ethanol]"]':                             'ethanol is not one of the components'
		'["A[methanol water]"]':                                'a pair is written A[a, b]'
	} {
		changed := text.replace('parameters = ["A[methanol, water]", "A[water, methanol]"]\n', 'parameters = ${bad}\n')
		if _ := flowsheet.load_fits_text(changed, dir: dir) {
			assert false, bad
		} else {
			assert err.msg().contains(want), err.msg()
		}
	}
	ethanol := os.read_file(os.join_path(dir, 'ethanol_water.toml'))!
	for liquid, want in {
		'UNIFAC':  "NRTL's or UNIQUAC's pair energy, and the model is UNIFAC"
		'UNIQUAC': "NRTL's α, and the model is UNIQUAC"
	} {
		changed := ethanol.replace('liquid = "NRTL"', 'liquid = "${liquid}"')
		if _ := flowsheet.load_fits_text(changed, dir: dir) {
			assert false, liquid
		} else {
			assert err.msg().contains(want), err.msg()
		}
	}
	eos := text.replace('format = 4', 'format = 3').replace('liquid = "NRTL"', 'equation = "PR"')
	if _ := flowsheet.load_fits_text(eos, dir: dir) {
		assert false, 'an equation of state'
	} else {
		assert err.msg().contains('needs a [model] with liquid = "NRTL" or "UNIQUAC"'), err.msg()
	}
}
