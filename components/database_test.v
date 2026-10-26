module components

import math
import os

// The built-in data, loaded once. Loading reads only index.toml; each find
// parses one compound's file.
const builtin = load_default() or { panic(err) }

fn lookup_error(d Database, q string) LookupError {
	r := d.find(q) or {
		// `err is LookupError` checks the error's concrete type and, in the
		// branch, lets you use it as one: Go's errors.As, built in.
		if err is LookupError {
			return err
		}
		panic(err)
	}
	panic('`${q}` found ${r.name}')
}

fn test_the_built_in_data_loads() {
	assert builtin.count() == 642
	assert builtin.dirs == [default_dir]
}

fn test_find_by_cas_name_or_inchikey() {
	for q in ['106-97-8', 'N-butane', 'n-butane', '  N-BUTANE  ', 'IJDNQMDRQITEOD-UHFFFAOYSA-N',
		'ijdnqmdrqiteod-uhfffaoysa-n'] {
		r := builtin.find(q)!
		assert r.ids.cas == '106-97-8', q
	}
}

fn test_an_unknown_name_suggests_the_closest() {
	// Names containing the query come first, shortest first.
	e := lookup_error(builtin, 'butane')
	assert !e.ambiguous
	assert e.suggestions.len == 5
	assert e.suggestions[..2] == ['N-butane (106-97-8)', 'Isobutane (75-28-5)']
	assert e.msg().starts_with('no compound `butane` in the database. Did you mean N-butane (106-97-8), Isobutane (75-28-5),')
	// Otherwise the smallest edit distance wins.
	assert lookup_error(builtin, 'watr').suggestions[0] == 'Water (7732-18-5)'
	assert lookup_error(builtin, 'tolune').suggestions[0] == 'Toluene (108-88-3)'
}

fn test_an_ambiguous_query_lists_every_match() {
	// ChemSep gives carbon monoxide formaldehyde's InChIKey (CO's own is
	// UGFAIRIUMAVXCW-UHFFFAOYSA-N), so this key names two compounds.
	e := lookup_error(builtin, 'WSFSSNUMVMOOMR-UHFFFAOYSA-N')
	assert e.ambiguous
	assert e.suggestions == ['Formaldehyde (50-00-0)', 'Carbon monoxide (630-08-0)']
	assert e.msg().contains('ask for one by its CAS number')
	// Its CAS number still finds each one.
	assert builtin.find('630-08-0')!.name == 'Carbon monoxide'
}

// A directory of the user's own, removed again by the caller.
fn user_dir(name string, files map[string]string) string {
	dir := os.join_path(os.temp_dir(), 'reflux_db_test_${os.getpid()}_${name}')
	os.rmdir_all(dir) or {}
	os.mkdir_all(dir) or { panic(err) }
	for file, text in files {
		os.write_file(os.join_path(dir, file), text) or { panic(err) }
	}
	return dir
}

// 12345-67-4 is a valid CAS number that ChemSep doesn't have.
const made_up = 'format = 1\nname = "Unobtainium"\nsource = "a test"\n\n[ids]\ncas = "12345-67-4"\n\n[constants]\ntc = 500.0\n'

fn test_a_users_directory_adds_and_overrides() {
	butane := os.read_file(os.join_path(default_dir, '106-97-8.toml'))!
	assert butane.contains('\ntc = 425.12\n')
	dir := user_dir('override', {
		'106-97-8.toml':   butane.replace('\ntc = 425.12\n', '\ntc = 425.125\n')
		'12345-67-4.toml': made_up
	})
	defer {
		os.rmdir_all(dir) or {}
	}
	if _ := builtin.find('Unobtainium') {
		assert false, 'the test compound is in ChemSep after all'
	}
	mine := load([default_dir, dir])!
	assert mine.count() == 643
	// The whole record is replaced: a later directory wins.
	assert mine.find('n-butane')!.constants.tc or { 0.0 } == 425.125
	assert mine.find('unobtainium')!.constants.tc or { 0.0 } == 500.0
	// The built-in database is untouched.
	assert builtin.find('n-butane')!.constants.tc or { 0.0 } == 425.12
}

fn test_a_name_used_twice_is_ambiguous() {
	dir := user_dir('clash', {
		'12345-67-4.toml': made_up.replace('Unobtainium', 'Water')
	})
	defer {
		os.rmdir_all(dir) or {}
	}
	e := lookup_error(load([default_dir, dir])!, 'water')
	assert e.ambiguous
	assert e.suggestions == ['Water (12345-67-4)', 'Water (7732-18-5)']
}

fn test_bad_directories_are_errors() {
	bad := user_dir('bad', {
		'12345-67-4.toml': made_up.replace('tc = 500.0', 'tc = "hot"')
	})
	defer {
		os.rmdir_all(bad) or {}
	}
	if _ := load([bad]) {
		assert false, 'a bad file loaded'
	} else {
		assert err.msg().contains('12345-67-4.toml: `constants.tc` should be a number')
	}
	// An index that lists a file that isn't there is stale.
	stale := user_dir('stale', {
		'index.toml': 'format = 1\n\n[compounds."106-97-8"]\nname = "N-butane"\n'
	})
	defer {
		os.rmdir_all(stale) or {}
	}
	if _ := load([stale]) {
		assert false, 'a stale index loaded'
	} else {
		assert err.msg().contains('lists 106-97-8, but') && err.msg().contains('is missing')
	}
	if _ := load([os.join_path(stale, 'nowhere')]) {
		assert false, 'a missing directory loaded'
	}
}

fn test_n_butane_as_thermo_sees_it() {
	r := builtin.find('n-butane')!
	c := r.component()!
	assert c.name == 'N-butane' && c.tc == 425.12 && c.pc == 3.796e6 && c.omega == 0.199
	assert c.mw == 0.0581222
	ig := r.ideal_gas()!
	assert ig.dhf or { 0.0 } == -125790.0
	assert ig.s0 or { 0.0 } == 309.91
	// ChemSep's eq. 16 in J/(kmol·K), converted: Cp = A + exp(B/T + C + DT + ET²).
	cp := fn (t f64) f64 {
		return 1e-3 * (44749.9 + math.exp(-338.141 / t + 11.8145 + 0.00097744 * t - 3.35913e-07 * t * t))
	}
	for t in [200.0, 300.0, 1000.0] {
		assert math.tolerance(ig.cp.cp(t), cp(t), 1e-14), 'Cp at ${t}'
	}
	// H(400) − H(300) by quadrature. The reference is M2 step 2's figure.
	assert math.tolerance(ig.cp.int_cp(300.0, 400.0), 11414.705, 1e-7)
}

fn test_a_compound_without_critical_constants() {
	r := builtin.find('3967-55-3')!
	if _ := r.component() {
		assert false
	} else {
		assert err.msg() == '1,2-Dichloroethylene Carbonate (3967-55-3) has no Tc, Pc, ω, so no equation of state can use it'
	}
}

fn test_properties_come_in_si_per_mol() {
	r := builtin.find('n-butane')!
	// Vapour pressure, eq. 101, already in Pa.
	vp := r.property('vapor_pressure')!
	assert vp.key == 'vapor_pressure' && vp.note == '' && vp.unit == 'Pa'
	want := math.exp(68.5773 - 4444.92 / 350.0 - 7.39584 * math.log(350.0) + 9.85743e-06 * 350.0 * 350.0)
	v := vp.at(350.0)
	assert math.tolerance(v.value, want, 1e-14) && v.in_range
	// Outside the fitted range [134.86, 425.16] K it extrapolates, and says so.
	assert !vp.at(100.0).in_range && vp.at(100.0).value > 0.0
	// Liquid density, eq. 105 in kmol/m³, becomes mol/m³.
	rho := r.property('liquid_density')!
	assert rho.unit == 'mol/m³'
	assert math.tolerance(rho.at(300.0).value, 1e3 * 1.0023 / math.pow(0.26457, 1.0 +
		math.pow(1.0 - 300.0 / 425.17, 0.27138)), 1e-14)
	// The heat of vaporization, eq. 106, takes Tc from the constants: zero above it.
	hv := r.property('heat_of_vaporization')!
	assert hv.unit == 'J/mol'
	assert hv.at(300.0).value > 1e4 && hv.at(425.13).value == 0.0
}

fn test_an_unusable_vapour_pressure_falls_back_to_antoine() {
	// T-cinnamaldehyde's vapour pressure is ChemSep's eq. 200.
	r := builtin.find('14371-10-9')!
	vp := r.property('vapor_pressure')!
	assert vp.key == 'antoine_vapor_pressure'
	assert vp.note == "the Antoine vapour pressure, because T-Cinnamaldehyde (14371-10-9) vapor_pressure is equation 200, which reflux can't evaluate"
	assert math.tolerance(vp.at(500.0).value, math.exp(20.9798 - 3757.48 / (500.0 - 92.5645)),
		1e-14)
	// With neither usable, the error gives both reasons.
	mut both := r
	both.correlations.antoine_vapor_pressure = none
	if _ := both.property('vapor_pressure') {
		assert false
	} else {
		assert err.msg().contains('equation 200') && err.msg().contains('no `antoine_vapor_pressure` data')
	}
}

fn test_defective_correlations_are_refused() {
	refused := fn (cas string, key string, message string) {
		r := builtin.find(cas) or { panic(err) }
		if p := r.property(key) {
			assert false, '${cas} ${key} gave ${p.key}'
		} else {
			assert err.msg().contains(message), err.msg()
		}
	}
	refused('67-66-3', 'solid_cp', 'Chloroform (67-66-3) solid_cp has an impossible range [283.0, 0.0] K (see data/components/CHANGES.md)')
	refused('1333-74-0p', 'rpp_cp', 'rpp_cp has no coefficients')
	refused('106-97-8', 'liquid_viscosity_rps', 'N-butane (106-97-8) has no `liquid_viscosity_rps` data')
	refused('106-97-8', 'vapour_pressure', 'no correlation called `vapour_pressure`')
	refused('106-97-8', 'relative_permittivity', 'is equation 120')
}

fn test_every_compound_builds() {
	// What a flowsheet needs from each compound: an ideal gas for all 642,
	// a component for all but the one that lacks Tc, Pc and ω, and a vapour
	// pressure, directly or by the Antoine fallback, for every one.
	// The counts were checked independently, with Python's tomllib over the
	// same files.
	mut built, mut fallbacks := 0, 0
	mut no_ideal_gas, mut no_psat := []string{}, []string{}
	for cas, _ in builtin.entries {
		r := builtin.find(cas)!
		if ig := r.ideal_gas() {
			cp := ig.cp.cp(298.15)
			assert math.is_finite(cp) && cp > 0.0, '${cas} Cp(298.15 K) = ${cp}'
		} else {
			no_ideal_gas << cas
		}
		if _ := r.component() {
			built++
		}
		if vp := r.property('vapor_pressure') {
			if vp.key != 'vapor_pressure' {
				fallbacks++
			}
		} else {
			no_psat << cas
		}
	}
	no_ideal_gas.sort()
	// The stub record, and R-1234yf, which ChemSep gives no ideal-gas Cp.
	assert no_ideal_gas == ['3967-55-3', '754-12-1']
	assert built == 641
	assert fallbacks == 26
	assert no_psat == ['3967-55-3']
}
