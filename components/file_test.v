module components

// The component-file writer and reader, on small hand-made records. The
// round trip of the whole ChemSep database is in chemsep/generate_test.v.

fn minimal() Record {
	return Record{
		name:   'Test'
		source: 'test'
		ids:    Ids{
			cas: '106-97-8'
		}
	}
}

const minimal_text = 'format = 1\nname = "Test"\nsource = "test"\n\n[ids]\ncas = "106-97-8"\n'

fn round_trip(r Record) !Record {
	return parse_toml(write_toml(r)!, '${r.ids.cas}.toml')!
}

fn test_a_minimal_record_and_nothing_more() {
	// Empty tables and absent values aren't written at all.
	assert write_toml(minimal())! == minimal_text
	assert parse_toml(minimal_text, '106-97-8.toml')! == minimal()
}

fn test_the_text_depends_only_on_the_record() {
	mut a := minimal()
	a.constants.set('tc', 425.12)!
	a.constants.set('pc', 3.796e6)!
	mut b := minimal()
	b.constants.set('pc', 3.796e6)!
	b.constants.set('tc', 425.12)!
	assert write_toml(a)! == write_toml(b)!
	// Keys come in the struct's order, not the order they were set in.
	assert write_toml(a)!.contains('tc = 425.12\npc = 3.796e+06\n')
}

fn test_floats_always_read_back_as_the_same_floats() {
	mut r := minimal()
	values := [2.0, 0.0, -0.0, 1e21, 1e-300, 0.1 + 0.2, 0.0581222, -125790.0, 5e-324]
	for i, v in values {
		r.correlations.set('vapor_pressure', CorrelationData{
			eq:     100
			coeffs: [v]
			tmin:   1.0
			tmax:   2.0
			units:  'Pa'
		})!
		text := write_toml(r)!
		back := round_trip(r)!
		got := back.correlations.vapor_pressure or { panic('lost') }
		assert got.coeffs[0] == v, 'value ${i}: ${v} came back as ${got.coeffs[0]}'
		// A whole number is written 2.0, never 2, or TOML would read an integer.
		if v == 2.0 {
			assert text.contains('coeffs = [2.0]')
		}
	}
}

fn test_values_toml_cannot_carry_are_refused() {
	for bad in [f64_nan(), f64_inf()] {
		mut r := minimal()
		r.constants.set('tc', bad)!
		if _ := write_toml(r) {
			assert false, '${bad} was written'
		}
	}
	mut r := minimal()
	r.name = 'tab\there'
	if _ := write_toml(r) {
		assert false, 'a control character was written'
	}
}

fn f64_nan() f64 {
	return 0.0 / zero()
}

fn f64_inf() f64 {
	return 1.0 / zero()
}

fn zero() f64 {
	return 0.0
}

fn test_strings_are_escaped() {
	mut r := minimal()
	// SMILES use backslashes for stereo bonds.
	r.ids.smiles = 'F/C=C\\F'
	r.name = 'say "hi"'
	r.ids.aliases = ['a', 'b\\c']
	text := write_toml(r)!
	assert text.contains('smiles = "F/C=C\\\\F"')
	assert text.contains('name = "say \\"hi\\""')
	assert round_trip(r)! == r
}

fn test_sources_override_the_file_source() {
	mut r := minimal()
	r.constants.set('tc', 425.125)!
	r.sources['constants.tc'] = 'CoolProp 8.0.0'
	r.correlations.set('ideal_gas_cp', CorrelationData{
		eq:     16
		coeffs: [1.0]
		tmin:   100.0
		tmax:   200.0
		units:  'J/kmol/K'
		tc:     425.0
	})!
	r.sources['correlations.ideal_gas_cp'] = 'a fit of my own'
	text := write_toml(r)!
	assert text.contains('tc = { value = 425.125, source = "CoolProp 8.0.0" }')
	assert text.contains('tc = 425.0\nsource = "a fit of my own"')
	assert round_trip(r)! == r
	// A source for a value that isn't there is a mistake.
	r.sources['constants.pc'] = 'nowhere'
	if _ := write_toml(r) {
		assert false, 'a source without a value was written'
	}
}

fn test_every_kind_of_field_round_trips() {
	mut r := minimal()
	r.ids = Ids{
		cas:           '7732-18-5'
		inchikey:      'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
		pubchem:       962
		smiles:        'O'
		structure:     'H2O'
		chemsep_index: 1
		family:        12
		aliases:       ['water', 'H2O']
	}
	r.constants.set('tc', 647.096)!
	r.constants.charge = -1
	r.parameters.set('mathias_copeman_c2', -0.1)!
	r.correlations.set('liquid_density', CorrelationData{
		eq:     105
		coeffs: [5.459, 0.30542, 647.13, 0.081]
		tmin:   273.16
		tmax:   333.15
		units:  'kmol/m3'
	})!
	r.groups.set('unifac_vle', {
		7: 1
	})!
	r.groups.set('asog', {
		19: 3
		2:  1
	})!
	assert round_trip(r)! == r
}

// The reader, on text that's wrong in one way at a time.

fn refused(text string, message string) {
	if r := parse_toml(text, '106-97-8.toml') {
		assert false, 'read ${r.name}; expected `${message}`'
	} else {
		assert err.msg().contains(message), err.msg()
	}
}

fn test_unknown_keys_are_refused() {
	refused(minimal_text + 'family = 1\n[extra]\nx = 1\n', 'unknown key `extra`')
	refused(minimal_text + '\n[constants]\nomgea = 0.2\n', 'unknown key `constants.omgea`')
	refused(minimal_text + '\n[parameters]\ntc = 1.0\n', 'unknown key `parameters.tc`')
	refused(minimal_text + 'colour = "blue"\n', 'unknown key `ids.colour`')
	refused(minimal_text +
		'\n[correlations.vapour_pressure]\neq = 101\ncoeffs = [1.0]\nt_range = [1.0, 2.0]\nunits = "Pa"\n',
		'unknown correlation `vapour_pressure`')
	refused(minimal_text +
		'\n[correlations.vapor_pressure]\neq = 101\ncoeffs = [1.0]\nt_range = [1.0, 2.0]\nunits = "Pa"\nform = 3\n',
		'unknown key `correlations.vapor_pressure.form`')
	refused(minimal_text + '\n[groups.unifac]\n"1" = 2\n', 'unknown group method `unifac`')
	refused(minimal_text + '\n[constants]\ntc = { value = 1.0, source = "x", note = "y" }\n',
		'unknown key `note`')
}

fn test_wrong_types_are_refused() {
	// V's toml would read the string as 0.0 without complaint.
	refused(minimal_text + '\n[constants]\ntc = "425.12"\n',
		'`constants.tc` should be a number, not a string')
	refused(minimal_text + 'pubchem = 7843.0\n', '`ids.pubchem` should be an integer, not a float')
	refused(minimal_text + '\n[constants]\ncharge = 0.5\n', 'should be an integer')
	refused(minimal_text + '\n[constants]\ncharge = { value = 0, source = "x" }\n',
		'should be an integer, not a table')
	refused(minimal_text + 'aliases = ["a", 1]\n', 'should be a string, not an integer')
	refused(minimal_text.replace('name = "Test"', 'name = 7'), '`name` should be a string')
	refused(minimal_text + '\n[constants]\ntc = { value = 1.0 }\n', 'no `source`')
}

fn test_numbers_are_checked() {
	// An integer is fine where a float is expected…
	r := parse_toml(minimal_text + '\n[constants]\ntc = 425\n', '106-97-8.toml')!
	assert r.constants.tc or { 0.0 } == 425.0
	// …unless it's too big to be exact. V's toml reads `inf` as such an
	// integer, and `nan` as a string; both are refused.
	refused(minimal_text + '\n[constants]\ntc = 9007199254740993\n', 'too large to be exact')
	refused(minimal_text + '\n[constants]\ntc = inf\n', 'too large to be exact')
	refused(minimal_text + '\n[constants]\ntc = nan\n', 'should be a number, not a string')
}

fn test_required_keys() {
	refused(minimal_text.replace('format = 1\n', ''), 'no `format` key')
	refused(minimal_text.replace('format = 1', 'format = 2'), 'format 2; this reflux reads format 1')
	refused(minimal_text.replace('name = "Test"\n', ''), 'no `name`')
	refused(minimal_text.replace('source = "test"\n', ''), 'no `source`')
	refused('format = 1\nname = "Test"\nsource = "test"\n', 'no [ids] table')
	refused(minimal_text.replace('cas = "106-97-8"', 'inchikey = "X"'), 'no `ids.cas`')
	corr := '\n[correlations.vapor_pressure]\neq = 101\ncoeffs = [1.0]\nt_range = [1.0, 2.0]\nunits = "Pa"\n'
	parse_toml(minimal_text + corr, '106-97-8.toml')!
	for key in ['eq = 101\n', 'coeffs = [1.0]\n', 't_range = [1.0, 2.0]\n', 'units = "Pa"\n'] {
		refused(minimal_text + corr.replace(key, ''), 'no `correlations.vapor_pressure.${key.all_before(' ')}`')
	}
	refused(minimal_text + corr.replace('[1.0, 2.0]', '[1.0, 2.0, 3.0]'), 'should be [Tmin, Tmax]')
}

fn test_the_file_must_be_named_after_its_cas_number() {
	parse_toml(minimal_text, 'some/dir/106-97-8.toml')!
	if _ := parse_toml(minimal_text, 'butane.toml') {
		assert false
	} else {
		assert err.msg() == 'butane.toml: holds 106-97-8, so it should be named 106-97-8.toml'
	}
}

fn test_groups_are_checked() {
	parse_toml(minimal_text + '\n[groups.asog]\n"19" = 3\n', '106-97-8.toml')!
	refused(minimal_text + '\n[groups.asog]\n"019" = 3\n', '`019` is not a group number')
	refused(minimal_text + '\n[groups.asog]\nx = 3\n', '`x` is not a group number')
	refused(minimal_text + '\n[groups.asog]\n"19" = 0\n', 'a count must be positive')
}

fn test_toml_errors_name_the_file() {
	refused(minimal_text + 'pubchem = 1\npubchem = 2\n', '106-97-8.toml: ')
	refused('format = = 1', '106-97-8.toml: ')
}

fn test_the_index_round_trips_sorted() {
	entries := [
		IndexEntry{
			cas:     '7732-18-5'
			name:    'Water'
			aliases: ['H2O']
		},
		IndexEntry{
			cas:      '106-97-8'
			name:     'N-butane'
			inchikey: 'IJDNQMDRQITEOD-UHFFFAOYSA-N'
		},
	]
	text := write_index(entries)!
	assert (text.index('106-97-8') or { -1 }) < (text.index('7732-18-5') or { -1 })
	back := parse_index(text)!
	assert back == [entries[1], entries[0]]
	// Each refusal for its own reason. (A key appended at the end of the
	// text would land in the last [compounds."…"] table, so the top-level
	// one goes straight after `format`.)
	bad := {
		text.replace('format = 1', 'format = 2'):             'format 2; this reflux reads format 1'
		text.replace('format = 1', 'format = 1\nextra = 1'):  'unknown key `extra`'
		text.replace('name = "Water"', 'title = "Water"'):    'unknown key `compounds."7732-18-5".title`'
		text.replace('name = "Water"\n', ''):                 'no `compounds."7732-18-5".name`'
		text.replace('aliases = ["H2O"]', 'aliases = "H2O"'): 'should be an array'
	}
	for t, message in bad {
		if _ := parse_index(t) {
			assert false, 'read an index that should fail with `${message}`'
		} else {
			assert err.msg().contains(message), err.msg()
		}
	}
}
