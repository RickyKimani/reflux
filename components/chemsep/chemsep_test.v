module chemsep

import components

// The vendored database, imported once for all the tests below. A const
// with a function call is computed when the test program starts.
const vendored = import_files(@VMODROOT, ['data/vendor/chemsep/chemsep1.xml',
	'data/vendor/chemsep/chemsep2.xml']) or { panic(err) }

fn find(cas string) components.Record {
	for r in vendored.records {
		if r.ids.cas == cas {
			return r
		}
	}
	panic('${cas} not imported')
}

// The vendored files: what M2 exit criterion 1 asks for.

fn test_every_record_imports() {
	assert vendored.sources.len == 2
	s1, s2 := vendored.sources[0], vendored.sources[1]
	assert s1.path == 'data/vendor/chemsep/chemsep1.xml' && s1.version == 'v8.49' && s1.records == 433
	assert s2.path == 'data/vendor/chemsep/chemsep2.xml' && s2.version == 'v8.50' && s2.records == 211
	// 644 records, two of them second copies.
	assert vendored.records.len == 642
	mut seen := map[string]bool{}
	for r in vendored.records {
		assert !seen[r.ids.cas], 'CAS ${r.ids.cas} twice'
		seen[r.ids.cas] = true
		assert r.name != ''
		assert r.source in ['ChemSep v8.49 (data/vendor/chemsep/chemsep1.xml)',
			'ChemSep v8.50 (data/vendor/chemsep/chemsep2.xml)']
	}
}

fn test_duplicates_keep_the_first_copy() {
	d := vendored.duplicates
	assert d.len == 2
	assert d[0].cas == '75-46-7' && d[0].name == 'Trifluoromethane'
	assert d[0].kept == 'chemsep2.xml record 179' && d[0].dropped == 'chemsep2.xml record 182'
	assert d[1].cas == '107-05-1' && d[1].differing == ['LiquidViscosity']
	// Only the first copies have a liquid viscosity, and the first
	// trifluoromethane has the refitted vapour pressure.
	assert find('107-05-1').correlations.liquid_viscosity != none
	vp := find('75-46-7').correlations.vapor_pressure or { panic('no Psat') }
	assert vp.coeffs[0] == 67.1034
}

fn test_repairs_are_the_known_defects() {
	// The files have 43 stray `<`s, but two are in the dropped duplicate
	// copies, and repairs are reported only for records kept.
	smiles := vendored.repairs.filter(it.key == 'ids.smiles')
	assert smiles.len == 41
	for r in smiles {
		assert r.before == r.after + '<'
		assert find(r.cas).ids.smiles or { '' } == r.after
	}
	assert smiles[0].cas == '3178-22-1' && smiles[0].after == 'CC(C)(C)C1CCCCC1'
	sums := vendored.repairs.filter(it.key == 'groups.asog')
	assert sums.len == 6
	assert vendored.repairs.len == 47
	// Diglyme's three ether oxygens: ASOG group 19, listed as 2 and 1.
	assert find('111-96-6').groups.asog[19] == 3
	assert sums.any(it.cas == '111-96-6' && it.before == 'group 19: 2 + 1'
		&& it.after == 'group 19: 3')
	// No SMILES anywhere still has markup in it.
	for r in vendored.records {
		assert !(r.ids.smiles or { '' }).contains_any('<>&')
	}
}

fn test_findings() {
	count := fn (k FindingKind) int {
		return vendored.findings.filter(it.kind == k).len
	}
	// Forms reflux can't evaluate, kept verbatim: 322 + 5 permittivities,
	// 5 vapour pressures (eq. 200), 1 liquid viscosity (118), 1 liquid Cp (124).
	assert count(.unsupported_form) == 334
	assert vendored.findings.filter(it.kind == .unsupported_form
		&& it.what == 'vapor_pressure: equation 200').len == 5
	// Inverted or empty ranges, e.g. chloroform's solid Cp over [283, 0] K.
	assert count(.bad_range) == 42
	assert vendored.findings.any(it.kind == .bad_range && it.name == 'Chloroform'
		&& it.what == 'solid_cp: range [283.0, 0.0] K')
	// Four eq.-2 solid densities with more than A and B, four eq.-4 amine
	// liquid Cps with an E.
	assert count(.unused_coefficients) == 8
	assert vendored.findings.filter(it.kind == .unused_coefficients
		&& it.what.starts_with('liquid_cp: equation 4')).len == 4
	assert count(.no_coefficients) == 1
	assert vendored.findings.any(it.kind == .no_coefficients && it.cas == '1333-74-0p')
	assert count(.partial_set) == 5
	assert count(.incomplete) == 1
	assert vendored.findings.any(it.kind == .incomplete && it.cas == '3967-55-3'
		&& it.what == 'no tc, pc, omega')
}

fn test_n_butane_against_hand_conversions() {
	r := find('106-97-8')
	assert r.name == 'N-butane'
	assert r.source == 'ChemSep v8.49 (data/vendor/chemsep/chemsep1.xml)'
	assert r.ids.smiles or { '' } == 'CCCC'
	assert r.ids.inchikey or { '' } == 'IJDNQMDRQITEOD-UHFFFAOYSA-N'
	assert r.ids.structure or { '' } == 'CH3(CH2)2CH3'
	assert r.ids.pubchem or { 0 } == 7843
	assert r.ids.chemsep_index or { 0 } == 5
	assert r.ids.family or { 0 } == 1
	c := r.constants
	// Values not per kmol come through unchanged…
	assert c.tc or { 0.0 } == 425.12
	assert c.pc or { 0.0 } == 3.796e6
	assert c.omega or { 0.0 } == 0.199
	assert c.zc or { 0.0 } == 0.274
	assert c.dipole or { 1.0 } == 0.0
	assert c.charge or { 1 } == 0
	// …and per-kmol ones are divided by 1000 exactly: 58.1222 kg/kmol is the
	// double nearest 0.0581222 kg/mol, and −1.2579E+08 J/kmol is −125790 J/mol.
	assert c.mw or { 0.0 } == 0.0581222
	assert c.vc or { 0.0 } == 0.000255
	assert c.v_liquid_tb or { 0.0 } == 0.00010048
	assert c.dhf or { 0.0 } == -125790.0
	assert c.dgf or { 0.0 } == -16700.0
	assert c.s0 or { 0.0 } == 309.91
	assert c.dh_fusion or { 0.0 } == 4661.0
	assert c.dh_combustion or { 0.0 } == -2657320.0
	p := r.parameters
	assert p.rackett or { 0.0 } == 0.2728
	assert p.costald_volume or { 0.0 } == 0.00025439
	assert p.vdw_area or { 0.0 } == 6.94e5
	assert p.parachor or { 0.0 } == 3.39e-5
	assert p.mathias_copeman_c1 or { 0.0 } == 0.84209
	assert p.mathias_copeman_c2 or { 0.0 } == -0.46406
	assert p.mathias_copeman_c3 or { 0.0 } == 0.84619
	assert p.uniquac_q_prime or { 0.0 } == 2.776
	// Correlations keep ChemSep's coefficients and units, untouched.
	vp := r.correlations.vapor_pressure or { panic('no Psat') }
	assert vp.eq == 101
	assert vp.coeffs == [68.5773, -4444.92, -7.39584, 9.85743e-06, 2.0]
	assert vp.tmin == 134.86 && vp.tmax == 425.16 && vp.units == 'Pa'
	assert vp.tc == none
	cp := r.correlations.ideal_gas_cp or { panic('no Cp') }
	assert cp.eq == 16 && cp.units == 'J/kmol/K' && cp.coeffs.len == 5
	visc := r.correlations.vapor_viscosity or { panic('no viscosity') }
	assert visc.eq == 102 && visc.coeffs == [2.7078e-08, 0.97147, -51.16, 6431.0]
	assert r.groups.unifac_vle == {
		1: 2
		2: 2
	}
	assert r.groups.gc_method.len == 2
}

fn test_a_chemsep2_compound() {
	// Cyrene, from the newer file.
	r := find('53716-82-8')
	assert r.source == 'ChemSep v8.50 (data/vendor/chemsep/chemsep2.xml)'
	assert r.name == 'Cyrene'
}

fn test_per_kmol_values_are_exact_decimal_shifts() {
	assert shift_decimal('58.1222', -3)! == '58.1222e-3'
	assert shift_decimal('3.796E+06', -3)! == '3.796e3'
	assert shift_decimal('-1.2579E+08', -3)! == '-1.2579e5'
	assert shift_decimal('2.886E-10', -3)! == '2.886e-13'
	// ChemSep's Vc of 0.0417 m³/kmol. Dividing by 1000 would give
	// 4.1700000000000004e-05; 220 of its 2637 distinct per-kmol values
	// would come out one unit in the last place away like that.
	assert number('0.0417', true)! == 4.17e-05
	assert 0.0417 / 1000.0 != 4.17e-05
	assert number('58.1222', true)! == 0.0581222
	for bad in ['', 'abc', '1e', '1.2.3', 'nan', 'inf', '1e+x'] {
		if v := number(bad, true) {
			assert false, '`${bad}` gave ${v}'
		}
	}
}

fn test_cas_numbers_are_checked() {
	for good in ['106-97-8', '7732-18-5', '1333-74-0p', '1333-74-0o', '50-00-0'] {
		check_cas(good)!
	}
	// A wrong check digit, wrong shapes, a capital suffix, and a path.
	for bad in ['106-97-7', '106-978', '1-97-8', '12345678-97-8', '106-97-8P', '106-97-8pp', '../97-8',
		'106-9a-8', ''] {
		if _ := check_cas(bad) {
			assert false, '`${bad}` passed'
		}
	}
}

fn test_every_table_key_is_a_record_field() {
	// A typo in fields.v would otherwise surface only when some file used
	// that element.
	mut r := components.Record{}
	for tag, f in fields {
		match f.kind {
			.constant { r.constants.set(f.key, 1.0) or { assert false, '${tag}: ${err}' } }
			.parameter { r.parameters.set(f.key, 1.0) or { assert false, '${tag}: ${err}' } }
			.correlation {
				r.correlations.set(f.key, components.CorrelationData{}) or {
					assert false, '${tag}: ${err}'
				}
			}
			.groups { r.groups.set(f.key, {}) or { assert false, '${tag}: ${err}' } }
			.text_id {
				assert f.key in ['cas', 'smiles', 'inchikey', 'structure'], tag
			}
			.int_id {
				assert f.key in ['chemsep_index', 'family', 'pubchem'], tag
			}
			.name, .int_constant {}
		}
		// Scalars are converted exactly when their unit is per kmol.
		// Correlations never are: their coefficients stay verbatim (D5).
		if f.kind in [.constant, .parameter] {
			assert f.per_kmol == f.units.contains('/kmol'), tag
		} else {
			assert !f.per_kmol, tag
		}
	}
}

// Small synthetic files, for the checks the vendored data never trips.

const minimal = '<CompoundID name="Name" value="Test" />\n<CAS name="CAS number" value="106-97-8" />'

fn file_with(body string) string {
	return '<compounds library="ChemSep v9.99 test">\n<compound>\n${body}\n</compound>\n</compounds>\n'
}

fn import_one(body string) !Import {
	mut imp := Import{}
	imp.add_file(file_with(body), 'dir/t.xml')!
	return imp
}

fn refused(body string, message string) {
	if imp := import_one(body) {
		assert false, 'imported ${imp.records.len} record(s), expected `${message}`'
	} else {
		assert err.msg().contains(message), err.msg()
	}
}

fn test_a_minimal_record() {
	imp := import_one(minimal)!
	assert imp.records.len == 1
	r := imp.records[0]
	assert r.name == 'Test' && r.ids.cas == '106-97-8'
	assert r.source == 'ChemSep v9.99 (dir/t.xml)'
	assert r.constants.tc == none && r.correlations.vapor_pressure == none
	assert imp.findings.len == 1 && imp.findings[0].kind == .incomplete
	assert imp.sources[0].path == 'dir/t.xml' && imp.sources[0].records == 1
}

fn test_the_smiles_repair() {
	imp := import_one(minimal + '\n<SMILES name="SMILES" value="CC&lt;" />')!
	// Written with a literal `<` as ChemSep does, the same result:
	imp2 := import_one(minimal + '\n<SMILES name="SMILES" value="CC<" />')!
	for i in [imp, imp2] {
		assert i.records[0].ids.smiles or { '' } == 'CC'
		assert i.repairs.len == 1
		assert i.repairs[0].cas == '106-97-8' && i.repairs[0].name == 'Test'
		assert i.repairs[0].before == 'CC<' && i.repairs[0].after == 'CC'
	}
	// Anywhere else, the same defect is refused.
	refused('<CompoundID name="Name" value="Test<" />\n<CAS name="CAS number" value="106-97-8" />',
		'contains markup')
	refused(minimal + '\n<SMILES name="SMILES" value="C<C" />', 'Malformed XML')
}

fn test_unknown_or_mislabelled_elements_are_refused() {
	refused(minimal + '\n<Colour name="Colour" value="blue" />', 'unknown element <Colour>')
	refused(minimal + '\n<CriticalVolume name="Vc" units="m3/mol" value="0.000255" />',
		'expected `m3/kmol`')
	refused(minimal + '\n<CriticalTemperature name="Tc" value="425.12" />', 'expected `K`')
	refused(minimal + '\n<CriticalTemperature name="Tc" units="K" value="425.12" extra="1" />',
		'unknown attribute `extra`')
	refused(minimal + '\n<CriticalTemperature name="Tc" units="K" value="hot" />',
		'`hot` is not a number')
	refused(minimal + '\n<CriticalTemperature name="Tc" units="K" value="nan" />',
		'not a finite number')
	refused(minimal + '\n<Family name="Family" value="1.5" />', 'not an integer')
	refused(minimal + '\n<CAS name="CAS number" value="106-97-8" />', '<CAS> appears twice')
	refused('<CompoundID name="Name" value="Test" />', 'no CAS number')
	refused('<CAS name="CAS number" value="106-97-8" />', 'has no name')
	refused('<CompoundID name="Name" value="Test" />\n<CAS name="CAS number" value="106-97-7" />',
		'fails its check digit')
	refused(minimal + '\nloose text', 'contains text')
}

fn correlation_body(children string) string {
	return minimal +
		'\n<VaporPressure name="Vapour pressure" units="Pa">\n${children}\n</VaporPressure>'
}

fn test_correlations_are_read_strictly() {
	ok := '<eqno value="101" />\n<A value="1" />\n<B value="2" />\n<Tmin units="K" value="100" />\n<Tmax units="K" value="200" />'
	imp := import_one(correlation_body(ok))!
	vp := imp.records[0].correlations.vapor_pressure or { panic('not read') }
	assert vp.eq == 101 && vp.coeffs == [1.0, 2.0] && vp.tmin == 100.0 && vp.tmax == 200.0
	assert vp.units == 'Pa'
	// A gap in the coefficients would shift their meaning.
	refused(correlation_body(ok.replace('<B value="2" />', '<C value="2" />')),
		'<C> where <B> should be')
	refused(correlation_body(ok.replace('<A value="1" />', '<A value="1" />\n<A value="1" />')),
		'<A> twice')
	refused(correlation_body(ok.replace('<Tmax units="K" value="200" />', '')), 'no <Tmax>')
	refused(correlation_body(ok.replace('<eqno value="101" />', '')), 'no <eqno>')
	refused(correlation_body(ok.replace('units="K" value="100"', 'units="C" value="100"')),
		'not in K')
	refused(correlation_body(ok + '\n<F value="1" />'), 'unknown child <F>')
	refused(correlation_body(ok).replace('units="Pa">', 'units="Pa" eq="101">'),
		'unknown attribute `eq`')
	// An inverted range is the data's problem: kept, and reported.
	inverted := import_one(correlation_body(ok.replace('value="200"', 'value="5"')))!
	assert inverted.findings.any(it.kind == .bad_range
		&& it.what == 'vapor_pressure: range [100.0, 5.0] K')
	// So is a form tdep can't evaluate.
	eq200 := import_one(correlation_body(ok.replace('"101"', '"200"')))!
	assert eq200.records[0].correlations.vapor_pressure != none
	assert eq200.findings.any(it.kind == .unsupported_form)
}

fn test_groups() {
	body := minimal +
		'\n<Asog name="ASOG">\n<group id="2" value="3" />\n<group id="19" value="2" />\n<group id="19" value="1" />\n</Asog>'
	imp := import_one(body)!
	assert imp.records[0].groups.asog == {
		2:  3
		19: 3
	}
	assert imp.repairs.len == 1 && imp.repairs[0].key == 'groups.asog'
	refused(body.replace('id="2"', 'id="two"'), '`two` is not an integer')
	refused(body.replace('<group id="2"', '<item id="2"'), 'unknown child <item>')
}

fn test_duplicates_across_files() {
	mut imp := Import{}
	imp.add_file(file_with(minimal + '\n<CriticalTemperature name="Tc" units="K" value="100" />'),
		'a.xml')!
	imp.add_file(file_with(minimal + '\n<CriticalTemperature name="Tc" units="K" value="200" />\n<SMILES name="SMILES" value="CC<" />'),
		'b.xml')!
	assert imp.records.len == 1
	assert imp.records[0].constants.tc or { 0.0 } == 100.0
	assert imp.duplicates.len == 1
	d := imp.duplicates[0]
	assert d.kept == 'a.xml record 1' && d.dropped == 'b.xml record 1'
	assert d.differing == ['CriticalTemperature', 'SMILES']
	// The dropped copy's repair isn't reported: its record isn't kept.
	assert imp.repairs.len == 0
}

fn test_files_that_are_not_chemsep() {
	mut imp := Import{}
	if _ := imp.add_file('<things library="ChemSep v8.49"></things>', 'x.xml') {
		assert false
	} else {
		assert err.msg().contains('not <compounds>')
	}
	if _ := imp.add_file('<compounds library="Somebody else"></compounds>', 'x.xml') {
		assert false
	} else {
		assert err.msg().contains('no ChemSep version')
	}
	// Errors say where.
	two := '<compounds library="ChemSep v1.0 x">\n<compound>\n${minimal}\n</compound>\n<compound>\n<CompoundID name="Name" value="B" />\n</compound>\n</compounds>'
	mut imp2 := Import{}
	if _ := imp2.add_file(two, 'y.xml') {
		assert false
	} else {
		assert err.msg() == 'y.xml record 2: no CAS number'
	}
}
