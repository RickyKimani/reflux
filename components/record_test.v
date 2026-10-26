module components

fn test_set_reaches_every_kind_of_field() {
	mut r := Record{}
	r.constants.set('tc', 425.12)!
	r.parameters.set('mathias_copeman_c3', 0.84619)!
	r.correlations.set('vapor_pressure', CorrelationData{ eq: 101, coeffs: [1.0, 2.0] })!
	r.groups.set('unifac_vle', {
		1: 2
		2: 2
	})!
	assert r.constants.tc or { 0.0 } == 425.12
	assert r.constants.pc == none
	assert r.parameters.mathias_copeman_c3 or { 0.0 } == 0.84619
	vp := r.correlations.vapor_pressure or { panic('not set') }
	assert vp.eq == 101 && vp.coeffs == [1.0, 2.0]
	assert r.groups.unifac_vle == {
		1: 2
		2: 2
	}
}

fn test_unknown_keys_are_errors() {
	mut r := Record{}
	// A typo, a key from the wrong table, and an int field set as a float.
	for key in ['omgea', 'rackett', 'charge'] {
		if _ := r.constants.set(key, 1.0) {
			assert false, '`${key}` should be refused as a constant'
		}
	}
	if _ := r.parameters.set('tc', 1.0) {
		assert false, '`tc` is not a parameter'
	}
	if _ := r.correlations.set('vapour_pressure', CorrelationData{}) {
		assert false, 'British spelling is not a key'
	}
	if _ := r.groups.set('unifac', {}) {
		assert false, 'unifac is not a method key'
	}
}

fn test_groups_are_copied_not_shared() {
	mut counts := {
		1: 2
	}
	mut g := Groups{}
	g.set('asog', counts)!
	counts[1] = 99
	assert g.asog[1] == 2
}
