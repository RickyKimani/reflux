module kinetics

import chem

fn species(names ...string) []Species {
	return names.map(Species{
		name:    it
		formula: chem.parse(it) or { panic(err) }
	})
}

fn reactions(sp []Species, eqs ...string) []Reaction {
	mut out := []Reaction{}
	for i, e in eqs {
		name := ['reforming', 'shift', 'overall', 'r4', 'r5'][i]
		out << parse_equation(name, e, sp) or { panic(err) }
	}
	return out
}

fn described(sp []Species, eqs ...string) string {
	rs := reactions(sp, ...eqs)
	d := dependence(rs, sp.len) or { panic(err) }
	return '${d.rank}: ' + d.combinations.map(it.describe(rs)).join('; ')
}

// The owner's set isn't independent: the overall reaction is steam
// reforming plus the shift, whichever order they come in.
fn test_the_owners_set() {
	sp := species('CH4', 'H2O', 'CO', 'H2', 'CO2')
	assert described(sp, 'CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2', 'CH4 + 2 H2O = CO2 + 4 H2') == '2: overall = reforming + shift'
	assert described(sp, 'CH4 + 2 H2O = CO2 + 4 H2', 'CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2') == '2: overall = reforming − shift'
	assert described(sp, 'CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2') == '2: '
	// Combinations name only reactions before them, with their factors.
	rs := reactions(sp, 'CH4 + 2 H2O = CO2 + 4 H2', 'CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2')
	d := dependence(rs, sp.len)!
	assert d.combinations.len == 1
	c := d.combinations[0]
	assert c.reaction == 2
	assert c.of.len == 2
	assert c.of[0].reaction == 0 && c.of[0].factor == one
	assert c.of[1].reaction == 1 && c.of[1].factor == one.neg()
}

fn test_factors() {
	sp := species('CO', 'O2', 'CO2', 'H2', 'H2O', 'CH4O')
	assert described(sp, 'CO + 0.5 O2 = CO2', '2 CO + O2 = 2 CO2') == '1: shift = 2 × reforming'
	assert described(sp, 'CO + 0.5 O2 = CO2', 'CO2 = CO + 1/2 O2') == '1: shift = − reforming'
	assert described(sp, '3 CO + 1.5 O2 = 3 CO2', 'CO + 0.5 O2 = CO2') == '1: shift = 1/3 × reforming'
	// Methanol synthesis: from CO is from CO₂ plus the shift.
	assert described(sp, 'CO2 + 3 H2 = CH4O + H2O', 'CO + H2O = CO2 + H2', 'CO + 2 H2 = CH4O') == '2: overall = reforming + shift'
	// The same reaction twice.
	assert described(sp, 'CO + H2O = CO2 + H2', 'CO + H2O = CO2 + H2') == '1: shift = reforming'
	// Independent, with species no reaction shares.
	assert described(sp, 'CO + 0.5 O2 = CO2', 'H2 + 0.5 O2 = H2O', 'CO + 2 H2 = CH4O') == '3: '
	// Three reactions among four species with three elements: the third is
	// a combination of the first two.
	assert described(sp, 'CO + 0.5 O2 = CO2', 'H2 + 0.5 O2 = H2O', 'CO2 + H2 = CO + H2O') == '2: overall = − reforming + shift'
}

fn test_charge_balances() {
	sp := [
		Species{
			name:    'H2O'
			formula: chem.parse('H2O')!
		},
		Species{
			name:    'H+'
			formula: chem.parse('H')!
			charge:  1
		},
		Species{
			name:    'OH-'
			formula: chem.parse('HO')!
			charge:  -1
		},
		Species{
			name:    'OH'
			formula: chem.parse('HO')!
		},
	]
	r := parse_equation('ionisation', "H2O = 'H+' + 'OH-'", sp)!
	assert r.nu == [-1.0, 1, 1]
	if _ := parse_equation('radical', "H2O = 'H+' + OH", sp) {
		assert false
	} else {
		assert err.msg() == "radical: H2O = 'H+' + OH doesn't balance: charge: 0 on the left, 1 on the right; no coefficients would balance it, as charge is only on the right"
		e := err as ReactionError
		assert e.imbalances.len == 1 && e.imbalances[0].element == 'charge'
	}
	// A formula two species share, at different charges, is refused as
	// ambiguous: name one.
	if _ := parse_equation('x', 'H2O = H + HO', sp) {
		assert false
	} else {
		assert err.msg() == 'x: column 11: HO is the formula of both OH- and OH: name one'
	}
}

fn test_find() {
	sp := [
		Species{
			name:    'Methane'
			cas:     '74-82-8'
			aliases: ['marsh gas']
			formula: chem.parse('CH4')!
		},
		Species{
			name:    'Ethanol'
			cas:     '64-17-5'
			aliases: ['alcohol']
			formula: chem.parse('C2H6O')!
		},
		Species{
			name:    'Dimethyl ether'
			cas:     '115-10-6'
			aliases: ['DME', 'alcohol']
			formula: chem.parse('C2H6O')!
		},
	]
	assert find(sp, 'methane')! == 0
	assert find(sp, '  METHANE ')! == 0
	assert find(sp, 'Marsh Gas')! == 0
	assert find(sp, '74-82-8')! == 0
	assert find(sp, 'CH4')! == 0
	assert find(sp, 'dme')! == 2
	assert why(sp, 'C2H5OH') == 'C2H5OH is the formula of both Ethanol and Dimethyl ether: name one'
	assert why(sp, 'alcohol') == 'alcohol names both Ethanol and Dimethyl ether'
	assert why(sp, 'CH3') == 'no species here is CH3'
	assert why(sp, 'ethanoll') == 'no species ethanoll'
	if _ := find(sp, 'ethanoll') {
		assert false
	} else {
		e := err as NameError
		assert e.suggestions == ['Ethanol']
	}
}

fn why(sp []Species, name string) string {
	k := find(sp, name) or { return err.msg() }
	return 'found ${k}'
}

// A count of species too small for the reactions is refused.
fn test_species_out_of_range() {
	sp := species('CO', 'O2', 'CO2')
	rs := reactions(sp, 'CO + 0.5 O2 = CO2')
	if _ := dependence(rs, 2) {
		assert false
	} else {
		assert err.msg() == 'reforming has a species beyond the 2 given'
	}
}
