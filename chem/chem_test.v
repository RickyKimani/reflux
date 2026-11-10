module chem

fn listed(f Formula) string {
	mut parts := []string{}
	for c in f.counts {
		parts << '${symbol(c.z)}${c.n}'
	}
	return parts.join(' ')
}

// Formulas as ChemSep writes them, and as a user types one in an equation.
fn test_parse() {
	cases := {
		'CH4':                     'CH4'
		'HOH':                     'H2O'
		'H2O':                     'H2O'
		'OCO':                     'CO2'
		'CH3OH':                   'CH4O'
		'C2H5OH':                  'C2H6O'
		'CH3(CH2)2CH3':            'C4H10'
		'Ca[OH]2':                 'CaH2O2'
		'((CH3)2SiO)6':            'C12H36O6Si6'
		'NH3':                     'H3N'
		'Co':                      'Co'
		'CO':                      'CO'
		'(C6H5)SH':                'C6H6S'
		'CH3CH2CHC(CH3)CHO':       'C6H10O'
		'C(CH3)3CH2C(CH3)2CH2CH3': 'C10H22'
		'UF6':                     'F6U'
	}
	for text, hill in cases {
		f := parse(text) or { panic('${text}: ${err}') }
		assert f.str() == hill, '${text}: ${f}'
	}
	assert listed(parse('CH3OH')!) == 'H4 C1 O1'
	assert parse('C6H12O6')!.unicode() == 'C₆H₁₂O₆'
	assert parse('CH4')!.unicode() == 'CH₄'
	// Without carbon, hydrogen takes its place among the rest.
	assert parse('HCl')!.str() == 'ClH'
	assert parse('NaOH')!.str() == 'HNaO'
	assert parse('C10H22')! == parse('CH3(CH2)8CH3')!
	assert parse('C2H6O')! != parse('C2H4O')!
}

fn reason(text string) string {
	parse(text) or {
		e := err as FormulaError
		return '${e.column}: ${e.reason}'
	}
	return 'parsed'
}

fn test_refused() {
	assert reason('') == '1: an empty formula'
	assert reason('Xx') == '1: no element Xx'
	assert reason('Cx') == '1: no element Cx'
	assert reason('ch4') == '1: no element c: a symbol starts with a capital'
	assert reason('2H2') == '1: a count with nothing before it'
	assert reason('(N2)0.781') == '6: a count must be whole, as no compound has a fraction of an atom'
	assert reason('C0') == '2: a count of 0'
	assert reason('CH3(CH2') == '4: a ( never closed'
	assert reason('CH3)2') == '4: a ) with nothing open'
	assert reason('C(H]') == '4: a ] where ) closes what column 2 opened'
	assert reason('C()') == '2: an empty group'
	assert reason('C-H') == '2: a - in a formula'
	assert reason('C99999999999') == '2: a count out of range'
	assert reason('(C99999)99999') == '8: a count out of range'
}

fn test_weights() {
	assert weight(6)? == 12.011 * 1e-3
	if _ := weight(43) {
		assert false
	}
	if _ := weight(118) {
		assert false
	}
	m := parse('CH4')!.mass()!
	assert m == 4 * (1.008 * 1e-3) + 12.011 * 1e-3
	if _ := parse('Tc2O7')!.mass() {
		assert false
	}
	assert element('Cl')? == 17
	if _ := element('CL') {
		assert false
	}
	assert symbol(118) == 'Og'
	assert elements == 118
}

fn smiles_of(text string) string {
	s := smiles(text) or {
		e := err as FormulaError
		return '${e.column}: ${e.reason}'
	}
	q := if s.charge != 0 { ' ${s.charge}' } else { '' }
	return s.formula.str() + q
}

// Compositions RDKit gives for the same strings.
fn test_smiles() {
	cases := {
		'C':                            'CH4'
		'O':                            'H2O'
		'[H][H]':                       'H2'
		'[C+]#[O-]':                    'CO'
		'[C-]#[O+]':                    'CO'
		'O=C=O':                        'CO2'
		'CO':                           'CH4O'
		'N':                            'H3N'
		'N#N':                          'N2'
		'c1ccccc1':                     'C6H6'
		'C1=CC=CC=C1':                  'C6H6'
		'n1ccccc1':                     'C5H5N'
		'c1cc[nH]c1':                   'C4H5N'
		'Cn1cccc1':                     'C5H7N'
		'c1ccoc1':                      'C4H4O'
		'c1ccsc1':                      'C4H4S'
		'c1ccc2ccccc2c1':               'C10H8'
		'c1ccc2sccc2c1':                'C8H6S'
		'CC(=O)O':                      'C2H4O2'
		'C[C@@H](O)C(=O)O':             'C3H6O3'
		'F/C=C/F':                      'C2H2F2'
		'[13CH4]':                      'CH4'
		'[2H]O[2H]':                    'H2O'
		'[NH4+].[Cl-]':                 'ClH4N'
		'[O-][N+](=O)c1ccccc1':         'C6H5NO2'
		'CS(=O)(=O)C':                  'C2H6O2S'
		'OP(=O)(O)O':                   'H3O4P'
		'C%10CC%10':                    'C3H6'
		'C1CC1':                        'C3H6'
		'C=1CCC1':                      'C4H6'
		'C1CCC=1':                      'C4H6'
		'C=1CC=1':                      'C3H4'
		'BrC(Cl)(F)I':                  'CBrClFI'
		'[Fe+2]':                       'Fe 2'
		'[O--]':                        'O -2'
		'[O-2]':                        'O -2'
		'[S-2]':                        'S -2'
		'[Fe+3]':                       'Fe 3'
		'[Si](C)(C)(C)C':               'C4H12Si'
		'C[Si](C)(C)O[Si](C)(C)C':      'C6H18OSi2'
		'[C@H](F)(Cl)Br':               'CHBrClF'
		'[CH3:1][OH:2]':                'CH4O'
		'Cc1ccc(C)nc1':                 'C7H9N'
		'O=c1cc[nH]cc1':                'C5H5NO'
		'[se]1cccc1':                   'C4H4Se'
		'[He]':                         'He'
		'C(C(C(C(C(C(C(C(C(C)))))))))': 'C10H22'
		'OCC(O)CO':                     'C3H8O3'
		'C1CCC2(CC1)CCCC2':             'C10H18'
		'N1CCNCC1':                     'C4H10N2'
		'[O-]C(=O)C':                   'C2H3O2 -1'
		'C[N+](C)(C)C':                 'C4H12N 1'
		'S':                            'H2S'
		'S(=O)(=O)(O)O':                'H2O4S'
		'B(O)(O)O':                     'BH3O3'
		'c1ccccc1-c1ccccc1':            'C12H10'
		'C(=O)=O':                      'CO2'
		'[O]':                          'O'
		'[CH2]':                        'CH2'
		'Cl':                           'ClH'
		'ClCl':                         'Cl2'
		'C#C':                          'C2H2'
		'C$C':                          'C2'
		'CC.CC':                        'C4H12'
		'[OH3+]':                       'H3O 1'
		'[Na+].[O-]Cl(=O)(=O)=O':       'ClNaO4'
		'P':                            'H3P'
		'P(Cl)(Cl)(Cl)(Cl)Cl':          'Cl5P'
		'O=S(=O)=O':                    'O3S'
	}
	for text, want in cases {
		got := smiles_of(text)
		assert got == want, '${text}: ${got}, not ${want}'
	}
}

fn test_smiles_refused() {
	assert smiles_of('') == '1: an empty SMILES string'
	assert smiles_of('*C') == '1: a wildcard atom, *, has no element'
	assert smiles_of('C[*]') == '3: a wildcard atom, *, has no element'
	assert smiles_of('C(') == '2: a branch never closed'
	assert smiles_of('C)') == '2: a ) with no branch open'
	assert smiles_of('C()C') == '2: an empty branch'
	assert smiles_of('(C)') == '1: a branch with no atom before it'
	assert smiles_of('C=') == '2: a bond with no atom after it'
	assert smiles_of('=C') == '1: a bond with no atom before it'
	assert smiles_of('C==C') == '3: two bonds in a row'
	assert smiles_of('C=(C)') == '2: a bond before a branch: it goes inside'
	assert smiles_of('C1CC') == '2: ring bond 1 never closed'
	assert smiles_of('C11') == '3: ring bond 1 closes on the atom it opened from'
	assert smiles_of('C=1CC#1') == '7: ring bond 1 has one order where it opens, column 3, and another here'
	assert smiles_of('1CC') == '1: a ring bond with no atom before it'
	assert smiles_of('C%1C') == '2: a % without two digits after it'
	assert smiles_of('C.') == '2: a . with no atom after it'
	assert smiles_of('.C') == '1: a . with no atom before it'
	assert smiles_of('C=.C') == '2: a bond with no atom after it'
	assert smiles_of('[C') == '1: a [ never closed'
	assert smiles_of('[13') == '1: a [ never closed'
	assert smiles_of('[Xx]') == '2: no element X'
	assert smiles_of('[x]') == '2: no aromatic atom x'
	assert smiles_of('[]') == '2: a bracket atom with no element'
	assert smiles_of('[C@X]') == '4: a X in a bracket atom'
	assert smiles_of('[C:]') == '3: a class with no number'
	assert smiles_of('Xe') == '1: no atom X outside brackets'
	assert smiles_of('K') == '1: no atom K outside brackets'
	assert smiles_of('[K+]') == 'K 1'
	assert smiles_of('C C') == '2: a space in SMILES'
	assert smiles_of('o[H][H]') == 'H2O'
}
