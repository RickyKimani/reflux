// Every compound's formula, from its structure, its SMILES and its molar
// mass (M9 step 4), against RDKit's reading of the same SMILES strings and
// chemicals' formulas by CAS number, from testdata/m9_formulas.json.
module components

import chem
import json2
import math
import os

// The built-in data, loaded once.
const everything = load_default() or { panic(err) }

struct FormulaCompound {
	cas             string
	name            string
	rdkit           ?string
	rdkit_charge    ?int
	rdkit_sanitized bool
	chemicals       ?string
}

struct FormulaFixture {
	weights   map[string]f64
	compounds []FormulaCompound
}

fn formula_fixture() !FormulaFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_formulas.json')
	return json2.decode[FormulaFixture](os.read_file(path)!)!
}

// The compounds whose files disagree with themselves, and which two of
// their three sources the formula comes from. ChemSep gives thiophenol
// benzothiophene's SMILES, ortho- and para-hydrogen SMILES that spell
// their names (o[H][H], p[H][H]), 2-methyl-2-pentenal its alcohol's SMILES,
// and dehydroabietic acid a structure with two hydrogens too many.
const disagreements = {
	'Phenyl mercaptan':    'structure'
	'Ortho-Hydrogen':      'structure'
	'Para-Hydrogen':       'structure'
	'2-Methyl-2-pentenal': 'structure'
	'Dehydroabietic acid': 'SMILES'
}

// Mixtures' pseudo-compounds, which have no formula.
const mixtures = ['Air', 'Dowtherm A']

// Each file's SMILES string reads to RDKit's composition, hydrogens and
// charge, except ortho- and para-hydrogen's, whose aromatic o and p RDKit
// gives one hydrogen and two that no valence asks for: they are a name, not
// a molecule.
fn test_smiles_read_as_rdkit_reads_them() {
	f := formula_fixture()!
	mut n := 0
	for c in f.compounds {
		want := c.rdkit or { continue }
		r := everything.find(c.cas)!
		s := r.sources()
		got := s.smiles or { panic('${c.name}: ${s.smiles_why}') }
		if c.name in ['Ortho-Hydrogen', 'Para-Hydrogen'] {
			assert got.formula.str() != want
			continue
		}
		assert got.formula.str() == want, '${c.name}: ${got.formula} vs ${want}'
		assert got.charge == c.rdkit_charge or { 0 }, c.name
		n++
	}
	assert n == 639
}

// Every compound but the mixtures has a formula, two of its sources
// agreeing, and it's the one chemicals' database gives for its CAS number
// wherever that has one. The files that disagree with themselves are
// exactly the ones listed, with the source taken.
fn test_formulas_are_two_sources_and_chemicals() {
	f := formula_fixture()!
	mut found := map[string]string{}
	mut checked := 0
	for c in f.compounds {
		r := everything.find(c.cas)!
		formula := r.formula() or {
			assert c.name in mixtures, '${err}'
			continue
		}
		assert c.name !in mixtures
		s := r.sources()
		st := s.structure or { panic(c.name) }
		sm := s.smiles or { panic(c.name) }
		if st != sm.formula {
			found[c.name] = if formula == st { 'structure' } else { 'SMILES' }
		}
		if want := c.chemicals {
			assert formula.str() == want, '${c.name}: ${formula} vs chemicals ${want}'
			checked++
		}
		assert s.weighs(formula), '${c.name}: ${formula.mass()!} vs ${s.mw}'
	}
	assert found == disagreements
	assert checked == 638
}

// The structure and the SMILES string agreeing are two sources, whatever the
// molar mass says; and the molar mass decides between them at 2e-4, where a
// hydrogen in 5,000 g/mol is beyond it.
fn test_two_of_three() {
	methane := everything.find('74-82-8')!
	mut wrong := methane
	wrong.constants.mw = 0.030
	assert wrong.formula()!.str() == 'CH4'
	f := chem.parse('CS2')!
	m := f.mass()!
	assert Sources{
		mw: m * (1 + 1.9e-4)
	}.weighs(f)
	assert !Sources{
		mw: m * (1 + 2.1e-4)
	}.weighs(f)
	assert !Sources{
		mw: m * (1 - 2.1e-4)
	}.weighs(f)
	// A structure a hydrogen short of its SMILES, in 563 g/mol: 1.8e-3 out,
	// so the SMILES string, which weighs, is taken.
	mut heavy := methane
	heavy.ids.structure = 'C40H81'
	heavy.ids.smiles = 'C'.repeat(40)
	heavy.constants.mw = chem.parse('C40H82')!.mass()!
	assert heavy.formula()!.str() == 'C40H82'
	heavy.ids.smiles = none
	if _ := heavy.formula() {
		assert false
	}
}

// A mixture's file, and the reasons each source gives.
fn test_no_formula_says_why() {
	air := everything.find('Air')!
	if _ := air.formula() {
		assert false
	} else {
		assert err.msg() == "Air (132259-10-0) has no formula that two of its structure, SMILES and molar mass agree on: its structure, (N2)0.781 (O2)0.209 (Ar)0.01, is a mixture's; it has no SMILES; its molar mass is 28.96 g/mol"
	}
	dow := everything.find('Dowtherm A')!
	if _ := dow.formula() {
		assert false
	} else {
		assert err.msg().contains('its SMILES gives C12H10O (170.211 g/mol)')
	}
	// Hydrogen sulfide's SMILES has a space after it, which reflux reads
	// past, as RDKit does.
	h2s := everything.find('7783-06-4')!
	assert h2s.ids.smiles or { '' } == 'S '
	assert h2s.formula()!.str() == 'H2S'
}

// The standard atomic weights reflux holds, the 84 elements below uranium
// that have one, against chemicals', to 2e-4: reflux's are CIAAW's of
// 2021, abridged, and chemicals' IUPAC's older ones, and the two differ by
// at most the 2021 table's conventional value for an interval (sulfur's
// 32.06 against 32.065, 1.6e-4).
fn test_atomic_weights() {
	f := formula_fixture()!
	mut worst := 0.0
	mut n := 0
	for sym, w in f.weights {
		z := chem.element(sym) or { panic(sym) }
		mine := chem.weight(z) or { continue }
		d := math.abs(mine * 1e3 - w) / w
		worst = math.max(worst, d)
		assert d <= 2e-4, '${sym}: ${mine * 1e3} vs ${w}'
		n++
	}
	println('atomic weights: ${n}, worst ${worst}')
	assert n == 84
}

// ΔG°f against ΔH°f − T°(S° − Σ nₑS°ₑ), by hand from the files and CODATA's
// entropies of the elements. Water: ΔS°f = 188.724 − 130.680 − 102.576 =
// −44.532 J/(mol·K), so −241814 + 298.15·44.532 = −228536.7842 J/mol, and
// the file's −228590 is 53.2158 below. Methane: ΔS°f = 186.27 − 5.74 −
// 261.36 = −80.83, so −74520 + 298.15·80.83 = −50420.5355, and the file's
// −50490 is 69.4645 below.
fn test_formation_gap() {
	water := everything.find('7732-18-5')!.formation_gap()!
	assert math.abs(water - -53.2158) < 1e-9
	methane := everything.find('74-82-8')!.formation_gap()!
	assert math.abs(methane - -69.4645) < 1e-9
	if _ := everything.find('Air')!.formation_gap() {
		assert false
	}
}
