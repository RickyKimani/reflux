module components

import chem
import math

// mass_tolerance is how far a formula's molar mass, from the standard
// atomic weights, may be from a component file's. ChemSep's masses use
// IUPAC's older weights (H 1.00794 against 1.008, S 32.065 against 32.06),
// which differ from CIAAW's of 2021 by at most sulfur's 1.6e-4, and a
// compound's mass by no more than its worst element's: carbon disulfide's
// by 1.3e-4, the most in the data. A hydrogen more or less is beyond it
// for any compound under 5,000 g/mol; CO and N2, 1.4e-4 apart, aren't.
pub const mass_tolerance = 2e-4

// Sources are what a compound's file says of its composition: the formula
// read from its structure and from its SMILES, or why each couldn't be,
// and its molar mass.
pub struct Sources {
pub:
	structure     ?chem.Formula
	structure_why string // why there's no structure formula
	smiles        ?chem.Smiles
	smiles_why    string
	mw            f64 // kg/mol, 0 if the file has none
}

// structure_text is ChemSep's structure as a formula: without the - that
// marks a ring's ends or a bond, or the spaces between groups.
fn structure_text(s string) string {
	return s.replace('-', '').replace(' ', '')
}

// sources reads the compound's formula from its structure and from its
// SMILES string.
pub fn (r Record) sources() Sources {
	mut structure := ?chem.Formula(none)
	mut structure_why := 'it has no structure'
	if st := r.ids.structure {
		if f := chem.parse(structure_text(st)) {
			structure, structure_why = f, ''
		} else {
			e := err as chem.FormulaError
			structure_why = if e.reason.starts_with('a count must be whole') {
				"its structure, ${st}, is a mixture's"
			} else {
				"its structure, ${st}, doesn't read as a formula: ${e.reason}"
			}
		}
	}
	mut smiles := ?chem.Smiles(none)
	mut smiles_why := 'it has no SMILES'
	if raw := r.ids.smiles {
		text := raw.trim_space()
		if s := chem.smiles(text) {
			smiles, smiles_why = s, ''
		} else {
			e := err as chem.FormulaError
			smiles_why = "its SMILES, ${text}, doesn't read: ${e.reason}"
		}
	}
	return Sources{
		structure:     structure
		structure_why: structure_why
		smiles:        smiles
		smiles_why:    smiles_why
		mw:            r.constants.mw or { 0.0 }
	}
}

// weighs is whether a formula's molar mass is the file's, within
// mass_tolerance.
pub fn (s Sources) weighs(f chem.Formula) bool {
	m := f.mass() or { return false }
	return s.mw > 0 && math.abs(m - s.mw) <= mass_tolerance * s.mw
}

// formula is the compound's formula, taken where two of its three sources
// agree: the structure and the SMILES string, or either of them and the
// molar mass. A compound with none that two agree on has no formula: a
// mixture's pseudo-compound, or a file with two of its three wrong.
pub fn (r Record) formula() !chem.Formula {
	s := r.sources()
	st := s.structure
	sm := s.smiles
	if a := st {
		if b := sm {
			if a == b.formula {
				return a
			}
		}
		if s.weighs(a) {
			return a
		}
	}
	if b := sm {
		if s.weighs(b.formula) {
			return b.formula
		}
	}
	mut said := []string{}
	if a := st {
		said << 'its structure gives ${a} (${grams(a)})'
	} else {
		said << s.structure_why
	}
	if b := sm {
		said << 'its SMILES gives ${b.formula} (${grams(b.formula)})'
	} else {
		said << s.smiles_why
	}
	said << if s.mw > 0 { 'its molar mass is ${s.mw * 1e3} g/mol' } else { 'it has no molar mass' }
	return error('${r.name} (${r.ids.cas}) has no formula that two of its structure, SMILES and molar mass agree on: ${said.join('; ')}')
}

fn grams(f chem.Formula) string {
	m := f.mass() or { return 'no molar mass' }
	return '${math.round(m * 1e7) / 1e4} g/mol'
}
