// Module kinetics is M9's chemistry of reactions: equations read with
// their species by name, formula or CAS number and balanced by counting
// elements, a set's independence, and each reaction's standard properties
// and equilibrium constant from the ideal gas (D7). Rate laws and sites
// follow in step 5.
module kinetics

import chem
import quantity
import thermo

// Species is a compound as reactions know it: the names it answers to,
// its formula and charge, and its ideal-gas data for standard properties.
// formula is none for a mixture's pseudo-compound, which can't react.
pub struct Species {
pub:
	name    string @[required]
	cas     string
	aliases []string
	formula ?chem.Formula
	charge  int
	gas     ?thermo.IdealGas
	dgf     ?f64 // ΔG°f at 298.15 K and 1 bar, J/mol: a check on ΔH°f and S°
}

// NameError is why a species couldn't be found, and what was close.
struct NameError {
	Error
	reason      string
	suggestions []string
}

fn (e NameError) msg() string {
	return e.reason
}

// find is the species a name means: its name or an alias, ignoring case,
// its CAS number, or its formula, CH3OH or C2H5OH as well as CH4O. A name
// two species answer to is refused, as is a formula two share.
pub fn find(species []Species, name string) !int {
	key := name.trim_space()
	lower := key.to_lower()
	mut found := []int{}
	for i, s in species {
		if s.name.to_lower() == lower || s.cas == key || s.aliases.any(it.to_lower() == lower) {
			found << i
		}
	}
	if found.len == 1 {
		return found[0]
	}
	if found.len > 1 {
		return NameError{
			reason: '${key} names ${both(species, found)}'
		}
	}
	if f := chem.parse(key) {
		for i, s in species {
			if g := s.formula {
				if g == f {
					found << i
				}
			}
		}
		if found.len == 1 {
			return found[0]
		}
		if found.len > 1 {
			return NameError{
				reason: '${key} is the formula of ${both(species, found)}: name one'
			}
		}
		return NameError{
			reason:      'no species here is ${f}'
			suggestions: near(species, key)
		}
	}
	return NameError{
		reason:      'no species ${key}'
		suggestions: near(species, key)
	}
}

// both lists the species found, "A and B" or "A, B and C".
fn both(species []Species, found []int) string {
	names := found.map(species[it].name)
	return if names.len == 2 {
		'both ${names[0]} and ${names[1]}'
	} else {
		'${names[..names.len - 1].join(', ')} and ${names.last()}'
	}
}

// near is the names and formulas closest to a name not found.
fn near(species []Species, name string) []string {
	mut accepted := []string{}
	for s in species {
		accepted << s.name
		accepted << s.aliases
		if f := s.formula {
			accepted << f.str()
		}
	}
	return quantity.closest(name, accepted)
}
