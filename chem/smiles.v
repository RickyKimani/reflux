module chem

// Smiles is what a SMILES string says of a compound's composition: its
// element counts, hydrogens included, and its net charge.
pub struct Smiles {
pub:
	formula Formula
	charge  int
}

// SmAtom is one atom of a SMILES string as read: its element, whether it
// was written aromatic or in brackets, a bracket atom's hydrogens and
// charge, and the sum of its bonds' orders, an aromatic bond counting 1.
struct SmAtom {
	z         int
	aromatic  bool
	bracket   bool
	hydrogens int
	charge    int
mut:
	bonds int
}

// Ring is a ring bond left open: its atom, its order if one was written
// (0 if not), and where.
struct Ring {
	atom   int
	order  int
	column int
}

// valences are each organic-subset element's normal valences, lowest
// first, by atomic number: B, C, N, O, P, S and the halogens.
fn valences(z int) []int {
	return match z {
		5 { [3] }
		6 { [4] }
		7 { [3, 5] }
		8 { [2] }
		15 { [3, 5] }
		16 { [2, 4, 6] }
		9, 17, 35, 53 { [1] }
		else { []int{} }
	}
}

// implicit is how many hydrogens an organic-subset atom has that SMILES
// leaves unwritten, by OpenSMILES' rule: enough to bring its bonds to the
// lowest normal valence that holds them, and none if none does. An
// aromatic atom gives one bond to its ring, so it takes one hydrogen
// fewer than that, as benzene's c takes one and pyridine's n none.
fn (a SmAtom) implicit() int {
	vs := valences(a.z)
	if a.aromatic {
		h := vs[0] - a.bonds - 1
		return if h > 0 { h } else { 0 }
	}
	for v in vs {
		if v >= a.bonds {
			return v - a.bonds
		}
	}
	return 0
}

// smiles reads a SMILES string's composition: atoms in the organic subset
// (B C N O P S F Cl Br I, and aromatic b c n o p s) with their implicit
// hydrogens, atoms in brackets with an isotope, chirality, hydrogens,
// charge and class, bonds, branches, ring bonds and dots. It checks the
// string's syntax and its ring bonds, not its chemistry.
pub fn smiles(text string) !Smiles {
	runes := text.runes()
	mut atoms := []SmAtom{}
	mut prev := -1
	mut order := 0 // a bond written before the next atom, 0 if none
	mut order_at := 0
	mut branches := []int{} // the atom each open branch hangs from
	mut branch_atoms := []int{} // how many atoms there were when it opened
	mut branch_at := []int{}
	mut rings := map[int]Ring{}
	mut dot := -1 // a . with no atom after it yet
	mut i := 0
	for i < runes.len {
		c := runes[i]
		if c == `(` {
			if prev < 0 {
				return fail(text, i + 1, 'a branch with no atom before it')
			}
			if order != 0 {
				return fail(text, order_at + 1, 'a bond before a branch: it goes inside')
			}
			branches << prev
			branch_atoms << atoms.len
			branch_at << i
			i++
			continue
		}
		if c == `)` {
			if branches.len == 0 {
				return fail(text, i + 1, 'a ) with no branch open')
			}
			if order != 0 {
				return fail(text, order_at + 1, 'a bond with no atom after it')
			}
			if atoms.len == branch_atoms.last() {
				return fail(text, branch_at.last() + 1, 'an empty branch')
			}
			prev = branches.pop()
			branch_atoms.delete_last()
			branch_at.delete_last()
			i++
			continue
		}
		bond := match c {
			`-`, `/`, `\\`, `:` { 1 }
			`=` { 2 }
			`#` { 3 }
			`$` { 4 }
			else { 0 }
		}
		if bond > 0 {
			if prev < 0 {
				return fail(text, i + 1, 'a bond with no atom before it')
			}
			if order != 0 {
				return fail(text, i + 1, 'two bonds in a row')
			}
			order, order_at = bond, i
			i++
			continue
		}
		if (c >= `0` && c <= `9`) || c == `%` {
			mut num := 0
			start := i
			if c == `%` {
				if i + 2 >= runes.len || !digit(runes[i + 1]) || !digit(runes[i + 2]) {
					return fail(text, i + 1, 'a % without two digits after it')
				}
				num = int(runes[i + 1] - `0`) * 10 + int(runes[i + 2] - `0`)
				i += 3
			} else {
				num = int(c - `0`)
				i++
			}
			if prev < 0 {
				return fail(text, start + 1, 'a ring bond with no atom before it')
			}
			if open := rings[num] {
				if open.atom == prev {
					return fail(text, start + 1, 'ring bond ${num} closes on the atom it opened from')
				}
				if open.order != 0 && order != 0 && open.order != order {
					return fail(text, start + 1, 'ring bond ${num} has one order where it opens, column ${open.column + 1}, and another here')
				}
				o := if order != 0 {
					order
				} else if open.order != 0 {
					open.order
				} else {
					1
				}
				atoms[open.atom].bonds += o
				atoms[prev].bonds += o
				rings.delete(num)
			} else {
				rings[num] = Ring{prev, order, start}
			}
			order = 0
			continue
		}
		if c == `.` {
			if prev < 0 {
				return fail(text, i + 1, 'a . with no atom before it')
			}
			if order != 0 {
				return fail(text, order_at + 1, 'a bond with no atom after it')
			}
			prev, dot = -1, i
			i++
			continue
		}
		atom, next := read_atom(text, runes, i)!
		dot = -1
		atoms << atom
		k := atoms.len - 1
		if prev >= 0 {
			o := if order != 0 { order } else { 1 }
			atoms[prev].bonds += o
			atoms[k].bonds += o
		}
		prev, order = k, 0
		i = next
	}
	if atoms.len == 0 {
		return fail(text, 1, 'an empty SMILES string')
	}
	if order != 0 {
		return fail(text, order_at + 1, 'a bond with no atom after it')
	}
	if dot >= 0 {
		return fail(text, dot + 1, 'a . with no atom after it')
	}
	if branches.len > 0 {
		return fail(text, branch_at.last() + 1, 'a branch never closed')
	}
	for num, open in rings {
		return fail(text, open.column + 1, 'ring bond ${num} never closed')
	}
	mut counts := []int{len: symbols.len + 1}
	mut charge := 0
	for a in atoms {
		if i64(counts[a.z]) + 1 > count_limit || i64(counts[1]) + 9 > count_limit {
			return fail(text, 1, 'too many atoms')
		}
		counts[a.z]++
		counts[1] += if a.bracket { a.hydrogens } else { a.implicit() }
		charge += a.charge
	}
	return Smiles{formula(counts), charge}
}

// read_atom reads the atom at runes[i], in the organic subset or in
// brackets, and where the string goes on.
fn read_atom(text string, runes []rune, i int) !(SmAtom, int) {
	c := runes[i]
	if c == `[` {
		return read_bracket(text, runes, i)
	}
	if c == `*` {
		return fail(text, i + 1, 'a wildcard atom, *, has no element')
	}
	if c == `C` && i + 1 < runes.len && runes[i + 1] == `l` {
		return SmAtom{
			z: 17
		}, i + 2
	}
	if c == `B` && i + 1 < runes.len && runes[i + 1] == `r` {
		return SmAtom{
			z: 35
		}, i + 2
	}
	z := match c {
		`B`, `b` { 5 }
		`C`, `c` { 6 }
		`N`, `n` { 7 }
		`O`, `o` { 8 }
		`P`, `p` { 15 }
		`S`, `s` { 16 }
		`F` { 9 }
		`I` { 53 }
		else { 0 }
	}
	if z == 0 {
		if c == ` ` || c == `\t` || c == `\n` || c == `\r` {
			return fail(text, i + 1, 'a space in SMILES')
		}
		if (c >= `a` && c <= `z`) || (c >= `A` && c <= `Z`) {
			return fail(text, i + 1, 'no atom ${c} outside brackets')
		}
		return fail(text, i + 1, 'a ${c} in SMILES')
	}
	return SmAtom{
		z:        z
		aromatic: c >= `a` && c <= `z`
	}, i + 1
}

// read_bracket reads an atom in brackets, [13CH4] or [O-] or [C@@H]:
// isotope, element, chirality, hydrogens, charge and class, in that order.
fn read_bracket(text string, runes []rune, start int) !(SmAtom, int) {
	mut i := start + 1
	for i < runes.len && digit(runes[i]) {
		i++
	}
	if i >= runes.len {
		return fail(text, start + 1, 'a [ never closed')
	}
	mut z := 0
	mut aromatic := false
	at := i
	if runes[i] == `*` {
		return fail(text, i + 1, 'a wildcard atom, *, has no element')
	}
	if runes[i] >= `a` && runes[i] <= `z` {
		aromatic = true
		two := if i + 1 < runes.len { runes[i..i + 2].string() } else { '' }
		if two in ['se', 'as', 'te'] {
			z = element(two.capitalize()) or { 0 }
			i += 2
		} else if runes[i] in [`b`, `c`, `n`, `o`, `p`, `s`] {
			z = element(runes[i].str().to_upper()) or { 0 }
			i++
		} else {
			return fail(text, at + 1, 'no aromatic atom ${runes[i]}')
		}
	} else if runes[i] >= `A` && runes[i] <= `Z` {
		if i + 1 < runes.len && runes[i + 1] >= `a` && runes[i + 1] <= `z` {
			if two := element(runes[i..i + 2].string()) {
				z = two
				i += 2
			}
		}
		if z == 0 {
			z = element(runes[i].str()) or { return fail(text, at + 1, 'no element ${runes[i]}') }
			i++
		}
	} else {
		return fail(text, at + 1, 'a bracket atom with no element')
	}
	// Chirality: @ or @@, or @ and a class with its number, as @TH1.
	if i < runes.len && runes[i] == `@` {
		i++
		if i < runes.len && runes[i] == `@` {
			i++
		} else if i + 2 < runes.len && runes[i..i + 2].string() in ['TH', 'AL', 'SP', 'TB', 'OH']
			&& digit(runes[i + 2]) {
			i += 2
			for i < runes.len && digit(runes[i]) {
				i++
			}
		}
	}
	mut hydrogens := 0
	if i < runes.len && runes[i] == `H` {
		hydrogens = 1
		i++
		if i < runes.len && digit(runes[i]) {
			hydrogens = int(runes[i] - `0`)
			i++
		}
	}
	mut charge := 0
	if i < runes.len && (runes[i] == `+` || runes[i] == `-`) {
		sign := if runes[i] == `+` { 1 } else { -1 }
		i++
		if i < runes.len && digit(runes[i]) {
			mut n := int(runes[i] - `0`)
			i++
			if i < runes.len && digit(runes[i]) {
				n = n * 10 + int(runes[i] - `0`)
				i++
			}
			charge = sign * n
		} else {
			charge = sign
			for i < runes.len && runes[i] == runes[i - 1] && charge > -9 && charge < 9 {
				charge += sign
				i++
			}
		}
	}
	if i < runes.len && runes[i] == `:` {
		i++
		if i >= runes.len || !digit(runes[i]) {
			return fail(text, i, 'a class with no number')
		}
		for i < runes.len && digit(runes[i]) {
			i++
		}
	}
	if i >= runes.len || runes[i] != `]` {
		if i >= runes.len {
			return fail(text, start + 1, 'a [ never closed')
		}
		return fail(text, i + 1, 'a ${runes[i]} in a bracket atom')
	}
	return SmAtom{
		z:         z
		aromatic:  aromatic
		bracket:   true
		hydrogens: hydrogens
		charge:    charge
	}, i + 1
}

fn digit(r rune) bool {
	return r >= `0` && r <= `9`
}
