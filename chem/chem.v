// Module chem knows the elements: their symbols and standard atomic
// weights, and a compound's element counts read from a formula, CH3OH or
// CH3(CH2)2CH3, or from a SMILES string. It depends on nothing, so that
// the component data can check a compound's formula with it and kinetics
// can read one typed in an equation.
module chem

// Count is how many atoms of one element a formula has.
pub struct Count {
pub:
	z int // the atomic number
	n int
}

// Formula is a compound's element counts, by atomic number, each above
// zero. Two formulas are the same compound's when they're equal, ==.
pub struct Formula {
pub:
	counts []Count
}

// FormulaError is every failure here: the text, the column (from 1, in
// characters) where it went wrong, and why.
pub struct FormulaError {
	Error
pub:
	text   string
	column int
	reason string
}

// msg is the reason, with where it is.
pub fn (e FormulaError) msg() string {
	return '${e.text}, column ${e.column}: ${e.reason}'
}

// symbols are the elements' symbols by atomic number, from 1.
const symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si',
	'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
	'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
	'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
	'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt',
	'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu',
	'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
	'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

// weights are the standard atomic weights in g/mol, by atomic number from
// 1, as CIAAW's 2021 table abridges them (Prohaska et al., Pure Appl.
// Chem. 94, 573, 2022), with the conventional value for an element whose
// weight is an interval (H, C, N, O and others). 0 is an element with no
// standard atomic weight, having no stable isotope.
const weights = [1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990,
	24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.95, 39.098, 40.078, 44.956, 47.867, 50.942,
	51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904,
	83.798, 85.468, 87.62, 88.906, 91.222, 92.906, 95.95, 0, 101.07, 102.91, 106.42, 107.87, 112.41,
	114.82, 118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24,
	0, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49,
	180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 0, 0,
	0, 0, 0, 0, 232.04, 231.04, 238.03]!

// count_limit bounds an element's count, far beyond any molecule's, so
// that no count overflows.
const count_limit = 1_000_000_000

// elements is how many elements there are.
pub const elements = symbols.len

// element is the atomic number of a symbol, or none.
pub fn element(symbol string) ?int {
	for i, s in symbols {
		if s == symbol {
			return i + 1
		}
	}
	return none
}

// symbol is element z's symbol.
pub fn symbol(z int) string {
	return symbols[z - 1]
}

// weight is element z's standard atomic weight in kg/mol, or none.
pub fn weight(z int) ?f64 {
	if z > weights.len || weights[z - 1] == 0 {
		return none
	}
	return weights[z - 1] * 1e-3
}

// count is how many atoms of element z the formula has.
pub fn (f Formula) count(z int) int {
	for c in f.counts {
		if c.z == z {
			return c.n
		}
	}
	return 0
}

// hill orders the counts as Hill's system writes them: carbon first, then
// hydrogen, then the rest by symbol; without carbon, every element by
// symbol, hydrogen among them.
pub fn (f Formula) hill() []Count {
	carbon := f.count(6) > 0
	mut first := []Count{}
	mut rest := []Count{}
	for c in f.counts {
		if carbon && (c.z == 6 || c.z == 1) {
			first << c
		} else {
			rest << c
		}
	}
	if first.len == 2 && first[0].z == 1 {
		first[0], first[1] = first[1], first[0]
	}
	rest.sort_with_compare(by_symbol)
	first << rest
	return first
}

fn by_symbol(a &Count, b &Count) int {
	return compare_strings(symbols[a.z - 1], symbols[b.z - 1])
}

// str writes the formula in Hill's order, CH4O.
pub fn (f Formula) str() string {
	mut s := ''
	for c in f.hill() {
		s += symbols[c.z - 1] + if c.n > 1 { '${c.n}' } else { '' }
	}
	return s
}

const subscripts = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']!

// unicode writes the formula in Hill's order with subscripts, CH₄O.
pub fn (f Formula) unicode() string {
	mut s := ''
	for c in f.hill() {
		s += symbols[c.z - 1]
		if c.n > 1 {
			for d in '${c.n}' {
				s += subscripts[d - `0`]
			}
		}
	}
	return s
}

// mass is the formula's molar mass from the standard atomic weights,
// kg/mol, or an error naming an element that has none.
pub fn (f Formula) mass() !f64 {
	mut m := 0.0
	for c in f.counts {
		w := weight(c.z) or {
			return error('${symbols[c.z - 1]} has no standard atomic weight')
		}
		m += f64(c.n) * w
	}
	return m
}

// formula is a Formula from counts by atomic number, zeros left out.
fn formula(counts []int) Formula {
	mut out := []Count{}
	for i, n in counts {
		if n > 0 {
			out << Count{i, n}
		}
	}
	return Formula{out}
}

// parse reads a formula: element symbols, each with a whole count after it
// if not 1, and groups in ( ) or [ ] with a count after them, as CH3OH,
// C2H5OH, CH3(CH2)2CH3 or Ca[OH]2. Case matters: Co is cobalt, CO carbon
// monoxide.
pub fn parse(text string) !Formula {
	runes := text.runes()
	if runes.len == 0 {
		return fail(text, 1, 'an empty formula')
	}
	// stack holds each open group's counts, and opens where it opened.
	mut stack := [][]int{len: 1, init: []int{len: symbols.len + 1}}
	mut opens := []int{}
	mut closers := []rune{}
	mut i := 0
	for i < runes.len {
		c := runes[i]
		if c == `(` || c == `[` {
			stack << []int{len: symbols.len + 1}
			opens << i
			closers << if c == `(` { `)` } else { `]` }
			i++
			continue
		}
		if c == `)` || c == `]` {
			if opens.len == 0 {
				return fail(text, i + 1, 'a ${c} with nothing open')
			}
			if c != closers.last() {
				return fail(text, i + 1, 'a ${c} where ${closers.last()} closes what column ${opens.last() + 1} opened')
			}
			start := opens.pop()
			closers.delete_last()
			inner := stack.pop()
			if !inner.any(it > 0) {
				return fail(text, start + 1, 'an empty group')
			}
			n, next := multiplier(text, runes, i + 1)!
			for z, k in inner {
				stack[stack.len - 1][z] = add(text, i + 1, stack.last()[z], i64(k) * i64(n))!
			}
			i = next
			continue
		}
		if c >= `A` && c <= `Z` {
			mut sym := c.str()
			mut next := i + 1
			if next < runes.len && runes[next] >= `a` && runes[next] <= `z` {
				two := sym + runes[next].str()
				if _ := element(two) {
					sym = two
					next++
				}
			}
			if next < runes.len && runes[next] >= `a` && runes[next] <= `z` {
				return fail(text, i + 1, 'no element ${sym}${runes[next]}')
			}
			z := element(sym) or { return fail(text, i + 1, 'no element ${sym}') }
			n, after := multiplier(text, runes, next)!
			stack[stack.len - 1][z] = add(text, i + 1, stack.last()[z], i64(n))!
			i = after
			continue
		}
		if c >= `0` && c <= `9` {
			return fail(text, i + 1, 'a count with nothing before it')
		}
		if c >= `a` && c <= `z` {
			return fail(text, i + 1, 'no element ${c}: a symbol starts with a capital')
		}
		return fail(text, i + 1, 'a ${c} in a formula')
	}
	if opens.len > 0 {
		return fail(text, opens.last() + 1, 'a ${runes[opens.last()]} never closed')
	}
	return formula(stack[0])
}

// multiplier reads the whole count at runes[i], 1 if there's none, and
// where the formula goes on.
fn multiplier(text string, runes []rune, i int) !(int, int) {
	mut j := i
	mut n := i64(0)
	for j < runes.len && runes[j] >= `0` && runes[j] <= `9` {
		n = n * 10 + i64(runes[j] - `0`)
		if n > count_limit {
			return fail(text, i + 1, 'a count out of range')
		}
		j++
	}
	if j == i {
		return 1, i
	}
	if j < runes.len && runes[j] == `.` {
		return fail(text, j + 1, 'a count must be whole, as no compound has a fraction of an atom')
	}
	if n == 0 {
		return fail(text, i + 1, 'a count of 0')
	}
	return int(n), j
}

// add is a + b, refused out of range.
fn add(text string, column int, a int, b i64) !int {
	if i64(a) + b > count_limit {
		return fail(text, column, 'a count out of range')
	}
	return a + int(b)
}

fn fail(text string, column int, reason string) IError {
	return FormulaError{
		text:   text
		column: column
		reason: reason
	}
}
