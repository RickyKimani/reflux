module kinetics

import chem
import quantity

// Reaction is an equation read and balanced: its species, as positions in
// the species list, in the order written, each as written, and their
// coefficients, below zero on the left, exactly and as doubles.
pub struct Reaction {
pub:
	name    string
	text    string
	species []int
	written []string
	exact   []Ratio
	nu      []f64
}

// Imbalance is one element's count on each side of an equation that
// doesn't balance, or the charge's, element 'charge'.
pub struct Imbalance {
pub:
	element string
	left    Ratio
	right   Ratio
}

// ReactionError is every failure here: the reaction, its equation, the
// column (from 1, in characters; 0 where the whole equation is at fault),
// why, and names it may have meant. An equation that doesn't balance has
// each element's imbalance, and, where its species allow, the equation
// balanced.
pub struct ReactionError {
	Error
pub:
	reaction    string
	text        string
	column      int
	reason      string
	suggestions []string
	imbalances  []Imbalance
	balanced    string
	hint        string // why there's no balanced equation to suggest
}

// msg is the reason, with where it is, the imbalances and the suggestions.
pub fn (e ReactionError) msg() string {
	mut s := if e.reaction != '' { '${e.reaction}: ' } else { '' }
	if e.column > 0 {
		s += 'column ${e.column}: '
	}
	s += e.reason
	if e.suggestions.len > 0 {
		s += ' (did you mean ${e.suggestions.join(', ')}?)'
	}
	if e.imbalances.len > 0 {
		s += ': ' + e.imbalances.map('${it.element}: ${decimal(it.left)} on the left, ${decimal(it.right)} on the right').join('; ')
	}
	if e.balanced != '' {
		s += '; balanced, it would be ${e.balanced}'
	} else if e.hint != '' {
		s += '; ${e.hint}'
	}
	return s
}

// Term is one term of an equation as read: its coefficient, its species as
// written, where it starts, and its side.
struct Term {
	coeff  Ratio
	text   string
	column int
	right  bool
}

// arrows are what an equation might be written with instead of =.
const arrows = ['<=>', '<->', '->', '=>', '<-', '→', '⇌', '⇄', '↔', '⟶', '⟷']

// parse_equation reads a reaction's equation, CH4 + H2O = CO + 3 H2: terms
// joined by + on each side of one =, each a species with a coefficient
// before it if not 1. A coefficient is a whole number, a decimal or a
// fraction, 1/2, and may touch its species (3H2) where that isn't itself a
// species' name (1-butene). A species is named as find takes it, in quotes
// if its name has a + or an = in it. Each species appears once. The
// equation must balance every element, and the charge.
pub fn parse_equation(name string, text string, species []Species) !Reaction {
	raw := split(name, text)!
	mut terms := []Term{}
	mut idx := []int{}
	for t in raw {
		term, k := read_term(name, text, t, species)!
		p := idx.index(k)
		if p >= 0 {
			reason := if terms[p].right == term.right {
				'${term.text} appears twice on the ${side(term.right)}'
			} else {
				'${term.text} is on both sides: write its net coefficient on one'
			}
			return fail(name, text, term.column, reason)
		}
		terms << term
		idx << k
	}
	mut exact := []Ratio{}
	for t in terms {
		exact << if t.right { t.coeff } else { t.coeff.neg() }
	}
	r := Reaction{
		name:    name
		text:    text
		species: idx
		written: terms.map(it.text)
		exact:   exact
		nu:      exact.map(it.value())
	}
	r.balance(species)!
	return r
}

fn side(right bool) string {
	return if right { 'right' } else { 'left' }
}

fn fail(name string, text string, column int, reason string) IError {
	return ReactionError{
		reaction: name
		text:     text
		column:   column
		reason:   reason
	}
}

// RawTerm is a term's text between separators, and where it starts.
struct RawTerm {
	text   string
	column int
	right  bool
}

// split cuts an equation into its terms at each + and the =, outside
// quotes, and refuses an arrow.
fn split(name string, text string) ![]RawTerm {
	runes := text.runes()
	for a in arrows {
		if at := index_outside_quotes(runes, a.runes()) {
			return fail(name, text, at + 1, 'an equation is written with =, not ${a}: whether a reaction runs both ways is for its rate law to say')
		}
	}
	mut terms := []RawTerm{}
	mut start := 0
	mut right := false
	mut quote := rune(0)
	mut quote_at := 0
	for i := 0; i <= runes.len; i++ {
		if i < runes.len && quote != 0 {
			if runes[i] == quote {
				quote = 0
			}
			continue
		}
		if i < runes.len && (runes[i] == `'` || runes[i] == `"`) {
			quote, quote_at = runes[i], i
			continue
		}
		if i == runes.len || runes[i] == `+` || runes[i] == `=` {
			terms << RawTerm{runes[start..i].string(), start + 1, right}
			if i < runes.len && runes[i] == `=` {
				if right {
					return fail(name, text, i + 1, 'a second =')
				}
				right = true
			}
			start = i + 1
		}
	}
	if quote != 0 {
		return fail(name, text, quote_at + 1, 'a quote never closed')
	}
	if !right {
		return fail(name, text, 0, 'no = between the two sides')
	}
	return terms
}

fn index_outside_quotes(runes []rune, pattern []rune) ?int {
	mut quote := rune(0)
	for i in 0 .. runes.len {
		if quote != 0 {
			if runes[i] == quote {
				quote = 0
			}
			continue
		}
		if runes[i] == `'` || runes[i] == `"` {
			quote = runes[i]
			continue
		}
		if i + pattern.len <= runes.len && runes[i..i + pattern.len] == pattern {
			return i
		}
	}
	return none
}

// read_term reads one term: its coefficient and its species.
fn read_term(name string, text string, t RawTerm, species []Species) !(Term, int) {
	runes := t.text.runes()
	mut i := 0
	for i < runes.len && runes[i] == ` ` {
		i++
	}
	column := t.column + i
	if i == runes.len {
		return fail(name, text, column, 'a term with nothing in it')
	}
	mut coeff := one
	mut rest := runes[i..].string().trim_space()
	if runes[i] >= `0` && runes[i] <= `9` {
		mut j := i
		for j < runes.len && ((runes[j] >= `0` && runes[j] <= `9`) || runes[j] == `.` || runes[j] == `/`) {
			j++
		}
		number := runes[i..j].string()
		after := runes[j..].string()
		if after.trim_space() == '' {
			return fail(name, text, column, 'a coefficient, ${number}, with no species after it')
		}
		if after.starts_with(' ') {
			coeff = read_coefficient(name, text, column, number)!
			rest = after.trim_space()
		} else if _ := find(species, rest) {
			// The whole term is a species' name, 1-butene or 2,2-dimethylbutane.
		} else if _ := find(species, after) {
			coeff = read_coefficient(name, text, column, number)!
			rest = after
		}
	}
	mut written := rest
	if rest.starts_with("'") || rest.starts_with('"') {
		q := rest[0..1]
		end := (rest[1..].index(q) or { return fail(name, text, column, 'a quote never closed') }) + 1
		if rest[end + 1..].trim_space() != '' {
			return fail(name, text, column, 'something after the quoted name ${rest[..end + 1]}')
		}
		written = rest[1..end]
	}
	k := find(species, written) or {
		e := err as NameError
		return ReactionError{
			reaction:    name
			text:        text
			column:      column
			reason:      e.reason
			suggestions: e.suggestions
		}
	}
	return Term{coeff, written, column, t.right}, k
}

// read_coefficient reads a whole number, a decimal or a fraction, exactly.
fn read_coefficient(name string, text string, column int, s string) !Ratio {
	bad := fail(name, text, column, "a coefficient, ${s}, that isn't a number: write 2, 0.5 or 1/2")
	parts := s.split('/')
	if parts.len > 2 {
		return bad
	}
	for part in parts {
		if part == '' || part.count('.') > 1 || part.starts_with('.') || part.ends_with('.') {
			return bad
		}
	}
	if parts.len == 2 && parts[1].contains('.') {
		return bad
	}
	num := decimal_ratio(parts[0]) or {
		return fail(name, text, column, 'a coefficient, ${s}, out of range')
	}
	mut c := num
	if parts.len == 2 {
		den := decimal_ratio(parts[1]) or {
			return fail(name, text, column, 'a coefficient, ${s}, out of range')
		}
		if is_zero(den) {
			return fail(name, text, column, 'a coefficient, ${s}, over zero')
		}
		c = div(num, den) or { return fail(name, text, column, 'a coefficient, ${s}, out of range') }
	}
	if is_zero(c) {
		return fail(name, text, column, 'a coefficient of 0: leave the species out')
	}
	return c
}

// decimal_ratio reads digits with at most one point in them, exactly.
fn decimal_ratio(s string) !Ratio {
	dot := s.index('.') or { s.len }
	digits := s.replace('.', '')
	if digits.len > 7 {
		return error('out of range')
	}
	places := if dot < s.len { s.len - dot - 1 } else { 0 }
	mut den := i64(1)
	for _ in 0 .. places {
		den *= 10
	}
	return quantity.ratio(digits.i64(), den)
}

// balance checks that the equation balances every element and the charge,
// and says how it doesn't if it doesn't.
fn (r Reaction) balance(species []Species) ! {
	mut formulas := []chem.Formula{}
	for j, k in r.species {
		f := species[k].formula or {
			return fail(r.name, r.text, 0, "${r.written[j]} has no formula, as a mixture's pseudo-compound, so it can't react")
		}
		formulas << f
	}
	mut zs := []int{}
	for f in formulas {
		for c in f.counts {
			if c.z !in zs {
				zs << c.z
			}
		}
	}
	zs = order_hill(zs)
	charged := r.species.any(species[it].charge != 0)
	mut rows := [][]Ratio{}
	mut labels := []string{}
	for z in zs {
		rows << formulas.map(whole(it.count(z)))
		labels << chem.symbol(z)
	}
	if charged {
		rows << r.species.map(whole(species[it].charge))
		labels << 'charge'
	}
	mut off := []Imbalance{}
	for i, row in rows {
		mut left := zero
		mut right := zero
		for j, n in row {
			term := n.mul(r.exact[j]) or { return range_error(r) }
			if r.exact[j].num < 0 {
				left = left.add(term.neg()) or { return range_error(r) }
			} else {
				right = right.add(term) or { return range_error(r) }
			}
		}
		if left != right {
			off << Imbalance{labels[i], left, right}
		}
	}
	if off.len == 0 {
		return
	}
	balanced, hint := r.suggest(rows, labels)
	return ReactionError{
		reaction:   r.name
		text:       r.text
		reason:     "${r.text} doesn't balance"
		imbalances: off
		balanced:   balanced
		hint:       hint
	}
}

fn range_error(r Reaction) IError {
	return fail(r.name, r.text, 0, 'coefficients too large to balance exactly')
}

// order_hill orders atomic numbers as Hill's system does: carbon, then
// hydrogen if there's carbon, then the rest by symbol.
fn order_hill(zs []int) []int {
	mut counts := []chem.Count{}
	for z in zs {
		counts << chem.Count{z, 1}
	}
	counts.sort(a.z < b.z)
	return chem.Formula{counts}.hill().map(it.z)
}

// suggest finds the equation's balanced form from the null space of its
// element matrix, where the species allow just one, or says why there's
// none.
fn (r Reaction) suggest(rows [][]Ratio, labels []string) (string, string) {
	basis := null_space(rows, r.species.len) or { return '', '' }
	if basis.len == 0 {
		mut lonely := []string{}
		for i, row in rows {
			mut on_left, mut on_right := false, false
			for j, n in row {
				if !is_zero(n) {
					if r.exact[j].num < 0 {
						on_left = true
					} else {
						on_right = true
					}
				}
			}
			if on_left != on_right {
				lonely << '${labels[i]} is only on the ${side(on_right)}'
			}
		}
		if lonely.len > 0 {
			return '', 'no coefficients would balance it, as ${lonely.join(' and ')}'
		}
		return '', 'no coefficients would balance these species'
	}
	if basis.len > 1 {
		return '', "these species balance in ${basis.len} independent ways, so there's no one equation to suggest"
	}
	mut v := smallest_whole(basis[0]) or { return '', '' }
	// Keep the first species that takes part on the side it was written on.
	for j, x in v {
		if x != 0 {
			if (x < 0) != (r.exact[j].num < 0) {
				v = v.map(-it)
			}
			break
		}
	}
	mut left := []string{}
	mut right := []string{}
	mut dropped := []string{}
	mut moved := []string{}
	for j, x in v {
		w := if x == 1 || x == -1 {
			r.written[j]
		} else {
			'${if x < 0 { -x } else { x }} ${r.written[j]}'
		}
		if x < 0 {
			left << w
		} else if x > 0 {
			right << w
		}
		if x == 0 {
			dropped << r.written[j]
		} else if (x < 0) != (r.exact[j].num < 0) {
			moved << r.written[j]
		}
	}
	if left.len == 0 || right.len == 0 {
		return '', 'no coefficients would balance these species'
	}
	mut balanced := '${left.join(' + ')} = ${right.join(' + ')}'
	if dropped.len > 0 {
		balanced += ', leaving out ${dropped.join(' and ')}'
	}
	if moved.len > 0 {
		balanced += ', with ${moved.join(' and ')} on the other side'
	}
	return balanced, ''
}

// str writes the equation back with each species as written and each
// coefficient as a decimal or a fraction.
pub fn (r Reaction) str() string {
	mut left := []string{}
	mut right := []string{}
	for j, c in r.exact {
		a := if c.num < 0 { c.neg() } else { c }
		w := if a == one { r.written[j] } else { '${decimal(a)} ${r.written[j]}' }
		if c.num < 0 {
			left << w
		} else {
			right << w
		}
	}
	return '${left.join(' + ')} = ${right.join(' + ')}'
}
