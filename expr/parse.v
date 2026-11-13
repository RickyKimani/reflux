module expr

import math
import quantity
import strconv

// The grammar, by precedence, loosest first. ^ binds tighter than a minus
// sign and to the right, as in mathematics and Python: -x^2 is −(x²), x^-2
// is x⁻², and a^b^c is a^(b^c). Division and multiplication go left to
// right, so a/b*c is (a/b)·c, as SymPy reads it.
//
//     expr   := term (("+" | "-") term)*
//     term   := factor (("*" | "·" | "×" | "/") factor)*
//     factor := ("-" | "+") factor | power
//     power  := atom ["in" unit] ["^" factor]
//     atom   := number [unit] | name | name "[" key "]" | function "(" expr ")"
//               | "(" expr ")" ["[" unit "]"]
//
// A number's unit follows it after a space and runs to the next space, or
// to a ) or , that isn't its own: "240.1 kJ/mol / (R*T)". Or it's in
// brackets, "240.1 [kJ/mol]", which is the way to write a unit that is
// also a name, such as R. A key is letters, digits and _, or quoted, as
// p['n-butane'].
//
// "in" reads the value just before it on a unit's scale, a pure number:
// (T in degC) is 320 at 593.15 K, binding as a unit binds to its number, so
// 0.0018*T in degC is 0.0018*(T in degC). A unit with an offset (°C, °F,
// barg, psig) reads only a context variable, T or P, which is absolute: a
// difference reads in K or Pa, and a definition may be one. A group may
// take a unit in brackets after it, (85 - 0.18*(T in degC)) [%], its value
// in that unit.

// The deepest nesting the parser follows, so that no text can exhaust the
// stack.
const max_nesting = 200

struct Parser {
	ctx   Context
	defs  []string
	text  string
	runes []rune
mut:
	pos   int
	nodes []Node
	depth int
	lines []int // where each line starts
}

fn new_parser(ctx Context, defs []string, text string) Parser {
	runes := text.runes()
	mut lines := [0]
	for i, c in runes {
		if c == `\n` {
			lines << i + 1
		}
	}
	return Parser{
		ctx:   ctx
		defs:  defs
		text:  text
		runes: runes
		lines: lines
	}
}

// where is the line and column, from 1, of a position.
fn (p &Parser) where(pos int) (int, int) {
	mut line := 0
	for line + 1 < p.lines.len && p.lines[line + 1] <= pos {
		line++
	}
	return line + 1, pos - p.lines[line] + 1
}

fn (p &Parser) fail(pos int, reason string) IError {
	return p.fail_with(pos, reason, []string{})
}

fn (p &Parser) fail_with(pos int, reason string, suggestions []string) IError {
	line, column := p.where(pos)
	return ExprError{
		text:        p.text
		line:        line
		column:      column
		reason:      reason
		suggestions: suggestions
	}
}

fn (mut p Parser) add(n Node, pos int) int {
	line, column := p.where(pos)
	p.nodes << Node{
		...n
		line:   line
		column: column
	}
	return p.nodes.len - 1
}

fn is_space(c rune) bool {
	return c == ` ` || c == `\t` || c == `\n` || c == `\r`
}

fn is_digit(c rune) bool {
	return c >= `0` && c <= `9`
}

fn is_letter(c rune) bool {
	return (c >= `a` && c <= `z`) || (c >= `A` && c <= `Z`)
}

fn is_name_char(c rune) bool {
	return is_letter(c) || is_digit(c) || c == `_`
}

// starts_unit is whether c can begin a unit written after a number.
fn starts_unit(c rune) bool {
	return is_letter(c) || c in [`°`, `µ`, `μ`, `%`]
}

fn (mut p Parser) skip() {
	for p.pos < p.runes.len && is_space(p.runes[p.pos]) {
		p.pos++
	}
}

fn (p &Parser) peek() rune {
	return if p.pos < p.runes.len { p.runes[p.pos] } else { 0 }
}

fn (mut p Parser) parse() !Expr {
	p.skip()
	if p.pos >= p.runes.len {
		return p.fail(0, 'the expression is empty')
	}
	root := p.expr()!
	p.skip()
	if p.pos < p.runes.len {
		return p.after_value()
	}
	return Expr{
		text:  p.text
		nodes: p.nodes
		root:  root
	}
}

fn (mut p Parser) expr() !int {
	mut left := p.term()!
	for {
		p.skip()
		c := p.peek()
		if c != `+` && c != `-` && c != `−` {
			return left
		}
		at := p.pos
		p.pos++
		right := p.operand(at, c)!
		right_term := p.term_from(right)!
		left = p.add(Node{
			kind: if c == `+` { Kind.add } else { Kind.sub }
			a:    left
			b:    right_term
		}, at)
	}
	return left
}

// operand reads what an operator needs after it, refusing nothing there.
fn (mut p Parser) operand(at int, op rune) !int {
	p.skip()
	if p.pos >= p.runes.len {
		return p.fail(at, 'a ${op} with nothing after it')
	}
	c := p.peek()
	if c in [`*`, `/`, `^`, `·`, `×`, `)`, `]`, `,`] {
		return p.fail(p.pos, 'a ${op} with nothing after it, before ${c}')
	}
	return p.factor()
}

fn (mut p Parser) term() !int {
	first := p.factor()!
	return p.term_from(first)
}

// term_from continues a term whose first factor has been read.
fn (mut p Parser) term_from(first int) !int {
	mut left := first
	for {
		p.skip()
		c := p.peek()
		if c !in [`*`, `·`, `×`, `/`] {
			return left
		}
		at := p.pos
		p.pos++
		right := p.operand(at, c)!
		left = p.add(Node{
			kind: if c == `/` { Kind.div } else { Kind.mul }
			a:    left
			b:    right
		}, at)
	}
	return left
}

fn (mut p Parser) factor() !int {
	p.skip()
	if p.depth >= max_nesting {
		return p.fail(p.pos, 'the expression is nested too deeply')
	}
	c := p.peek()
	if c == `-` || c == `−` || c == `+` {
		at := p.pos
		p.pos++
		p.depth++
		inner := p.operand(at, c)!
		p.depth--
		if c == `+` {
			return inner
		}
		return p.add(Node{
			kind: .neg
			a:    inner
		}, at)
	}
	return p.power()
}

fn (mut p Parser) power() !int {
	p.skip()
	start := p.pos
	atom := p.atom()!
	base := p.reading(atom, start)!
	p.skip()
	if p.peek() != `^` {
		return base
	}
	at := p.pos
	p.pos++
	p.depth++
	exponent := p.operand(at, `^`)!
	p.depth--
	return p.add(Node{
		kind: .pow
		a:    base
		b:    exponent
	}, at)
}

fn (mut p Parser) atom() !int {
	p.skip()
	if p.pos >= p.runes.len {
		return p.fail(p.pos, 'the expression ends where a value should be')
	}
	c := p.peek()
	if c == `(` {
		open := p.pos
		p.pos++
		p.skip()
		if p.peek() == `)` {
			return p.fail(open, 'empty parentheses')
		}
		p.depth++
		inner := p.expr()!
		p.depth--
		p.skip()
		if p.peek() != `)` {
			if p.pos >= p.runes.len {
				return p.fail(open, 'a ( with no ) to close it')
			}
			return p.after_value()
		}
		p.pos++
		return p.unit_after(inner, open)
	}
	if is_digit(c) || (c == `.` && p.pos + 1 < p.runes.len && is_digit(p.runes[p.pos + 1])) {
		return p.number()
	}
	if is_letter(c) || c == `_` {
		return p.name()
	}
	return match c {
		`)` { p.fail(p.pos, 'a ) with no ( before it, or nothing inside') }
		`]` { p.fail(p.pos, 'a ] with no [ before it') }
		`,` { p.fail(p.pos, 'a , where a value should be: a function here takes one argument') }
		`*`, `/`, `^`, `·`, `×` { p.fail(p.pos, 'a ${c} with nothing before it') }
		`[` {
			p.fail(p.pos, 'a [ with no name before it: a unit in brackets follows its number, as 2 [K]')
		}
		else { p.fail(p.pos, 'unexpected ${c}') }
	}
}

// after_value explains what follows a complete value and can't: a ), a
// second value side by side, or a character with no meaning here.
fn (p &Parser) after_value() IError {
	c := p.peek()
	if c == `)` {
		return p.fail(p.pos, 'a ) with no ( before it')
	}
	if c == `,` {
		return p.fail(p.pos, 'a , here: a function takes one argument, and there are no lists')
	}
	if is_digit(c) || is_letter(c) || c == `(` || c == `_` {
		return p.fail(p.pos, 'two values side by side: put an operator between them, such as *')
	}
	if c == `=` {
		return p.fail(p.pos, 'an = here: an expression gives a value, and assigns nothing')
	}
	return p.fail(p.pos, 'unexpected ${c}')
}

// number reads a number and any unit after it.
fn (mut p Parser) number() !int {
	start := p.pos
	for p.pos < p.runes.len && is_digit(p.runes[p.pos]) {
		p.pos++
	}
	if p.peek() == `.` {
		p.pos++
		for p.pos < p.runes.len && is_digit(p.runes[p.pos]) {
			p.pos++
		}
	}
	if p.peek() == `e` || p.peek() == `E` {
		mut k := p.pos + 1
		if k < p.runes.len && (p.runes[k] == `+` || p.runes[k] == `-`) {
			k++
		}
		if k < p.runes.len && is_digit(p.runes[k]) {
			p.pos = k
			for p.pos < p.runes.len && is_digit(p.runes[p.pos]) {
				p.pos++
			}
		}
	}
	text := p.runes[start..p.pos].string()
	// The text is a number's by now, so strconv refuses it only out of range.
	value := strconv.atof64(text) or { return p.fail(start, '${text} is beyond the range of a number') }
	if !math.is_finite(value) {
		return p.fail(start, '${text} is beyond the range of a number')
	}
	// A unit: in brackets, or after a space.
	gap := p.pos
	p.skip()
	spaced := p.pos > gap
	c := p.peek()
	if c == `[` {
		open := p.pos
		p.pos++
		from := p.pos
		for p.pos < p.runes.len && p.runes[p.pos] != `]` {
			p.pos++
		}
		if p.pos >= p.runes.len {
			return p.fail(open, 'a [ with no ] to close its unit')
		}
		unit_text := p.runes[from..p.pos].string().trim_space()
		p.pos++
		return p.with_unit(start, text, value, unit_text, from, true)
	}
	if !spaced && (is_letter(c) || c == `_`) {
		word := p.word()
		if _ := quantity.unit(word) {
			return p.fail(gap, 'put a space between a number and its unit: ${text} ${word}')
		}
		return p.fail(gap, 'a number and a name run together: write ${text}*${word.all_before('(')}')
	}
	if spaced && p.at_in() && p.unit_follows(p.pos + 2) {
		return p.fail(p.pos, "a number before in: in reads a variable on a unit's scale, as (T in degC), and a number takes its unit after it, as ${text} K (for inches, ${text} [in])")
	}
	if spaced && starts_unit(c) {
		from := p.pos
		word := p.word()
		lead := leading_name(word)
		if lead != '' && (lead == word || word[lead.len..].starts_with('(')
			|| word[lead.len..].starts_with('[')) && p.known(lead) {
			if lead == word {
				if _ := quantity.unit(word) {
					return p.fail(from, '${text} ${word}: ${word} is ${p.meaning(word)} and also a unit: write ${text}*${word} to multiply, or ${text} [${word}] for the unit')
				}
			}
			return p.fail(from, 'a number and a name side by side: write ${text}*${lead} to multiply')
		}
		return p.with_unit(start, text, value, word, from, false)
	}
	p.pos = gap
	return p.add(Node{
		kind:  .num
		value: value
		text:  text
	}, start)
}

// word reads a unit written after a number: up to a space, a comma, or a )
// or ] that closes something opened before it.
fn (mut p Parser) word() string {
	from := p.pos
	mut depth := 0
	for p.pos < p.runes.len {
		c := p.runes[p.pos]
		if is_space(c) || c == `,` || c == `]` {
			break
		}
		if c == `(` {
			depth++
		} else if c == `)` {
			if depth == 0 {
				break
			}
			depth--
		}
		p.pos++
	}
	return p.runes[from..p.pos].string()
}

fn leading_name(word string) string {
	mut i := 0
	for i < word.len && (word[i].is_letter() || word[i].is_digit() || word[i] == `_`) {
		i++
	}
	if i == 0 || word[0].is_digit() {
		return ''
	}
	return word[..i]
}

// known is whether a name means something here: a function, a constant, a
// variable, a family or a definition.
fn (p &Parser) known(name string) bool {
	return reserved(p.ctx, name) != '' || name in p.defs
}

fn (p &Parser) meaning(name string) string {
	what := reserved(p.ctx, name)
	return if what != '' { what } else { 'a definition' }
}

// read_unit reads a unit's text, from its place in the expression.
fn (p &Parser) read_unit(unit_text string, from int) !quantity.Unit {
	if unit_text == '' {
		return p.fail(from, 'empty brackets where a unit should be')
	}
	return quantity.unit(unit_text) or {
		if err is quantity.UnitError {
			e := err as quantity.UnitError
			return p.fail_with(from + e.column - 1, 'in the unit ${unit_text}: ${e.reason}', e.suggestions)
		}
		return p.fail(from, err.msg())
	}
}

fn has_offset(u quantity.Unit) bool {
	return u.pre != 0.0 || u.post != 0.0 || u.gauge
}

fn (p &Parser) offset_refused(unit_text string, from int) IError {
	return p.fail(from, "${unit_text} has an offset, so it can't stand in an expression: give the temperature in K or the pressure in Pa or bar, or read a variable on its scale, as (T in ${unit_text})")
}

// at_in is whether the word in starts here, alone.
fn (p &Parser) at_in() bool {
	return p.pos + 1 < p.runes.len && p.runes[p.pos] == `i` && p.runes[p.pos + 1] == `n`
		&& (p.pos + 2 >= p.runes.len || !is_name_char(p.runes[p.pos + 2]))
}

// unit_follows is whether a unit follows a space from k on.
fn (p &Parser) unit_follows(k int) bool {
	mut j := k
	for j < p.runes.len && is_space(p.runes[j]) {
		j++
	}
	return j > k && j < p.runes.len && (starts_unit(p.runes[j]) || p.runes[j] == `[`)
}

// bracketed reads a unit in brackets, the [ at p.pos: its text and where
// the text starts.
fn (mut p Parser) bracketed() !(string, int) {
	open := p.pos
	p.pos++
	from := p.pos
	for p.pos < p.runes.len && p.runes[p.pos] != `]` {
		p.pos++
	}
	if p.pos >= p.runes.len {
		return p.fail(open, 'a [ with no ] to close its unit')
	}
	text := p.runes[from..p.pos].string().trim_space()
	p.pos++
	return text, from
}

// reading reads " in unit" after a value, if it's there: the value's number
// on the unit's scale (D2, as amended before step 10). A unit with an
// offset reads only a context variable: T and P are absolute, where a
// difference, a member or a definition may not be.
fn (mut p Parser) reading(base int, start int) !int {
	save := p.pos
	p.skip()
	if p.pos == save || !p.at_in() {
		p.pos = save
		return base
	}
	at := p.pos
	if !p.unit_follows(p.pos + 2) {
		return p.fail(at, 'in needs a unit after it, as (T in degC)')
	}
	p.pos += 2
	p.skip()
	braced := p.peek() == `[`
	mut unit_text := ''
	mut from := p.pos
	if braced {
		unit_text, from = p.bracketed()!
	} else {
		unit_text = p.word()
	}
	u := p.read_unit(unit_text, from)!
	if has_offset(u) {
		what := match p.nodes[base].kind {
			.var { '' }
			.def { 'a definition, which may be a difference' }
			.member { 'a member of a family' }
			.num { 'a number' }
			else { 'an expression, which may be a difference' }
		}
		if what != '' {
			return p.fail(start, '${unit_text} has an offset, so it reads only a temperature or a pressure the context gives, such as T or P, not ${what}: read a difference in K or Pa')
		}
	}
	return p.add(Node{
		kind:   .reading
		a:      base
		value:  u.scale
		pre:    u.pre
		post:   u.post
		unit:   unit_text
		braced: braced
		dim:    u.dim
	}, start)
}

// unit_after reads a unit in brackets after a group's closing bracket, if
// one is there: (85 - 0.18*t) [%] is the group's value in that unit.
fn (mut p Parser) unit_after(inner int, open int) !int {
	save := p.pos
	p.skip()
	if p.peek() != `[` {
		p.pos = save
		return inner
	}
	unit_text, from := p.bracketed()!
	u := p.read_unit(unit_text, from)!
	if has_offset(u) {
		return p.offset_refused(unit_text, from)
	}
	return p.add(Node{
		kind:   .scaled
		a:      inner
		value:  u.scale
		unit:   unit_text
		braced: true
		dim:    u.dim
	}, open)
}

// with_unit makes a number with its unit, in SI.
fn (mut p Parser) with_unit(start int, text string, value f64, unit_text string, from int, braced bool) !int {
	u := p.read_unit(unit_text, from)!
	if has_offset(u) {
		return p.offset_refused(unit_text, from)
	}
	si := value * u.scale
	if !math.is_finite(si) {
		return p.fail(start, '${text} ${unit_text} is beyond the range of a number')
	}
	return p.add(Node{
		kind:   .num
		value:  si
		text:   text
		unit:   unit_text
		braced: braced
		dim:    u.dim
	}, start)
}

// name reads a name: a variable, a family member, a constant, a definition
// or a function's call.
fn (mut p Parser) name() !int {
	start := p.pos
	for p.pos < p.runes.len && is_name_char(p.runes[p.pos]) {
		p.pos++
	}
	name := p.runes[start..p.pos].string()
	after := p.pos
	p.skip()
	next := p.peek()
	if name == 'log' {
		return p.fail(start, 'log is ambiguous: write ln for the natural logarithm, or log10')
	}
	if name in functions {
		if next != `(` {
			return p.fail(start, '${name} needs its argument in parentheses: ${name}(…)')
		}
		return p.call(start, name)
	}
	if next == `(` {
		if p.known(name) {
			return p.fail(p.pos, "${name} isn't a function: write ${name}*(…) to multiply")
		}
		return p.fail_with(start, 'no function ${name}', quantity.closest(name, functions[..]))
	}
	for f, family in p.ctx.families {
		if family.name == name {
			if next != `[` {
				return p.fail(start, '${name} is a family of variables: say which, as ${name}[${family.members[0] or {
					'…'
				}}]')
			}
			return p.member(start, f)
		}
	}
	if next == `[` {
		families := p.ctx.families.map(it.name)
		if families.len == 0 {
			return p.fail(p.pos, "nothing here is a family of variables, so there's nothing to choose with [ ]")
		}
		return p.fail_with(start, "${name} isn't a family of variables: only ${families.join(', ')} take [ ]",
			quantity.closest(name, families))
	}
	p.pos = after
	for i, v in p.ctx.variables {
		if v.name == name {
			return p.add(Node{
				kind:  .var
				index: i
				text:  name
			}, start)
		}
	}
	for i, c in constants {
		if c.name == name {
			return p.add(Node{
				kind:  .constant
				index: i
				text:  name
			}, start)
		}
	}
	d := p.defs.index(name)
	if d >= 0 {
		return p.add(Node{
			kind:  .def
			index: d
			text:  name
		}, start)
	}
	mut every := []string{}
	every << p.ctx.variables.map(it.name)
	every << p.ctx.families.map(it.name)
	every << constants[..].map(it.name)
	every << p.defs
	every << functions[..]
	return p.fail_with(start, 'no name ${name} here', quantity.closest(name, every))
}

// call reads a function's argument, in parentheses.
fn (mut p Parser) call(start int, name string) !int {
	open := p.pos
	p.pos++
	p.skip()
	if p.peek() == `)` {
		return p.fail(open, '${name}() needs an argument')
	}
	p.depth++
	arg := p.expr()!
	p.depth--
	p.skip()
	if p.peek() == `,` {
		return p.fail(p.pos, '${name} takes one argument')
	}
	if p.peek() != `)` {
		if p.pos >= p.runes.len {
			return p.fail(open, 'a ( with no ) to close it')
		}
		return p.after_value()
	}
	p.pos++
	return p.add(Node{
		kind:  .call
		a:     arg
		index: functions[..].index(name)
		text:  name
	}, start)
}

// member reads a family member's key, in brackets.
fn (mut p Parser) member(start int, f int) !int {
	family := p.ctx.families[f]
	open := p.pos
	p.pos++
	p.skip()
	key_at := p.pos
	mut key := ''
	c := p.peek()
	if c == `'` || c == `"` {
		p.pos++
		from := p.pos
		for p.pos < p.runes.len && p.runes[p.pos] != c {
			p.pos++
		}
		if p.pos >= p.runes.len {
			return p.fail(key_at, 'a quote with no quote to close it')
		}
		key = p.runes[from..p.pos].string()
		p.pos++
	} else {
		for p.pos < p.runes.len && is_name_char(p.runes[p.pos]) {
			p.pos++
		}
		key = p.runes[key_at..p.pos].string()
	}
	p.skip()
	if p.peek() != `]` {
		if p.pos >= p.runes.len {
			return p.fail(open, 'a [ with no ] to close it')
		}
		if key == '' {
			return p.fail(p.pos, '${family.name}[ ] needs a name inside, as ${family.name}[${family.members[0] or {
				'…'
			}}]')
		}
		return p.fail(p.pos, "a key with ${p.peek()} in it must be quoted, as ${family.name}['${key}${p.peek()}…']")
	}
	p.pos++
	if key == '' {
		return p.fail(key_at, '${family.name}[ ] needs a name inside')
	}
	m := family.keys[key] or {
		return p.fail_with(key_at, 'no ${key} in ${family.name}[ ]', quantity.closest(key, family.keys.keys()))
	}
	return p.add(Node{
		kind:  .member
		a:     f
		b:     m
		index: p.ctx.slot(f, m)
		text:  family.members[m]
	}, start)
}
