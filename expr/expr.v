// Package expr is reflux's expression language (M9, decisions D1, D2, D4
// and D5): every equation a user gives, a rate law, a rate constant, an
// adsorption term, a pressure drop, is parsed here into a tree, never
// evaluated as text.
//
// - Every name is bound: to a variable the caller's Context has (T, P, a
//   species' partial pressure p[CH4]), a definition, or a built-in constant
//   at the 2019 SI's exact value (R, not a typed 8.314). A name that isn't
//   there is refused with suggestions.
// - Every number may carry a unit, "240.1 kJ/mol", converted to SI as it's
//   read, and every node has a dimension, checked as the tree is built: an
//   exp of something with units, or a pressure added to a temperature, is
//   refused with its reason.
// - Every tree renders back as text that reads back to the same tree, and
//   as Unicode and LaTeX, to hold against the paper.
//
// A tree is an arena of nodes, each referring to its operands by index, in
// an order where every operand comes before its use: data, not closures,
// as step 3's tape will be. expr imports quantity and nothing else of
// reflux's (D1).
module expr

import math
import quantity { Dim }

// Kind is what a node is.
pub enum Kind as u8 {
	num      // a number, with its unit: value is in SI
	var      // a context variable: index is its slot
	member   // a family member, p[CH4]: a is the family, b the member, index the slot
	constant // a built-in constant: index into constants
	def      // a definition: index into the Definitions
	neg      // −a
	add      // a + b
	sub      // a − b
	mul      // a·b
	div      // a/b
	pow      // a^b
	call     // a function of a: index into functions
}

// Node is one node of a tree. line and column are where it starts in the
// text, from 1, in characters.
pub struct Node {
pub:
	kind   Kind
	a      int = -1
	b      int = -1
	index  int
	value  f64
	text   string // a number or a name, as written
	unit   string // a number's unit, as written
	braced bool   // the unit was written in brackets, 2 [K]
	line   int
	column int
pub mut:
	dim Dim
}

// Expr is a parsed, bound and checked expression: its nodes, its root, and
// the text it came from. dim is the root's dimension, and families the
// context's family names, for writing a member back.
pub struct Expr {
pub:
	text     string
	nodes    []Node
	root     int
	dim      Dim
	families []string
}

// Variable is a scalar a context gives, such as T or P.
pub struct Variable {
pub:
	name string
	dim  Dim
}

// Family is a set of variables named by a key, such as each species'
// partial pressure, p[CH4]. members are each member's name as rendered (a
// formula, say), and keys every name a member is known by (its name, its
// formula, its CAS number), each giving the member's number.
pub struct Family {
pub:
	name    string
	dim     Dim
	members []string
	keys    map[string]int
}

// Context is what an expression's names may be bound to: its variables
// and families. Their slots, an evaluation's inputs, are the variables in
// order, then each family's members in order.
pub struct Context {
pub:
	variables []Variable
	families  []Family
}

// slots is how many inputs an evaluation takes.
pub fn (c Context) slots() int {
	mut n := c.variables.len
	for f in c.families {
		n += f.members.len
	}
	return n
}

// slot is the input of family f's member m.
pub fn (c Context) slot(f int, m int) int {
	mut n := c.variables.len
	for i in 0 .. f {
		n += c.families[i].members.len
	}
	return n + m
}

// Constant is a built-in constant, with its unit as written for a reader.
pub struct Constant {
pub:
	name  string
	value f64
	unit  string
	dim   Dim
	about string
}

// constants are the built-in constants (D2): the 2019 SI's exact values,
// and standard gravity and the standard pressure, exact by definition. R is
// N_A·k_B and F is N_A·e, each the double nearest the exact product.
pub const constants = [
	Constant{'R', 8.31446261815324, 'J/(mol*K)', quantity.dim(1, 2, -2, -1, -1, 0), 'the molar gas constant, N_A·k_B'},
	Constant{'N_A', 6.02214076e23, '1/mol', quantity.dim(0, 0, 0, 0, -1, 0), 'the Avogadro constant'},
	Constant{'k_B', 1.380649e-23, 'J/K', quantity.dim(1, 2, -2, -1, 0, 0), 'the Boltzmann constant'},
	Constant{'h', 6.62607015e-34, 'J*s', quantity.dim(1, 2, -1, 0, 0, 0), 'the Planck constant'},
	Constant{'F', 96485.33212331001, 'C/mol', quantity.dim(0, 0, 1, 0, -1, 1), 'the Faraday constant, N_A·e'},
	Constant{'c', 299792458.0, 'm/s', quantity.dim(0, 1, -1, 0, 0, 0), 'the speed of light in vacuum'},
	Constant{'g_n', 9.80665, 'm/s^2', quantity.dim(0, 1, -2, 0, 0, 0), 'standard gravity'},
	Constant{'pi', math.pi, '', quantity.Dim{}, 'π'},
	Constant{'p_std', 1e5, 'Pa', quantity.dim(1, -1, -2, 0, 0, 0), 'the standard pressure, 1 bar'},
]!

// functions are the functions an expression may call, each of one
// argument. exp, ln and log10 take a pure number; sqrt halves its
// argument's dimension and abs keeps it.
pub const functions = ['exp', 'ln', 'log10', 'sqrt', 'abs']!

// ExprError is every failure here: the text, the line and column (from 1,
// in characters) where it went wrong, why, and any names it may have
// meant. definition is the definition the text is, if it is one.
pub struct ExprError {
	Error
pub:
	text        string
	line        int
	column      int
	reason      string
	suggestions []string
	definition  string
}

// msg is the reason, with where it is and any suggestions.
pub fn (e ExprError) msg() string {
	hint := if e.suggestions.len > 0 { ' (did you mean ${e.suggestions.join(', ')}?)' } else { '' }
	what := if e.definition != '' { '${e.definition}: ' } else { '' }
	return '${what}line ${e.line}, column ${e.column}: ${e.reason}${hint}'
}

// Source is a definition as given: its name and its text.
pub struct Source {
pub:
	name string
	text string
}

// Definitions are named expressions, bound and checked, which other
// expressions may use by name. order lists them so that each comes after
// every definition it uses.
pub struct Definitions {
pub:
	names []string
	exprs []Expr
	order []int
}

// find is the definition of this name, or -1.
pub fn (d &Definitions) find(name string) int {
	for i, n in d.names {
		if n == name {
			return i
		}
	}
	return -1
}

// reserved is what a name means already, if anything: a function, a
// constant, a variable or a family.
fn reserved(ctx Context, name string) string {
	if name in functions {
		return 'a function'
	}
	for c in constants {
		if c.name == name {
			return 'a built-in constant, ${c.about}'
		}
	}
	for v in ctx.variables {
		if v.name == name {
			return 'a variable'
		}
	}
	for f in ctx.families {
		if f.name == name {
			return 'a family of variables, as in ${f.name}[${f.members[0] or { '…' }}]'
		}
	}
	return ''
}

// define parses, binds and checks a set of definitions. Each may use the
// others, in any order, but not itself through any chain of them; and none
// may take the name of a function, a constant or one of ctx's variables.
pub fn define(ctx Context, sources []Source) !Definitions {
	mut names := []string{}
	for s in sources {
		what := reserved(ctx, s.name)
		if what != '' {
			return ExprError{
				text:       s.text
				line:       1
				column:     1
				reason:     "${s.name} is ${what}, so it can't name a definition"
				definition: s.name
			}
		}
		if !is_identifier(s.name) {
			return ExprError{
				text:       s.text
				line:       1
				column:     1
				reason:     '"${s.name}" can\'t name a definition: a name is letters, digits and _, starting with a letter'
				definition: s.name
			}
		}
		if s.name in names {
			return ExprError{
				text:       s.text
				line:       1
				column:     1
				reason:     '${s.name} is defined twice'
				definition: s.name
			}
		}
		names << s.name
	}
	// Parse each, with every definition's name known, then check them in an
	// order where each comes after what it uses, refusing a cycle.
	mut trees := []Expr{}
	for s in sources {
		mut p := new_parser(ctx, names, s.text)
		trees << p.parse() or { return in_definition(err, s.name) }
	}
	order := ordered(names, trees, sources)!
	mut checked := []Expr{len: trees.len}
	for i in order {
		partial := Definitions{
			names: names
			exprs: checked
		}
		checked[i] = check(ctx, &partial, trees[i]) or { return in_definition(err, names[i]) }
	}
	return Definitions{
		names: names
		exprs: checked
		order: order
	}
}

fn in_definition(err IError, name string) IError {
	if err is ExprError {
		e := err as ExprError
		return ExprError{
			...e
			definition: name
		}
	}
	return err
}

// ordered lists the definitions so that each follows those it uses, by a
// depth-first walk, or refuses the first cycle found, saying its chain.
fn ordered(names []string, trees []Expr, sources []Source) ![]int {
	mut state := []int{len: names.len} // 0 not seen, 1 on the path, 2 done
	mut order := []int{}
	mut path := []int{}
	for i in 0 .. names.len {
		visit(i, names, trees, sources, mut state, mut order, mut path)!
	}
	return order
}

fn visit(i int, names []string, trees []Expr, sources []Source, mut state []int, mut order []int, mut path []int) ! {
	if state[i] == 2 {
		return
	}
	if state[i] == 1 {
		start := path.index(i)
		chain := path[start..].map(names[it])
		next := path[(start + 1) % path.len] or { i }
		at := first_use(trees[i], next)
		return ExprError{
			text:       sources[i].text
			line:       at.line
			column:     at.column
			reason:     '${names[i]} uses itself: ${chain.join(' → ')} → ${names[i]}'
			definition: names[i]
		}
	}
	state[i] = 1
	path << i
	for n in trees[i].nodes {
		if n.kind == .def {
			visit(n.index, names, trees, sources, mut state, mut order, mut path)!
		}
	}
	path.pop()
	state[i] = 2
	order << i
}

// first_use is the first node of e that uses definition d, or the root.
fn first_use(e Expr, d int) Node {
	for n in e.nodes {
		if n.kind == .def && n.index == d {
			return n
		}
	}
	return e.nodes[e.root]
}

// parse reads an expression in a context, with these definitions, binds
// its names and checks its units.
pub fn parse(ctx Context, defs &Definitions, text string) !Expr {
	mut p := new_parser(ctx, defs.names, text)
	tree := p.parse()!
	return check(ctx, defs, tree)
}

// must_be refuses an expression whose dimension isn't d: what is what the
// expression is for, as a message says it ("a rate per mass of catalyst").
pub fn (e Expr) must_be(d Dim, what string) ! {
	if e.dim.eq(d) {
		return
	}
	root := e.nodes[e.root]
	return ExprError{
		text:   e.text
		line:   root.line
		column: root.column
		reason: '${what} must be a ${d.name()}, but this is a ${e.dim.name()}'
	}
}

// uses lists the definitions an expression uses directly.
pub fn (e Expr) uses() []int {
	mut out := []int{}
	for n in e.nodes {
		if n.kind == .def && n.index !in out {
			out << n.index
		}
	}
	return out
}

fn is_identifier(s string) bool {
	if s == '' || !(s[0].is_letter()) {
		return false
	}
	for c in s {
		if !(c.is_letter() || c.is_digit() || c == `_`) {
			return false
		}
	}
	return true
}
