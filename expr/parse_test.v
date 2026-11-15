// The expression language's syntax and names (M9 step 2, criterion 1):
// what reads, to which tree, and every refusal with its line, its column
// and a suggestion. Seeded random text never panics.
import expr { Kind }
import math
import quantity

// rate is a rate law's context: T and P, and each species' partial
// pressure, mole fraction and concentration, by formula, name or CAS
// number.
fn rate() expr.Context {
	species := ['CH4', 'H2O', 'H2', 'CO', 'CO2']
	names := ['Methane', 'Water', 'Hydrogen', 'Carbon monoxide', 'Carbon dioxide']
	cas := ['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9']
	mut keys := map[string]int{}
	for i, s in species {
		keys[s] = i
		keys[names[i]] = i
		keys[cas[i]] = i
	}
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, species, keys},
			expr.Family{'y', quantity.Dim{}, species, keys},
			expr.Family{'c', quantity.dim(0, -3, 0, 0, 1, 0), species, keys},
		]
	}
}

const no_defs = expr.Definitions{}

fn read(text string) !expr.Expr {
	return expr.parse(rate(), &no_defs, text)
}

// kinds lists a tree's kinds in node order, its operands before their use.
fn kinds(text string) []Kind {
	e := read(text) or { panic('${text}: ${err.msg()}') }
	return e.nodes.map(it.kind)
}

fn test_precedence_and_association() {
	// ^ binds tighter than a minus sign and to the right; · and / to the left.
	assert kinds('-T^2') == [Kind.var, .num, .pow, .neg]
	assert kinds('T^-2') == [Kind.var, .num, .neg, .pow]
	assert kinds('T^2^3') == [Kind.var, .num, .num, .pow, .pow]
	e := read('T^2^3')!
	assert e.nodes[e.root].b == 3 // the exponent is 2^3
	assert kinds('P/P/P') == [Kind.var, .var, .div, .var, .div]
	assert kinds('P/P*P') == [Kind.var, .var, .div, .var, .mul]
	assert kinds('P-P-P') == [Kind.var, .var, .sub, .var, .sub]
	assert kinds('-2*T') == [Kind.num, .neg, .var, .mul]
	assert kinds('+T') == [Kind.var]
	assert kinds('(T)') == [Kind.var]
	assert kinds('exp(T/T)') == [Kind.var, .var, .div, .call]
}

fn test_names_bind() {
	// A member by formula, name or CAS number is the same member.
	for key in ['CH4', "'Methane'", '"74-82-8"', "'CH4'"] {
		e := read('p[${key}]')!
		n := e.nodes[e.root]
		assert n.kind == .member && n.a == 0 && n.b == 0 && n.index == 2, key
		assert n.text == 'CH4', key // as the context names it
	}
	y := read('y[CO2]')!
	assert y.nodes[y.root].index == 2 + 5 + 4 // after T, P and p's five
	// Constants, by name.
	r := read('R')!
	assert r.nodes[r.root].kind == .constant && expr.constants[r.nodes[r.root].index].name == 'R'
	// Spaces anywhere between tokens, and a line break.
	assert read(' p [ CH4 ] ')!.nodes.len == 1
	assert kinds('T +\n  T') == [Kind.var, .var, .add]
}

fn test_numbers_and_units() {
	// Every way to write a number, and a unit after a space or in brackets.
	for text, want in {
		'2':               2.0
		'2.5':             2.5
		'.5':              0.5
		'4.22e15':         4.22e15
		'1E+3':            1e3
		'3.5e-2':          3.5e-2
		'240.1 kJ/mol':    240100.0
		'240.1 [kJ/mol]':  240100.0
		'1 bar':           1e5
		'6.65e-4 bar^-1':  6.65e-9
		'1 kmol/(kg*h)':   1e3 / 3600.0
		'2 [R]':           2.0 * 5.0 / 9.0
		'300 K':           300.0
		'1 kJ/(kg*°C)':    1e3
		'8.314 J/(mol*K)': 8.314
	} {
		e := read(text) or { panic('${text}: ${err.msg()}') }
		assert e.nodes.len == 1 && math.abs(e.nodes[0].value - want) <= 1e-15 * math.abs(want), '${text}: ${e.nodes[0].value}'
	}
	// A unit runs to the next space, or to a ) that isn't its own.
	assert kinds('exp(-240.1 kJ/mol / (R*T))') == [Kind.num, .neg, .constant, .var, .mul, .div,
		.call]
	assert kinds('(1 bar)') == [Kind.num]
	assert kinds('1 J/(mol*K) * T') == [Kind.num, .var, .mul]
	e := read('6.65e-4 bar^-1 * p[CH4]')!
	assert e.nodes[0].unit == 'bar^-1' && !e.nodes[0].braced
	assert read('2 [K]')!.nodes[0].braced
}

// refused asserts that text is refused at a line and column with a reason
// mentioning want, and returns the error.
fn refused(text string, line int, column int, want string) expr.ExprError {
	read(text) or {
		assert err is expr.ExprError, err.msg()
		e := err as expr.ExprError
		assert e.text == text, err.msg()
		assert e.reason.contains(want), '${text}: ${err.msg()}'
		assert e.line == line && e.column == column, '${text}: line ${e.line}, column ${e.column}, not ${line}, ${column}: ${err.msg()}'
		return e
	}
	assert false, '"${text}" was read'
	return expr.ExprError{}
}

fn test_every_syntax_error_has_its_place() {
	refused('', 1, 1, 'empty')
	refused('   ', 1, 1, 'empty')
	refused('T +', 1, 3, 'a + with nothing after it')
	refused('T + * T', 1, 5, 'a + with nothing after it, before *')
	refused('* T', 1, 1, 'a * with nothing before it')
	refused('T ^', 1, 3, 'a ^ with nothing after it')
	refused('-', 1, 1, 'a - with nothing after it')
	refused('(T + P', 1, 1, 'a ( with no ) to close it')
	refused('T + P)', 1, 6, 'a ) with no ( before it')
	refused('()', 1, 1, 'empty parentheses')
	refused('T T', 1, 3, 'two values side by side')
	refused('T (P)', 1, 3, "T isn't a function: write T*(…) to multiply")
	refused('2 (P)', 1, 3, 'two values side by side')
	refused('T = 2', 1, 3, 'an = here')
	refused('T $ P', 1, 3, 'unexpected $')
	refused('2T', 1, 2, 'run together: write 2*T')
	refused('2kg', 1, 2, 'put a space between a number and its unit: 2 kg')
	refused('2 T', 1, 3, 'a number and a name side by side: write 2*T')
	refused('2 exp(T)', 1, 3, 'write 2*exp')
	refused('2 p[CH4]', 1, 3, 'write 2*p')
	refused('2 R', 1, 3, 'R is a built-in constant, the molar gas constant, N_A·k_B and also a unit: write 2*R to multiply, or 2 [R] for the unit')
	refused('2 h', 1, 3, 'also a unit')
	refused('2 P', 1, 3, 'P is a variable and also a unit')
	refused('25 °C', 1, 4, '°C has an offset')
	refused('2 barg', 1, 3, 'barg has an offset')
	refused('2 [kmol/kg*h]', 1, 11, 'in the unit kmol/kg*h: a * after a / is ambiguous')
	refused('2 [kmol/[kg.h]', 1, 9, 'group units with ( and ), not [ and ]')
	refused('2 [', 1, 3, 'a [ with no ] to close its unit')
	refused('2 []', 1, 4, 'empty brackets')
	refused('1e999', 1, 1, 'beyond the range of a number')
	refused('1e300 GPa', 1, 1, 'beyond the range of a number')
	// 200 levels of nesting read; the 201st is refused, so that no text can
	// exhaust the stack.
	assert read('('.repeat(150) + 'T' + ')'.repeat(150))!.nodes.len == 1
	refused('('.repeat(300) + 'T' + ')'.repeat(300), 1, 201, 'nested too deeply')
	refused('-'.repeat(300) + 'T', 1, 201, 'nested too deeply')
	// A line break counts lines and resets the column.
	refused('T +\n  * P', 2, 3, 'a + with nothing after it, before *')
	refused('T\n+ (P', 2, 3, 'a ( with no ) to close it')
}

fn test_every_name_error_has_a_suggestion() {
	assert 'T' in refused('Tt', 1, 1, 'no name Tt here').suggestions
	assert 'R' in refused('r', 1, 1, 'no name r here').suggestions
	assert 'exp' in refused('expp(T)', 1, 1, 'no function expp').suggestions
	refused('log(T)', 1, 1, 'log is ambiguous: write ln for the natural logarithm, or log10')
	refused('log', 1, 1, 'log is ambiguous')
	refused('exp T', 1, 1, 'exp needs its argument in parentheses')
	refused('exp()', 1, 4, 'exp() needs an argument')
	refused('exp(T, P)', 1, 6, 'exp takes one argument')
	// A number's unit ends at a comma, which is then the call's.
	refused('sqrt(2 bar, P)', 1, 11, 'sqrt takes one argument')
	refused('T(P)', 1, 2, "T isn't a function: write T*(…) to multiply")
	refused('p', 1, 1, 'p is a family of variables: say which, as p[CH4]')
	assert 'p' in refused('q[CH4]', 1, 1, "q isn't a family of variables: only p, y, c take [ ]").suggestions
	assert 'CH4' in refused('p[CH5]', 1, 3, 'no CH5 in p[ ]').suggestions
	assert 'Methane' in refused("p['methane']", 1, 3, 'no methane in p[ ]').suggestions
	refused('p[]', 1, 3, 'p[ ] needs a name inside')
	refused('p[n-butane]', 1, 4, "a key with - in it must be quoted, as p['n-…']")
	refused('p[CH4', 1, 2, 'a [ with no ] to close it')
	refused("p['CH4]", 1, 3, 'a quote with no quote to close it')
}

fn test_eval_reads_the_inputs() {
	e := read('p[CH4] * y[H2O] / T * 300 K + 2 bar')!
	mut x := []f64{len: rate().slots()}
	x[0] = 500.0 // T
	x[2] = 3e5 // p[CH4]
	x[2 + 5 + 1] = 0.25 // y[H2O]
	assert e.eval(&no_defs, x) == 3e5 * 0.25 / 500.0 * 300.0 + 2e5
}

// Rng is a small seeded generator for the fuzzing: splitmix64.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) below(n int) int {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return int((z ^ (z >> 31)) % u64(n))
}

// same_tree is whether two trees are the same, node for node, but for
// where each was written.
fn same_tree(a expr.Expr, b expr.Expr) bool {
	if a.root != b.root || a.nodes.len != b.nodes.len {
		return false
	}
	for i, x in a.nodes {
		y := b.nodes[i]
		if x.kind != y.kind || x.a != y.a || x.b != y.b || x.index != y.index
			|| math.f64_bits(x.value) != math.f64_bits(y.value) || x.text != y.text
			|| x.unit != y.unit || x.braced != y.braced || !x.dim.eq(y.dim) {
			return false
		}
	}
	return true
}

// round_trips asserts that an expression's plain form reads back to the
// same tree, and that its other renderings are made.
fn round_trips(t string, e expr.Expr) {
	back := read(e.plain()) or { panic('"${t}" wrote "${e.plain()}", which doesn\'t read: ${err.msg()}') }
	assert same_tree(e, back), '"${t}" wrote "${e.plain()}", a different tree'
	assert e.unicode() != '' && e.latex() != ''
}

fn test_random_text_never_panics() {
	// 20,000 seeded strings of the language's pieces, and of stray ones.
	// Each is refused with an ExprError with a place and a reason, or reads
	// to a tree that round-trips.
	pieces := ['T', 'P', 'R', 'p[CH4]', 'y[H2O]', "p['Water']", 'p[', ']', '[K]', ' bar', ' kJ/mol',
		'2', '0.5', '1e3', '.', 'e', '(', ')', '+', '-', '*', '/', '^', '·', '−', 'exp(', 'ln(',
		'sqrt(', 'abs(', 'log10(', ',', ' ', '\n', 'Ω', '°', 'pi', 'k_B', 'x', '=', "'", '1/(mol*K)',
		'!']
	mut rng := Rng{20261110}
	mut accepted := 0
	for _ in 0 .. 20000 {
		mut t := ''
		for _ in 0 .. 1 + rng.below(10) {
			t += pieces[rng.below(pieces.len)]
		}
		e := read(t) or {
			assert err is expr.ExprError, '"${t}": ${err.msg()}'
			e := err as expr.ExprError
			assert e.line >= 1 && e.column >= 1 && e.reason != '', '"${t}": ${err.msg()}'
			continue
		}
		accepted++
		round_trips(t, e)
	}
	assert accepted > 100
}

// random_text is a seeded, well-formed expression of pure numbers, so that
// every one reads: a term with units appears only as a ratio to its own
// unit, as (T/(300 K)). Spacing and parentheses vary.
fn random_text(mut rng Rng, depth int) string {
	leaves := ['y[H2O]', "y['Methane']", '2', '0.5', '1e3', '3.5e-2', 'pi', '(T / (300 K))',
		'(p[CH4]/(1 bar))', '(P / p_std)', '(240.1 kJ/mol / (R*T))', '(c[CO] * 1 [m^3/mol])']
	if depth == 0 || rng.below(4) == 0 {
		return leaves[rng.below(leaves.len)]
	}
	sp := [' ', '', '  '][rng.below(3)]
	a := random_text(mut rng, depth - 1)
	return match rng.below(9) {
		0 { '${a}${sp}+${sp}${random_text(mut rng, depth - 1)}' }
		1 { '${a}${sp}-${sp}${random_text(mut rng, depth - 1)}' }
		2 { '${a}${sp}*${sp}${random_text(mut rng, depth - 1)}' }
		3 { '${a}${sp}/${sp}(${random_text(mut rng, depth - 1)})' }
		4 { '(${a})^${['2', '0.5', '-1', '(1/3)', 'y[H2]'][rng.below(5)]}' }
		5 { '-${sp}(${a})' }
		6 { '${['exp', 'ln', 'sqrt', 'abs', 'log10'][rng.below(5)]}(${a})' }
		7 { '(${a})' }
		else { '${a}${sp}·${sp}${random_text(mut rng, depth - 1)}' }
	}
}

fn test_random_expressions_read_and_round_trip() {
	// 10,000 seeded well-formed expressions: each reads, is a pure number,
	// and round-trips.
	mut rng := Rng{20261111}
	for _ in 0 .. 10000 {
		t := random_text(mut rng, 5)
		e := read(t) or { panic('"${t}" was refused: ${err.msg()}') }
		assert e.dim.is_dimensionless(), '"${t}": ${e.dim}'
		round_trips(t, e)
	}
}
