module flowsheet

import expr
import math
import quantity
import toml.ast

// Section is one table of a flowsheet file as it's read. Reading a key marks
// it used, finish refuses any key left over, and every failure names its key
// and line. It's what a Maker reads a unit's specification from.
pub struct Section {
pub:
	path string // 'unit.cooler'
	text string // the whole file
pub mut:
	table      map[string]ast.Value
	name       string   // the unit's name, for a unit
	components []string // the file's, for what goes by component
	used       []string
	context    Context // the file's species and reactions, for a unit's maker
}

// Sized is a quantity read as one of two dimensions.
pub struct Sized {
pub:
	value f64
	dim   quantity.Dimension
}

// has is whether a key is there. It doesn't count as reading it.
pub fn (s Section) has(key string) bool {
	return key in s.table
}

fn (s Section) full(key string) string {
	return if s.path == '' {
		key
	} else if key == '' {
		s.path
	} else {
		'${s.path}.${key}'
	}
}

// line_of is the line a value starts on, or 0 for an empty table or array.
fn line_of(v ast.Value) int {
	return match v {
		ast.Quoted, ast.Number, ast.Bool, ast.Date, ast.Time, ast.DateTime, ast.Null {
			v.pos.line_nr
		}
		[]ast.Value {
			if v.len > 0 { line_of(v[0]) } else { 0 }
		}
		map[string]ast.Value {
			mut line := 0
			for _, x in v {
				line = line_of(x)
				if line > 0 {
					break
				}
			}
			line
		}
	}
}

// header_line is the line of a table's header, if it has one of its own.
fn header_line(text string, path string) int {
	for i, l in text.split_into_lines() {
		t := l.trim_space()
		if t == '[${path}]' || t.starts_with('[${path}] ') {
			return i + 1
		}
	}
	return 0
}

// line is the line to report for a key, or for the table with ''.
fn (s Section) line(key string) int {
	if v := s.table[key] {
		line := line_of(v)
		if line > 0 {
			return line
		}
	}
	return s.line_of_root(s.path)
}

// line_of_root is the line of the table at a path. The file as a whole has
// none.
fn (s Section) line_of_root(path string) int {
	if path == '' {
		return 0
	}
	h := header_line(s.text, path)
	return if h > 0 { h } else { line_of(ast.Value(s.table)) }
}

// fail is a failure at a key, with any suggestions.
pub fn (s Section) fail(key string, reason string, suggestions []string) IError {
	return file_error(s.line(key), s.full(key), reason, suggestions)
}

// located gives an error this section's place, if it has none.
fn (s Section) located(err IError) IError {
	if err is FlowsheetError {
		if err.kind == .bad_file {
			return err
		}
		return FlowsheetError{
			kind:   err.kind
			unit:   err.unit
			reason: err.reason
			cause:  err.cause
			line:   s.line('')
			key:    s.path
		}
	}
	return s.fail('', err.msg(), [])
}

// require is a key's value, marked used.
fn (mut s Section) require(key string) !ast.Value {
	s.used << key
	return s.table[key] or { return s.fail('', '${key} is missing', []) }
}

// text is a string.
pub fn (mut s Section) text(key string) !string {
	v := s.require(key)!
	if v is ast.Quoted {
		return v.text
	}
	return s.fail(key, 'it should be a string', [])
}

// boolean is true or false.
pub fn (mut s Section) boolean(key string) !bool {
	v := s.require(key)!
	if v is ast.Bool {
		return v.text == 'true'
	}
	return s.fail(key, 'it should be true or false', [])
}

// number is a TOML number.
pub fn (mut s Section) number(key string) !f64 {
	v := s.require(key)!
	if v is ast.Number {
		return v.f64()
	}
	return s.fail(key, 'it should be a number', [])
}

// integer is a whole number.
pub fn (mut s Section) integer(key string) !i64 {
	v := s.require(key)!
	if v is ast.Number {
		if v.f64() == f64(v.i64()) {
			return v.i64()
		}
	}
	return s.fail(key, 'it should be a whole number', [])
}

// quantity is a number with its unit (D10), in SI.
pub fn (mut s Section) quantity(key string, d quantity.Dimension) !f64 {
	v := s.require(key)!
	return s.quantity_of(key, v, d)
}

fn (s Section) quantity_of(key string, v ast.Value, d quantity.Dimension) !f64 {
	if v is ast.Quoted {
		return quantity.parse(v.text, d) or {
			if err is quantity.QuantityError {
				return file_error(v.pos.line_nr, s.full(key), err.reason, err.suggestions)
			}
			return s.fail(key, err.msg(), [])
		}
	}
	return s.fail(key, 'it needs a unit, as a string, such as "${quantity.units(d)[0]}"', [])
}

// measure is a number with a unit of any dimension, read as an expression
// reads one (M9): a volume, a length or a density, which M4's named
// dimensions haven't. what is what it's for, as a message says it.
pub fn (mut s Section) measure(key string, d quantity.Dim, what string) !f64 {
	text := s.text(key)!
	defs := expr.Definitions{}
	e := expr.parse(expr.Context{}, &defs, text) or { return s.fail(key, reason_of(err), []) }
	e.must_be(d, what) or { return s.fail(key, reason_of(err), []) }
	v := e.eval(&defs, []f64{})
	if !math.is_finite(v) {
		return s.fail(key, '${text} is ${v}', [])
	}
	return v
}

// expression reads an expression in a context, with its definitions (M9
// step 10): a number, or a string the expression language reads, such as
// "0.85 - 0.0018*(T in degC)", of dimension d. what is what it's for, as a
// message says it.
pub fn (mut s Section) expression(key string, ctx expr.Context, defs &expr.Definitions, d quantity.Dim, what string) !expr.Expr {
	v := s.require(key)!
	text := match v {
		ast.Quoted { v.text }
		ast.Number { v.text }
		else { return s.fail(key, 'it should be a number, or an expression in quotes', []) }
	}
	e := expr.parse(ctx, defs, text) or { return s.fail(key, reason_of(err), []) }
	e.must_be(d, what) or { return s.fail(key, reason_of(err), []) }
	return e
}

// reason_of is an expression's failure without its place in the string,
// which the key's line already gives.
fn reason_of(err IError) string {
	if err is expr.ExprError {
		return err.reason
	}
	return err.msg()
}

// quantity_either is a number with a unit of either dimension, such as a flow
// by moles or by mass.
pub fn (mut s Section) quantity_either(key string, a quantity.Dimension, b quantity.Dimension) !Sized {
	v := s.require(key)!
	if x := s.quantity_of(key, v, a) {
		return Sized{x, a}
	}
	x := s.quantity_of(key, v, b)!
	return Sized{x, b}
}

// fraction is a plain number or a percentage.
pub fn (mut s Section) fraction(key string) !f64 {
	v := s.require(key)!
	return s.fraction_of(key, v)
}

fn (s Section) fraction_of(key string, v ast.Value) !f64 {
	if v is ast.Number {
		return v.f64()
	}
	return s.quantity_of(key, v, .fraction)
}

// fractions is a list of fractions.
pub fn (mut s Section) fractions(key string) ![]f64 {
	v := s.require(key)!
	if v is []ast.Value {
		mut out := []f64{cap: v.len}
		for x in v {
			out << s.fraction_of(key, x)!
		}
		return out
	}
	return s.fail(key, 'it should be a list of fractions', [])
}

// strings is a list of strings. A single string is a list of one.
pub fn (mut s Section) strings(key string) ![]string {
	v := s.require(key)!
	if v is ast.Quoted {
		return [v.text]
	}
	if v is []ast.Value {
		mut out := []string{cap: v.len}
		for x in v {
			if x is ast.Quoted {
				out << x.text
			} else {
				return s.fail(key, 'it should be a list of strings', [])
			}
		}
		return out
	}
	return s.fail(key, 'it should be a list of strings', [])
}

// one_of is which of these keys is given: exactly one must be (D3).
pub fn (s Section) one_of(keys []string) !string {
	given := keys.filter(s.has(it))
	if given.len == 0 {
		return s.fail('', 'a specification is missing: give one of ${keys.join(', ')}', [])
	}
	if given.len > 1 {
		return s.fail(given[1], 'give only one of ${keys.join(', ')}, not ${given.join(' and ')}',
			[])
	}
	return given[0]
}

// child is a table within this one, as a section of its own.
pub fn (mut s Section) child(key string) !Section {
	v := s.require(key)!
	if v is map[string]ast.Value {
		return Section{
			path:  s.full(key)
			text:  s.text
			table: v
		}
	}
	return s.fail(key, 'it should be a table', [])
}

// tables is a list of tables, as sections.
pub fn (mut s Section) tables(key string) ![]Section {
	v := s.require(key)!
	if v is []ast.Value {
		mut out := []Section{cap: v.len}
		for i, x in v {
			if x is map[string]ast.Value {
				out << Section{
					path:  '${s.full(key)}[${i}]'
					text:  s.text
					table: x
				}
			} else {
				return s.fail(key, 'it should be a list of tables', [])
			}
		}
		return out
	}
	return s.fail(key, 'it should be a list of tables', [])
}

// by_component reads a table with a value per component, in the components'
// order. A component left out is 0. Fractions may be plain numbers.
pub fn (mut s Section) by_component(key string, d quantity.Dimension) ![]f64 {
	v := s.require(key)!
	if v is map[string]ast.Value {
		mut out := []f64{len: s.components.len}
		for name, x in v {
			i := s.components.index(name)
			if i < 0 {
				return file_error(line_of(x), s.full(key), '${name} is not one of the components',
					quantity.closest(name, s.components))
			}
			out[i] = if d == .fraction {
				s.fraction_of(key, x)!
			} else {
				s.quantity_of(key, x, d)!
			}
		}
		return out
	}
	return s.fail(key, 'it should be a table with a value for each component', [])
}

// ports maps a unit's ports to the streams its table names: `in` and `out`
// list them in port order, and a key named for a port gives that port's.
fn (mut s Section) ports(u Unit) !map[string]string {
	mut out := map[string]string{}
	s.bind(mut out, 'in', u.inlet_ports())!
	s.bind(mut out, 'out', u.outlet_ports())!
	return out
}

fn (mut s Section) bind(mut out map[string]string, key string, ports []string) ! {
	if s.has(key) {
		streams := s.strings(key)!
		if streams.len != ports.len {
			return s.fail(key, 'it needs ${ports.len} streams, for ${ports.join(', ')}, not ${streams.len}',
				[])
		}
		for i, port in ports {
			out[port] = streams[i]
		}
		return
	}
	for port in ports {
		if !s.has(port) {
			return s.fail('', '${key} is missing', [])
		}
		out[port] = s.text(port)!
	}
}

// finish refuses a key nothing read (D9). A misspelt key is suggested the keys
// that were read.
pub fn (s Section) finish() ! {
	for key, _ in s.table {
		if key !in s.used {
			return s.fail(key, 'no key ${key} here', quantity.closest(key, s.used))
		}
	}
}
