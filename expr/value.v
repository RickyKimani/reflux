module expr

import math

// value is the expression's value at x, as eval gives it, or why it has none
// (M8's D9). Where eval would give ±∞ or NaN, value refuses, at the node
// where it happened and in words:
// - an input with no value, NaN, and why[slot] if it's given;
// - a division by zero, which is how a product's conversion fails when the
//   feed has none of it;
// - a logarithm of a number not above zero, or a square root of a negative
//   one;
// - zero to a negative power, or a negative number to a power that isn't
//   whole;
// - a result too large for a double.
//
// It does eval's arithmetic in eval's order, so where there is a value it's
// eval's, to the bit. The failure is an ExprError at the node's line and
// column, saying which part of the expression had no value.
pub fn (e &Expr) value(defs &Definitions, x []f64, why []string) !f64 {
	return e.checked(defs, x, why, e.root)
}

fn (e &Expr) checked(defs &Definitions, x []f64, why []string, i int) !f64 {
	n := e.nodes[i]
	match n.kind {
		.num {
			return n.value
		}
		.var, .member {
			v := x[n.index]
			if !math.is_finite(v) {
				reason := why[n.index] or { '' }
				tail := if reason != '' { ': ${reason}' } else { '' }
				return e.no_value(n, '${e.form(i, .plain)} has no value here${tail}')
			}
			return v
		}
		.constant {
			return constants[n.index].value
		}
		.def {
			d := defs.exprs[n.index]
			return d.checked(defs, x, why, d.root) or {
				return e.no_value(n, '${n.text} has no value here: ${reason_in(err)}')
			}
		}
		.neg {
			return -e.checked(defs, x, why, n.a)!
		}
		.add, .sub, .mul {
			a := e.checked(defs, x, why, n.a)!
			b := e.checked(defs, x, why, n.b)!
			r := match n.kind {
				.add { a + b }
				.sub { a - b }
				else { a * b }
			}
			if !math.is_finite(r) {
				return e.no_value(n, '${e.form(i, .plain)} is too large for a number here')
			}
			return r
		}
		.div {
			a := e.checked(defs, x, why, n.a)!
			b := e.checked(defs, x, why, n.b)!
			if b == 0.0 {
				return e.no_value(n, 'it divides by ${e.form(n.b, .plain)}, which is zero here')
			}
			r := a / b
			if !math.is_finite(r) {
				return e.no_value(n, '${e.form(i, .plain)} is too large for a number here')
			}
			return r
		}
		.pow {
			a := e.checked(defs, x, why, n.a)!
			b := e.checked(defs, x, why, n.b)!
			r := raise(a, b)
			if !math.is_finite(r) {
				if a == 0.0 && b < 0.0 {
					return e.no_value(n, '${e.form(n.a, .plain)} is zero here, and its power ${e.form(n.b, .plain)} is negative')
				}
				if a < 0.0 && b != math.trunc(b) {
					return e.no_value(n, '${e.form(n.a, .plain)} is ${number_text(a)} here, and its power ${e.form(n.b, .plain)} is ${number_text(b)}, not a whole number')
				}
				return e.no_value(n, '${e.form(i, .plain)} is too large for a number here')
			}
			return r
		}
		.call {
			a := e.checked(defs, x, why, n.a)!
			r := apply(n.index, a)
			if !math.is_finite(r) {
				arg := e.form(n.a, .plain)
				return match functions[n.index] {
					'ln', 'log10' {
						e.no_value(n, 'it takes ${n.text} of ${arg}, which is ${number_text(a)} here: a logarithm needs a number above zero')
					}
					'sqrt' {
						e.no_value(n, 'it takes the square root of ${arg}, which is ${number_text(a)} here')
					}
					else {
						e.no_value(n, '${n.text}(${arg}) is too large for a number here: its argument is ${number_text(a)}')
					}
				}
			}
			return r
		}
		.reading {
			r := read_on(n, e.checked(defs, x, why, n.a)!)
			if !math.is_finite(r) {
				return e.no_value(n, '${e.form(i, .plain)} is too large for a number here')
			}
			return r
		}
		.scaled {
			r := e.checked(defs, x, why, n.a)! * n.value
			if !math.is_finite(r) {
				return e.no_value(n, '${e.form(i, .plain)} is too large for a number here')
			}
			return r
		}
	}
}

fn (e &Expr) no_value(n Node, reason string) IError {
	return ExprError{
		text:   e.text
		line:   n.line
		column: n.column
		reason: reason
	}
}

// number_text is a number for a message: a whole one without its .0.
fn number_text(v f64) string {
	if v == math.trunc(v) && math.abs(v) < 1e15 {
		return i64(v).str()
	}
	return v.str()
}

// reason_in is an expression's failure without its place.
fn reason_in(err IError) string {
	if err is ExprError {
		return err.reason
	}
	return err.msg()
}
