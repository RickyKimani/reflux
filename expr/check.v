module expr

import math
import quantity { Dim, Ratio }

// check gives every node of a tree its dimension, refusing what has none:
// a sum of two dimensions, an exp of a quantity with units, a power with
// units, or a quantity with units raised to a power that isn't a fixed
// fraction. Every operand comes before its use, so one pass in order does.
fn check(ctx Context, defs &Definitions, tree Expr) !Expr {
	mut nodes := tree.nodes.clone()
	for i in 0 .. nodes.len {
		n := nodes[i]
		nodes[i].dim = match n.kind {
			.num {
				n.dim
			}
			.var {
				ctx.variables[n.index].dim
			}
			.member {
				ctx.families[n.a].dim
			}
			.constant {
				constants[n.index].dim
			}
			.def {
				defs.exprs[n.index].dim
			}
			.neg {
				nodes[n.a].dim
			}
			.add, .sub {
				a, b := nodes[n.a].dim, nodes[n.b].dim
				if !a.eq(b) {
					op := if n.kind == .add { '+' } else { '−' }
					return fail_at(tree, n, 'the two sides of ${op} must have the same units, but the left is ${describe(a)} and the right ${describe(b)}')
				}
				a
			}
			.mul {
				nodes[n.a].dim.mul(nodes[n.b].dim) or { return fail_at(tree, n, err.msg()) }
			}
			.div {
				nodes[n.a].dim.div(nodes[n.b].dim) or { return fail_at(tree, n, err.msg()) }
			}
			.pow {
				power_dim(defs, nodes, tree, n)!
			}
			.call {
				arg := nodes[n.a].dim
				match functions[n.index] {
					'sqrt' {
						arg.pow(Ratio{1, 2}) or { return fail_at(tree, n, err.msg()) }
					}
					'abs' {
						arg
					}
					else {
						if !arg.is_dimensionless() {
							return fail_at(tree, n, 'the argument of ${n.text} must be a pure number, but it is ${describe(arg)}')
						}
						Dim{}
					}
				}
			}
		}
	}
	return Expr{
		text:     tree.text
		nodes:    nodes
		root:     tree.root
		dim:      nodes[tree.root].dim
		families: ctx.families.map(it.name)
	}
}

// power_dim is the dimension of a^b: a pure number's power is a pure
// number, and a quantity with units can be raised only to a fixed power,
// which must be a fraction p/q with q at most 10⁶.
fn power_dim(defs &Definitions, nodes []Node, tree Expr, n Node) !Dim {
	base, exponent := nodes[n.a].dim, nodes[n.b].dim
	if !exponent.is_dimensionless() {
		return fail_at(tree, leftmost(nodes, n.b), 'a power must be a pure number, but this one is ${describe(exponent)}')
	}
	if base.is_dimensionless() {
		return Dim{}
	}
	x := fold(defs, nodes, n.b) or {
		return fail_at(tree, leftmost(nodes, n.b), 'a quantity with units, here ${describe(base)}, can only be raised to a fixed power, and this power varies')
	}
	r := to_ratio(x) or {
		return fail_at(tree, leftmost(nodes, n.b), "a quantity with units, here ${describe(base)}, can only be raised to a fraction, and ${x} isn't one with a denominator up to 10⁶")
	}
	return base.pow(r) or { fail_at(tree, n, err.msg()) }
}

// fold is a subtree's value if it uses no variable, directly or through a
// definition.
fn fold(defs &Definitions, nodes []Node, i int) ?f64 {
	n := nodes[i]
	return match n.kind {
		.num {
			n.value
		}
		.var, .member {
			none
		}
		.constant {
			constants[n.index].value
		}
		.def {
			d := defs.exprs[n.index]
			fold(defs, d.nodes, d.root)?
		}
		.neg {
			-fold(defs, nodes, n.a)?
		}
		.add {
			fold(defs, nodes, n.a)? + fold(defs, nodes, n.b)?
		}
		.sub {
			fold(defs, nodes, n.a)? - fold(defs, nodes, n.b)?
		}
		.mul {
			fold(defs, nodes, n.a)? * fold(defs, nodes, n.b)?
		}
		.div {
			fold(defs, nodes, n.a)? / fold(defs, nodes, n.b)?
		}
		.pow {
			math.pow(fold(defs, nodes, n.a)?, fold(defs, nodes, n.b)?)
		}
		.call {
			apply(n.index, fold(defs, nodes, n.a)?)
		}
	}
}

// to_ratio is the simplest fraction within 1e-12 of x, its denominator at
// most 10⁶, by continued fractions: 2.5 is 5/2, 1/3 is 1/3, 1.37 is 137/100.
fn to_ratio(x f64) ?Ratio {
	if !math.is_finite(x) || math.abs(x) > 1e6 {
		return none
	}
	tol := 1e-12 * math.max(1.0, math.abs(x))
	mut h0, mut h1 := i64(0), i64(1)
	mut k0, mut k1 := i64(1), i64(0)
	mut y := x
	for _ in 0 .. 64 {
		a := i64(math.floor(y))
		h0, h1 = h1, a * h1 + h0
		k0, k1 = k1, a * k1 + k0
		// quantity.ratio holds the same bound, so this changes no answer; it
		// stops the next step's a·k1 overflowing an i64.
		if k1 > 1_000_000 {
			return none
		}
		if math.abs(f64(h1) / f64(k1) - x) <= tol {
			return quantity.ratio(h1, k1) or { return none }
		}
		f := y - f64(a)
		if f == 0.0 {
			break
		}
		y = 1.0 / f
	}
	return none
}

// describe is a dimension for a message: its name and its base units, or
// just its base units, or "a pure number".
fn describe(d Dim) string {
	if d.is_dimensionless() {
		return 'a pure number'
	}
	name := d.name()
	if name == d.str() {
		return 'in ${name}'
	}
	article := if name[0] in [`a`, `e`, `i`, `o`, `u`] { 'an' } else { 'a' }
	return '${article} ${name}, ${d.str()}'
}

// leftmost is where a subtree begins in the text: its leftmost operand's
// place, since a binary node's own place is its operator's. A message
// about a whole operand points there.
fn leftmost(nodes []Node, i int) Node {
	n := nodes[i]
	return match n.kind {
		.add, .sub, .mul, .div, .pow { leftmost(nodes, n.a) }
		else { n }
	}
}

fn fail_at(tree Expr, n Node, reason string) IError {
	return ExprError{
		text:   tree.text
		line:   n.line
		column: n.column
		reason: reason
	}
}
