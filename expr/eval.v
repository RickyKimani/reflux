module expr

import math

// eval is the expression's value, in SI, with x its inputs in the context's
// slots (Context.slots), also in SI. It walks the tree, which is the
// reference that step 3's tape is held to; it's not for the inner loops.
pub fn (e &Expr) eval(defs &Definitions, x []f64) f64 {
	return e.value_of(defs, x, e.root)
}

fn (e &Expr) value_of(defs &Definitions, x []f64, i int) f64 {
	n := e.nodes[i]
	return match n.kind {
		.num { n.value }
		.var, .member { x[n.index] }
		.constant { constants[n.index].value }
		.def { defs.exprs[n.index].eval(defs, x) }
		.neg { -e.value_of(defs, x, n.a) }
		.add { e.value_of(defs, x, n.a) + e.value_of(defs, x, n.b) }
		.sub { e.value_of(defs, x, n.a) - e.value_of(defs, x, n.b) }
		.mul { e.value_of(defs, x, n.a) * e.value_of(defs, x, n.b) }
		.div { e.value_of(defs, x, n.a) / e.value_of(defs, x, n.b) }
		.pow { raise(e.value_of(defs, x, n.a), e.value_of(defs, x, n.b)) }
		.call { apply(n.index, e.value_of(defs, x, n.a)) }
		.reading { read_on(n, e.value_of(defs, x, n.a)) }
		.scaled { e.value_of(defs, x, n.a) * n.value }
	}
}

// raise is a^b: a whole power up to 16 by repeated multiplication, a
// negative one as 1 over that, whose bits IEEE arithmetic makes the same
// on every platform; any other by math.pow, the C library's. The tree and
// the tape both raise this way, so they give the same bits.
@[inline]
fn raise(a f64, b f64) f64 {
	if b == math.trunc(b) && math.abs(b) <= 16.0 {
		n := int(math.abs(b))
		mut y := 1.0
		for _ in 0 .. n {
			y *= a
		}
		return if b < 0.0 { 1.0 / y } else { y }
	}
	return math.pow(a, b)
}

// apply is function f of x.
fn apply(f int, x f64) f64 {
	return match functions[f] {
		'exp' { math.exp(x) }
		'ln' { math.log(x) }
		'log10' { math.log10(x) }
		'sqrt' { math.sqrt(x) }
		else { math.abs(x) }
	}
}
