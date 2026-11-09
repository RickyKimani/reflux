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
		.pow { math.pow(e.value_of(defs, x, n.a), e.value_of(defs, x, n.b)) }
		.call { apply(n.index, e.value_of(defs, x, n.a)) }
	}
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
