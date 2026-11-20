module expr

import math
import numeric

// jet is the expression's value with its first two derivatives along one
// direction, by forward-mode differentiation over the tree: x is each
// slot's value as a jet, seeded with the direction (numeric.along for the
// variable differentiated by, numeric.constant for the rest). Its value is
// eval's, to the bit, since every node computes it as eval does.
//
// The tape gives first derivatives, by a reverse sweep for every output at
// once. A second derivative in one variable is what an activity model's
// τ(T) needs for its excess heat capacity (M6 step 2), and this gives it
// without a second tape. Like eval, it walks the tree: it's for values
// kept until T moves, not for inner loops.
pub fn (e &Expr) jet(defs &Definitions, x []numeric.Jet) numeric.Jet {
	return e.jet_of(defs, x, e.root)
}

fn (e &Expr) jet_of(defs &Definitions, x []numeric.Jet, i int) numeric.Jet {
	n := e.nodes[i]
	return match n.kind {
		.num { numeric.constant(n.value) }
		.var, .member { x[n.index] }
		.constant { numeric.constant(constants[n.index].value) }
		.def { defs.exprs[n.index].jet(defs, x) }
		.neg { e.jet_of(defs, x, n.a).neg() }
		.add { e.jet_of(defs, x, n.a) + e.jet_of(defs, x, n.b) }
		.sub { e.jet_of(defs, x, n.a) - e.jet_of(defs, x, n.b) }
		.mul { e.jet_of(defs, x, n.a) * e.jet_of(defs, x, n.b) }
		.div { e.jet_of(defs, x, n.a) / e.jet_of(defs, x, n.b) }
		.pow { raise_jet(e.jet_of(defs, x, n.a), e.jet_of(defs, x, n.b)) }
		.call { apply_jet(n.index, e.jet_of(defs, x, n.a)) }
		.reading {
			(e.jet_of(defs, x, n.a) - numeric.constant(n.post)) / numeric.constant(n.value) - numeric.constant(n.pre)
		}
		.scaled { e.jet_of(defs, x, n.a) * numeric.constant(n.value) }
	}
}

// raise_jet is aᵇ as raise computes it: a whole power up to 16 that doesn't
// move by repeated multiplication, any other by math.pow, so that its value
// is raise's bits.
fn raise_jet(a numeric.Jet, b numeric.Jet) numeric.Jet {
	if b.d != 0.0 || b.dd != 0.0 {
		return a.pow_jet(b)
	}
	p := b.v
	if p == math.trunc(p) && math.abs(p) <= 16.0 {
		n := int(math.abs(p))
		mut y := numeric.constant(1.0)
		for _ in 0 .. n {
			y = y * a
		}
		return if p < 0.0 { numeric.constant(1.0) / y } else { y }
	}
	return a.pow(p)
}

// apply_jet is function f of x, as apply.
fn apply_jet(f int, x numeric.Jet) numeric.Jet {
	return match functions[f] {
		'exp' { x.exp() }
		'ln' { x.log() }
		'log10' { x.log10() }
		'sqrt' { x.sqrt() }
		else { x.abs() }
	}
}
