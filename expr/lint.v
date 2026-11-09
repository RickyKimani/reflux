module expr

import math

// LookAlike is a number typed where a built-in constant was likely meant:
// one with the constant's dimension, within 1% of it but not it, such as
// 8.314 J/(mol*K) for R. relative is (number − constant)/constant. A paper
// fitted with 8.314 may still write it (D2); this says how far it is.
pub struct LookAlike {
pub:
	node     int
	constant string
	value    f64
	relative f64
}

// look_alike_limit is how near a number must be to a constant to be one.
const look_alike_limit = 0.01

// look_alikes lists the expression's numbers that look like a built-in
// constant. The standard pressure is left out: a pressure near 1 bar, such
// as 1 atm, is usually meant to be what it is.
pub fn (e &Expr) look_alikes() []LookAlike {
	mut out := []LookAlike{}
	for i, n in e.nodes {
		if n.kind != .num {
			continue
		}
		for c in constants {
			if c.name == 'p_std' || !n.dim.eq(c.dim) {
				continue
			}
			rel := (n.value - c.value) / c.value
			if rel != 0.0 && math.abs(rel) <= look_alike_limit {
				out << LookAlike{
					node:     i
					constant: c.name
					value:    n.value
					relative: rel
				}
			}
		}
	}
	return out
}
