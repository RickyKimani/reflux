// expr_render writes reflux's renderings of expressions, for
// tools/gen_m9_fixtures.py to hold against SymPy:
//
//     v run tools/expr_render.v cases.json [readings]
//
// cases.json is a JSON array of expression texts. Each is read in the
// fixture's context: the pure numbers a, b, x and y, and a family p of pure
// numbers with members CH4 and H2O; or with readings, T in K and P in Pa.
// The output, on stdout, is a JSON array with each case's plain, Unicode and
// LaTeX renderings, or its error.
module main

import expr
import json2
import os
import quantity

fn context() expr.Context {
	pure := quantity.Dim{}
	return expr.Context{
		variables: ['a', 'b', 'x', 'y'].map(expr.Variable{it, pure})
		families:  [
			expr.Family{'p', pure, ['CH4', 'H2O'], {
				'CH4': 0
				'H2O': 1
			}},
		]
	}
}

// The readings' context: T and P, for readings on a unit's scale.
fn readings_context() expr.Context {
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', quantity.Dimension.pressure.dim()},
		]
	}
}

struct Rendered {
	plain   string
	unicode string
	latex   string
	error   string
}

fn main() {
	if os.args.len !in [2, 3] || (os.args.len == 3 && os.args[2] != 'readings') {
		eprintln('usage: v run tools/expr_render.v cases.json [readings]')
		exit(2)
	}
	texts := json2.decode[[]string](os.read_file(os.args[1])!)!
	ctx := if os.args.len == 3 { readings_context() } else { context() }
	mut out := []Rendered{}
	for text in texts {
		e := expr.parse(ctx, &expr.Definitions{}, text) or {
			out << Rendered{
				error: err.msg()
			}
			continue
		}
		out << Rendered{
			plain:   e.plain()
			unicode: e.unicode()
			latex:   e.latex()
		}
	}
	println(json2.encode(out))
}
