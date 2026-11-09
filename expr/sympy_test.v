// The expression language against SymPy and NIST (M9 step 2, criteria 1
// and 2), from testdata/m9_expressions.json: each case read as SymPy reads
// it, its renderings the ones SymPy's reading confirmed, and the built-in
// constants the SI's exact values.
import expr
import json2
import math
import os
import quantity

struct Case {
	text    string
	points  [][]f64
	values  []f64
	scales  []f64
	plain   string
	unicode string
	latex   string
}

struct ConstantRow {
	name  string
	exact string
	value f64
	unit  string
}

struct Fixture {
	cases     []Case
	constants []ConstantRow
}

fn load() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_expressions.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

// The fixture's context, as tools/expr_render.v has it.
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

fn test_every_case_reads_as_sympy_reads_it() {
	// SymPy's value at every point, to 1e-12 of the point's scale, the
	// largest magnitude of any subexpression: rounding moves a value by a
	// few ulps of that, and a wrong precedence or association by far more.
	f := load()!
	assert f.cases.len >= 250
	defs := expr.Definitions{}
	mut worst := 0.0
	for c in f.cases {
		e := expr.parse(context(), &defs, c.text) or { panic('${c.text}: ${err.msg()}') }
		for k, x in c.points {
			got := e.eval(&defs, x)
			err := math.abs(got - c.values[k]) / math.max(1.0, c.scales[k])
			worst = math.max(worst, err)
			assert err <= 1e-12, "${c.text} at ${x}: ${got} vs SymPy's ${c.values[k]}"
		}
	}
	println("worst error against SymPy, over each point's scale: ${worst:.2e}")
}

fn test_every_rendering_is_the_one_sympy_confirmed() {
	// tools/gen_m9_fixtures.py decoded each rendering and found it to mean
	// the text, by SymPy, at every point. Any change to how reflux writes
	// an expression shows here, and means running the generator again.
	f := load()!
	for c in f.cases {
		e := expr.parse(context(), &expr.Definitions{}, c.text)!
		assert e.plain() == c.plain, '${c.text}: plain ${e.plain()}, not ${c.plain}'
		assert e.unicode() == c.unicode, '${c.text}: unicode ${e.unicode()}, not ${c.unicode}'
		assert e.latex() == c.latex, '${c.text}: latex ${e.latex()}, not ${c.latex}'
		// And the plain form reads back to the same tree.
		back := expr.parse(context(), &expr.Definitions{}, c.plain)!
		assert back.nodes.len == e.nodes.len && back.root == e.root, c.text
		for i, n in e.nodes {
			m := back.nodes[i]
			assert n.kind == m.kind && n.a == m.a && n.b == m.b && n.index == m.index
				&& math.f64_bits(n.value) == math.f64_bits(m.value) && n.text == m.text, '${c.text}: node ${i}'
		}
	}
}

fn test_constants_are_the_si_values() {
	// Each constant is the double nearest its exact value, as the SI
	// Brochure defines it, R and F computed exactly; and its unit is its
	// dimension's, at scale 1.
	f := load()!
	assert f.constants.len == expr.constants.len
	for i, row in f.constants {
		c := expr.constants[i]
		assert c.name == row.name
		assert math.f64_bits(c.value) == math.f64_bits(row.value), '${c.name}: ${c.value} vs ${row.exact}'
		assert c.unit == row.unit, c.name
		if c.unit != '' {
			u := quantity.unit(c.unit)!
			assert u.dim.eq(c.dim) && u.scale == 1.0, c.name
		}
	}
}
