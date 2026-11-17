module flowsheet

import expr
import math
import quantity
import toml.ast

// The most points a study may have, all its inputs' values together.
const most_points = 1_000_000

// StudyPlan is a [study.NAME] as read (M8's D16): the inputs it varies and
// their values, the outputs it reports, and its chunk (D7). The study
// module runs it.
pub struct StudyPlan {
pub:
	name    string
	axes    []Axis // one or two
	outputs []Output
	chunk   int = 8
	line    int
}

// Axis is an input a study varies, and its values in SI, in order.
pub struct Axis {
pub:
	input  Input
	values []f64
	unit   string // the values' SI unit, '' for a pure number
}

// Output is an expression a study reports at every point, by name.
pub struct Output {
pub:
	name string
	expr expr.Expr
	unit string // its value's SI unit, '' for a pure number
}

// read_studies reads [study.NAME] (D16): vary, one or two inputs, each a
// key with from, to and points, or values; report, named expressions; and
// chunk, how many points run from one cold start, 8 by default (D7).
fn read_studies(mut studies Section, l Loaded) ![]StudyPlan {
	mut out := []StudyPlan{}
	for name in studies.table.keys() {
		mut t := studies.child(name)!
		mut axes := []Axis{}
		mut varying := t.tables('vary')!
		if varying.len == 0 || varying.len > 2 {
			return t.fail('vary', 'a study varies one input or two, each { key = "...", from = "...", to = "...", points = N } or { key = "...", values = [...] }',
				[])
		}
		for mut a in varying {
			key := a.text('key')!
			input := l.input(key) or { return at_key(err, a.line('key'), a.full('key')) }
			for spec in l.specs {
				if spec.input.key == key {
					return a.fail('key', "spec ${spec.name} varies ${key}, to meet its target, so a study can't",
						[])
				}
			}
			if axes.any(it.input.key == key) {
				return a.fail('key', '${key} is varied twice', [])
			}
			values := axis_values(mut a, input)!
			a.finish()!
			axes << Axis{
				input:  input
				values: values
				unit:   unit_of(input.dim)
			}
		}
		mut total := 1
		for a in axes {
			total *= a.values.len
		}
		if total > most_points {
			return t.fail('vary', '${total} points are more than a study takes, ${most_points}',
				[])
		}
		mut outputs := []Output{}
		mut report := t.child('report')!
		for out_name in report.table.keys() {
			text := report.text(out_name)!
			e := l.expression(text) or {
				if err is expr.ExprError {
					x := err as expr.ExprError
					return file_error(report.line(out_name), report.full(out_name), 'column ${x.column}: ${x.reason}',
						x.suggestions)
				}
				return report.fail(out_name, err.msg(), [])
			}
			outputs << Output{
				name: out_name
				expr: e
				unit: unit_of(e.dim)
			}
		}
		if outputs.len == 0 {
			return t.fail('report', 'a study reports at least one expression, as T_out = "product.T"',
				[])
		}
		report.finish()!
		chunk := if t.has('chunk') { t.integer('chunk')! } else { 8 }
		if chunk < 1 || chunk > most_points {
			return t.fail('chunk', 'a chunk of ${chunk} points', [])
		}
		t.finish()!
		out << StudyPlan{
			name:    name
			axes:    axes
			outputs: outputs
			chunk:   int(chunk)
			line:    studies.line_of_root('study.${name}')
		}
	}
	return out
}

// axis_values are an input's values: from, to and points, evenly spaced and
// ending at to exactly, or a list.
fn axis_values(mut a Section, input Input) ![]f64 {
	if a.has('values') {
		v := a.require('values')!
		if v !is []ast.Value || (v as []ast.Value).len == 0 {
			return a.fail('values', 'it needs a list of values, as ["700 degC", "800 degC"]',
				[])
		}
		mut out := []f64{}
		for x in v as []ast.Value {
			out << input.value_in(a, 'values', x)!
		}
		return out
	}
	if !a.has('from') || !a.has('to') || !a.has('points') {
		return a.fail('', 'give from, to and points, or values', [])
	}
	from := input.value_in(a, 'from', a.require('from')!)!
	to := input.value_in(a, 'to', a.require('to')!)!
	n := a.integer('points')!
	if n < 2 || n > most_points {
		return a.fail('points', 'from and to need 2 points or more, not ${n}', [])
	}
	mut out := []f64{len: int(n)}
	for i in 0 .. n {
		out[i] = from + (to - from) * f64(i) / f64(n - 1)
	}
	out[n - 1] = to
	return out
}

// value_in reads one of an input's values from a file, as its maker reads
// the input: a pure number, or a quantity in quotes with its unit.
fn (i Input) value_in(t Section, key string, v ast.Value) !f64 {
	if v is ast.Number {
		if i.dim.is_dimensionless() {
			return v.f64()
		}
		return t.fail(key, 'a value of ${i.key} needs its unit, in quotes, as "${v.text} ${si_text(i.dim)}"',
			[])
	}
	if v is ast.Quoted {
		x := i.parse(v.text) or {
			return file_error(v.pos.line_nr, t.full(key), err.msg(), []string{})
		}
		if !math.is_finite(x) {
			return file_error(v.pos.line_nr, t.full(key), '${v.text} is not a finite value',
				[]string{})
		}
		return x
	}
	return t.fail(key, 'it needs a number, or a quantity in quotes', [])
}

// unit_of is a dimension's SI unit for a report, '' for a pure number.
fn unit_of(d quantity.Dim) string {
	return if d.is_dimensionless() { '' } else { si_text(d) }
}
