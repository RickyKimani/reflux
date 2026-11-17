module study

import encoding.csv
import flowsheet
import json2
import math
import strings

// The study report's own format (M8's D17), which docs/report-format.md
// describes beside the solve's: it changes only with a new number.
pub const report_format = 1

// Report is a study's results, or why it couldn't run: what `reflux study`
// writes as JSON and CSV, and prints as a table. Every number is SI.
pub struct Report {
pub:
	format  int
	reflux  string
	file    string
	study   string
	status  string // 'ran', even with failed points, or 'failed'
	error   ?flowsheet.ErrorReport @[json_null]
	inputs  []InputReport
	outputs []OutputReport
	specs   []string // the file's specs, whose values each point gives
	chunk   int
	points  []PointReport
}

// InputReport is an input a study varies.
pub struct InputReport {
pub:
	key  string
	unit string // its values', in SI: '' for a pure number
}

// OutputReport is an expression a study reports.
pub struct OutputReport {
pub:
	name       string
	expression string // as reflux writes it back
	unit       string // its values', in SI: '' for a pure number
}

// PointReport is a point of a study.
pub struct PointReport {
pub:
	index            int
	inputs           []f64
	status           string // 'solved', 'no_value' or 'failed'
	outputs          []Value
	error            string // why a failed point failed, else ''
	specs            []f64
	passes           int
	spec_evaluations int
	steps            int
	warm             bool
}

// Value is an output at a point: its value, or null and why it has none.
pub struct Value {
pub:
	value  ?f64 @[json_null]
	reason string
}

// report is a run study's report.
pub fn (r Result) report(file string, specs []string) Report {
	plan := r.plan
	return Report{
		format:  report_format
		reflux:  flowsheet.version()
		file:    file
		study:   plan.name
		status:  'ran'
		inputs:  plan.axes.map(InputReport{
			key:  it.input.key
			unit: it.unit
		})
		outputs: plan.outputs.map(OutputReport{
			name:       it.name
			expression: it.expr.plain()
			unit:       it.unit
		})
		specs:   specs
		chunk:   plan.chunk
		points:  r.points.map(point_report(it))
	}
}

fn point_report(p Point) PointReport {
	mut cells := []Value{}
	for k, x in p.outputs {
		cells << if math.is_finite(x) {
			Value{
				value: x
			}
		} else {
			Value{
				value:  none
				reason: p.reasons[k]
			}
		}
	}
	return PointReport{
		index:            p.index
		inputs:           p.inputs
		status:           p.status.str()
		outputs:          cells
		error:            p.error
		specs:            p.specs
		passes:           p.passes
		spec_evaluations: p.evaluations
		steps:            p.steps
		warm:             p.warm
	}
}

// load_failure is the report of a study whose file didn't load.
pub fn load_failure(file string, err IError) Report {
	return Report{
		format: report_format
		reflux: flowsheet.version()
		file:   file
		status: 'failed'
		error:  flowsheet.error_report(err)
	}
}

// json is the report as JSON, a value to a line.
pub fn (r Report) json() string {
	return json2.encode(r, prettify: true, indent_string: '  ') + '\n'
}

// with_unit is a column's name with its unit, as the CSV and the table
// write it.
fn with_unit(name string, unit string) string {
	return if unit == '' { name } else { '${name} [${unit}]' }
}

// csv is a row for each point, in the grid's order: its index, its inputs,
// its outputs, NaN where one has no value, and its status. Each number
// reads back as the same double.
pub fn (r Report) csv() string {
	mut head := ['point']
	head << r.inputs.map(with_unit(it.key, it.unit))
	head << r.outputs.map(with_unit(it.name, it.unit))
	head << 'status'
	mut w := csv.new_writer()
	w.write(head) or {}
	for p in r.points {
		mut row := [p.index.str()]
		row << p.inputs.map(it.str())
		for v in p.outputs {
			row << if x := v.value { x.str() } else { 'NaN' }
		}
		row << p.status
		w.write(row) or {}
	}
	return w.str()
}

// table is the study for a person to read: what ran, a row for each point
// to 5 significant figures, and why any output had no value or any point
// failed.
pub fn (r Report) table() string {
	mut b := strings.new_builder(4096)
	n_solved := r.points.filter(it.status == 'solved').len
	none_ := r.points.filter(it.status == 'no_value').len
	failed_ := r.points.filter(it.status == 'failed').len
	warm := r.points.filter(it.warm).len
	b.writeln('${r.file}, study ${r.study}: ${r.points.len} points, ${n_solved} solved, ${none_} with an output that has no value, ${failed_} failed; in chunks of ${r.chunk}, ${warm} from a neighbour and ${r.points.len - warm} cold.')
	mut head := r.inputs.map(with_unit(it.key, it.unit))
	head << r.outputs.map(with_unit(it.name, it.unit))
	mut rows := [][]string{}
	for p in r.points {
		mut row := p.inputs.map(flowsheet.figures(it))
		for k, v in p.outputs {
			row << if x := v.value {
				flowsheet.figures(x)
			} else if p.status == 'failed' && k == 0 {
				'failed'
			} else {
				'—'
			}
		}
		rows << row
	}
	mut widths := head.map(math.max(10, it.runes().len))
	for row in rows {
		for i, cell in row {
			widths[i] = math.max(widths[i], cell.runes().len)
		}
	}
	b.writeln('')
	mut titles := []string{}
	for i, h in head {
		titles << pad(h, widths[i])
	}
	b.writeln(titles.join('  '))
	for row in rows {
		mut cells := []string{}
		for i, cell in row {
			cells << pad(cell, widths[i])
		}
		b.writeln(cells.join('  '))
	}
	// Each output's reasons, once each, with how many points had it.
	mut said := map[string]int{}
	mut seen := []string{}
	for p in r.points {
		for k, v in p.outputs {
			if v.reason != '' {
				key := '${r.outputs[k].name}: ${v.reason}'
				if key !in said {
					seen << key
				}
				said[key]++
			}
		}
	}
	if seen.len > 0 {
		b.writeln('')
		b.writeln('Outputs with no value:')
		for key in seen {
			where := if said[key] == r.points.len {
				'at every point'
			} else {
				'at ${said[key]} points'
			}
			b.writeln('  ${key.all_before(': ')}, ${where}: ${key.all_after(': ')}')
		}
	}
	failures := r.points.filter(it.status == 'failed')
	if failures.len > 0 {
		b.writeln('')
		b.writeln('Failed points:')
		for p in failures {
			mut at := []string{}
			for i, x in p.inputs {
				at << '${r.inputs[i].key} = ${with_unit(flowsheet.figures(x), r.inputs[i].unit).replace('[', '').replace(']', '')}'
			}
			b.writeln('  point ${p.index} (${at.join(', ')}): ${p.error}')
		}
	}
	return b.str()
}

// pad right-aligns a cell.
fn pad(s string, width int) string {
	n := s.runes().len
	return if n >= width { s } else { ' '.repeat(width - n) + s }
}
