module study

import encoding.csv
import expr
import flowsheet
import json2
import math
import strings

// The optimisation report's own format (M8's D17), which
// docs/report-format.md describes beside the study's: it changes only with
// a new number.
pub const optimize_report_format = 1

// OptimizeReport is an optimisation's outcome, or why it couldn't run: what
// `reflux optimize` writes as JSON and CSV, and prints as a table. Every
// number is SI.
pub struct OptimizeReport {
pub:
	format      int
	reflux      string
	file        string
	optimize    string
	status      string // 'optimal', or 'failed' with its reason or error
	error       ?flowsheet.ErrorReport @[json_null] // the file's, if it didn't load
	reason      string // why there's no optimum, else ''
	goal        string // 'minimize' or 'maximize'
	objective   ObjectiveReport
	inputs      []VariedReport
	constraints []ConstraintReport
	specs       []SpecValue
	iterations  []IterationReport
	evaluations int // solves at points the optimiser chose
	solves      int // every solve, the gradients' included
}

// ObjectiveReport is what's minimised or maximised, and its value at the
// optimum, or null.
pub struct ObjectiveReport {
pub:
	expression string
	unit       string
	value      ?f64 @[json_null]
}

// VariedReport is an input an optimisation varies: its bounds, its start,
// its value at the optimum (the last iterate's, if it failed), which bound
// it's on, if either, and that bound's multiplier, how fast the objective
// improves as the bound is relaxed, per unit of the input.
pub struct VariedReport {
pub:
	key        string
	unit       string
	lower      f64
	upper      f64
	start      f64
	value      f64
	bound      string // 'lower', 'upper' or ''
	multiplier f64
}

// ConstraintReport is a constraint at the optimum: its two sides, whether
// it's active, and its multiplier, as ConstraintValue says; null for each
// number if there's no optimum.
pub struct ConstraintReport {
pub:
	equation   string
	unit       string
	lhs        ?f64 @[json_null]
	rhs        ?f64 @[json_null]
	active     bool
	multiplier ?f64 @[json_null]
}

// SpecValue is a spec's value at the optimum.
pub struct SpecValue {
pub:
	name  string
	value f64
}

// IterationReport is an iterate: its inputs and objective, in SI, its worst
// constraint's violation and the Lagrangian's gradient, both as the
// optimiser scales them, and the step taken from it.
pub struct IterationReport {
pub:
	inputs       []f64
	objective    f64
	violation    f64
	stationarity f64
	step         f64
}

// report is an optimisation's report.
pub fn (o Optimum) report(file string, specs []string) OptimizeReport {
	plan := o.plan
	mut inputs := []VariedReport{}
	for k, v in plan.vary {
		x := o.inputs[k]
		bound := if x == v.lo {
			'lower'
		} else if x == v.hi {
			'upper'
		} else {
			''
		}
		inputs << VariedReport{
			key:        v.input.key
			unit:       v.unit
			lower:      v.lo
			upper:      v.hi
			start:      v.start
			value:      x
			bound:      bound
			multiplier: if o.ok { (o.lower[k] or { 0.0 }) + (o.upper[k] or { 0.0 }) } else { 0.0 }
		}
	}
	mut constraints := []ConstraintReport{}
	for i, q in plan.constraints {
		constraints << if o.ok {
			c := o.constraints[i]
			ConstraintReport{
				equation:   q.plain()
				unit:       unit_of_dim(q.lhs)
				lhs:        c.lhs
				rhs:        c.rhs
				active:     c.active
				multiplier: c.multiplier
			}
		} else {
			ConstraintReport{
				equation:   q.plain()
				unit:       unit_of_dim(q.lhs)
				lhs:        none
				rhs:        none
				multiplier: none
			}
		}
	}
	mut spec_values := []SpecValue{}
	for k, name in specs {
		if v := o.specs[k] {
			spec_values << SpecValue{
				name:  name
				value: v
			}
		}
	}
	return OptimizeReport{
		format:      optimize_report_format
		reflux:      flowsheet.version()
		file:        file
		optimize:    plan.name
		status:      if o.ok { 'optimal' } else { 'failed' }
		reason:      o.error
		goal:        plan.objective.name
		objective:   ObjectiveReport{
			expression: plan.objective.expr.plain()
			unit:       plan.objective.unit
			value:      if o.ok { o.objective } else { none }
		}
		inputs:      inputs
		constraints: constraints
		specs:       spec_values
		iterations:  o.iterations.map(IterationReport{
			inputs:       it.inputs
			objective:    it.objective
			violation:    it.violation
			stationarity: it.stationarity
			step:         it.step
		})
		evaluations: o.evaluations
		solves:      o.solves
	}
}

// unit_of_dim is an expression's SI unit, '' for a pure number.
fn unit_of_dim(e expr.Expr) string {
	return if e.dim.is_dimensionless() { '' } else { flowsheet.si_text(e.dim) }
}

// optimize_load_failure is the report of an optimisation whose file didn't
// load.
pub fn optimize_load_failure(file string, err IError) OptimizeReport {
	return OptimizeReport{
		format: optimize_report_format
		reflux: flowsheet.version()
		file:   file
		status: 'failed'
		error:  flowsheet.error_report(err)
	}
}

// json is the report as JSON, a value to a line.
pub fn (r OptimizeReport) json() string {
	return json2.encode(r, prettify: true, indent_string: '  ') + '\n'
}

// csv is a row for each iterate: its number from 0, its inputs, its
// objective, its violation and stationarity, and the step taken from it.
// Each number reads back as the same double.
pub fn (r OptimizeReport) csv() string {
	mut head := ['iteration']
	head << r.inputs.map(with_unit(it.key, it.unit))
	head << with_unit(r.goal, r.objective.unit)
	head << ['violation', 'stationarity', 'step']
	mut w := csv.new_writer()
	w.write(head) or {}
	for k, it in r.iterations {
		mut row := [k.str()]
		row << it.inputs.map(it.str())
		row << [it.objective.str(), it.violation.str(), it.stationarity.str(), it.step.str()]
		w.write(row) or {}
	}
	return w.str()
}

// table is the optimisation for a person to read: the optimum or why
// there's none, each input with its bounds, each constraint, and each
// iterate, to 5 significant figures.
pub fn (r OptimizeReport) table() string {
	mut b := strings.new_builder(4096)
	if r.status == 'optimal' {
		b.writeln('${r.file}, optimisation ${r.optimize}: an optimum in ${r.iterations.len} iterations, ${r.evaluations} evaluations and ${r.solves} solves.')
	} else {
		b.writeln('${r.file}, optimisation ${r.optimize}: no optimum, after ${r.iterations.len} iterations, ${r.evaluations} evaluations and ${r.solves} solves: ${r.reason}')
	}
	value := if v := r.objective.value {
		' = ${with_unit_value(v, r.objective.unit)}'
	} else {
		''
	}
	b.writeln('${r.goal} ${r.objective.expression}${value}')
	b.writeln('')
	mut rows := [['input', 'least', 'most', 'start', 'value', 'bound', 'multiplier']]
	for v in r.inputs {
		rows << [with_unit(v.key, v.unit), flowsheet.figures(v.lower), flowsheet.figures(v.upper),
			flowsheet.figures(v.start), flowsheet.figures(v.value),
			if v.bound == '' { '—' } else { v.bound },
			if v.bound == '' { '—' } else { flowsheet.figures(v.multiplier) }]
	}
	write_rows(mut b, rows)
	if r.constraints.len > 0 {
		b.writeln('')
		mut c_rows := [['constraint', 'left', 'right', 'active', 'multiplier']]
		for c in r.constraints {
			c_rows << [c.equation, if x := c.lhs { flowsheet.figures(x) } else { '—' },
				if x := c.rhs { flowsheet.figures(x) } else { '—' }, if c.lhs == none {
					'—'
				} else if c.active {
					'yes'
				} else {
					'no'
				}, if x := c.multiplier { flowsheet.figures(x) } else { '—' }]
		}
		write_rows(mut b, c_rows)
	}
	if r.iterations.len > 0 {
		b.writeln('')
		mut i_rows := [['iteration']]
		i_rows[0] << r.inputs.map(it.key)
		i_rows[0] << [r.goal, 'violation', 'stationarity', 'step']
		for k, it in r.iterations {
			mut row := [k.str()]
			row << it.inputs.map(flowsheet.figures(it))
			row << [flowsheet.figures(it.objective), flowsheet.figures(it.violation),
				flowsheet.figures(it.stationarity), flowsheet.figures(it.step)]
			i_rows << row
		}
		write_rows(mut b, i_rows)
	}
	return b.str()
}

// with_unit_value is a value to 5 significant figures, with its unit.
fn with_unit_value(v f64, unit string) string {
	return if unit == '' { flowsheet.figures(v) } else { '${flowsheet.figures(v)} ${unit}' }
}

// write_rows writes rows as columns, the first row their titles.
fn write_rows(mut b strings.Builder, rows [][]string) {
	mut widths := []int{len: rows[0].len}
	for row in rows {
		for i, cell in row {
			widths[i] = math.max(widths[i], cell.runes().len)
		}
	}
	for row in rows {
		mut cells := []string{}
		for i, cell in row {
			cells << pad(cell, widths[i])
		}
		b.writeln(cells.join('  ').trim_right(' '))
	}
}
