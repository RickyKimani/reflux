module flowsheet

import expr
import toml.ast

// OptimizePlan is an [optimize.NAME] as read (M8's D13): the inputs it
// varies within their bounds, the expression it minimises or maximises, and
// the equations and inequalities it's subject to. The study module runs it.
pub struct OptimizePlan {
pub:
	name        string
	vary        []Varied
	objective   Output // its name is minimize or maximize
	maximize    bool
	constraints []expr.Equation // as the file lists them
	line        int
}

// Varied is an input an optimisation varies: its bounds and its start, in
// SI.
pub struct Varied {
pub:
	input Input
	lo    f64
	hi    f64
	start f64
	unit  string // the values' SI unit, '' for a pure number
}

// read_optimizations reads [optimize.NAME] (D13): vary, a list of inputs,
// each { key = "...", between = [least, most] } with an optional start, the
// file's value if not given; minimize or maximize, an expression; and
// subject_to, an optional list of equations and inequalities over the
// flowsheet's variables.
fn read_optimizations(mut tables Section, l Loaded) ![]OptimizePlan {
	mut out := []OptimizePlan{}
	for name in tables.table.keys() {
		mut t := tables.child(name)!
		mut varying := t.tables('vary')!
		if varying.len == 0 {
			return t.fail('vary', 'an optimisation varies at least one input, as [{ key = "unit.tube.volume", between = ["0.1 m^3", "10 m^3"] }]',
				[])
		}
		mut vary := []Varied{}
		for mut a in varying {
			key := a.text('key')!
			input := l.input(key) or { return at_key(err, a.line('key'), a.full('key')) }
			for spec in l.specs {
				if spec.input.key == key {
					return a.fail('key', "spec ${spec.name} varies ${key}, to meet its target, so an optimisation can't",
						[])
				}
			}
			if vary.any(it.input.key == key) {
				return a.fail('key', '${key} is varied twice', [])
			}
			v := a.require('between')!
			if v !is []ast.Value || (v as []ast.Value).len != 2 {
				return a.fail('between', 'it needs two values, the least and the most ${key} may take, as ["0.1 m^3", "10 m^3"]',
					[])
			}
			pair := v as []ast.Value
			lo := input.value_in(a, 'between', pair[0])!
			hi := input.value_in(a, 'between', pair[1])!
			if !(lo < hi) {
				return a.fail('between', 'its least value must be below its most', [])
			}
			start := if a.has('start') {
				input.value_in(a, 'start', a.require('start')!)!
			} else {
				input.value
			}
			if !(start >= lo && start <= hi) {
				where := if a.has('start') {
					'its start'
				} else {
					"the file's ${key}, ${input.text},"
				}
				return a.fail(if a.has('start') { 'start' } else { 'between' }, '${where} is outside it: an optimisation starts within its bounds',
					[])
			}
			a.finish()!
			vary << Varied{
				input: input
				lo:    lo
				hi:    hi
				start: start
				unit:  unit_of(input.dim)
			}
		}
		minimize, maximize := t.has('minimize'), t.has('maximize')
		if minimize && maximize {
			return t.fail('maximize', 'an optimisation has one objective, minimize or maximize, not both',
				[])
		}
		if !minimize && !maximize {
			return t.fail('', 'an optimisation needs an objective, minimize = "..." or maximize = "..."',
				[])
		}
		goal := if maximize { 'maximize' } else { 'minimize' }
		text := t.text(goal)!
		objective := l.expression(text) or {
			if err is expr.ExprError {
				x := err as expr.ExprError
				return file_error(t.line(goal), t.full(goal), 'column ${x.column}: ${x.reason}',
					x.suggestions)
			}
			return t.fail(goal, err.msg(), [])
		}
		mut constraints := []expr.Equation{}
		if t.has('subject_to') {
			for c in t.strings('subject_to')! {
				q := expr.equation(l.variables.context, &no_definitions, c) or {
					if err is expr.ExprError {
						e := err as expr.ExprError
						return file_error(t.line('subject_to'), t.full('subject_to'), '"${c}", column ${e.column}: ${e.reason}',
							e.suggestions)
					}
					return t.fail('subject_to', err.msg(), [])
				}
				constraints << q
			}
		}
		t.finish()!
		out << OptimizePlan{
			name:        name
			vary:        vary
			objective:   Output{
				name: goal
				expr: objective
				unit: unit_of(objective.dim)
			}
			maximize:    maximize
			constraints: constraints
			line:        tables.line_of_root('optimize.${name}')
		}
	}
	return out
}
