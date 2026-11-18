module flowsheet

import expr
import kinetics
import math
import os
import quantity
import thermo
import toml
import toml.ast

// A file's fits (M8's D14, D15): [fit.NAME], in a format-3 file with the
// reactions they fit, the model whose kᵢⱼ they fit, or both. The fit module
// runs them. docs/flowsheet-format.md describes the keys.

// FitFile is a file's fits with what they fit: its components, its
// reactions if it has any, and its fluid if it has a [model].
pub struct FitFile {
pub:
	names    []string
	kinetics ?kinetics.Kinetics
	fluid    ?thermo.Fluid
	fits     []FitPlan
}

// FitKind is what a fit fits: rate constants to rates, or kᵢⱼ to bubble
// points.
pub enum FitKind {
	rates
	vle
}

// Part is which number of a named constant a fit varies: an Arrhenius
// constant's A or Eₐ, a van 't Hoff constant's A or ΔH, a constant's whole
// value, or a pair's kᵢⱼ.
pub enum Part {
	a
	ea
	dh
	value
	kij
}

// FitParameter is a number a fit varies, as the file names it, k1.A,
// K_E.dH, K_H or kij[ethane, propane], with its start in SI, its dimension,
// and the unit the file wrote it in, for the fitted section ('' for a pure
// number or SI).
pub struct FitParameter {
pub:
	name     string
	constant string // the named constant, k1; '' for a kᵢⱼ
	part     Part
	start    f64
	dim      quantity.Dim
	unit     string
	i        int // a kᵢⱼ's components
	j        int
}

// Weight is a column's measurement error: relative, a fraction of each
// measurement, or absolute, in SI.
pub struct Weight {
pub:
	column   string
	relative bool
	value    f64
	dim      quantity.Dim // an absolute weight's, which must be its column's
}

// FitPlan is a [fit.NAME] as read (D14, D15).
pub struct FitPlan {
pub:
	name       string
	kind       FitKind
	data       string // the data's file, as given
	path       string // and as found, beside the file
	parameters []FitParameter
	weights    []Weight
	line       int
}

// The sections a fits file may have that its reader leaves to others.
const fit_others = ['feed', 'unit', 'solver', 'estimate', 'spec', 'study', 'optimize']

// load_fits reads a file's fits. See load_fits_text.
pub fn load_fits(path string, params LoadParams) !FitFile {
	text := os.read_file(path) or { return file_error(0, path, "it can't be read: ${err.msg()}", []) }
	return load_fits_text(text, LoadParams{
		...params
		dir: os.dir(path)
	})
}

// load_fits_text reads a format-3 file's components, its [model] and its
// reactions if it has them, and its [fit.NAME]: each fit's parameters
// found among the reactions' named constants or the model's pairs, every
// key checked and every refusal at its line. A flowsheet's own sections are
// left alone, so a flowsheet file's fits read too.
pub fn load_fits_text(text string, params LoadParams) !FitFile {
	check_text(text)!
	doc := toml.parse_text(text) or {
		line, reason := syntax_error(err.msg())
		return file_error(line, '', reason, [])
	}
	root := doc.ast.table as map[string]ast.Value
	mut top := Section{
		path:  ''
		table: root
		text:  text
	}
	given := top.integer('format')!
	if given != 3 {
		return top.fail('format', 'fits need format 3, not ${given}', [])
	}
	mut names := top.child('components')!
	component_names := names.strings('names')!
	if component_names.len == 0 {
		return names.fail('names', 'a file needs at least one component', [])
	}
	names_line := names.line('names')
	data := if names.has('data') { names.strings('data')! } else { []string{} }
	names.finish()!
	db := database(params, data) or { return top.fail('components', err.msg(), []) }
	mut species := []kinetics.Species{}
	for i, name in component_names {
		if name in component_names[..i] {
			return names.fail('names', '${name} is listed twice', [])
		}
		r := db.find(name) or { return names.fail('names', err.msg(), []) }
		species << species_of(r)
	}
	mut fluid := ?thermo.Fluid(none)
	if top.has('model') {
		mut m := top.child('model')!
		fluid = read_model(mut m, db, component_names, names_line)!
	}
	mut k := ?kinetics.Kinetics(none)
	if top.has('reaction') {
		k = read_reactions(mut top, species)!
	}
	for s in fit_others {
		if top.has(s) {
			top.used << s
		}
	}
	mut fits := []FitPlan{}
	if top.has('fit') {
		mut all := top.child('fit')!
		for name in all.table.keys() {
			mut t := all.child(name)!
			fits << read_fit(mut t, name, component_names, k, fluid, params.dir)!
		}
	}
	top.finish()!
	return FitFile{
		names:    component_names
		kinetics: k
		fluid:    fluid
		fits:     fits
	}
}

// read_fit reads one [fit.NAME]: data, a CSV file beside this one;
// parameters, the numbers it varies; start, their first values if not the
// file's; and weights, each measured column's error.
fn read_fit(mut t Section, name string, names []string, k ?kinetics.Kinetics, fluid ?thermo.Fluid, dir string) !FitPlan {
	data := t.text('data')!
	path := if os.is_abs_path(data) { data } else { os.join_path(dir, data) }
	wanted := t.strings('parameters')!
	if wanted.len == 0 {
		return t.fail('parameters', 'a fit varies at least one parameter, as ["k1.A", "k1.Ea"]',
			[])
	}
	vle := wanted.any(it.starts_with('kij'))
	if vle && wanted.any(!it.starts_with('kij')) {
		return t.fail('parameters', 'a fit is of rates or of kij, not both: fit them in two [fit.NAME]',
			[])
	}
	mut starts := map[string]string{}
	if t.has('start') {
		mut s := t.child('start')!
		for key in s.table.keys() {
			if key !in wanted {
				return s.fail(key, '${key} is not one of the parameters', quantity.closest(key,
					wanted))
			}
			starts[key] = s.text(key)!
		}
		s.finish()!
	}
	mut parameters := []FitParameter{}
	for w in wanted {
		if parameters.any(it.name == w) {
			return t.fail('parameters', '${w} is listed twice', [])
		}
		p := if vle {
			kij_parameter(w, names, fluid) or { return t.fail('parameters', err.msg(), []) }
		} else {
			rate_parameter(w, k) or {
				hint := if err is FlowsheetError { err.suggestions } else { []string{} }
				return t.fail('parameters', err.msg(), hint)
			}
		}
		parameters << if w in starts {
			value := quantity_value(starts[w], p.dim) or {
				mut s := t.child('start')!
				return s.fail(w, err.msg(), [])
			}
			FitParameter{
				...p
				start: value
			}
		} else {
			p
		}
	}
	mut weights := []Weight{}
	if t.has('weights') {
		mut s := t.child('weights')!
		for column in s.table.keys() {
			text := s.text(column)!
			weights << weight_of(column, text) or { return s.fail(column, err.msg(), []) }
		}
		s.finish()!
	}
	t.finish()!
	return FitPlan{
		name:       name
		kind:       if vle { FitKind.vle } else { FitKind.rates }
		data:       data
		path:       path
		parameters: parameters
		weights:    weights
		line:       t.line_of_root(t.path)
	}
}

// rate_parameter is one of the reactions' named constants' numbers:
// NAME.A and NAME.Ea of an Arrhenius constant, NAME.A and NAME.dH of a van
// 't Hoff one, or NAME of a constant that doesn't depend on T.
fn rate_parameter(w string, kin ?kinetics.Kinetics) !FitParameter {
	k := kin or { return error('${w}: the file has no reactions whose constants it could be') }
	constants := k.constants()
	name, part := if w.contains('.') {
		w.all_before_last('.'), w.all_after_last('.')
	} else {
		w, ''
	}
	mut found := ?kinetics.NamedConstant(none)
	for c in constants {
		if c.name == name {
			found = c
		}
	}
	c := found or {
		mut every := []string{}
		for c in constants {
			every << c.name
		}
		return FlowsheetError{
			kind:        .bad_file
			reason:      'no constant ${name}: the reactions have ${every.join(', ')}'
			suggestions: quantity.closest(name, every)
		}
	}
	match c.value {
		kinetics.Arrhenius {
			v := c.value as kinetics.Arrhenius
			if part == 'A' {
				return parameter(w, name, .a, v.a)
			}
			if part == 'Ea' {
				return parameter(w, name, .ea, v.ea)
			}
			return error('${name} is ${c.what}, A·exp(−Ea/RT): fit ${name}.A or ${name}.Ea')
		}
		kinetics.VantHoff {
			v := c.value as kinetics.VantHoff
			if part == 'A' {
				return parameter(w, name, .a, v.a)
			}
			if part == 'dH' {
				return parameter(w, name, .dh, v.dh)
			}
			return error('${name} is ${c.what}, A·exp(−dH/RT): fit ${name}.A or ${name}.dH')
		}
		kinetics.Fixed, kinetics.Written {
			// A number with its unit is a constant that doesn't depend on T,
			// which a file writes in quotes; anything with T or a name in it
			// isn't.
			text := c.value.text()
			if _ := expr.parse(expr.Context{}, &expr.Definitions{}, text) {
				if part != '' {
					return error('${name} is ${c.what}, a number that does not depend on T: fit ${name}')
				}
				return parameter(w, name, .value, text)
			}
			return error("${name} is ${c.what}, written as an expression, which a fit can't take apart: give it as { A, Ea }, { A, dH } or a number")
		}
	}
}

// parameter is a constant's number, from its text as the file writes it,
// with its value in SI and its unit.
fn parameter(w string, name string, part Part, text string) !FitParameter {
	e := expr.parse(expr.Context{}, &expr.Definitions{}, text) or {
		return error('${w}: ${text} is not a number with its unit')
	}
	fields := text.trim_space().fields()
	return FitParameter{
		name:     w
		constant: name
		part:     part
		start:    e.eval(&expr.Definitions{}, [])
		dim:      e.dim
		unit:     if fields.len > 1 { fields[1..].join(' ') } else { '' }
	}
}

// kij_parameter is a pair's kᵢⱼ, kij[a, b], from the file's model.
fn kij_parameter(w string, names []string, fl ?thermo.Fluid) !FitParameter {
	f := fl or { return error('${w}: kij needs a [model], whose pairs it fits') }
	if !w.starts_with('kij[') || !w.ends_with(']') || w.count(',') != 1 {
		return error('${w}: a pair is written kij[a, b], two of the components')
	}
	inner := w.all_after('kij[').all_before_last(']')
	a, b := inner.all_before(',').trim_space(), inner.all_after(',').trim_space()
	i, j := names.index(a), names.index(b)
	if i < 0 || j < 0 {
		bad := if i < 0 { a } else { b }
		return error('${w}: ${bad} is not one of the components')
	}
	if i == j {
		return error('${w}: a pair needs two components')
	}
	if f.model !is thermo.Cubic {
		return error("${w}: kij is a cubic's, and the model is not one")
	}
	c := f.model as thermo.Cubic
	return FitParameter{
		name:  w
		part:  .kij
		start: c.kij[i * names.len + j]
		dim:   quantity.Dim{}
		i:     i
		j:     j
	}
}

// quantity_value reads a number with its unit, in SI, of a dimension.
fn quantity_value(text string, d quantity.Dim) !f64 {
	e := expr.parse(expr.Context{}, &expr.Definitions{}, text) or {
		return error('"${text}" is not a number with its unit')
	}
	if !e.dim.eq(d) {
		return error('"${text}" is ${e.dim.name()}, not ${d.name()}')
	}
	return e.eval(&expr.Definitions{}, [])
}

// weight_of reads a column's error: "5 %" of each measurement, or a number
// with the column's unit, which the fit checks against the column.
fn weight_of(column string, text string) !Weight {
	t := text.trim_space()
	if t.ends_with('%') {
		v := t.all_before_last('%').trim_space().f64()
		if !(v > 0.0) {
			return error('"${text}" must be above 0 %')
		}
		return Weight{
			column:   column
			relative: true
			value:    v / 100.0
		}
	}
	e := expr.parse(expr.Context{}, &expr.Definitions{}, t) or {
		return error('"${text}" is neither a percentage nor a number with its unit')
	}
	v := e.eval(&expr.Definitions{}, [])
	if !(v > 0.0) || !math.is_finite(v) {
		return error('"${text}" must be above 0')
	}
	return Weight{
		column: column
		value:  v
		dim:    e.dim
	}
}
