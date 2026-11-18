module fit

import encoding.csv
import flowsheet
import json2
import kinetics
import math
import quantity
import strings

// The fit report's own format (M8's D17), which docs/report-format.md
// describes beside the study's and the optimisation's: it changes only
// with a new number.
pub const report_format = 1

// Report is a fit's outcome, or why it couldn't run: what `reflux fit`
// writes as JSON and CSV, and prints as a table. Every number is SI, each
// parameter in the file's own form.
pub struct Report {
pub:
	format       int
	reflux       string
	file         string
	fit          string
	status       string // 'fitted', or 'failed' with its reason or error
	error        ?flowsheet.ErrorReport @[json_null]
	reason       string
	kind         string // 'rates' or 'vle'
	data         string // the data's file, as the fit names it
	parameters   []ParameterReport
	correlation  [][]f64
	centred      []Centred
	measurements int // the residuals
	dof          int
	s            ?f64 @[json_null] // the residuals' standard deviation, √(Σr²/dof)
	t            ?f64 @[json_null] // Student's t at 0.975 on dof
	residuals    []Residual
	section      string // the fitted constants as the file writes them
	iterations   int
	evaluations  int
}

// ParameterReport is a parameter: its name, its SI unit and the file's,
// its start, and its value, standard error and 95% interval, in SI.
pub struct ParameterReport {
pub:
	name           string
	unit           string
	file_unit      string
	start          f64
	value          ?f64 @[json_null]
	standard_error ?f64 @[json_null]
	lower          ?f64 @[json_null]
	upper          ?f64 @[json_null]
}

// report is a fit's report, from the file it fitted.
pub fn (f Fitted) report(file string, fits flowsheet.FitFile) Report {
	plan := f.plan
	mut parameters := []ParameterReport{}
	for q, p in plan.parameters {
		parameters << if f.ok {
			ParameterReport{
				name:           p.name
				unit:           si_unit(p.dim)
				file_unit:      p.unit
				start:          p.start
				value:          f.values[q]
				standard_error: f.standard_errors[q]
				lower:          f.lower[q]
				upper:          f.upper[q]
			}
		} else {
			ParameterReport{
				name:           p.name
				unit:           si_unit(p.dim)
				file_unit:      p.unit
				start:          p.start
				value:          none
				standard_error: none
				lower:          none
				upper:          none
			}
		}
	}
	return Report{
		format:       report_format
		reflux:       flowsheet.version()
		file:         file
		fit:          plan.name
		status:       if f.ok { 'fitted' } else { 'failed' }
		reason:       f.error
		kind:         plan.kind.str()
		data:         plan.data
		parameters:   parameters
		correlation:  f.correlation
		centred:      f.centred
		measurements: f.residuals.len
		dof:          f.dof
		s:            if f.ok { math.sqrt(f.ssr / f64(f.dof)) } else { none }
		t:            if f.ok { f.t } else { none }
		residuals:    f.residuals
		section:      if f.ok { section(f, fits) } else { '' }
		iterations:   f.iterations
		evaluations:  f.evaluations
	}
}

// si_unit is a dimension's SI unit, '' for a pure number.
fn si_unit(d quantity.Dim) string {
	return if d.is_dimensionless() { '' } else { flowsheet.si_text(d) }
}

// load_failure is the report of a fit whose file didn't load.
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

// csv is a row for each residual: its row in the data, its column, the
// measured and the model's value, and the residual over its weight. Each
// number reads back as the same double.
pub fn (r Report) csv() string {
	mut w := csv.new_writer()
	w.write(['row', 'column', 'measured', 'model', 'weighted']) or {}
	for x in r.residuals {
		w.write([x.row.str(), x.column, x.measured.str(), x.model.str(), x.weighted.str()]) or {}
	}
	return w.str()
}

// table is the fit for a person to read: each parameter with its interval,
// their correlations, the residuals' spread, and the fitted constants as
// the file writes them.
pub fn (r Report) table() string {
	mut b := strings.new_builder(4096)
	if r.status != 'fitted' {
		b.writeln('${r.file}, fit ${r.fit}: no fit: ${r.reason}')
		return b.str()
	}
	s := r.s or { 0.0 }
	b.writeln('${r.file}, fit ${r.fit}: ${r.parameters.len} parameters to ${r.measurements} measurements in ${r.data}, in ${r.iterations} iterations; s = ${flowsheet.figures(s)} on ${r.dof} degrees of freedom.')
	b.writeln('')
	mut rows := [['parameter', 'value', 'standard error', '95% interval', 'unit']]
	for p in r.parameters {
		// In the file's own unit, as it was written.
		scale := file_scale(p.file_unit)
		v, e, lo, hi := p.value or { 0.0 }, p.standard_error or { 0.0 }, p.lower or { 0.0 }, p.upper or {
			0.0
		}
		rows << [p.name, flowsheet.figures(v / scale), flowsheet.figures(e / scale),
			'${flowsheet.figures(lo / scale)} to ${flowsheet.figures(hi / scale)}', p.file_unit]
	}
	write_table(mut b, rows)
	if r.parameters.len > 1 {
		b.writeln('')
		b.writeln('Correlations:')
		mut c_rows := [[''].clone()]
		c_rows[0] << r.parameters.map(it.name)
		for p, row in r.correlation {
			mut cells := [r.parameters[p].name]
			for q, v in row {
				rounded := math.round(v * 1000.0) / 1000.0
				cells << if q > p {
					''
				} else if rounded == 0.0 {
					'0'
				} else {
					'${rounded}'
				}
			}
			c_rows << cells
		}
		write_table(mut b, c_rows)
	}
	if r.centred.len > 0 {
		b.writeln('')
		b.writeln("Fitted centred on ${flowsheet.figures(r.centred[0].t_ref)} K, the data's mean T, and written back as the file gives them: ${r.centred.map(it.constant).join(', ')}.")
	}
	b.writeln('')
	b.writeln('Fitted, as the file writes them:')
	b.write_string(r.section)
	return b.str()
}

// write_table writes rows as columns, the first row their titles.
fn write_table(mut b strings.Builder, rows [][]string) {
	mut widths := []int{len: rows[0].len}
	for row in rows {
		for i, cell in row {
			widths[i] = math.max(widths[i], cell.runes().len)
		}
	}
	for row in rows {
		mut cells := []string{}
		for i, cell in row {
			cells << cell + ' '.repeat(widths[i] - cell.runes().len)
		}
		b.writeln(cells.join('  ').trim_right(' '))
	}
}

// file_scale is a file's unit's SI scale, 1 for a pure number.
fn file_scale(unit string) f64 {
	if unit == '' {
		return 1.0
	}
	u := quantity.unit(unit) or { return 1.0 }
	return u.scale
}

// figures10 is a value to 10 significant figures, as the shortest decimal
// that reads back to that double: by math.round, so that every compiler
// writes the same.
fn figures10(v f64) string {
	if v == 0.0 || !math.is_finite(v) {
		return v.str()
	}
	e := int(math.floor(math.log10(math.abs(v))))
	k := 9 - e
	r := if k >= 0 {
		math.round(v * math.pow(10, k)) / math.pow(10, k)
	} else {
		math.round(v / math.pow(10, -k)) * math.pow(10, -k)
	}
	text := r.str()
	return if text.ends_with('.0') { text[..text.len - 2] } else { text }
}

// in_unit is an SI value in the file's unit, with it, to 10 significant
// figures: "45480.71595 mol/(kg*s*bar^2)".
fn in_unit(v f64, unit string) string {
	if unit == '' {
		return figures10(v)
	}
	u := quantity.unit(unit) or { return figures10(v) }
	return '${figures10(v / u.scale)} ${unit}'
}

// value_of is a constant's fitted number in the file's unit, or the
// original where it wasn't fitted.
fn value_of(f Fitted, constant string, part flowsheet.Part, original string) string {
	for q, p in f.plan.parameters {
		if p.constant == constant && p.part == part {
			return in_unit(f.values[q], p.unit)
		}
	}
	return original
}

// section is the fitted constants as the file writes them: each constant's
// table, in the file's units, its other number as it was where only one
// was fitted; or each pair's kij.
fn section(f Fitted, fits flowsheet.FitFile) string {
	mut b := strings.new_builder(1024)
	plan := f.plan
	if plan.kind == .vle {
		for q, p in plan.parameters {
			b.writeln('${p.name} = ${figures10(f.values[q])}')
		}
		return b.str()
	}
	k := fits.kinetics or { return '' }
	mut done := []string{}
	for p in plan.parameters {
		if p.constant in done {
			continue
		}
		done << p.constant
		for c in k.constants() {
			if c.name != p.constant {
				continue
			}
			b.writeln('# ${c.name}, ${c.what}')
			match c.value {
				kinetics.Arrhenius {
					a := c.value as kinetics.Arrhenius
					centred := if a.t_ref != '' { ', T_ref = "${a.t_ref}"' } else { '' }
					key := if a.t_ref != '' { 'k_ref' } else { 'A' }
					b.writeln('{ ${key} = "${value_of(f, c.name, .a, a.a)}", Ea = "${value_of(f, c.name, .ea, a.ea)}"${centred} }')
				}
				kinetics.VantHoff {
					v := c.value as kinetics.VantHoff
					centred := if v.t_ref != '' { ', T_ref = "${v.t_ref}"' } else { '' }
					key := if v.t_ref != '' { 'K_ref' } else { 'A' }
					b.writeln('{ ${key} = "${value_of(f, c.name, .a, v.a)}", dH = "${value_of(f, c.name, .dh, v.dh)}"${centred} }')
				}
				kinetics.Fixed, kinetics.Written {
					b.writeln('"${value_of(f, c.name, .value, '')}"')
				}
			}
		}
	}
	return b.str()
}
