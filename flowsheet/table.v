module flowsheet

import math
import strings

// How wide a stream table may grow before its other streams go into another
// below it.
const table_width = 100

// table is the report for a person to read: how it was solved, a stream table
// in engineering units, and each unit's duty or work.
pub fn (r Report) table() string {
	mut b := strings.new_builder(4096)
	r.summary(mut b)
	r.stream_table(mut b)
	r.unit_table(mut b)
	return b.str()
}

// trace is every pass of every loop, each tear's residual (D7) in a column.
pub fn (r Report) trace() string {
	mut b := strings.new_builder(1024)
	loops := r.blocks.filter(it.tears.len > 0)
	if loops.len == 0 {
		b.writeln('No loops, so nothing was iterated.')
	}
	for blk in loops {
		b.writeln('The loop of ${blk.units.join(', ')}, torn at ${blk.tears.join(', ')}:')
		widths := blk.tears.map(math.max(10, it.runes().len))
		mut head := '  ${'pass':6}'
		for i, name in blk.tears {
			head += '  ${pad_left(name, widths[i])}'
		}
		b.writeln(head)
		for k, row in blk.residuals {
			mut line := '  ${k + 1:6}'
			for i, x in row {
				line += '  ${pad_left(scientific(x, 3), widths[i])}'
			}
			b.writeln(line)
		}
		if blk.residuals.len == 0 {
			b.writeln('  (no pass finished)')
		}
	}
	return b.str()
}

fn method_name(m Method) string {
	return match m {
		.direct { 'direct substitution' }
		.wegstein { 'Wegstein' }
		.broyden { 'Broyden' }
	}
}

fn (r Report) summary(mut b strings.Builder) {
	what := if r.file == '' { 'The flowsheet' } else { r.file }
	solver := r.solver or { SolverReport{} }
	loops := r.blocks.filter(it.tears.len > 0)
	if loops.len == 0 {
		b.writeln('${what}: solved, with no loops, so each unit once.')
	} else {
		b.writeln("${what}: solved by ${method_name(solver.method)}, to ${scientific(solver.tol, 2)} of each tear variable's scale.")
	}
	for blk in loops {
		b.writeln('  The loop of ${blk.units.join(', ')}: torn at ${blk.tears.join(', ')}, ${blk.residuals.len} passes.')
	}
	for note in r.notes {
		b.writeln('  ${note}')
	}
	if all := r.balance {
		worst := r.worst_unit or { Balance{} }
		b.writeln('  The whole closes its balances to ${scientific(all.component, 2)} of the flow and ${scientific(all.energy, 2)} of F·RT,')
		b.writeln('  and the worst unit to ${scientific(worst.component, 2)} and ${scientific(worst.energy, 2)}.')
	}
}

// Column is a stream's column of the stream table, top to bottom.
struct Column {
	cells []string
	width int
}

fn (r Report) stream_table(mut b strings.Builder) {
	mut labels := ['Stream', 'From', 'To', 'T (K)', 'P (bar)', 'Vapour fraction', 'Flow (mol/s)',
		'Flow (kg/s)', 'H (kJ/mol)', 'S (J/mol/K)', 'Mole fractions']
	for c in r.components {
		labels << '  ${c.name}'
	}
	mut label_width := 0
	for l in labels {
		label_width = math.max(label_width, l.runes().len)
	}
	mut columns := []Column{}
	for s in r.streams {
		mut cells := [s.name, s.from or { '-' }, s.to or { '-' }, figures(s.t), figures(s.p / 1e5),
			figures(s.vapour_fraction), figures(s.molar_flow), figures(s.mass_flow),
			figures(s.h / 1e3), figures(s.s), '']
		for x in s.mole_fractions {
			cells << figures(x)
		}
		mut width := 10
		for c in cells {
			width = math.max(width, c.runes().len)
		}
		columns << Column{cells, width}
	}
	// As many streams to a table as fit, and at least one.
	mut start := 0
	for start < columns.len {
		mut end_at := start + 1
		mut used := label_width + 2 + columns[start].width
		for end_at < columns.len && used + 2 + columns[end_at].width <= table_width {
			used += 2 + columns[end_at].width
			end_at++
		}
		b.writeln('')
		for i, l in labels {
			mut line := pad_right(l, label_width)
			for c in columns[start..end_at] {
				line += '  ${pad_left(c.cells[i], c.width)}'
			}
			b.writeln(line.trim_right(' '))
		}
		start = end_at
	}
}

fn (r Report) unit_table(mut b strings.Builder) {
	mut kinds := []string{} // of energy, in the order they first appear
	for u in r.units {
		for e in u.energy {
			if e.name !in kinds {
				kinds << e.name
			}
		}
	}
	heads := kinds.map('${it.capitalize()} (kW)')
	mut name_width, mut type_width := 'Unit'.len, 'Type'.len
	mut widths := heads.map(math.max(10, it.len))
	// Each unit's cell under each kind, and each column as wide as its widest:
	// an adiabatic drum's duty can be a rounding's -9.3132e-13 kW.
	mut cells := [][]string{}
	for u in r.units {
		name_width = math.max(name_width, u.name.runes().len)
		type_width = math.max(type_width, u.kind.runes().len)
		mut row := []string{len: kinds.len}
		for e in u.energy {
			row[kinds.index(e.name)] = figures(e.rate / 1e3)
		}
		for i, cell in row {
			widths[i] = math.max(widths[i], cell.runes().len)
		}
		cells << row
	}
	b.writeln('')
	mut head := '${pad_right('Unit', name_width)}  ${pad_right('Type', type_width)}'
	for i, h in heads {
		head += '  ${pad_left(h, widths[i])}'
	}
	b.writeln(head.trim_right(' '))
	for k, u in r.units {
		mut line := '${pad_right(u.name, name_width)}  ${pad_right(u.kind, type_width)}'
		for i, cell in cells[k] {
			line += '  ${pad_left(cell, widths[i])}'
		}
		b.writeln(line.trim_right(' '))
	}
	warned := r.units.filter(it.warnings.len > 0)
	if warned.len > 0 {
		b.writeln('')
		b.writeln('Warnings:')
		for u in warned {
			for w in u.warnings {
				b.writeln('  ${u.name}: ${w}')
			}
		}
	}
}

fn pad_left(s string, width int) string {
	n := s.runes().len
	return if n >= width { s } else { ' '.repeat(width - n) + s }
}

fn pad_right(s string, width int) string {
	n := s.runes().len
	return if n >= width { s } else { s + ' '.repeat(width - n) }
}

// figures is x to 5 significant figures: as a decimal from a thousandth up to a
// million, and as a power of ten outside that.
fn figures(x f64) string {
	if x == 0.0 || !math.is_finite(x) {
		return if x == 0.0 { '0' } else { '${x}' }
	}
	d, e := digits(x, 5)
	sign := if x < 0.0 { '-' } else { '' }
	return match true {
		e < -3 || e > 5 { '${sign}${d[..1]}.${d[1..]}${exponent(e)}' }
		e >= 4 { sign + d + '0'.repeat(e - 4) }
		e >= 0 { '${sign}${d[..e + 1]}.${d[e + 1..]}' }
		else { '${sign}0.${'0'.repeat(-e - 1)}${d}' }
	}
}

// scientific is x to n significant figures as a power of ten, such as 9.5e-14.
fn scientific(x f64, n int) string {
	if x == 0.0 || !math.is_finite(x) {
		return if x == 0.0 { '0' } else { '${x}' }
	}
	d, e := digits(x, n)
	sign := if x < 0.0 { '-' } else { '' }
	rest := if n > 1 { '.${d[1..]}' } else { '' }
	return '${sign}${d[..1]}${rest}${exponent(e)}'
}

fn exponent(e int) string {
	return if e < 0 { 'e-${-e:02}' } else { 'e+${e:02}' }
}

// digits is |x| to n significant figures: those n digits, and the power of ten
// of the first. It rounds by math.round, half away from zero, and prints only
// the integer, because printf's rounding of a tie differs between tcc's C
// library and gcc's (docs/v-for-go.md, §14).
fn digits(x f64, n int) (string, int) {
	a := math.abs(x)
	mut e := int(math.floor(math.log10(a)))
	mut m := scaled(a, n - 1 - e)
	// Rounding can carry into a new digit, and log10 can fall short of a power of
	// ten: one more either way. Where it rounds up just short of one, so does m.
	if m >= math.pow(10.0, n) {
		e++
		m = scaled(a, n - 1 - e)
	}
	return '${i64(m)}', e
}

// scaled is a·10^k rounded to a whole number. For k below 0 it divides, since
// 10^k isn't exact then: 846064.9999999999 · 0.1 is 84606.5, which rounds up,
// but divided by 10 it's 84606.49999999999.
fn scaled(a f64, k int) f64 {
	return if k >= 0 {
		math.round(a * math.pow(10.0, k))
	} else {
		math.round(a / math.pow(10.0, -k))
	}
}
