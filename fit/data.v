module fit

import encoding.csv
import math
import os
import quantity
import strconv

// Data is a fit's measurements (M8's D14): a CSV file whose header names
// each column with its unit, `T [K]`, `P [bar]`, `y[CH4]` or
// `r[reforming] [mol/(kg*s)]`, and a row for each state, every cell in SI
// (NaN where it's empty: not measured).
pub struct Data {
pub:
	file    string
	columns []Column
	rows    [][]f64
}

// Column is a column as the header gives it: its name, `T` or `y[CH4]`,
// and its unit, '' for a pure number.
pub struct Column {
pub:
	name string
	unit string
}

// DataError is data that can't be read: the file, and where in it.
pub struct DataError {
	Error
pub:
	file   string
	line   int // from 1, the header's 1; 0 for the file as a whole
	column string
	reason string
}

// msg is the reason, with the file, its line and its column.
pub fn (e DataError) msg() string {
	at := if e.line > 0 { ':${e.line}' } else { '' }
	col := if e.column != '' { ', ${e.column}' } else { '' }
	return '${e.file}${at}${col}: ${e.reason}'
}

fn data_fail(file string, line int, column string, reason string) IError {
	return DataError{
		file:   file
		line:   line
		column: column
		reason: reason
	}
}

// read_data reads a fit's CSV file: its header's columns, and each row as
// text, converted to SI by kind_of each column. Blank lines are skipped.
pub fn read_data(path string) !Data {
	text := os.read_file(path) or { return data_fail(path, 0, '', "it can't be read: ${err.msg()}") }
	mut reader := csv.new_reader(text)
	mut header := []string{}
	mut cells := [][]string{}
	for {
		row := reader.read() or { break }
		if row.len == 1 && row[0].trim_space() == '' {
			continue
		}
		if header.len == 0 {
			header = row.map(it.trim_space())
			continue
		}
		cells << row.map(it.trim_space())
	}
	if header.len == 0 {
		return data_fail(path, 0, '', 'it has no header')
	}
	if cells.len == 0 {
		return data_fail(path, 0, '', 'it has no rows')
	}
	mut columns := []Column{}
	for h in header {
		c := column_of(h)
		if c.name == '' {
			return data_fail(path, 1, h, 'a column needs a name, as `T [K]`')
		}
		if columns.any(it.name == c.name) {
			return data_fail(path, 1, h, '${c.name} is a column twice')
		}
		columns << c
	}
	mut rows := [][]f64{}
	for r, row in cells {
		if row.len != columns.len {
			return data_fail(path, r + 2, '', 'it has ${row.len} cells, not ${columns.len}')
		}
		mut values := []f64{len: columns.len}
		for k, cell in row {
			values[k] = si(cell, columns[k]) or { return data_fail(path, r + 2, header[k], err.msg()) }
		}
		rows << values
	}
	return Data{
		file:    path
		columns: columns
		rows:    rows
	}
}

// column_of splits a header's cell into its name and unit: `T [K]` is T in
// K, `y[CH4]` a pure number, `r[a] [mol/(kg*s)]` r[a] in mol/(kg·s).
fn column_of(h string) Column {
	if h.ends_with(']') && h.contains(' [') {
		return Column{
			name: h.all_before_last(' [').trim_space()
			unit: h.all_after_last(' [').all_before_last(']').trim_space()
		}
	}
	return Column{
		name: h
	}
}

// si is a cell in SI: T and P by quantity's own reading, which knows °C
// and bar(g) and refuses an absolute value at or below zero; a fraction
// as it is, or in %; anything else times its unit's scale. An empty cell is
// NaN, not measured.
fn si(cell string, c Column) !f64 {
	if cell == '' {
		return math.nan()
	}
	v := strconv.atof64(cell) or { return error('"${cell}" is not a number') }
	if !math.is_finite(v) {
		return error('"${cell}" is not a finite number')
	}
	if c.name == 'T' {
		return quantity.parse('${cell} ${c.unit}', .temperature)!
	}
	if c.name == 'P' {
		return quantity.parse('${cell} ${c.unit}', .pressure)!
	}
	if c.unit == '' {
		return v
	}
	u := quantity.unit(c.unit)!
	return v * u.scale
}

// index is a column's place, or -1.
pub fn (d &Data) index(name string) int {
	for k, c in d.columns {
		if c.name == name {
			return k
		}
	}
	return -1
}

// dim is a column's unit's dimension.
fn (c Column) dim() !quantity.Dim {
	if c.unit == '' {
		return quantity.Dim{}
	}
	return quantity.unit(c.unit)!.dim
}
