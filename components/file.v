module components

import math
import os
import toml

// The component-file format this code reads and writes (docs/data-format.md).
pub const format_version = 1

// Integers up to 2⁵³ convert to f64 exactly. A larger one in a float field is
// refused: among other things, V's toml parser reads `inf` as one.
const max_exact_int = i64(9007199254740992)

// ---------------------------------------------------------------- writing

// write_toml writes a record as a component file. The same record always
// gives the same text: keys in the order of the struct fields, floats in
// V's shortest round-trip form, LF line endings. It fails on values TOML
// can't carry (NaN, infinities, control characters in strings) and on an
// override in `sources` that names no value.
pub fn write_toml(r Record) !string {
	mut used := map[string]bool{}
	mut out := []string{}
	out << 'format = ${format_version}'
	out << 'name = ${quote(r.name)!}'
	out << 'source = ${quote(r.source)!}'
	out << section('ids', ids_lines(r.ids)!)
	out << section('constants', scalar_lines(r.constants, 'constants', r.sources, mut used)!)
	out << section('parameters', scalar_lines(r.parameters, 'parameters', r.sources, mut
		used)!)
	$for f in Correlations.fields {
		$if f.typ is ?CorrelationData {
			key := 'correlations.${f.name}'
			// Unwrapped in a typed function: inside `$for`, `if d := …` on an
			// Option struct field leaves d an Option (a V 0.5.2 bug).
			lines := correlation_lines(r.correlations.$(f.name), r.sources[key] or { '' })!
			if lines.len > 0 {
				out << section(key, lines)
				used[key] = key in r.sources
			}
		}
	}
	$for f in Groups.fields {
		$if f.typ is map[int]int {
			counts := r.groups.$(f.name).clone()
			out << section('groups.${f.name}', counts.keys().map('"${it}" = ${counts[it]}'))
		}
	}
	for key, _ in r.sources {
		if !used[key] {
			return error('${r.ids.cas}: a source for `${key}`, which has no value')
		}
	}
	return out.filter(it != '').join('\n') + '\n'
}

// correlation_lines writes one correlation's table, or nothing if it's
// absent; src is its source if it differs from the file's, or ''.
fn correlation_lines(opt ?CorrelationData, src string) ![]string {
	d := opt or { return [] }
	mut lines := ['eq = ${d.eq}', 'coeffs = ${float_array(d.coeffs)!}',
		't_range = ${float_array([d.tmin, d.tmax])!}', 'units = ${quote(d.units)!}']
	if tc := d.tc {
		lines << 'tc = ${float(tc)!}'
	}
	if src != '' {
		lines << 'source = ${quote(src)!}'
	}
	return lines
}

// section is a table header and its lines after a blank line, or nothing if
// the table would be empty.
fn section(name string, lines []string) string {
	if lines.len == 0 {
		return ''
	}
	return '\n[${name}]\n' + lines.join('\n')
}

fn ids_lines(ids Ids) ![]string {
	mut lines := []string{}
	$for f in Ids.fields {
		$if f.typ is string {
			// `ids.$(f.name)` can't go inside an interpolation; V's parser
			// doesn't accept it there.
			v := ids.$(f.name)
			lines << '${f.name} = ${quote(v)!}'
		} $else $if f.typ is ?string {
			if v := ids.$(f.name) {
				lines << '${f.name} = ${quote(v)!}'
			}
		} $else $if f.typ is ?int {
			if v := ids.$(f.name) {
				lines << '${f.name} = ${v}'
			}
		} $else $if f.typ is []string {
			v := ids.$(f.name)
			if v.len > 0 {
				mut quoted := []string{}
				for s in v {
					quoted << quote(s)!
				}
				lines << '${f.name} = [${quoted.join(', ')}]'
			}
		}
	}
	return lines
}

// scalar_lines writes the ?f64 and ?int fields of Constants or Parameters.
// A value whose source differs from the file's becomes an inline table.
fn scalar_lines[T](s T, table string, sources map[string]string, mut used map[string]bool) ![]string {
	mut lines := []string{}
	$for f in T.fields {
		$if f.typ is ?f64 {
			if v := s.$(f.name) {
				key := '${table}.${f.name}'
				if src := sources[key] {
					lines << '${f.name} = { value = ${float(v)!}, source = ${quote(src)!} }'
					used[key] = true
				} else {
					lines << '${f.name} = ${float(v)!}'
				}
			}
		} $else $if f.typ is ?int {
			if v := s.$(f.name) {
				lines << '${f.name} = ${v}'
			}
		}
	}
	return lines
}

// float prints x so that TOML reads it back as the same float: V's
// shortest round-trip form always has a point or an exponent ('2.0', not
// '2'), which TOML needs to see a float rather than an integer.
fn float(x f64) !string {
	if !math.is_finite(x) {
		return error('${x} cannot be written: TOML has no portable NaN or infinity')
	}
	return x.str()
}

fn float_array(xs []f64) !string {
	mut parts := []string{}
	for x in xs {
		parts << float(x)!
	}
	return '[${parts.join(', ')}]'
}

// quote writes s as a TOML basic string. SMILES use `\` for stereo bonds, so
// backslashes and quotes are escaped. Control characters are refused rather
// than escaped: nothing in the data should have one.
fn quote(s string) !string {
	for c in s {
		if c < 0x20 || c == 0x7f {
			return error('control character in `${s}`')
		}
	}
	return '"' + s.replace('\\', '\\\\').replace('"', '\\"') + '"'
}

// ---------------------------------------------------------------- reading

// read_toml reads a component file. See parse_toml.
pub fn read_toml(path string) !Record {
	text := os.read_file(path)!
	return parse_toml(text, path)!
}

// parse_toml reads the text of a component file, strictly. Errors name the
// file and the key:
// - a key the format doesn't have, so a typo can't silently drop a value;
// - a value of the wrong type (V's toml would quietly read a string as 0.0);
// - a missing required key (format, name, source, ids.cas; a correlation's
//   eq, coeffs, t_range and units);
// - a CAS number that isn't the file's name.
pub fn parse_toml(text string, file string) !Record {
	name := os.file_name(file)
	r := parse_doc(text) or { return error('${name}: ${err}') }
	if name != '${r.ids.cas}.toml' {
		return error('${name}: holds ${r.ids.cas}, so it should be named ${r.ids.cas}.toml')
	}
	return r
}

fn parse_doc(text string) !Record {
	doc := toml.parse_text(text)!
	top := doc.to_any().as_map()
	allowed := ['format', 'name', 'source', 'ids', 'constants', 'parameters', 'correlations', 'groups']
	for k, _ in top {
		if k !in allowed {
			return error('unknown key `${k}`')
		}
	}
	format := want_int(top['format'] or { return error('no `format` key') }, 'format')!
	if format != format_version {
		return error('format ${format}; this reflux reads format ${format_version}')
	}
	mut r := Record{
		name:   want_string(top['name'] or { return error('no `name`') }, 'name')!
		source: want_string(top['source'] or { return error('no `source`') }, 'source')!
	}
	read_fields(mut r.ids, want_table(top['ids'] or { return error('no [ids] table') }, 'ids')!,
		'ids', mut r.sources)!
	if r.ids.cas == '' {
		return error('no `ids.cas`')
	}
	if v := top['constants'] {
		read_fields(mut r.constants, want_table(v, 'constants')!, 'constants', mut r.sources)!
	}
	if v := top['parameters'] {
		read_fields(mut r.parameters, want_table(v, 'parameters')!, 'parameters', mut r.sources)!
	}
	if v := top['correlations'] {
		for key, t in want_table(v, 'correlations')! {
			d, src := read_correlation(want_table(t, 'correlations.${key}')!, 'correlations.${key}')!
			r.correlations.set(key, d) or { return error('unknown correlation `${key}`') }
			if s := src {
				r.sources['correlations.${key}'] = s
			}
		}
	}
	if v := top['groups'] {
		for key, t in want_table(v, 'groups')! {
			counts := read_groups(want_table(t, 'groups.${key}')!, 'groups.${key}')!
			r.groups.set(key, counts) or { return error('unknown group method `${key}`') }
		}
	}
	return r
}

// read_fields fills the fields of Ids, Constants or Parameters from a
// table, by name. A float field also takes an integer, or an inline table
// { value = …, source = … } whose source goes into `sources`.
fn read_fields[T](mut s T, table map[string]toml.Any, prefix string, mut sources map[string]string) ! {
	for key, v in table {
		what := '${prefix}.${key}'
		mut found := false
		$for f in T.fields {
			if f.name == key {
				found = true
				$if f.typ is string {
					s.$(f.name) = want_string(v, what)!
				} $else $if f.typ is ?string {
					s.$(f.name) = want_string(v, what)!
				} $else $if f.typ is ?int {
					s.$(f.name) = want_int(v, what)!
				} $else $if f.typ is []string {
					mut list := []string{}
					for item in want_array(v, what)! {
						list << want_string(item, '${what} item')!
					}
					s.$(f.name) = list
				} $else $if f.typ is ?f64 {
					if v is map[string]toml.Any {
						for k, _ in v {
							if k !in ['value', 'source'] {
								return error('${what}: unknown key `${k}`')
							}
						}
						s.$(f.name) = want_number(v['value'] or {
							return error('${what}: no `value`')
						}, what)!
						sources[what] = want_string(v['source'] or {
							return error('${what}: no `source`')
						}, '${what}.source')!
					} else {
						s.$(f.name) = want_number(v, what)!
					}
				}
			}
		}
		if !found {
			return error('unknown key `${what}`')
		}
	}
}

fn read_correlation(t map[string]toml.Any, what string) !(CorrelationData, ?string) {
	mut d := CorrelationData{}
	mut src := ?string(none)
	for k, v in t {
		match k {
			'eq' {
				d.eq = want_int(v, '${what}.eq')!
			}
			'coeffs' {
				for x in want_array(v, '${what}.coeffs')! {
					d.coeffs << want_number(x, '${what}.coeffs')!
				}
			}
			't_range' {
				r := want_array(v, '${what}.t_range')!
				if r.len != 2 {
					return error('${what}.t_range should be [Tmin, Tmax], not ${r.len} numbers')
				}
				d.tmin = want_number(r[0], '${what}.t_range')!
				d.tmax = want_number(r[1], '${what}.t_range')!
			}
			'units' {
				d.units = want_string(v, '${what}.units')!
			}
			'tc' {
				d.tc = want_number(v, '${what}.tc')!
			}
			'source' {
				src = want_string(v, '${what}.source')!
			}
			else {
				return error('unknown key `${what}.${k}`')
			}
		}
	}
	for need in ['eq', 'coeffs', 't_range', 'units'] {
		if need !in t {
			return error('no `${what}.${need}`')
		}
	}
	return d, src
}

// read_groups reads { "19" = 3, … }: group numbers as keys, written as
// plain decimals, with positive counts.
fn read_groups(t map[string]toml.Any, what string) !map[int]int {
	mut counts := map[int]int{}
	for k, v in t {
		id := k.int()
		if id.str() != k || id <= 0 {
			return error('${what}: `${k}` is not a group number')
		}
		n := want_int(v, '${what}."${k}"')!
		if n <= 0 {
			return error('${what}."${k}": a count must be positive, not ${n}')
		}
		counts[id] = n
	}
	return counts
}

// The typed getters. toml.Any's own conversions turn a value of the wrong
// type into a zero without complaint, so these check first.

fn want_string(v toml.Any, what string) !string {
	if v is string {
		return v
	}
	return error('`${what}` should be a string, not ${kind(v)}')
}

fn want_int(v toml.Any, what string) !int {
	n := integer_of(v) or { return error('`${what}` should be an integer, not ${kind(v)}') }
	if n < i64(min_i32) || n > i64(max_i32) {
		return error('`${what}` = ${n} is out of range')
	}
	return int(n)
}

fn want_number(v toml.Any, what string) !f64 {
	if v is f64 {
		// V 0.5.2's toml never hands back a non-finite f64 (it reads `nan` as
		// a string and `inf` as an integer), so no test can reach this. It
		// stays in case the parser is fixed to read them as floats.
		if !math.is_finite(v) {
			return error('`${what}` is not finite')
		}
		return v
	}
	n := integer_of(v) or { return error('`${what}` should be a number, not ${kind(v)}') }
	if n > max_exact_int || n < -max_exact_int {
		return error('`${what}` = ${n} is too large to be exact as a float')
	}
	return f64(n)
}

// integer_of returns an integer value as i64, whichever of V's integer
// variants toml chose (it uses more than one), or none for a non-integer.
fn integer_of(v toml.Any) ?i64 {
	return match v {
		i64 {
			v
		}
		int {
			i64(v)
		}
		u64 {
			if v > u64(max_i64) {
				// Beyond i64: certainly too large; report it as the largest.
				max_i64
			} else {
				i64(v)
			}
		}
		else {
			none
		}
	}
}

fn want_table(v toml.Any, what string) !map[string]toml.Any {
	if v is map[string]toml.Any {
		return v
	}
	return error('`${what}` should be a table, not ${kind(v)}')
}

fn want_array(v toml.Any, what string) ![]toml.Any {
	if v is []toml.Any {
		return v
	}
	return error('`${what}` should be an array, not ${kind(v)}')
}

// kind names a value's TOML type, for error messages. (toml.Any's
// type_name() says just "Any" for tables and arrays.)
fn kind(v toml.Any) string {
	return match v {
		string { 'a string ("${v}")' }
		i64, int, u64 { 'an integer' }
		f64, f32 { 'a float' }
		bool { 'a boolean' }
		[]toml.Any { 'an array' }
		map[string]toml.Any { 'a table' }
		else { 'a date or time' }
	}
}
