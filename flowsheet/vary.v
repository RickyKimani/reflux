module flowsheet

import quantity
import toml.ast

// Input is a number the file gives a feed or a unit, by its key (M8's D4),
// which a spec, a study or an optimiser may vary: 'feed.syngas.T',
// 'unit.cooler.outlet_T', 'unit.reactor.volume',
// 'feed.syngas.molar_flows.methane', 'unit.split.fractions[0]'. Its value
// and dimension are as the feed's or unit's maker read them.
pub struct Input {
pub:
	key   string
	table string // 'feed' or 'unit'
	name  string // the feed's or unit's
	index int    // into the flowsheet's feeds or units
	path  string // the key within its table: 'outlet_T', 'fractions[0]'
	dim   quantity.Dim
	value f64    // the file's, in SI
	text  string // the file's, as written
	line  int
	// named is the dimension its maker read it as, if it's one of M4's.
	named ?quantity.Dimension
}

// input is the number at a key, which must be one the file gives and its
// feed's or unit's maker reads as a number (D4). Anything else is refused
// at its line, with the keys like it that can be varied.
pub fn (l Loaded) input(key string) !Input {
	table := key.all_before('.')
	if table !in ['feed', 'unit'] || !key.contains('.') {
		return file_error(0, key, 'a key to vary starts with feed. or unit., as unit.cooler.outlet_T',
			[])
	}
	rest := key.all_after('.')
	names := if table == 'feed' {
		l.sheet.feeds.map(it.name)
	} else {
		l.sheet.units.map(it.unit.name)
	}
	mut name := ''
	for n in names {
		if rest.starts_with(n + '.') && n.len > name.len {
			name = n
		}
	}
	if name == '' {
		what := if table == 'feed' { 'feed' } else { 'unit' }
		return file_error(0, key, 'no ${what} ${rest.all_before('.')}', quantity.closest(rest.all_before('.'),
			names).map('${table}.${it}.${rest.all_after('.')}'))
	}
	path := rest[name.len + 1..]
	index := names.index(name)
	tree := l.table_of(table, name)
	readings := l.readings(table, name, tree)!
	for r in readings {
		if r.key == key {
			v := value_at(tree, steps(path)) or { ast.Value(ast.Null{}) }
			return Input{
				key:   key
				table: table
				name:  name
				index: index
				path:  path
				dim:   r.dim
				named: r.named
				value: r.value
				text:  text_of(v)
				line:  line_of(v)
			}
		}
	}
	keys := readings.map(it.key)
	if v := value_at(tree, steps(path)) {
		return file_error(line_of(v), key, "it isn't a number the ${table} reads as one, so it can't be varied: only its numbers can be",
			quantity.closest(key, keys))
	}
	listed := if keys.len > 0 {
		': it gives ${keys.map(it.all_after('${table}.${name}.')).join(', ')}'
	} else {
		''
	}
	return file_error(0, key, '${table}.${name} has no ${path} in the file${listed}', quantity.closest(key,
		keys))
}

// with is the file with these inputs at these values, in SI, each feed and
// unit they touch made again by its maker, which checks the values as it
// checks the file (D4). The rest is shared, so remaking costs a maker's
// call and a flash for a feed. Each value is checked to be what the maker
// read, to the bit.
pub fn (l Loaded) with(inputs []Input, values []f64) !Loaded {
	if inputs.len != values.len {
		return bad_spec('', 'got ${values.len} values for ${inputs.len} inputs')
	}
	mut root := l.root.clone()
	mut sheet := Flowsheet{
		units: l.sheet.units.clone()
		feeds: l.sheet.feeds.clone()
	}
	mut done := []string{}
	for i, input in inputs {
		whose := '${input.table}.${input.name}'
		if whose in done {
			continue
		}
		done << whose
		mut tree := ast.Value(l.table_of(input.table, input.name))
		for j in i .. inputs.len {
			other := inputs[j]
			if other.table == input.table && other.name == input.name {
				old := value_at(tree, steps(other.path)) or {
					return file_error(0, other.key, 'it is not in the file', [])
				}
				tree = replaced(tree, steps(other.path), written(old, values[j], other.dim))!
			}
		}
		table := tree as map[string]ast.Value
		mut s := l.section(input.table, input.name, table)
		if input.table == 'feed' {
			stream := read_stream(mut s, l.fluid) or { return s.located(err) }
			sheet.feeds[input.index] = Feed{input.name, stream}
		} else {
			u, _, _ := make_unit(mut s, input.name, l.makers)!
			sheet.units[input.index] = Placed{
				...sheet.units[input.index]
				unit: u
			}
		}
		for j in i .. inputs.len {
			other := inputs[j]
			if other.table == input.table && other.name == input.name {
				got := s.log.list.filter(it.key == other.key)
				if got.len == 0 || got.any(it.value != values[j]) {
					return bad_spec(input.name, '${other.key} was made ${got.map(it.value)} from ${values[j]}')
				}
			}
		}
		mut all := (root[input.table] or { ast.Value(map[string]ast.Value{}) }) as map[string]ast.Value
		all = all.clone()
		all[input.name] = ast.Value(table)
		root[input.table] = ast.Value(all)
	}
	return Loaded{
		...l
		root:  root
		sheet: sheet
	}
}

// table_of is a feed's or unit's table in the file.
fn (l Loaded) table_of(table string, name string) map[string]ast.Value {
	all := l.root[table] or { return map[string]ast.Value{} }
	if all is map[string]ast.Value {
		t := all[name] or { return map[string]ast.Value{} }
		if t is map[string]ast.Value {
			return t
		}
	}
	return map[string]ast.Value{}
}

// section is a feed's or unit's table as its maker reads it.
fn (l Loaded) section(table string, name string, tree map[string]ast.Value) Section {
	return Section{
		path:       '${table}.${name}'
		text:       l.text
		table:      tree
		name:       if table == 'unit' { name } else { '' }
		components: l.names
		context:    l.context
		log:        &Readings{}
	}
}

// readings are the numbers a feed's or unit's maker reads from a table.
fn (l Loaded) readings(table string, name string, tree map[string]ast.Value) ![]Reading {
	mut s := l.section(table, name, tree)
	if table == 'feed' {
		read_stream(mut s, l.fluid) or { return s.located(err) }
	} else {
		make_unit(mut s, name, l.makers)!
	}
	return s.log.list
}

// steps splits a key within a table: 'molar_flows.methane' is two steps,
// 'fractions[0]' one with an index.
fn steps(path string) []string {
	mut out := []string{}
	for part in path.split('.') {
		if part.ends_with(']') && part.contains('[') {
			out << part.all_before('[')
			out << '[' + part.all_after('[')
		} else {
			out << part
		}
	}
	return out
}

// index_of is a step's list index, if it is one.
fn index_of(step string) ?int {
	if step.starts_with('[') && step.ends_with(']') {
		inside := step[1..step.len - 1]
		if inside.len > 0 && inside.bytes().all(it.is_digit()) {
			return inside.int()
		}
	}
	return none
}

// value_at is the value at a path in a table.
fn value_at(v ast.Value, path []string) ?ast.Value {
	if path.len == 0 {
		return v
	}
	if i := index_of(path[0]) {
		if v is []ast.Value {
			if i < v.len {
				return value_at(v[i], path[1..])
			}
		}
		return none
	}
	if v is map[string]ast.Value {
		next := v[path[0]] or { return none }
		return value_at(next, path[1..])
	}
	return none
}

// replaced is a copy of a table with the value at a path replaced; only the
// tables and lists on the path are copied.
fn replaced(v ast.Value, path []string, new ast.Value) !ast.Value {
	if path.len == 0 {
		return new
	}
	if i := index_of(path[0]) {
		if v is []ast.Value {
			mut list := v.clone()
			list[i] = replaced(v[i], path[1..], new)!
			return ast.Value(list)
		}
		return error('no list at ${path[0]}')
	}
	if v is map[string]ast.Value {
		mut t := v.clone()
		t[path[0]] = replaced(v[path[0]] or { return error('no ${path[0]}') }, path[1..],
			new)!
		return ast.Value(t)
	}
	return error('no table at ${path[0]}')
}

// written is a value as the file would give it, where the old one was: a
// pure number as a number, anything else with its SI unit, so that the
// maker reads it back to the bit.
fn written(old ast.Value, value f64, d quantity.Dim) ast.Value {
	if d.is_dimensionless() {
		return match old {
			ast.Quoted, ast.Number {
				ast.Value(ast.Number{
					pos:  old.pos
					text: value.str()
				})
			}
			else {
				ast.Value(ast.Number{
					text: value.str()
				})
			}
		}
	}
	text := '${value.str()} ${si_text(d)}'
	return match old {
		ast.Quoted, ast.Number {
			ast.Value(ast.Quoted{
				pos:   old.pos
				text:  text
				quote: `"`
			})
		}
		else {
			ast.Value(ast.Quoted{
				text:  text
				quote: `"`
			})
		}
	}
}

// si_text is a dimension's unit in SI, as a file writes it, scale 1: a
// named dimension's own (K, Pa, mol/s, W, J/mol), or else its base units,
// kg/m^3.
pub fn si_text(d quantity.Dim) string {
	named := {
		'K':     quantity.Dimension.temperature
		'Pa':    .pressure
		'mol/s': .molar_flow
		'kg/s':  .mass_flow
		'm^3/s': .volume_flow
		'W':     .power
		'J/mol': .molar_enthalpy
		'J/kg':  .specific_enthalpy
	}
	for text, n in named {
		if d.eq(n.dim()) {
			return text
		}
	}
	symbols := ['kg', 'm', 's', 'K', 'mol', 'A']
	mut top := []string{}
	mut bottom := []string{}
	for i, e in d.exponents() {
		if e.num == 0 {
			continue
		}
		mag := if e.num < 0 { e.neg() } else { e }
		part := if mag.den != 1 {
			'${symbols[i]}^(${mag})'
		} else if mag.num == 1 {
			symbols[i]
		} else {
			'${symbols[i]}^${mag.num}'
		}
		if e.num > 0 {
			top << part
		} else {
			bottom << part
		}
	}
	num := if top.len == 0 { '1' } else { top.join('*') }
	if bottom.len == 0 {
		return num
	}
	return if bottom.len == 1 { '${num}/${bottom[0]}' } else { '${num}/(${bottom.join('*')})' }
}

// text_of is a value as the file writes it, for a message.
fn text_of(v ast.Value) string {
	return match v {
		ast.Quoted, ast.Number { v.text }
		else { '' }
	}
}

// readings_of is every number a feed's or unit's maker reads from the file,
// each at its key: what can be varied there.
pub fn (l Loaded) readings_of(table string, name string) []Reading {
	return l.readings(table, name, l.table_of(table, name)) or { []Reading{} }
}
