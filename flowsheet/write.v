module flowsheet

import strings
import toml.ast

// write is the file again, in a fixed order, with its numbers and units as
// they were given, so that loading it gives the same flowsheet (criterion 6).
// Its comments are gone: V's parser doesn't tie them to anything.
pub fn (l Loaded) write() string {
	mut b := strings.new_builder(1024)
	b.write_string('format = ${value_text(l.root['format'] or { ast.Null{} })}\n')
	for name in ['components', 'model', 'solver'] {
		if v := l.root[name] {
			write_table(mut b, name, v as map[string]ast.Value)
		}
	}
	for name in ['feed', 'estimate', 'unit'] {
		if v := l.root[name] {
			for key, t in v as map[string]ast.Value {
				write_table(mut b, '${name}.${key_text(key)}', t as map[string]ast.Value)
			}
		}
	}
	return b.str()
}

// write_table writes a table under its header. Anything nested goes inline:
// nothing in the format nests deeper than one line can hold.
fn write_table(mut b strings.Builder, header string, table map[string]ast.Value) {
	b.write_string('\n[${header}]\n')
	for key, v in table {
		b.write_string('${key_text(key)} = ${value_text(v)}\n')
	}
}

// basic_string is text as a TOML basic string. V's parser hands strings over
// decoded, so what it undid has to be done again.
fn basic_string(text string) string {
	mut b := strings.new_builder(text.len + 2)
	b.write_u8(`"`)
	for c in text {
		match c {
			`\\` { b.write_string(r'\\') }
			`"` { b.write_string(r'\"') }
			`\n` { b.write_string(r'\n') }
			`\r` { b.write_string(r'\r') }
			`\t` { b.write_string(r'\t') }
			`\f` { b.write_string(r'\f') }
			`\b` { b.write_string(r'\b') }
			else { b.write_u8(c) }
		}
	}
	b.write_u8(`"`)
	return b.str()
}

// key_text is a key, quoted if it has to be.
fn key_text(key string) string {
	bare := key != '' && key.bytes().all(it.is_letter() || it.is_digit() || it == `_` || it == `-`)
	return if bare { key } else { basic_string(key) }
}

// value_text is a value as TOML.
fn value_text(v ast.Value) string {
	return match v {
		ast.Quoted {
			if v.quote == `'` && !v.is_multiline {
				"'${v.text}'"
			} else {
				basic_string(v.text)
			}
		}
		ast.Number, ast.Bool {
			v.text
		}
		[]ast.Value {
			'[' + v.map(value_text(it)).join(', ') + ']'
		}
		map[string]ast.Value {
			mut parts := []string{cap: v.len}
			for key, x in v {
				parts << '${key_text(key)} = ${value_text(x)}'
			}
			'{ ' + parts.join(', ') + ' }'
		}
		else {
			''
		}
	}
}
