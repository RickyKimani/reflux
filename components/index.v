module components

import toml

// IndexEntry is what index.toml says about one compound: enough to find it
// by CAS number, name, alias or InChIKey without parsing its file.
pub struct IndexEntry {
pub:
	cas      string
	name     string
	inchikey ?string
	aliases  []string
}

// IndexEntry.of takes a record's index entry.
pub fn IndexEntry.of(r Record) IndexEntry {
	return IndexEntry{
		cas:      r.ids.cas
		name:     r.name
		inchikey: r.ids.inchikey
		aliases:  r.ids.aliases
	}
}

// write_index writes index.toml for a set of compounds, sorted by CAS
// number, so the file doesn't change when the source reorders its records.
pub fn write_index(entries []IndexEntry) !string {
	mut sorted := entries.clone()
	sorted.sort(a.cas < b.cas)
	mut out := ['format = ${format_version}']
	for e in sorted {
		out << ''
		out << '[compounds.${quote(e.cas)!}]'
		out << 'name = ${quote(e.name)!}'
		if k := e.inchikey {
			out << 'inchikey = ${quote(k)!}'
		}
		if e.aliases.len > 0 {
			mut quoted := []string{}
			for a in e.aliases {
				quoted << quote(a)!
			}
			out << 'aliases = [${quoted.join(', ')}]'
		}
	}
	return out.join('\n') + '\n'
}

// parse_index reads index.toml, as strictly as parse_toml reads a component
// file.
pub fn parse_index(text string) ![]IndexEntry {
	doc := toml.parse_text(text) or { return error('index.toml: ${err}') }
	top := doc.to_any().as_map()
	for k, _ in top {
		if k !in ['format', 'compounds'] {
			return error('index.toml: unknown key `${k}`')
		}
	}
	format := want_int(top['format'] or { return error('index.toml: no `format`') }, 'format')!
	if format != format_version {
		return error('index.toml: format ${format}; this reflux reads format ${format_version}')
	}
	mut entries := []IndexEntry{}
	compounds := want_table(top['compounds'] or { toml.Any(map[string]toml.Any{}) }, 'compounds')!
	for cas, v in compounds {
		what := 'compounds."${cas}"'
		t := want_table(v, what)!
		mut name := ''
		mut inchikey := ?string(none)
		mut aliases := []string{}
		for k, x in t {
			match k {
				'name' {
					name = want_string(x, '${what}.name')!
				}
				'inchikey' {
					inchikey = want_string(x, '${what}.inchikey')!
				}
				'aliases' {
					for a in want_array(x, '${what}.aliases')! {
						aliases << want_string(a, '${what}.aliases')!
					}
				}
				else {
					return error('index.toml: unknown key `${what}.${k}`')
				}
			}
		}
		if name == '' {
			return error('index.toml: no `${what}.name`')
		}
		entries << IndexEntry{
			cas:      cas
			name:     name
			inchikey: inchikey
			aliases:  aliases
		}
	}
	return entries
}
