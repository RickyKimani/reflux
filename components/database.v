module components

import os
import strings

// default_root is the data reflux comes with: this repository's data/.
// @VMODROOT is the directory holding v.mod, filled in when reflux is
// compiled, so this works for anything built from the repository, the C
// API's library among them. A data root holds:
// - components/, the compounds: ChemSep converted, one file each;
// - vendor/chemsep/, ChemSep's kᵢⱼ banks (pr.ipd, srk.ipd) and
//   Mathias–Copeman constants (mc-PR.lib).
// Set REFLUX_DATA to a root elsewhere, a copy of this one for instance.
pub const default_root = os.join_path(@VMODROOT, 'data')

// default_dir is the default root's compounds.
pub const default_dir = os.join_path(default_root, 'components')

// Database finds compounds in one or more directories of component files,
// and holds where its root's kᵢⱼ banks and Mathias–Copeman constants are.
// A compound whose CAS number appears in a later directory replaces the
// earlier one entirely, which is how a directory of your own overrides the
// built-in data. Files are parsed only when a compound is asked for.
pub struct Database {
pub:
	root   string   // the data root
	dirs   []string // the root's components/, then any others, in order
	vendor string   // the root's vendor/chemsep/
mut:
	entries map[string]Entry // CAS number → the latest entry
}

struct Entry {
	IndexEntry
	path string // the component file
}

// data_root is where a database opens unless told otherwise: $REFLUX_DATA if
// that's set, and default_root if not.
pub fn data_root() string {
	return os.getenv_opt('REFLUX_DATA') or { default_root }
}

// load_default opens data_root().
pub fn load_default() !Database {
	return open(data_root(), [])
}

// open opens a data root: its compounds, then those in any extra
// directories, which replace the root's compounds of the same CAS number;
// and its kᵢⱼ banks and Mathias–Copeman constants, so that a copy of the data
// moved elsewhere is used whole.
pub fn open(root string, extra []string) !Database {
	comps := os.join_path(root, 'components')
	vendor := os.join_path(root, 'vendor', 'chemsep')
	if !os.is_dir(comps) {
		hint := if os.exists(os.join_path(root, 'index.toml')) {
			'; it looks like a components directory itself, so its parent is the root'
		} else {
			''
		}
		return error('${root} is not a data root: it has no components directory${hint}')
	}
	if !os.is_dir(vendor) {
		return error('${root} is not a data root: it has no vendor/chemsep directory')
	}
	mut db := load_dirs([comps, ...extra])!
	return Database{
		...db
		root:   root
		vendor: vendor
	}
}

// load_dirs reads directories of component files, in order. A directory with
// an index.toml is read through it, as the built-in one is. A directory
// without one, such as a directory of your own, has every *.toml file in it
// read, which is quick for a handful of files.
fn load_dirs(dirs []string) !Database {
	mut db := Database{
		dirs: dirs
	}
	for dir in dirs {
		index := os.join_path(dir, 'index.toml')
		if os.exists(index) {
			for e in parse_index(os.read_file(index)!)! {
				path := os.join_path(dir, '${e.cas}.toml')
				if !os.exists(path) {
					return error('${index} lists ${e.cas}, but ${path} is missing')
				}
				db.entries[e.cas] = Entry{
					IndexEntry: e
					path:       path
				}
			}
		} else {
			mut names := os.ls(dir) or { return error('cannot read component directory ${dir}: ${err}') }
			names.sort()
			for name in names.filter(it.ends_with('.toml')) {
				path := os.join_path(dir, name)
				r := read_toml(path)!
				db.entries[r.ids.cas] = Entry{
					IndexEntry: IndexEntry.of(r)
					path:       path
				}
			}
		}
	}
	return db
}

// count returns how many compounds the database holds.
pub fn (db Database) count() int {
	return db.entries.len
}

// list returns every compound's index entry, sorted by CAS number.
pub fn (db Database) list() []IndexEntry {
	mut out := []IndexEntry{}
	for _, e in db.entries {
		out << e.IndexEntry
	}
	out.sort(a.cas < b.cas)
	return out
}

// LookupError is find's error. It carries the closest matches, so a caller
// can offer them; its message lists them too.
pub struct LookupError {
	Error
pub:
	query       string
	suggestions []string // 'name (CAS)', closest first
	ambiguous   bool     // the query matched more than one compound
}

// msg makes LookupError an IError, with the suggestions in the message.
pub fn (e LookupError) msg() string {
	if e.ambiguous {
		return '`${e.query}` matches ${e.suggestions.join(' and ')}; ask for one by its CAS number'
	}
	if e.suggestions.len == 0 {
		return 'no compound `${e.query}` in the database'
	}
	return 'no compound `${e.query}` in the database. Did you mean ${e.suggestions.join(', ')}?'
}

// find returns the compound with this CAS number, InChIKey, name or alias.
// Names and aliases are matched ignoring case and surrounding spaces. A
// query that matches nothing, or more than one compound, is a LookupError.
pub fn (db Database) find(query string) !Record {
	q := query.trim_space()
	if q == '' {
		// Not ambiguous: a compound without an InChIKey would match it below.
		return LookupError{
			query: q
		}
	}
	if e := db.entries[q] {
		return read_toml(e.path)!
	}
	lower := q.to_lower()
	mut matches := []Entry{}
	for _, e in db.entries {
		if e.name.to_lower() == lower || e.aliases.any(it.to_lower() == lower)
			|| (e.inchikey or { '' }) == q.to_upper() {
			matches << e
		}
	}
	if matches.len == 1 {
		return read_toml(matches[0].path)!
	}
	if matches.len > 1 {
		matches.sort(a.cas < b.cas)
		return LookupError{
			query:       q
			suggestions: matches.map('${it.name} (${it.cas})')
			ambiguous:   true
		}
	}
	return LookupError{
		query:       q
		suggestions: db.closest(lower, 5)
	}
}

// closest returns the n compounds whose name or an alias is nearest the
// query: first those that contain it ("butane" finds N-butane and
// Isobutane), then by edit distance.
fn (db Database) closest(lower string, n int) []string {
	mut scored := []Scored{}
	for _, e in db.entries {
		mut best := max_int
		for name in [e.name, ...e.aliases] {
			l := name.to_lower()
			// Containing the query outranks any edit distance; among those,
			// shorter names are closer.
			score := if l.contains(lower) {
				l.len - lower.len
			} else {
				1000 + strings.levenshtein_distance(l, lower)
			}
			if score < best {
				best = score
			}
		}
		scored << Scored{
			score: best
			label: '${e.name} (${e.cas})'
		}
	}
	scored.sort_with_compare(fn (a &Scored, b &Scored) int {
		if a.score != b.score {
			return a.score - b.score
		}
		return compare_strings(a.label, b.label)
	})
	return scored#[..n].map(it.label)
}

struct Scored {
	score int
	label string
}
