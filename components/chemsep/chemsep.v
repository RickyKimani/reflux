// Module chemsep imports ChemSep's pure-component database, the XML files
// vendored in data/vendor/chemsep/, into components.Record values.
//
// It is strict on purpose. Every element must be one it knows (fields.v),
// carrying the unit it expects, and every number must parse. A future
// ChemSep release that adds or changes something then fails the import
// loudly instead of being misread. The known defect of ChemSep 8.50 is
// repaired, and every repair is reported, as are duplicates and anything
// else a person should look at.
module chemsep

import encoding.xml
import math
import os
import strconv
import components
import tdep

// Source describes one imported file.
pub struct Source {
pub:
	path    string // relative to the root given to import_files, with forward slashes
	library string // the root element's `library` attribute
	version string // e.g. 'v8.49'
	records int    // compound records in the file, duplicates included
}

// Repair is a defect in the source that the importer fixed.
pub struct Repair {
pub:
	cas    string
	name   string
	key    string // e.g. 'ids.smiles'
	before string
	after  string
}

// Duplicate is a compound whose CAS number appeared again. The first
// occurrence in reading order is kept (M2 decision D8).
pub struct Duplicate {
pub:
	cas       string
	name      string
	kept      string   // where, e.g. 'chemsep2.xml record 179'
	dropped   string   // where the dropped copy was
	differing []string // the elements that differ between the two copies
}

pub enum FindingKind {
	bad_range           // a correlation range that is inverted, or ends at or below 0 K
	unsupported_form    // an equation form tdep can't evaluate
	unused_coefficients // non-zero coefficients the form doesn't use
	no_coefficients     // a correlation with no coefficients at all
	partial_set         // some of the Mathias–Copeman constants missing
	incomplete          // lacks a constant that a thermo.Component needs
}

// Finding is something the import report shows, though nothing was changed.
pub struct Finding {
pub:
	kind FindingKind
	cas  string
	name string
	what string
}

// Import is the result of importing one or more ChemSep files. Repairs and
// findings are for the records kept, not for dropped duplicates.
pub struct Import {
pub mut:
	sources    []Source
	records    []components.Record // one per CAS number, in reading order
	repairs    []Repair
	duplicates []Duplicate
	findings   []Finding
mut:
	kept map[string]Kept // CAS number → where its kept record came from
}

struct Kept {
	where    string
	elements map[string]string // element name → canonical form
}

// Parsed is one compound record, with what was found while reading it.
struct Parsed {
	record   components.Record
	repairs  []Repair
	findings []Finding
}

// import_files imports ChemSep XML files in order. The paths are relative
// to root, and the records' provenance names them that way, e.g.
// 'data/vendor/chemsep/chemsep1.xml'. A compound whose CAS number was
// already seen, in an earlier file or earlier in the same one, is dropped
// and listed as a duplicate.
pub fn import_files(root string, paths []string) !Import {
	mut imp := Import{}
	for path in paths {
		text := os.read_file(os.join_path(root, path))!
		imp.add_file(text, path.replace('\\', '/'))!
	}
	return imp
}

// add_file imports the text of one file. `path` is how the records'
// provenance names it.
fn (mut imp Import) add_file(text string, path string) ! {
	file := os.file_name(path)
	// ChemSep 8.50's defect is a stray `<` ending some SMILES values, as in
	// value="CC(C)(C)C1CCCCC1<". `<"` can't occur in well-formed XML, so it's
	// escaped here, and the `<` is removed (and reported) when the value is
	// read. Any other misplaced `<` still fails the parse.
	doc := xml.XMLDocument.from_string(text.replace('<"', '&lt;"')) or {
		return error('${file}: ${err}')
	}
	root := doc.root
	if root.name != 'compounds' {
		return error('${file}: the root element is <${root.name}>, not <compounds>')
	}
	library := root.attributes['library'] or { '' }
	version := library_version(library) or { return error('${file}: ${err}') }
	source := 'ChemSep ${version} (${path})'
	mut n := 0
	for node in elements(root) or { return error('${file}: ${err}') } {
		n++
		where := '${file} record ${n}'
		if node.name != 'compound' {
			return error('${where}: <${node.name}>, not <compound>')
		}
		p := compound(node, source) or { return error('${where}: ${err}') }
		cas := p.record.ids.cas
		elems := canonical_elements(node)
		if kept := imp.kept[cas] {
			imp.duplicates << Duplicate{
				cas:       cas
				name:      p.record.name
				kept:      kept.where
				dropped:   where
				differing: differing(kept.elements, elems)
			}
			continue
		}
		imp.kept[cas] = Kept{
			where:    where
			elements: elems
		}
		imp.records << p.record
		imp.repairs << p.repairs
		imp.findings << p.findings
	}
	imp.sources << Source{
		path:    path
		library: library
		version: version
		records: n
	}
}

// library_version finds the version in a `library` attribute such as
// 'ChemSep v8.49 pure component data 1 - Copyright …'.
fn library_version(library string) !string {
	words := library.fields()
	if words.len >= 2 && words[0] == 'ChemSep' && words[1].starts_with('v') {
		return words[1]
	}
	return error('no ChemSep version in the library attribute `${library}`')
}

// compound converts one <compound> element.
fn compound(node xml.XMLNode, source string) !Parsed {
	mut rec := components.Record{
		source: source
	}
	mut repairs := []Repair{}
	mut notes := []Finding{} // cas and name are filled in at the end
	mut seen := map[string]bool{}
	for el in elements(node)! {
		f := fields[el.name] or { return error('unknown element <${el.name}>') }
		if seen[el.name] {
			return error('<${el.name}> appears twice')
		}
		seen[el.name] = true
		units := el.attributes['units'] or { '' }
		if units != f.units {
			return error('<${el.name}> has units `${units}`, expected `${f.units}`')
		}
		match f.kind {
			.name {
				rec.name = scalar(el)!
			}
			.text_id {
				mut v := scalar_raw(el)!
				if f.key == 'smiles' && v.ends_with('<') {
					repairs << Repair{
						key:    'ids.smiles'
						before: v
						after:  v#[..-1]
					}
					v = v#[..-1]
				}
				check_text(el.name, v)!
				match f.key {
					'cas' { rec.ids.cas = v }
					'smiles' { rec.ids.smiles = v }
					'inchikey' { rec.ids.inchikey = v }
					'structure' { rec.ids.structure = v }
					else { return error('no text id `${f.key}`') }
				}
			}
			.int_id {
				v := integer(scalar(el)!)!
				match f.key {
					'chemsep_index' { rec.ids.chemsep_index = v }
					'family' { rec.ids.family = v }
					'pubchem' { rec.ids.pubchem = v }
					else { return error('no integer id `${f.key}`') }
				}
			}
			.constant {
				rec.constants.set(f.key, number(scalar(el)!, f.per_kmol)!)!
			}
			.int_constant {
				rec.constants.charge = integer(scalar(el)!)!
			}
			.parameter {
				rec.parameters.set(f.key, number(scalar(el)!, f.per_kmol)!)!
			}
			.correlation {
				d := correlation(el, units)!
				rec.correlations.set(f.key, d)!
				notes << correlation_findings(f.key, d)
			}
			.groups {
				counts, sums := groups(el, f.key)!
				rec.groups.set(f.key, counts)!
				repairs << sums
			}
		}
	}
	if rec.ids.cas == '' {
		return error('no CAS number')
	}
	components.check_cas(rec.ids.cas)!
	if rec.name == '' {
		return error('${rec.ids.cas} has no name')
	}
	notes << record_findings(rec)
	return Parsed{
		record:   rec
		repairs:  repairs.map(Repair{
			...it
			cas:  rec.ids.cas
			name: rec.name
		})
		findings: notes.map(Finding{
			...it
			cas:  rec.ids.cas
			name: rec.name
		})
	}
}

// correlation reads a correlation element: <eqno>, the coefficients A, B, …
// in order with none skipped, and <Tmin> and <Tmax> in K.
fn correlation(el xml.XMLNode, units string) !components.CorrelationData {
	check_attributes(el, ['name', 'units'])!
	mut d := components.CorrelationData{
		units: units
	}
	mut have := map[string]bool{}
	for c in elements(el)! {
		if have[c.name] {
			return error('<${el.name}> has <${c.name}> twice')
		}
		have[c.name] = true
		match c.name {
			'eqno' {
				check_attributes(c, ['value'])!
				d.eq = integer(c.attributes['value'] or { '' })!
			}
			'A', 'B', 'C', 'D', 'E' {
				check_attributes(c, ['value'])!
				expected := 'ABCDE'[d.coeffs.len..d.coeffs.len + 1]
				if c.name != expected {
					return error('<${el.name}> has <${c.name}> where <${expected}> should be')
				}
				d.coeffs << number(c.attributes['value'] or { '' }, false)!
			}
			'Tmin', 'Tmax' {
				check_attributes(c, ['units', 'value'])!
				if c.attributes['units'] or { '' } != 'K' {
					return error('<${el.name}> <${c.name}> is not in K')
				}
				t := number(c.attributes['value'] or { '' }, false)!
				if c.name == 'Tmin' {
					d.tmin = t
				} else {
					d.tmax = t
				}
			}
			else {
				return error('<${el.name}> has an unknown child <${c.name}>')
			}
		}
	}
	for need in ['eqno', 'Tmin', 'Tmax'] {
		if !have[need] {
			return error('<${el.name}> has no <${need}>')
		}
	}
	// An inverted or empty range is a defect in the data, not in the
	// file: it's kept as it is and reported (see correlation_findings).
	return d
}

// groups reads a group-count element's <group id=… value=…> children.
//
// Six ChemSep 8.50 compounds list an ASOG group twice, for instance
// diglyme's ether group 19 as 2 and then 1. Each time the sum is the
// right count (diglyme has three ether oxygens), so repeats are added up,
// and each sum is returned as a repair to report. One compound, diphenyl
// disulfide, has a UNIFAC-LLE list that is just "group 10, count 0": a
// zero count contributes nothing, so it's dropped, and reported, keeping
// the component format's rule that counts are positive.
fn groups(el xml.XMLNode, key string) !(map[int]int, []Repair) {
	check_attributes(el, ['name'])!
	mut counts := map[int]int{}
	mut listed := map[int][]int{}
	for c in elements(el)! {
		if c.name != 'group' {
			return error('<${el.name}> has an unknown child <${c.name}>')
		}
		check_attributes(c, ['id', 'value'])!
		id := integer(c.attributes['id'] or { '' })!
		n := integer(c.attributes['value'] or { '' })!
		if id <= 0 || n < 0 {
			return error('<${el.name}> has group ${id} with count ${n}')
		}
		counts[id] += n
		listed[id] << n
	}
	mut repairs := []Repair{}
	for id, ns in listed {
		if counts[id] == 0 {
			counts.delete(id)
			repairs << Repair{
				key:    'groups.${key}'
				before: 'group ${id}: ' + ns.map(it.str()).join(' + ')
				after:  'dropped'
			}
		} else if ns.len > 1 {
			repairs << Repair{
				key:    'groups.${key}'
				before: 'group ${id}: ' + ns.map(it.str()).join(' + ')
				after:  'group ${id}: ${counts[id]}'
			}
		}
	}
	return counts, repairs
}

// correlation_findings notes what the report should show about a
// correlation: a range that can't be right, a form reflux can't evaluate,
// unused non-zero coefficients, or no coefficients at all.
fn correlation_findings(key string, d components.CorrelationData) []Finding {
	mut found := []Finding{}
	if d.tmin > d.tmax || d.tmin < 0.0 || d.tmax <= 0.0 {
		found << Finding{
			kind: .bad_range
			what: '${key}: range [${d.tmin}, ${d.tmax}] K'
		}
	}
	if !tdep.is_supported(d.eq) {
		found << Finding{
			kind: .unsupported_form
			what: '${key}: equation ${d.eq}'
		}
	} else if d.coeffs.len == 0 {
		found << Finding{
			kind: .no_coefficients
			what: '${key}: equation ${d.eq} with no coefficients'
		}
	} else if tdep.unused_nonzero(d.eq, d.coeffs) {
		found << Finding{
			kind: .unused_coefficients
			what: '${key}: equation ${d.eq} with coefficients ${d.coeffs}'
		}
	}
	return found
}

// record_findings notes gaps in a whole record.
fn record_findings(rec components.Record) []Finding {
	mut found := []Finding{}
	p := rec.parameters
	mc := [p.mathias_copeman_c1 != none, p.mathias_copeman_c2 != none, p.mathias_copeman_c3 != none]
	if mc.any(it) && !mc.all(it) {
		mut have := []string{}
		for i, name in ['C1', 'C2', 'C3'] {
			if mc[i] {
				have << name
			}
		}
		found << Finding{
			kind: .partial_set
			what: 'Mathias–Copeman: only ${have.join(', ')}'
		}
	}
	c := rec.constants
	missing := ['tc', 'pc', 'omega', 'mw'].filter(match it {
		'tc' { c.tc == none }
		'pc' { c.pc == none }
		'omega' { c.omega == none }
		else { c.mw == none }
	})
	if missing.len > 0 {
		found << Finding{
			kind: .incomplete
			what: 'no ${missing.join(', ')}'
		}
	}
	return found
}

// elements returns a node's child elements, skipping formatting whitespace
// and comments. Text or CDATA where elements belong is an error.
fn elements(node xml.XMLNode) ![]xml.XMLNode {
	mut out := []xml.XMLNode{}
	for c in node.children {
		match c {
			xml.XMLNode {
				out << c
			}
			xml.XMLComment {}
			string {
				if c.trim_space() != '' {
					return error('<${node.name}> contains text `${c}`')
				}
			}
			xml.XMLCData {
				return error('<${node.name}> contains CDATA')
			}
		}
	}
	return out
}

// scalar returns a value element's `value`, checked: no children, only the
// expected attributes, and no stray markup.
fn scalar(el xml.XMLNode) !string {
	v := scalar_raw(el)!
	check_text(el.name, v)!
	return v
}

// scalar_raw is scalar without the markup check, for values that may need
// repairing first.
fn scalar_raw(el xml.XMLNode) !string {
	if elements(el)!.len > 0 {
		return error('<${el.name}> has children, expected a value')
	}
	check_attributes(el, ['name', 'units', 'value'])!
	return el.attributes['value'] or { return error('<${el.name}> has no value') }
}

fn check_text(name string, v string) ! {
	if v.contains_any('<>&') {
		return error('<${name}> value `${v}` contains markup')
	}
}

// check_attributes refuses attributes outside `allowed`.
fn check_attributes(el xml.XMLNode, allowed []string) ! {
	for k, _ in el.attributes {
		if k !in allowed {
			return error('<${el.name}> has an unknown attribute `${k}`')
		}
	}
}

// number parses a decimal number, which must be finite. per_kmol values
// are divided by 1000 exactly (see shift_decimal).
fn number(s string, per_kmol bool) !f64 {
	text := if per_kmol { shift_decimal(s, -3)! } else { s }
	v := strconv.atof64(text) or { return error('`${s}` is not a number') }
	if !math.is_finite(v) {
		return error('`${s}` is not a finite number')
	}
	return v
}

// shift_decimal moves a decimal number's point by `places`: '58.1222' with
// −3 becomes '58.1222e-3'. Parsing that gives the double nearest to
// 0.0581222 itself. Dividing the parsed 58.1222 by 1000 instead rounds
// twice, and can land a unit in the last place away.
fn shift_decimal(s string, places int) !string {
	i := s.index_any('eE')
	if i < 0 {
		return '${s}e${places}'
	}
	e := strconv.atoi(s[i + 1..].trim_left('+')) or { return error('`${s}` is not a number') }
	return '${s[..i]}e${e + places}'
}

fn integer(s string) !int {
	return strconv.atoi(s) or { return error('`${s}` is not an integer') }
}

// canonical_elements maps each child element's name to a canonical text
// form, for telling which elements two copies of a compound differ in.
fn canonical_elements(node xml.XMLNode) map[string]string {
	mut out := map[string]string{}
	for el in elements(node) or { [] } {
		out[el.name] = canonical(el)
	}
	return out
}

fn canonical(n xml.XMLNode) string {
	mut keys := n.attributes.keys()
	keys.sort()
	mut s := '<${n.name}'
	for k in keys {
		s += ' ${k}="${n.attributes[k]}"'
	}
	s += '>'
	for c in elements(n) or { [] } {
		s += canonical(c)
	}
	return s + '</${n.name}>'
}

// differing lists the element names whose canonical forms differ between
// a and b, including those only one has, in a's order and then b's.
fn differing(a map[string]string, b map[string]string) []string {
	mut out := []string{}
	for k, v in a {
		if (b[k] or { '' }) != v {
			out << k
		}
	}
	for k, _ in b {
		if k !in a {
			out << k
		}
	}
	return out
}
