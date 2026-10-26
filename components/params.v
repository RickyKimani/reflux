module components

import math
import os
import strconv
import thermo

// ChemSep's model-parameter files, read straight from data/vendor/chemsep/:
// the kᵢⱼ libraries pr.ipd and srk.ipd, and mc-PR.lib, Mathias–Copeman
// constants fitted for Peng–Robinson. They're small, and nothing in them
// needs converting, so unlike the XML they aren't imported.
//
// Both are text files with a free-form header, then a section ([IPD] or
// [LIB]) holding a `Comment=` line, `#` comments and one entry per line.
// The shipped files have defects (NOTICE.md), so a data line that can't be
// read is skipped and kept as a Defect, never dropped silently. When a pair
// or a compound is listed more than once, the first listing is the one used
// and the rest are reported.

// vendor_dir holds the vendored ChemSep files.
pub const vendor_dir = os.join_path(@VMODROOT, 'data', 'vendor', 'chemsep')

// Defect is a data line that couldn't be read.
pub struct Defect {
pub:
	line int // 1-based
	text string
	why  string
}

// Kij is one binary interaction parameter as a file lists it.
pub struct Kij {
pub:
	cas1    string
	cas2    string
	k       f64
	comment string // e.g. 'Hydrogen/Nitrogen T=90-113K p210' (the DECHEMA page)
	line    int
}

// KijBank is one equation of state's kᵢⱼ file.
pub struct KijBank {
pub:
	file    string
	comment string // its Comment= line
	entries []Kij  // in file order
	defects []Defect
}

// load_kij loads ChemSep's kᵢⱼ for a cubic equation of state: DECHEMA's,
// for Peng–Robinson and Soave–Redlich–Kwong. ChemSep has none for vdW or RK.
pub fn load_kij(kind thermo.CubicKind) !KijBank {
	name := match kind {
		.pr { 'pr.ipd' }
		.srk { 'srk.ipd' }
		else { return error('ChemSep has no kᵢⱼ for ${kind}') }
	}
	return read_ipd(os.join_path(vendor_dir, name))!
}

// read_ipd reads a ChemSep .ipd file of kᵢⱼ.
pub fn read_ipd(path string) !KijBank {
	return parse_ipd(os.read_file(path)!, os.file_name(path))!
}

// parse_ipd reads the text of an .ipd file. Each data line is
// `CAS1 CAS2 k12 comment…`.
pub fn parse_ipd(text string, file string) !KijBank {
	comment, lines := lib_section(text, 'IPD') or { return error('${file}: ${err}') }
	mut entries := []Kij{}
	mut defects := []Defect{}
	for l in lines {
		f, rest := leading_fields(l.text, 3)
		if f.len < 3 {
			defects << Defect{l.n, l.text, 'fewer than three fields'}
			continue
		}
		check_cas(f[0]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		check_cas(f[1]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		if f[0] == f[1] {
			defects << Defect{l.n, l.text, 'a compound paired with itself'}
			continue
		}
		k := finite_number(f[2]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		entries << Kij{
			cas1:    f[0]
			cas2:    f[1]
			k:       k
			comment: rest
			line:    l.n
		}
	}
	return KijBank{
		file:    file
		comment: comment
		entries: entries
		defects: defects
	}
}

// for_pair returns every listing of the pair a, b, in file order, whichever
// way round each line names it.
pub fn (b KijBank) for_pair(a string, c string) []Kij {
	return b.entries.filter((it.cas1 == a && it.cas2 == c) || (it.cas1 == c && it.cas2 == a))
}

// McEntry is one compound's Mathias–Copeman constants as mc-PR.lib lists
// them. A `*` in the file, an unknown value, is none.
pub struct McEntry {
pub:
	cas     string
	c1      ?f64
	c2      ?f64
	c3      ?f64
	comment string // the name, and the fit's error on vapour pressure
	line    int
}

// McBank is mc-PR.lib.
pub struct McBank {
pub:
	file    string
	comment string
	entries []McEntry // in file order
	defects []Defect
}

// load_mathias_copeman_pr loads ChemSep's Mathias–Copeman constants for
// Peng–Robinson (mc-PR.lib). See M2 decision D10: the XML's own constants
// don't say which equation they were fitted for, and agree with none of
// these.
pub fn load_mathias_copeman_pr() !McBank {
	return read_lib(os.join_path(vendor_dir, 'mc-PR.lib'))!
}

// read_lib reads a ChemSep Mathias–Copeman .lib file.
pub fn read_lib(path string) !McBank {
	return parse_lib(os.read_file(path)!, os.file_name(path))!
}

// parse_lib reads the text of a Mathias–Copeman .lib file. Each data line is
// `CAS c1 c2 c3`, then either `# name` or tab-separated name and fit error.
pub fn parse_lib(text string, file string) !McBank {
	comment, lines := lib_section(text, 'LIB') or { return error('${file}: ${err}') }
	mut entries := []McEntry{}
	mut defects := []Defect{}
	for l in lines {
		f, rest := leading_fields(l.text, 4)
		if f.len < 4 {
			defects << Defect{l.n, l.text, 'fewer than four fields'}
			continue
		}
		check_cas(f[0]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		mut c := []?f64{}
		for s in f[1..4] {
			if s == '*' {
				c << none
				continue
			}
			c << finite_number(s) or {
				defects << Defect{l.n, l.text, err.msg()}
				break
			}
		}
		if c.len < 3 {
			continue
		}
		entries << McEntry{
			cas:     f[0]
			c1:      c[0]
			c2:      c[1]
			c3:      c[2]
			comment: rest.trim_left('#').trim_space().replace('\t', ', ')
			line:    l.n
		}
	}
	return McBank{
		file:    file
		comment: comment
		entries: entries
		defects: defects
	}
}

// find returns a compound's constants: its first listing, when it has more.
pub fn (b McBank) find(cas string) ?McEntry {
	for e in b.entries {
		if e.cas == cas {
			return e
		}
	}
	return none
}

// alpha returns the entry as a thermo α function. It fails if any of the
// three constants is unknown (`*` in the file).
pub fn (e McEntry) alpha() !thermo.Alpha {
	c1 := e.c1 or { return error('${e.cas}: Mathias–Copeman c1 is unknown') }
	c2 := e.c2 or { return error('${e.cas}: Mathias–Copeman c2 is unknown') }
	c3 := e.c3 or { return error('${e.cas}: Mathias–Copeman c3 is unknown') }
	return thermo.MathiasCopeman{
		c1: c1
		c2: c2
		c3: c3
	}
}

// What the two formats share.

struct Line {
	n    int
	text string
}

// lib_section returns a section's Comment= text and its data lines: everything
// after the `[name]` header that isn't blank, a `#` comment, or the Comment=
// line itself.
fn lib_section(text string, name string) !(string, []Line) {
	lines := text.split_into_lines()
	header := lines.index('[${name}]')
	if header < 0 {
		return error('no [${name}] section')
	}
	mut comment := ''
	mut data := []Line{}
	for i in header + 1 .. lines.len {
		l := lines[i]
		t := l.trim_space()
		if t == '' || t.starts_with('#') {
			continue
		}
		if t.starts_with('Comment=') {
			comment = t.all_after('Comment=').trim_space()
			continue
		}
		if t.starts_with('[') {
			break // the next section
		}
		data << Line{i + 1, l}
	}
	return comment, data
}

// leading_fields splits off the first n whitespace-separated fields (tabs
// and spaces alike: the files mix them) and returns the rest of the line,
// trimmed.
fn leading_fields(s string, n int) ([]string, string) {
	mut fields := []string{}
	mut i := 0
	for fields.len < n {
		for i < s.len && s[i] in [` `, `\t`] {
			i++
		}
		if i >= s.len {
			break
		}
		start := i
		for i < s.len && s[i] !in [` `, `\t`] {
			i++
		}
		fields << s[start..i]
	}
	return fields, s[i..].trim_space()
}

fn finite_number(s string) !f64 {
	v := strconv.atof64(s) or { return error('`${s}` is not a number') }
	if !math.is_finite(v) {
		return error('`${s}` is not finite')
	}
	return v
}
