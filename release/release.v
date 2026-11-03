// Package release is what has to be true before a version of reflux is
// tagged (M4, decision D15; roadmap A8): its CHANGELOG entry, dated, for
// v.mod's version; no tag of it yet; and every golden result changed since the
// last release named in that entry. tools/release.vsh gathers the facts from
// git and the files, asks check, and runs the full check before it prints the
// tag command. The logic is here, apart from git, so that it's tested.
module release

import strconv

// Entry is one release in the CHANGELOG: its version, its date (or
// 'unreleased'), and everything under its heading.
pub struct Entry {
pub:
	version string
	date    string
	text    string
}

// entries reads the CHANGELOG's releases, newest first: each heading of the
// form `## [x.y.z] — 2026-11-03`, and what follows it up to the next. A heading
// with nothing after its brackets is an entry still, with the date '', so that
// check can say what it lacks.
pub fn entries(changelog string) []Entry {
	mut out := []Entry{}
	mut current := ?Entry(none)
	mut lines := []string{}
	for line in changelog.split_into_lines() {
		if line.starts_with('## [') {
			if e := current {
				out << Entry{
					...e
					text: lines.join('\n')
				}
			}
			current = Entry{
				version: line.all_after('## [').all_before(']')
				date:    line.all_after(']').trim_space().trim_string_left('—').trim_space()
			}
			lines = []
		} else {
			lines << line // the preamble's too, dropped at the first heading
		}
	}
	if e := current {
		out << Entry{
			...e
			text: lines.join('\n')
		}
	}
	return out
}

// Version is x.y.z, as numbers to compare.
pub struct Version {
pub:
	major int
	minor int
	patch int
}

// parse_version reads x.y.z, each a whole number with no leading zeros, as
// semantic versioning has them. A tag's v isn't part of it.
pub fn parse_version(text string) !Version {
	parts := text.split('.')
	if parts.len != 3 {
		return error('${text} is not x.y.z')
	}
	mut n := []int{}
	for p in parts {
		// atoi alone would take a sign; it refuses '', and a number too large.
		if !p.bytes().all(it.is_digit()) || (p.len > 1 && p[0] == `0`) {
			return error('${text} is not x.y.z')
		}
		n << strconv.atoi(p) or { return error('${text} is not x.y.z') }
	}
	return Version{n[0], n[1], n[2]}
}

// before is whether a comes before b.
pub fn (a Version) before(b Version) bool {
	if a.major != b.major {
		return a.major < b.major
	}
	if a.minor != b.minor {
		return a.minor < b.minor
	}
	return a.patch < b.patch
}

// str is x.y.z.
pub fn (v Version) str() string {
	return '${v.major}.${v.minor}.${v.patch}'
}

// latest is the newest release among git tags, a v and x.y.z, or none.
pub fn latest(tags []string) ?Version {
	mut best := ?Version(none)
	for t in tags {
		if !t.starts_with('v') {
			continue
		}
		v := parse_version(t[1..]) or { continue }
		if b := best {
			if b.before(v) {
				best = v
			}
		} else {
			best = v
		}
	}
	return best
}

// Candidate is what a release would be made from.
pub struct Candidate {
pub:
	version   string   // v.mod's
	changelog string   // CHANGELOG.md
	tags      []string // git's tags
	clean     bool     // whether the working tree has nothing to commit
	branch    string   // the branch checked out, or '' for none
	// golden are the golden flowsheets' reports changed since the latest
	// release, by file name.
	golden []string
}

// check is what stands in the way of releasing a candidate, one reason a
// line: nothing, if it may be tagged.
pub fn check(c Candidate) []string {
	mut out := []string{}
	if !c.clean {
		out << 'the working tree has changes: commit them, or set them aside'
	}
	if c.branch != 'main' {
		out << 'a release is tagged on main, and this is ${if c.branch == '' {
			'no branch'
		} else {
			c.branch
		}}'
	}
	version := parse_version(c.version) or {
		out << "v.mod's version, ${c.version}, is not x.y.z"
		return out
	}
	all := entries(c.changelog)
	newest := all[0] or {
		out << 'the CHANGELOG has no release entries'
		return out
	}
	if newest.version != c.version {
		out << "the CHANGELOG's newest entry is ${newest.version}, not v.mod's ${version}"
	}
	if !is_date(newest.date) {
		said := if newest.date == '' { '' } else { ' (it says ${newest.date})' }
		out << "the CHANGELOG's entry for ${newest.version} has no date${said}: give it today's"
	}
	if 'v${version}' in c.tags {
		out << 'v${version} is tagged already'
	} else if last := latest(c.tags) {
		if !last.before(version) {
			out << 'v${version} is not after the latest release, v${last}'
		}
	}
	for f in c.golden {
		if !newest.text.contains(f) {
			out << "the golden result ${f} changed since the latest release, and the CHANGELOG entry for ${newest.version} doesn't say why"
		}
	}
	return out
}

// is_date is whether text is a date, yyyy-mm-dd.
fn is_date(text string) bool {
	if text.len != 10 || text[4] != `-` || text[7] != `-` {
		return false
	}
	for i, c in text {
		if i != 4 && i != 7 && !c.is_digit() {
			return false
		}
	}
	month := text[5..7].int()
	day := text[8..10].int()
	return month >= 1 && month <= 12 && day >= 1 && day <= 31
}
