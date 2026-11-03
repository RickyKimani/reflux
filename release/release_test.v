// The release gate (D15): what check refuses, and that the CHANGELOG and
// v.mod agree, as the next release will need them to.
import os
import release
import v.vmod

const changelog = '# Changelog

Notes.

## [0.2.0] — 2027-01-15

### Changed
- `reference.json` moved: the flash is more exact.

## [0.1.0] — 2026-11-03

The first release.
'

fn ok_candidate() release.Candidate {
	return release.Candidate{
		version:   '0.2.0'
		changelog: changelog
		tags:      ['v0.1.0']
		clean:     true
		branch:    'main'
		golden:    ['reference.json']
	}
}

fn test_the_changelog_is_read_newest_first() {
	es := release.entries(changelog)
	assert es.len == 2
	assert es[0].version == '0.2.0' && es[0].date == '2027-01-15'
	assert es[0].text.contains('`reference.json` moved') && !es[0].text.contains('first release')
	assert es[1].version == '0.1.0' && es[1].text.contains('The first release.')
	// Each entry's text is its own: not the preamble's, nor a newer entry's.
	assert !es[0].text.contains('Notes.') && !es[1].text.contains('reference.json')
	assert release.entries('# Changelog\n\nNothing yet.\n').len == 0
	// A heading with no date is an entry still, and spaces after one aren't in it.
	undated := release.entries('## [0.3.0]\n\nSoon.\n\n## [0.2.0] — 2027-01-15  \n')
	assert undated.len == 2 && undated[0].version == '0.3.0' && undated[0].date == ''
	assert undated[0].text.contains('Soon.') && undated[1].date == '2027-01-15'
}

fn test_a_ready_release_passes() {
	assert release.check(ok_candidate()) == []
	// The first release has no tags, and so no golden results to account for.
	first := release.Candidate{
		version:   '0.1.0'
		changelog: changelog.all_after('### Changed\n- `reference.json` moved: the flash is more exact.\n\n')
		clean:     true
		branch:    'main'
	}
	assert release.check(first) == [], '${release.check(first)}'
}

fn test_each_reason_to_wait_is_given() {
	c := ok_candidate()
	// Each wrong in turn, and only it, gives its own reason.
	expect(release.Candidate{ ...c, clean: false }, 'the working tree has changes')
	expect(release.Candidate{ ...c, branch: 'm4-flowsheet' }, 'a release is tagged on main, and this is m4-flowsheet')
	expect(release.Candidate{ ...c, branch: '' }, 'a release is tagged on main, and this is no branch')
	expect(release.Candidate{ ...c, version: '0.3.0' }, "the CHANGELOG's newest entry is 0.2.0, not v.mod's 0.3.0")
	expect(release.Candidate{ ...c, changelog: changelog.replace('2027-01-15', 'unreleased') },
		"the CHANGELOG's entry for 0.2.0 has no date (it says unreleased): give it today's")
	expect(release.Candidate{ ...c, changelog: changelog.replace(' — 2027-01-15', '') },
		"the CHANGELOG's entry for 0.2.0 has no date: give it today's")
	expect(release.Candidate{ ...c, tags: ['v0.1.0', 'v0.2.0'] }, 'v0.2.0 is tagged already')
	expect(release.Candidate{ ...c, tags: ['v0.1.0', 'v1.0.0'] }, 'v0.2.0 is not after the latest release, v1.0.0')
	expect(release.Candidate{ ...c, golden: ['reference.json', 'saturated_loop.json'] },
		'the golden result saturated_loop.json changed since the latest release')
	expect(release.Candidate{ ...c, version: '0.2' }, "v.mod's version, 0.2, is not x.y.z")
	expect(release.Candidate{ ...c, version: 'v0.2.0' }, "v.mod's version, v0.2.0, is not x.y.z")
	expect(release.Candidate{ ...c, changelog: '# Changelog\n' }, 'the CHANGELOG has no release entries')
}

// expect is that check gives exactly one reason, starting so.
fn expect(c release.Candidate, reason string) {
	got := release.check(c)
	assert got.len == 1 && got[0].starts_with(reason), '${got}'
}

fn test_versions_are_semantic() {
	for text in ['0.1.0', '0.0.0', '10.20.30'] {
		release.parse_version(text) or { assert false, text }
	}
	for text in ['0.1', '0.1.0.0', '01.1.0', '0.a.0', '', 'v', '0..1', 'v0.1.0', '0.-1.0', '0.+1.0',
		' 0.1.0', '0.1.0 '] {
		if v := release.parse_version(text) {
			assert false, '${text} read as ${v}'
		}
	}
	a, b := release.parse_version('0.9.9')!, release.parse_version('0.10.0')!
	assert a.before(b) && !b.before(a) && !a.before(a)
	assert release.parse_version('1.0.0')!.before(release.parse_version('1.0.1')!)
	assert release.parse_version('1.2.0')!.before(release.parse_version('2.0.0')!)
	// A release's tag is v and x.y.z: not x.y.z alone, nor vv, nor another letter.
	assert release.latest(['v0.9.9', 'v0.10.0', 'notes', 'v2', 'v0.2.0', '3.0.0', 'vv4.0.0', 'r9.0.0']) or {
		panic('none')
	} == b
	if v := release.latest(['notes']) {
		assert false, '${v}'
	}
}

fn test_dates_are_checked() {
	for date, ok in {
		'2026-11-03':  true
		'2026-13-03':  false
		'2026-11-3':   false
		'unreleased':  false
		'2026/11/03':  false
		'2026-00-10':  false
		'2026-11-32':  false
		'2026-11-00':  false
		'2o26-11-03':  false
		'2026-11-033': false
		'2026/11-03':  false
		'2026-11/03':  false
	} {
		c := release.Candidate{
			...ok_candidate()
			changelog: changelog.replace('2027-01-15', date)
		}
		assert (release.check(c) == []) == ok, date
	}
}

fn test_the_changelog_is_for_this_version() {
	// The committed CHANGELOG's newest entry is v.mod's version, so that the
	// next release only has to date it.
	root := @VMODROOT
	text := os.read_file(os.join_path(root, 'CHANGELOG.md'))!
	version := vmod.decode(os.read_file(os.join_path(root, 'v.mod'))!)!.version
	newest := release.entries(text)[0] or { panic('no entries') }
	assert newest.version == version, '${newest.version} against ${version}'
}
