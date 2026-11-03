#!/usr/bin/env -S v

// release.vsh is the release gate (M4, decision D15; roadmap A8): whether
// this commit may be tagged as v.mod's version, and if so, how. CHANGELOG.md's
// "Making a release" is the procedure:
//
//     v run tools/release.vsh
//
// It refuses, with every reason, unless release.check passes: on main, with a
// clean tree; v.mod's version the CHANGELOG's newest entry, and dated; not
// tagged yet, and after the latest release; and every golden result changed
// since that release named in the entry. Then it runs the full check, the
// golden flowsheets among it, and prints the tag command. It tags nothing.
import release
import v.vmod

// lines is a git command's output, a line each, or the script stops. Commands
// are argument arrays, as in check.vsh: no shell parses them.
fn lines(cmd []string) []string {
	r := exec(cmd)
	if r.exit_code != 0 {
		eprintln('${cmd.join(' ')} failed: ${r.output}')
		exit(1)
	}
	return r.output.split_into_lines().map(it.trim_space()).filter(it != '')
}

chdir(dir(dir(@FILE)))!
version := vmod.decode(read_file('v.mod')!)!.version
tags := lines(['git', 'tag', '--list'])
mut golden := []string{}
if last := release.latest(tags) {
	golden = lines(['git', 'diff', '--name-only', 'v${last}', 'HEAD', '--', 'testdata/flowsheets/golden']).map(file_name(it))
}
problems := release.check(release.Candidate{
	version:   version
	changelog: read_file('CHANGELOG.md')!
	tags:      tags
	clean:     lines(['git', 'status', '--porcelain']).len == 0
	branch:    lines(['git', 'branch', '--show-current'])[0] or { '' }
	golden:    golden
})
if problems.len > 0 {
	eprintln("v${version} can't be released yet:")
	for p in problems {
		eprintln('  - ${p}')
	}
	exit(1)
}
println('> v run tools/check.vsh')
if system_args(['v', 'run', 'tools/check.vsh']) != 0 {
	eprintln("\nthe check failed, so v${version} can't be released")
	exit(1)
}
println('\nv${version} is ready. Tag it, and push the tag:\n')
println('    git tag -a v${version} -m "reflux ${version}"')
println('    git push origin v${version}')
