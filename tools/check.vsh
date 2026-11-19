#!/usr/bin/env -S v

// check.vsh runs the same checks as CI, so "it passes locally" means the
// same thing as "CI is green".
//
//     v run tools/check.vsh            # tcc
//     v run tools/check.vsh -cc gcc    # gcc, plus an optimised -prod build
//
// A .vsh file is a V script: top-level statements run in order, and the
// `os` module is imported implicitly, so `system_args`, `chdir` and friends
// need no `os.` prefix. It's the V answer to a Makefile or a shell script,
// and it runs the same on Windows, Linux and any shell, nushell included.

// run executes one command and stops the script if it fails. Commands are
// argument arrays, not strings: no shell parses them, so paths with spaces
// need no quoting and nothing can be injected.
fn run(cmd []string) {
	println('\n> ${cmd.join(' ')}')
	code := system_args(cmd)
	if code != 0 {
		eprintln('\nFAILED (exit ${code}): ${cmd.join(' ')}')
		exit(code)
	}
}

// expect runs a command that must exit with `want`, such as a failure.
fn expect(cmd []string, want int) {
	println('\n> ${cmd.join(' ')}    (should exit ${want})')
	code := system_args(cmd)
	if code != want {
		eprintln('\nFAILED (exit ${code}, not ${want}): ${cmd.join(' ')}')
		exit(1)
	}
}

// smoke runs a built reflux as someone would: the demo, the reference
// flowsheet with both its reports, and a flowsheet that isn't there, which must
// fail with exit status 1 (criterion 8); and M9's kinetics check and explain
// on Xu and Froment's kinetics.
fn smoke(exe string, reports string) {
	run([exe, 'demo'])
	run([exe, 'solve', 'testdata/flowsheets/reference.toml', '--json', '${reports}.json', '--csv',
		'${reports}.csv'])
	expect([exe, 'solve', 'testdata/flowsheets/no such file.toml'], 1)
	run([exe, 'kinetics', 'check', 'testdata/kinetics/xu_froment.toml'])
	run([exe, 'kinetics', 'explain', 'testdata/kinetics/xu_froment.toml', '900 K', '15 bar', 'CH4=0.2',
		'H2O=0.45', 'H2=0.2', 'CO=0.1', 'CO2=0.05'])
	rm('${reports}.json') or {}
	rm('${reports}.csv') or {}
}

// vcc is `v -cc <cc>` followed by args. With gcc, V 0.5.2's module cache
// makes every build about three times slower: 40 s against 12 s for a small
// test on Windows, most of it spent scanning cached objects' dependencies
// and replaying a build that V's cached path gives up on. So gcc builds get
// -nocache; tcc is as fast either way. V has no `+` for arrays: `<<` appends
// one element, or every element of another array, like Go's append(a, b...).
fn vcc(cc string, args []string) []string {
	mut cmd := ['v', '-cc', cc]
	if cc != 'tcc' {
		cmd << '-nocache'
	}
	cmd << args
	return cmd
}

args := arguments()
i := args.index('-cc') // -1 when absent; index does not return an Option
cc := if i < 0 {
	'tcc'
} else {
	args[i + 1] or {
		eprintln('-cc needs a value, e.g. -cc gcc')
		exit(2)
	}
}

// Run from the project root, wherever the script was started from.
chdir(dir(dir(@FILE)))!

// CI builds the V commit pinned in .v-version. @VCURRENTHASH is the commit
// of the V that compiled this script, filled in at compile time. Warn
// rather than fail: a local V that differs from CI's is fine while
// experimenting, but it explains "passes here, fails in CI".
pinned := read_file('.v-version')!.trim_space()
if !pinned.starts_with(@VCURRENTHASH) {
	eprintln('warning: local V is ${@VCURRENTHASH} but CI pins ${pinned[..7]} (.v-version)')
}

// `$if` is evaluated at compile time, like a Go build tag but inline.
ext := $if windows { '.exe' } $else { '' }
// A fresh name per run, removed at the end: on Windows the previous run's
// binary can still be locked (antivirus scanning it, typically) when the
// next build tries to overwrite it.
tag := '${cc}_${getpid()}'
exe := join_path(temp_dir(), 'reflux_check_${tag}${ext}')
exe_prod := join_path(temp_dir(), 'reflux_check_${tag}_prod${ext}')

// `v test` compiles each _test.v file into a program of its own, and its
// runner allows one job per 8 GiB of RAM: one at a time on a 14 GiB
// machine. A build peaks near 1 GB, so run 4, the runner's own cap. (Not
// one per CPU: on Windows, `import runtime` defeats V 0.5.2's default
// compiler, through tcc's psapi.h, and only its fallback builds the
// script.) VJOBS set by hand still wins.
if getenv('VJOBS') == '' {
	setenv('VJOBS', '4', true)
}
run(['v', 'fmt', '-verify', '.'])
run(['v', 'vet', '.'])
run(vcc(cc, ['-W', 'test', '.']))
run(vcc(cc, ['-W', '-o', exe, '.']))
smoke(exe, join_path(temp_dir(), 'reflux_check_${tag}'))
if cc != 'tcc' {
	// tcc does no optimisation, so -prod only means something with gcc or clang.
	run(vcc(cc, ['-prod', '-o', exe_prod, '.']))
	smoke(exe_prod, join_path(temp_dir(), 'reflux_check_${tag}_prod'))
}
rm(exe) or {}
rm(exe_prod) or {}
println('\nAll checks passed with ${cc}.')
