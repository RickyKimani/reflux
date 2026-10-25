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

run(['v', 'fmt', '-verify', '.'])
run(['v', 'vet', '.'])
run(['v', '-cc', cc, '-W', 'test', '.'])
run(['v', '-cc', cc, '-W', '-o', exe, '.'])
run([exe])
if cc != 'tcc' {
	// tcc does no optimisation, so -prod only means something with gcc or clang.
	run(['v', '-cc', cc, '-prod', '-o', exe_prod, '.'])
	run([exe_prod])
}
rm(exe) or {}
rm(exe_prod) or {}
println('\nAll checks passed with ${cc}.')
