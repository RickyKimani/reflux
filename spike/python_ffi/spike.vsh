#!/usr/bin/env -S v

// spike.vsh builds the ffi module as a shared library and runs the Python
// scenarios against it.
//
//     v run spike/python_ffi/spike.vsh            # gcc
//     v run spike/python_ffi/spike.vsh -cc tcc    # to see what tcc does
//
// Needs uv (for Python and NumPy) on PATH.

fn run(cmd []string) int {
	println('\n> ${cmd.join(' ')}')
	return system_args(cmd)
}

args := arguments()
i := args.index('-cc')
cc := if i < 0 { 'gcc' } else { args[i + 1] or { 'gcc' } }

here := dir(@FILE)
out_dir := join_path(temp_dir(), 'reflux_ffi_spike')
mkdir_all(out_dir)!
// V adds the extension itself: .dll on Windows, .so elsewhere.
base := join_path(out_dir, 'reflux_ffi_${cc}')
lib := base + $if windows { '.dll' } $else { '.so' }

mut build := ['v', '-cc', cc, '-W', '-shared', '-enable-globals']
$if windows {
	// Without this the DLL needs libgcc_s_seh-1.dll, which Python won't find:
	// since 3.8 it doesn't search PATH for a DLL's dependencies.
	if cc == 'gcc' {
		build << ['-cflags', '-static-libgcc']
	}
}
build << ['-o', base, join_path(here, 'ffi')]
if run(build) != 0 {
	eprintln('build failed')
	exit(1)
}
exit(run(['uv', 'run', join_path(here, 'run.py'), lib]))
