#!/usr/bin/env -S v

// python.vsh builds reflux's C API into the Python package, and runs the
// package's tests (M10a, decisions D1, D3, D10 and D11).
//
//     v run tools/python.vsh build             # python/reflux/_lib/, with gcc -prod
//     v run tools/python.vsh build -debug      # the same without -prod: quicker
//     v run tools/python.vsh test              # build, then pytest with uv
//     v run tools/python.vsh test --python <path to python>
//
// The library is built with gcc only: tcc's Windows DLLs export nothing, and
// its Linux .so files ask for an executable stack that glibc refuses (M1's
// spike). It exports exactly the functions capi/reflux.h declares: on Windows
// by a .def file, and elsewhere by a linker version script, both written here
// from the header. -fvisibility=hidden wouldn't do on Linux: it can't hide the
// collector's symbols, which come from an object V built beforehand.

// The Python tests' dependencies, pinned (D11).
const pytest_version = '9.2.0'
const numpy_version = '2.5.4'

// run executes one command and stops the script if it fails. Commands are
// argument arrays, as in check.vsh: no shell parses them.
fn run(cmd []string) {
	println('\n> ${cmd.join(' ')}')
	code := system_args(cmd)
	if code != 0 {
		eprintln('\nFAILED (exit ${code}): ${cmd.join(' ')}')
		exit(code)
	}
}

// header_functions are the functions capi/reflux.h declares, in its order.
fn header_functions(header string) []string {
	mut names := []string{}
	for line in header.split_into_lines() {
		t := line.trim_space()
		if t.starts_with('/*') || t.starts_with('*') || t.starts_with('#') || !t.contains('(') {
			continue
		}
		name := t.all_before('(').all_after_last(' ').trim_left('*')
		if name.starts_with('reflux_') {
			names << name
		}
	}
	return names
}

// export_flags writes the export list for this OS into dir, and returns the
// linker flags that apply it.
fn export_flags(names []string, dir string) ![]string {
	mkdir_all(dir)!
	$if windows {
		path := join_path(dir, 'reflux.def')
		write_file(path, 'EXPORTS\n' + names.map('    ${it}\n').join(''))!
		return ['-ldflags', path]
	} $else {
		path := join_path(dir, 'reflux.map')
		write_file(path, '{\n  global:\n' + names.map('    ${it};\n').join('') + '  local: *;\n};\n')!
		return ['-ldflags', '-Wl,--version-script=${path}']
	}
}

// build makes the library: python/reflux/_lib/reflux.dll, or reflux.so.
fn build(prod bool) ! {
	names := header_functions(read_file('capi/reflux.h')!)
	if names.len == 0 {
		return error('capi/reflux.h declares no functions')
	}
	// CAPI_LIBRARY tells capi/gc_threads.h it's in the library, whose first call
	// starts the collector in the caller's thread.
	mut cmd := ['v', '-cc', 'gcc', '-nocache', '-W', '-shared', '-cflags', '-DCAPI_LIBRARY']
	if prod {
		cmd << '-prod'
	}
	$if windows {
		// Without it a DLL that calls one of libgcc's helpers needs
		// libgcc_s_seh-1.dll, which Python won't find: since 3.8 it doesn't
		// search PATH for a DLL's dependencies. (gcc 16.2's builds call none,
		// so far: M10a step 1.)
		cmd << ['-cflags', '-static-libgcc']
	}
	cmd << export_flags(names, join_path(temp_dir(), 'reflux_capi'))!
	// V adds the extension: .dll on Windows, .so elsewhere.
	mkdir_all(join_path('python', 'reflux', '_lib'))!
	cmd << ['-o', join_path('python', 'reflux', '_lib', 'reflux'), 'capi']
	run(cmd)
}

args := arguments()[1..]
command := args[0] or { '' }
prod := '-debug' !in args
// Run from the project root, wherever the script was started from.
chdir(dir(dir(@FILE)))!
match command {
	'build' {
		build(prod)!
	}
	'test' {
		build(prod)!
		// D12: the V API's answers for the reference states, from a program
		// built with the library's compiler and flags, for test_reference.py
		// to compare with bit for bit.
		ref_dir := join_path(temp_dir(), 'reflux_capi')
		mkdir_all(ref_dir)!
		exe := join_path(ref_dir, 'capi_reference')
		mut ref_cmd := ['v', '-cc', 'gcc', '-nocache', '-W']
		if prod {
			ref_cmd << '-prod'
		}
		ref_cmd << ['-o', exe, join_path('tools', 'capi_reference.v')]
		run(ref_cmd)
		reference := join_path(ref_dir, 'reference.json')
		run([$if windows { exe + '.exe' } $else { exe }, reference])
		setenv('REFLUX_REFERENCE', reference, true)
		mut cmd := ['uv', 'run', '--no-project']
		i := args.index('--python')
		if i >= 0 {
			cmd << ['--python', args[i + 1] or {
				eprintln('--python needs a path')
				exit(2)
			}]
		}
		cmd << ['--with', 'pytest==${pytest_version}', '--with', 'numpy==${numpy_version}', 'python',
			'-m', 'pytest', join_path('python', 'tests'), '-q']
		run(cmd)
	}
	else {
		eprintln('usage: v run tools/python.vsh build|test [-debug] [--python <path>]')
		exit(2)
	}
}
