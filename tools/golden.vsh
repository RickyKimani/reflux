#!/usr/bin/env -S v

// golden.vsh writes the golden flowsheets' results (M4, D13): each flowsheet's
// report, as `reflux solve --json` writes it, in testdata/flowsheets/golden/,
// and the flowsheets with reactors (M9 step 10) in its reactors/ folder.
// flowsheet/golden_test.v and reactors/golden_test.v solve each again and
// compare. Run it only on purpose, when a change is meant to move the
// answers, and say why in the commit:
//
//     v run tools/golden.vsh [name ...]
//
// With names, it writes only those flowsheets; with none, every one.
//
// The numbers' last bits can differ between C compilers, so it's run with
// V's default, tcc, as the test is; the test allows 1e-9 of each value's scale.
import flowsheet
import reactors

const golden = ['reference', 'every_unit', 'no_recycle', 'interlocking', 'saturated_loop']
const with_reactors = ['reformer', 'methanol']

chdir(dir(dir(@FILE)))!
asked := arguments()[1..]
for name in asked {
	if name !in golden && name !in with_reactors {
		eprintln('no golden flowsheet ${name}: ${golden.join(', ')}, ${with_reactors.join(', ')}')
		exit(2)
	}
}
mut every := golden[..].clone()
every << with_reactors[..]
for name in every {
	if asked.len > 0 && name !in asked {
		continue
	}
	// Forward slashes, so that the report reads the same on every system.
	path := 'testdata/flowsheets/${name}.toml'
	reacts := name in with_reactors
	loaded := if reacts {
		flowsheet.load(path, makers: reactors.makers())!
	} else {
		flowsheet.load(path)!
	}
	report := flowsheet.Report{
		...loaded.report(loaded.solve()!)
		file: path
	}
	dir_ := if reacts {
		join_path('testdata', 'flowsheets', 'golden', 'reactors')
	} else {
		join_path('testdata', 'flowsheets', 'golden')
	}
	mkdir_all(dir_)!
	out := join_path(dir_, '${name}.json')
	write_file(out, report.json())!
	println('wrote ${out}')
}
