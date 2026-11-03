#!/usr/bin/env -S v

// golden.vsh writes the golden flowsheets' results (M4, D13): each flowsheet's
// report, as `reflux solve --json` writes it, in testdata/flowsheets/golden/.
// flowsheet/golden_test.v solves each again and compares. Run it only on
// purpose, when a change is meant to move the answers, and say why in the
// commit:
//
//     v run tools/golden.vsh
//
// The numbers' last bits can differ between C compilers, so it's run with
// V's default, tcc, as the test is; the test allows 1e-9 of each value's scale.
import flowsheet

const golden = ['reference', 'every_unit', 'no_recycle', 'interlocking', 'saturated_loop']

chdir(dir(dir(@FILE)))!
mkdir_all(join_path('testdata', 'flowsheets', 'golden'))!
for name in golden {
	// Forward slashes, so that the report reads the same on every system.
	path := 'testdata/flowsheets/${name}.toml'
	loaded := flowsheet.load(path)!
	report := flowsheet.Report{
		...loaded.report(loaded.solve()!)
		file: path
	}
	out := join_path('testdata', 'flowsheets', 'golden', '${name}.json')
	write_file(out, report.json())!
	println('wrote ${out}')
}
