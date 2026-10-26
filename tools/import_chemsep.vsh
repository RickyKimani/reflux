#!/usr/bin/env -S v

// import_chemsep.vsh regenerates data/components/ from the vendored ChemSep
// XML in data/vendor/chemsep/.
//
//     v run tools/import_chemsep.vsh
//
// It writes every file components.chemsep.generate() produces and removes
// any file in the directory that it no longer produces, so the directory is
// exactly the generated set. Nothing else writes there. The test in
// components/chemsep/generate_test.v fails while the committed directory
// differs from what this would write, so CI catches a forgotten rerun.
import components.chemsep

chdir(dir(dir(@FILE)))!

imp := chemsep.import_files('.', chemsep.vendored)!
files := chemsep.generate(imp)!

out := chemsep.output_dir
mkdir_all(out)!
mut written, mut unchanged := 0, 0
for name, text in files {
	path := join_path(out, name)
	if exists(path) && read_file(path)! == text {
		unchanged++
		continue
	}
	write_file(path, text)!
	written++
}
mut removed := 0
for name in ls(out)! {
	if name !in files {
		rm(join_path(out, name))!
		removed++
		println('removed ${name}')
	}
}

println('${imp.records.len} compounds from ${chemsep.vendored.len} files: ${written} files written, ${unchanged} unchanged, ${removed} removed')
println('${imp.repairs.len} repairs, ${imp.duplicates.len} duplicates dropped, ${imp.findings.len} findings: see ${out}/CHANGES.md')
