#!/usr/bin/env -S v

// data_quality.vsh regenerates docs/data-quality.md: the sweep over every
// compound's data, the PR vapour-pressure comparison, and the outliers
// against CoolProp (M2 exit criteria 5 and 6).
//
//     v run tools/data_quality.vsh
//
// Run it after changing the data, the importer, or anything the sweep
// evaluates. The test in components/quality/quality_test.v fails until the
// committed report matches a fresh run.
import components.quality

chdir(dir(dir(@FILE)))!
text := quality.generate()!
write_file(quality.report_path, text)!
println('wrote ${quality.report_path}')
