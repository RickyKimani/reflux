#!/usr/bin/env -S v

// flash_sweep.vsh runs M3's robustness sweep (decision D8, exit criterion
// 5) and writes docs/flash-robustness.md:
//
//     v -prod -cc gcc run tools/flash_sweep.vsh        # all 10,000 points
//     v -prod -cc gcc run tools/flash_sweep.vsh 1000   # just the first 1,000
//
// With gcc and -prod it takes about three minutes, the build included, and
// under tcc about twenty. sweep/sweep_test.v runs the first 1,000 points in
// CI.
import components
import sweep

chdir(dir(dir(@FILE)))!
args := arguments()
count := if args.len > 1 { args[1].int() } else { sweep.full_count }
db := components.load_default()!
pool := sweep.candidates(db)!
outcomes := sweep.sweep(db, count, fn (i int, o sweep.Outcome) {
	if (i + 1) % 500 == 0 {
		println('${i + 1} points')
	}
	if o.class == .unexplained {
		println('point ${i}: unexplained: ${o.detail}')
	}
})!
write_file(sweep.report_path, sweep.report(pool.len, outcomes))!
t := sweep.tally(outcomes)
println('wrote ${sweep.report_path}: ${t}')
