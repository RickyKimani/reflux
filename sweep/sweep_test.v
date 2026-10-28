// The robustness sweep's CI part (exit criterion 5): the first 1,000 of the
// seeded points that tools/flash_sweep.vsh runs 10,000 of, with no
// unexplained failure; and the committed report of the full sweep, which
// must show none either.
import components
import math
import numeric
import os
import sweep

fn test_the_first_thousand_points() {
	db := components.load_default()!
	outcomes := sweep.sweep(db, sweep.ci_count, fn (i int, o sweep.Outcome) {
		if o.class == .unexplained {
			eprintln('point ${i}: ${o.point.names} at ${o.point.t} K, ${o.point.p} Pa: ${o.detail}')
		}
	})!
	t := sweep.tally(outcomes)
	println(t)
	assert t.unexplained == 0
	assert t.single + t.split + t.three_phase + t.near_critical == sweep.ci_count
}

// The points of the full sweep, beyond the first 1,000, that found the
// flash's bugs in step 9, among them: K-values that became 0 or +inf in
// Rachford–Rice (1523 and 1879), a split that ran out of steps from one
// start (7301), a third phase hidden behind a later start's failure (1269
// and 1989), and a Newton step that never ended (3430, when DEM's K went
// unclamped). Points 322 and 449 are in the CI subset.
const found_bugs = [1269, 1523, 1879, 1989, 2423, 3119, 3430, 3613, 5116, 7301, 8792, 9175, 9497,
	9961]

fn test_the_points_that_found_bugs() {
	db := components.load_default()!
	all := sweep.points(sweep.candidates(db)!, sweep.full_count)
	for i in found_bugs {
		o := sweep.run(db, all[i])
		assert o.class != .unexplained, 'point ${i}: ${o.detail}'
	}
}

// The points where the first start's split has an unstable phase, a local
// minimum of G, and a later start finds the verified split. A flash that
// refused at the first start would call each of them three phases.
const later_start_splits = [209, 554, 1822, 2261, 2584, 5333, 5433, 5772, 7642, 8089, 9311]

fn test_a_later_start_finds_the_split() {
	db := components.load_default()!
	all := sweep.points(sweep.candidates(db)!, sweep.full_count)
	for i in later_start_splits {
		o := sweep.run(db, all[i])
		assert o.class == .split, 'point ${i}: ${o.class} ${o.detail}'
	}
}

fn test_the_committed_report_shows_no_unexplained_failure() {
	// Regenerating 10,000 points takes too long for CI, so this reads what
	// the last full run recorded. Its counts row for unexplained failures
	// must be zero for all of them and for the CI subset.
	path := os.join_path(os.dir(@FILE), '..', sweep.report_path)
	text := os.read_file(path) or { panic('${path}: ${err}. Run `v -prod -cc gcc run tools/flash_sweep.vsh`.') }
	assert text.contains('| Outcome | All ${sweep.full_count} |'), 'the report is not of the full sweep'
	assert text.contains('| **Unexplained** | **0** | **0** |'), 'the full sweep has unexplained failures'
}

fn test_the_points_are_seeded() {
	// splitmix64 from seed 0 gives 0xe220a8397b1dcdaf first (Vigna's
	// reference), and the same seed gives the same points every time.
	mut r := sweep.Rng{0}
	assert r.u64() == u64(0xe220a8397b1dcdaf)
	assert r.u64() == u64(0x6e789e6aa1b965f4)
	db := components.load_default()!
	pool := sweep.candidates(db)!
	a := sweep.points(pool, 20)
	b := sweep.points(pool, 20)
	for k in 0 .. 20 {
		assert a[k].cas == b[k].cas && a[k].z == b[k].z && a[k].t == b[k].t && a[k].p == b[k].p
		assert a[k].t >= sweep.t_min && a[k].t <= sweep.t_max
		assert a[k].p >= sweep.p_min && a[k].p <= sweep.p_max
		assert a[k].cas.len >= sweep.nc_min && a[k].cas.len <= sweep.nc_max
		// Ten fractions summed carry about 1e-15 of rounding.
		assert a[k].z.all(it > 0.0) && math.abs(numeric.sum(a[k].z) - 1.0) <= 1e-14
	}
}
