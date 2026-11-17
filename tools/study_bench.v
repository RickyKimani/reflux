// study_bench times a study on 1, 2, 4 and 8 workers (M8 step 7, D8): the
// owner's 50-point study of the reforming PFR, 700 to 1000 °C, in chunks of
// 8 and of 2. Beside it, the same number of threads running a loop that
// allocates nothing, which says how far this machine's cores scale on
// their own, and what the study allocates a point, which the collector's
// lock serialises. Each time is the best of three.
//
//     v -prod -cc gcc -o study_bench.exe tools/study_bench.v
//     study_bench.exe
module main

import flowsheet
import math
import os
import reactors
import study
import time

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

// spin is arithmetic that allocates nothing, about as long as a point.
fn spin(n int) f64 {
	mut x := 0.0
	for i in 0 .. n {
		x += math.sin(f64(i) * 1e-3) * math.exp(-f64(i % 100) * 1e-2)
	}
	return x
}

fn spun(threads int, n int) f64 {
	t0 := time.sys_mono_now()
	mut ts := []thread f64{}
	for _ in 0 .. threads {
		ts << spawn spin(n)
	}
	ts.wait()
	return ms_since(t0)
}

fn main() {
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reformer_dwsim.toml'))!
	for chunk in [8, 2] {
		file := text.replace('format = 2', 'format = 3') +
			'\n[study.feed]\nvary = [{ key = "feed.feed.T", from = "700 degC", to = "1000 degC", points = 50 }]\nreport = { methane = "product.n[CH4]" }\nchunk = ${chunk}\n'
		l := flowsheet.load_text(file, makers: reactors.makers())!
		mut one := 0.0
		for w in [1, 2, 4, 8] {
			mut best := 1e300
			mut bytes := u64(0)
			for _ in 0 .. 3 {
				b := u64(gc_heap_usage().total_bytes)
				t0 := time.sys_mono_now()
				r := study.run(l, l.studies[0], workers: w)
				best = math.min(best, ms_since(t0))
				bytes = u64(gc_heap_usage().total_bytes) - b
				assert r.points.all(it.status == .solved)
			}
			if w == 1 {
				one = best
			}
			println('chunks of ${chunk}, ${w} workers: ${best:.0f} ms, ${one / best:.2f} times one; ${f64(bytes) / 50e6:.1f} MB a point')
		}
	}
	n := 20_000_000
	base := math.min(spun(1, n), spun(1, n))
	for w in [2, 4, 8] {
		t := math.min(spun(w, n), spun(w, n))
		println('a loop that allocates nothing on ${w} threads: ${f64(w) * base / t:.2f} times one')
	}
}
