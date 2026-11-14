// bench_reactor times reflux on a flowsheet with a reactor, for M9 step 11
// (criterion 9). tools/gen_m9_bench.py starts it once and hands it one feed
// temperature at a time, between DWSIM's solves of the same, so that the
// two programs are timed within seconds of each other:
//
//     v -prod -cc gcc -o bench_reactor.exe tools/bench_reactor.v
//     bench_reactor.exe --build          # how it was built
//     bench_reactor.exe file.toml        # then "T repeats" on each line of stdin
//
// It answers the file's loading first, then each line with a JSON object:
// each solve's time, the feed flashed at T and the flowsheet solved, as a
// study's point costs; the feed's and the product's flows; and the steps
// the PFR's integrator took. It's a .v program, not a .vsh script, because
// a script runs whenever V compiles it (docs/v-for-go.md, §1).
module main

import benchmarks
import flowsheet
import json2
import os
import reactors
import time

struct Load {
	load_ms f64
}

struct Point {
	t       f64
	ms      []f64
	feed    []f64
	product []f64
	steps   int
	error   string
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

// solve_at solves the flowsheet with every feed at T, repeats times.
fn solve_at(l flowsheet.Loaded, t f64, repeats int) Point {
	mut ms := []f64{}
	mut last := flowsheet.Solution{}
	mut feeds := []flowsheet.Feed{}
	for _ in 0 .. repeats {
		t0 := time.sys_mono_now()
		feeds = []flowsheet.Feed{}
		for f in l.sheet.feeds {
			s := flowsheet.feed_tp(l.fluid, f.stream.flows, t, f.stream.p) or {
				return Point{
					t:     t
					error: err.msg()
				}
			}
			feeds << flowsheet.Feed{f.name, s}
		}
		sheet := flowsheet.Flowsheet{
			units: l.sheet.units
			feeds: feeds
		}
		last = sheet.solve(l.fluid, l.params) or {
			return Point{
				t:     t
				error: err.msg()
			}
		}
		ms << ms_since(t0)
	}
	mut steps := 0
	for u in last.units {
		if p := u.profile {
			steps += p.points.len - 1
		}
	}
	return Point{
		t:       t
		ms:      ms
		feed:    feeds[0].stream.flows
		product: last.streams[last.products[0]].flows
		steps:   steps
	}
}

fn main() {
	$if !prod {
		eprintln('build with -prod: v -prod -cc gcc -o bench_reactor.exe tools/bench_reactor.v')
		exit(1)
	}
	if os.args.len == 2 && os.args[1] == '--build' {
		println(json2.encode(benchmarks.Build{
			compiler: benchmarks.compiler()
			v:        @VCURRENTHASH
		}))
		return
	}
	if os.args.len != 2 {
		eprintln('usage: bench_reactor --build, or bench_reactor file.toml with "T repeats" lines on stdin')
		exit(2)
	}
	t0 := time.sys_mono_now()
	loaded := flowsheet.load(os.args[1], makers: reactors.makers()) or {
		eprintln(err.msg())
		exit(1)
	}
	println(json2.encode(Load{ms_since(t0)}))
	for {
		line := os.get_line().trim_space()
		if line == '' {
			break
		}
		fields := line.split(' ')
		repeats := if fields.len > 1 { fields[1].int() } else { 1 }
		println(json2.encode(solve_at(loaded, fields[0].f64(), repeats)))
	}
}
