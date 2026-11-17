// M8 step 0 spike: what a study's warm starts save, whether threads sharing
// one loaded file give one thread's answers, and what a solve allocates.
//
//     v -prod -cc gcc run spike/study
//
// 1. Warm starts by hand. The reference flowsheet's cooler is swept, and each
//    point is solved three ways: cold; its tear started from the last point's
//    answer, through SolveParams.estimates; and from a secant predictor
//    through the last two. The passes each takes, and how far the warm
//    answers are from the cold one, set criterion 4's expectations.
// 2. Workers. The owner's study (testdata/flowsheets/reformer_dwsim.toml, 48
//    feed temperatures) on 1, 2, 4 and 8 threads sharing one Loaded: every
//    answer's bits against one thread's, and the time. Then a checksum of
//    every number of that study, of the published model's (reformer.toml)
//    and of both with their PFRs adiabatic, where T moves at every state:
//    two builds of reflux with the same checksums give the same answers,
//    which is how M8 step 1 showed that what it keeps changes none.
// 3. What a solve allocates, by the collector's count, for the owner's case
//    and the published model.
//
// It's run from the repository root. Not part of the simulator.
module main

import flowsheet
import math
import os
import reactors
import time

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

// with_cooler is the reference flowsheet with its cooler's outlet at t.
fn with_cooler(l flowsheet.Loaded, t f64) !flowsheet.Flowsheet {
	mut units := l.sheet.units.clone()
	for i, p in units {
		if p.unit.name == 'cooler' {
			h := flowsheet.new_heater('cooler', mode: .outlet_t, value: t)!
			units[i] = flowsheet.Placed{
				unit:    h
				inlets:  p.inlets
				outlets: p.outlets
			}
		}
	}
	return flowsheet.Flowsheet{
		units: units
		feeds: l.sheet.feeds
	}
}

// passes is how many passes a solve's blocks with tears took.
fn passes(s flowsheet.Solution) int {
	mut n := 0
	for h in s.history {
		n += h.residuals.len
	}
	return n
}

// apart is how far two solutions are: the largest difference in any stream's
// flows as a fraction of its total, or in its T or P as a fraction of them.
fn apart(a flowsheet.Solution, b flowsheet.Solution) f64 {
	mut worst := 0.0
	for name, x in a.streams {
		y := b.streams[name]
		total := math.max(x.total(), 1e-300)
		for i, f in x.flows {
			worst = math.max(worst, math.abs(f - y.flows[i]) / total)
		}
		worst = math.max(worst, math.abs(x.t - y.t) / x.t)
		worst = math.max(worst, math.abs(x.p - y.p) / x.p)
	}
	return worst
}

// predicted is a tear's first value from the last two answers, b after a:
// 2b − a in its flows, P and molar H, a flow held at 0.
fn predicted(l flowsheet.Loaded, a flowsheet.Stream, b flowsheet.Stream) !flowsheet.Stream {
	mut flows := []f64{len: b.flows.len}
	for i in 0 .. flows.len {
		flows[i] = math.max(2 * b.flows[i] - a.flows[i], 0)
	}
	return flowsheet.feed_ph(l.fluid, flows, 2 * b.p - a.p, 2 * b.h - a.h)
}

struct Timed {
	sol flowsheet.Solution
	ms  f64
}

// timed solves three times and keeps the quickest.
fn timed(sheet flowsheet.Flowsheet, l flowsheet.Loaded, params flowsheet.SolveParams) !Timed {
	mut best := 1e300
	mut sol := flowsheet.Solution{}
	for _ in 0 .. 3 {
		t0 := time.sys_mono_now()
		sol = sheet.solve(l.fluid, params)!
		best = math.min(best, ms_since(t0))
	}
	return Timed{sol, best}
}

// sweep solves n points from t0, dt apart, each cold, warm and predicted.
fn sweep(l flowsheet.Loaded, method flowsheet.Method, tol f64, t0 f64, dt f64, n int) ! {
	mut cold_n, mut warm_n, mut pred_n := 0, 0, 0
	mut cold_ms, mut warm_ms, mut pred_ms := 0.0, 0.0, 0.0
	mut warm_apart, mut pred_apart := 0.0, 0.0
	mut warm_most, mut pred_most := 0, 0
	mut last := []flowsheet.Solution{}
	mut tear := ''
	for k in 0 .. n {
		sheet := with_cooler(l, t0 + k * dt)!
		params := flowsheet.SolveParams{
			...l.params
			method: method
			tol:    tol
		}
		cold := timed(sheet, l, params)!
		tear = cold.sol.history[0].tears[0]
		if k >= 2 {
			// The first two points have no neighbours to start from.
			cold_n += passes(cold.sol)
			cold_ms += cold.ms
			warm := timed(sheet, l, flowsheet.SolveParams{
				...params
				estimates: {
					tear: last[k - 1].streams[tear]
				}
			})!
			warm_n += passes(warm.sol)
			warm_most = math.max(warm_most, passes(warm.sol))
			warm_ms += warm.ms
			warm_apart = math.max(warm_apart, apart(warm.sol, cold.sol))
			start := predicted(l, last[k - 2].streams[tear], last[k - 1].streams[tear])!
			pred := timed(sheet, l, flowsheet.SolveParams{
				...params
				estimates: {
					tear: start
				}
			})!
			pred_n += passes(pred.sol)
			pred_most = math.max(pred_most, passes(pred.sol))
			pred_ms += pred.ms
			pred_apart = math.max(pred_apart, apart(pred.sol, cold.sol))
		}
		last << cold.sol
	}
	m := f64(n - 2)
	println('  ${method:-8} tol ${tol:.0e}, ${n} points ${dt:+.1f} K apart: passes a point, cold ${cold_n / m:.1f}, warm ${warm_n / m:.1f} (at most ${warm_most}), predicted ${pred_n / m:.1f} (at most ${pred_most}); ms a point ${cold_ms / m:.2f}, ${warm_ms / m:.2f}, ${pred_ms / m:.2f}; warm answers within ${warm_apart:.1e} of cold, predicted ${pred_apart:.1e} (tear ${tear})')
}

// Out is a point's answer as bits, so that two can be compared exactly.
struct Out {
	bits  []u64
	steps int
	err   string
}

// solve_points solves the flowsheet with its feeds at each T in turn.
fn solve_points(l flowsheet.Loaded, ts []f64) []Out {
	mut out := []Out{}
	for t in ts {
		out << solve_at(l, t)
	}
	return out
}

fn solve_at(l flowsheet.Loaded, t f64) Out {
	mut feeds := []flowsheet.Feed{}
	for f in l.sheet.feeds {
		s := flowsheet.feed_tp(l.fluid, f.stream.flows, t, f.stream.p) or {
			return Out{
				err: err.msg()
			}
		}
		feeds << flowsheet.Feed{f.name, s}
	}
	sheet := flowsheet.Flowsheet{
		units: l.sheet.units
		feeds: feeds
	}
	sol := sheet.solve(l.fluid, l.params) or {
		return Out{
			err: err.msg()
		}
	}
	mut bits := []u64{}
	mut steps := 0
	for name in sol.products {
		s := sol.streams[name]
		for f in s.flows {
			bits << math.f64_bits(f)
		}
		bits << math.f64_bits(s.t)
		bits << math.f64_bits(s.p)
		bits << math.f64_bits(s.h)
	}
	for u in sol.units {
		for e in u.energy {
			bits << math.f64_bits(e.rate)
		}
		if p := u.profile {
			steps += p.points.len - 1
			for pt in p.points {
				bits << math.f64_bits(pt.volume)
				bits << math.f64_bits(pt.t)
				bits << math.f64_bits(pt.p)
				for f in pt.flows {
					bits << math.f64_bits(f)
				}
			}
		}
	}
	return Out{
		bits:  bits
		steps: steps
	}
}

// same says whether two runs gave the same bits at every point, and how many
// points and numbers that was.
fn same(a []Out, b []Out) (bool, int) {
	mut numbers := 0
	if a.len != b.len {
		return false, 0
	}
	for k in 0 .. a.len {
		if a[k].err != b[k].err || a[k].steps != b[k].steps || a[k].bits != b[k].bits {
			return false, numbers
		}
		numbers += a[k].bits.len
	}
	return true, numbers
}

// adiabatic is a loaded flowsheet with its PFRs adiabatic, so that T moves
// at every state along them and K is computed at each.
fn adiabatic(l flowsheet.Loaded) flowsheet.Loaded {
	mut units := l.sheet.units.clone()
	for i, p in units {
		u := p.unit
		if u is reactors.Pfr {
			// A struct update on the smart-cast u would copy nothing (V 0.5.2).
			pfr := u as reactors.Pfr
			units[i] = flowsheet.Placed{
				unit:    reactors.Pfr{
					...pfr
					thermal: .adiabatic
				}
				inlets:  p.inlets
				outlets: p.outlets
			}
		}
	}
	return flowsheet.Loaded{
		...l
		sheet: flowsheet.Flowsheet{
			units: units
			feeds: l.sheet.feeds
		}
	}
}

// checksum is one number from every bit of a run's answers (FNV-1a), to
// compare two builds of reflux by: the same checksum, the same answers.
fn checksum(outs []Out) u64 {
	mut h := u64(14695981039346656037)
	for o in outs {
		for b in o.bits {
			for k in 0 .. 8 {
				h = (h ^ ((b >> (8 * k)) & 0xff)) * 1099511628211
			}
		}
		h = (h ^ u64(o.steps)) * 1099511628211
	}
	return h
}

fn workers(l flowsheet.Loaded, ts []f64) {
	mut t0 := time.sys_mono_now()
	serial := solve_points(l, ts)
	mut one := ms_since(t0)
	for _ in 0 .. 2 {
		t0 = time.sys_mono_now()
		again := solve_points(l, ts)
		one = math.min(one, ms_since(t0))
		ok, _ := same(again, serial)
		assert ok, 'one thread, run twice, gave different bits'
	}
	failed := serial.filter(it.err != '').len
	mut steps := 0
	for o in serial {
		steps += o.steps
	}
	println('  one thread: ${ts.len} points in ${one:.0f} ms, ${failed} failed, ${steps} integrator steps, checksum ${checksum(serial):016x}')
	for w in [2, 4, 8] {
		mut best := 1e300
		mut all_same := true
		mut numbers := 0
		for _ in 0 .. 3 {
			t0 = time.sys_mono_now()
			mut threads := []thread []Out{}
			for k in 0 .. w {
				// Thread k takes every w-th point from the k-th.
				mut mine := []f64{}
				for i := k; i < ts.len; i += w {
					mine << ts[i]
				}
				threads << spawn solve_points(l, mine)
			}
			parts := threads.wait()
			best = math.min(best, ms_since(t0))
			mut together := []Out{len: ts.len}
			for k in 0 .. w {
				for j, o in parts[k] {
					together[k + j * w] = o
				}
			}
			ok, n := same(together, serial)
			all_same = all_same && ok
			numbers = n
		}
		println('  ${w} threads: ${best:.0f} ms (${one / best:.2f} times one), every bit of ${numbers} numbers the same as one thread, three runs: ${all_same}')
	}
}

// allocated solves the file's flowsheet and counts what the collector handed
// out for it.
fn allocated(path string) ! {
	l := flowsheet.load(path, makers: reactors.makers())!
	l.solve()!
	before := gc_heap_usage().total_bytes
	t0 := time.sys_mono_now()
	sol := l.solve()!
	ms := ms_since(t0)
	bytes := gc_heap_usage().total_bytes - before
	mut steps := 0
	for u in sol.units {
		if p := u.profile {
			steps += p.points.len - 1
		}
	}
	println('  ${os.file_name(path)}: a solve allocates ${f64(bytes) / 1e6:.2f} MB in ${ms:.1f} ms, over ${steps} integrator steps: ${f64(bytes) / 1e3 / steps:.1f} kB a step')
}

fn run() ! {
	$if !prod {
		println('not built with -prod: the times below mean little')
	}
	println('1. Warm starts: the reference flowsheet, its cooler swept from 230 K')
	reference := flowsheet.load('testdata/flowsheets/reference.toml', makers: reactors.makers())!
	for method in [flowsheet.Method.wegstein, .direct, .broyden] {
		sweep(reference, method, reference.params.tol, 230, -1, 21)!
	}
	sweep(reference, .wegstein, reference.params.tol, 230, -0.1, 21)!
	sweep(reference, .wegstein, reference.params.tol, 230, -5, 9)!
	sweep(reference, .broyden, reference.params.tol, 230, -5, 9)!
	sweep(reference, .wegstein, 1e-12, 230, -1, 21)!

	println("2. Workers: the owner's study, 48 feed temperatures from 700 to 1000 °C")
	owner := flowsheet.load('testdata/flowsheets/reformer_dwsim.toml', makers: reactors.makers())!
	ts := []f64{len: 48, init: 973.15 + 300.0 * index / 47}
	workers(owner, ts)

	published := flowsheet.load('testdata/flowsheets/reformer.toml', makers: reactors.makers())!
	t0 := time.sys_mono_now()
	model := solve_points(published, ts)
	println('  the published model (reformer.toml), the same 48 points on one thread: ${ms_since(t0):.0f} ms, ${model.filter(it.err != '').len} failed, checksum ${checksum(model):016x}')
	for name, sheet in {
		'the published model, adiabatic': adiabatic(published)
		"the owner's case, adiabatic":    adiabatic(owner)
	} {
		t1 := time.sys_mono_now()
		outs := solve_points(sheet, ts)
		mut steps := 0
		for o in outs {
			steps += o.steps
		}
		println('  ${name}: ${ms_since(t1):.0f} ms, ${outs.filter(it.err != '').len} failed, ${steps} integrator steps, checksum ${checksum(outs):016x}')
	}
	println('3. What a solve allocates')
	allocated('testdata/flowsheets/reformer_dwsim.toml')!
	allocated('testdata/flowsheets/reformer.toml')!
}

fn main() {
	run() or {
		eprintln(err.msg())
		exit(1)
	}
}
