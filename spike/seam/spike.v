// M6 step 0 spike: the flash behind a phase-model seam (D3). seamflash is
// flash's PT flash and stability test with thermo.Model replaced by an
// interface, Phases, and a cubic behind it by its volume roots (make.py
// writes it). This runs the flash robustness sweep's seeded mixtures
// through both and compares every answer bit for bit, then times both on
// M3's benchmark cases, each the least of many rounds, taking turns.
//
//     python spike/seam/make.py
//     v -prod -cc gcc run spike/seam [points]
//
// It's run from the repository root. Not part of the simulator.
module main

import benchmarks
import components
import flash
import math
import os
import spike.seam.seamflash as sf
import sweep
import thermo
import time

struct Tally {
mut:
	same_split  int
	same_single int
	same_error  int
	differ      []string
}

fn bits(a []f64) []u64 {
	return a.map(math.f64_bits(it))
}

fn same_state(a thermo.State, b thermo.State) bool {
	return math.f64_bits(a.t) == math.f64_bits(b.t) && math.f64_bits(a.p) == math.f64_bits(b.p)
		&& bits(a.n) == bits(b.n) && math.f64_bits(a.v) == math.f64_bits(b.v)
		&& math.f64_bits(a.z) == math.f64_bits(b.z) && bits(a.ln_phi) == bits(b.ln_phi)
}

// differs says how two answers differ, or '' if they're the same bits.
fn differs(a flash.FlashResult, b sf.FlashResult) string {
	if math.f64_bits(a.beta) != math.f64_bits(b.beta) {
		return 'beta ${a.beta} against ${b.beta}'
	}
	if a.phases.len != b.phases.len {
		return '${a.phases.len} phases against ${b.phases.len}'
	}
	for i, pa in a.phases {
		pb := b.phases[i]
		if int(pa.label) != int(pb.label) || math.f64_bits(pa.fraction) != math.f64_bits(pb.fraction)
			|| !same_state(pa.state, pb.state) {
			return 'phase ${i} differs'
		}
	}
	if a.iterations != b.iterations || int(a.path) != int(b.path) {
		return 'iterations ${a.iterations} against ${b.iterations}'
	}
	if a.stability.trials.len != b.stability.trials.len {
		return 'trials ${a.stability.trials.len} against ${b.stability.trials.len}'
	}
	for i, ta in a.stability.trials {
		tb := b.stability.trials[i]
		if math.f64_bits(ta.tm) != math.f64_bits(tb.tm) || bits(ta.w) != bits(tb.w)
			|| ta.iterations != tb.iterations {
			return 'trial ${i} differs'
		}
	}
	return ''
}

fn compare(db components.Database, count int) !Tally {
	pool := sweep.candidates(db)!
	mut tally := Tally{}
	for pt in sweep.points(pool, count) {
		fl := db.fluid(.pr, pt.cas) or { continue }
		m := fl.model
		eos := sf.Phases(sf.Eos{m})
		ra := flash.pt(m, pt.t, pt.p, pt.z) or {
			theirs := seam_error(eos, pt)
			if theirs == '' {
				tally.differ << "${pt.index}: flash failed (${err.msg()}) and seamflash didn't"
			} else if theirs != err.msg() {
				tally.differ << '${pt.index}: different errors'
			} else {
				tally.same_error++
			}
			continue
		}
		rb := sf.pt(eos, pt.t, pt.p, pt.z) or {
			tally.differ << '${pt.index}: seamflash failed: ${err.msg()}'
			continue
		}
		d := differs(ra, rb)
		if d != '' {
			tally.differ << '${pt.index}: ${d}'
		} else if ra.phases.len == 1 {
			tally.same_single++
		} else {
			tally.same_split++
		}
	}
	return tally
}

// seam_error is seamflash's error at a point, or '' if it gave an answer.
fn seam_error(eos sf.Phases, pt sweep.Point) string {
	sf.pt(eos, pt.t, pt.p, pt.z) or { return err.msg() }
	return ''
}

// least_ns is the least time over rounds of one pass over the states, by
// f, in ns per flash.
fn least_ns(states [][]f64, z []f64, rounds int, f fn (f64, f64, []f64) !int) !f64 {
	mut best := math.inf(1)
	for _ in 0 .. rounds {
		t0 := time.sys_mono_now()
		for s in states {
			f(s[0], s[1], z)!
		}
		best = math.min(best, f64(time.sys_mono_now() - t0) / f64(states.len))
	}
	return best
}

fn time_cases() ! {
	fx := benchmarks.load()!
	println('\nM3 benchmark cases, the least of 40 rounds each, in µs a flash:')
	println('case                         flash   seamflash  ratio')
	for c in fx.cases {
		if c.kind !in ['pt2', 'pt1'] {
			continue
		}
		s := fx.systems[c.system]
		fl := benchmarks.fluid(s)!
		m := fl.model
		eos := sf.Phases(sf.Eos{m})
		states := [][]f64{len: c.p.len, init: [c.t[index], c.p[index]]}
		direct := fn [m] (t f64, p f64, z []f64) !int {
			return flash.pt(m, t, p, z)!.phases.len
		}
		seamed := fn [eos] (t f64, p f64, z []f64) !int {
			return sf.pt(eos, t, p, z)!.phases.len
		}
		mut a, mut b := math.inf(1), math.inf(1)
		for _ in 0 .. 4 { // taking turns, ten rounds each
			a = math.min(a, least_ns(states, s.z, 10, direct)!)
			b = math.min(b, least_ns(states, s.z, 10, seamed)!)
		}
		println('${s.name + ' ' + c.kind:-26} ${a / 1e3:7.2f}  ${b / 1e3:10.2f}  ${b / a:5.3f}')
	}
}

fn main() {
	count := if os.args.len > 1 { os.args[1].int() } else { 10000 }
	db := components.load_default()!
	t0 := time.sys_mono_now()
	tally := compare(db, count)!
	println('${count} sweep points: ${tally.same_split} splits and ${tally.same_single} single phases the same bits, ${tally.same_error} failures the same error, ${tally.differ.len} different (${f64(time.sys_mono_now() - t0) / 1e9:.1f} s)')
	for d in tally.differ#[..20] {
		println('  ${d}')
	}
	time_cases()!
}
