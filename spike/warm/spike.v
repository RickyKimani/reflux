// M8 step 5 spike: what a warm start must carry (D6, as step 0 amended it).
// Each case is a sweep of one input, every point solved cold and then from
// its neighbour, carrying in turn the units' last answers, the tears, a
// secant predictor and the specs' values. It prints, a point: the passes of
// the blocks with tears, the specs' evaluations, the PFR's steps and the
// milliseconds, and how far each warm answer is from its cold one. A case's
// label, or part of it, as the argument runs that case alone.
//
// Step 5 also measured carrying what the tears' method learned, Broyden's
// inverse Jacobian and Wegstein's slopes, with code since removed: it took
// more passes, not fewer (docs/m8-plan.md, step 5). The solver now starts
// tears warm only where that pays, so on Wegstein and Broyden "tears" here
// is "cold".
//
//     v -prod -cc gcc run spike/warm
//
// It's run from the repository root. Not part of the simulator.
module main

import flowsheet
import math
import os
import reactors
import time

struct Variant {
	name    string
	units   bool
	tears   bool
	predict bool
	specs   bool
}

const variants = [
	Variant{'units', true, false, false, false},
	Variant{'tears', true, true, false, false},
	Variant{'tears, predicted', true, true, true, false},
	Variant{'specs', true, true, false, true},
	Variant{'everything', true, true, true, true},
	Variant{'specs, tears cold', true, false, true, true},
]

struct Cost {
mut:
	passes      int
	most        int
	evaluations int
	steps       int
	ms          f64
	apart       f64
}

fn ms_since(t0 u64) f64 {
	return f64(time.sys_mono_now() - t0) / 1e6
}

fn passes(s flowsheet.Solution) int {
	mut n := 0
	for h in s.history {
		n += h.residuals.len
	}
	return n
}

fn evaluations(s flowsheet.Solution) int {
	mut n := 0
	for h in s.history {
		n += h.spec_residuals.len
	}
	return n
}

fn steps(s flowsheet.Solution) int {
	mut n := 0
	for u in s.units {
		if p := u.profile {
			n += p.points.len - 1
		}
	}
	return n
}

// apart is how far two solutions are: the largest difference in a stream's
// flows as a fraction of its total, or in its T or P as a fraction of them.
fn apart(a flowsheet.Solution, b flowsheet.Solution) f64 {
	mut worst := 0.0
	for name, x in a.streams {
		y := b.streams[name] or { continue }
		total := math.max(x.total(), 1e-300)
		for i, f in x.flows {
			worst = math.max(worst, math.abs(f - y.flows[i]) / total)
		}
		worst = math.max(worst, math.abs(x.t - y.t) / x.t)
		worst = math.max(worst, math.abs(x.p - y.p) / x.p)
	}
	return worst
}

fn add(mut c Cost, s flowsheet.Solution, ms f64) {
	c.passes += passes(s)
	c.most = math.max(c.most, passes(s))
	c.evaluations += evaluations(s)
	c.steps += steps(s)
	c.ms += ms
}

fn line(name string, c Cost, n int) string {
	m := f64(n)
	return '    ${name:-26} passes ${c.passes / m:5.1f} (most ${c.most:2}), spec evaluations ${c.evaluations / m:5.1f}, steps ${c.steps / m:6.1f}, ${c.ms / m:7.2f} ms; from cold ${c.apart:.1e}'
}

// sweep solves each value of a key cold, then from its neighbour by each
// variant, the first point of each cold.
fn sweep(label string, text string, key string, values []f64) ! {
	if os.args.len > 1 && !label.contains(os.args[1]) {
		return
	}
	l := flowsheet.load_text(text, makers: reactors.makers())!
	input := l.input(key)!
	mut colds := []flowsheet.Solution{}
	mut cold := Cost{}
	for v in values {
		lv := l.with([input], [v])!
		t0 := time.sys_mono_now()
		sol := lv.solve()!
		add(mut cold, sol, ms_since(t0))
		colds << sol
	}
	n := values.len - 1
	println('  ${label}: ${values.len} points of ${key}')
	println(line('cold', Cost{
		...cold
		passes:      cold.passes - passes(colds[0])
		evaluations: cold.evaluations - evaluations(colds[0])
		steps:       cold.steps - steps(colds[0])
		ms:          cold.ms * f64(n) / f64(values.len)
	}, n))
	for w in variants {
		mut c := Cost{}
		mut sols := [colds[0]]
		for k in 1 .. values.len {
			lv := l.with([input], [values[k]])!
			start := flowsheet.Start{
				last:    sols[k - 1]
				before:  if k >= 2 { ?flowsheet.Solution(sols[k - 2]) } else { none }
				units:   w.units
				tears:   w.tears
				predict: w.predict
				specs:   w.specs
			}
			t0 := time.sys_mono_now()
			mut failure := ''
			sol := lv.solve_from(start) or {
				failure = err.msg()
				flowsheet.Solution{}
			}
			if failure != '' {
				println('    ${w.name}: point ${k} failed: ${failure}')
				break
			}
			add(mut c, sol, ms_since(t0))
			c.apart = math.max(c.apart, apart(sol, colds[k]))
			sols << sol
		}
		println(line(w.name, c, n))
	}
}

fn spaced(from f64, step f64, n int) []f64 {
	return []f64{len: n, init: from + f64(index) * step}
}

fn run() ! {
	root := 'testdata/flowsheets'
	reference := os.read_file('${root}/reference.toml')!
	for method in ['wegstein', 'broyden', 'direct'] {
		text := reference.replace('method = "wegstein"', 'method = "${method}"')
		sweep('reference, ${method}, 1 K apart', text, 'unit.cooler.outlet_T', spaced(230,
			-1, 21))!
		if method != 'direct' {
			sweep('reference, ${method}, 0.1 K apart', text, 'unit.cooler.outlet_T',
				spaced(230, -0.1, 21))!
			sweep('reference, ${method}, 5 K apart', text, 'unit.cooler.outlet_T', spaced(230,
				-5, 9))!
		}
	}
	reformer := os.read_file('${root}/reformer_dwsim.toml')!
	sweep("the owner's PFR", reformer, 'feed.feed.T', spaced(973.15, 300.0 / 47, 48))!
	pfr := reformer[reformer.index('[unit.reformer]') or { 0 }..reformer.index('[site.Ni]') or {
		0
	}]
	cstr := reformer.replace(pfr, '[unit.reformer]\ntype = "cstr"\nin = "feed"\nout = "product"\nvolume = "1 m^3"\ncatalyst_loading = "0.386 kg/m^3"\n\n')
	sweep('a CSTR on its kinetics', cstr, 'feed.feed.T', spaced(973.15, 300.0 / 23, 24))!
	gibbs := reformer.replace(pfr, '[unit.reformer]\ntype = "gibbs_reactor"\nin = "feed"\nout = "product"\nadiabatic = true\n\n')
	sweep('an adiabatic Gibbs reactor', gibbs, 'feed.feed.T', spaced(973.15, 300.0 / 23,
		24))!
	heater := 'format = 3\n\n[components]\nnames = ["methane", "propane", "n-pentane"]\n\n[model]\nequation = "PR"\n\n[feed.f]\nT = "300 K"\nP = "10 bar"\nflow = "100 mol/s"\nmole_fractions = { methane = 0.6, propane = 0.3, n-pentane = 0.1 }\n\n[unit.heat]\ntype = "heater"\nin = "f"\nout = "out"\nduty = "0 W"\n\n[spec.cool]\nvary = "unit.heat.duty"\nbetween = ["-2 MW", "1 MW"]\nuntil = "out.T = 250 K"\n'
	sweep("a heater's spec", heater, 'feed.f.T', spaced(280, 2, 21))!
	spec := reference.replace('format = 1', 'format = 3') +
		'\n[spec.vapour]\nvary = "unit.cooler.outlet_T"\nbetween = ["200 K", "260 K"]\nuntil = "gas.F = 60 mol/s"\n'
	sweep('the reference flowsheet, with a spec', spec, 'feed.feed.P', spaced(30e5, 0.1e5,
		11))!
}

fn main() {
	run() or {
		eprintln(err.msg())
		exit(1)
	}
}
