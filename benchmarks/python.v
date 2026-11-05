module benchmarks

import components
import flash
import json2
import numeric
import os
import strings
import thermo
import time

// M10a's speed (exit criterion 10): what a call from Python costs over V's
// own loop doing the same work, for properties and a PT flash at 1, 100 and
// 10,000 points. tools/gen_m10a_bench.py times V's loop through
// tools/bench_case.v and the C API through the Python package, one after the
// other, and writes testdata/m10a_bench.json; tools/benchmarks.vsh adds it to
// docs/benchmarks.md.

pub const python_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm10a_bench.json')

// The work timed, and the numbers of points.
pub const python_work = ['props', 'flash_pt']
pub const python_points = [1, 100, 10000]

// CallTiming is one piece of work at one number of points: each call's time,
// in µs, from V's own loop and from Python through the C API.
pub struct CallTiming {
pub:
	what      string // props or flash_pt
	points    int
	v_us      []f64
	python_us []f64
}

pub struct PythonSource {
pub:
	generated string
	generator string
	python    string
	numpy     string
}

pub struct PythonFixture {
pub:
	source       PythonSource
	machine      Machine
	reflux_build Build
	timings      []CallTiming
}

// load_python reads the fixture.
pub fn load_python() !PythonFixture {
	text := os.read_file(python_fixture_path) or {
		return error('${python_fixture_path}: ${err}. Regenerate it: `uv run --python <3.11 or later> tools/gen_m10a_bench.py`.')
	}
	return json2.decode[PythonFixture](text)!
}

// python_states are the points timed: methane, ethane and n-butane, 30/30/40,
// at 20 bar, from 220 to 400 K, two phases and one. One point is at 300 K.
pub fn python_states(n int) []f64 {
	if n == 1 {
		return [300.0]
	}
	return []f64{len: n, init: 220.0 + 180.0 * f64(index) / f64(n - 1)}
}

pub const python_p = 20e5
pub const python_z = [0.3, 0.3, 0.4]

// time_v times V's own loop over n points, reps times, doing the work the C
// API does at each: for props, the stable volume root, the properties per
// mole and the phase's label; for flash_pt, the PT flash and each phase's
// properties. Each call's time is in µs.
pub fn time_v(what string, n int, reps int) ![]f64 {
	if what !in python_work {
		return error('no work ${what}: ${python_work}')
	}
	db := components.load_default()!
	fl := db.fluid(.pr, ['methane', 'ethane', 'n-butane'])!
	ts := python_states(n)
	mut sink := 0.0
	mut us := []f64{}
	for _ in 0 .. reps + 1 {
		t0 := time.sys_mono_now()
		for t in ts {
			if what == 'props' {
				st := thermo.solve_tp(fl.model, t, python_p, python_z, .stable)!
				pr := fl.props(st)
				sink += pr.h / numeric.sum(python_z)
				if thermo.phase_label(fl.model, st) == .liquid {
					sink += 1.0
				}
			} else {
				r := flash.pt(fl.model, t, python_p, python_z)!
				for ph in r.phases {
					sink += ph.fraction * fl.props(ph.state).h
				}
			}
		}
		us << f64(time.sys_mono_now() - t0) / 1e3
	}
	if !(sink == sink) {
		return error('NaN in the work')
	}
	// The first call warms the caches; it isn't one of the reps.
	return us[1..]
}

// long_duration prints µs as duration does, and in seconds from a second up.
fn long_duration(us f64) string {
	return if round_sig(us, 3) < 1e6 { duration(us) } else { '${sig(us / 1e6, 3)} s' }
}

// signed_duration prints a difference of µs, which noise can make negative.
fn signed_duration(us f64) string {
	return if us < 0.0 { '−' + long_duration(-us) } else { long_duration(us) }
}

fn work_label(what string) string {
	return if what == 'props' { 'Properties' } else { 'PT flash' }
}

// python_section is M10a's part of docs/benchmarks.md.
pub fn python_section(pf PythonFixture) string {
	mut b := strings.new_builder(2048)
	b.writeln('')
	b.writeln('## From Python')
	b.writeln('')
	b.writeln("M10a step 7, exit criterion 10: what a call from Python costs over V's own loop doing the same work. Methane, ethane and n-butane, 30/30/40, at 20 bar, from 220 to 400 K (one point at 300 K), with Peng–Robinson and the database's kᵢⱼ. `tools/gen_m10a_bench.py` writes `testdata/m10a_bench.json`, timing V and Python alternately, on ${pf.machine.cpu}, ${pf.machine.os}, on ${pf.source.generated}, with Python ${pf.source.python} and NumPy ${pf.source.numpy}.")
	b.writeln('')
	b.writeln("- **V's loop** is `tools/bench_case.v`, built with ${pf.reflux_build.compiler}, doing at each point what the C API does: for properties, the stable volume root, the properties per mole and the phase's label; for the PT flash, the flash and each phase's properties.")
	b.writeln("- **Each pair is timed in eight rounds,** V then Python, about four seconds in all, so that the laptop's power state moves both alike.")
	b.writeln('- **From Python** is one call of `Fluid.props` or `Fluid.flash_pt` with every point, through the library built with the same compiler, its arrays allocated and its results shaped as the package does.')
	b.writeln('- **The overhead** is the difference of the medians: per call at one point, and per point at more.')
	b.writeln('')
	b.writeln("| | Points | V's loop | From Python | Overhead | Per point |")
	b.writeln('|---|---|---|---|---|---|')
	for t in pf.timings {
		v, py := median_of(t.v_us), median_of(t.python_us)
		b.writeln('| ${work_label(t.what)} | ${t.points} | ${long_duration(v)} | ${long_duration(py)} | ${signed_duration(py - v)} | ${signed_duration((py - v) / f64(t.points))} |')
	}
	b.writeln('')
	b.writeln('Medians of each call; per point is the overhead over the points.')
	b.writeln('')
	b.writeln("What the overhead is (M10a step 7). A call costs about 0.1 to 0.2 ms more from Python, at any size: the package's arrays and result objects. Past that, the C API adds a few µs a point to properties (the inputs' checks, and copying the results out), and to the PT flash nothing measurable: V's own loop, run inside the library from Python, takes as long as the call. But V's code runs slower inside the library, in a Python process, than in a V program of its own, the flash 1.3 to 1.8 times as long. A heap made large enough that the collector rarely runs doesn't remove it, so it's more than the collections, and its cause isn't known yet: M8's question, with parallel speed.")
	return b.str()
}
