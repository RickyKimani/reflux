// Package benchmarks is M3's benchmark harness (decision D10, exit criterion
// 8), roadmap A5's `benchmark` module, renamed because V's standard library
// has a `benchmark` that V's own test runner imports. It times PT, PH and
// bubble-point flashes of a 2-, a 5- and a 10-component mixture in reflux.
// tools/gen_m3_bench.py runs it (through tools/bench_case.v) between
// thermo's and DWSIM's timings of each case, so that all three are timed
// within seconds of each other, and writes testdata/m3_bench.json;
// tools/benchmarks.vsh writes docs/benchmarks.md from that.
module benchmarks

import flash
import json2
import math
import os
import tdep
import thermo
import time

pub const fixture_path = os.join_path(@VMODROOT, 'testdata', 'm3_bench.json')
pub const report_path = os.join_path(@VMODROOT, 'docs', 'benchmarks.md')

// The kinds of flash, in the report's order.
pub const kinds = ['pt2', 'pt1', 'ph', 'bubble']!

// A case is timed for at least this long, and at least min_calls times.
const min_seconds = 1.0
const min_calls = 200

// Component is one compound with DWSIM's constants: Tc in K, Pc in Pa, M in
// kg/mol, and ChemSep's eq. 16 ideal-gas Cp in J/(kmol·K).
pub struct Component {
pub:
	name  string
	cas   string
	tc    f64
	pc    f64
	omega f64
	mw    f64
	cp    []f64
}

pub struct System {
pub:
	name       string
	components []Component
	kij        []f64
	z          []f64
}

// Answer is a program's result at a case's first state: β for PT, T for PH
// and the bubble point, and for PH the H it computed at the state, in J/mol.
pub struct Answer {
pub:
	beta f64
	t    f64
	h    f64
}

// Timing is one program's timing of one case, in microseconds per flash, or
// the error it gave instead.
pub struct Timing {
pub:
	median_us f64
	min_us    f64
	calls     int
	answer    Answer
	error     string
}

// Case is one flash at a case's states: (t[k], p[k]), or p[k] alone for the
// bubble point, with each program's timing of it.
pub struct Case {
pub:
	system        int
	kind          string
	t             []f64
	p             []f64
	reflux        Timing
	thermo        Timing
	dwsim_default Timing
	dwsim_tight   Timing
}

// Build is how reflux's timing program was built.
pub struct Build {
pub:
	compiler string
	v        string
}

pub struct Source {
pub:
	generated string
	python    string
	thermo    string
	dwsim     string
	pythonnet string
}

pub struct Machine {
pub:
	cpu          string
	logical_cpus int
	os           string
}

pub struct Fixture {
pub:
	source       Source
	machine      Machine
	reflux_build Build
	states       int
	systems      []System
	cases        []Case
}

// load reads the fixture.
pub fn load() !Fixture {
	text := os.read_file(fixture_path) or {
		return error('${fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m3_bench.py`.')
	}
	return json2.decode[Fixture](text)!
}

// fluid is Peng–Robinson with the system's constants and kᵢⱼ, and its
// ideal-gas Cp, with no formation data: H is relative to the ideal gas at
// 298.15 K, as DWSIM's is.
pub fn fluid(s System) !thermo.Fluid {
	comps := s.components.map(thermo.Component{
		name:  it.name
		tc:    it.tc
		pc:    it.pc
		omega: it.omega
		mw:    it.mw
	})
	m := thermo.new_cubic(.pr, comps, kij: s.kij)!
	return thermo.new_fluid(m, s.components.map(thermo.IdealGas{
		cp: thermo.Correlated{
			c:     tdep.Correlation{
				eq:   16
				a:    it.cp[0]
				b:    it.cp[1]
				c:    it.cp[2]
				d:    it.cp[3]
				e:    it.cp[4]
				tmax: math.inf(1)
			}
			scale: 1e-3
		}
	}))!
}

// molar_h is a flash result's H for one mole of feed.
pub fn molar_h(fl thermo.Fluid, r flash.FlashResult) f64 {
	mut h := 0.0
	for ph in r.phases {
		h += ph.fraction * fl.props(ph.state).h
	}
	return h
}

// answer_of is what the report compares: β for PT, T otherwise.
fn answer_of(kind string, r flash.FlashResult) Answer {
	return if kind in ['pt2', 'pt1'] {
		Answer{
			beta: r.beta
		}
	} else {
		Answer{
			t: r.t
		}
	}
}

// prepare turns a case into one closure per state, each running reflux's
// flash there. PH is by state: H is computed at the PT state first,
// untimed.
pub fn prepare(fl thermo.Fluid, s System, c Case) ![]fn () !flash.FlashResult {
	mut out := []fn () !flash.FlashResult{}
	z := s.z
	m := fl.model
	for k, p in c.p {
		match c.kind {
			'pt2', 'pt1' {
				t := c.t[k]
				out << fn [m, t, p, z] () !flash.FlashResult {
					return flash.pt(m, t, p, z)
				}
			}
			'ph' {
				h := molar_h(fl, flash.pt(m, c.t[k], p, z)!)
				out << fn [fl, p, h, z] () !flash.FlashResult {
					return flash.ph(fl, p, h, z)
				}
			}
			'bubble' {
				out << fn [m, p, z] () !flash.FlashResult {
					return flash.pvf(m, p, 0.0, z)
				}
			}
			else {
				return error('unknown kind of case: ${c.kind}')
			}
		}
	}
	return out
}

// time_calls times each call on its own, after one untimed pass over them,
// cycling through them for at least min_seconds and min_calls. The answer
// is the first call's.
pub fn time_calls(kind string, calls []fn () !flash.FlashResult) !Timing {
	if calls.len == 0 {
		return error('no states to time')
	}
	first := calls[0]()!
	for f in calls[1..] {
		f()!
	}
	mut ns := []i64{}
	sw := time.new_stopwatch()
	for i := 0; ns.len < min_calls || sw.elapsed().seconds() < min_seconds; i++ {
		f := calls[i % calls.len]
		t0 := time.sys_mono_now()
		f()!
		ns << i64(time.sys_mono_now() - t0)
	}
	return summary(ns, answer_of(kind, first))
}

// compiler is the C compiler this program was built with, and whether with
// -prod.
pub fn compiler() string {
	mut cc := 'an unknown C compiler'
	$if gcc {
		cc = 'gcc'
		version := os.exec(['gcc', '-dumpfullversion'])
		if version.exit_code == 0 {
			cc += ' ' + version.output.trim_space()
		}
	} $else $if clang {
		cc = 'clang'
	} $else $if tinyc {
		cc = 'tcc'
	} $else $if msvc {
		cc = 'MSVC'
	}
	$if prod {
		cc += ' -prod'
	}
	return cc
}

// cpu_name is the processor's name, as Windows or Linux gives it.
pub fn cpu_name() string {
	$if windows {
		r := os.exec(['powershell', '-NoProfile', '-Command', '(Get-CimInstance Win32_Processor).Name'])
		if r.exit_code == 0 && r.output.trim_space() != '' {
			return r.output.trim_space()
		}
	} $else {
		info := os.read_file('/proc/cpuinfo') or { '' }
		for line in info.split_into_lines() {
			if line.starts_with('model name') {
				return line.all_after(':').trim_space()
			}
		}
	}
	return 'an unknown processor'
}

// summary is the median and minimum of ns, in microseconds.
fn summary(ns []i64, answer Answer) Timing {
	mut sorted := ns.clone()
	sorted.sort()
	n := sorted.len
	median := if n % 2 == 1 {
		f64(sorted[n / 2])
	} else {
		0.5 * f64(sorted[n / 2 - 1] + sorted[n / 2])
	}
	return Timing{
		median_us: median / 1e3
		min_us:    f64(sorted[0]) / 1e3
		calls:     n
		answer:    answer
	}
}
