module benchmarks

import json2
import os
import strings

// M5b's speed (step 11, exit criterion 10): a phase's transport in reflux,
// DWSIM 9 and thermo, and IF97 behind reflux's flash beside IAPWS-95's,
// water.if97's functions and DWSIM's Steam Tables package, in one session;
// and what each phase's transport costs a solve's report.
// tools/gen_m5b_bench.py times them all, reflux through tools/bench_m5b.v,
// and writes testdata/m5b_bench.json; tools/benchmarks.vsh adds it to
// docs/benchmarks.md.

pub const m5b_fixture_path = os.join_path(@VMODROOT, 'testdata', 'm5b_bench.json')

// The transport cases, in the report's order, with the phase each is.
pub const m5b_transport_cases = ['benzene-toluene liquid', 'reforming gas', 'methane-propane dense gas']!

// The water cases' programs, by reflux's names.
pub const m5b_water = ['IF97 flash', 'IAPWS-95 flash', 'IF97 functions']!

pub const m5b_kinds = ['pt', 'ph']!

pub const m5b_reports = ['reference', 'rankine_if97']!

// M5bTimed is one program's case: each round's µs (per call for the water
// cases), and its answers: a phase's values (−1 for none), a PT call's
// molar volume or h, a PH call's T, a report's phases with a viscosity and
// with none; or the error its first call raised.
pub struct M5bTimed {
pub:
	name    string
	kind    string
	us      []f64
	answers []f64
	error   string
}

// M5bRun is one of reflux's runs.
pub struct M5bRun {
pub:
	transport []M5bTimed
	water     []M5bTimed
	report    []M5bTimed
}

pub struct M5bSource {
pub:
	generated string
	python    string
	dwsim     string
	thermo    string
	pythonnet string
}

pub struct M5bFixture {
pub:
	source             M5bSource
	machine            Machine
	reflux_build       Build
	reflux_before      M5bRun
	reflux_after       M5bRun
	dwsim              []M5bTimed
	dwsim_steam_tables []M5bTimed
	thermo             []M5bTimed
	seconds            f64
}

// load_m5b reads M5b's benchmark's fixture.
pub fn load_m5b() !M5bFixture {
	text := os.read_file(m5b_fixture_path) or {
		return error('${m5b_fixture_path}: ${err}. Regenerate it on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m5b_bench.py`.')
	}
	return json2.decode[M5bFixture](text)!
}

// m5b_case is a program's case by name and kind, or one with nothing.
pub fn m5b_case(cases []M5bTimed, name string, kind string) M5bTimed {
	for c in cases {
		if c.name == name && c.kind == kind {
			return c
		}
	}
	return M5bTimed{}
}

// m5b_reflux_us is reflux's median µs in a case: the faster of its runs
// before the others' and after.
pub fn m5b_reflux_us(before []M5bTimed, after []M5bTimed, name string, kind string) f64 {
	a := median_of(m5b_case(before, name, kind).us)
	b := median_of(m5b_case(after, name, kind).us)
	return if a < b { a } else { b }
}

// m5b_section is M5b's benchmark's part of docs/benchmarks.md.
pub fn m5b_section(f M5bFixture) string {
	mut b := strings.new_builder(4096)
	b.writeln('')
	b.writeln('## Transport and IF97')
	b.writeln('')
	b.writeln('M5b step 11, exit criterion 10. `tools/gen_m5b_bench.py` writes `testdata/m5b_bench.json` on ${f.machine.cpu}, ${f.machine.os}, on ${f.source.generated}, with DWSIM ${dwsim_version(f.source.dwsim)}, thermo ${f.source.thermo} from Python, and reflux built with ${f.reflux_build.compiler}. Each time is the median of its rounds; reflux ran before the others and again after, and the faster is shown.')
	b.writeln('')
	b.writeln("A phase's transport: reflux's `transport.Mixture.phase`, every property of the phase in one call (a liquid's viscosity, conductivity and surface tension, a vapour's viscosity and conductivity), on Peng–Robinson's phase; DWSIM's Peng–Robinson package's own functions for the same (`AUX_LIQVISCm`, `AUX_CONDTL` and `AUX_SURFTM`, or `AUX_VAPVISCm` and `AUX_CONDTG`) on a stream it calculated, through pythonnet; thermo's mixture properties at their default methods, on its own pure data, each property's `mixture_property`. Their values differ by method, as M5b step 10 found.")
	b.writeln('')
	b.writeln('| Phase | reflux | thermo | DWSIM | DWSIM over reflux |')
	b.writeln('|---|---|---|---|---|')
	for name in m5b_transport_cases {
		r := m5b_reflux_us(f.reflux_before.transport, f.reflux_after.transport, name, 'phase')
		th := m5b_case(f.thermo, name, 'phase')
		d := median_of(m5b_case(f.dwsim, name, 'phase').us)
		tcell := if th.error != '' { '—' } else { duration(median_of(th.us)) }
		b.writeln('| ${name} | ${duration(r)} | ${tcell} | ${duration(d)} | ${times(d, r)} |')
	}
	b.writeln('')
	b.writeln("Water at M5a's three states (300 K and 3 MPa, 500 K and 1 MPa, 700 K and 30 MPa), per call: IF97 behind reflux's flash (`flash.pt` and `flash.ph` on its phase model, M5b's D12), IAPWS-95 behind it, IF97's own functions (`water.if97.state_pt` and `state_ph`), and DWSIM's Steam Tables (IAPWS-IF97) package's flashes. A PH call goes back from the h of each one's own PT call.")
	b.writeln('')
	b.writeln('| | IF97 flash | IAPWS-95 flash | IF97 functions | DWSIM Steam Tables |')
	b.writeln('|---|---|---|---|---|')
	for kind in m5b_kinds {
		mut cells := []string{}
		for name in m5b_water {
			cells << duration(m5b_reflux_us(f.reflux_before.water, f.reflux_after.water,
				name, kind))
		}
		st := m5b_case(f.dwsim_steam_tables, 'Steam Tables', kind)
		cells << if st.error != '' { '—' } else { duration(median_of(st.us)) }
		row := if kind == 'pt' { 'PT' } else { 'PH' }
		b.writeln('| ${row} | ${cells.join(' | ')} |')
	}
	b.writeln('')
	b.writeln("What each phase's transport costs a report (D11: worked out for the report, after the solve, so a solve costs nothing more): `Loaded.report` of a solved flowsheet, with each phase's transport and without it.")
	b.writeln('')
	b.writeln('| Flowsheet | phases | with transport | without | added per phase |')
	b.writeln('|---|---|---|---|---|')
	for name in m5b_reports {
		with := m5b_reflux_us(f.reflux_before.report, f.reflux_after.report, name, 'with transport')
		without := m5b_reflux_us(f.reflux_before.report, f.reflux_after.report, name,
			'without transport')
		counts := m5b_case(f.reflux_before.report, name, 'with transport').answers
		n := counts[0] + counts[1]
		b.writeln('| `${name}.toml` | ${int(n)} | ${duration(with)} | ${duration(without)} | ${duration((with - without) / n)} |')
	}
	return b.str()
}
