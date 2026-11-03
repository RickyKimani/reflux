// The reports (D11, exit criterion 8): the JSON is the solution, number for
// number, read without the type that wrote it; it balances on its own numbers;
// it's what docs/report-format.md says; and the CSV, the stream table and the
// trace show the same values.
import components
import encoding.csv
import flowsheet { FlowsheetError, Loaded, Report, Solution }
import json2
import math
import os
import strconv
import thermo

const db = components.load_default() or { panic(err) }

fn read_flowsheet(name string) string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', name)) or { panic(err) }
}

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

// solved is a flowsheet's report, and what it came from.
fn solved(text string) !(Loaded, Solution, Report) {
	l := load(text)!
	sol := l.solve()!
	return l, sol, l.report(sol)
}

// The JSON, read as JSON: none of the Report type's own fields or attributes.
type Object = map[string]json2.Any

fn parse(text string) Object {
	return (json2.decode[json2.Any](text) or { panic(err) }).as_map()
}

fn list(a json2.Any) []Object {
	return a.as_array().map(it.as_map())
}

fn numbers(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

fn sum(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

// two_loops is a unit, then a loop, then another, each torn once: the history
// has an entry for each loop, and none for the unit before them.
const two_loops = '
format = 1

[components]
names = ["methane", "n-butane"]

[model]
equation = "PR"

[feed.feed]
T = "300 K"
P = "30 bar"
flow = "10 mol/s"
mole_fractions = { methane = 0.5, n-butane = 0.5 }

[unit.pre]
type = "heater"
in = "feed"
out = "fed"
outlet_T = "290 K"

[unit.first]
type = "mixer"
in = ["fed", "back1"]
out = "m1"

[unit.cool]
type = "heater"
in = "m1"
out = "c1"
outlet_T = "250 K"

[unit.split1]
type = "splitter"
in = "c1"
out = ["back1", "on"]
fractions = [0.4, 0.6]

[unit.second]
type = "mixer"
in = ["on", "back2"]
out = "m2"

[unit.warm]
type = "heater"
in = "m2"
out = "w2"
outlet_T = "260 K"

[unit.split2]
type = "splitter"
in = "w2"
out = ["back2", "product"]
fractions = [0.3, 0.7]
'

// empty_stream is a gas through a drum, which has no liquid to give the
// splitter after it.
const empty_stream = '
format = 1

[components]
names = ["methane", "ethane"]

[model]
equation = "PR"

[feed.in]
T = "300 K"
P = "10 bar"
flow = "1 mol/s"
mole_fractions = { methane = 0.9, ethane = 0.1 }

[unit.drum]
type = "flash_drum"
in = "in"
vapour = "gas"
liquid = "none"
duty = "0 W"

[unit.share]
type = "splitter"
in = "none"
out = ["a", "b"]
fractions = [0.5, 0.5]
'

// interlocking is the reference flowsheet's loop with a second recycle through
// a pump of its own: two cycles in one block.
const interlocking = read_flowsheet('reference.toml').all_before('[unit.mixer]') +
	'
[unit.mixer]
type = "mixer"
in = ["feed", "recycle", "recycle2"]
out = "mixed"

[unit.cooler]
type = "heater"
in = "mixed"
out = "cold"
outlet_T = "230 K"

[unit.drum]
type = "flash_drum"
in = "cold"
vapour = "gas"
liquid = "liquid"
duty = "0 W"

[unit.split]
type = "splitter"
in = "liquid"
out = ["back", "back2", "on"]
fractions = [0.3, 0.2, 0.5]

[unit.pump]
type = "pump"
in = "back"
out = "recycle"
outlet_P = "31 bar"
efficiency = 0.75

[unit.pump2]
type = "pump"
in = "back2"
out = "recycle2"
outlet_P = "32 bar"
efficiency = 0.75
'

struct Case {
	text  string
	types []string // each unit's type, as its file has it
}

// The flowsheets the tests report on.
const cases = {
	'reference':    Case{read_flowsheet('reference.toml'), ['mixer', 'heater', 'flash_drum', 'splitter',
		'pump', 'valve', 'flash_drum']}
	'every unit':   Case{read_flowsheet('every_unit.toml'), ['mixer', 'compressor', 'heater', 'expander',
		'component_separator', 'valve', 'splitter', 'flash_drum', 'pump']}
	'two loops':    Case{two_loops, ['heater', 'mixer', 'heater', 'splitter', 'mixer', 'heater',
		'splitter']}
	'interlocking': Case{interlocking, ['mixer', 'heater', 'flash_drum', 'splitter', 'pump', 'pump']}
	'empty':        Case{empty_stream, ['flash_drum', 'splitter']}
}

fn test_the_json_reads_back_as_the_report() {
	// json2 writes each double as the shortest decimal that reads back as it, so
	// a report comes back from its JSON to the last bit, a failure's too.
	for name, sheet in cases {
		l, _, r := solved(sheet.text)!
		assert json2.decode[Report](r.json())! == r, name
		failed := l.failure(error('a failure from outside'))
		assert json2.decode[Report](failed.json())! == failed
	}
	load('format = 2') or {
		r := flowsheet.load_failure(err)
		assert json2.decode[Report](r.json())! == r
		return
	}
	assert false, 'format 2 loaded'
}

fn test_the_json_is_laid_out_a_value_to_a_line() {
	_, _, r := solved(cases['reference'].text)!
	assert r.json().starts_with('{\n  "format": 1,\n  "reflux": "')
	assert r.json().contains('\n      "molar_flows": [\n        70,\n        10,\n')
	assert r.json().ends_with('\n}\n')
}

fn test_the_json_is_the_solution() {
	for name, sheet in cases {
		types := sheet.types
		l, sol, r := solved(sheet.text)!
		j := parse(r.json())
		assert j['format']!.int() == 1 && j['status']!.str() == 'solved'
		assert j['reflux']!.str() == flowsheet.version()
		assert j['error']! is json2.Null
		// The components, from the database itself.
		comps := j['components']!.as_array().map(it.as_map())
		assert comps.len == l.names.len
		mut mw := []f64{}
		for i, c in comps {
			rec := db.find(l.names[i])!.component()!
			assert c['name']!.str() == l.names[i] && c['compound']!.str() == rec.name
			assert c['molar_mass']!.f64() == rec.mw
			mw << rec.mw
		}
		// Every stream once, each where the flowsheet's ports say.
		mut from, mut to := map[string]string{}, map[string]string{}
		for p in l.sheet.units {
			for s in p.outlets {
				from[s] = p.unit.name
			}
			for s in p.inlets {
				to[s] = p.unit.name
			}
		}
		streams := list(j['streams']!)
		assert streams.len == sol.streams.len, name
		mut seen := []string{}
		for s in streams {
			at := '${name} ${s['name']!.str()}'
			st := sol.streams[s['name']!.str()] or { panic(at) }
			assert s['name']!.str() !in seen, at
			seen << s['name']!.str()
			for key, ends in {
				'from': from
				'to':   to
			} {
				if unit := ends[s['name']!.str()] {
					assert s[key]!.str() == unit, at
				} else {
					assert s[key]! is json2.Null, at
				}
			}
			assert s['T']!.f64() == st.t && s['P']!.f64() == st.p, at
			assert s['H']!.f64() == st.h && s['S']!.f64() == st.s && s['V']!.f64() == st.v, at
			assert s['vapour_fraction']!.f64() == st.beta, at
			flows := numbers(s['molar_flows']!)
			assert flows == st.flows && s['molar_flow']!.f64() == st.total(), at
			mass := numbers(s['mass_flows']!)
			for i, n in flows {
				assert mass[i] == n * mw[i], at
			}
			assert math.abs(s['mass_flow']!.f64() - sum(mass)) <= 1e-15 * sum(mass), at
			assert math.abs(s['volume_flow']!.f64() - st.total() * st.v) <= 1e-15 * st.total() * st.v, at
			x := numbers(s['mole_fractions']!)
			w := numbers(s['mass_fractions']!)
			assert math.abs(sum(x) - 1.0) <= 1e-15 && math.abs(sum(w) - 1.0) <= 1e-15, at
			mut per_mole := 0.0
			for i, xi in x {
				per_mole += xi * mw[i]
				if st.total() > 0.0 {
					assert math.abs(xi - flows[i] / st.total()) <= 1e-15, at
				}
			}
			for i, xi in x {
				assert math.abs(w[i] - xi * mw[i] / per_mole) <= 1e-15, at
			}
			// Its phases make up the stream, the less dense first.
			phases := list(s['phases']!)
			assert phases.len == st.result.phases.len, at
			mut beta, mut h, mut entropy, mut volume := 0.0, 0.0, 0.0, 0.0
			mut z := []f64{len: x.len}
			mut density := []f64{}
			for ph in phases {
				f := ph['fraction']!.f64()
				y := numbers(ph['mole_fractions']!)
				beta += f
				h += f * ph['H']!.f64()
				entropy += f * ph['S']!.f64()
				volume += f * ph['V']!.f64()
				mut m := 0.0
				for i, yi in y {
					z[i] += f * yi
					m += yi * mw[i]
				}
				density << m / ph['V']!.f64()
				assert ph['phase']!.str() in ['vapour', 'liquid'], at
			}
			assert math.abs(beta - 1.0) <= 1e-14 && math.abs(h - st.h) <= 1e-12 * math.abs(st.h), at
			assert math.abs(entropy - st.s) <= 1e-12 * math.abs(st.s), at
			assert math.abs(volume - st.v) <= 1e-12 * st.v, at
			for i in 0 .. x.len {
				assert math.abs(z[i] - x[i]) <= 1e-12, at
			}
			if phases.len == 2 {
				assert density[0] < density[1], at
				assert phases[0]['fraction']!.f64() == st.beta, at
			} else {
				// One phase: its label is what the vapour fraction says.
				label := if st.beta == 1.0 { 'vapour' } else { 'liquid' }
				assert phases[0]['phase']!.str() == label, at
			}
		}
		// Every unit, its energy, and its balance worked out again here.
		units := list(j['units']!)
		assert units.len == l.sheet.units.len
		mut worst_comp, mut worst_heat := 0.0, 0.0
		for i, u in units {
			p := l.sheet.units[i]
			run := sol.units[i]
			at := '${name} ${p.unit.name}'
			assert u['name']!.str() == p.unit.name && u['type']!.str() == types[i], at
			assert u['inlets']!.as_array().map(it.str()) == p.inlets, at
			assert u['outlets']!.as_array().map(it.str()) == p.outlets, at
			energy := list(u['energy']!)
			assert energy.len == run.energy.len, at
			mut q := 0.0
			for k, e in energy {
				assert e['name']!.str() == run.energy[k].name, at
				assert e['rate']!.f64() == run.energy[k].rate, at
				q += run.energy[k].rate
			}
			comp, heat := imbalance(run.inlets, run.outlets, q)
			b := u['balance']!.as_map()
			assert b['component']!.f64() == comp && b['energy']!.f64() == heat, at
			assert comp <= 1e-10 && heat <= 1e-10, at
			worst_comp = math.max(worst_comp, comp)
			worst_heat = math.max(worst_heat, heat)
		}
		// The plan, with each loop's passes.
		blocks := list(j['blocks']!)
		assert blocks.len == sol.plan.blocks.len
		mut loop := 0
		for k, b in blocks {
			plan := sol.plan.blocks[k]
			assert b['units']!.as_array().map(it.str()) == plan.units.map(l.sheet.units[it].unit.name)
			assert b['tears']!.as_array().map(it.str()) == plan.tears
			assert b['cycles']!.int() == plan.cycles
			rows := b['residuals']!.as_array().map(numbers(it))
			if plan.tears.len > 0 {
				assert rows == sol.history[loop].residuals
				loop++
			} else {
				assert rows.len == 0
			}
		}
		assert loop == sol.history.len
		solver := j['solver']!.as_map()
		assert solver['method']!.str() == l.params.method.str()
		assert solver['tol']!.f64() == l.params.tol
		assert solver['max_iterations']!.int() == l.params.max_iterations
		whole := j['balance']!.as_map()
		worst := j['worst_unit']!.as_map()
		assert whole['component']!.f64() == sol.balance.component
		assert whole['energy']!.f64() == sol.balance.energy
		assert worst['component']!.f64() == worst_comp && worst['energy']!.f64() == worst_heat
	}
}

// imbalance is a unit's balances from scratch: its worst component against its
// flow, and its energy against F·RT at its hottest inlet.
fn imbalance(inlets []flowsheet.Stream, outlets []flowsheet.Stream, q f64) (f64, f64) {
	nc := inlets[0].flows.len
	mut comp := []f64{len: nc}
	mut energy, mut flow, mut t := q, 0.0, 0.0
	for s in inlets {
		for i, n in s.flows {
			comp[i] += n
		}
		energy += s.total() * s.h
		flow += s.total()
		t = math.max(t, s.t)
	}
	for s in outlets {
		for i, n in s.flows {
			comp[i] -= n
		}
		energy -= s.total() * s.h
	}
	if flow == 0.0 {
		return 0.0, 0.0
	}
	mut worst := 0.0
	for c in comp {
		worst = math.max(worst, math.abs(c))
	}
	return worst / flow, math.abs(energy) / (flow * thermo.r_gas * t)
}

fn test_the_streams_are_in_the_order_of_the_process() {
	_, _, r := solved(read_flowsheet('reference.toml'))!
	assert r.streams.map(it.name) == ['feed', 'mixed', 'cold', 'gas', 'liquid', 'back', 'on', 'recycle',
		'let down', 'flash gas', 'product']
}

fn test_the_report_balances_on_its_own_numbers() {
	// From the JSON alone, as a reader would check it: what comes in goes out,
	// in moles of each component, in mass and in energy. A unit closes to the
	// tear tolerance on the streams reported, since a tear's consumer was
	// calculated a pass before its value, and the whole closes to what the tears
	// leave open.
	for name, sheet in cases {
		_, _, r := solved(sheet.text)!
		j := parse(r.json())
		tol := j['solver']!.as_map()['tol']!.f64()
		mw := list(j['components']!).map(it['molar_mass']!.f64())
		mut streams := map[string]Object{}
		for s in list(j['streams']!) {
			streams[s['name']!.str()] = s
		}
		mut torn := []string{}
		for b in list(j['blocks']!) {
			torn << b['tears']!.as_array().map(it.str())
		}
		// What a pass can leave open in a tear: tol of its flow in each component,
		// which carries its enthalpy, and tol of RT in its H.
		mut open_flow, mut open_energy := 0.0, 0.0
		for t in torn {
			s := streams[t]
			f := s['molar_flow']!.f64()
			open_flow += tol * f
			open_energy += tol * f * (f64(mw.len) * math.abs(s['H']!.f64()) +
				thermo.r_gas * s['T']!.f64())
		}
		mut comp := []f64{len: mw.len}
		mut mass, mut energy, mut flow, mut t_flow := 0.0, 0.0, 0.0, 0.0
		for _, s in streams {
			f := numbers(s['molar_flows']!)
			fh := s['molar_flow']!.f64() * s['H']!.f64()
			if s['from']! is json2.Null {
				for i, n in f {
					comp[i] += n
				}
				mass += s['mass_flow']!.f64()
				energy += fh
				flow += s['molar_flow']!.f64()
				t_flow += s['molar_flow']!.f64() * s['T']!.f64()
			}
			if s['to']! is json2.Null {
				for i, n in f {
					comp[i] -= n
				}
				mass -= s['mass_flow']!.f64()
				energy -= fh
			}
		}
		for u in list(j['units']!) {
			for e in list(u['energy']!) {
				energy += e['rate']!.f64()
			}
		}
		rounding := 1e-12 * flow
		for c in comp {
			assert math.abs(c) <= open_flow + rounding, name
		}
		mut heaviest := 0.0
		for m in mw {
			heaviest = math.max(heaviest, m)
		}
		assert math.abs(mass) <= (open_flow + rounding) * heaviest * f64(mw.len), name
		rt := thermo.r_gas * t_flow / flow
		assert math.abs(energy) <= open_energy + 1e-10 * flow * rt, name
		for u in list(j['units']!) {
			ins := u['inlets']!.as_array().map(streams[it.str()])
			outs := u['outlets']!.as_array().map(streams[it.str()])
			mut c := []f64{len: mw.len}
			mut q, mut f_in := 0.0, 0.0
			for s in ins {
				for i, n in numbers(s['molar_flows']!) {
					c[i] += n
				}
				q += s['molar_flow']!.f64() * s['H']!.f64()
				f_in += s['molar_flow']!.f64()
			}
			for s in outs {
				for i, n in numbers(s['molar_flows']!) {
					c[i] -= n
				}
				q -= s['molar_flow']!.f64() * s['H']!.f64()
			}
			for e in list(u['energy']!) {
				q += e['rate']!.f64()
			}
			at := '${name} ${u['name']!.str()}'
			for x in c {
				assert math.abs(x) <= open_flow + 1e-10 * f_in, at
			}
			assert math.abs(q) <= open_energy + 1e-10 * f_in * rt, at
		}
	}
}

// documented is every key docs/report-format.md lists, by the object it's in:
// '' for the report itself, `streams[]` for a stream, and so on.
fn documented() map[string][]string {
	text := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md')) or { panic(err) }
	mut out := map[string][]string{}
	mut paths := []string{}
	for line in text.split_into_lines() {
		if line.starts_with('## ') {
			paths = []
		} else if line == '### The report' {
			paths = ['']
			out[''] = []
		} else if line.starts_with('### ') {
			paths = []
			for i, part in line.split('`') {
				if i % 2 == 1 {
					paths << part
					out[part] = []
				}
			}
		} else if line.starts_with('| `') {
			for p in paths {
				out[p] << line.all_after('| `').all_before('`')
			}
		}
	}
	return out
}

fn test_the_json_is_as_documented() {
	// Every key of every object, in order, and nothing more: what M10 and any
	// script that reads a report rely on.
	doc := documented()
	_, _, r := solved(read_flowsheet('every_unit.toml'))!
	j := parse(r.json())
	assert j.keys() == doc['']
	assert list(j['components']!)[0].keys() == doc['components[]']
	stream := list(j['streams']!)[0]
	assert stream.keys() == doc['streams[]']
	assert list(stream['phases']!)[0].keys() == doc['streams[].phases[]']
	unit := list(j['units']!).filter(it['energy']!.as_array().len > 0)[0]
	assert unit.keys() == doc['units[]']
	assert list(unit['energy']!)[0].keys() == doc['units[].energy[]']
	assert unit['balance']!.as_map().keys() == doc['units[].balance']
	assert j['solver']!.as_map().keys() == doc['solver']
	assert list(j['blocks']!)[0].keys() == doc['blocks[]']
	assert j['balance']!.as_map().keys() == doc['balance']
	assert j['worst_unit']!.as_map().keys() == doc['worst_unit']
	load('format = 2') or {
		f := parse(flowsheet.load_failure(err).json())
		assert f.keys() == doc['']
		assert f['error']!.as_map().keys() == doc['error']
		mut all := doc.keys()
		all.sort()
		assert all == ['', 'balance', 'blocks[]', 'components[]', 'error', 'solver', 'streams[]',
			'streams[].phases[]', 'units[]', 'units[].balance', 'units[].energy[]', 'worst_unit']
		return
	}
	assert false, 'format 2 loaded'
}

fn test_a_failure_is_reported_with_how_far_it_got() {
	// Out of passes: the error and the passes, and no streams, units or
	// balances, since nothing unproven is reported (D8).
	text := read_flowsheet('reference.toml').replace('method = "wegstein"', 'method = "direct"\nmax_iterations = 3')
	l := load(text)!
	l.solve() or {
		r := l.failure(err)
		history := (err as FlowsheetError).history
		e := r.error or { panic('no error') }
		assert r.status == 'failed' && e.kind == 'no_convergence' && e.message == err.msg()
		assert e.reason == (err as FlowsheetError).reason && e.unit == '' && e.line == 0
		assert r.streams.len == 0 && r.units.len == 0
		assert r.components.len == 5 && r.solver or { panic('no solver') }.max_iterations == 3
		assert r.blocks.len == 3 && r.blocks[0].tears == ['recycle']
		assert r.blocks[0].residuals.len == 3 && r.blocks[0].residuals == history[0].residuals
		assert r.blocks[1].residuals.len == 0 && r.blocks[2].residuals.len == 0
		j := parse(r.json())
		assert j['status']!.str() == 'failed'
		assert j['balance']! is json2.Null && j['worst_unit']! is json2.Null
		assert j['streams']!.as_array().len == 0
		assert j['error']!.as_map()['cause']! is json2.Null
		return
	}
	assert false, 'it converged in 3 passes'
}

fn test_a_failure_after_a_loop_keeps_the_loop() {
	// A valve can't raise P, which it only finds once the loop before it has
	// converged. The loop's passes are in the failure all the same.
	text := read_flowsheet('reference.toml').replace('outlet_P = "5 bar"', 'outlet_P = "40 bar"')
	l := load(text)!
	l.solve() or {
		e := err as FlowsheetError
		assert e.kind == .bad_specification && e.unit == 'valve', err.msg()
		assert e.history.len == 1 && e.history[0].residuals.len == 4
		r := l.failure(err)
		assert r.blocks[0].residuals == e.history[0].residuals
		assert r.blocks[1].residuals.len == 0
		assert (r.error or { panic('no error') }).unit == 'valve'
		return
	}
	assert false, 'a valve raised P'
}

fn test_a_flash_failure_has_its_cause() {
	// The cooler's outlet below the flash's range: the unit, and the flash's own
	// failure as the cause.
	text := read_flowsheet('reference.toml').replace('outlet_T = "230 K"', 'outlet_T = "1 K"')
	l := load(text)!
	l.solve() or {
		e := err as FlowsheetError
		r := l.failure(err)
		got := r.error or { panic('no error') }
		assert got.kind == 'flash_failed' && got.unit == 'cooler', err.msg()
		cause := e.cause or { panic('no cause: ${err.msg()}') }
		assert got.cause or { '' } == cause.str()
		assert parse(r.json())['error']!.as_map()['cause']!.str() == cause.str()
		return
	}
	assert false, 'it flashed at 1 K'
}

fn test_a_file_that_does_not_load_is_reported() {
	text := read_flowsheet('reference.toml').replace('type = "heater"', 'type = "heatr"')
	load(text) or {
		r := flowsheet.load_failure(err)
		e := r.error or { panic('no error') }
		assert r.status == 'failed' && e.kind == 'bad_file'
		assert e.line == 26 && e.key == 'unit.cooler.type' && e.reason == 'no unit type heatr'
		assert e.suggestions == ['heater'] && e.message == err.msg()
		assert r.components.len == 0 && r.streams.len == 0 && r.blocks.len == 0
		if _ := r.solver {
			assert false, 'a solver for a file that did not load'
		}
		assert parse(r.json())['solver']! is json2.Null
		return
	}
	assert false, 'heatr loaded'
}

fn test_an_error_from_outside_is_reported_as_other() {
	l := load(read_flowsheet('reference.toml'))!
	r := l.failure(error('from a unit of some other module'))
	e := r.error or { panic('no error') }
	assert e.kind == 'other' && e.message == 'from a unit of some other module'
	assert e.reason == e.message && e.unit == '' && e.suggestions.len == 0
}

fn test_an_empty_stream_keeps_its_composition() {
	// A gas through a drum: there's no liquid, so the liquid outlet's flows are
	// 0. It keeps the gas's composition, and with it fractions and a phase. The
	// splitter it goes to has nothing to balance.
	_, _, r := solved(cases['empty'].text)!
	s := r.streams.filter(it.name == 'none')[0]
	assert s.molar_flow == 0.0 && s.mass_flow == 0.0 && s.volume_flow == 0.0
	assert s.molar_flows == [0.0, 0.0] && s.mass_flows == [0.0, 0.0]
	assert math.abs(s.mole_fractions[0] - 0.9) <= 1e-15 && math.abs(sum(s.mass_fractions) - 1.0) <= 1e-15
	assert s.phases.len == 1 && s.phases[0].phase == 'vapour' && s.to or { '' } == 'share'
	assert r.units[1].balance == flowsheet.Balance{}
	assert r.blocks.len == 2 && r.blocks.all(it.tears.len == 0)
	assert r.trace() == 'No loops, so nothing was iterated.\n'
	assert r.table().starts_with('The flowsheet: solved, with no loops, so each unit once.\n')
}

// csv_names is the CSV's columns as docs/report-format.md lists them.
fn csv_names(names []string) []string {
	text := os.read_file(os.join_path(@VMODROOT, 'docs', 'report-format.md')) or { panic(err) }
	block := text.all_after('## The CSV').all_before('then for each').split_into_lines().filter(it.starts_with('    ')).join(' ')
	mut out := block.split(',').map(it.trim_space())
	for key in ['molar_flows', 'mass_flows', 'mole_fractions', 'mass_fractions'] {
		unit := match key {
			'molar_flows' { ' [mol/s]' }
			'mass_flows' { ' [kg/s]' }
			else { '' }
		}
		for c in names {
			out << '${key}.${c}${unit}'
		}
	}
	return out
}

fn test_the_csv_holds_the_streams() {
	// A compound with a comma in its name, and a stream with a quote and a
	// comma, which have to be quoted. Every number reads back as the same double.
	text := read_flowsheet('reference.toml').replace('"n-pentane"', '"1,3-butadiene"').replace('n-pentane = 0.05',
		'"1,3-butadiene" = 0.05').replace('"let down"', r'"let \"down\", gently"')
	l, _, r := solved(text)!
	mut reader := csv.new_reader(r.csv())
	mut rows := [][]string{}
	for {
		rows << reader.read() or { break }
	}
	assert rows[0] == csv_names(l.names)
	assert rows.len == r.streams.len + 1
	for k, s in r.streams {
		row := rows[k + 1]
		assert row.len == rows[0].len
		assert row[0] == s.name && row[1] == s.from or { '' } && row[2] == s.to or { '' }
		mut want := [s.t, s.p, s.vapour_fraction, s.molar_flow, s.mass_flow, s.volume_flow, s.h,
			s.s, s.v]
		want << s.molar_flows
		want << s.mass_flows
		want << s.mole_fractions
		want << s.mass_fractions
		for i, x in want {
			assert strconv.atof64(row[i + 3])! == x, '${s.name} ${rows[0][i + 3]}'
		}
	}
	assert r.streams.any(it.name == 'let "down", gently')
}

// fields is a line of the stream table cut into cells, which two spaces or
// more separate.
fn fields(line string) []string {
	return line.split('  ').map(it.trim_space()).filter(it != '')
}

// shown is whether a cell is x to 5 significant figures.
fn shown(cell string, x f64) bool {
	if x == 0.0 {
		return cell == '0'
	}
	got := strconv.atof64(cell) or { return false }
	e := math.floor(math.log10(math.abs(x)))
	return math.abs(got - x) <= 0.5 * math.pow(10.0, e - 4.0) * (1.0 + 1e-9)
}

fn test_the_stream_table_shows_every_value() {
	_, sol, solved_report := solved(read_flowsheet('reference.toml'))!
	r := Report{
		...solved_report
		file: 'reference.toml'
	}
	t := r.table()
	lines := t.split_into_lines()
	assert lines[0] == "reference.toml: solved by Wegstein, to 1.0e-09 of each tear variable's scale."
	assert lines[1] == '  The loop of mixer, cooler, drum, split, pump: torn at recycle, 4 passes.'
	// The balances, the whole's and the worst unit's, to 2 figures.
	mut said := []f64{}
	for line in lines[2..4] {
		for word in line.split(' ') {
			if word.contains('e-') {
				said << strconv.atof64(word.trim_right(',.'))!
			}
		}
	}
	want := [sol.balance.component, sol.balance.energy, sol.worst.component, sol.worst.energy]
	assert said.len == 4
	for i, x in said {
		e := math.floor(math.log10(want[i]))
		assert math.abs(x - want[i]) <= 0.5 * math.pow(10.0, e - 1.0) * (1.0 + 1e-9), lines[2..4].str()
	}
	// As many streams to a table as fit: another column wouldn't, in 12 or more.
	// Each is right-aligned, 10 wide or its widest, 2 apart, after a blank line.
	heads := lines.filter(it.starts_with('Stream  '))
	assert heads.len == 2 && heads[0].len + 12 > 100
	mut header := 'Stream         '
	for name in r.streams[..7].map(it.name) {
		header += '  ${' '.repeat(10 - name.len)}${name}'
	}
	assert heads[0] == header
	for i, line in lines {
		if line.starts_with('Stream  ') || line.starts_with('Unit  ') {
			assert lines[i - 1] == '', line
		}
	}
	// Each stream's column, as label → cell, from tables that fit in 100 columns.
	mut columns := map[string]map[string]string{}
	mut names := []string{}
	for line in lines {
		f := fields(line)
		if names.len > 0 || line.starts_with('Stream  ') {
			assert line.runes().len <= 100, line
		}
		assert !line.ends_with(' '), 'a trailing space: "${line}"'
		if line.starts_with('Stream  ') {
			names = f[1..]
			for n in names {
				columns[n] = map[string]string{}
			}
		} else if line == '' || line.starts_with('Unit  ') {
			names = []
		} else if names.len > 0 && f.len == names.len + 1 {
			for k, n in names {
				columns[n][f[0]] = f[k + 1]
			}
		}
	}
	assert columns.len == sol.streams.len
	comps := r.components.map(it.name)
	for s in r.streams {
		c := columns[s.name]
		assert c['From'] == s.from or { '-' } && c['To'] == s.to or { '-' }, s.name
		assert shown(c['T (K)'], s.t) && shown(c['P (bar)'], s.p / 1e5), s.name
		assert shown(c['Vapour fraction'], s.vapour_fraction), s.name
		assert shown(c['Flow (mol/s)'], s.molar_flow) && shown(c['Flow (kg/s)'], s.mass_flow), s.name
		assert shown(c['H (kJ/mol)'], s.h / 1e3) && shown(c['S (J/mol/K)'], s.s), s.name
		for i, name in comps {
			assert shown(c[name], s.mole_fractions[i]), '${s.name} ${name}'
		}
	}
	// Each kind of energy in a column of its own, in kW.
	head := lines.filter(it.starts_with('Unit  '))[0]
	assert head.count('Duty (kW)') == 1 && head.count('Work (kW)') == 1, head
	for u in r.units {
		line := lines.filter(it.starts_with('${u.name}  '))[0]
		assert fields(line)[..2] == [u.name, u.kind]
		for e in u.energy {
			column := '${e.name.capitalize()} (kW)'
			assert line.len == (head.index(column) or { -1 }) + column.len, line
			assert shown(fields(line).last(), e.rate / 1e3), line
		}
	}
}

fn test_a_wide_energy_widens_its_column() {
	// gcc's build gives the second drum's adiabatic duty as -9.3132e-13 kW, a
	// rounding where tcc's gives 0: 11 characters, wider than the column's 10.
	_, _, r := solved(read_flowsheet('reference.toml'))!
	mut units := []flowsheet.UnitReport{}
	for u in r.units {
		mut energy := []flowsheet.EnergyStream{}
		for e in u.energy {
			energy << flowsheet.EnergyStream{e.name, if e.name == 'duty' {
				-9.3132e-10
			} else {
				e.rate
			}}
		}
		units << flowsheet.UnitReport{
			...u
			energy: energy
		}
	}
	wide := Report{
		...r
		units: units
	}
	lines := wide.table().split_into_lines()
	head := lines.filter(it.starts_with('Unit  '))[0]
	for u in wide.units {
		line := lines.filter(it.starts_with('${u.name}  '))[0]
		for e in u.energy {
			column := '${e.name.capitalize()} (kW)'
			assert line.len == (head.index(column) or { -1 }) + column.len, line
		}
		if u.energy.any(it.name == 'duty') {
			assert fields(line).last() == '-9.3132e-13', line
		}
	}
}

fn test_the_table_says_how_it_was_solved() {
	for method, words in {
		'direct':  'direct substitution'
		'broyden': 'Broyden'
	} {
		text := read_flowsheet('reference.toml').replace('"wegstein"', '"${method}"')
		_, _, r := solved(text)!
		assert r.table().starts_with('The flowsheet: solved by ${words}, to 1.0e-09'), method
	}
	// A plan's notes, after the loops they're about. Only a search of hundreds
	// of cycles makes one, so this plan is given one.
	l, sol, r := solved(interlocking)!
	assert r.blocks[0].cycles == 2 && r.notes.len == 0
	noted := l.report(flowsheet.Solution{
		...sol
		plan: flowsheet.Plan{
			...sol.plan
			notes: ['a note from the plan']
		}
	})
	assert noted.notes == ['a note from the plan']
	assert noted.table().contains(' passes.\n  a note from the plan\n  The whole closes')
}

fn test_a_table_can_be_exactly_100_wide() {
	// The widest label is 15, and 6 streams of 10 and one of 11 are 100 with
	// their gaps: all 7 fit, as they do in 99.
	text := read_flowsheet('reference.toml').replace('"mixed"', '"mixed fluid"')
	_, _, r := solved(text)!
	heads := r.table().split_into_lines().filter(it.starts_with('Stream  '))
	assert heads[0].len == 100 && fields(heads[0]).len == 8, heads[0]
}

fn test_the_trace_has_every_pass() {
	_, sol, r := solved(read_flowsheet('reference.toml'))!
	lines := r.trace().split_into_lines()
	assert lines[0] == 'The loop of mixer, cooler, drum, split, pump, torn at recycle:'
	assert fields(lines[1]) == ['pass', 'recycle']
	rows := sol.history[0].residuals
	assert lines.len == rows.len + 2
	for k, row in rows {
		f := fields(lines[k + 2])
		assert f[0] == '${k + 1}'
		got := strconv.atof64(f[1])!
		e := math.floor(math.log10(row[0]))
		assert math.abs(got - row[0]) <= 0.5 * math.pow(10.0, e - 2.0) * (1.0 + 1e-9), f[1]
	}
}

fn test_a_loop_that_fails_in_its_first_pass_says_so() {
	// The tear's estimate is a gas, which the pump can't take, so no pass ends.
	text := read_flowsheet('reference.toml').replace('[unit.mixer]', '[estimate.back]\nT = "300 K"\nP = "30 bar"\nmolar_flows = { methane = "30 mol/s" }\n\n[unit.mixer]')
	l := load(text)!
	l.solve() or {
		assert (err as FlowsheetError).unit == 'pump', err.msg()
		r := l.failure(err)
		assert r.blocks[0].tears == ['back'] && r.blocks[0].residuals.len == 0
		assert r.trace() == 'The loop of pump, mixer, cooler, drum, split, torn at back:\n    pass        back\n  (no pass finished)\n'
		return
	}
	assert false, 'a pump took a gas'
}

fn test_the_version_is_the_one_in_v_mod() {
	text := os.read_file(os.join_path(@VMODROOT, 'v.mod'))!
	assert text.contains("version: '${flowsheet.version()}'")
	assert flowsheet.version() != 'unknown'
}
