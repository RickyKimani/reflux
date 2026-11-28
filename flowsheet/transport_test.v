// A phase's transport in a flowsheet (M5b step 7): format 5's
// [[transport.pairs]] and its refusals at their lines, and report format
// 5's values, each the transport model's own bits, with model.transport.
import components
import flowsheet { FlowsheetError, Loaded, Report }
import math
import os
import thermo
import transport
import water

const db = components.load_default() or { panic(err) }

// A liquid methanol–water stream by NRTL, warmed: one liquid phase, each
// compound below its Tc, so every value is the liquid's rules'.
const liquid = 'format = 5

[components]
names = ["methanol", "water"]

[model]
liquid = "NRTL"

[feed.cold]
T = "300 K"
P = "1 bar"
flow = "10 mol/s"
mole_fractions = { methanol = 0.3, water = 0.7 }

[unit.warm]
type = "heater"
in = "cold"
out = "warm"
outlet_T = "320 K"
'

const pair = '
[[transport.pairs]]
a = "water"
b = "methanol"
G = "1.2 - 150 K / T"
'

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

fn solved(text string) !(Loaded, Report) {
	l := load(text)!
	return l, l.report(l.solve()!)
}

fn refused(text string, line int, key string, want string, suggestion string) {
	load(text) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		assert e.key == key, 'key ${e.key}, not ${key}: ${err.msg()}'
		assert e.reason.contains(want), err.msg()
		if suggestion != '' {
			assert suggestion in e.suggestions, err.msg()
		}
		return
	}
	assert false, 'accepted: ${want}'
}

// lined is the 1-based line of a text's first occurrence in text.
fn lined(text string, s string) int {
	assert text.contains(s), s
	return text.all_before(s).count('\n') + 1
}

// Line is G = a − b/T, the file's pair as a transport.Parameter.
struct Line {
	a f64
	b f64
}

fn (l Line) value(t f64) !f64 {
	return l.a - l.b / t
}

fn value_of(v ?f64) f64 {
	return v or { math.nan() }
}

// The file's G is the liquid's Grunberg–Nissan Gᵢⱼ (D6, D15): each liquid's
// viscosity is the transport model's with that pair, bit for bit, and
// ln μ moves from the logarithmic rule's by exactly x₁x₂G(T), to rounding;
// its conductivity and surface tension don't move. The report lists the
// pair as the file wrote it, with its source.
fn test_a_files_pair_is_the_liquids_g() {
	l, rep := solved(liquid + pair)!
	_, plain := solved(liquid)!
	records := ['methanol', 'water'].map(db.find(it)!)
	base := components.mixture(records)!
	with := transport.Mixture{
		...base
		pairs: [transport.Pair{
			i: 1
			j: 0
			g: Line{1.2, 150.0}
		}]
	}
	for k, st in rep.streams {
		assert st.phases.map(it.phase) == ['liquid'], st.name
		ph := st.phases[0]
		want := with.phase(.liquid, st.t, st.p, ph.mole_fractions, ph.v)
		assert value_of(ph.viscosity) == want.viscosity.value, st.name
		assert value_of(ph.thermal_conductivity) == want.conductivity.value, st.name
		assert value_of(ph.surface_tension) == want.surface_tension.value, st.name
		assert ph.transport_notes == want.notes, st.name
		x := ph.mole_fractions
		g := 1.2 - 150.0 / st.t
		p0 := plain.streams[k].phases[0]
		moved := math.log(value_of(ph.viscosity) / value_of(p0.viscosity))
		assert math.abs(moved - x[0] * x[1] * g) <= 1e-14, '${st.name}: ln μ moved ${moved}, not ${x[0] * x[1] * g}'
		assert ph.thermal_conductivity == p0.thermal_conductivity
		assert ph.surface_tension == p0.surface_tension
	}
	m := rep.model or { panic('no model') }
	tm := m.transport or { panic('no transport') }
	assert tm.pairs.len == 1
	assert tm.pairs[0].a == 'water' && tm.pairs[0].b == 'methanol'
	assert tm.pairs[0].parameters.map('${it.name} = ${it.text}') == ['G = 1.2 - 150 K / T']
	assert tm.pairs[0].source == 'the file'
	pm := plain.model or { panic('no model') }
	ptm := pm.transport or { panic('no transport') }
	assert ptm.pairs.map('${it.a} ${it.b} ${it.parameters[0].text} ${it.source}') == [
		'methanol water 0 none listed: 0',
	]
	assert tm.methods.map('${it.property} ${it.phase}') == ['viscosity vapour',
		'thermal_conductivity vapour', 'viscosity liquid', 'thermal_conductivity liquid',
		'surface_tension liquid']
	// Written back, the file keeps its pair, and reports the same.
	again := load(l.write())!
	assert again.report(again.solve()!).json() == rep.json()
}

// A pair's G with no value at a stream's T leaves its liquid without a
// viscosity, saying why; the rest stands.
fn test_a_g_with_no_value_leaves_the_viscosity_none() {
	text := liquid + pair.replace('1.2 - 150 K / T', 'sqrt(1 - T / 310 K)')
	_, rep := solved(text)!
	cold := rep.streams[0].phases[0]
	assert cold.viscosity != none
	warm := rep.streams[1].phases[0]
	assert warm.viscosity == none
	assert warm.thermal_conductivity != none && warm.surface_tension != none
	assert warm.transport_notes.len == 1
	assert warm.transport_notes[0].starts_with('no viscosity: the Grunberg–Nissan G of Water and Methanol at 320.0 K: '), warm.transport_notes[0]
	assert warm.transport_notes[0].contains('square root'), warm.transport_notes[0]
}

// A G large enough to overflow Grunberg and Nissan's exp (0.21 G ≈ 840 > 709)
// leaves no viscosity, not an infinite one.
fn test_a_g_that_overflows_leaves_the_viscosity_none() {
	_, rep := solved(liquid + pair.replace('1.2 - 150 K / T', '4000'))!
	for st in rep.streams {
		ph := st.phases[0]
		assert ph.viscosity == none, st.name
		assert ph.transport_notes == ['no viscosity: the mixing rule gives +inf, not a positive finite number'], '${ph.transport_notes}'
	}
}

// What [transport] can't be is refused at its line.
fn test_transport_is_refused_at_its_line() {
	text := liquid + pair
	at := lined(text, '[[transport.pairs]]')
	// A table is placed at its first key's line, as a section always is.
	refused(text.replace('format = 5', 'format = 4'), at + 1, 'transport', 'transport needs format 5 (M5b), and this file is format 4',
		'')
	refused(text.replace('a = "water"', 'a = "watr"'), at + 1, 'transport.pairs[0].a', 'watr is not one of the components',
		'water')
	refused(text.replace('b = "methanol"', 'b = "water"'), at + 2, 'transport.pairs[0].b', 'a pair is of two components, and this one is water twice',
		'')
	twice := text + pair.replace('a = "water"\nb = "methanol"', 'a = "methanol"\nb = "water"')
	refused(twice, lined(twice, 'a = "methanol"') + 1, 'transport.pairs[1].b', 'methanol and water are a pair already',
		'')
	refused(text.replace('1.2 - 150 K / T', '150 K'), at + 3, 'transport.pairs[0].G', 'G',
		'')
	refused(text.replace('G = "1.2 - 150 K / T"', 'G = 1.2'), at + 3, 'transport.pairs[0].G', 'it should be a string',
		'')
	refused(text.replace('G = "1.2 - 150 K / T"', 'g = "1.2"'), at + 1, 'transport.pairs[0]', 'G is missing',
		'')
	refused(text + 'alpha = "0.3"\n', at + 4, 'transport.pairs[0].alpha', 'alpha', '')
	refused(text.replace('[[transport.pairs]]', '[transport]\nmethod = "Wilke"\n\n[[transport.pairs]]'),
		at + 1, 'transport.method', 'method', '')
}

// A fits file reads format 5 with its [transport] left to the flowsheet.
fn test_a_format_5_files_fits_read() {
	f := flowsheet.load_fits_text(liquid + pair, db: db)!
	assert f.fits.len == 0 && f.names == ['methanol', 'water']
}

// Every phase of the golden flowsheets reports its transport model's bits
// (criterion 4): a cubic's and a γ–φ fluid's by the compounds'
// correlations (components.mixture, built again here), IAPWS-95 water's by
// M5a's functions at its T and density, IF97 water's by IAPWS's industrial
// forms on its state; none without a reason but a vapour's surface tension.
fn test_every_phase_reports_its_models_bits() {
	for name in ['every_unit', 'decanter', 'rankine', 'rankine_if97'] {
		l := flowsheet.load(os.join_path(@VMODROOT, 'testdata', 'flowsheets', '${name}.toml'),
			db: db
		)!
		rep := l.report(l.solve()!)
		records := l.names.map(db.find(it)!)
		mix := components.mixture(records)!
		mut values := 0
		for st in rep.streams {
			for ph in st.phases {
				label := match ph.phase {
					'vapour' { thermo.PhaseLabel.vapor }
					'liquid2' { .liquid2 }
					else { .liquid }
				}
				at := '${name} ${st.name} ${ph.phase}'
				mut want := []f64{}
				if name == 'rankine' {
					rho := water.molar_mass / ph.v
					want << water.viscosity(st.t, rho) or { math.nan() }
					want << water.thermal_conductivity(st.t, rho) or { math.nan() }
					want << if ph.phase == 'vapour' {
						math.nan()
					} else {
						water.surface_tension(st.t) or { math.nan() }
					}
				} else {
					tr := if name == 'rankine_if97' {
						water.TransportIf97{}.phase(label, st.t, st.p, ph.mole_fractions,
							ph.v)
					} else {
						mix.phase(label, st.t, st.p, ph.mole_fractions, ph.v)
					}
					want << tr.viscosity.value
					want << tr.conductivity.value
					want << tr.surface_tension.value
				}
				got := [ph.viscosity, ph.thermal_conductivity, ph.surface_tension]
				for k, g in got {
					if math.is_nan(want[k]) {
						assert g == none, '${at} ${k}: ${g}'
					} else {
						assert value_of(g) == want[k], '${at} ${k}: ${g} against ${want[k]}'
						values++
					}
				}
				// Each none has its reason, but a vapour's surface tension.
				nones := got.filter(it == none).len - if ph.phase == 'vapour' { 1 } else { 0 }
				if nones > 0 {
					assert ph.transport_notes.any(it.starts_with('no ')), '${at}: ${ph.transport_notes}'
				}
			}
		}
		assert values > 0, name
	}
}

// A Loaded made by hand without a transport model reports none, saying so.
fn test_a_loaded_without_transport_reports_none() {
	l := load(liquid)!
	hand := Loaded{
		...l
		transport: none
	}
	rep := hand.report(hand.solve()!)
	ph := rep.streams[0].phases[0]
	assert ph.viscosity == none && ph.thermal_conductivity == none && ph.surface_tension == none
	assert ph.transport_notes == ['no transport model']
}
