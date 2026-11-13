// Reactors in a flowsheet file (M9 step 9, D11, criteria 7 and 10): every
// type read from format 2 and solved, the whole checked on its elements,
// the report's profile, a long PFR reaching the Gibbs reactor's
// equilibrium, and what the file refuses, at its line.
import flowsheet
import json2
import math
import os
import reactors

const xu_froment = os.read_file(os.join_path(@VMODROOT, 'testdata', 'kinetics', 'xu_froment.toml')) or {
	panic(err)
}

// reactions is Xu and Froment's file from its sites on: its sites and
// reactions, for a flowsheet file.
fn reactions() string {
	return xu_froment[xu_froment.index('[site.Ni]') or { panic('no site') }..]
}

const head = 'format = 2

[components]
names = ["methane", "water", "hydrogen", "carbon monoxide", "carbon dioxide"]

[model]
equation = "PR"

[feed.feed]
T = "900 K"
P = "5 bar"
flow = "10 mol/s"
mole_fractions = { methane = 0.25, water = 0.7, hydrogen = 0.05, "carbon monoxide" = 0, "carbon dioxide" = 0 }
'

// chain is a reactor of each type in series: a conversion reactor, a PFR, a
// CSTR, an equilibrium reactor over the shift alone, and a Gibbs reactor.
const chain = head + '
[unit.convert]
type = "conversion_reactor"
in = "feed"
out = "s1"
conversions = [{ reaction = "reforming", key = "methane", conversion = 0.3 }]

[unit.tube]
type = "pfr"
in = "s1"
out = "s2"
volume = "0.5 m^3"
length = "2 m"
catalyst = { loading = "10 kg/m^3" }
pressure_drop = "0.2 bar"

[unit.tank]
type = "cstr"
in = "s2"
out = "s3"
volume = "0.2 m^3"
catalyst_loading = "10 kg/m^3"
T = "850 K"

[unit.shift]
type = "equilibrium_reactor"
in = "s3"
out = "s4"
reactions = ["shift"]
T = "600 K"

[unit.gibbs]
type = "gibbs_reactor"
in = "s4"
out = "product"
adiabatic = true
'

fn load(text string) !flowsheet.Loaded {
	return flowsheet.load_text(text, makers: reactors.makers())
}

// Every reactor solves in its chain, each closing its elements to 1e-10 of
// its flow and its energy to 1e-10 of F·RT (the solver's checks, on the
// elements where something reacts), and the whole too.
fn test_a_chain_of_every_reactor() {
	l := load(chain + reactions())!
	sol := l.solve()!
	assert l.types == ['conversion_reactor', 'pfr', 'cstr', 'equilibrium_reactor', 'gibbs_reactor']
	for r in sol.units {
		b := r.balance()
		assert r.conserved.len == 3 // C, H and O
		assert b.component <= 1e-10 && b.energy <= 1e-10, '${r.name}: ${b}'
	}
	assert sol.balance.component <= 1e-12 && sol.balance.energy <= 1e-10, '${sol.balance}'
	// The equilibrium reactor ran the shift alone: methane passes through.
	assert sol.streams['s4'].flows[0] == sol.streams['s3'].flows[0]
	// Adiabatic, from 600 K, the Gibbs reactor methanates, which is
	// exothermic, and warms.
	assert sol.streams['product'].t > sol.streams['s4'].t
	assert sol.streams['product'].flows[0] > sol.streams['s4'].flows[0]
}

// The report: format 2, the PFR's profile from its inlet to its outlet,
// its points' keys docs/report-format.md's, and every other unit's null.
fn test_the_profile_in_the_report() {
	l := load(chain + reactions())!
	sol := l.solve()!
	j := json2.decode[json2.Any](json2.encode(l.report(sol)))!.as_map()
	assert j['format']!.int() == 2
	units := j['units']!.as_array()
	for u in units {
		m := u.as_map()
		p := m['profile']!
		if m['type']!.str() != 'pfr' {
			assert p is json2.Null
			continue
		}
		points := p.as_array()
		first := points[0].as_map()
		last := points.last().as_map()
		assert first.keys() == ['volume', 'length', 'T', 'P', 'molar_flows']
		assert first['volume']!.f64() == 0 && last['volume']!.f64() == 0.5
		assert last['length']!.f64() == 2
		assert math.abs(last['P']!.f64() - (5e5 - 0.2e5)) <= 1e-6
		out := sol.streams['s2']
		flows := last['molar_flows']!.as_array().map(it.f64())
		for i, f in flows {
			assert math.abs(f - out.flows[i]) <= 1e-15 * 10, '${i}'
		}
	}
}

// A long PFR reaches equilibrium: at 900 K and 5 bar, with a hundred times
// the catalyst, its outlet is the Gibbs reactor's, reforming and the shift
// being the independent reactions of these species and Xu and Froment's K
// the Gibbs energies'. Their rates run on fugacities here: on partial
// pressures, in Peng–Robinson's gas, they'd stop where Q on pressures is K,
// Σν ln φ short of equilibrium (step 7). Within 1e-10 in each mole fraction
// (measured: 2.8e-14).
fn test_a_long_pfr_reaches_equilibrium() {
	text := head + '
[feed.again]
T = "900 K"
P = "5 bar"
flow = "10 mol/s"
mole_fractions = { methane = 0.25, water = 0.7, hydrogen = 0.05, "carbon monoxide" = 0, "carbon dioxide" = 0 }

[unit.tube]
type = "pfr"
in = "feed"
out = "near"
volume = "10 m^3"
catalyst = { loading = "1000 kg/m^3" }

[unit.gibbs]
type = "gibbs_reactor"
in = "again"
out = "at"
' + reactions().replace('law = "lhhw"', 'law = "lhhw"
basis = "fugacity"')
	l := load(text)!
	sol := l.solve()!
	near := sol.streams['near']
	at := sol.streams['at']
	mut worst := 0.0
	for i in 0 .. 5 {
		d := math.abs(near.flows[i] / near.total() - at.flows[i] / at.total())
		worst = math.max(worst, d)
	}
	println('a long PFR is within ${worst} of equilibrium')
	assert worst <= 1e-10, '${worst}'
}

fn refused(text string) string {
	load(text) or { return err.msg() }
	return 'loaded'
}

// What the file refuses, each at its line: a reactor in format 1, a
// reaction it doesn't have, a reactor that needs reactions in a file
// without, two thermal modes, Ergun without a length, an integrator it
// doesn't know, a conversion's key that isn't a reactant, and a volume
// that isn't one.
fn test_refused() {
	pfr := '
[unit.tube]
type = "pfr"
in = "feed"
out = "out"
volume = "1 m^3"
catalyst = { loading = "10 kg/m^3" }
'
	ok := head + pfr + reactions()
	assert refused(ok) == 'loaded'
	assert refused(ok.replace('format = 2', 'format = 1')).contains('a reactor needs format 2 (M9), and this file is format 1')
	assert refused(head + pfr) == "flowsheet: line 16: unit.tube.type: a reactor of this type needs the file's reactions, and it has none"
	assert refused(ok.replace('volume = "1 m^3"', 'volume = "1 m^3"\nreactions = ["reform"]')).ends_with('unit.tube.reactions: no reaction reform')
	assert refused(ok.replace('volume = "1 m^3"', 'volume = "1 m^3"\nT = "900 K"\nadiabatic = true')).contains('give only one of T, duty and adiabatic, not T and adiabatic')
	assert refused(ok.replace('volume = "1 m^3"', 'volume = "1 m^3"\npressure_drop = "ergun"')).ends_with("unit.tube.pressure_drop: Ergun's pressure drop needs the reactor's length")
	assert refused(ok.replace('volume = "1 m^3"', 'volume = "1 m^3"\nintegrator = "rodas5"')).contains('no integrator rodas5: rodas4 or dopri5 (did you mean rodas4?)')
	assert refused(ok.replace('volume = "1 m^3"', 'volume = "1 m^2"')).ends_with("unit.tube.volume: a reactor's volume must be a volume, but this is an area")
	conv := head + '
[unit.convert]
type = "conversion_reactor"
in = "feed"
out = "out"
conversions = [{ reaction = "reforming", key = "hydrogen", conversion = 0.3 }]
' + reactions()
	assert refused(conv).contains('Hydrogen is not a reactant of reforming')
}

// A reaction without a rate law is left out of a PFR by default, and
// refused where it's named: methanation, here only for an equilibrium.
fn test_a_reaction_without_a_rate() {
	more := reactions() + '
[reaction.methanation]
equation = "CO + 3 H2 = CH4 + H2O"
'
	pfr := head + '
[unit.tube]
type = "pfr"
in = "feed"
out = "out"
volume = "0.1 m^3"
catalyst = { loading = "10 kg/m^3" }
'
	l := load(pfr + more)!
	l.solve()!
	got := refused(pfr.replace('volume = "0.1 m^3"', 'volume = "0.1 m^3"\nreactions = ["methanation"]') +
		more)
	assert got.ends_with('unit.tube: methanation has no rate law'), got
}

// A Gibbs reactor in a recycle loop: the tears converge, and the whole
// closes its elements to what they left open.
fn test_a_reactor_in_a_loop() {
	text := head + '
[unit.mix]
type = "mixer"
in = ["feed", "back"]
out = "mixed"

[unit.gibbs]
type = "gibbs_reactor"
in = "mixed"
out = "hot"
T = "900 K"

[unit.split]
type = "splitter"
in = "hot"
out = ["back", "product"]
fractions = [0.5, 0.5]
' + reactions()
	l := load(text)!
	sol := l.solve()!
	assert sol.history.len == 1
	assert sol.balance.component <= 1e-9 && sol.balance.energy <= 1e-9, '${sol.balance}'
}

// The owner's Tutorial 2 as a file (M9 step 10): the conversion typed as
// the tutorial prints it, with T read in °C, or as DWSIM's percentage;
// solved at 30 bar it warns that the outlet is past equilibrium, in the
// report's JSON and its table, and at 70 bar it doesn't. What the key
// refuses, it refuses at its line.
const tutorial2 = 'format = 2

[components]
names = ["carbon monoxide", "hydrogen", "methanol"]

[model]
equation = "PR"

[feed.feed]
T = "320 degC"
P = "30 bar"
flow = "1500 kg/h"
mole_fractions = { "carbon monoxide" = 0.2, hydrogen = 0.8 }

[definitions]
t = "T in degC"

[unit.reactor]
type = "conversion_reactor"
in = "feed"
out = "product"
conversions = [{ reaction = "synthesis", key = "carbon monoxide", conversion = "0.85 - 0.0018*(T in degC)" }]

[reaction.synthesis]
equation = "CO + 2 H2 = CH3OH"
'

fn test_the_tutorials_conversion_as_printed() {
	l := load(tutorial2)!
	sol := l.solve()!
	feed, out := sol.streams['feed'], sol.streams['product']
	assert math.abs((feed.flows[0] - out.flows[0]) / feed.flows[0] - 0.274) <= 1e-12
	report := l.report(sol)
	assert report.units[0].warnings.len == 1
	assert report.units[0].warnings[0].starts_with('synthesis is past equilibrium at the outlet: Q/K = 4.55')
	j := json2.decode[json2.Any](report.json())!.as_map()
	w := j['units']!.as_array()[0].as_map()['warnings']!.as_array()
	assert w.len == 1 && w[0].str().contains('Q/K')
	assert report.table().contains('Warnings:\n  reactor: synthesis is past equilibrium')
	// The same through a definition, and as a percentage; and no warning
	// at 70 bar.
	for conv in ['"0.85 - 0.0018*t"', '"(85 - 0.18*(T in degC)) [%]"'] {
		l2 := load(tutorial2.replace('"0.85 - 0.0018*(T in degC)"', conv))!
		s2 := l2.solve()!
		x := (s2.streams['feed'].flows[0] - s2.streams['product'].flows[0]) / s2.streams['feed'].flows[0]
		assert math.abs(x - 0.274) <= 1e-12, conv
	}
	l3 := load(tutorial2.replace('P = "30 bar"', 'P = "70 bar"'))!
	s3 := l3.solve()!
	assert l3.report(s3).units[0].warnings.len == 0
	assert !l3.report(s3).table().contains('Warnings')
	// Refusals, at the conversion's line.
	assert refused(tutorial2.replace('(T in degC)"', 'T"')).contains("line 22: unit.reactor.conversions[0].conversion: the two sides of − must have the same units, but the left is a pure number and the right a temperature, K: to use a temperature's number in a unit, read it on that unit's scale, as (T in degC)")
	assert refused(tutorial2.replace('(T in degC)"', '(P in bar)"')).contains('line 22: unit.reactor.conversions[0].conversion: a conversion may depend on T alone, and this uses P')
	// P through a definition is refused too.
	assert refused(tutorial2.replace('t = "T in degC"', 't = "T in degC"
pb = "P in bar"').replace('"0.85 - 0.0018*(T in degC)"',
		'"0.003*pb"')).contains('a conversion may depend on T alone, and this uses P')
	assert refused(tutorial2.replace('"0.85 - 0.0018*(T in degC)"', '"0.0018 1/K * T"')).contains('two values side by side')
	why := refused(tutorial2.replace('"0.85 - 0.0018*(T in degC)"', '"0.5 m"'))
	assert why.contains('a conversion must be a fraction, but this is a length'), why
	// A conversion beyond 1 at the reactor's T fails, naming the T.
	hot := load(tutorial2.replace('"0.85 - 0.0018*(T in degC)"', '"1.6 - 0.0018*(T in degC)"'))!
	hot.solve() or {
		assert err.msg().contains("synthesis's conversion is 1.024 at 593.15 K, not from 0 to 1"), err.msg()
		return
	}
	assert false, 'a conversion above 1'
}
