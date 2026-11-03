// Units of measure against NIST SP 811 (exit criterion 7), from
// testdata/m4_units.json: every factor and every affine conversion, every
// accepted unit covered, the round trip back from SI, each refusal for its
// own reason, and seeded random input that never panics.
import json2
import math
import os
import quantity { Dimension }

struct Factor {
	unit      string
	dimension string
	factor    f64
	exact     bool
}

struct Point {
	text      string
	dimension string
	si        f64
}

struct Fixture {
	factors []Factor
	points  []Point
}

const dimensions = [Dimension.temperature, .temperature_difference, .pressure, .pressure_difference,
	.molar_flow, .mass_flow, .volume_flow, .power, .molar_enthalpy, .specific_enthalpy, .fraction]

fn load() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm4_units.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

fn dimension(name string) !Dimension {
	for d in dimensions {
		if d.str() == name {
			return d
		}
	}
	return error('no dimension ${name}')
}

fn test_every_factor_matches_nist_sp811() {
	// An exact factor (an SI prefix, the minute or hour, or one SP 811 calls
	// exact) must match to rounding: 4 ulps, since reflux derives it from the
	// definitions in a few operations. The rest SP 811 prints to 7
	// significant figures, so reflux's exact definition must round to them:
	// within half a unit in the 7th figure.
	f := load()!
	for c in f.factors {
		d := dimension(c.dimension)!
		got := quantity.parse('1 ${c.unit}', d)!
		if c.exact {
			assert math.abs(got - c.factor) <= 4.0 * 2.220446049250313e-16 * math.abs(c.factor), '${c.unit} (${c.dimension}): ${got} vs ${c.factor}'
		} else {
			half := 0.5 * math.pow(10.0, math.floor(math.log10(math.abs(c.factor))) - 6.0)
			assert math.abs(got - c.factor) <= half, '${c.unit} (${c.dimension}): ${got} vs ${c.factor}'
		}
	}
}

fn test_every_affine_conversion_matches_nist_sp811() {
	// Celsius, Fahrenheit and Rankine by SP 811's formulas, and gauge
	// pressures above a standard atmosphere. A shift and a scale round a few
	// times, so 1e-12 relative.
	f := load()!
	for c in f.points {
		got := quantity.parse(c.text, dimension(c.dimension)!)!
		assert math.abs(got - c.si) <= 1e-12 * c.si, '${c.text}: ${got} vs ${c.si}'
	}
}

fn test_every_accepted_unit_is_checked() {
	// No unit escapes the two tests above: each is in the fixture's factors,
	// or, for an absolute temperature or a gauge pressure, in its points.
	f := load()!
	for d in dimensions {
		for u in quantity.units(d) {
			in_factors := f.factors.any(it.unit == u && it.dimension == d.str())
			in_points := f.points.any(it.dimension == d.str() && it.text.fields()[1] == u)
			assert in_factors || in_points, '${u} (${d.name()}) is not in testdata/m4_units.json'
		}
	}
}

fn test_to_unit_inverts_parse() {
	// Back from SI to every unit, for values of both signs where the
	// dimension allows them: to 1e-12 of the value.
	for d in dimensions {
		for u in quantity.units(d) {
			values := if d in [.temperature, .pressure] { [2.5, 300.0] } else { [2.5, -40.0, 0.0] }
			for x in values {
				si := quantity.parse('${x} ${u}', d)!
				back := quantity.to_unit(si, u, d)!
				assert math.abs(back - x) <= 1e-12 * math.max(1.0, math.abs(x)), '${x} ${u} (${d.name()}): ${back}'
			}
		}
	}
}

// refused asserts that text is refused as a d, with reason mentioning
// want, and returns the error's suggestions.
fn refused(text string, d Dimension, want string) []string {
	quantity.parse(text, d) or {
		assert err is quantity.QuantityError, err.msg()
		e := err as quantity.QuantityError
		assert e.text == text && e.dimension == d, err.msg()
		assert e.reason.contains(want), '${text}: ${err.msg()}'
		return e.suggestions
	}
	assert false, '"${text}" was accepted as a ${d.name()}'
	return []
}

fn test_each_refusal_has_its_reason() {
	refused('', .pressure, 'empty')
	refused('   ', .pressure, 'empty')
	refused('30', .pressure, 'needs a unit')
	refused('1 2 3', .pressure, 'a number and one unit')
	refused('abc K', .temperature, 'not a number')
	refused('inf bar', .pressure, 'not a finite number')
	refused('nan K', .temperature, 'not a finite number')
	refused('1e400 Pa', .pressure, 'not a number')
	// A number and its unit run together: the fix is in the message.
	refused('30bar', .pressure, '"30 bar"')
	refused('1e5Pa', .pressure, '"1e5 Pa"')
	// Physically impossible.
	refused('-300 degC', .temperature, 'not above absolute zero')
	refused('0 K', .temperature, 'not above absolute zero')
	refused('-1 bar', .pressure, 'not above zero')
	refused('-2 barg', .pressure, 'not above zero')
	// A unit of another dimension is named as one.
	assert 'bar' !in refused('30 K', .pressure, 'is a unit of temperature, not pressure')
	refused('2 barg', .pressure_difference, 'is a unit of pressure, not pressure difference')
	refused('5 kmol/h', .mass_flow, 'is a unit of molar flow, not mass flow')
	// A near miss gets suggestions, nearest first, and capitals are forgiven.
	assert refused('30 Bar', .pressure, 'no unit Bar')[0] == 'bar'
	assert refused('30 bars', .pressure, 'no unit bars')[0] == 'bar'
	assert 'kmol/h' in refused('5 kmole/h', .molar_flow, 'no unit kmole/h')
	assert refused('25 degc', .temperature, 'no unit degc')[0] == 'degC'
	// The masculine ordinal º, a common slip for the degree sign °.
	assert refused('25 ºC', .temperature, 'no unit ºC')[0] == '°C'
	// Only ignoring case puts kPa first: case-sensitively, Pa and kPa are
	// both two edits from KPA, and Pa comes first in the table.
	assert refused('5 KPA', .pressure, 'no unit KPA')[0] == 'kPa'
	// At most three: six pressure units are within two edits of ps.
	assert refused('30 ps', .pressure, 'no unit ps') == ['Pa', 'psi', 'kPa']
	// Nothing near: the units there are.
	assert refused('3 parsecs', .pressure, 'the units of pressure are Pa, kPa').len == 0
}

fn test_differences_and_fractions_take_any_sign() {
	// ΔT and ΔP are differences, so a negative one is fine and has no
	// offset; a fraction may be a bare number or a percentage.
	assert quantity.parse('-5 K', .temperature_difference)! == -5.0
	assert quantity.parse('10 degC', .temperature_difference)! == 10.0
	assert math.abs(quantity.parse('9 degF', .temperature_difference)! - 5.0) <= 1e-15
	assert quantity.parse('-2 bar', .pressure_difference)! == -2e5
	assert quantity.parse('0.75', .fraction)! == 0.75
	assert quantity.parse('75 %', .fraction)! == 0.75
	assert quantity.parse('  30   bar ', .pressure)! == 30e5
}

// Rng is a small seeded generator for the fuzzing: splitmix64.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) below(n int) int {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return int((z ^ (z >> 31)) % u64(n))
}

fn test_random_text_never_panics() {
	// 20,000 seeded strings: pieces of numbers and units, spaces, stray
	// characters and Unicode, and every valid input with one character
	// deleted, inserted or changed. Each is read as every dimension. Any
	// answer is fine, but it's finite, and an absolute T or P is above
	// zero; a panic would end the test.
	pieces := ['1', '0', '-', '+', '.', 'e', '5', '_', ' ', '  ', 'K', 'degC', '°', 'bar', 'barg',
		'k', 'mol', '/', 'h', 's', 'Pa', '%', 'nan', 'inf', 'Ω', '\t', 'x', '9999999999', 'e999']
	mut rng := Rng{20261029}
	mut texts := []string{}
	for _ in 0 .. 15000 {
		mut t := ''
		for _ in 0 .. 1 + rng.below(6) {
			t += pieces[rng.below(pieces.len)]
		}
		texts << t
	}
	f := load()!
	valid := f.points.map(it.text)
	for _ in 0 .. 5000 {
		v := valid[rng.below(valid.len)]
		i := rng.below(v.len)
		c := pieces[rng.below(pieces.len)]
		texts << match rng.below(3) {
			0 { v[..i] + v[i + 1..] }
			1 { v[..i] + c + v[i..] }
			else { v[..i] + c + v[i + 1..] }
		}
	}
	mut accepted := 0
	for t in texts {
		for d in dimensions {
			x := quantity.parse(t, d) or { continue }
			accepted++
			assert math.is_finite(x), '"${t}" as ${d.name()}: ${x}'
			if d in [.temperature, .pressure] {
				assert x > 0.0, '"${t}" as ${d.name()}: ${x}'
			}
		}
	}
	assert accepted > 0
}
