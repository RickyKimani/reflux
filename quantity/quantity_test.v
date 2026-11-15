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

// M9 step 1: the dimension algebra (D3). Everything below is new; the
// tests above are M4's, unchanged.

// M4's table as it was at v0.2.0, frozen: every M4 unit must read through
// the algebra to exactly these bits, since the golden reports were made
// with them. SI = (x + pre)·scale + post.
struct M4Row {
	name  string
	dim   Dimension
	scale f64
	pre   f64
	post  f64
}

const m4_pound = 0.45359237
const m4_inch = 0.0254
const m4_g0 = 9.80665
const m4_lbf = m4_pound * m4_g0
const m4_psi = m4_lbf / (m4_inch * m4_inch)
const m4_foot = 12.0 * m4_inch
const m4_atm = 101325.0
const m4_btu = 1055.05585262
const m4_kcal = 4186.8
const m4_mmhg = 13595.1 * m4_g0 * 1e-3
const m4_torr = m4_atm / 760.0
const m4_gallon = 231.0 * m4_inch * m4_inch * m4_inch
const m4_hp = 550.0 * m4_foot * m4_lbf
const m4_minute = 60.0
const m4_hour = 3600.0
const m4_rankine = 5.0 / 9.0

const m4_table = [
	M4Row{'K', .temperature, 1.0, 0.0, 0.0},
	M4Row{'degC', .temperature, 1.0, 273.15, 0.0},
	M4Row{'°C', .temperature, 1.0, 273.15, 0.0},
	M4Row{'degF', .temperature, m4_rankine, 459.67, 0.0},
	M4Row{'°F', .temperature, m4_rankine, 459.67, 0.0},
	M4Row{'R', .temperature, m4_rankine, 0.0, 0.0},
	M4Row{'K', .temperature_difference, 1.0, 0.0, 0.0},
	M4Row{'degC', .temperature_difference, 1.0, 0.0, 0.0},
	M4Row{'°C', .temperature_difference, 1.0, 0.0, 0.0},
	M4Row{'degF', .temperature_difference, m4_rankine, 0.0, 0.0},
	M4Row{'°F', .temperature_difference, m4_rankine, 0.0, 0.0},
	M4Row{'R', .temperature_difference, m4_rankine, 0.0, 0.0},
	M4Row{'Pa', .pressure, 1.0, 0.0, 0.0},
	M4Row{'kPa', .pressure, 1e3, 0.0, 0.0},
	M4Row{'MPa', .pressure, 1e6, 0.0, 0.0},
	M4Row{'bar', .pressure, 1e5, 0.0, 0.0},
	M4Row{'mbar', .pressure, 1e2, 0.0, 0.0},
	M4Row{'atm', .pressure, m4_atm, 0.0, 0.0},
	M4Row{'psi', .pressure, m4_psi, 0.0, 0.0},
	M4Row{'psia', .pressure, m4_psi, 0.0, 0.0},
	M4Row{'mmHg', .pressure, m4_mmhg, 0.0, 0.0},
	M4Row{'torr', .pressure, m4_torr, 0.0, 0.0},
	M4Row{'barg', .pressure, 1e5, 0.0, m4_atm},
	M4Row{'psig', .pressure, m4_psi, 0.0, m4_atm},
	M4Row{'Pa', .pressure_difference, 1.0, 0.0, 0.0},
	M4Row{'kPa', .pressure_difference, 1e3, 0.0, 0.0},
	M4Row{'MPa', .pressure_difference, 1e6, 0.0, 0.0},
	M4Row{'bar', .pressure_difference, 1e5, 0.0, 0.0},
	M4Row{'mbar', .pressure_difference, 1e2, 0.0, 0.0},
	M4Row{'atm', .pressure_difference, m4_atm, 0.0, 0.0},
	M4Row{'psi', .pressure_difference, m4_psi, 0.0, 0.0},
	M4Row{'mmHg', .pressure_difference, m4_mmhg, 0.0, 0.0},
	M4Row{'torr', .pressure_difference, m4_torr, 0.0, 0.0},
	M4Row{'mol/s', .molar_flow, 1.0, 0.0, 0.0},
	M4Row{'mol/min', .molar_flow, 1.0 / m4_minute, 0.0, 0.0},
	M4Row{'mol/h', .molar_flow, 1.0 / m4_hour, 0.0, 0.0},
	M4Row{'kmol/s', .molar_flow, 1e3, 0.0, 0.0},
	M4Row{'kmol/min', .molar_flow, 1e3 / m4_minute, 0.0, 0.0},
	M4Row{'kmol/h', .molar_flow, 1e3 / m4_hour, 0.0, 0.0},
	M4Row{'lbmol/h', .molar_flow, 1e3 * m4_pound / m4_hour, 0.0, 0.0},
	M4Row{'kg/s', .mass_flow, 1.0, 0.0, 0.0},
	M4Row{'kg/min', .mass_flow, 1.0 / m4_minute, 0.0, 0.0},
	M4Row{'kg/h', .mass_flow, 1.0 / m4_hour, 0.0, 0.0},
	M4Row{'g/s', .mass_flow, 1e-3, 0.0, 0.0},
	M4Row{'t/h', .mass_flow, 1e3 / m4_hour, 0.0, 0.0},
	M4Row{'lb/s', .mass_flow, m4_pound, 0.0, 0.0},
	M4Row{'lb/h', .mass_flow, m4_pound / m4_hour, 0.0, 0.0},
	M4Row{'m3/s', .volume_flow, 1.0, 0.0, 0.0},
	M4Row{'m³/s', .volume_flow, 1.0, 0.0, 0.0},
	M4Row{'m3/h', .volume_flow, 1.0 / m4_hour, 0.0, 0.0},
	M4Row{'m³/h', .volume_flow, 1.0 / m4_hour, 0.0, 0.0},
	M4Row{'L/s', .volume_flow, 1e-3, 0.0, 0.0},
	M4Row{'L/min', .volume_flow, 1e-3 / m4_minute, 0.0, 0.0},
	M4Row{'L/h', .volume_flow, 1e-3 / m4_hour, 0.0, 0.0},
	M4Row{'gal/min', .volume_flow, m4_gallon / m4_minute, 0.0, 0.0},
	M4Row{'W', .power, 1.0, 0.0, 0.0},
	M4Row{'kW', .power, 1e3, 0.0, 0.0},
	M4Row{'MW', .power, 1e6, 0.0, 0.0},
	M4Row{'kJ/s', .power, 1e3, 0.0, 0.0},
	M4Row{'kJ/h', .power, 1e3 / m4_hour, 0.0, 0.0},
	M4Row{'MJ/h', .power, 1e6 / m4_hour, 0.0, 0.0},
	M4Row{'GJ/h', .power, 1e9 / m4_hour, 0.0, 0.0},
	M4Row{'kcal/h', .power, m4_kcal / m4_hour, 0.0, 0.0},
	M4Row{'Btu/h', .power, m4_btu / m4_hour, 0.0, 0.0},
	M4Row{'hp', .power, m4_hp, 0.0, 0.0},
	M4Row{'J/mol', .molar_enthalpy, 1.0, 0.0, 0.0},
	M4Row{'kJ/mol', .molar_enthalpy, 1e3, 0.0, 0.0},
	M4Row{'kJ/kmol', .molar_enthalpy, 1.0, 0.0, 0.0},
	M4Row{'MJ/kmol', .molar_enthalpy, 1e3, 0.0, 0.0},
	M4Row{'Btu/lbmol', .molar_enthalpy, m4_btu / (1e3 * m4_pound), 0.0, 0.0},
	M4Row{'J/kg', .specific_enthalpy, 1.0, 0.0, 0.0},
	M4Row{'kJ/kg', .specific_enthalpy, 1e3, 0.0, 0.0},
	M4Row{'MJ/kg', .specific_enthalpy, 1e6, 0.0, 0.0},
	M4Row{'Btu/lb', .specific_enthalpy, m4_btu / m4_pound, 0.0, 0.0},
	M4Row{'%', .fraction, 1e-2, 0.0, 0.0},
]

fn test_every_m4_unit_reads_to_the_same_bits() {
	// Through the algebra, each of M4's 75 units has M4's scale to the bit,
	// and each value reads to the bit M4's table gave, for 1, -40.5 (where
	// the dimension takes it) and 1234.5678, and back.
	assert m4_table.len == 75
	for row in m4_table {
		u := quantity.unit(row.name)!
		assert u.scale == row.scale, '${row.name}: scale ${u.scale} against ${row.scale}'
		assert u.dim.eq(row.dim.dim()), '${row.name}: ${u.dim}'
		for x in [1.0, -40.5, 1234.5678] {
			want := (x + row.pre) * row.scale + row.post
			if row.dim in [.temperature, .pressure] && want <= 0.0 {
				continue
			}
			got := quantity.parse('${x} ${row.name}', row.dim)!
			assert got == want, '${x} ${row.name}: ${got} against ${want}'
			back := quantity.to_unit(want, row.name, row.dim)!
			assert back == (want - row.post) / row.scale - row.pre, '${row.name}: back to ${back}'
		}
	}
	// And units(d) still lists exactly M4's units, in M4's order.
	for d in dimensions {
		assert quantity.units(d) == m4_table.filter(it.dim == d).map(it.name), d.name()
	}
}

struct NewFactor {
	unit   string
	factor f64
	exact  bool
	dim    string
}

struct NewFixture {
	factors []NewFactor
}

fn load_new() !NewFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_units.json')
	return json2.decode[NewFixture](os.read_file(path)!)!
}

fn test_every_new_unit_matches_nist_sp811() {
	// As M4's: an exact factor to 4 ulps, the rest to SP 811's 7 printed
	// figures. And each unit's dimension, as quantity writes it.
	f := load_new()!
	for c in f.factors {
		u := quantity.unit(c.unit)!
		if c.exact {
			assert math.abs(u.scale - c.factor) <= 4.0 * 2.220446049250313e-16 * math.abs(c.factor), '${c.unit}: ${u.scale} vs ${c.factor}'
		} else {
			half := 0.5 * math.pow(10.0, math.floor(math.log10(math.abs(c.factor))) - 6.0)
			assert math.abs(u.scale - c.factor) <= half, '${c.unit}: ${u.scale} vs ${c.factor}'
		}
		assert u.dim.str() == c.dim, '${c.unit}: ${u.dim.str()} vs ${c.dim}'
	}
}

fn test_every_symbol_and_prefix_is_checked() {
	// No symbol escapes the test above: each stands alone in the fixture.
	// And each prefix is there on some symbol that takes it.
	f := load_new()!
	names := f.factors.map(it.unit)
	for s in quantity.symbol_names() {
		assert s in names, '${s} is not in testdata/m9_units.json'
	}
	for p in quantity.prefixes() {
		assert names.any(it.starts_with(p) && it[p.len..] in quantity.prefixed()), 'no unit with the prefix ${p} in testdata/m9_units.json'
	}
}

fn ratio_refused(num i64, den i64, want string) {
	if r := quantity.ratio(num, den) {
		assert false, '${num}/${den} gave ${r}'
	} else {
		assert err.msg().contains(want), err.msg()
	}
}

fn test_ratios_are_exact_and_bounded() {
	assert quantity.ratio(2, 4)! == quantity.Ratio{1, 2}
	assert quantity.ratio(3, -6)! == quantity.Ratio{-1, 2}
	assert quantity.ratio(0, 7)! == quantity.Ratio{0, 1}
	assert quantity.ratio(1, 2)!.add(quantity.ratio(1, 3)!)! == quantity.Ratio{5, 6}
	assert quantity.ratio(-3, 4)!.mul(quantity.ratio(2, 3)!)! == quantity.Ratio{-1, 2}
	assert quantity.Ratio{3, 4}.str() == '3/4' && quantity.Ratio{-2, 1}.str() == '-2'
	ratio_refused(1, 0, 'over zero')
	ratio_refused(1_000_001, 1, 'out of range')
	ratio_refused(1, 1_000_001, 'out of range')
	// A sum whose denominator outgrows the limit is refused, not wrapped.
	big := quantity.ratio(1, 999_983)!
	if sum := big.add(quantity.ratio(1, 999_979)!) {
		assert false, 'a sum over ${sum.den}'
	} else {
		assert err.msg().contains('out of range'), err.msg()
	}
}

fn test_the_dimension_algebra() {
	pressure := Dimension.pressure.dim()
	energy := quantity.dim(1, 2, -2, 0, 0, 0)
	volume := quantity.dim(0, 3, 0, 0, 0, 0)
	// A pressure is an energy per volume; a power an energy per time.
	assert energy.div(volume)!.eq(pressure)
	assert energy.div(quantity.dim(0, 0, 1, 0, 0, 0))!.eq(Dimension.power.dim())
	// A square root, squared, is exact.
	half := pressure.pow(quantity.ratio(1, 2)!)!
	assert half.mul(half)!.eq(pressure)
	assert half.pow(quantity.ratio(2, 1)!)!.eq(pressure)
	assert !half.eq(pressure)
	// Every base counts: a coulomb is an ampere-second, not a second.
	assert quantity.unit('C')!.dim.eq(quantity.unit('A*s')!.dim)
	assert !quantity.unit('C')!.dim.eq(quantity.unit('s')!.dim)
	assert !quantity.dim(0, 0, 0, 0, 1, 0).eq(quantity.dim(0, 0, 0, 0, 1, 1))
	// x/x is a pure number.
	assert pressure.div(pressure)!.is_dimensionless()
	assert !pressure.is_dimensionless()
	assert quantity.Dim{}.str() == '1'
	// Names: M4's dimensions, then the common ones, then SI base units. A
	// difference shares its absolute quantity's dimension, and that is the
	// name a message gives.
	assert Dimension.temperature_difference.dim().name() == 'temperature'
	assert pressure.name() == 'pressure'
	assert energy.name() == 'energy'
	assert quantity.dim(1, -1, -1, 0, 0, 0).name() == 'viscosity'
	assert quantity.dim(-1, 0, -1, 0, 1, 0).name() == 'mol/(kg·s)'
	assert half.name() == 'kg^(1/2)/(m^(1/2)·s)'
	// Exponents in Base's order.
	e := quantity.dim(1, -1, -2, 0, 0, 0).exponents()
	assert e[int(quantity.Base.mass)] == quantity.Ratio{1, 1}
	assert e[int(quantity.Base.time)] == quantity.Ratio{-2, 1}
}

// same asserts that two units read to the same dimension and, within 4
// ulps, the same scale.
fn same(a string, b string) {
	x := quantity.unit(a) or { panic('${a}: ${err.msg()}') }
	y := quantity.unit(b) or { panic('${b}: ${err.msg()}') }
	assert x.dim.eq(y.dim), '${a} is ${x.dim}, ${b} is ${y.dim}'
	assert math.abs(x.scale - y.scale) <= 4.0 * 2.220446049250313e-16 * x.scale, '${a}: ${x.scale}, ${b}: ${y.scale}'
}

fn test_compound_units() {
	// The same unit written several ways.
	same('J/mol/K', 'J/(mol*K)')
	same('J/(mol·K)', 'J/(mol⋅K)')
	same('m3', 'm^3')
	same('m³', 'm^(3/1)')
	same('s^-1', '1/s')
	same('s⁻¹', '1/s')
	same('s⁻²', 's^-2')
	same('bar^0.5', 'bar^(1/2)')
	same('bar^-0.5', '1/bar^0.5')
	same('bar^0.5*bar^0.5', 'bar')
	same('(bar^(1/3))^3', 'bar')
	same('cP', 'mPa*s')
	same('kmol/(kg*h)', 'kmol*kg^-1*h^-1')
	same('kmol/(kg*h)', 'mol/(g*h)')
	same('kJ/kmol', 'J/mol')
	same('L', 'dm3')
	same('mol/L', 'kmol/m3')
	same('(kmol/h)^2', 'kmol^2/h^2')
	same('Btu/(lb*degF)', 'kcal/(kg*K)')
	same('1/(1/s)', 's')
	// A Dim written as quantity writes it reads back as SI base units.
	for text in ['kmol*bar^0.5/(kg*h)', 'J/(mol*K)', 'C/mol', 'bar^-2', 'cP', '%', 'lb/ft^3'] {
		u := quantity.unit(text)!
		back := quantity.unit(u.dim.str())!
		assert back.dim.eq(u.dim) && back.scale == 1.0, '${text}: ${u.dim.str()} reads as ${back.dim}, ${back.scale}'
	}
	// A whole power is a product, to the bit, the same on every platform:
	// cm⁴ is four centimetres, where this machine's math.pow gives 1e-8.
	assert quantity.unit('cm3')!.scale == 1e-2 * 1e-2 * 1e-2
	assert quantity.unit('cm4')!.scale == 1e-2 * 1e-2 * 1e-2 * 1e-2
	assert quantity.unit('km^-2')!.scale == 1.0 / (1e3 * 1e3)
	// A fraction's sign, on either part.
	same('bar^(-1/2)', 'bar^-0.5')
	same('bar^(1/-2)', '1/bar^(1/2)')
	assert !quantity.unit('bar^(-1/2)')!.dim.eq(quantity.unit('bar^(1/2)')!.dim)
	// Celsius and gauge pressures keep their offsets only standing alone.
	assert quantity.unit('°C')!.pre == 273.15
	assert quantity.unit('kJ/(kg*°C)')!.pre == 0.0
	assert quantity.unit('barg')!.post == 101325.0 && quantity.unit('barg')!.gauge
	assert !quantity.unit('bar')!.gauge
}

// bad asserts that a unit is refused at a column, with reason mentioning
// want, and returns the error.
fn bad(text string, column int, want string) quantity.UnitError {
	quantity.unit(text) or {
		assert err is quantity.UnitError, err.msg()
		e := err as quantity.UnitError
		assert e.text == text, err.msg()
		assert e.reason.contains(want), '${text}: ${err.msg()}'
		assert e.column == column, '${text}: column ${e.column}, not ${column}: ${err.msg()}'
		return e
	}
	assert false, '"${text}" was read as a unit'
	return quantity.UnitError{}
}

fn test_each_unit_refusal_has_its_reason() {
	bad('', 1, 'empty')
	// kmol/kg*h might be kmol/(kg·h) or kmol·h/kg.
	bad('kmol/kg*h', 8, 'ambiguous')
	bad('J/mol·K', 6, 'ambiguous')
	// DWSIM writes kmol/[kg.h].
	bad('kmol/[kg.h]', 6, 'with ( and ), not [ and ]')
	bad('kg.h', 3, 'put * or · between units, not .')
	bad('kmol/(kg*h', 6, 'a ( with no ) to close it')
	bad('kmol/h)', 7, 'a ) with no ( before it')
	bad('kmol/(kg.h)', 9, 'put * or ·')
	bad('kmol/', 6, 'missing at the end')
	bad('*h', 1, 'missing before *')
	bad('kmol//h', 6, 'missing before /')
	bad('^2', 1, 'nothing to raise')
	bad('m^', 2, 'a ^ with no power after it')
	bad('m^1.', 2, 'digits after its decimal point')
	bad('m^12345678', 2, 'too many digits')
	bad('m^1000001', 2, 'out of range')
	bad('m^0.1234567', 2, 'too many digits')
	bad('m3^2', 3, 'two exponents')
	bad('m²³⁴⁵⁶⁷⁸⁹', 2, 'too many digits')
	bad('m⁻', 2, 'a ⁻ with no power after it')
	bad('m^(1/0)', 2, 'over zero')
	bad('m^(1/2', 2, 'a fraction')
	bad('m^(a/2)', 4, 'a fraction of whole numbers')
	bad('2/s', 1, 'the 1 of 1/s')
	// Digits are a power only right after a symbol, as in m3.
	bad('(m)3', 4, 'a number can only be an exponent')
	bad('kmol/2', 6, 'the 1 of 1/s')
	// km to the 400th is 1e1200 m^400: no double holds its scale.
	bad('km^400', 1, 'out of range')
	bad('J/barg', 3, 'gauge pressure, which can only stand alone')
	bad('('.repeat(17) + 'm' + ')'.repeat(17), 17, 'too many parentheses')
	// An unknown name, with suggestions among every symbol and prefix.
	e := bad('kmole/h', 1, 'no unit kmole')
	assert e.unknown == 'kmole' && 'kmol' in e.suggestions, e.msg()
	assert 'Pa' in bad('kmol/(kg*h*pa)', 12, 'no unit pa').suggestions
	// A symbol that takes no prefix, and no prefix of pico or smaller.
	bad('katm', 1, 'no unit katm')
	bad('kkmol', 1, 'no unit kkmol')
	bad('ps', 1, 'no unit ps')
	bad('ºC', 1, 'no unit ºC') // the masculine ordinal, not the degree sign
	// A whole unit is under 17 nested parentheses.
	same('('.repeat(16) + 'm' + ')'.repeat(16), 'm')
}

fn test_parse_reads_any_unit_of_the_dimension() {
	assert quantity.parse('2 lbmol/min', .molar_flow)! == 2.0 * ((1e3 * m4_pound) / 60.0)
	assert quantity.parse('1 mol/ms', .molar_flow)! == 1e3
	assert quantity.parse('5 m^3/h', .volume_flow)! == quantity.parse('5 m3/h', .volume_flow)!
	assert quantity.parse('3 kJ/s', .power)! == 3e3
	assert quantity.parse('1 kcal/s', .power)! == m4_kcal
	assert quantity.parse('100 mol/mol', .fraction)! == 100.0
	assert quantity.parse('10 mbar', .pressure_difference)! == 1e3
	assert math.abs(quantity.to_unit(1.0, 'kmol/min', .molar_flow)! - 0.06) <= 1e-17
	// A unit of another dimension is named by what it is.
	refused('100 kmol/(kg*h)', .molar_flow, 'kmol/(kg*h) is a unit of mol/(kg·s), not molar flow')
	refused('1 J/(mol*K)', .molar_enthalpy, 'is a unit of molar entropy, not molar enthalpy')
	refused('1 m/s', .volume_flow, 'is a unit of velocity, not volume flow')
	refused('1 J/barg', .pressure, 'gauge pressure, which can only stand alone')
	// A unit that doesn't parse says why and where.
	refused('5 kmol/kg*h', .molar_flow, 'kmol/kg*h: a * after a / is ambiguous')
	refused('5 kmol/kg*h', .molar_flow, 'at character 8')
	refused('5 kmol/[kg.h]', .molar_flow, 'not [ and ]')
	// A finite number and a finite scale can still overflow.
	refused('1e300 GPa', .pressure, 'beyond the range of a number')
	refused('1e300 GPa', .pressure_difference, 'beyond the range of a number')
}

fn test_random_units_never_panic() {
	// 30,000 seeded strings of unit syntax: symbols, prefixes, operators,
	// exponents of every form, parentheses, superscripts and stray Unicode.
	// Any answer is fine, but an accepted unit has a finite scale above
	// zero, and its dimension, written out, reads back as itself.
	pieces := ['m', 'k', 'mol', 'g', 's', 'h', 'K', 'Pa', 'bar', 'barg', '°C', 'J', 'L', 'lb',
		'µ', 'μ', 'c', 'P', '%', '*', '·', '/', '(', ')', '^', '-', '.', '0', '1', '2', '3',
		'9', '²', '³', '⁻', '¹', '[', 'Ω', ' ', 'e', '999999', '(1/2)', '^-1', '^0.5']
	mut rng := Rng{20261109}
	mut accepted := 0
	for _ in 0 .. 30000 {
		mut t := ''
		for _ in 0 .. 1 + rng.below(8) {
			t += pieces[rng.below(pieces.len)]
		}
		u := quantity.unit(t) or {
			assert err is quantity.UnitError, '"${t}": ${err.msg()}'
			continue
		}
		accepted++
		assert math.is_finite(u.scale) && u.scale > 0.0, '"${t}": ${u.scale}'
		back := quantity.unit(u.dim.str()) or { panic('"${t}": its dimension ${u.dim} does not read: ${err.msg()}') }
		assert back.dim.eq(u.dim) && back.scale == 1.0, '"${t}": ${u.dim}'
	}
	assert accepted > 1000
}
