// Package quantity reads numbers with units, such as "230 K", "30 bar" or
// "100 kmol/h", and converts them to SI, for the flowsheet files (M4,
// decision D10). Inside reflux everything is SI; this is the boundary.
//
// Each unit is a row of a fixed table: SI = (x + pre)·scale + post. Celsius
// and Fahrenheit shift before they scale, and gauge pressures add a
// standard atmosphere after. The factors that aren't SI prefixes or time
// come from the exact definitions NIST SP 811 gives: the pound, the inch,
// standard gravity, the IT calorie and Btu, the conventional millimetre of
// mercury. quantity_test.v checks every one against SP 811's table.
module quantity

import math
import strconv
import strings

// Dimension is what kind of quantity a number is.
pub enum Dimension {
	temperature
	temperature_difference // ΔT: no offset, so 1 degC is 1 K
	pressure               // absolute, or gauge only where the unit says so
	pressure_difference    // ΔP: no gauge units
	molar_flow
	mass_flow
	volume_flow
	power // also a heat duty
	molar_enthalpy
	specific_enthalpy
	fraction // a bare number, or a percentage
}

// name is the dimension as a message says it.
pub fn (d Dimension) name() string {
	return match d {
		.temperature { 'temperature' }
		.temperature_difference { 'temperature difference' }
		.pressure { 'pressure' }
		.pressure_difference { 'pressure difference' }
		.molar_flow { 'molar flow' }
		.mass_flow { 'mass flow' }
		.volume_flow { 'volume flow' }
		.power { 'power' }
		.molar_enthalpy { 'molar enthalpy' }
		.specific_enthalpy { 'specific enthalpy' }
		.fraction { 'fraction' }
	}
}

// The exact definitions (NIST SP 811, appendix B).
const pound = 0.45359237 // kg
const inch = 0.0254 // m
const standard_gravity = 9.80665 // m/s²
const pound_force = pound * standard_gravity // N
const psi = pound_force / (inch * inch) // Pa
const foot = 12.0 * inch // m
const atmosphere = 101325.0 // Pa
const btu = 1055.05585262 // J, the IT Btu
const kilocalorie = 4186.8 // J, the IT kilocalorie
// The conventional millimetre of mercury: 13.5951 g/cm³ under standard
// gravity, so 13595.1 kg/m³ · 9.80665 m/s² · 1e-3 m.
const mm_hg = 13595.1 * standard_gravity * 1e-3 // Pa
const torr = atmosphere / 760.0 // Pa
const us_gallon = 231.0 * inch * inch * inch // m³
const horsepower = 550.0 * foot * pound_force // W, 550 ft·lbf/s
const minute = 60.0 // s
const hour = 3600.0 // s
const rankine = 5.0 / 9.0 // K per degree Fahrenheit or Rankine

// Unit is one row: SI = (x + pre)·scale + post.
struct Unit {
	name  string
	dim   Dimension
	scale f64
	pre   f64
	post  f64
}

const table = [
	Unit{'K', .temperature, 1.0, 0.0, 0.0},
	Unit{'degC', .temperature, 1.0, 273.15, 0.0},
	Unit{'°C', .temperature, 1.0, 273.15, 0.0},
	Unit{'degF', .temperature, rankine, 459.67, 0.0},
	Unit{'°F', .temperature, rankine, 459.67, 0.0},
	Unit{'R', .temperature, rankine, 0.0, 0.0},
	Unit{'K', .temperature_difference, 1.0, 0.0, 0.0},
	Unit{'degC', .temperature_difference, 1.0, 0.0, 0.0},
	Unit{'°C', .temperature_difference, 1.0, 0.0, 0.0},
	Unit{'degF', .temperature_difference, rankine, 0.0, 0.0},
	Unit{'°F', .temperature_difference, rankine, 0.0, 0.0},
	Unit{'R', .temperature_difference, rankine, 0.0, 0.0},
	Unit{'Pa', .pressure, 1.0, 0.0, 0.0},
	Unit{'kPa', .pressure, 1e3, 0.0, 0.0},
	Unit{'MPa', .pressure, 1e6, 0.0, 0.0},
	Unit{'bar', .pressure, 1e5, 0.0, 0.0},
	Unit{'mbar', .pressure, 1e2, 0.0, 0.0},
	Unit{'atm', .pressure, atmosphere, 0.0, 0.0},
	Unit{'psi', .pressure, psi, 0.0, 0.0},
	Unit{'psia', .pressure, psi, 0.0, 0.0},
	Unit{'mmHg', .pressure, mm_hg, 0.0, 0.0},
	Unit{'torr', .pressure, torr, 0.0, 0.0},
	Unit{'barg', .pressure, 1e5, 0.0, atmosphere},
	Unit{'psig', .pressure, psi, 0.0, atmosphere},
	Unit{'Pa', .pressure_difference, 1.0, 0.0, 0.0},
	Unit{'kPa', .pressure_difference, 1e3, 0.0, 0.0},
	Unit{'MPa', .pressure_difference, 1e6, 0.0, 0.0},
	Unit{'bar', .pressure_difference, 1e5, 0.0, 0.0},
	Unit{'mbar', .pressure_difference, 1e2, 0.0, 0.0},
	Unit{'atm', .pressure_difference, atmosphere, 0.0, 0.0},
	Unit{'psi', .pressure_difference, psi, 0.0, 0.0},
	Unit{'mmHg', .pressure_difference, mm_hg, 0.0, 0.0},
	Unit{'torr', .pressure_difference, torr, 0.0, 0.0},
	Unit{'mol/s', .molar_flow, 1.0, 0.0, 0.0},
	Unit{'mol/min', .molar_flow, 1.0 / minute, 0.0, 0.0},
	Unit{'mol/h', .molar_flow, 1.0 / hour, 0.0, 0.0},
	Unit{'kmol/s', .molar_flow, 1e3, 0.0, 0.0},
	Unit{'kmol/min', .molar_flow, 1e3 / minute, 0.0, 0.0},
	Unit{'kmol/h', .molar_flow, 1e3 / hour, 0.0, 0.0},
	// A pound-mole is a pound's worth of moles: 453.59237 mol.
	Unit{'lbmol/h', .molar_flow, 1e3 * pound / hour, 0.0, 0.0},
	Unit{'kg/s', .mass_flow, 1.0, 0.0, 0.0},
	Unit{'kg/min', .mass_flow, 1.0 / minute, 0.0, 0.0},
	Unit{'kg/h', .mass_flow, 1.0 / hour, 0.0, 0.0},
	Unit{'g/s', .mass_flow, 1e-3, 0.0, 0.0},
	Unit{'t/h', .mass_flow, 1e3 / hour, 0.0, 0.0},
	Unit{'lb/s', .mass_flow, pound, 0.0, 0.0},
	Unit{'lb/h', .mass_flow, pound / hour, 0.0, 0.0},
	Unit{'m3/s', .volume_flow, 1.0, 0.0, 0.0},
	Unit{'m³/s', .volume_flow, 1.0, 0.0, 0.0},
	Unit{'m3/h', .volume_flow, 1.0 / hour, 0.0, 0.0},
	Unit{'m³/h', .volume_flow, 1.0 / hour, 0.0, 0.0},
	Unit{'L/s', .volume_flow, 1e-3, 0.0, 0.0},
	Unit{'L/min', .volume_flow, 1e-3 / minute, 0.0, 0.0},
	Unit{'L/h', .volume_flow, 1e-3 / hour, 0.0, 0.0},
	Unit{'gal/min', .volume_flow, us_gallon / minute, 0.0, 0.0},
	Unit{'W', .power, 1.0, 0.0, 0.0},
	Unit{'kW', .power, 1e3, 0.0, 0.0},
	Unit{'MW', .power, 1e6, 0.0, 0.0},
	Unit{'kJ/s', .power, 1e3, 0.0, 0.0},
	Unit{'kJ/h', .power, 1e3 / hour, 0.0, 0.0},
	Unit{'MJ/h', .power, 1e6 / hour, 0.0, 0.0},
	Unit{'GJ/h', .power, 1e9 / hour, 0.0, 0.0},
	Unit{'kcal/h', .power, kilocalorie / hour, 0.0, 0.0},
	Unit{'Btu/h', .power, btu / hour, 0.0, 0.0},
	Unit{'hp', .power, horsepower, 0.0, 0.0},
	Unit{'J/mol', .molar_enthalpy, 1.0, 0.0, 0.0},
	Unit{'kJ/mol', .molar_enthalpy, 1e3, 0.0, 0.0},
	Unit{'kJ/kmol', .molar_enthalpy, 1.0, 0.0, 0.0},
	Unit{'MJ/kmol', .molar_enthalpy, 1e3, 0.0, 0.0},
	// A Btu per pound-mole: btu / (1000·pound) J/mol, which is 2.326.
	Unit{'Btu/lbmol', .molar_enthalpy, btu / (1e3 * pound), 0.0, 0.0},
	Unit{'J/kg', .specific_enthalpy, 1.0, 0.0, 0.0},
	Unit{'kJ/kg', .specific_enthalpy, 1e3, 0.0, 0.0},
	Unit{'MJ/kg', .specific_enthalpy, 1e6, 0.0, 0.0},
	Unit{'Btu/lb', .specific_enthalpy, btu / pound, 0.0, 0.0},
	Unit{'%', .fraction, 1e-2, 0.0, 0.0},
]!

// QuantityError is every failure here: the text that couldn't be read, the
// dimension it was read as, why, and any units it may have meant.
pub struct QuantityError {
	Error
pub:
	text        string
	dimension   Dimension
	reason      string
	suggestions []string
}

// msg is the reason, with the text and the dimension it was read as, and
// any suggestions.
pub fn (e QuantityError) msg() string {
	hint := if e.suggestions.len > 0 { ' (did you mean ${e.suggestions.join(', ')}?)' } else { '' }
	return '"${e.text}" as a ${e.dimension.name()}: ${e.reason}${hint}'
}

fn fail(text string, d Dimension, reason string, suggestions []string) IError {
	return QuantityError{
		text:        text
		dimension:   d
		reason:      reason
		suggestions: suggestions
	}
}

// units lists the units a dimension accepts, in the table's order.
pub fn units(d Dimension) []string {
	return table.filter(it.dim == d).map(it.name)
}

// find is the unit of this name and dimension.
fn find(name string, d Dimension) ?Unit {
	for u in table {
		if u.dim == d && u.name == name {
			return u
		}
	}
	return none
}

// parse reads "value unit" as a quantity of dimension d and returns it in
// SI: K, Pa, mol/s, kg/s, m³/s, W, J/mol, J/kg, or a plain fraction. A
// fraction may be a bare number; everything else needs its unit, after a
// space. An absolute temperature or pressure must be above zero.
pub fn parse(text string, d Dimension) !f64 {
	parts := text.fields()
	if parts.len == 0 {
		return fail(text, d, 'it is empty', []string{})
	}
	if parts.len > 2 {
		return fail(text, d, 'expected a number and one unit', []string{})
	}
	value := number(text, d, parts[0])!
	if parts.len == 1 {
		if d == .fraction {
			return value
		}
		return fail(text, d, 'it needs a unit, such as ${units(d)[..math.min(3, units(d).len)].join(' or ')}',
			[]string{})
	}
	u := find(parts[1], d) or { return unknown(text, d, parts[1]) }
	si := (value + u.pre) * u.scale + u.post
	if d == .temperature && !(si > 0.0) {
		return fail(text, d, 'that is ${round6(si)} K, not above absolute zero', []string{})
	}
	if d == .pressure && !(si > 0.0) {
		return fail(text, d, 'that is ${round6(si)} Pa, not above zero', []string{})
	}
	return si
}

// round6 rounds x to six decimals, for messages.
fn round6(x f64) f64 {
	return math.round(x * 1e6) / 1e6
}

// number reads the value, finite. V's strconv also takes "inf", "nan" and
// digit separators ("1_000"); the first two are refused here.
fn number(text string, d Dimension, s string) !f64 {
	v := strconv.atof64(s) or {
		// "30bar": a number and a unit run together.
		mut i := 0
		for i < s.len && in_number(s[i]) {
			i++
		}
		if i > 0 && i < s.len {
			return fail(text, d, 'put a space between the number and its unit: "${s[..i]} ${s[i..]}"',
				[]string{})
		}
		return fail(text, d, '"${s}" is not a number', []string{})
	}
	if !math.is_finite(v) {
		return fail(text, d, '"${s}" is not a finite number', []string{})
	}
	return v
}

// in_number is whether c can be part of a number as strconv reads it.
fn in_number(c u8) bool {
	return c.is_digit() || c == `.` || c == `-` || c == `+` || c == `e` || c == `E` || c == `_`
}

// unknown explains a unit this dimension doesn't have: a unit of another
// dimension, a near miss, or else the units it does have.
fn unknown(text string, d Dimension, name string) IError {
	accepted := units(d)
	for u in table {
		if u.name == name {
			return fail(text, d, '${name} is a unit of ${u.dim.name()}, not ${d.name()}', closest(name,
				accepted))
		}
	}
	near := closest(name, accepted)
	if near.len > 0 {
		return fail(text, d, 'no unit ${name}', near)
	}
	return fail(text, d, 'no unit ${name}; the units of ${d.name()} are ${accepted.join(', ')}',
		[]string{})
}

// closest is the accepted units within two edits of name, ignoring case:
// nearest first, and at most three.
fn closest(name string, accepted []string) []string {
	mut scored := []string{}
	for dist in 0 .. 3 {
		for a in accepted {
			if strings.levenshtein_distance(name.to_lower(), a.to_lower()) == dist && a !in scored {
				scored << a
			}
		}
	}
	return scored[..math.min(3, scored.len)]
}

// to_unit converts an SI value to the named unit of dimension d: the
// inverse of parse, for reports.
pub fn to_unit(si f64, unit string, d Dimension) !f64 {
	u := find(unit, d) or { return unknown('${si} ${unit}', d, unit) }
	return (si - u.post) / u.scale - u.pre
}
