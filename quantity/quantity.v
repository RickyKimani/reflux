// Package quantity reads numbers with units, such as "230 K", "30 bar" or
// "100 kmol/h", and converts them to SI, for the flowsheet files (M4,
// decision D10). Inside reflux everything is SI; this is the boundary.
//
// Since M9 (its D3) it has a dimension algebra: a dimension is a vector of
// exponents of the SI base units, and a unit is a scale on one. So a unit
// can be compound, such as kmol/(kg*h), J/(mol*K) or bar^-0.5, and a unit
// of the wrong dimension is refused with what it is a unit of. M4's named
// dimensions (pressure, molar flow, ...) are dimensions of the algebra, and
// read every unit they read before, to the same bits.
//
// Each unit symbol is a row of a fixed table: SI = (x + pre)·scale + post.
// Celsius and Fahrenheit shift before they scale and gauge pressures add a
// standard atmosphere after, but only standing alone: in a compound unit,
// °C is a kelvin of difference and a gauge pressure is refused. The factors
// that aren't SI prefixes or time come from the exact definitions NIST
// SP 811 gives. quantity_test.v checks every one against SP 811's table.
module quantity

import math
import strconv
import strings

// Base is an SI base quantity, in the order a dimension's exponents are
// written.
pub enum Base {
	mass        // kg
	length      // m
	time        // s
	temperature // K
	amount      // mol
	current     // A, for the coulomb and the Faraday constant
}

// The base units' symbols, in Base's order.
const base_symbols = ['kg', 'm', 's', 'K', 'mol', 'A']!

// Ratio is an exact rational number num/den in lowest terms, den > 0. A
// dimension's exponents are ratios, so that bar^0.5 squared is bar again,
// exactly.
pub struct Ratio {
pub:
	num i64
	den i64 = 1
}

// ratio_limit bounds a ratio's numerator and denominator, far beyond any
// unit's, so that no sum or product of two in range overflows an i64.
const ratio_limit = i64(1_000_000)

// ratio is num/den in lowest terms, refused over zero or out of range.
pub fn ratio(num i64, den i64) !Ratio {
	if den == 0 {
		return error('a ratio over zero')
	}
	mut n, mut d := num, den
	if d < 0 {
		n, d = -n, -d
	}
	g := gcd(if n < 0 { -n } else { n }, d)
	n, d = n / g, d / g
	if n > ratio_limit || n < -ratio_limit || d > ratio_limit {
		return error('an exponent of ${Ratio{n, d}} is out of range')
	}
	return Ratio{n, d}
}

fn gcd(a i64, b i64) i64 {
	mut x, mut y := a, b
	for y != 0 {
		x, y = y, x % y
	}
	return x
}

// add is a + b, exact, or refused out of range.
pub fn (a Ratio) add(b Ratio) !Ratio {
	return ratio(a.num * b.den + b.num * a.den, a.den * b.den)
}

// mul is a·b, exact, or refused out of range.
pub fn (a Ratio) mul(b Ratio) !Ratio {
	return ratio(a.num * b.num, a.den * b.den)
}

// neg is −a.
pub fn (a Ratio) neg() Ratio {
	return Ratio{-a.num, a.den}
}

// value is the ratio as the nearest double.
pub fn (a Ratio) value() f64 {
	return f64(a.num) / f64(a.den)
}

// str writes a whole ratio as an integer, and any other as p/q.
pub fn (a Ratio) str() string {
	return if a.den == 1 { '${a.num}' } else { '${a.num}/${a.den}' }
}

// Dim is a dimension: the exponent of each SI base unit. It has named
// fields, not an array, because V 0.5.2 leaves a fixed array of structs
// zeroed, its field defaults unapplied, inside another struct's zero value
// (docs/v-for-go.md, §7).
pub struct Dim {
pub:
	mass        Ratio
	length      Ratio
	time        Ratio
	temperature Ratio
	amount      Ratio
	current     Ratio
}

// dim is the dimension with these whole exponents of kg, m, s, K, mol and A.
pub fn dim(mass int, length int, time int, temperature int, amount int, current int) Dim {
	return Dim{Ratio{i64(mass), 1}, Ratio{i64(length), 1}, Ratio{i64(time), 1}, Ratio{i64(temperature), 1}, Ratio{i64(amount), 1}, Ratio{i64(current), 1}}
}

// exponents is the dimension's exponents in Base's order.
pub fn (d Dim) exponents() [6]Ratio {
	return [d.mass, d.length, d.time, d.temperature, d.amount, d.current]!
}

fn from_exponents(e [6]Ratio) Dim {
	return Dim{e[0], e[1], e[2], e[3], e[4], e[5]}
}

// mul is the dimension of a product.
pub fn (a Dim) mul(b Dim) !Dim {
	x, y := a.exponents(), b.exponents()
	mut e := [6]Ratio{}
	for i in 0 .. 6 {
		e[i] = x[i].add(y[i])!
	}
	return from_exponents(e)
}

// div is the dimension of a quotient.
pub fn (a Dim) div(b Dim) !Dim {
	x, y := a.exponents(), b.exponents()
	mut e := [6]Ratio{}
	for i in 0 .. 6 {
		e[i] = x[i].add(y[i].neg())!
	}
	return from_exponents(e)
}

// pow is the dimension raised to p.
pub fn (a Dim) pow(p Ratio) !Dim {
	x := a.exponents()
	mut e := [6]Ratio{}
	for i in 0 .. 6 {
		e[i] = x[i].mul(p)!
	}
	return from_exponents(e)
}

// eq is whether two dimensions are the same.
pub fn (a Dim) eq(b Dim) bool {
	x, y := a.exponents(), b.exponents()
	for i in 0 .. 6 {
		if x[i].num != y[i].num || x[i].den != y[i].den {
			return false
		}
	}
	return true
}

// is_dimensionless is whether the dimension is a pure number's.
pub fn (d Dim) is_dimensionless() bool {
	return d.eq(Dim{})
}

// str writes the dimension in SI base units, the positive powers over the
// negative ones: kg/(m·s²) for a pressure, mol/(kg·s) for a rate per mass
// of catalyst, 1 for a pure number. A fractional power is written ^(p/q).
pub fn (d Dim) str() string {
	mut top := []string{}
	mut bottom := []string{}
	for i, e in d.exponents() {
		if e.num > 0 {
			top << power_of(base_symbols[i], e)
		} else if e.num < 0 {
			bottom << power_of(base_symbols[i], e.neg())
		}
	}
	num := if top.len == 0 { '1' } else { top.join('·') }
	if bottom.len == 0 {
		return num
	}
	den := if bottom.len == 1 { bottom[0] } else { '(' + bottom.join('·') + ')' }
	return '${num}/${den}'
}

fn power_of(symbol string, e Ratio) string {
	if e.den != 1 {
		return '${symbol}^(${e})'
	}
	if e.num == 1 {
		return symbol
	}
	return symbol + superscript(e.num)
}

const superscripts = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']!

// superscript writes a power above 1 as superscript digits. A negative power
// never comes here: str writes it below the line, as a positive one.
fn superscript(n i64) string {
	mut out := ''
	for c in n.str().bytes() {
		out += superscripts[int(c - u8(`0`))]
	}
	return out
}

// Known is a dimension's name, for messages, beyond the named dimensions.
struct Known {
	name string
	dim  Dim
}

const known = [
	Known{'length', dim(0, 1, 0, 0, 0, 0)},
	Known{'area', dim(0, 2, 0, 0, 0, 0)},
	Known{'volume', dim(0, 3, 0, 0, 0, 0)},
	Known{'mass', dim(1, 0, 0, 0, 0, 0)},
	Known{'time', dim(0, 0, 1, 0, 0, 0)},
	Known{'amount of substance', dim(0, 0, 0, 0, 1, 0)},
	Known{'electric current', dim(0, 0, 0, 0, 0, 1)},
	Known{'electric charge', dim(0, 0, 1, 0, 0, 1)},
	Known{'force', dim(1, 1, -2, 0, 0, 0)},
	Known{'energy', dim(1, 2, -2, 0, 0, 0)},
	Known{'density', dim(1, -3, 0, 0, 0, 0)},
	Known{'viscosity', dim(1, -1, -1, 0, 0, 0)},
	Known{'concentration', dim(0, -3, 0, 0, 1, 0)},
	Known{'molar volume', dim(0, 3, 0, 0, -1, 0)},
	Known{'molar entropy', dim(1, 2, -2, -1, -1, 0)},
	Known{'velocity', dim(0, 1, -1, 0, 0, 0)},
]!

// name is the dimension as a message says it: a named dimension's name,
// another common one's, or else its SI base units.
pub fn (d Dim) name() string {
	for n in [Dimension.temperature, .pressure, .molar_flow, .mass_flow, .volume_flow, .power,
		.molar_enthalpy, .specific_enthalpy, .fraction] {
		if d.eq(n.dim()) {
			return n.name()
		}
	}
	for k in known {
		if d.eq(k.dim) {
			return k.name
		}
	}
	return d.str()
}

// Dimension is a named kind of quantity, as M4's flowsheet files ask for
// one. Each is a Dim; a difference has its absolute quantity's.
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

// dim is the named dimension's exponents.
pub fn (d Dimension) dim() Dim {
	return match d {
		.temperature, .temperature_difference { dim(0, 0, 0, 1, 0, 0) }
		.pressure, .pressure_difference { dim(1, -1, -2, 0, 0, 0) }
		.molar_flow { dim(0, 0, -1, 0, 1, 0) }
		.mass_flow { dim(1, 0, -1, 0, 0, 0) }
		.volume_flow { dim(0, 3, -1, 0, 0, 0) }
		.power { dim(1, 2, -3, 0, 0, 0) }
		.molar_enthalpy { dim(1, 2, -2, 0, -1, 0) }
		.specific_enthalpy { dim(0, 2, -2, 0, 0, 0) }
		.fraction { Dim{} }
	}
}

// absolute is whether an offset unit (°C, °F, barg, psig) keeps its offset:
// only for an absolute temperature or pressure.
fn (d Dimension) absolute() bool {
	return d == .temperature || d == .pressure
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
const calorie = 4.1868 // J, the IT calorie
const kilocalorie = 4186.8 // J, the IT kilocalorie
// The conventional millimetre of mercury: 13.5951 g/cm³ under standard
// gravity, so 13595.1 kg/m³ · 9.80665 m/s² · 1e-3 m.
const mm_hg = 13595.1 * standard_gravity * 1e-3 // Pa
const torr = atmosphere / 760.0 // Pa
const us_gallon = 231.0 * inch * inch * inch // m³
const horsepower = 550.0 * foot * pound_force // W, 550 ft·lbf/s
const minute = 60.0 // s
const hour = 3600.0 // s
const day = 86400.0 // s
const rankine = 5.0 / 9.0 // K per degree Fahrenheit or Rankine

// Dimensions the table uses.
const dim_length = dim(0, 1, 0, 0, 0, 0)
const dim_mass = dim(1, 0, 0, 0, 0, 0)
const dim_time = dim(0, 0, 1, 0, 0, 0)
const dim_temperature = dim(0, 0, 0, 1, 0, 0)
const dim_amount = dim(0, 0, 0, 0, 1, 0)
const dim_current = dim(0, 0, 0, 0, 0, 1)
const dim_charge = dim(0, 0, 1, 0, 0, 1)
const dim_volume = dim(0, 3, 0, 0, 0, 0)
const dim_force = dim(1, 1, -2, 0, 0, 0)
const dim_energy = dim(1, 2, -2, 0, 0, 0)
const dim_power = dim(1, 2, -3, 0, 0, 0)
const dim_pressure = dim(1, -1, -2, 0, 0, 0)
const dim_viscosity = dim(1, -1, -1, 0, 0, 0)

// Symbol is one row: SI = (x + pre)·scale + post. pre and post apply only
// when the symbol stands alone as an absolute temperature or pressure.
// prefixed is whether it takes SI prefixes (km, mmol, kPa, cP), and gauge
// whether it is a gauge pressure, which can only stand alone.
struct Symbol {
	name     string
	scale    f64
	dim      Dim
	pre      f64
	post     f64
	prefixed bool
	gauge    bool
}

// The kilogram stands alone, so that it is exactly 1, and the gram takes
// the other prefixes. A pound-mole is a pound's worth of moles, 453.59237
// mol, and the poise, 0.1 Pa·s, is there for the centipoise.
const symbols = [
	Symbol{'m', 1.0, dim_length, 0.0, 0.0, true, false},
	Symbol{'in', inch, dim_length, 0.0, 0.0, false, false},
	Symbol{'ft', foot, dim_length, 0.0, 0.0, false, false},
	Symbol{'kg', 1.0, dim_mass, 0.0, 0.0, false, false},
	Symbol{'g', 1e-3, dim_mass, 0.0, 0.0, true, false},
	Symbol{'t', 1e3, dim_mass, 0.0, 0.0, false, false},
	Symbol{'lb', pound, dim_mass, 0.0, 0.0, false, false},
	Symbol{'s', 1.0, dim_time, 0.0, 0.0, true, false},
	Symbol{'min', minute, dim_time, 0.0, 0.0, false, false},
	Symbol{'h', hour, dim_time, 0.0, 0.0, false, false},
	Symbol{'d', day, dim_time, 0.0, 0.0, false, false},
	Symbol{'K', 1.0, dim_temperature, 0.0, 0.0, true, false},
	Symbol{'degC', 1.0, dim_temperature, 273.15, 0.0, false, false},
	Symbol{'°C', 1.0, dim_temperature, 273.15, 0.0, false, false},
	Symbol{'degF', rankine, dim_temperature, 459.67, 0.0, false, false},
	Symbol{'°F', rankine, dim_temperature, 459.67, 0.0, false, false},
	Symbol{'R', rankine, dim_temperature, 0.0, 0.0, false, false},
	Symbol{'mol', 1.0, dim_amount, 0.0, 0.0, true, false},
	Symbol{'lbmol', 1e3 * pound, dim_amount, 0.0, 0.0, false, false},
	Symbol{'A', 1.0, dim_current, 0.0, 0.0, true, false},
	Symbol{'C', 1.0, dim_charge, 0.0, 0.0, true, false},
	Symbol{'L', 1e-3, dim_volume, 0.0, 0.0, true, false},
	Symbol{'gal', us_gallon, dim_volume, 0.0, 0.0, false, false},
	Symbol{'N', 1.0, dim_force, 0.0, 0.0, true, false},
	Symbol{'J', 1.0, dim_energy, 0.0, 0.0, true, false},
	Symbol{'cal', calorie, dim_energy, 0.0, 0.0, false, false},
	Symbol{'kcal', kilocalorie, dim_energy, 0.0, 0.0, false, false},
	Symbol{'Btu', btu, dim_energy, 0.0, 0.0, false, false},
	Symbol{'W', 1.0, dim_power, 0.0, 0.0, true, false},
	Symbol{'hp', horsepower, dim_power, 0.0, 0.0, false, false},
	Symbol{'Pa', 1.0, dim_pressure, 0.0, 0.0, true, false},
	Symbol{'bar', 1e5, dim_pressure, 0.0, 0.0, true, false},
	Symbol{'atm', atmosphere, dim_pressure, 0.0, 0.0, false, false},
	Symbol{'psi', psi, dim_pressure, 0.0, 0.0, false, false},
	Symbol{'psia', psi, dim_pressure, 0.0, 0.0, false, false},
	Symbol{'mmHg', mm_hg, dim_pressure, 0.0, 0.0, false, false},
	Symbol{'torr', torr, dim_pressure, 0.0, 0.0, false, false},
	Symbol{'barg', 1e5, dim_pressure, 0.0, atmosphere, false, true},
	Symbol{'psig', psi, dim_pressure, 0.0, atmosphere, false, true},
	Symbol{'P', 0.1, dim_viscosity, 0.0, 0.0, true, false},
	Symbol{'%', 1e-2, Dim{}, 0.0, 0.0, false, false},
]!

// Prefix is an SI prefix. pico and smaller, and tera and larger, are left
// out: no unit here needs them, and "ps" would otherwise read as a
// picosecond where "psi" was meant.
struct Prefix {
	name  string
	scale f64
}

const prefix_table = [
	Prefix{'G', 1e9},
	Prefix{'M', 1e6},
	Prefix{'k', 1e3},
	Prefix{'h', 1e2},
	Prefix{'d', 1e-1},
	Prefix{'c', 1e-2},
	Prefix{'m', 1e-3},
	Prefix{'µ', 1e-6}, // the micro sign, U+00B5
	Prefix{'μ', 1e-6}, // the Greek mu, U+03BC
	Prefix{'u', 1e-6},
	Prefix{'n', 1e-9},
]!

// symbol_names lists every unit symbol, before prefixes.
pub fn symbol_names() []string {
	return symbols.map(it.name)
}

// prefixed lists the symbols that take SI prefixes.
pub fn prefixed() []string {
	return symbols.filter(it.prefixed).map(it.name)
}

// prefixes lists the SI prefixes a symbol may take.
pub fn prefixes() []string {
	return prefix_table.map(it.name)
}

// lookup finds a symbol, or a prefix and a symbol that takes one: kmol, cP.
fn lookup(name string) ?Symbol {
	for s in symbols {
		if s.name == name {
			return s
		}
	}
	for p in prefix_table {
		if name.starts_with(p.name) {
			rest := name[p.name.len..]
			for s in symbols {
				if s.prefixed && s.name == rest {
					return Symbol{
						...s
						name:  name
						scale: p.scale * s.scale
					}
				}
			}
		}
	}
	return none
}

// every_name is every symbol, with every prefix where it takes one, for
// suggestions.
fn every_name() []string {
	mut out := symbol_names()
	for s in symbols.filter(it.prefixed) {
		for p in prefix_table {
			out << p.name + s.name
		}
	}
	return out
}

// Unit is a unit read: SI = (x + pre)·scale + post. pre and post are only
// set for a symbol standing alone, such as °C or barg.
pub struct Unit {
pub:
	text  string
	scale f64
	dim   Dim
	pre   f64
	post  f64
	gauge bool
}

// UnitError is a unit that couldn't be read: the text, the column (in
// characters, from 1) where it went wrong, why, and any units it may have
// meant. unknown is the name that isn't a unit, when that's why.
pub struct UnitError {
	Error
pub:
	text        string
	column      int
	reason      string
	suggestions []string
	unknown     string
}

// msg is the reason, with the text, the column and any suggestions.
pub fn (e UnitError) msg() string {
	hint := if e.suggestions.len > 0 { ' (did you mean ${e.suggestions.join(', ')}?)' } else { '' }
	return '"${e.text}", column ${e.column}: ${e.reason}${hint}'
}

// unit reads a unit: a symbol, with or without a prefix, or a compound of
// them.
//
//     unit   := factor (("*" | "·") factor | "/" factor)*
//     factor := ("1" | name | "(" unit ")") [exponent]
//     exponent := digits after a name (m3), superscripts (m³, s⁻¹), or "^"
//                 and an integer, a decimal or "(p/q)" (bar^-2, bar^0.5)
//
// Division is left to right, so J/mol/K is J/(mol·K), but a * after a / is
// refused as ambiguous: kmol/kg*h might mean kmol/(kg·h) or kmol·h/kg.
pub fn unit(text string) !Unit {
	return read_unit(text) or {
		// The fields one by one: V 0.5.2's ...err on a smart-cast error
		// copies nothing (docs/v-for-go.md, §4).
		if err is UnitError && err.unknown != '' {
			return UnitError{
				text:        err.text
				column:      err.column
				reason:      err.reason
				unknown:     err.unknown
				suggestions: closest(err.unknown, every_name())
			}
		}
		return err
	}
}

// read_unit is unit without suggestions, which parse makes from the
// dimension's own units instead.
fn read_unit(text string) !Unit {
	if text == '' {
		return UnitError{
			text:   text
			column: 1
			reason: 'it is empty'
		}
	}
	if s := lookup(text) {
		return Unit{text, s.scale, s.dim, s.pre, s.post, s.gauge}
	}
	mut r := Reader{
		text:  text
		runes: text.runes()
	}
	p := r.expr()!
	if r.pos < r.runes.len {
		return r.unexpected()
	}
	scale := p.num / p.den
	if !math.is_finite(scale) || !(scale > 0.0) {
		return r.fail(0, 'its scale, ${scale}, is out of range')
	}
	return Unit{
		text:  text
		scale: scale
		dim:   p.dim
	}
}

// Part is a piece of a unit being read: SI = num/den per unit. The scale is
// kept as a numerator and a denominator, each a product, so that kmol/h is
// 1e3/3600 in one division, as M4 computed it.
struct Part {
	num f64 = 1.0
	den f64 = 1.0
	dim Dim
}

fn (a Part) times(b Part) !Part {
	return Part{a.num * b.num, a.den * b.den, a.dim.mul(b.dim)!}
}

fn (a Part) over(b Part) !Part {
	return Part{a.num * b.den, a.den * b.num, a.dim.div(b.dim)!}
}

fn (p Part) pow(e Ratio) !Part {
	mut num, mut den := power(p.num, e), power(p.den, e)
	if e.num < 0 {
		num, den = den, num
	}
	return Part{num, den, p.dim.pow(e)!}
}

// power is x^|e|: by repeated multiplication for a whole power, whose bits
// IEEE arithmetic makes the same on every platform, else by math.pow, which
// is the C library's. math.pow can be an ulp nearer the exact value (1e-8
// for cm⁴, where the product is an ulp above), but CI runs on Windows and
// Linux, whose libraries differ.
fn power(x f64, e Ratio) f64 {
	n := if e.num < 0 { -e.num } else { e.num }
	if e.den == 1 && n <= 16 {
		mut y := 1.0
		for _ in 0 .. n {
			y *= x
		}
		return y
	}
	return math.pow(x, f64(n) / f64(e.den))
}

// The deepest nesting of parentheses read, so that no text can exhaust the
// stack.
const max_depth = 16

struct Reader {
	text  string
	runes []rune
mut:
	pos   int
	depth int
}

fn (r Reader) fail(at int, reason string) IError {
	return UnitError{
		text:   r.text
		column: at + 1
		reason: reason
	}
}

fn is_operator(c rune) bool {
	return c in [`*`, `·`, `⋅`, `/`, `(`, `)`, `^`, `[`, `]`, `.`] || (c >= `0` && c <= `9`)
		|| superscript_digit(c) >= 0 || c == `⁻`
}

fn superscript_digit(c rune) int {
	return match c {
		`⁰` { 0 }
		`¹` { 1 }
		`²` { 2 }
		`³` { 3 }
		`⁴` { 4 }
		`⁵` { 5 }
		`⁶` { 6 }
		`⁷` { 7 }
		`⁸` { 8 }
		`⁹` { 9 }
		else { -1 }
	}
}

fn (mut r Reader) expr() !Part {
	mut acc := r.factor()!
	mut divided := false
	for r.pos < r.runes.len {
		c := r.runes[r.pos]
		if c in [`*`, `·`, `⋅`] {
			if divided {
				return r.fail(r.pos, 'a * after a / is ambiguous: put what divides in parentheses, as in J/(mol*K)')
			}
			at := r.pos
			r.pos++
			acc = acc.times(r.factor()!) or { return r.fail(at, err.msg()) }
		} else if c == `/` {
			at := r.pos
			r.pos++
			acc = acc.over(r.factor()!) or { return r.fail(at, err.msg()) }
			divided = true
		} else {
			break
		}
	}
	return acc
}

fn (mut r Reader) factor() !Part {
	start := r.pos
	p, named := r.primary()!
	has, e := r.exponent(named)!
	if !has {
		return p
	}
	if r.pos < r.runes.len {
		c := r.runes[r.pos]
		if c == `^` || c == `⁻` || superscript_digit(c) >= 0 || (c >= `0` && c <= `9`) {
			return r.fail(r.pos, 'two exponents in a row')
		}
	}
	return p.pow(e) or { return r.fail(start, err.msg()) }
}

// primary reads a symbol, a parenthesised unit or the 1 of 1/s, and says
// whether it was a symbol, after which digits are an exponent (m3).
fn (mut r Reader) primary() !(Part, bool) {
	if r.pos >= r.runes.len {
		return r.fail(r.pos, 'a unit is missing at the end')
	}
	c := r.runes[r.pos]
	if c == `(` {
		if r.depth >= max_depth {
			return r.fail(r.pos, 'too many parentheses inside one another')
		}
		open := r.pos
		r.pos++
		r.depth++
		p := r.expr()!
		r.depth--
		if r.pos >= r.runes.len || r.runes[r.pos] != `)` {
			if r.pos < r.runes.len && r.runes[r.pos] != `)` {
				return r.unexpected()
			}
			return r.fail(open, 'a ( with no ) to close it')
		}
		r.pos++
		return p, false
	}
	if c == `1` && (r.pos + 1 >= r.runes.len || !(r.runes[r.pos + 1] >= `0` && r.runes[r.pos + 1] <= `9`)) {
		r.pos++
		return Part{}, false
	}
	if is_operator(c) {
		return r.unexpected()
	}
	start := r.pos
	for r.pos < r.runes.len && !is_operator(r.runes[r.pos]) {
		r.pos++
	}
	name := r.runes[start..r.pos].string()
	s := lookup(name) or {
		return UnitError{
			text:    r.text
			column:  start + 1
			reason:  'no unit ${name}'
			unknown: name
		}
	}
	if s.gauge {
		return r.fail(start, '${name} is a gauge pressure, which can only stand alone')
	}
	return Part{
		num: s.scale
		dim: s.dim
	}, true
}

// unexpected explains a character where a unit or an operator should be.
fn (r Reader) unexpected() IError {
	c := r.runes[r.pos]
	return match c {
		`[`, `]` { r.fail(r.pos, 'group units with ( and ), not [ and ]') }
		`.` { r.fail(r.pos, 'put * or · between units, not .') }
		`)` { r.fail(r.pos, 'a ) with no ( before it') }
		`^` { r.fail(r.pos, 'an exponent with nothing to raise') }
		`*`, `·`, `⋅`, `/` { r.fail(r.pos, 'a unit is missing before ${c}') }
		else {
			if (c >= `0` && c <= `9`) || superscript_digit(c) >= 0 || c == `⁻` {
				r.fail(r.pos, 'a number can only be an exponent, or the 1 of 1/s')
			} else {
				r.fail(r.pos, 'unexpected ${c}')
			}
		}
	}
}

// exponent reads an exponent, if one follows: digits right after a symbol
// (m3), superscripts (m³, s⁻¹), or ^ and an integer, a decimal or a
// fraction in parentheses.
fn (mut r Reader) exponent(named bool) !(bool, Ratio) {
	if r.pos >= r.runes.len {
		return false, Ratio{}
	}
	c := r.runes[r.pos]
	if named && c >= `0` && c <= `9` {
		start := r.pos
		for r.pos < r.runes.len && r.runes[r.pos] >= `0` && r.runes[r.pos] <= `9` {
			r.pos++
		}
		e := r.integer(start, r.runes[start..r.pos].string(), false)!
		return true, e
	}
	if c == `⁻` || superscript_digit(c) >= 0 {
		start := r.pos
		negative := c == `⁻`
		if negative {
			r.pos++
		}
		mut digits := ''
		for r.pos < r.runes.len && superscript_digit(r.runes[r.pos]) >= 0 {
			digits += superscript_digit(r.runes[r.pos]).str()
			r.pos++
		}
		if digits == '' {
			return r.fail(start, 'a ⁻ with no power after it')
		}
		e := r.integer(start, digits, negative)!
		return true, e
	}
	if c != `^` {
		return false, Ratio{}
	}
	caret := r.pos
	r.pos++
	if r.pos < r.runes.len && r.runes[r.pos] == `(` {
		r.pos++
		num := r.signed_digits()!
		if r.pos >= r.runes.len || r.runes[r.pos] != `/` {
			return r.fail(caret, 'a power in parentheses is a fraction, such as ^(1/2)')
		}
		r.pos++
		den := r.signed_digits()!
		if r.pos >= r.runes.len || r.runes[r.pos] != `)` {
			return r.fail(caret, 'a power in parentheses is a fraction, such as ^(1/2)')
		}
		r.pos++
		e := ratio(num, den) or { return r.fail(caret, err.msg()) }
		return true, e
	}
	negative := r.pos < r.runes.len && r.runes[r.pos] == `-`
	if r.pos < r.runes.len && (r.runes[r.pos] == `-` || r.runes[r.pos] == `+`) {
		r.pos++
	}
	start := r.pos
	for r.pos < r.runes.len && r.runes[r.pos] >= `0` && r.runes[r.pos] <= `9` {
		r.pos++
	}
	whole := r.runes[start..r.pos].string()
	mut fraction := ''
	if r.pos < r.runes.len && r.runes[r.pos] == `.` {
		r.pos++
		from := r.pos
		for r.pos < r.runes.len && r.runes[r.pos] >= `0` && r.runes[r.pos] <= `9` {
			r.pos++
		}
		fraction = r.runes[from..r.pos].string()
		if fraction == '' {
			return r.fail(caret, 'a power needs digits after its decimal point')
		}
	}
	if whole == '' && fraction == '' {
		return r.fail(caret, 'a ^ with no power after it')
	}
	if whole.len > 7 || fraction.len > 6 {
		return r.fail(caret, 'a power with too many digits')
	}
	scale := i64(math.powi(10, fraction.len))
	mut num := (if whole == '' { i64(0) } else { strconv.parse_int(whole, 10, 64)! }) * scale
	if fraction != '' {
		num += strconv.parse_int(fraction, 10, 64)!
	}
	if negative {
		num = -num
	}
	e := ratio(num, scale) or { return r.fail(caret, err.msg()) }
	return true, e
}

fn (mut r Reader) signed_digits() !i64 {
	start := r.pos
	negative := r.pos < r.runes.len && r.runes[r.pos] == `-`
	if negative {
		r.pos++
	}
	from := r.pos
	for r.pos < r.runes.len && r.runes[r.pos] >= `0` && r.runes[r.pos] <= `9` {
		r.pos++
	}
	digits := r.runes[from..r.pos].string()
	if digits == '' || digits.len > 7 {
		return r.fail(start, 'a power in parentheses is a fraction of whole numbers, such as ^(1/2)')
	}
	n := strconv.parse_int(digits, 10, 64)!
	return if negative { -n } else { n }
}

fn (r Reader) integer(start int, digits string, negative bool) !Ratio {
	if digits.len > 7 {
		return r.fail(start, 'a power with too many digits')
	}
	n := strconv.parse_int(digits, 10, 64)!
	return ratio(if negative { -n } else { n }, 1) or { return r.fail(start, err.msg()) }
}

// QuantityError is every failure of a number with its unit: the text that
// couldn't be read, the dimension it was read as, why, and any units it may
// have meant.
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

// units lists the common units of a named dimension: M4's, which a message
// offers and suggests from. Any unit of the dimension reads.
pub fn units(d Dimension) []string {
	return named_units.filter(it.dim == d).map(it.name)
}

struct Named {
	name string
	dim  Dimension
}

const named_units = [
	Named{'K', .temperature},
	Named{'degC', .temperature},
	Named{'°C', .temperature},
	Named{'degF', .temperature},
	Named{'°F', .temperature},
	Named{'R', .temperature},
	Named{'K', .temperature_difference},
	Named{'degC', .temperature_difference},
	Named{'°C', .temperature_difference},
	Named{'degF', .temperature_difference},
	Named{'°F', .temperature_difference},
	Named{'R', .temperature_difference},
	Named{'Pa', .pressure},
	Named{'kPa', .pressure},
	Named{'MPa', .pressure},
	Named{'bar', .pressure},
	Named{'mbar', .pressure},
	Named{'atm', .pressure},
	Named{'psi', .pressure},
	Named{'psia', .pressure},
	Named{'mmHg', .pressure},
	Named{'torr', .pressure},
	Named{'barg', .pressure},
	Named{'psig', .pressure},
	Named{'Pa', .pressure_difference},
	Named{'kPa', .pressure_difference},
	Named{'MPa', .pressure_difference},
	Named{'bar', .pressure_difference},
	Named{'mbar', .pressure_difference},
	Named{'atm', .pressure_difference},
	Named{'psi', .pressure_difference},
	Named{'mmHg', .pressure_difference},
	Named{'torr', .pressure_difference},
	Named{'mol/s', .molar_flow},
	Named{'mol/min', .molar_flow},
	Named{'mol/h', .molar_flow},
	Named{'kmol/s', .molar_flow},
	Named{'kmol/min', .molar_flow},
	Named{'kmol/h', .molar_flow},
	Named{'lbmol/h', .molar_flow},
	Named{'kg/s', .mass_flow},
	Named{'kg/min', .mass_flow},
	Named{'kg/h', .mass_flow},
	Named{'g/s', .mass_flow},
	Named{'t/h', .mass_flow},
	Named{'lb/s', .mass_flow},
	Named{'lb/h', .mass_flow},
	Named{'m3/s', .volume_flow},
	Named{'m³/s', .volume_flow},
	Named{'m3/h', .volume_flow},
	Named{'m³/h', .volume_flow},
	Named{'L/s', .volume_flow},
	Named{'L/min', .volume_flow},
	Named{'L/h', .volume_flow},
	Named{'gal/min', .volume_flow},
	Named{'W', .power},
	Named{'kW', .power},
	Named{'MW', .power},
	Named{'kJ/s', .power},
	Named{'kJ/h', .power},
	Named{'MJ/h', .power},
	Named{'GJ/h', .power},
	Named{'kcal/h', .power},
	Named{'Btu/h', .power},
	Named{'hp', .power},
	Named{'J/mol', .molar_enthalpy},
	Named{'kJ/mol', .molar_enthalpy},
	Named{'kJ/kmol', .molar_enthalpy},
	Named{'MJ/kmol', .molar_enthalpy},
	Named{'Btu/lbmol', .molar_enthalpy},
	Named{'J/kg', .specific_enthalpy},
	Named{'kJ/kg', .specific_enthalpy},
	Named{'MJ/kg', .specific_enthalpy},
	Named{'Btu/lb', .specific_enthalpy},
	Named{'%', .fraction},
]!

// parse reads "value unit" as a quantity of dimension d and returns it in
// SI: K, Pa, mol/s, kg/s, m³/s, W, J/mol, J/kg, or a plain fraction. A
// fraction may be a bare number; everything else needs its unit, after a
// space. Any unit of d's dimension reads, compound or not. An absolute
// temperature or pressure must be above zero.
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
	u := read(text, d, parts[1])!
	pre, post := if d.absolute() { u.pre, u.post } else { 0.0, 0.0 }
	si := (value + pre) * u.scale + post
	if !math.is_finite(si) {
		return fail(text, d, 'that is beyond the range of a number', []string{})
	}
	if d == .temperature && !(si > 0.0) {
		return fail(text, d, 'that is ${round6(si)} K, not above absolute zero', []string{})
	}
	if d == .pressure && !(si > 0.0) {
		return fail(text, d, 'that is ${round6(si)} Pa, not above zero', []string{})
	}
	return si
}

// read reads a unit of dimension d, with the reason if it isn't one.
fn read(text string, d Dimension, name string) !Unit {
	u := read_unit(name) or {
		if err is UnitError {
			if err.unknown != '' {
				return unknown(text, d, name)
			}
			return fail(text, d, '${name}: ${err.reason}, at character ${err.column}', []string{})
		}
		return err
	}
	if !u.dim.eq(d.dim()) || (u.gauge && d != .pressure) {
		what := if u.gauge { 'pressure' } else { u.dim.name() }
		return fail(text, d, '${name} is a unit of ${what}, not ${d.name()}', closest(name, units(d)))
	}
	return u
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

// unknown explains a unit that isn't one: a near miss of the dimension's
// common units, or else what they are.
fn unknown(text string, d Dimension, name string) IError {
	accepted := units(d)
	near := closest(name, accepted)
	if near.len > 0 {
		return fail(text, d, 'no unit ${name}', near)
	}
	return fail(text, d, 'no unit ${name}; the units of ${d.name()} are ${accepted.join(', ')}',
		[]string{})
}

// closest is the accepted units within two edits of name, ignoring case:
// nearest first, and at most three.
pub fn closest(name string, accepted []string) []string {
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

// to_unit converts an SI value to a unit of dimension d: the inverse of
// parse, for reports.
pub fn to_unit(si f64, name string, d Dimension) !f64 {
	u := read('${si} ${name}', d, name)!
	pre, post := if d.absolute() { u.pre, u.post } else { 0.0, 0.0 }
	return (si - post) / u.scale - pre
}
