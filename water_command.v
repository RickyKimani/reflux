// `reflux water` (M5a, D13): a steam table on the command line.
//
//     reflux water SPEC SPEC [--if97] [--json OUT]
module main

import flash
import flowsheet
import math
import quantity
import strconv
import thermo
import water
import water.if97
import json2

const water_usage = "usage:
  reflux water SPEC SPEC [--if97] [--json OUT]
      Water's properties at the state two SPECs fix, by IAPWS-95, or by
      IAPWS-IF97 with --if97, per kg and per mole, on IAPWS's reference
      (u = s = 0 for the liquid at the triple point). --json writes them
      (docs/report-format.md). A SPEC is a name, = and a number with its
      unit, with or without a space between them:
        T=450degC  P=40bar  h=3000kJ/kg  s=6.5kJ/(kg*K)  x=0.9  rho=500kg/m3
      h, s and rho may be molar: h=54kJ/mol. The pairs that fix one state
      are T and P, T and rho, T and x, P and x, P and h, and P and s; IF97
      takes each but T and rho.

The exit status is 0 when it worked, 1 when water has no state there (or
the formulation none) or the report could not be written, and 2 when the
command is wrong, a pair that doesn't fix one state among them."

// Given is a property a state is asked from, in the order a pair is named.
enum Given {
	t
	p
	h
	s
	x
	rho
}

fn (g Given) str() string {
	return match g {
		.t { 'T' }
		.p { 'P' }
		.h { 'h' }
		.s { 's' }
		.x { 'x' }
		.rho { 'rho' }
	}
}

// Spec is one SPEC: what it gives, as typed, and its value in SI per kg
// (K, Pa, J/kg, J/(kg·K), a fraction, kg/m³).
struct Spec {
	given Given
	text  string
	value f64
}

struct WaterOptions {
mut:
	specs []Spec
	if97  bool
	json  string
}

// WrongPair is a pair of SPECs that doesn't fix one state: a wrong command,
// not a state water lacks.
struct WrongPair {
	Error
	reason string
}

fn (e WrongPair) msg() string {
	return e.reason
}

// WaterPhase is one phase's properties per kg, on IAPWS's reference; a
// property it has none of is NaN (a mixture's cp, a viscosity out of its
// correlation's range).
struct WaterPhase {
	rho    f64 // kg/m³
	v      f64 // m³/kg
	h      f64 // J/kg
	s      f64 // J/(kg·K)
	u      f64 // J/kg
	cv     f64 // J/(kg·K)
	cp     f64 // J/(kg·K)
	w      f64 // m/s
	mu     f64 // Pa·s
	lambda f64 // W/(m·K)
}

// WaterState is the state two SPECs fix: one phase, or two (the saturated
// liquid and vapour, with the mixture as whole).
struct WaterState {
	formulation string // IAPWS-95 or IAPWS-IF97
	specs       []Spec
	t           f64 // K
	p           f64 // Pa
	phase       string
	x           f64 // the vapour's mass fraction, NaN off saturation
	region      int // IF97's, 0 by IAPWS-95
	whole       WaterPhase
	liquid      ?WaterPhase
	vapour      ?WaterPhase
	sigma       f64 // N/m at a saturated state, else NaN
	notes       []string
}

// water_command runs `reflux water`, and returns its exit status.
fn water_command(args []string, mut c Console) int {
	if args.len > 0 && args[0] in ['help', '--help', '-h'] {
		c.say(water_usage)
		return 0
	}
	o := water_options(args) or {
		c.complain('reflux water: ${err.msg()}\n\n${water_usage}')
		return 2
	}
	st := steam(o.specs[0], o.specs[1], o.if97) or {
		c.complain('reflux water: ${err.msg()}')
		return if err is WrongPair { 2 } else { 1 }
	}
	c.say(st.table())
	if o.json != '' && !write(o.json, st.json(), mut c) {
		return 1
	}
	return 0
}

fn water_options(args []string) !WaterOptions {
	mut o := WaterOptions{}
	mut i := 0
	for i < args.len {
		arg := args[i]
		i++
		name := arg.all_before('=')
		if name == '--json' {
			mut value := arg.all_after_first('=')
			if !arg.contains('=') {
				if i >= args.len || args[i].starts_with('-') {
					return error('--json needs a file name')
				}
				value = args[i]
				i++
			}
			if value == '' {
				return error('--json needs a file name')
			}
			if o.json != '' {
				return error('--json is given twice')
			}
			o.json = value
		} else if arg == '--if97' {
			o.if97 = true
		} else if arg.starts_with('--') {
			return error('no option ${arg}${did_you_mean(name, ['--json', '--if97'])}')
		} else {
			// A negative number's minus can't start a SPEC: each has its name.
			o.specs << read_spec(arg)!
		}
	}
	if o.specs.len != 2 {
		return error('it needs two SPECs, not ${o.specs.len}')
	}
	if o.specs[0].given == o.specs[1].given {
		return error('${o.specs[0].given} is given twice')
	}
	return o
}

// Each SPEC's dimensions, per kg and per mole.
const specific_entropy = quantity.dim(0, 2, -2, -1, 0, 0)
const molar_entropy = quantity.dim(1, 2, -2, -1, -1, 0)
const mass_density = quantity.dim(1, -3, 0, 0, 0, 0)
const molar_density = quantity.dim(0, -3, 0, 0, 1, 0)

// read_spec reads NAME=VALUE, the value a number and its unit, with or
// without a space between them, and x a bare fraction or a percentage.
fn read_spec(arg string) !Spec {
	if !arg.contains('=') {
		return error('${arg}: a SPEC is a name, = and a value, as T=450degC')
	}
	name := arg.all_before('=').trim_space()
	text := arg.all_after_first('=').trim_space()
	given := match name.to_lower() {
		't' {
			Given.t
		}
		'p' {
			Given.p
		}
		'h' {
			Given.h
		}
		's' {
			Given.s
		}
		'x' {
			Given.x
		}
		'rho', 'ρ' {
			Given.rho
		}
		else {
			return error('${arg}: no property ${name}${did_you_mean(name, ['T', 'P', 'h', 's',
				'x', 'rho'])}; a state is two of T, P, h, s, x and rho')
		}
	}
	mut i := 0
	for i < text.len && (text[i].is_digit() || text[i] in [`.`, `-`, `+`, `e`, `E`]) {
		i++
	}
	// An e that begins a unit isn't the number's: 2e5Pa is 2e5 Pa, but 2eV
	// would be 2 eV.
	for i > 0 && text[i - 1] in [`e`, `E`, `+`, `-`] {
		i--
	}
	num := text[..i]
	u := text[i..].trim_space()
	v := strconv.atof64(num) or { return error('${arg}: "${num}" is not a number') }
	if !math.is_finite(v) {
		return error('${arg}: "${num}" is not a finite number')
	}
	value := match given {
		.t {
			quantity.parse('${num} ${u}', .temperature) or { return error('${arg}: ${err.msg()}') }
		}
		.p {
			quantity.parse('${num} ${u}', .pressure) or { return error('${arg}: ${err.msg()}') }
		}
		.x {
			x := if u == '' {
				v
			} else if u == '%' {
				v / 100.0
			} else {
				return error('${arg}: x is a fraction, a number from 0 to 1 or a percentage')
			}
			if !(x >= 0.0 && x <= 1.0) {
				return error('${arg}: x is from 0 to 1')
			}
			x
		}
		.h {
			per(arg, v, u, quantity.Dimension.specific_enthalpy.dim(), quantity.Dimension.molar_enthalpy.dim(),
				1.0 / water.molar_mass, 'kJ/kg or kJ/mol')!
		}
		.s {
			per(arg, v, u, specific_entropy, molar_entropy, 1.0 / water.molar_mass, 'kJ/(kg*K) or J/(mol*K)')!
		}
		.rho {
			r := per(arg, v, u, mass_density, molar_density, water.molar_mass, 'kg/m3 or mol/m3')!
			if !(r > 0.0) {
				return error('${arg}: a density is above zero')
			}
			r
		}
	}
	return Spec{
		given: given
		text:  text
		value: value
	}
}

// per is a value in a unit of either dimension, per kg in SI: a molar one
// times to_kg, 1/M for an energy per mole, M for moles per m³.
fn per(arg string, v f64, u string, specific quantity.Dim, molar quantity.Dim, to_kg f64, units string) !f64 {
	if u == '' {
		return error('${arg}: it needs a unit, such as ${units}')
	}
	unit := quantity.unit(u) or { return error('${arg}: ${err.msg()}') }
	if unit.dim.eq(specific) {
		return v * unit.scale
	}
	if unit.dim.eq(molar) {
		return v * unit.scale * to_kg
	}
	return error('${arg}: ${u} is a unit of ${unit.dim.name()}, not of ${units}')
}

// steam is the state two SPECs fix, by IAPWS-95 or IF97. A pair that
// doesn't fix one state is a WrongPair, whatever the values.
fn steam(a Spec, b Spec, use_if97 bool) !WaterState {
	first, second := if int(a.given) < int(b.given) { a, b } else { b, a }
	check_pair(first.given, second.given, use_if97)!
	mut st := if use_if97 { by_if97(first, second)! } else { by_iapws95(first, second)! }
	st = WaterState{
		...st
		specs: [a, b]
	}
	return with_ice(st)
}

// check_pair refuses a pair that can fix two states, or none, or that
// reflux doesn't solve, saying why. The pairs are in Given's order.
fn check_pair(a Given, b Given, use_if97 bool) ! {
	match true {
		a == .t && b == .rho && use_if97 {
			return WrongPair{
				reason: 'IF97 has no state from T and rho: its regions 1, 2 and 5 are equations in P and T. Give P, or leave out --if97 for IAPWS-95.'
			}
		}
		(a == .t && b in [.p, .x, .rho]) || (a == .p && b in [.h, .s, .x]) {
			return
		}
		a == .t && b == .h {
			return WrongPair{
				reason: 'T and h fix no state: on an isotherm h hardly moves with P (in the ideal gas not at all). Give P, x or rho with T.'
			}
		}
		a == .t && b == .s {
			return WrongPair{
				reason: "T and s can fix two states: below 4 °C liquid water's s rises with P to its density maximum and falls beyond it. Give P, x or rho with T."
			}
		}
		a == .p && b == .rho {
			return WrongPair{
				reason: "P and rho can fix two states: liquid water's density has a maximum, near 4 °C at low pressure. Give T with rho, or h or s with P."
			}
		}
		b == .x {
			// h, s or rho with x: along the saturation line each has an
			// extremum at some x (the saturated vapour's h at about 2803
			// kJ/kg, near 235 °C), so two temperatures can share it.
			return WrongPair{
				reason: '${a} and x can fix two states: along the saturation line ${a} at a given x can rise and fall. Give T or P with x.'
			}
		}
		else {
			return WrongPair{
				reason: "reflux water doesn't solve for T and P from ${a} and ${b} yet. Give T or P with one of them."
			}
		}
	}
}

// IAPWS-95's range in reflux: from 251.165 K, where its model has no
// state below (M5a's D7), to 1273.15 K and 1000 MPa (R6-95, section 6).
const t_max_95 = 1273.15
const p_max_95 = 1e9

fn in_range_95(t f64, p f64) ! {
	if !(t >= water.liquid_t_min) {
		return error('IAPWS-95 has no state below ${water.liquid_t_min} K in reflux (the lowest T of liquid water, at ice III, and the edge of what its model holds): not ${shown(t)} K')
	}
	if !(t <= t_max_95) {
		return error('IAPWS-95 holds to ${t_max_95} K (R6-95): not ${shown(t)} K')
	}
	if !(p <= p_max_95) {
		return error('IAPWS-95 holds to 1000 MPa (R6-95): not ${shown(p * 1e-6)} MPa')
	}
}

// by_iapws95 solves a pair on IAPWS-95: (T, P) by the PT flash, (T, ρ)
// directly or by the lever rule inside the dome, T or P with x by the
// saturated phases, P with h or s by the PH and PS flashes. Each phase's
// properties are then the release's Table 3 at its (T, ρ).
fn by_iapws95(a Spec, b Spec) !WaterState {
	m := thermo.by_roots(water.iapws95())
	z := [1.0]
	mut notes := []string{}
	// What's given is checked first, so that the flash isn't asked.
	if a.given == .t {
		in_range_95(a.value, 0.0)!
	}
	if a.given == .p {
		in_range_95(water.liquid_t_min, a.value)!
	}
	match b.given {
		.p {
			t, p := a.value, b.value
			in_range_95(t, p)!
			r := flash.pt(m, t, p, z) or { return error('no state at ${shown(t)} K and ${shown(p * 1e-6)} MPa: ${err.msg()}') }
			return one95(t, p, 1.0 / r.phases[0].state.v, mut notes)
		}
		.rho {
			t, rho := a.value, b.value / water.molar_mass
			if !(rho <= 4.0 * water.critical_rho) {
				return error("IAPWS-95 has no state denser than 4 ρc (${4.0 * water.critical_rho_mass} kg/m³) in reflux (M5a's D7): not ${shown(b.value)} kg/m³")
			}
			if t < water.critical_t {
				sat := saturated95(m, t, true)!
				if rho < 1.0 / sat.liquid.v && rho > 1.0 / sat.vapour.v {
					v := 1.0 / rho
					x := (v - sat.liquid.v) / (sat.vapour.v - sat.liquid.v)
					return two95(t, sat.p, sat.liquid.v, sat.vapour.v, x, mut notes)
				}
			}
			p := water.props95(t, rho).p
			if !(p > 0.0) {
				return error('no state at ${shown(t)} K and ${shown(b.value)} kg/m³: its pressure would be ${shown(p * 1e-6)} MPa')
			}
			in_range_95(t, p)!
			return one95(t, p, rho, mut notes)
		}
		.x {
			sat := saturated95(m, a.value, a.given == .t)!
			return two95(sat.t, sat.p, sat.liquid.v, sat.vapour.v, b.value, mut notes)
		}
		.h, .s {
			fl := water.fluid(none, none)!
			dh, ds := water.reference_offsets(fl)
			p := a.value
			r := if b.given == .h {
				flash.ph(fl, p, b.value * water.molar_mass + dh, z) or {
					return error('no state at ${shown(p * 1e-6)} MPa and h = ${shown(b.value * 1e-3)} kJ/kg: ${err.msg()}')
				}
			} else {
				flash.ps(fl, p, b.value * water.molar_mass + ds, z) or {
					return error('no state at ${shown(p * 1e-6)} MPa and s = ${shown(b.value * 1e-3)} kJ/(kg K): ${err.msg()}')
				}
			}
			in_range_95(r.t, p)!
			if r.phases.len == 2 {
				vapour := r.phases.filter(it.label == .vapor)
				liquid := r.phases.filter(it.label != .vapor)
				if vapour.len == 1 && liquid.len == 1 {
					return two95(r.t, p, liquid[0].state.v, vapour[0].state.v, r.beta, mut notes)
				}
			}
			return one95(r.t, p, 1.0 / r.phases[0].state.v, mut notes)
		}
		else {
			return error('no pair ${a.given} and ${b.given}') // check_pair refused it
		}
	}
}

struct Saturated95 {
	t      f64
	p      f64
	liquid thermo.State
	vapour thermo.State
}

// saturated95 is the saturated liquid and vapour at T (by_t) or at P.
fn saturated95(m thermo.PhaseModel, value f64, by_t bool) !Saturated95 {
	if by_t && !(value < water.critical_t) {
		return error('water has no saturated states at or above its critical point, ${water.critical_t} K: not ${shown(value)} K')
	}
	if !by_t && !(value < water.critical_p) {
		return error('water has no saturated states at or above its critical point, ${shown(water.critical_p * 1e-6)} MPa: not ${shown(value * 1e-6)} MPa')
	}
	r := if by_t {
		flash.tvf(m, value, 0.0, [1.0]) or {
			return error('no saturated states at ${shown(value)} K: ${err.msg()}')
		}
	} else {
		flash.pvf(m, value, 0.0, [1.0]) or {
			return error('no saturated states at ${shown(value * 1e-6)} MPa: ${err.msg()}')
		}
	}
	vapour := r.phases.filter(it.label == .vapor)
	liquid := r.phases.filter(it.label != .vapor)
	if vapour.len != 1 || liquid.len != 1 {
		return error('no saturated phases at ${shown(r.t)} K')
	}
	in_range_95(r.t, r.p)!
	return Saturated95{
		t:      r.t
		p:      r.p
		liquid: liquid[0].state
		vapour: vapour[0].state
	}
}

// phase95 is one phase's properties at T and the molar density rho, by
// IAPWS-95's Table 3 and the transport releases (R12-08, R15-11).
fn phase95(t f64, rho f64, mut notes []string) WaterPhase {
	pr := water.props95(t, rho)
	mm := water.molar_mass
	rho_kg := rho * mm
	mu := water.viscosity(t, rho_kg) or {
		notes << 'mu: ${err.msg()}'
		math.nan()
	}
	lambda := water.thermal_conductivity(t, rho_kg) or {
		notes << 'lambda: ${err.msg()}'
		math.nan()
	}
	return WaterPhase{
		rho:    rho_kg
		v:      1.0 / rho_kg
		h:      pr.h / mm
		s:      pr.s / mm
		u:      pr.u / mm
		cv:     pr.cv / mm
		cp:     pr.cp / mm
		w:      pr.w
		mu:     mu
		lambda: lambda
	}
}

fn one95(t f64, p f64, rho f64, mut notes []string) WaterState {
	ph := phase95(t, rho, mut notes)
	return WaterState{
		formulation: 'IAPWS-95'
		t:           t
		p:           p
		phase:       label(t, p, ph.rho)
		x:           math.nan()
		whole:       ph
		sigma:       math.nan()
		notes:       notes
	}
}

// two95 is the mixture of the saturated phases of molar volumes vl and vv
// at T and P, x the vapour's fraction; x = 0 or 1 is one saturated phase.
fn two95(t f64, p f64, vl f64, vv f64, x f64, mut notes []string) WaterState {
	liquid := phase95(t, 1.0 / vl, mut notes)
	vapour := phase95(t, 1.0 / vv, mut notes)
	return saturated_state('IAPWS-95', t, p, x, 0, liquid, vapour, mut notes)
}

// saturated_state is a saturated state of quality x from its phases: one
// of them at x = 0 or 1, the lever rule's mixture between.
fn saturated_state(formulation string, t f64, p f64, x f64, region int, liquid WaterPhase, vapour WaterPhase, mut notes []string) WaterState {
	sigma := water.surface_tension(t) or {
		notes << 'sigma: ${err.msg()}'
		math.nan()
	}
	if x == 0.0 || x == 1.0 {
		return WaterState{
			formulation: formulation
			t:           t
			p:           p
			phase:       if x == 0.0 { 'saturated liquid' } else { 'saturated vapour' }
			x:           x
			region:      region
			whole:       if x == 0.0 { liquid } else { vapour }
			sigma:       sigma
			notes:       notes
		}
	}
	nan := math.nan()
	v := liquid.v + x * (vapour.v - liquid.v)
	return WaterState{
		formulation: formulation
		t:           t
		p:           p
		phase:       'liquid and vapour'
		x:           x
		region:      region
		whole:       WaterPhase{
			rho:    1.0 / v
			v:      v
			h:      liquid.h + x * (vapour.h - liquid.h)
			s:      liquid.s + x * (vapour.s - liquid.s)
			u:      liquid.u + x * (vapour.u - liquid.u)
			cv:     nan
			cp:     nan
			w:      nan
			mu:     nan
			lambda: nan
		}
		liquid:      liquid
		vapour:      vapour
		sigma:       sigma
		notes:       notes
	}
}

// label is a single phase's name: supercritical above both critical T and
// P, vapour above Tc below pc, and below Tc liquid or vapour by its side
// of ρc, as IAPWS-95's branches are split (M5a's D7).
fn label(t f64, p f64, rho f64) string {
	if t >= water.critical_t {
		return if p >= water.critical_p { 'supercritical' } else { 'vapour' }
	}
	return if rho > water.critical_rho_mass { 'liquid' } else { 'vapour' }
}

// by_if97 solves a pair by IF97's functions, which refuse what's outside
// its range; a wet state's phases are its saturated states at x = 0 and 1.
fn by_if97(a Spec, b Spec) !WaterState {
	st := state97(a, b) or { return error('no state by IF97: ${err.msg()}') }
	mut notes := []string{}
	if st.region == 4 {
		liquid := if a.given == .t { if97.state_tx(st.t, 0.0)! } else { if97.state_px(st.p, 0.0)! }
		vapour := if a.given == .t { if97.state_tx(st.t, 1.0)! } else { if97.state_px(st.p, 1.0)! }
		return saturated_state('IAPWS-IF97', st.t, st.p, st.x, 4, phase97(liquid, mut notes),
			phase97(vapour, mut notes), mut notes)
	}
	ph := phase97(st, mut notes)
	return WaterState{
		formulation: 'IAPWS-IF97'
		t:           st.t
		p:           st.p
		phase:       label(st.t, st.p, ph.rho)
		x:           math.nan()
		region:      st.region
		whole:       ph
		sigma:       math.nan()
		notes:       notes
	}
}

// state97 is IF97's state at a pair, in Given's order.
fn state97(a Spec, b Spec) !if97.State {
	if a.given == .t && b.given == .p {
		return if97.state_pt(b.value, a.value)!
	}
	if a.given == .t && b.given == .x {
		return if97.state_tx(a.value, b.value)!
	}
	if a.given == .p && b.given == .x {
		return if97.state_px(a.value, b.value)!
	}
	if a.given == .p && b.given == .h {
		return if97.state_ph(a.value, b.value)!
	}
	if a.given == .p && b.given == .s {
		return if97.state_ps(a.value, b.value)!
	}
	return error('no pair ${a.given} and ${b.given}') // check_pair refused it
}

// phase97 is one IF97 state's properties, with the transport releases'
// industrial forms on its density (R12-08's μ₂ = 1, R15-11's Eq. 25).
fn phase97(st if97.State, mut notes []string) WaterPhase {
	rho := st.rho()
	mu := water.viscosity_industrial(st.t, rho) or {
		notes << 'mu: ${err.msg()}'
		math.nan()
	}
	lambda := if c := water.conductivity_industrial(st) {
		c.lambda
	} else {
		notes << 'lambda: ${err.msg()}'
		math.nan()
	}
	return WaterPhase{
		rho:    rho
		v:      st.v
		h:      st.h
		s:      st.s
		u:      st.u
		cv:     st.cv
		cp:     st.cp
		w:      st.w
		mu:     mu
		lambda: lambda
	}
}

// with_ice refuses a state where ice is the only stable phase above
// liquid water's melting pressure (ice III, V or VI by R14-08), and notes
// one below 273.16 K that ice Ih makes metastable.
fn with_ice(st WaterState) !WaterState {
	if top, ice := melting_top(st.t) {
		if st.p > top {
			return error('ice ${ice} at ${shown(st.t)} K and ${shown(st.p * 1e-6)} MPa: above ${shown(top * 1e-6)} MPa, its melting pressure there (R14-08), water is no fluid')
		}
	}
	if st.t < triple_t {
		mut notes := st.notes.clone()
		if !math.is_nan(st.x) {
			notes << 'below the triple point, 273.16 K, ice Ih is the stable phase: this liquid and vapour are metastable'
		} else if st.p > sublimation_pressure(st.t) && st.p < ih_melting_pressure(st.t) {
			notes << 'ice Ih is the stable phase here, between its sublimation pressure (${shown(sublimation_pressure(st.t))} Pa) and its melting pressure (${shown(ih_melting_pressure(st.t) * 1e-6)} MPa) at this T (R14-08): this ${st.phase} is metastable'
		}
		return WaterState{
			...st
			notes: notes
		}
	}
	return st
}

// R14-08's melting and sublimation curves of ice: each p(T) as its Eqs. 1
// to 6 have it, in Pa.
const triple_t = 273.16
const triple_p = 611.657

fn ih_melting_pressure(t f64) f64 {
	theta := t / triple_t
	return triple_p * (1.0 + 0.119539337e7 * (1.0 - math.pow(theta, 3.0)) +
		0.808183159e5 * (1.0 - math.pow(theta, 25.75)) +
		0.333826860e4 * (1.0 - math.pow(theta, 103.75)))
}

fn sublimation_pressure(t f64) f64 {
	theta := t / triple_t
	return triple_p * math.exp((-0.212144006e2 * math.pow(theta, 0.333333333e-2) +
		0.273203819e2 * math.pow(theta, 0.120666667e1) +
		-0.610598130e1 * math.pow(theta, 0.170333333e1)) / theta)
}

// melting_top is the melting pressure above which ice bounds liquid water
// at T, and that ice's name, from 251.165 K to 355 K; above 355 K it's ice
// VII's, beyond 2216 MPa and IAPWS-95's range, so none.
fn melting_top(t f64) ?(f64, string) {
	if t < 251.165 || t > 355.0 {
		return none
	}
	if t <= 256.164 {
		return 208.566e6 * (1.0 - 0.299948 * (1.0 - math.pow(t / 251.165, 60.0))), 'III'
	}
	if t <= 273.31 {
		return 350.100e6 * (1.0 - 1.18721 * (1.0 - math.pow(t / 256.164, 8.0))), 'V'
	}
	return 632.400e6 * (1.0 - 1.07476 * (1.0 - math.pow(t / 273.31, 4.6))), 'VI'
}

// shown is x to 7 significant figures, rounded by math.round (printf's
// ties differ between compilers): plain from 1e-4 to 1e7, else with an
// exponent, and no trailing zeros.
fn shown(x f64) string {
	if math.is_nan(x) {
		return 'NaN'
	}
	if math.is_inf(x, 0) {
		return if x > 0 { 'inf' } else { '-inf' }
	}
	if x == 0.0 {
		return '0'
	}
	a := math.abs(x)
	mut e := int(math.floor(math.log10(a)))
	mut m := mantissa(a, e)
	if m >= 10_000_000 {
		e++
		m = mantissa(a, e)
	} else if m < 1_000_000 {
		e--
		m = mantissa(a, e)
	}
	digits := m.str()
	sign := if x < 0 { '-' } else { '' }
	if e >= 7 || e < -4 {
		frac := digits[1..].trim_right('0')
		return '${sign}${digits[..1]}${if frac == '' { '' } else { '.' + frac }}e${e}'
	}
	if e >= 6 {
		return sign + digits
	}
	if e >= 0 {
		frac := digits[e + 1..].trim_right('0')
		return '${sign}${digits[..e + 1]}${if frac == '' { '' } else { '.' + frac }}'
	}
	return '${sign}0.${'0'.repeat(-e - 1)}${digits.trim_right('0')}'
}

// mantissa is a's seven digits from its 10^e place, as an integer.
fn mantissa(a f64, e int) i64 {
	k := 6 - e
	return i64(math.round(if k >= 0 { a * math.pow(10, k) } else { a / math.pow(10, -k) }))
}

// table is the state as `reflux water` prints it: one phase per kg and per
// mole, or a saturated mixture's two phases and the whole.
fn (st WaterState) table() string {
	mut lines := []string{}
	mut head := 'Water by ${st.formulation}: ${st.phase}'
	if st.phase == 'liquid and vapour' {
		head += ", x = ${shown(st.x)} (the vapour's mass fraction)"
	}
	head += ', at ${shown(st.t)} K and ${shown(st.p * 1e-6)} MPa'
	if st.region > 0 {
		head += ', in its region ${st.region}'
	}
	lines << head + '.'
	lines << "h, s and u on IAPWS's reference: u = s = 0 for the liquid at the triple point."
	lines << ''
	mm := water.molar_mass
	mut columns := []WaterPhase{}
	mut titles := []string{}
	if liquid := st.liquid {
		vapour := st.vapour or { WaterPhase{} }
		titles = ['liquid', 'vapour', 'mixture', 'per mole']
		columns = [liquid, vapour, st.whole]
	} else {
		titles = ['per kg', 'per mole']
		columns = [st.whole]
	}
	// Each column is a number, right-aligned, and its unit.
	mut header := ['']
	for title in titles {
		header << [title, '']
	}
	mut rows := [header]
	rows << row('T', shown(st.t), 'K', titles.len)
	rows << row('P', shown(st.p * 1e-6), 'MPa', titles.len)
	for name in ['rho', 'v', 'h', 's', 'u', 'cv', 'cp', 'w', 'mu', 'lambda'] {
		mut cells := [name]
		for col in columns {
			cells << cell(name, col, 1.0)
		}
		cells << if name in ['w', 'mu', 'lambda'] { ['', ''] } else { cell(name, st.whole, mm) }
		rows << cells
	}
	lines << aligned(rows)
	if !math.is_nan(st.sigma) {
		lines << ''
		lines << 'sigma ${shown(st.sigma * 1e3)} mN/m, the surface tension at ${shown(st.t)} K (R1-76).'
	}
	for note in st.notes {
		lines << 'Note: ${note}'
	}
	return lines.join('\n')
}

// row is a T or P row: its value in the first column and blanks after.
fn row(name string, value string, unit string, columns int) []string {
	mut r := [name, value, unit]
	for _ in 1 .. columns {
		r << ['', '']
	}
	return r
}

// cell is one property of a phase and its unit, per kg (mm = 1) or per
// mole (mm the molar mass, kg/mol); blank where it has none.
fn cell(name string, ph WaterPhase, mm f64) []string {
	molar := mm != 1.0
	x, unit := match name {
		'rho' { ph.rho / mm, if molar { 'mol/m3' } else { 'kg/m3' } }
		'v' { ph.v * mm, if molar { 'm3/mol' } else { 'm3/kg' } }
		'h' { ph.h * mm * 1e-3, if molar { 'kJ/mol' } else { 'kJ/kg' } }
		'u' { ph.u * mm * 1e-3, if molar { 'kJ/mol' } else { 'kJ/kg' } }
		's' {
			ph.s * mm * (if molar { 1.0 } else { 1e-3 }), if molar {
				'J/(mol K)'
			} else {
				'kJ/(kg K)'
			}
		}
		'cv' {
			ph.cv * mm * (if molar { 1.0 } else { 1e-3 }), if molar {
				'J/(mol K)'
			} else {
				'kJ/(kg K)'
			}
		}
		'cp' {
			ph.cp * mm * (if molar { 1.0 } else { 1e-3 }), if molar {
				'J/(mol K)'
			} else {
				'kJ/(kg K)'
			}
		}
		'w' { ph.w, 'm/s' }
		'mu' { ph.mu * 1e6, 'uPa s' }
		else { ph.lambda * 1e3, 'mW/(m K)' }
	}
	return if math.is_finite(x) { [shown(x), unit] } else { ['', ''] }
}

// aligned lays rows out in columns: the names, then each number
// right-aligned and its unit left-aligned after it.
fn aligned(rows [][]string) string {
	mut widths := []int{len: rows[0].len}
	for r in rows {
		for i, c in r {
			widths[i] = math.max(widths[i], c.runes().len)
		}
	}
	mut out := []string{}
	for r in rows {
		mut line := r[0] + ' '.repeat(widths[0] - r[0].runes().len)
		for i in 1 .. r.len {
			pad := ' '.repeat(widths[i] - r[i].runes().len)
			line += if i % 2 == 1 { '   ' + pad + r[i] } else { ' ' + r[i] + pad }
		}
		out << line.trim_right(' ')
	}
	return out.join('\n')
}

// The water report's JSON (docs/report-format.md, "The water report"):
// SI throughout, per kg and per mole, a property it has none of left out.
struct WaterJson {
	format      int
	reflux      string
	formulation string
	reference   string
	given       []GivenJson
	t           f64 @[json: 'T']
	p           f64 @[json: 'P']
	phase       string
	x           ?f64
	region      ?int
	per_kg      PhaseJson
	per_mole    PhaseJson
	liquid      ?PhaseJson
	vapour      ?PhaseJson
	sigma       ?f64
	notes       []string
}

struct GivenJson {
	name  string
	text  string
	value f64
}

struct PhaseJson {
	rho    ?f64
	v      ?f64
	h      ?f64
	s      ?f64
	u      ?f64
	cv     ?f64
	cp     ?f64
	w      ?f64
	mu     ?f64
	lambda ?f64
}

// json_number is x, or none where it isn't finite, which the report leaves
// out. (Not `finite`: module main's functions keep their C names, and
// mingw's math.h declares an int finite(double), so gcc on Windows refuses.)
fn json_number(x f64) ?f64 {
	return if math.is_finite(x) { x } else { none }
}

// per_json is a phase's properties per kg (mm = 1) or per mole (mm the
// molar mass); a mole has no transport properties of its own.
fn per_json(ph WaterPhase, mm f64) PhaseJson {
	molar := mm != 1.0
	return PhaseJson{
		rho:    json_number(ph.rho / mm)
		v:      json_number(ph.v * mm)
		h:      json_number(ph.h * mm)
		s:      json_number(ph.s * mm)
		u:      json_number(ph.u * mm)
		cv:     json_number(ph.cv * mm)
		cp:     json_number(ph.cp * mm)
		w:      if molar { none } else { json_number(ph.w) }
		mu:     if molar { none } else { json_number(ph.mu) }
		lambda: if molar { none } else { json_number(ph.lambda) }
	}
}

fn (st WaterState) json() string {
	mut liquid := ?PhaseJson(none)
	mut vapour := ?PhaseJson(none)
	if l := st.liquid {
		liquid = per_json(l, 1.0)
	}
	if v := st.vapour {
		vapour = per_json(v, 1.0)
	}
	r := WaterJson{
		format:      1
		reflux:      flowsheet.version()
		formulation: st.formulation
		reference:   "IAPWS's: u = s = 0 for the liquid at the triple point"
		given:       st.specs.map(GivenJson{
			name:  it.given.str()
			text:  it.text
			value: it.value
		})
		t:           st.t
		p:           st.p
		phase:       st.phase
		x:           json_number(st.x)
		region:      if st.region > 0 { st.region } else { none }
		per_kg:      per_json(st.whole, 1.0)
		per_mole:    per_json(st.whole, water.molar_mass)
		liquid:      liquid
		vapour:      vapour
		sigma:       json_number(st.sigma)
		notes:       st.notes
	}
	return json2.encode(r, prettify: true, indent_string: '  ') + '\n'
}
