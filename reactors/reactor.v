module reactors

import flowsheet
import kinetics
import math
import numeric
import thermo

// What every reactor shares: its errors, the elements it conserves, its
// outlet's state, and the temperature an energy balance sets.

// spec_error is a reactor's specification it can't use.
fn spec_error(unit string, reason string) IError {
	return flowsheet.FlowsheetError{
		kind:   .bad_specification
		unit:   unit
		reason: reason
	}
}

// failed is a reactor's equilibrium, integration or Newton failing.
fn failed(unit string, reason string) IError {
	return flowsheet.FlowsheetError{
		kind:   .reaction_failed
		unit:   unit
		reason: reason
	}
}

// named gives an error from flowsheet's own functions, which name no unit,
// the reactor's name.
fn named(unit string, err IError) IError {
	if err is flowsheet.FlowsheetError {
		e := err as flowsheet.FlowsheetError
		return flowsheet.FlowsheetError{
			kind:   e.kind
			unit:   unit
			reason: e.reason
			cause:  e.cause
		}
	}
	return failed(unit, err.msg())
}

// element_rows is what any reaction conserves: a row per element of the
// species, and a row of its own for a species without a formula, which
// can't react. Every reactor reports these, so a flowsheet with reactors is
// checked on its elements.
pub fn element_rows(species []kinetics.Species) [][]f64 {
	mut zs := []int{}
	for sp in species {
		f := sp.formula or { continue }
		for c in f.counts {
			if c.z !in zs {
				zs << c.z
			}
		}
	}
	zs.sort()
	mut rows := [][]f64{}
	for z in zs {
		mut row := []f64{len: species.len}
		for i, sp in species {
			if f := sp.formula {
				row[i] = f64(f.count(z))
			}
		}
		rows << row
	}
	for i, sp in species {
		if sp.formula == none {
			mut row := []f64{len: species.len}
			row[i] = 1
			rows << row
		}
	}
	return rows
}

// Thermal is how a reactor's temperature is set.
pub enum Thermal {
	isothermal // at a given T, or a PFR's inlet's: the duty is what that takes
	adiabatic  // no heat: the outlet's enthalpy is the inlet's
	duty       // a given heat, W, into the reactor
}

// at_tp is the outlet stream with these flows at T and P.
fn at_tp(unit string, fl thermo.Fluid, flows []f64, t f64, p f64) !flowsheet.Stream {
	return flowsheet.feed_tp(fl, flows, t, p) or { return named(unit, err) }
}

// enthalpy_flow is a stream's H·F, W.
fn enthalpy_flow(s flowsheet.Stream) f64 {
	return s.total() * s.h
}

// Outlet is a reactor's outlet at a temperature: what it makes there, at
// its outlet's pressure. An interface, not a closure: V never frees one.
interface Outlet {
	outlet_at(t f64) !flowsheet.Stream
}

// Balance is H_out(T) − H_target, W, for Brent.
struct Balance {
	o      Outlet
	target f64 // W
}

fn (b Balance) at(t f64) !f64 {
	s := b.o.outlet_at(t)!
	return s.total() * s.h - b.target
}

// solve_t finds the outlet whose enthalpy flow is target, W: Brent on
// H_out(T) − target, bracketed outwards from a guess, between 0.5 and
// 10,000 K. H rises with T, so the root is above the guess where H falls
// short.
fn solve_t(unit string, o Outlet, target f64, guess f64) !flowsheet.Stream {
	b := Balance{o, target}
	g0 := b.at(guess) or { return failed(unit, 'no outlet at ${guess} K: ${err.msg()}') }
	if g0 == 0 {
		return o.outlet_at(guess)!
	}
	up := g0 < 0
	mut lo, mut hi := guess, guess
	mut step := 50.0
	mut found := false
	for _ in 0 .. 40 {
		next := if up { math.min(1e4, hi + step) } else { math.max(0.5, lo - step) }
		g := b.at(next) or { return failed(unit, 'no outlet at ${next} K: ${err.msg()}') }
		if up {
			lo, hi = hi, next
		} else {
			lo, hi = next, lo
		}
		if (g > 0) == up || g == 0 {
			found = true
			break
		}
		if (!up && next <= 0.5) || (up && next >= 1e4) {
			break
		}
		step *= 2
	}
	if !found {
		return failed(unit, 'no temperature from 0.5 to 10,000 K balances its energy')
	}
	r := numeric.brent(b, lo, hi, xtol: 1e-9) or {
		return failed(unit, 'its energy balance did not converge: ${err.msg()}')
	}
	return o.outlet_at(r.x)!
}

// one_phase refuses an outlet that would split: a reactor here works in one
// phase, a second being M6's.
fn one_phase(unit string, s flowsheet.Stream) ! {
	if s.result.phases.len > 1 {
		return failed(unit, 'its outlet would split into two phases at ${s.t} K and ${s.p} Pa, which needs M6')
	}
}

// pressure_out is a reactor's outlet pressure: the inlet's less a drop.
fn pressure_out(unit string, inlet flowsheet.Stream, drop f64) !f64 {
	p := inlet.p - drop
	if !(p > 0) {
		return spec_error(unit, 'a drop of ${drop} Pa from ${inlet.p} Pa leaves no pressure')
	}
	return p
}
