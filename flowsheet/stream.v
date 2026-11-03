// Module flowsheet solves steady-state processes: streams, unit operations,
// and, as M4 goes on, the graph, the solver and the file format.
module flowsheet

import flash
import math
import numeric
import thermo

// Stream is a material stream (D2): each component's molar flow, and the
// flash that gave its state. A stream is always flashed, so it can't be half
// specified, and its H, S, V and phases are always at hand.
pub struct Stream {
pub:
	flows []f64 // mol/s of each component
	t     f64   // K
	p     f64   // Pa
	beta  f64   // vapour fraction
	h     f64   // molar enthalpy, J/mol, on the Fluid's reference state
	s     f64   // molar entropy, J/(mol·K)
	v     f64   // molar volume, m³/mol
	// result is the flash that produced it, per mole of stream. Its z is the
	// composition, which an empty stream keeps though its flows are all 0.
	result flash.FlashResult
}

// total is the stream's molar flow, mol/s.
pub fn (s Stream) total() f64 {
	return numeric.sum(s.flows)
}

// of_flash is the stream that a flash produced, with the given flows.
fn of_flash(fl thermo.Fluid, flows []f64, r flash.FlashResult) Stream {
	mut h, mut s, mut v := 0.0, 0.0, 0.0
	for phase in r.phases {
		pr := fl.props(phase.state)
		h += phase.fraction * pr.h
		s += phase.fraction * pr.s
		v += phase.fraction * phase.state.v
	}
	return Stream{
		flows:  flows.clone()
		t:      r.t
		p:      r.p
		beta:   r.beta
		h:      h
		s:      s
		v:      v
		result: r
	}
}

// check_flows checks a feed's flows. The flash would refuse them too, but
// as mole fractions.
fn check_flows(fl thermo.Fluid, flows []f64) ! {
	nc := fl.model.components().len
	if flows.len != nc {
		return bad_spec('', 'got ${flows.len} flows for ${nc} components')
	}
	for i, n in flows {
		if !math.is_finite(n) || n < 0.0 {
			return bad_spec('', 'the flow of component ${i} is ${n} mol/s')
		}
	}
	total := numeric.sum(flows)
	if !(total > 0.0) || !math.is_finite(total) {
		return bad_spec('', 'a feed needs a finite total flow above zero, not ${total} mol/s')
	}
}

// feed_tp is a feed at T (K) and P (Pa), with its flows in mol/s.
pub fn feed_tp(fl thermo.Fluid, flows []f64, t f64, p f64) !Stream {
	check_flows(fl, flows)!
	r := flash.pt(fl.model, t, p, flows) or { return flash_failed('', err) }
	return of_flash(fl, flows, r)
}

// feed_ph is a feed at P (Pa) and molar enthalpy H (J/mol).
pub fn feed_ph(fl thermo.Fluid, flows []f64, p f64, h f64) !Stream {
	check_flows(fl, flows)!
	r := flash.ph(fl, p, h, flows) or { return flash_failed('', err) }
	return of_flash(fl, flows, r)
}

// feed_tvf is a feed at T (K) and vapour fraction β.
pub fn feed_tvf(fl thermo.Fluid, flows []f64, t f64, beta f64) !Stream {
	check_flows(fl, flows)!
	r := flash.tvf(fl.model, t, beta, flows) or { return flash_failed('', err) }
	return of_flash(fl, flows, r)
}

// feed_pvf is a feed at P (Pa) and vapour fraction β.
pub fn feed_pvf(fl thermo.Fluid, flows []f64, p f64, beta f64) !Stream {
	check_flows(fl, flows)!
	r := flash.pvf(fl.model, p, beta, flows) or { return flash_failed('', err) }
	return of_flash(fl, flows, r)
}
