module water

import math
import thermo
import transport

// Transport95 is a stream of IAPWS-95 water's transport (M5b's D10): at
// each phase's T and density, R12-08's viscosity and R15-11's conductivity
// with their critical enhancements, and R1-76's surface tension for a
// liquid, M5a's functions; outside their ranges, none, with why. It's a
// transport.Model, so a flowsheet asks it as it asks a mixture's.
pub struct Transport95 {}

// phase implements transport.Model: one compound, water, at molar volume v
// (m³/mol).
pub fn (_ Transport95) phase(label thermo.PhaseLabel, t f64, p f64, x []f64, v f64) transport.PhaseTransport {
	if x.len != 1 {
		why := '${x.len} mole fractions for water alone'
		return transport.PhaseTransport{
			viscosity:       transport.Result{
				why: why
			}
			conductivity:    transport.Result{
				why: why
			}
			surface_tension: transport.Result{
				why: why
			}
		}
	}
	rho := if v > 0 && !math.is_inf(v, 1) { molar_mass / v } else { math.nan() }
	mu := if m := viscosity(t, rho) {
		transport.Result{
			value: m
		}
	} else {
		transport.Result{
			why: err.msg()
		}
	}
	lambda := if l := thermal_conductivity(t, rho) {
		transport.Result{
			value: l
		}
	} else {
		transport.Result{
			why: err.msg()
		}
	}
	mut sigma := transport.Result{
		why: 'a gas has no surface tension'
	}
	if label != .vapor {
		sigma = if s := surface_tension(t) {
			transport.Result{
				value: s
			}
		} else {
			transport.Result{
				why: err.msg()
			}
		}
	}
	return transport.PhaseTransport{
		viscosity:       mu
		conductivity:    lambda
		surface_tension: sigma
	}
}

// methods implements transport.Model: M5a's functions (D10).
pub fn (_ Transport95) methods() []transport.MethodUsed {
	return [
		transport.MethodUsed{'viscosity', 'any', "IAPWS R12-08 with its critical enhancement, at the phase's T and IAPWS-95 density"},
		transport.MethodUsed{'thermal_conductivity', 'any', "IAPWS R15-11 with its critical enhancement, at the phase's T and IAPWS-95 density"},
		transport.MethodUsed{'surface_tension', 'liquid', 'IAPWS R1-76(2014), at T'},
	]
}
