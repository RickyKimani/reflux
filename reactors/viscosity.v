// Module reactors is M9's reactors (D10): conversion, equilibrium over
// chosen reactions, Gibbs, CSTR and PFR, each a flowsheet unit with a maker
// for the file (D1), and the gas viscosity a packed bed's pressure drop
// needs.
module reactors

import components
import math

// Viscosity is a gas mixture's viscosity, Pa·s, at T (K), P (Pa) and mole
// fractions y, for Ergun's pressure drop (D10, as amended in step 0). The
// PFR takes it as an interface, so a test can give it DWSIM's.
pub interface Viscosity {
	viscosity(t f64, p f64, y []f64) !f64
}

// Wilke is a gas mixture's viscosity at low pressure: each compound's from
// its component file's vapour-viscosity correlation (ChemSep's equation 102,
// say), mixed by Wilke's rule,
//
//     μ = Σᵢ yᵢμᵢ / Σⱼ yⱼΦᵢⱼ,
//     Φᵢⱼ = [1 + (μᵢ/μⱼ)^½ (Mⱼ/Mᵢ)^¼]² / [8(1 + Mᵢ/Mⱼ)]^½.
//
// A mole-fraction mean of the μᵢ misses how a little heavy gas raises a
// light one's viscosity: in a reforming gas, mostly hydrogen, by 30% (step
// 0). Pressure is ignored, the gas being near ideal; a correlation is used
// outside its fitted range, as every property is (M2's D7).
pub struct Wilke {
	pure []components.Property
	mw   []f64 // kg/mol
}

// new_wilke takes each compound's vapour-viscosity correlation and molar
// mass from its component file.
pub fn new_wilke(records []components.Record) !Wilke {
	mut pure := []components.Property{}
	mut mw := []f64{}
	for r in records {
		pure << r.property('vapor_viscosity')!
		m := r.constants.mw or { return error('${r.name} has no molar mass') }
		if !(m > 0) {
			return error("${r.name}'s molar mass, ${m} kg/mol, isn't one")
		}
		mw << m
	}
	return Wilke{pure, mw}
}

// viscosity is the mixture's, Pa·s.
pub fn (w Wilke) viscosity(t f64, _ f64, y []f64) !f64 {
	n := w.pure.len
	if y.len != n {
		return error('${y.len} mole fractions for ${n} compounds')
	}
	mut mu := []f64{len: n}
	for i in 0 .. n {
		mu[i] = w.pure[i].at(t).value
		if !(mu[i] > 0) || math.is_inf(mu[i], 1) {
			return error('the viscosity of compound ${i} at ${t} K is ${mu[i]} Pa·s')
		}
	}
	mut out := 0.0
	for i in 0 .. n {
		if y[i] == 0 {
			continue
		}
		mut den := 0.0
		for j in 0 .. n {
			a := 1 + math.sqrt(mu[i] / mu[j]) * math.pow(w.mw[j] / w.mw[i], 0.25)
			den += y[j] * a * a / math.sqrt(8 * (1 + w.mw[i] / w.mw[j]))
		}
		out += y[i] * mu[i] / den
	}
	return out
}
