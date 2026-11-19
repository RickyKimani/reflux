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
	mut k := w.room()
	w.at(t, mut k)!
	return w.mixed(k, y)
}

// WilkeAt is what Wilke's rule needs at a T, whatever the composition: each
// compound's viscosity, and for each pair the two factors of Φᵢⱼ, its
// numerator's root 1 + (μᵢ/μⱼ)^½ (Mⱼ/Mᵢ)^¼ and its denominator
// [8(1 + Mᵢ/Mⱼ)]^½. A reactor keeps one, so that a T that hasn't moved
// costs it no correlation, power or square root (M8's D11): an isothermal
// reactor's is made once, and mixed at each composition along it.
pub struct WilkeAt {
pub mut:
	kept bool  // whether it holds a T's values, at t
	t    f64   // K
	mu   []f64 // Pa·s
	a    []f64 // n×n, row-major
	root []f64 // n×n, row-major
}

// room is a WilkeAt of this mixture's size, holding nothing yet.
pub fn (w Wilke) room() WilkeAt {
	n := w.pure.len
	return WilkeAt{
		mu:   []f64{len: n}
		a:    []f64{len: n * n}
		root: []f64{len: n * n}
	}
}

// at fills k for T, unless it holds that T's values already. A compound
// whose viscosity there isn't one is refused, and then k holds nothing.
pub fn (w Wilke) at(t f64, mut k WilkeAt) ! {
	n := w.pure.len
	if k.mu.len != n || k.a.len != n * n || k.root.len != n * n {
		return error('room for ${k.mu.len} viscosities for ${n} compounds')
	}
	if k.kept && t == k.t {
		return
	}
	k.kept = false
	for i in 0 .. n {
		k.mu[i] = w.pure[i].at(t).value
		if !(k.mu[i] > 0) || math.is_inf(k.mu[i], 1) {
			return error('the viscosity of compound ${i} at ${t} K is ${k.mu[i]} Pa·s')
		}
	}
	for i in 0 .. n {
		for j in 0 .. n {
			k.a[i * n + j] = 1 + math.sqrt(k.mu[i] / k.mu[j]) * math.pow(w.mw[j] / w.mw[i], 0.25)
			k.root[i * n + j] = math.sqrt(8 * (1 + w.mw[i] / w.mw[j]))
		}
	}
	k.t = t
	k.kept = true
}

// mixed is the mixture's viscosity, Pa·s, at mole fractions y, from what
// at filled, by Wilke's rule: the arithmetic the rule is written in above,
// in its order, so that it's the same number however k came to be filled.
pub fn (w Wilke) mixed(k &WilkeAt, y []f64) f64 {
	n := w.pure.len
	mut out := 0.0
	for i in 0 .. n {
		if y[i] == 0 {
			continue
		}
		mut den := 0.0
		for j in 0 .. n {
			a := k.a[i * n + j]
			den += y[j] * a * a / k.root[i * n + j]
		}
		out += y[i] * k.mu[i] / den
	}
	return out
}
