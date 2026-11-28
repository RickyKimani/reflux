module transport

import math

// Wilke is a gas mixture's viscosity at low pressure (M9's D10, moved here
// by M5b's D2 with its arithmetic): each compound's from its vapour-viscosity
// correlation (ChemSep's equation 102, say), mixed by Wilke's rule,
//
//     μ = Σᵢ yᵢμᵢ / Σⱼ yⱼΦᵢⱼ,
//     Φᵢⱼ = [1 + (μᵢ/μⱼ)^½ (Mⱼ/Mᵢ)^¼]² / [8(1 + Mᵢ/Mⱼ)]^½.
//
// A mole-fraction mean of the μᵢ misses how a little heavy gas raises a
// light one's viscosity: in a reforming gas, mostly hydrogen, by 30% (M9
// step 0). Pressure is ignored; a phase's viscosity takes Lucas's factor
// on it (D5), and a PFR's Ergun drop takes it as it is (D11). A correlation
// is used outside its fitted range, as every property is (M2's D7).
//
// Its Φᵢⱼ are also Mason and Saxena's Aᵢⱼ at ε = 1, so the same factors
// give a gas's conductivity (conductivity, D4).
pub struct Wilke {
pub:
	pure []Correlated
	mw   []f64 // kg/mol
}

// new_wilke takes each compound's vapour-viscosity correlation, and its
// molar mass, which must be positive.
pub fn new_wilke(pure []Correlated, mw []f64) !Wilke {
	if pure.len != mw.len {
		return error('${pure.len} viscosities for ${mw.len} molar masses')
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
		k.mu[i] = w.pure[i].at(t)
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

// conductivity is the gas's thermal conductivity at low pressure, W/(m·K),
// by Wassiljewa's form with Mason and Saxena's Aᵢⱼ at ε = 1 (D4),
//
//     λ = Σᵢ yᵢλᵢ / Σⱼ yⱼAᵢⱼ,  Aᵢⱼ = Φᵢⱼ,
//
// Wilke's Φᵢⱼ of the viscosities at that T, which at filled; lambda holds
// each compound's conductivity there.
pub fn (w Wilke) conductivity(k &WilkeAt, y []f64, lambda []f64) f64 {
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
		out += y[i] * lambda[i] / den
	}
	return out
}
