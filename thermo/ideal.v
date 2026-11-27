module thermo

import math
import tdep

// Ideal-gas heat capacity forms. Each one gives Cp(T) and the two integrals
// that H and S need, ∫Cp dT and ∫Cp/T dT. Everything is per mole, in
// J/(mol·K) and J/mol.
pub type CpForm = AlyLee | Correlated | PlanckEinstein | Polynomial

// PowerTerm is c·T^p.
pub struct PowerTerm {
pub:
	c f64
	p int
}

// Polynomial is Cp = Σ cₖ·T^pₖ. Negative powers are allowed, so the common
// A + B·T + C·T² + D·T⁻² form is four terms.
pub struct Polynomial {
pub:
	terms []PowerTerm
}

// AlyLee is DIPPR equation 107 (Aly & Lee, 1981), the form ChemSep and DIPPR
// use for most ideal-gas heat capacities:
//
//     Cp = A + B·[(C/T)/sinh(C/T)]² + D·[(E/T)/cosh(E/T)]²
//
// A, B and D in J/(mol·K); C and E in K.
//
// It's evaluated by tdep's eq. 107, which has the closed-form integrals and
// stays finite at low T, where C/T is large enough that sinh and cosh
// overflow a float.
pub struct AlyLee {
pub:
	a f64
	b f64
	c f64
	d f64
	e f64
}

fn (f AlyLee) correlation() tdep.Correlation {
	return tdep.Correlation{
		eq:   107
		a:    f.a
		b:    f.b
		c:    f.c
		d:    f.d
		e:    f.e
		tmax: math.inf(1)
	}
}

// Correlated is Cp from one of tdep's temperature-dependent correlations,
// such as ChemSep's eq. 16, the form of most of its ideal-gas heat
// capacities. The correlation keeps its source's units, and scale converts
// its value to J/(mol·K): 1e-3 for ChemSep's J/(kmol·K).
//
// Forms with no closed-form integral are integrated numerically. If that
// fails (it takes a correlation that overflows somewhere between the two
// temperatures), int_cp and int_cp_over_t return NaN. Fluid.props stays
// infallible that way, and the NaN shows in H and S instead of hiding as a
// wrong number.
pub struct Correlated {
pub:
	c     tdep.Correlation @[required]
	scale f64 = 1.0
}

// EinsteinTerm is one Planck–Einstein term of a PlanckEinstein Cp: c in
// J/(mol·K), θ in K.
pub struct EinsteinTerm {
pub:
	c     f64
	theta f64
}

// PlanckEinstein is a constant plus Planck–Einstein terms, each the heat
// capacity of one harmonic vibration of characteristic temperature θ:
//
//     Cp = a + Σ cᵢ·E(θᵢ/T),   E(u) = u²·eᵘ/(eᵘ − 1)²
//
// a in J/(mol·K). It's the form of the multiparameter reference equations'
// ideal parts: IAPWS-95's five terms (M5a's D5), and GERG-2008's and Span
// and Wagner's. Aly–Lee's sinh term is one such term: (x/sinh x)² = E(2x),
// so B·[(C/T)/sinh(C/T)]² is c = B, θ = 2C.
//
// Both integrals are closed forms. With u = θ/T,
//
//     ∫E dT   = θ/(eᵘ − 1)
//     ∫E/T dT = u/(eᵘ − 1) − ln(1 − e⁻ᵘ)
//
// written with expm1 and log1p, so that they keep their digits where u is
// small (a hot gas) and stay finite where it's large (u = 2·10⁴ at 1 K for
// IAPWS-95's largest θ, where eᵘ overflows a float).
pub struct PlanckEinstein {
pub:
	a     f64
	terms []EinsteinTerm
}

// einstein is E(u) = u²eᵘ/(eᵘ − 1)², written as u²e⁻ᵘ/(1 − e⁻ᵘ)², which
// underflows to 0 instead of overflowing where u is large.
fn einstein(u f64) f64 {
	e := math.exp(-u)
	d := -math.expm1(-u) // 1 − e⁻ᵘ, with its digits where u is small
	return u * u * e / (d * d)
}

// h_anti is an antiderivative of Cp: a·T + Σ cᵢθᵢ/(e^(θᵢ/T) − 1).
fn (f PlanckEinstein) h_anti(t f64) f64 {
	mut s := f.a * t
	for term in f.terms {
		s += term.c * term.theta / math.expm1(term.theta / t)
	}
	return s
}

// s_anti is an antiderivative of Cp/T: a·ln T + Σ cᵢ[u/(eᵘ − 1) − ln(1 − e⁻ᵘ)].
fn (f PlanckEinstein) s_anti(t f64) f64 {
	mut s := f.a * math.log(t)
	for term in f.terms {
		u := term.theta / t
		s += term.c * (u / math.expm1(u) - math.log1p(-math.exp(-u)))
	}
	return s
}

// cp returns the ideal-gas heat capacity at T, J/(mol·K).
pub fn (f CpForm) cp(t f64) f64 {
	return match f {
		Polynomial {
			mut s := 0.0
			for term in f.terms {
				s += term.c * math.pow(t, f64(term.p))
			}
			s
		}
		AlyLee {
			f.correlation().eval(t)
		}
		Correlated {
			f.scale * f.c.eval(t)
		}
		PlanckEinstein {
			mut s := f.a
			for term in f.terms {
				s += term.c * einstein(term.theta / t)
			}
			s
		}
	}
}

// int_cp returns ∫Cp dT from t1 to t2, J/mol.
pub fn (f CpForm) int_cp(t1 f64, t2 f64) f64 {
	return match f {
		Polynomial { f.h_anti(t2) - f.h_anti(t1) }
		AlyLee { f.correlation().integral(t1, t2) or { math.nan() } }
		Correlated { f.scale * (f.c.integral(t1, t2) or { math.nan() }) }
		PlanckEinstein { f.h_anti(t2) - f.h_anti(t1) }
	}
}

// int_cp_over_t returns ∫Cp/T dT from t1 to t2, J/(mol·K).
pub fn (f CpForm) int_cp_over_t(t1 f64, t2 f64) f64 {
	return match f {
		Polynomial { f.s_anti(t2) - f.s_anti(t1) }
		AlyLee { f.correlation().integral_over_t(t1, t2) or { math.nan() } }
		Correlated { f.scale * (f.c.integral_over_t(t1, t2) or { math.nan() }) }
		PlanckEinstein { f.s_anti(t2) - f.s_anti(t1) }
	}
}

// h_anti is an antiderivative of Cp.
fn (f Polynomial) h_anti(t f64) f64 {
	mut s := 0.0
	for term in f.terms {
		// ∫T^p dT = T^(p+1)/(p+1), except p = −1, which gives ln T.
		s += if term.p == -1 {
			term.c * math.log(t)
		} else {
			term.c * math.pow(t, f64(term.p + 1)) / f64(term.p + 1)
		}
	}
	return s
}

// s_anti is an antiderivative of Cp/T.
fn (f Polynomial) s_anti(t f64) f64 {
	mut s := 0.0
	for term in f.terms {
		// ∫T^(p−1) dT = T^p/p, except p = 0, which gives ln T.
		s += if term.p == 0 {
			term.c * math.log(t)
		} else {
			term.c * math.pow(t, f64(term.p)) / f64(term.p)
		}
	}
	return s
}

// IdealGas is one component's ideal-gas data. The formation values put H
// and S on the reference state that reflux uses everywhere: elements in
// their standard states at 298.15 K and P° = 1 bar. With them, reactors
// need no separate heat-of-reaction bookkeeping.
//
// They are Options (`?f64`) because data sets often lack them. Without them
// H, S, G, U and A are still correct *differences*, but not absolute, and
// Props.absolute says so.
pub struct IdealGas {
pub:
	cp  CpForm @[required]
	dhf ?f64 // ΔH°f at 298.15 K, J/mol
	s0  ?f64 // S° at 298.15 K and 1 bar, J/(mol·K)
}

// Reference temperature and pressure for formation properties.
pub const t_ref = 298.15

pub const p_ref = 1e5
