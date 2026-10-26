module thermo

import math
import tdep

// Ideal-gas heat capacity forms. Each one gives Cp(T) and the two integrals
// that H and S need, ∫Cp dT and ∫Cp/T dT. Everything is per mole, in
// J/(mol·K) and J/mol.
pub type CpForm = AlyLee | Correlated | Polynomial

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
	}
}

// int_cp returns ∫Cp dT from t1 to t2, J/mol.
pub fn (f CpForm) int_cp(t1 f64, t2 f64) f64 {
	return match f {
		Polynomial { f.h_anti(t2) - f.h_anti(t1) }
		AlyLee { f.correlation().integral(t1, t2) or { math.nan() } }
		Correlated { f.scale * (f.c.integral(t1, t2) or { math.nan() }) }
	}
}

// int_cp_over_t returns ∫Cp/T dT from t1 to t2, J/(mol·K).
pub fn (f CpForm) int_cp_over_t(t1 f64, t2 f64) f64 {
	return match f {
		Polynomial { f.s_anti(t2) - f.s_anti(t1) }
		AlyLee { f.correlation().integral_over_t(t1, t2) or { math.nan() } }
		Correlated { f.scale * (f.c.integral_over_t(t1, t2) or { math.nan() }) }
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
