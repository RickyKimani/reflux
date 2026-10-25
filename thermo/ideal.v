module thermo

import math

// Ideal-gas heat capacity forms. Each one gives Cp(T) and the two integrals
// that H and S need, ∫Cp dT and ∫Cp/T dT, in closed form. Everything is
// per mole, in J/(mol·K) and J/mol. Database units (DIPPR uses J/(kmol·K))
// are converted when the data is loaded, never here.
pub type CpForm = AlyLee | Polynomial

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
pub struct AlyLee {
pub:
	a f64
	b f64
	c f64
	d f64
	e f64
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
			mut s := f.a
			if f.b != 0.0 {
				s += f.b * x_over_sinh_sq(f.c / t)
			}
			if f.d != 0.0 {
				s += f.d * x_over_cosh_sq(f.e / t)
			}
			s
		}
	}
}

// int_cp returns ∫Cp dT from t1 to t2, J/mol.
pub fn (f CpForm) int_cp(t1 f64, t2 f64) f64 {
	return f.h_anti(t2) - f.h_anti(t1)
}

// int_cp_over_t returns ∫Cp/T dT from t1 to t2, J/(mol·K).
pub fn (f CpForm) int_cp_over_t(t1 f64, t2 f64) f64 {
	return f.s_anti(t2) - f.s_anti(t1)
}

// h_anti is an antiderivative of Cp.
fn (f CpForm) h_anti(t f64) f64 {
	match f {
		Polynomial {
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
		AlyLee {
			// d/dT[C·coth(C/T)] = (C/T)²/sinh²(C/T), and
			// d/dT[−E·tanh(E/T)] = (E/T)²/cosh²(E/T).
			mut s := f.a * t
			if f.b != 0.0 {
				s += f.b * f.c / math.tanh(f.c / t)
			}
			if f.d != 0.0 {
				s -= f.d * f.e * math.tanh(f.e / t)
			}
			return s
		}
	}
}

// s_anti is an antiderivative of Cp/T.
fn (f CpForm) s_anti(t f64) f64 {
	match f {
		Polynomial {
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
		AlyLee {
			// With x = C/T: d/dT[x·coth x − ln sinh x] = x²/(T·sinh²x), and
			// with y = E/T: d/dT[−(y·tanh y − ln cosh y)] = y²/(T·cosh²y).
			mut s := f.a * math.log(t)
			if f.b != 0.0 {
				x := f.c / t
				s += f.b * (x / math.tanh(x) - ln_sinh(x))
			}
			if f.d != 0.0 {
				y := f.e / t
				s -= f.d * (y * math.tanh(y) - ln_cosh(y))
			}
			return s
		}
	}
}

// The Aly–Lee terms overflow if written naively: sinh and cosh exceed f64
// for arguments above about 710, which C/T reaches at low T. These forms
// rewrite them in terms of e^(−x), which only underflows harmlessly to 0.

// x_over_sinh_sq returns (x/sinh x)² for x ≥ 0.
fn x_over_sinh_sq(x f64) f64 {
	if x == 0.0 {
		return 1.0
	}
	if x < 1.0 {
		r := x / math.sinh(x)
		return r * r
	}
	e := math.exp(-x)
	r := 2.0 * x * e / (1.0 - e * e)
	return r * r
}

// x_over_cosh_sq returns (x/cosh x)² for x ≥ 0.
fn x_over_cosh_sq(x f64) f64 {
	e := math.exp(-x)
	r := 2.0 * x * e / (1.0 + e * e)
	return r * r
}

// ln_sinh returns ln(sinh x) for x > 0.
fn ln_sinh(x f64) f64 {
	if x < 1.0 {
		return math.log(math.sinh(x))
	}
	return x - math.ln2 + math.log1p(-math.exp(-2.0 * x))
}

// ln_cosh returns ln(cosh x) for x ≥ 0. For small x, cosh x is so close to
// 1 that ln(cosh x) loses relative precision; cosh x − 1 = 2·sinh²(x/2)
// keeps it.
fn ln_cosh(x f64) f64 {
	if x < 1.0 {
		s := math.sinh(0.5 * x)
		return math.log1p(2.0 * s * s)
	}
	return x - math.ln2 + math.log1p(math.exp(-2.0 * x))
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
