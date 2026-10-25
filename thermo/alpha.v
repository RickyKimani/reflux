module thermo

import math

// Alpha is the temperature dependence of a cubic's attraction parameter,
// a(T) = a_c·α(Tr). It is chosen per component, because the parameters of
// PRSV, Twu or Mathias–Copeman belong to the component, not to the
// equation.
//
// This is a sum type: the module owns the whole list, and `match` must
// handle every variant. Adding one (say `Twu{l f64, m f64, n f64}`) makes
// the compiler point at every `match` that needs a new branch.
pub type Alpha = InvSqrtTr | Prsv | Soave | Unity

// Unity is α = 1 (van der Waals).
pub struct Unity {}

// InvSqrtTr is α = 1/√Tr (the original Redlich–Kwong).
pub struct InvSqrtTr {}

// Soave is α = [1 + m(1 − √Tr)]², used by SRK and PR with m(ω).
pub struct Soave {
pub:
	m f64
}

// Prsv is Stryjek & Vera's (1986) α for Peng–Robinson:
//
//     α = [1 + κ(1 − √Tr)]²,   κ = κ₀(ω) + κ₁(1 + √Tr)(0.7 − Tr)
pub struct Prsv {
pub:
	kappa0 f64
	kappa1 f64
}

// soave_srk is the SRK (1972) m(ω) correlation.
pub fn soave_srk(omega f64) Alpha {
	return Soave{
		m: 0.480 + 1.574 * omega - 0.176 * omega * omega
	}
}

// soave_pr is the PR (1976) m(ω) correlation.
pub fn soave_pr(omega f64) Alpha {
	return Soave{
		m: 0.37464 + 1.54226 * omega - 0.26992 * omega * omega
	}
}

// prsv builds the PRSV α from the acentric factor and the fitted κ₁.
//
// Stryjek & Vera recommend κ₁ = 0 above Tr = 0.7. That is not done here.
// κ itself stays continuous at Tr = 0.7, but dκ/dT would jump, and since
// Hʳ depends on dα/dT the enthalpy would jump too. A step in H makes
// energy balances and PH flashes misbehave.
pub fn prsv(omega f64, kappa1 f64) Alpha {
	return Prsv{
		kappa0: 0.378893 + omega * (1.4897153 + omega * (-0.17131848 + omega * 0.0196554))
		kappa1: kappa1
	}
}

// eval returns α and dα/dTr.
pub fn (a Alpha) eval(tr f64) (f64, f64) {
	match a {
		Unity {
			return 1.0, 0.0
		}
		InvSqrtTr {
			sq := math.sqrt(tr)
			return 1.0 / sq, -0.5 / (tr * sq)
		}
		Soave {
			// Inside this branch `a` is a Soave, so a.m is accessible:
			// a smart cast, like a Go type switch variable.
			sq := math.sqrt(tr)
			k := 1.0 + a.m * (1.0 - sq)
			return k * k, -a.m * k / sq
		}
		Prsv {
			sq := math.sqrt(tr)
			kappa := a.kappa0 + a.kappa1 * (1.0 + sq) * (0.7 - tr)
			dkappa := a.kappa1 * ((0.7 - tr) / (2.0 * sq) - (1.0 + sq))
			k := 1.0 + kappa * (1.0 - sq)
			dk := dkappa * (1.0 - sq) - kappa / (2.0 * sq)
			return k * k, 2.0 * k * dk
		}
	}
}
