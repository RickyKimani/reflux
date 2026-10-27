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
pub type Alpha = InvSqrtTr | MathiasCopeman | Prsv | Soave | Twu | Unity

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

// Twu is the α of Twu et al. (1991):
//
//     α = Tr^(N(M−1)) · exp(L(1 − Tr^(NM)))
//
// L, M and N are fitted per component. Unlike the Soave form, α stays
// positive and keeps falling at high Tr, which is why it extrapolates to
// supercritical components better.
pub struct Twu {
pub:
	l f64
	m f64
	n f64
}

// MathiasCopeman is the α of Mathias & Copeman (1983), fitted per
// component, usually to polar compounds' vapour pressures:
//
//     α = [1 + c₁u + c₂u² + c₃u³]²,   u = 1 − √Tr
//
// Above Tc the original uses only the c₁ term, because the cubic in u
// misbehaves at high Tr. That is the default here, matching Caleb Bell's
// `thermo`. It keeps α and α′ continuous at Tc but makes α″, and so Cp, jump
// there by c₂. Set `truncated: false` for the smooth form at every Tr.
pub struct MathiasCopeman {
pub:
	c1        f64
	c2        f64
	c3        f64
	truncated bool = true
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

// eval returns α, dα/dTr and d²α/dTr². The formulas are written out in
// docs/design.md under "Second derivatives of the cubic".
pub fn (a Alpha) eval(tr f64) (f64, f64, f64) {
	match a {
		Unity {
			return 1.0, 0.0, 0.0
		}
		InvSqrtTr {
			sq := math.sqrt(tr)
			return 1.0 / sq, -0.5 / (tr * sq), 0.75 / (tr * tr * sq)
		}
		Twu {
			alpha, g1, g2 := a.ln_derivs(tr)
			return alpha, alpha * g1, alpha * (g2 + g1 * g1)
		}
		// One branch for three variants: `a` stays an Alpha here, with no
		// smart cast, which is all a.root needs.
		Soave, Prsv, MathiasCopeman {
			k, dk, d2k := a.root(tr)
			return k * k, 2.0 * k * dk, 2.0 * (dk * dk + k * d2k)
		}
	}
}

// root returns √α, d√α/dTr and d²√α/dTr². cubic.v builds √aᵢ from it, so
// nothing divides by √α.
//
// Soave, PRSV and Mathias–Copeman are written as α = k², and for them √α is
// k itself, signed. Soave's k = 1 + m(1 − √Tr) reaches zero at
// Tr = (1 + 1/m)², about 11 for nitrogen with PR (1,400 K) and 12.6 for
// methane (2,400 K). Through that point k passes smoothly through zero, and
// α = k² touches zero and rises again, as the published form does. Taking
// √α = |k| instead would put a kink in √aᵢ there, and so in a mixture's
// a(T) and a jump in its enthalpy. One consequence past the zero: a cross
// term (1 − kᵢⱼ)·√aᵢ·√aⱼ is negative while only one of the pair has passed
// its zero. Caleb Bell's `thermo` takes |k|, so the two differ there, and
// only there. Neither is physical: Soave's α isn't meant for such Tr, and
// Twu's α, which stays positive and keeps falling, is the one to use.
pub fn (a Alpha) root(tr f64) (f64, f64, f64) {
	match a {
		Unity {
			return 1.0, 0.0, 0.0
		}
		InvSqrtTr {
			// √α = Tr^(−1/4).
			q := 1.0 / math.sqrt(math.sqrt(tr))
			return q, -0.25 * q / tr, 0.3125 * q / (tr * tr)
		}
		Twu {
			// ln √α = ½ ln α, so (√α)′ = √α·½(ln α)′ and
			// (√α)″ = √α·[½(ln α)″ + ¼(ln α)′²].
			alpha, g1, g2 := a.ln_derivs(tr)
			r := math.sqrt(alpha)
			return r, 0.5 * r * g1, r * (0.5 * g2 + 0.25 * g1 * g1)
		}
		Soave {
			// Inside this branch `a` is a Soave, so a.m is accessible:
			// a smart cast, like a Go type switch variable.
			sq := math.sqrt(tr)
			return 1.0 + a.m * (1.0 - sq), -a.m / (2.0 * sq), a.m / (4.0 * tr * sq)
		}
		Prsv {
			sq := math.sqrt(tr)
			// κ = κ₀ + κ₁·h(Tr) and k = 1 + κ·u(Tr).
			h := (1.0 + sq) * (0.7 - tr)
			dh := (0.7 - tr) / (2.0 * sq) - (1.0 + sq)
			d2h := -1.0 / sq - (0.7 - tr) / (4.0 * tr * sq)
			kappa := a.kappa0 + a.kappa1 * h
			dkappa := a.kappa1 * dh
			d2kappa := a.kappa1 * d2h
			u := 1.0 - sq
			du := -1.0 / (2.0 * sq)
			d2u := 1.0 / (4.0 * tr * sq)
			return 1.0 + kappa * u, dkappa * u + kappa * du, d2kappa * u + 2.0 * dkappa * du +
				kappa * d2u
		}
		MathiasCopeman {
			sq := math.sqrt(tr)
			u := 1.0 - sq
			du := -1.0 / (2.0 * sq)
			d2u := 1.0 / (4.0 * tr * sq)
			// Above Tc the truncated form keeps only c₁.
			c2, c3 := if a.truncated && tr >= 1.0 { 0.0, 0.0 } else { a.c2, a.c3 }
			dk_du := a.c1 + u * (2.0 * c2 + 3.0 * c3 * u)
			d2k_du2 := 2.0 * c2 + 6.0 * c3 * u
			return 1.0 + u * (a.c1 + u * (c2 + u * c3)), dk_du * du, d2k_du2 * du * du + dk_du * d2u
		}
	}
}

// ln_derivs returns Twu's α with (ln α)′ and (ln α)″, from
// ln α = N(M−1)·ln Tr + L(1 − Tr^(NM)). eval and root both build on them.
fn (a Twu) ln_derivs(tr f64) (f64, f64, f64) {
	nm := a.n * a.m
	trp := math.pow(tr, nm)
	alpha := math.pow(tr, a.n * (a.m - 1.0)) * math.exp(a.l * (1.0 - trp))
	g1 := a.n * (a.m - 1.0) / tr - a.l * nm * trp / tr
	g2 := -a.n * (a.m - 1.0) / (tr * tr) - a.l * nm * (nm - 1.0) * trp / (tr * tr)
	return alpha, g1, g2
}
