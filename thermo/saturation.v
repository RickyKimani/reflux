module thermo

import math

const sat_max_iter = 500
const sat_tol = 1e-10

// saturation_pressure finds Psat(T) for a single-component model: the
// pressure at which the liquid and vapor roots have equal fugacity. It works
// for any Model, not only cubics. It uses successive substitution,
// P ← P·φL/φV, which converges steadily but slows down as T approaches Tc.
pub fn saturation_pressure(m Model, t f64) !f64 {
	comps := m.components()
	if comps.len != 1 {
		return error('saturation_pressure needs a single-component model, got ${comps.len} components')
	}
	c := comps[0]
	c.check()!
	if t <= 0.0 || t >= c.tc {
		return error('${c.name}: T = ${t} K is not between 0 and Tc = ${c.tc} K, so there is no saturation pressure')
	}
	n := [1.0]
	// Start from the definition of the acentric factor (Edmister):
	//   log10(Psat/Pc) ≈ (7/3)(1 + ω)(1 − Tc/T)
	mut p := c.pc * math.pow(10.0, 7.0 / 3.0 * (1.0 + c.omega) * (1.0 - c.tc / t))
	fx := fix(m, t, n)
	// lo is the last pressure with no liquid root, hi the last with no
	// vapour root. Once there are both, the window where both roots exist
	// lies between them, and the next try is their geometric mean (M5a's D7):
	// near Tc the window is narrower than a step of 1.1 or 0.9 (0.07% of P
	// for IAPWS-95 at 646 K), which those steps bounce across. A path that
	// never meets both is the one it always was.
	mut lo := 0.0
	mut hi := 0.0
	for _ in 0 .. sat_max_iter {
		// With no liquid root, P is below the liquid branch's minimum: too
		// low. With no vapor root, P is above the vapor branch's maximum.
		vl := liquid_volume(fx, p) or {
			lo = p
			p = if hi > 0.0 { math.sqrt(lo * hi) } else { p * 1.1 }
			continue
		}
		vv := vapor_volume(fx, p) or {
			hi = p
			p = if lo > 0.0 { math.sqrt(lo * hi) } else { p * 0.9 }
			continue
		}
		if math.abs(vl - vv) <= 1e-8 * vv {
			return error('${c.name}: lost the two-phase region at T = ${t} K (too close to Tc?)')
		}
		diff := gibbs_residual(fx, p, vl) - gibbs_residual(fx, p, vv)
		p *= math.exp(diff)
		if math.abs(diff) < sat_tol {
			return p
		}
	}
	return error('${c.name}: saturation pressure at ${t} K did not converge in ${sat_max_iter} iterations')
}
