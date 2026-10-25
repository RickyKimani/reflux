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
	for _ in 0 .. sat_max_iter {
		// With no liquid root, P is below the liquid branch's minimum: too
		// low. With no vapor root, P is above the vapor branch's maximum.
		vl := liquid_volume(m, t, p, n) or {
			p *= 1.1
			continue
		}
		vv := vapor_volume(m, t, p, n) or {
			p *= 0.9
			continue
		}
		if math.abs(vl - vv) <= 1e-8 * vv {
			return error('${c.name}: lost the two-phase region at T = ${t} K (too close to Tc?)')
		}
		diff := gibbs_residual(m, t, p, vl, n) - gibbs_residual(m, t, p, vv, n)
		p *= math.exp(diff)
		if math.abs(diff) < sat_tol {
			return p
		}
	}
	return error('${c.name}: saturation pressure at ${t} K did not converge in ${sat_max_iter} iterations')
}
