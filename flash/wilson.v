module flash

import math
import thermo

// wilson_k returns Wilson's (1968) estimate of the K-values,
//
//     Kᵢ = (Pcᵢ/P)·exp(5.373(1 + ωᵢ)(1 − Tcᵢ/T)),
//
// from each component's critical constants alone. It's the usual first guess
// for stability trials and for a flash.
pub fn wilson_k(comps []thermo.Component, t f64, p f64) []f64 {
	return comps.map(it.pc / p * math.exp(5.373 * (1.0 + it.omega) * (1.0 - it.tc / t)))
}
