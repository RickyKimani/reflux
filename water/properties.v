module water

import math
import thermo

// reference_offsets are D6's constants for a fluid made by `fluid`: its h
// and s on reflux's reference less IAPWS's, per mole, J/mol and J/(mol·K).
// D4 and D5 make the total Helmholtz energy IAPWS-95's but for a term linear
// in T, so they're the same at every state (iapws95_test.v holds them to
// 1e-13 over Table 7's). They're taken at Table 6's state, 500 K and
// 838.025 kg/m³.
pub fn reference_offsets(fl thermo.Fluid) (f64, f64) {
	t, rho := 500.0, 838.025 / molar_mass
	pr := fl.props(state_at_density(iapws95(), t, rho))
	ia := props95(t, rho)
	return pr.h - ia.h, pr.s - ia.s
}

// state_at_density is the thermo.State of one mole at (T, ρ), ρ in mol/m³.
fn state_at_density(m thermo.Model, t f64, rho f64) thermo.State {
	v := 1.0 / rho
	p := thermo.pressure(m, t, v, [1.0])
	return thermo.State{
		t:      t
		p:      p
		n:      [1.0]
		v:      v
		z:      p * v / (thermo.r_gas * t)
		ln_phi: [0.0]
	}
}

// Props95 are IAPWS-95's properties at one (T, ρ) on IAPWS's own reference
// (u = s = 0 for the liquid at the triple point), per mole, from the
// release's Table 3. They're the steam-table side of D6: reflux's fluid
// gives the same state's properties on its own reference, which differ by
// the constants of reference_offset.
pub struct Props95 {
pub:
	p  f64 // Pa
	u  f64 // J/mol
	h  f64 // J/mol
	s  f64 // J/(mol·K)
	cv f64 // J/(mol·K)
	cp f64 // J/(mol·K)
	w  f64 // m/s
}

// props95 evaluates Table 3 at T (K) and the molar density rho (mol/m³).
pub fn props95(t f64, rho f64) Props95 {
	tau := critical_t / t
	delta := rho / critical_rho
	o := phi_ideal(tau, delta)
	r := phi_residual(tau, delta)
	rt := r95 * t
	dpd := 1.0 + 2.0 * delta * r.d + delta * delta * r.dd // (∂p/∂ρ)_T/RT
	num := 1.0 + delta * r.d - delta * tau * r.dt
	cv := -r95 * tau * tau * (o.tt + r.tt)
	return Props95{
		p:  rho * rt * (1.0 + delta * r.d)
		u:  rt * tau * (o.t + r.t)
		h:  rt * (1.0 + tau * (o.t + r.t) + delta * r.d)
		s:  r95 * (tau * (o.t + r.t) - o.f - r.f)
		cv: cv
		cp: cv + r95 * num * num / dpd
		// w² = RT(…)/M per kg: Table 3 with the specific gas constant.
		w:  math.sqrt(rt / molar_mass * (dpd - num * num / (tau * tau * (o.tt + r.tt))))
	}
}
