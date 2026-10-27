module thermo

import numeric

// PhaseLabel is what a phase is called, for unit operations and for people.
pub enum PhaseLabel {
	liquid
	vapor
}

// pip returns the phase identification parameter of Venkatarathnam and
// Oellrich (2011) at (T, V, n):
//
//     Π = V·[(∂²P/∂T∂V)/(∂P/∂T) − (∂²P/∂V²)/(∂P/∂V)]
//
// It needs no saturation curve or critical point, so it labels a phase in
// any mixture, supercritical ones included: Π > 1 is liquid-like and Π < 1
// vapour-like (phase_label). Caleb Bell's `thermo` labels phases this way.
//
// ∂P/∂V and ∂P/∂T are analytic. The second derivatives would need third
// derivatives of F, which the Model contract doesn't have (D6), so they're
// central differences of the first ones: in V for ∂²P/∂V², and in T for
// ∂²P/∂T∂V. With steps of 1e-5 of T and of V − V_min, Π is good to about
// 1e-8 relative, plenty for a label.
pub fn pip(m Model, t f64, v f64, n []f64) f64 {
	d := pressure_derivs(m, t, v, n)
	hv := 1e-5 * (v - m.min_volume(n))
	ht := 1e-5 * t
	d2p_dv2 := (dp_dv_at(m, t, v + hv, n) - dp_dv_at(m, t, v - hv, n)) / (2.0 * hv)
	d2p_dtdv := (dp_dv_at(m, t + ht, v, n) - dp_dv_at(m, t - ht, v, n)) / (2.0 * ht)
	return v * (d2p_dtdv / d.dp_dt - d2p_dv2 / d.dp_dv)
}

// phase_label calls a single phase liquid when Π > 1 and vapour otherwise,
// as `thermo` does. In a two-phase split the flash labels by density
// instead, the less dense phase being the vapour (D6), since Π can call
// both phases liquid-like.
pub fn phase_label(m Model, st State) PhaseLabel {
	return if pip(m, st.t, st.v, st.n) > 1.0 { .liquid } else { .vapor }
}

// dp_dv_at is (∂P/∂V) at (T, V, n), from the volume solver's cheap derivatives.
fn dp_dv_at(m Model, t f64, v f64, n []f64) f64 {
	d := m.volume_derivs(t, v, n)
	return -r_gas * t * (d.f_vv + numeric.sum(n) / (v * v))
}
