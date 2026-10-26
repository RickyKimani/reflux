module thermo

import math
import numeric

// Everything in this file is written against Model, not Cubic, so any model
// that supplies residual2 (CPA, PC-SAFT, a reference equation) gets these
// properties for free. The formulas are tabulated in docs/design.md under
// "Properties from F".

// PressureDerivs are P and its derivatives at one (T, V, n).
pub struct PressureDerivs {
pub:
	p     f64   // Pa
	dp_dv f64   // (∂P/∂V)_T,n, Pa/m³
	dp_dt f64   // (∂P/∂T)_V,n, Pa/K
	dp_dn []f64 // (∂P/∂nᵢ)_T,V, Pa/mol
	v_bar []f64 // partial molar volumes V̄ᵢ = −(∂P/∂nᵢ)/(∂P/∂V), m³/mol
}

// pressure_derivs evaluates P and its derivatives at (T, V, n).
pub fn pressure_derivs(m Model, t f64, v f64, n []f64) PressureDerivs {
	return pressure_derivs_from(m.residual2(t, v, n), t, v, n)
}

fn pressure_derivs_from(r Residual2, t f64, v f64, n []f64) PressureDerivs {
	rt := r_gas * t
	p := rt * (numeric.sum(n) / v - r.f_v)
	dp_dv := -rt * (r.f_vv + numeric.sum(n) / (v * v))
	dp_dn := r.f_vn.map(rt * (1.0 / v - it))
	return PressureDerivs{
		p:     p
		dp_dv: dp_dv
		dp_dt: p / t - rt * r.f_tv
		dp_dn: dp_dn
		v_bar: dp_dn.map(-it / dp_dv)
	}
}

// ResidualProps are a State's residual properties on the (T, P) basis:
// Mʳ = M(T, P, n) − Mⁱᵍ(T, P, n). They are extensive, in J or J/K for the
// moles in the State; divide by Σnᵢ for molar values.
pub struct ResidualProps {
pub:
	h  f64 // Hʳ, J
	s  f64 // Sʳ, J/K
	g  f64 // Gʳ, J
	u  f64 // Uʳ, J
	a  f64 // Aʳ, J
	cv f64 // Cvʳ, J/K
	cp f64 // Cpʳ, J/K
}

// residual_props returns the residual properties of a solved State.
pub fn residual_props(m Model, st State) ResidualProps {
	r := m.residual2(st.t, st.v, st.n)
	return residual_props_from(r, pressure_derivs_from(r, st.t, st.v, st.n), st)
}

// residual_props_from does the work of residual_props from an already
// evaluated residual, so that Fluid.props needs only one residual2 call.
fn residual_props_from(r Residual2, pd PressureDerivs, st State) ResidualProps {
	t, v := st.t, st.v
	ntot := numeric.sum(st.n)
	rt := r_gas * t
	h := -rt * t * r.f_t + st.p * v - ntot * rt
	s := -rt * r.f_t - r_gas * r.f + ntot * r_gas * math.log(st.z)
	u := -rt * t * r.f_t
	cv := -r_gas * t * (t * r.f_tt + 2.0 * r.f_t)
	return ResidualProps{
		h:  h
		s:  s
		g:  h - t * s
		u:  u
		a:  u - t * s
		cv: cv
		// Cp − Cv = −T(∂P/∂T)²/(∂P/∂V) for the real fluid and nR for the
		// ideal gas, so the residual parts differ by the difference.
		cp: cv - t * pd.dp_dt * pd.dp_dt / pd.dp_dv - ntot * r_gas
	}
}

// LnPhiDerivs are the derivatives of ln φᵢ that a Newton flash needs.
pub struct LnPhiDerivs {
pub:
	dt []f64 // (∂ln φᵢ/∂T)_P,n, 1/K
	dp []f64 // (∂ln φᵢ/∂P)_T,n, 1/Pa
	// dn is n·(∂ln φᵢ/∂nⱼ)_T,P, nc×nc row-major, dimensionless. It is
	// symmetric, and Gibbs–Duhem makes Σᵢ nᵢ·dn[i][j] = 0 for every j.
	dn []f64
}

// ln_phi_derivs returns the T, P and composition derivatives of ln φᵢ at a
// solved State (Michelsen & Mollerup, chapter 2).
pub fn ln_phi_derivs(m Model, st State) LnPhiDerivs {
	t := st.t
	nc := st.n.len
	r := m.residual2(t, st.v, st.n)
	pd := pressure_derivs_from(r, t, st.v, st.n)
	ntot := numeric.sum(st.n)
	rt := r_gas * t
	mut dt := []f64{len: nc}
	mut dp := []f64{len: nc}
	mut dn := []f64{len: nc * nc}
	for i in 0 .. nc {
		dt[i] = r.f_tn[i] + 1.0 / t - pd.v_bar[i] * pd.dp_dt / rt
		dp[i] = pd.v_bar[i] / rt - 1.0 / st.p
		for j in 0 .. nc {
			dn[i * nc + j] = ntot * r.f_nn[i * nc + j] + 1.0 +
				ntot * pd.dp_dn[i] * pd.dp_dn[j] / (rt * pd.dp_dv)
		}
	}
	return LnPhiDerivs{
		dt: dt
		dp: dp
		dn: dn
	}
}
