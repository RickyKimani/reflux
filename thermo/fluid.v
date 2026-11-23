module thermo

import math
import numeric

// Fluid is what a stream needs: a PhaseModel for the phases and their
// non-ideality (M6's D3), plus ideal-gas data for each component. Keeping
// them apart means a model stays purely residual, so CPA or SAFT never deal
// with heat capacities. What needs F(T, V, n) itself asks helmholtz for it.
pub struct Fluid {
pub:
	model PhaseModel
	ideal []IdealGas
}

// new_fluid pairs an equation of state, whose phases are its volume roots,
// with ideal-gas data, one entry per component, in the model's component
// order.
pub fn new_fluid(model Model, ideal []IdealGas) !Fluid {
	return new_phase_fluid(by_roots(model), ideal)
}

// new_phase_fluid is new_fluid for any phase model, with or without an
// F(T, V, n) behind it.
pub fn new_phase_fluid(model PhaseModel, ideal []IdealGas) !Fluid {
	nc := model.components().len
	if ideal.len != nc {
		return error('got ideal-gas data for ${ideal.len} components, but the model has ${nc}')
	}
	return Fluid{
		model: model
		ideal: ideal.clone()
	}
}

// Props are the total properties of a State. They are extensive (J or J/K
// for the moles in the State) and use the reference state described at
// IdealGas.
pub struct Props {
pub:
	h     f64 // J
	s     f64 // J/K
	g     f64 // J
	u     f64 // J
	a     f64 // J
	cp    f64 // J/K
	cv    f64 // J/K
	w     f64 // speed of sound, m/s
	mu_jt f64 // Joule–Thomson coefficient (∂T/∂P)_H, K/Pa
	// absolute is false when some component lacks ΔH°f or S°. H, S, G, U and
	// A then sit on an arbitrary zero. Differences between states of the
	// same composition are still correct.
	absolute bool
}

// props returns the total properties of a solved State: ideal gas at the
// same T, P and n, plus the residual part from the model.
pub fn (f Fluid) props(st State) Props {
	t, p, v := st.t, st.p, st.v
	res := f.model.phase_residual(st)
	ntot := numeric.sum(st.n)
	comps := f.model.components()

	mut h_ig := 0.0
	mut s_ig := 0.0
	mut cp_ig := 0.0
	mut mass := 0.0
	mut absolute := true
	for i, ig in f.ideal {
		ni := st.n[i]
		if ni == 0.0 {
			continue
		}
		// An Option field is read with `or`, like a Result.
		dhf := ig.dhf or {
			absolute = false
			0.0
		}
		s0 := ig.s0 or {
			absolute = false
			0.0
		}
		h_ig += ni * (dhf + ig.cp.int_cp(t_ref, t))
		// −R ln(P/P°) for the pressure, −R ln xᵢ for ideal mixing.
		s_ig += ni * (s0 + ig.cp.int_cp_over_t(t_ref, t) - r_gas * math.log(p / p_ref) - r_gas * math.log(ni / ntot))
		cp_ig += ni * ig.cp.cp(t)
		mass += ni * comps[i].mw
	}

	h := h_ig + res.h
	s := s_ig + res.s
	cp := cp_ig + res.cp
	cv := cp_ig - ntot * r_gas + res.cv
	u := h - p * v
	return Props{
		h:        h
		s:        s
		g:        h - t * s
		u:        u
		a:        u - t * s
		cp:       cp
		cv:       cv
		// w² = (∂P/∂ρ)_S = (Cp/Cv)·(∂P/∂ρ)_T, with ρ the mass density m/V.
		w:        math.sqrt(-v * v / mass * cp / cv * res.dp_dv)
		mu_jt:    (t * res.dv_dt - v) / cp
		absolute: absolute
	}
}
