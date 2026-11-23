module thermo

// PhaseModel is the flash's seam (M6's D3): whatever gives a phase at (T, P,
// n). flash asks a model for nothing else, so it works alike on an equation
// of state, by its volume roots (VolumeRoots), and on a γ–φ fluid, whose
// liquid has no F(T, V, n) at all.
//
// A model holds nothing about a flash, so threads may share one. A State
// handed back to it is one it gave.
pub interface PhaseModel {
	// components lists the components, in the order used by n: their
	// critical constants for a flash's starts, their molar masses for a
	// phase's mass density, with the State's volume.
	components() []Component
	// phase is the phase at (T, P, n) that root asks for, the liquid-like,
	// the vapour-like or the one of lower G, with its volume and ln φᵢ.
	// Where a model has only one of them, that's the one it gives, whatever
	// was asked for.
	phase(t f64, p f64, n []f64, root Root) !State
	// ln_phi_derivs are a phase's derivatives of ln φᵢ in T, P and n. An
	// error says why the model has none there.
	ln_phi_derivs(st State) !LnPhiDerivs
	// phase_label is what a single phase is called.
	phase_label(st State) PhaseLabel
	// phase_residual is a phase's residual properties, for Fluid.props to
	// add to the ideal gas's. It can't fail, as props can't: what the model
	// can't give is NaN, which shows in what's made of it, as a Cp
	// correlation's does (ideal.v).
	phase_residual(st State) PhaseResidual
}

// Scoped is a phase model that can keep what depends on T alone between
// calls, in a workspace for one thread: the γ–φ fluid's Pˢᵃᵗ, Vᴸ and its
// activity model's parameters (M6's D4). A model holds nothing about a
// flash, so threads share one, and each flash asks for its own copy.
pub interface Scoped {
	scoped() PhaseModel
}

// for_one_thread is m with a workspace of its own, if it keeps one, or m:
// what a flash calls first, so that its many phases at one T share what
// depends on T alone.
pub fn for_one_thread(m PhaseModel) PhaseModel {
	if m is Scoped {
		return m.scoped()
	}
	return m
}

// PhaseResidual is what Fluid.props needs of a phase beyond the ideal gas:
// its residual properties on the (T, P) basis, extensive as in
// ResidualProps, and its volume's two derivatives, for the speed of sound
// and the Joule–Thomson coefficient.
pub struct PhaseResidual {
pub:
	h     f64 // Hʳ, J
	s     f64 // Sʳ, J/K
	cp    f64 // Cpʳ, J/K
	cv    f64 // Cvʳ, J/K
	dv_dt f64 // (∂V/∂T)_P,n, m³/K
	dp_dv f64 // (∂P/∂V)_T,n, Pa/m³
}

// VolumeRoots is a Helmholtz model as a PhaseModel, its phases being its
// volume roots. Each method is the function of the same name on a Model
// (solve_tp for phase), and gives its bits.
pub struct VolumeRoots {
pub:
	model Model
}

// by_roots puts a Helmholtz model behind the seam, for flash:
//
//     r := flash.pt(thermo.by_roots(eos), t, p, z)!
//
// A model only has to give F(T, V, n), so a new equation of state needs
// nothing more than this.
pub fn by_roots(m Model) PhaseModel {
	return VolumeRoots{
		model: m
	}
}

// components implements PhaseModel.
pub fn (r VolumeRoots) components() []Component {
	return r.model.components()
}

// phase implements PhaseModel: solve_tp's root.
pub fn (r VolumeRoots) phase(t f64, p f64, n []f64, root Root) !State {
	return solve_tp(r.model, t, p, n, root)
}

// ln_phi_derivs implements PhaseModel, from F's second derivatives, which
// always exist.
pub fn (r VolumeRoots) ln_phi_derivs(st State) !LnPhiDerivs {
	return ln_phi_derivs(r.model, st)
}

// phase_label implements PhaseModel, by Π and the four-covolume rule.
pub fn (r VolumeRoots) phase_label(st State) PhaseLabel {
	return phase_label(r.model, st)
}

// phase_residual implements PhaseModel, from one residual2.
pub fn (r VolumeRoots) phase_residual(st State) PhaseResidual {
	r2 := r.model.residual2(st.t, st.v, st.n)
	pd := pressure_derivs_from(r2, st.t, st.v, st.n)
	res := residual_props_from(r2, pd, st)
	return PhaseResidual{
		h:     res.h
		s:     res.s
		cp:    res.cp
		cv:    res.cv
		dv_dt: -pd.dp_dt / pd.dp_dv
		dp_dv: pd.dp_dv
	}
}

// NoHelmholtzError refuses work that needs F(T, V, n) on a fluid that
// hasn't one, as a γ–φ fluid hasn't: the TV and UV flashes and the reactors.
pub struct NoHelmholtzError {
	Error
}

// msg says what's missing.
pub fn (e NoHelmholtzError) msg() string {
	return "it needs an equation of state's Helmholtz energy, F(T, V, n), and the fluid's model has none: it gives its phases at T and P alone"
}

// helmholtz is the equation of state behind a fluid's phases, for what
// works in F(T, V, n) itself, or a NoHelmholtzError where there's none.
pub fn (f Fluid) helmholtz() !Model {
	m := f.model
	if m is VolumeRoots {
		return m.model
	}
	return NoHelmholtzError{}
}
