// The phase seam (M6's D3, criterion 4): a Helmholtz model behind it gives
// the bits of the functions on a Model, at every root, and Fluid is on it.
// An external test, through the pub API only.
import math
import thermo { Component, Root }

const propane = Component{
	name:  'propane'
	tc:    369.8
	pc:    42.48e5
	omega: 0.152
	mw:    0.044097
}

const hexane = Component{
	name:  'n-hexane'
	tc:    507.6
	pc:    30.25e5
	omega: 0.301
	mw:    0.086177
}

fn same(a []f64, b []f64) bool {
	return a.map(math.f64_bits(it)) == b.map(math.f64_bits(it))
}

fn ideal() []thermo.IdealGas {
	r := thermo.r_gas
	return [9.0, 17.0].map(thermo.IdealGas{
		cp: thermo.Polynomial{
			terms: [thermo.PowerTerm{r * it, 0}]
		}
	})
}

fn models() ![]thermo.Model {
	pr := thermo.new_cubic(.pr, [propane, hexane], kij: [0.0, 0.01, 0.01, 0.0])!
	srk := thermo.new_cubic(.srk, [propane, hexane])!
	shifted := thermo.translate(pr, [3e-6, 6e-6])!
	return [thermo.Model(pr), srk, shifted]
}

fn test_volume_roots_give_the_models_bits() {
	mut compared, mut two_roots, mut refused := 0, 0, 0
	for m in models()! {
		seam := thermo.by_roots(m)
		assert seam.components() == m.components()
		fl := thermo.new_fluid(m, ideal())!
		for t in [220.0, 300.0, 330.0, 360.0, 420.0, 700.0] {
			for p in [1e3, 2e5, 8e5, 12e5, 60e5, 400e5] {
				for n in [[1.0, 0.0], [0.7, 0.3], [0.2, 2.8]] {
					mut volumes := []f64{}
					for root in [Root.liquid, .vapor, .stable] {
						want := thermo.solve_tp(m, t, p, n, root) or {
							if _ := seam.phase(t, p, n, root) {
								assert false, 'a phase where solve_tp has none'
							}
							refused++
							continue
						}
						at := '${t} K, ${p} Pa, ${n}, ${root}'
						st := seam.phase(t, p, n, root)!
						assert same([st.t, st.p, st.v, st.z], [want.t, want.p, want.v, want.z]), at
						assert same(st.n, want.n) && same(st.ln_phi, want.ln_phi), at
						volumes << st.v

						d, dw := seam.ln_phi_derivs(st)!, thermo.ln_phi_derivs(m, want)
						assert same(d.dt, dw.dt) && same(d.dp, dw.dp) && same(d.dn, dw.dn), at
						assert seam.phase_label(st) == thermo.phase_label(m, want), at

						res, rw := seam.phase_residual(st), thermo.residual_props(m, want)
						assert same([res.h, res.s, res.cp, res.cv], [rw.h, rw.s, rw.cp, rw.cv]), at
						pd := thermo.pressure_derivs(m, t, want.v, n)
						assert same([res.dv_dt, res.dp_dv], [-pd.dp_dt / pd.dp_dv, pd.dp_dv]), at

						// And Fluid, on the seam, adds them to the ideal gas's: at
						// H's and Cv's own arithmetic.
						pr := fl.props(st)
						mut h_ig, mut cp_ig, mut total := 0.0, 0.0, 0.0
						for i, ig in fl.ideal {
							total += n[i]
							if n[i] == 0.0 {
								continue
							}
							h_ig += n[i] * (0.0 + ig.cp.int_cp(thermo.t_ref, t))
							cp_ig += n[i] * ig.cp.cp(t)
						}
						assert same([pr.h, pr.cp, pr.cv], [h_ig + rw.h, cp_ig + rw.cp,
							cp_ig - total * thermo.r_gas + rw.cv]), at
						compared++
					}
					if volumes.len == 3 && volumes[1] > 1.5 * volumes[0] {
						two_roots++
					}
				}
			}
		}
	}
	println('${compared} phases compared, ${two_roots} states with two roots, ${refused} refused')
	// Both roots are among them, or "the root asked for" wouldn't be tested.
	assert compared > 500 && two_roots > 20
}

// A bad state is refused through the seam as solve_tp refuses it, with the
// same typed error.
fn test_the_seam_passes_refusals_on() {
	seam := thermo.by_roots(models()![0])
	if _ := seam.phase(-1.0, 1e5, [0.5, 0.5], .stable) {
		assert false, 'a phase at −1 K'
	} else {
		assert err is thermo.InputError, err.msg()
	}
}

fn test_a_fluid_gives_back_its_equation_of_state() {
	for m in models()! {
		fl := thermo.new_fluid(m, ideal())!
		back := fl.helmholtz()!
		assert back is thermo.Cubic == (m is thermo.Cubic)
		a, b := back.residual(300.0, 2e-4, [0.5, 0.5]), m.residual(300.0, 2e-4, [0.5, 0.5])
		assert same([a.f, a.f_v, a.f_t], [b.f, b.f_v, b.f_t]) && same(a.f_n, b.f_n)
	}
}

// NoF is a phase model with no Helmholtz energy to give, whatever is inside.
struct NoF {
	inner thermo.PhaseModel
}

fn (m NoF) components() []Component {
	return m.inner.components()
}

fn (m NoF) phase(t f64, p f64, n []f64, root Root) !thermo.State {
	return m.inner.phase(t, p, n, root)
}

fn (m NoF) ln_phi_derivs(st thermo.State) !thermo.LnPhiDerivs {
	return m.inner.ln_phi_derivs(st)
}

fn (m NoF) phase_label(st thermo.State) thermo.PhaseLabel {
	return m.inner.phase_label(st)
}

fn (m NoF) phase_residual(st thermo.State) thermo.PhaseResidual {
	return m.inner.phase_residual(st)
}

fn test_a_fluid_without_one_says_so() {
	seam := thermo.PhaseModel(NoF{thermo.by_roots(models()![0])})
	fl := thermo.new_phase_fluid(seam, ideal())!
	if _ := fl.helmholtz() {
		assert false, 'an F(T, V, n) from a model with none'
	} else {
		assert err is thermo.NoHelmholtzError, err.msg()
		assert err.msg().contains('F(T, V, n)') && err.msg().contains('T and P alone'), err.msg()
	}
	if _ := thermo.new_phase_fluid(seam, ideal()[..1]) {
		assert false, 'one ideal gas for two components'
	} else {
		assert err.msg().contains('1 components'), err.msg()
	}
}
