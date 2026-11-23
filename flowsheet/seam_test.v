// A flowsheet needs nothing of its fluid but the phase seam (M6's D3): with
// the cubic hidden behind a model that has no F(T, V, n) to give, every
// golden flowsheet's report is the same text, every number's shortest form.
import flowsheet
import os
import thermo

struct Hidden {
	inner thermo.PhaseModel
}

fn (m Hidden) components() []thermo.Component {
	return m.inner.components()
}

fn (m Hidden) phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State {
	return m.inner.phase(t, p, n, root)
}

fn (m Hidden) ln_phi_derivs(st thermo.State) !thermo.LnPhiDerivs {
	return m.inner.ln_phi_derivs(st)
}

fn (m Hidden) phase_label(st thermo.State) thermo.PhaseLabel {
	return m.inner.phase_label(st)
}

fn (m Hidden) phase_residual(st thermo.State) thermo.PhaseResidual {
	return m.inner.phase_residual(st)
}

fn test_the_golden_flowsheets_need_only_the_seam() {
	for name in ['reference', 'every_unit', 'no_recycle', 'interlocking', 'saturated_loop'] {
		l := flowsheet.load(os.join_path(@VMODROOT, 'testdata', 'flowsheets', '${name}.toml'))!
		hidden := flowsheet.Loaded{
			...l
			fluid: thermo.new_phase_fluid(Hidden{l.fluid.model}, l.fluid.ideal)!
		}
		if _ := hidden.fluid.helmholtz() {
			assert false, "the model isn't hidden"
		}
		want := l.report(l.solve()!).json()
		got := hidden.report(hidden.solve()!).json()
		assert got == want, name
	}
}
