module flash

import math
import thermo

fn methane_butane_fluid() !thermo.Fluid {
	m := thermo.new_cubic(.pr, [
		thermo.Component{
			name:  'methane'
			tc:    190.6
			pc:    45.99e5
			omega: 0.012
			mw:    0.016043
		},
		thermo.Component{
			name:  'n-butane'
			tc:    425.1
			pc:    37.96e5
			omega: 0.200
			mw:    0.058123
		},
	], kij: [0.0, 0.02, 0.02, 0.0])!
	// Smith, Van Ness & Abbott's Cp/R for the two, times R.
	cp := [[1.702, 9.081e-3, -2.164e-6], [1.935, 36.915e-3, -11.402e-6]]
	return thermo.new_fluid(m, cp.map(thermo.IdealGas{
		cp: thermo.Polynomial{
			terms: [
				thermo.PowerTerm{thermo.r_gas * it[0], 0},
				thermo.PowerTerm{thermo.r_gas * it[1], 1},
				thermo.PowerTerm{thermo.r_gas * it[2], 2},
			]
		}
	}))!
}

fn test_the_q_function_jacobian_matches_differences() {
	// eval's Jacobian, Q's Hessian, against central differences of its
	// residuals, for PH and PS, at a two-phase point a little off the
	// solution.
	fl := methane_butane_fluid()!
	z := [0.4, 0.6]
	r := pt(fl.model, 250.0, 30e5, z)!
	assert r.phases.len == 2
	for spec in [Spec.h, .s] {
		s := StateFlash{
			fl:      fl
			spec:    spec
			value:   if spec == .h { -15000.0 } else { -60.0 }
			p:       30e5
			z:       z
			present: [0, 1]
			params:  StateParams{}
		}
		y := r.phases[0].state.n
		mut v := [0.9 * r.beta * y[0], 1.1 * r.beta * y[1], math.log(255.0)]
		n := v.len
		mut f := []f64{len: n}
		mut jac := []f64{len: n * n}
		s.eval(v, mut f, mut jac)!
		for col in 0 .. n {
			h := 1e-6 * math.max(1.0, math.abs(v[col]))
			mut up := v.clone()
			mut down := v.clone()
			up[col] += h
			down[col] -= h
			mut fu := []f64{len: n}
			mut fd := []f64{len: n}
			mut scratch := []f64{len: n * n}
			s.eval(up, mut fu, mut scratch)!
			s.eval(down, mut fd, mut scratch)!
			for row in 0 .. n {
				fdiff := (fu[row] - fd[row]) / (2.0 * h)
				err := math.abs(jac[row * n + col] - fdiff)
				assert err <= 1e-6 * math.max(1.0, math.abs(fdiff)), '${spec}: J[${row}][${col}] ${jac[row * n + col]} vs ${fdiff}'
			}
		}
	}
}

fn test_the_ideal_gas_start_aims_at_the_ideal_gas() {
	// ideal_t searches for where the feed as an ideal gas has the specified
	// H or S. From the state's own H and S, that gap at its T is minus its
	// residual H or S (61 J/mol and 0.10 J/(mol·K) at 1 bar and 400 K). A
	// wrong ideal mixing term, 2R·Σz ln z = 11 J/(mol·K) here, would show
	// at once. Only the start depends on it, which is why it needs its own
	// test: the bracket finds the answer from anywhere.
	fl := methane_butane_fluid()!
	z := [0.4, 0.6]
	r := pt(fl.model, 400.0, 1e5, z)!
	assert r.phases.len == 1
	st := r.phases[0].state
	pr := fl.props(st)
	res := thermo.residual_props(fl.model, st)
	for spec in [Spec.h, .s] {
		s := StateFlash{
			fl:      fl
			spec:    spec
			value:   if spec == .h { pr.h } else { pr.s }
			p:       1e5
			z:       z
			present: [0, 1]
			params:  StateParams{}
		}
		// Fluid.props computes the ideal gas on its own, so this checks
		// IdealGap's formula against it: measured, they agree to 5e-13.
		gap := IdealGap{s}.at(400.0)
		own := if spec == .h { -res.h } else { -res.s }
		assert math.abs(gap - own) <= 1e-9 * math.max(1.0, math.abs(s.value)), '${spec}: ${gap} vs ${own}'
	}
}
