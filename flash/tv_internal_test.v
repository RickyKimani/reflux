module flash

import math
import thermo

fn tv_fluid() !thermo.Fluid {
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

// differences compares an analytic Jacobian with central differences of
// its residuals.
fn differences(eval fn ([]f64, mut []f64, mut []f64) !, v []f64, what string) ! {
	n := v.len
	mut f := []f64{len: n}
	mut jac := []f64{len: n * n}
	eval(v, mut f, mut jac)!
	for col in 0 .. n {
		h := 1e-6 * math.max(1e-3, math.abs(v[col]))
		mut up := v.clone()
		mut down := v.clone()
		up[col] += h
		down[col] -= h
		mut fu := []f64{len: n}
		mut fd := []f64{len: n}
		mut scratch := []f64{len: n * n}
		eval(up, mut fu, mut scratch)!
		eval(down, mut fd, mut scratch)!
		for row in 0 .. n {
			fdiff := (fu[row] - fd[row]) / (2.0 * h)
			err := math.abs(jac[row * n + col] - fdiff)
			assert err <= 1e-5 * math.max(1.0, math.abs(fdiff)), '${what}: J[${row}][${col}] ${jac[row * n + col]} vs ${fdiff}'
		}
	}
}

// split_near is the methane/n-butane split at 250 K and 30 bar, its molar
// volume, and a point a little off it in (nⱽ, Vⱽ).
fn split_near() !(thermo.Fluid, FlashResult, f64, []f64) {
	fl := tv_fluid()!
	r := pt(fl.model, 250.0, 30e5, [0.4, 0.6])!
	mut v := 0.0
	for ph in r.phases {
		v += ph.fraction * ph.state.molar_volume()
	}
	y := r.phases[0].state.n
	x := [0.9 * r.beta * y[0], 1.1 * r.beta * y[1], 0.95 * r.beta * r.phases[0].state.molar_volume()]
	return fl, r, v, x
}

fn test_the_helmholtz_derivatives_match_differences() {
	// helmholtz_eval's gradient against central differences of A/RT, and
	// its Hessian against central differences of the gradient, a little off
	// the TV solution.
	fl, _, v, x := split_near()!
	s := volume_flash(fl.model, 250.0, v, [0.4, 0.6], VolumeParams{})!
	mut g := []f64{len: 3}
	mut h := []f64{len: 9}
	s.helmholtz_eval(x, mut g, mut h)!
	for k in 0 .. 3 {
		step := 1e-6 * math.abs(x[k])
		mut up := x.clone()
		mut down := x.clone()
		up[k] += step
		down[k] -= step
		mut scratch := []f64{len: 3}
		mut hs := []f64{len: 9}
		au := s.helmholtz_eval(up, mut scratch, mut hs)!
		ad := s.helmholtz_eval(down, mut scratch, mut hs)!
		fdiff := (au - ad) / (2.0 * step)
		assert math.abs(g[k] - fdiff) <= 1e-5 * math.max(1.0, math.abs(fdiff)), 'g[${k}] ${g[k]} vs ${fdiff}'
	}
	differences(fn [s] (x []f64, mut f []f64, mut jac []f64) ! {
		s.helmholtz_eval(x, mut f, mut jac)!
	}, x, 'Helmholtz')!
}

fn test_the_uv_jacobian_matches_differences() {
	fl, _, v, x := split_near()!
	e := EnergyFlash{
		fl:  fl
		u:   -15000.0
		tvs: volume_flash(fl.model, 250.0, v, [0.4, 0.6], VolumeParams{})!
	}
	mut xt := x.clone()
	xt << math.log(255.0)
	differences(fn [e] (x []f64, mut f []f64, mut jac []f64) ! {
		e.eval(x, mut f, mut jac)!
	}, xt, 'UV')!
}

fn test_the_one_phase_uv_jacobian_matches_differences() {
	// The vapour's volume at 350 K and 10 bar, evaluated at 340 K against a
	// U it doesn't have: r_U is far from 0 there, so the Jacobian's −r_U
	// term shows. Newton converges without it, only more slowly, which is
	// why a planted mistake there survived every other test.
	fl := tv_fluid()!
	r := pt(fl.model, 350.0, 10e5, [0.4, 0.6])!
	e := EnergyFlash{
		fl:  fl
		u:   -15000.0
		tvs: volume_flash(fl.model, 350.0, r.phases[0].state.molar_volume(), [0.4, 0.6],
			VolumeParams{})!
	}
	differences(fn [e] (x []f64, mut f []f64, mut jac []f64) ! {
		SingleU{e}.eval(x, mut f, mut jac)!
	}, [math.log(340.0)], 'one-phase UV')!
}
