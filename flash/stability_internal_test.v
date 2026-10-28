module flash

import math
import thermo

// The second-order finish is only as good as its gradient and Hessian.
// These check them against central differences of tm* itself, away from a
// stationary point, where every term of the Hessian matters.

fn feed_at(m thermo.Model, t f64, p f64, z []f64) !Feed {
	st := thermo.solve_tp(m, t, p, z, .stable)!
	present := []int{len: z.len, init: index}.filter(z[it] > 0.0)
	mut d := []f64{len: z.len}
	for i in present {
		d[i] = math.log(z[i]) + st.ln_phi[i]
	}
	return Feed{
		m:       m
		t:       t
		p:       p
		z:       z
		d:       d
		present: present
	}
}

fn check_objective(f Feed, a []f64) ! {
	k := a.len
	mut g := []f64{len: k}
	mut h := []f64{len: k * k}
	f.objective(a, mut g, mut h)!
	for j in 0 .. k {
		step := 1e-6 * a[j]
		mut ap := a.clone()
		mut am := a.clone()
		ap[j] += step
		am[j] -= step
		mut gp := []f64{len: k}
		mut gm := []f64{len: k}
		mut scratch := []f64{len: k * k}
		fp := f.objective(ap, mut gp, mut scratch)!
		fm := f.objective(am, mut gm, mut scratch)!
		fd := (fp - fm) / (2.0 * step)
		assert math.abs(fd - g[j]) <= 1e-6 * math.max(math.abs(g[j]), 1e-3), 'g[${j}] ${g[j]} vs ${fd}'
		for i in 0 .. k {
			fd2 := (gp[i] - gm[i]) / (2.0 * step)
			assert math.abs(fd2 - h[i * k + j]) <= 1e-6 * math.max(math.abs(h[i * k + j]), 1e-3), 'H[${i}][${j}] ${h[i * k + j]} vs ${fd2}'
		}
	}
	// H is symmetric.
	for i in 0 .. k {
		for j in 0 .. i {
			assert math.abs(h[i * k + j] - h[j * k + i]) <= 1e-12 * math.max(1.0, math.abs(h[i * k + j]))
		}
	}
}

fn test_the_newton_finish_has_the_right_derivatives() {
	methane := thermo.Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	butane := thermo.Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	co2 := thermo.Component{
		name:  'carbon dioxide'
		tc:    304.2
		pc:    73.83e5
		omega: 0.224
		mw:    0.04401
	}
	binary := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.02, 0.02, 0.0])!
	fb := feed_at(binary, 250.0, 30e5, [0.4, 0.6])!
	// W far from z and with ΣW ≠ 1, so that r, the 1/ΣW and the r/2 terms
	// all count.
	check_objective(fb, [2.0 * math.sqrt(0.9), 2.0 * math.sqrt(0.05)])!
	check_objective(fb, [2.0 * math.sqrt(0.3), 2.0 * math.sqrt(1.2)])!
	ternary := thermo.new_cubic(.pr, [methane, butane, co2],
		kij: [0.0, 0.02, 0.10, 0.02, 0.0, 0.13, 0.10, 0.13, 0.0]
	)!
	ft := feed_at(ternary, 300.0, 50e5, [0.5, 0.3, 0.2])!
	check_objective(ft, [2.0 * math.sqrt(0.7), 2.0 * math.sqrt(0.02), 2.0 * math.sqrt(0.4)])!
	// An absent component stays out: α covers the components present only.
	fa := feed_at(ternary, 300.0, 50e5, [0.5, 0.5, 0.0])!
	assert fa.present == [0, 1]
	check_objective(fa, [2.0 * math.sqrt(0.8), 2.0 * math.sqrt(0.3)])!
}
