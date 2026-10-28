module flash

import math
import thermo

fn methane_butane() !thermo.Cubic {
	return thermo.new_cubic(.pr, [
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
}

fn test_the_jacobian_matches_differences() {
	// eval's analytic Jacobian against central differences of F, for each
	// kind of vapour-fraction flash, away from the solution.
	m := methane_butane()!
	z := [0.4, 0.6]
	for given in [Given.t, Given.p] {
		for beta in [0.0, 0.3, 1.0] {
			value := if given == .t { 250.0 } else { 30e5 }
			s := Saturation{
				m:       m
				given:   given
				value:   value
				beta:    beta
				z:       z
				present: [0, 1]
				params:  VfParams{}
			}
			v := s.start()!
			n := v.len
			mut f := []f64{len: n}
			mut jac := []f64{len: n * n}
			s.eval(v, mut f, mut jac)!
			for col in 0 .. n {
				h := 1e-6
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
					assert err <= 1e-6 * math.max(1.0, math.abs(fdiff)), '${given} β ${beta}: J[${row}][${col}] ${jac[row * n + col]} vs ${fdiff}'
				}
			}
		}
	}
}

fn test_phases_the_wrong_way_round_are_refused() {
	// result is handed a real bubble point with its phases swapped, so that
	// the phase of fraction β is the denser: it must refuse, not relabel.
	m := methane_butane()!
	z := [0.4, 0.6]
	r := bubble_p(m, 250.0, z)!
	vap, liq := r.phases[0].state, r.phases[1].state
	s := Saturation{
		m:       m
		given:   .t
		value:   250.0
		beta:    0.0
		z:       z
		present: [0, 1]
		params:  VfParams{}
	}
	s.result(r.t, r.p, vap.n, vap, liq.n, liq, 0) or {
		assert (err as FlashError).kind == .phases_swapped, err.msg()
		return
	}
	assert false, 'swapped phases were accepted'
}

fn test_a_metastable_root_is_refused() {
	// Nearly pure n-butane at 300 K and 2 bar, below its vapour pressure,
	// has a liquid root of higher Gibbs energy than its vapour root. A
	// "boundary" with a phase there is between metastable phases.
	m := methane_butane()!
	w := [0.01, 0.99]
	liq := thermo.solve_tp(m, 300.0, 2e5, w, .liquid)!
	vap := thermo.solve_tp(m, 300.0, 2e5, w, .vapor)!
	stable := thermo.solve_tp(m, 300.0, 2e5, w, .stable)!
	assert liq.v < vap.v && stable.v == vap.v, 'the liquid root should be metastable here'
	s := Saturation{
		m:       m
		given:   .t
		value:   300.0
		beta:    0.0
		z:       w
		present: [0, 1]
		params:  VfParams{}
	}
	cur := SatPoint{
		t:     300.0
		p:     2e5
		d:     [1.0, 1.0]
		x:     w
		y:     w
		sum_x: 1.0
		sum_y: 1.0
		sx:    liq
		sy:    vap
	}
	s.on_stable_roots(cur) or {
		assert (err as FlashError).kind == .no_convergence && err.msg().contains('higher Gibbs'), err.msg()
		return
	}
	assert false, 'a metastable root was accepted'
}

fn test_a_collapse_onto_the_feed_is_refused() {
	// At 250 K and 200 bar the feed has one volume root, so K = 1 solves
	// the equations: Newton started there stops at once, and attempt must
	// say the iteration collapsed rather than call it a bubble point.
	m := methane_butane()!
	s := Saturation{
		m:       m
		given:   .t
		value:   250.0
		beta:    0.0
		z:       [0.4, 0.6]
		present: [0, 1]
		params:  VfParams{}
	}
	s.attempt([0.0, 0.0, math.log(200e5)]) or {
		assert (err as FlashError).kind == .no_convergence && err.msg().contains('collapsed'), err.msg()
		return
	}
	assert false, 'a collapse was accepted'
}
