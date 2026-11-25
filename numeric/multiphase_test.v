module numeric

import math

// A seeded generator of its own, so these tests don't move with another's.
struct Lcg {
mut:
	s u64
}

fn (mut r Lcg) f64() f64 {
	r.s = r.s * 6364136223846793005 + 1442695040888963407
	return f64(r.s >> 11) / f64(u64(1) << 53)
}

// With two phases, Q's minimum is Rachford–Rice's root where that lies in
// (0, 1), and a single phase at β = 0 or 1 beyond it: 200 seeded feeds and
// K-values, four to eight components, K from 1e-3 to 1e3.
fn test_two_phases_are_rachford_rices() {
	mut r := Lcg{7}
	mut inside, mut beyond := 0, 0
	mut worst := 0.0
	mut steps := 0
	for _ in 0 .. 200 {
		nc := 4 + int(r.f64() * 5.0)
		mut z := []f64{len: nc, init: 0.05 + r.f64()}
		total := sum(z)
		z = z.map(it / total)
		lnk := []f64{len: nc, init: (r.f64() * 2.0 - 1.0) * math.log(1e3)}
		rr := rachford_rice(z, lnk.map(math.exp(it))) or { continue }
		m := multiphase_rr(z, [[]f64{len: nc}, lnk], [0.5, 0.5], 1e-14)!
		steps += m.iterations
		if rr.beta > 0.0 && rr.beta < 1.0 {
			inside++
			worst = math.max(worst, math.abs(m.beta[1] - rr.beta))
			for i in 0 .. nc {
				worst = math.max(worst, math.abs(m.x[1][i] - rr.y[i]))
				worst = math.max(worst, math.abs(m.x[0][i] - rr.x[i]))
			}
		} else {
			beyond++
			// The feed itself is the one phase, and the other's Σx ≤ 1.
			one := if rr.beta >= 1.0 { 1 } else { 0 }
			assert m.beta[one] == 1.0 && m.beta[1 - one] == 0.0, '${rr.beta}: ${m.beta}'
			assert m.sum_x[1 - one] <= 1.0 + 1e-12, '${m.sum_x}'
			for i in 0 .. nc {
				assert math.abs(m.x[one][i] - z[i]) <= 1e-14
			}
		}
	}
	println('two phases: ${inside} inside (0, 1), ${beyond} beyond, worst ${worst:.1e}, ${steps} steps')
	// Measured: a few ulps; the gradient's tolerance, 1e-14, is in Σx.
	// The seeded draw: 171 inside (0, 1), 18 beyond.
	assert inside == 171 && beyond == 18, '${inside}, ${beyond}'
	assert worst <= 1e-12, '${worst}'
	// Newton's steps, recorded, which can't grow: a line search that took a
	// rise in Q would take more.
	assert steps <= 1103, '${steps} steps'
}

// Three phases from a known answer: with Kᵢⱼ = xᵢⱼ/cᵢ (so that each
// phase's fugacity, xᵢⱼ/Kᵢⱼ, is the same), Q's minimum is the β the feed
// was made from. And a fourth phase whose K-values give Σᵢxᵢ = 0.6 there
// stays absent, with that sum reported.
fn test_three_phases_from_a_known_answer() {
	mut r := Lcg{11}
	mut worst := 0.0
	for _ in 0 .. 100 {
		nc := 3 + int(r.f64() * 6.0)
		mut x := [][]f64{len: 3, init: []f64{len: nc}}
		for j in 0 .. 3 {
			for i in 0 .. nc {
				x[j][i] = math.exp((r.f64() * 2.0 - 1.0) * 6.0)
			}
			s := sum(x[j])
			x[j] = x[j].map(it / s)
		}
		mut beta := [0.1 + r.f64(), 0.1 + r.f64(), 0.1 + r.f64()]
		bs := sum(beta)
		beta = beta.map(it / bs)
		z := []f64{len: nc, init: beta[0] * x[0][index] + beta[1] * x[1][index] + beta[2] * x[2][index]}
		c := []f64{len: nc, init: math.exp(r.f64() * 4.0 - 2.0)}
		mut lnk := [][]f64{len: 4, init: []f64{len: nc}}
		for j in 0 .. 3 {
			for i in 0 .. nc {
				lnk[j][i] = math.log(x[j][i] / c[i])
			}
		}
		// The fourth: K proportional to a composition w, scaled so that
		// Σᵢ zᵢKᵢ/Eᵢ = 0.6 at the answer, where Eᵢ = zᵢcᵢ.
		mut w := []f64{len: nc, init: r.f64() + 0.01}
		ws := sum(w)
		w = w.map(it / ws)
		for i in 0 .. nc {
			lnk[3][i] = math.log(0.6 * w[i] / c[i])
		}
		m := multiphase_rr(z, lnk, [0.25, 0.25, 0.25, 0.25], 1e-14)!
		for j in 0 .. 3 {
			worst = math.max(worst, math.abs(m.beta[j] - beta[j]))
			for i in 0 .. nc {
				worst = math.max(worst, math.abs(m.x[j][i] - x[j][i]) / x[j][i])
			}
		}
		assert m.beta[3] == 0.0, '${m.beta}'
		assert math.abs(m.sum_x[3] - 0.6) <= 1e-12, '${m.sum_x}'
	}
	println('three phases: worst ${worst:.1e}')
	// Measured: the gradient's tolerance over the Hessian's conditioning.
	assert worst <= 1e-10, '${worst}'
}

// A phase that leaves reaches 0 exactly. Two components in three phases,
// water and 1-butanol by NRTL a step from their heteroazeotrope's T (the
// published data's test, after 0.6.0's R moved it): the Hessian is
// singular, Newton's steps are 1e5 long, and the middle phase, cut to its
// bound each step, kept a residue of β·ε, until at 1.6e-319 its step −β/d
// underflowed to 0 and nothing moved. The inputs are the doubles that
// failed. Its answer is checked from scratch: the phases present have
// Σx = 1, the one held has Σx ≤ 1, and the feed is their sum.
fn test_a_leaving_phase_reaches_zero() {
	z := [0.8, 0.2]
	lnk := [[0.25265486030507067, -2.126001781856114], [-0.0, -0.0],
		[-0.4854475926527495, 0.778552092635867]]
	m := multiphase_rr(z, lnk, [0.21045616834899977, 0.7795438316510003, 0.01], 1e-14)!
	assert m.beta.filter(it == 0.0).len == 1, '${m.beta}'
	for j in 0 .. 3 {
		if m.beta[j] > 0.0 {
			assert math.abs(m.sum_x[j] - 1.0) <= 1e-14, '${m.sum_x}'
		} else {
			assert m.sum_x[j] <= 1.0, '${m.sum_x}'
		}
	}
	for i in 0 .. 2 {
		mut back := 0.0
		for j in 0 .. 3 {
			back += m.beta[j] * m.x[j][i]
		}
		assert math.abs(back - z[i]) <= 1e-14, 'component ${i}: ${back}'
	}
	println('a leaving phase: β = ${m.beta} in ${m.iterations} steps')
	// Measured: 4 steps.
	assert m.iterations <= 10
}

// K-values from 1e-30 to 1e25 are scaled, never overflowed: a heavy
// component's mole fraction 1e-30 in the light phase keeps its digits.
fn test_extreme_k_values() {
	z := [0.3, 0.4, 0.3]
	lnk := [[0.0, 0.0, 0.0], [math.log(1e25), 0.0, math.log(1e-30)]]
	rr := rachford_rice(z, [1e25, 1.0, 1e-30])!
	m := multiphase_rr(z, lnk, [0.5, 0.5], 1e-14)!
	assert math.abs(m.beta[1] / rr.beta - 1.0) <= 1e-12, '${m.beta} ${rr.beta}'
	assert math.abs(m.x[1][2] / rr.y[2] - 1.0) <= 1e-12, '${m.x[1][2]} ${rr.y[2]}'
}

fn test_bad_input_is_refused() {
	for args in [
		[[0.5, 0.5], [0.5]],
		[[0.5, -0.5], [0.5, 0.5]],
		[[0.5, math.nan()], [0.5, 0.5]],
	] {
		if _ := multiphase_rr(args[0], [[0.0, 0.0], [1.0, -1.0]], args[1], 1e-12) {
			assert false, '${args}'
		}
	}
	if _ := multiphase_rr([0.5, 0.5], [[0.0, 0.0], [1.0, -1.0]], [0.0, 0.0], 1e-12) {
		assert false, 'a start with no phase'
	}
	if _ := multiphase_rr([0.5, 0.5], [[0.0, 0.0], [math.inf(1), -1.0]], [0.5, 0.5], 1e-12) {
		assert false, 'an infinite K'
	}
}
