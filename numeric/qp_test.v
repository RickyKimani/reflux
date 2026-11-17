module numeric

import math
import rand

// kkt_holds checks a QP's answer from scratch, not from the method's
// bookkeeping: every row met, H·x + g = Σ uᵢ·rowsᵢ, each inequality's
// multiplier at least zero and zero unless its row is met with equality.
// The tolerances are rounding in a few dozen products of numbers of order
// one.
fn kkt_holds(p Qp, s QpSolution) {
	n := p.g.len
	for i, row in p.rows {
		r := dot(row, s.x) - p.b[i]
		if i < p.eq {
			assert math.abs(r) <= 1e-10, 'row ${i}: ${r}'
		} else {
			assert r >= -1e-10, 'row ${i}: ${r}'
			assert s.u[i] >= 0.0, 'row ${i}: u ${s.u[i]}'
			assert math.abs(s.u[i] * r) <= 1e-10, 'row ${i}: u·r ${s.u[i] * r}'
		}
		assert s.active[i] == (s.u[i] != 0.0) || s.active[i]
	}
	for j in 0 .. n {
		mut r := p.g[j]
		for k in 0 .. n {
			r += p.h[j * n + k] * s.x[k]
		}
		for i, row in p.rows {
			r -= s.u[i] * row[j]
		}
		assert math.abs(r) <= 1e-10, 'stationarity ${j}: ${r}'
	}
}

// Nocedal and Wright's example 16.4: (x₁ − 1)² + (x₂ − 2.5)² over five
// inequalities, whose minimum (1.4, 1.7) has the first alone active, with
// multiplier 0.8.
fn test_nocedal_and_wrights_example() {
	p := Qp{
		h:    [2.0, 0, 0, 2]
		g:    [-2.0, -5]
		rows: [[1.0, -2], [-1.0, -2], [-1.0, 2], [1.0, 0], [0.0, 1]]
		b:    [-2.0, -6, -2, 0, 0]
	}
	s := qp(p)!
	kkt_holds(p, s)
	assert math.abs(s.x[0] - 1.4) <= 1e-15 && math.abs(s.x[1] - 1.7) <= 1e-15, '${s.x}'
	assert math.abs(s.u[0] - 0.8) <= 1e-15 && s.active == [true, false, false, false, false]
}

// Projection onto a box: x is p clamped, and each bound's multiplier how
// far p is beyond it.
fn test_a_box() {
	target := [3.0, -0.5, 0.25, -7.0]
	mut rows := [][]f64{}
	mut b := []f64{}
	for j in 0 .. 4 {
		mut lo := []f64{len: 4}
		lo[j] = 1.0
		rows << lo
		b << -1.0
		mut hi := []f64{len: 4}
		hi[j] = -1.0
		rows << hi
		b << -1.0
	}
	mut h := []f64{len: 16}
	for j in 0 .. 4 {
		h[j * 5] = 1.0
	}
	p := Qp{
		h:    h
		g:    target.map(-it)
		rows: rows
		b:    b
	}
	s := qp(p)!
	kkt_holds(p, s)
	assert s.x == [1.0, -0.5, 0.25, -1.0]
	assert s.u == [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0]
}

// Equalities: ½‖x‖² with Σ x = 1 is 1/n each; the same row twice is met
// once and the repeat left out; a third equality inconsistent with the
// first is infeasible.
fn test_equalities() {
	mut h := []f64{len: 9}
	for j in 0 .. 3 {
		h[j * 4] = 1.0
	}
	p := Qp{
		h:    h
		g:    [0.0, 0, 0]
		rows: [[1.0, 1, 1], [2.0, 2, 2]]
		b:    [1.0, 2]
		eq:   2
	}
	s := qp(p)!
	kkt_holds(p, s)
	for v in s.x {
		assert math.abs(v - 1.0 / 3.0) <= 1e-15
	}
	if _ := qp(Qp{
		...p
		b: [1.0, 3]
	})
	{
		assert false, 'an inconsistent equality met'
	} else {
		assert err is QpError && (err as QpError).kind == .infeasible, err.msg()
	}
}

// No point at all: x₁ ≥ 1 and −x₁ ≥ 0; a row of zeros above zero; an H
// that isn't positive definite; sizes that disagree.
fn test_refusals() {
	h := [1.0, 0, 0, 1]
	cases := [
		Qp{
			h:    h
			g:    [0.0, 0]
			rows: [[1.0, 0], [-1.0, 0]]
			b:    [1.0, 0]
		},
		Qp{
			h:    h
			g:    [0.0, 0]
			rows: [[0.0, 0]]
			b:    [1.0]
		},
		Qp{
			h: [1.0, 2, 2, 1]
			g: [0.0, 0]
		},
		Qp{
			h:    h
			g:    [0.0, 0]
			rows: [[1.0]]
			b:    [0.0]
		},
	]
	kinds := [QpFailure.infeasible, .infeasible, .not_convex, .bad_input]
	for k, p in cases {
		if _ := qp(p) {
			assert false, '${k}: solved'
		} else {
			assert err is QpError && (err as QpError).kind == kinds[k], '${k}: ${err.msg()}'
		}
	}
}

// Random QPs with a point known to meet every row, each answer checked
// from scratch: 2 to 6 variables, up to n − 1 equalities and 2n
// inequalities, many of them active.
fn test_random_qps() {
	rand.seed([u32(2026), 1117])
	mut most := 0
	for _ in 0 .. 500 {
		n := rand.intn(5) or { 0 } + 2
		eq := rand.intn(n) or { 0 }
		m := eq + rand.intn(2 * n + 1) or { 0 }
		mut l := []f64{len: n * n}
		for i in 0 .. n {
			for j in 0 .. i + 1 {
				l[i * n + j] = rand.f64() * 2.0 - 1.0
			}
			l[i * n + i] = 0.5 + rand.f64()
		}
		mut h := []f64{len: n * n}
		for i in 0 .. n {
			for j in 0 .. n {
				for k in 0 .. n {
					h[i * n + j] += l[i * n + k] * l[j * n + k]
				}
			}
		}
		feasible := []f64{len: n, init: rand.f64() * 4.0 - 2.0}
		mut rows := [][]f64{}
		mut b := []f64{}
		for i in 0 .. m {
			row := []f64{len: n, init: rand.f64() * 2.0 - 1.0}
			slack := if i < eq { 0.0 } else { rand.f64() * 0.5 }
			b << dot(row, feasible) - slack
			rows << row
		}
		p := Qp{
			h:    h
			g:    []f64{len: n, init: rand.f64() * 10.0 - 5.0}
			rows: rows
			b:    b
			eq:   eq
		}
		s := qp(p)!
		kkt_holds(p, s)
		most = math.max(most, s.iterations)
	}
	println('random QPs: ${most} steps at most')
	// Measured: 10 steps at most, which can't grow.
	assert most <= 10, '${most}'
}
