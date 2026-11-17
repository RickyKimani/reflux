module numeric

import math

// Qp is a strictly convex quadratic programme, the subproblem each step of
// sqp solves (M8's D12):
//
//     minimise ½·xᵀ·H·x + gᵀ·x
//     subject to rows[i]·x = b[i] for i < eq, and rows[i]·x ≥ b[i] for the rest,
//
// with H n×n, row-major and symmetric positive definite.
pub struct Qp {
pub:
	h    []f64
	g    []f64
	rows [][]f64
	b    []f64
	eq   int // how many of rows, from the first, are equalities
}

// QpSolution is a QP's minimum: x, each row's multiplier u (zero where the
// row isn't active, of either sign for an equality and at least zero for an
// inequality), so that H·x + g = Σ uᵢ·rowsᵢ, and which rows are active.
pub struct QpSolution {
pub:
	x          []f64
	u          []f64
	active     []bool
	iterations int // steps of the dual method
}

// QpFailure is what kind of failure a QpError is.
pub enum QpFailure {
	bad_input      // sizes that don't agree, or a value that isn't finite
	not_convex     // H isn't positive definite
	infeasible     // no x meets every row
	too_many_steps // the dual method's budget spent, which rounding can cause
}

// QpError is qp's failure.
pub struct QpError {
	Error
pub:
	kind       QpFailure
	reason     string
	iterations int
}

// msg is the reason.
pub fn (e QpError) msg() string {
	return 'qp: ${e.reason}'
}

fn qp_fail(kind QpFailure, reason string, iterations int) IError {
	return QpError{
		kind:       kind
		reason:     reason
		iterations: iterations
	}
}

// qp solves a strictly convex QP by Goldfarb and Idnani's dual active-set
// method (Math. Programming 27, 1983). It starts from the unconstrained
// minimum, −H⁻¹·g, and adds one row at a time: each equality not met, then
// the inequality most violated, by its distance s/‖row‖. Adding a row p moves x
// along z and the active rows' multipliers along −r, where
//
//     H·z + N·r = nₚ,   Nᵀ·z = 0,
//
// N the active rows as columns, so that the active rows stay met while p's
// multiplier grows by t. The step t is the smaller of
// - t₂ = −sₚ/(zᵀ·nₚ), which meets p, and
// - t₁, the least uⱼ/rⱼ over active inequalities with rⱼ > 0, which takes one
//   of them to zero: that row is dropped, and p is tried again.
//
// z = 0 means nₚ is a combination of the active rows; with no t₁ either,
// nothing can meet p, and the QP is infeasible. Every step keeps the
// multipliers of the inequalities at least zero and lowers no constraint
// already met, so the method ends, in exact arithmetic, at the minimum or at
// a proof that there's none.
//
// The problems here are small, a handful of variables and rows, so z and r
// come from the whole system, factorised afresh at each step, rather than
// the updated factorisations a large QP needs.
pub fn qp(p Qp) !QpSolution {
	n := p.g.len
	m := p.rows.len
	if n == 0 || p.h.len != n * n || p.b.len != m || p.eq < 0 || p.eq > m {
		return qp_fail(.bad_input, 'expected H n×n, g of n, b one a row, and eq from 0 to the rows',
			0)
	}
	for i, row in p.rows {
		if row.len != n {
			return qp_fail(.bad_input, 'row ${i} has ${row.len} entries, not ${n}', 0)
		}
		if !math.is_finite(p.b[i]) || row.any(!math.is_finite(it)) {
			return qp_fail(.bad_input, 'row ${i} is not finite', 0)
		}
	}
	if p.h.any(!math.is_finite(it)) || p.g.any(!math.is_finite(it)) {
		return qp_fail(.bad_input, 'H or g is not finite', 0)
	}
	ch := cholesky(p.h, n) or {
		return qp_fail(.not_convex, "H isn't positive definite: ${err.msg()}", 0)
	}
	mut x := ch.solve(p.g.map(-it))!
	// The active rows in the order added, each with its sign (an equality is
	// added the way it's violated, as an inequality, and never dropped) and
	// its multiplier.
	mut active := []int{}
	mut sign := []f64{len: m, init: 1.0}
	mut u := []f64{}
	mut added := []bool{len: m}
	mut norms := []f64{len: m}
	for i, row in p.rows {
		norms[i] = math.sqrt(dot(row, row))
	}
	mut steps := 0
	most := 50 + 10 * (n + m)
	for {
		// The next row: an equality not yet added and not met, else the most
		// violated inequality. An equality met though not active, as one
		// that repeats another, stays out unless a later step moves x off
		// it.
		mut next := -1
		for i in 0 .. p.eq {
			ax := dot(p.rows[i], x)
			if !added[i] && math.abs(ax - p.b[i]) > violation_tol(p.b[i], ax) {
				next = i
				break
			}
		}
		if next < 0 {
			mut worst := 0.0
			for i in p.eq .. m {
				if added[i] {
					continue
				}
				ax := dot(p.rows[i], x)
				s := ax - p.b[i]
				if s < -violation_tol(p.b[i], ax) && norms[i] > 0.0 && -s / norms[i] > worst {
					worst = -s / norms[i]
					next = i
				} else if s < -violation_tol(p.b[i], ax) && norms[i] == 0.0 {
					return qp_fail(.infeasible, 'row ${i} is 0 ≥ ${p.b[i]}', steps)
				}
			}
		}
		if next < 0 {
			break
		}
		if next < p.eq {
			sign[next] = if dot(p.rows[next], x) - p.b[next] > 0.0 { -1.0 } else { 1.0 }
		}
		np := p.rows[next].map(sign[next] * it)
		bp := sign[next] * p.b[next]
		// nₚᵀ·H⁻¹·nₚ, z's size were nothing active: the scale for z = 0.
		w := ch.solve(np)!
		reference := dot(np, w)
		mut up := 0.0
		for {
			steps++
			if steps > most {
				return qp_fail(.too_many_steps, 'no minimum in ${most} steps', steps)
			}
			z, r := qp_direction(p, active, sign, np) or {
				return qp_fail(.too_many_steps, 'the active rows became dependent: ${err.msg()}',
					steps)
			}
			mut t1 := math.inf(1)
			mut k := -1
			for j, row in active {
				if row >= p.eq && r[j] > 0.0 {
					t := u[j] / r[j]
					if t < t1 {
						t1, k = t, j
					}
				}
			}
			ztn := dot(z, np)
			sp := dot(np, x) - bp
			t2 := if ztn > 1e-12 * reference { math.max(-sp / ztn, 0.0) } else { math.inf(1) }
			t := math.min(t1, t2)
			if math.is_inf(t, 1) {
				return qp_fail(.infeasible, 'no point meets row ${next} with the rows active', steps)
			}
			for j in 0 .. active.len {
				u[j] -= t * r[j]
			}
			up += t
			if !math.is_inf(t2, 1) {
				for i in 0 .. n {
					x[i] += t * z[i]
				}
			}
			if !math.is_inf(t2, 1) && t2 <= t1 {
				active << next
				u << up
				added[next] = true
				break
			}
			added[active[k]] = false
			active.delete(k)
			u.delete(k)
		}
	}
	mut out_u := []f64{len: m}
	mut on := []bool{len: m}
	for j, row in active {
		out_u[row] = sign[row] * u[j]
		on[row] = true
	}
	return QpSolution{
		x:          x
		u:          out_u
		active:     on
		iterations: steps
	}
}

// violation_tol is how far below its bound a row may be and count as met:
// rounding in rows·x and b.
fn violation_tol(b f64, ax f64) f64 {
	return 1e-12 * math.max(1.0, math.max(math.abs(b), math.abs(ax)))
}

// qp_direction solves H·z + N·r = nₚ, Nᵀ·z = 0 for z and r, N the active
// rows, signed, as columns.
fn qp_direction(p Qp, active []int, sign []f64, np []f64) !([]f64, []f64) {
	n := p.g.len
	q := active.len
	size := n + q
	mut k := []f64{len: size * size}
	for i in 0 .. n {
		for j in 0 .. n {
			k[i * size + j] = p.h[i * n + j]
		}
	}
	for j, row in active {
		for i in 0 .. n {
			v := sign[row] * p.rows[row][i]
			k[i * size + n + j] = v
			k[(n + j) * size + i] = v
		}
	}
	mut rhs := []f64{len: size}
	for i in 0 .. n {
		rhs[i] = np[i]
	}
	sol := lu_factor(k, size)!.solve(rhs)!
	return sol[..n], sol[n..]
}
