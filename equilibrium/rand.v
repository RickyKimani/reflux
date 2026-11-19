module equilibrium

import math
import numeric

// The RAND method (White, Johnson and Dantzig, 1958; D9) in the form NASA's
// CEA uses (Gordon and McBride, NASA RP-1311, 1994, §2.5 and §3.3): Newton
// on the constraints' potentials π and on ln n, the total amount as a
// variable of its own, with every species' amount in log space, so that
// one that vanishes is a large negative ln nⱼ and never a negative nⱼ.
//
// It minimises
//
//     G/RT = Σⱼ nⱼ μⱼ/RT,   μⱼ/RT = gⱼ + ln(P/p°) + ln(nⱼ/n) + ln φⱼ,
//
// subject to C·n = b, at fixed φ: the caller updates φ between calls for a
// non-ideal gas (D9). C's rows are the elements, or for an equilibrium over
// chosen reactions, what those reactions conserve. Each step solves, for π
// and Δln n,
//
//     Σᵢ (Σⱼ cₖⱼcᵢⱼnⱼ) πᵢ + (Σⱼ cₖⱼnⱼ) Δln n = bₖ − Σⱼ cₖⱼnⱼ + Σⱼ cₖⱼnⱼ μⱼ/RT
//     Σᵢ (Σⱼ cᵢⱼnⱼ) πᵢ + (Σⱼ nⱼ − n) Δln n = n − Σⱼ nⱼ + Σⱼ nⱼ μⱼ/RT
//
// then Δln nⱼ = −μⱼ/RT + Σᵢ cᵢⱼπᵢ + Δln n, damped by CEA's control factor.

// max_steps bounds the RAND steps of one minimisation.
const max_steps = 500

// rand_tol is how far a full step may move every species' ln n for the
// minimisation to have converged: Newton's error after it is of the order of
// its square, rounding's. The answer is the species' amounts, so ln n, the
// total as a variable, needs no test of its own.
const rand_tol = 1e-10

// Problem is a minimisation: the constraints' rows over the species, their
// right-hand sides, and each species' standard potential μ°/RT, or the
// stand-in an equilibrium over chosen reactions uses.
struct Problem {
	c    [][]f64
	b    []f64
	g    []f64
	ln_p f64 // ln(P/p°)
}

// Point is where a minimisation is: each species' ln n, ln n itself, the
// constraints' potentials, and how many steps it has taken.
struct Point {
mut:
	ln_n     []f64
	ln_total f64
	pi       []f64
	steps    int
}

// minimise takes RAND steps from pt at fixed ln φ until a full step moves
// no species' ln n by more than rand_tol, or refuses after max_steps. A step
// that small is never damped, the control factor being 1 for it.
fn (pr Problem) minimise(ln_phi []f64, mut pt Point) ! {
	s := pr.g.len
	m := pr.c.len
	w := m + 1
	mut n := []f64{len: s}
	mut mu := []f64{len: s}
	mut a := []f64{len: w * w}
	mut rhs := []f64{len: w}
	mut dln := []f64{len: s}
	for pt.steps < max_steps {
		pt.steps++
		total := math.exp(pt.ln_total)
		mut sum := 0.0
		for j in 0 .. s {
			n[j] = math.exp(pt.ln_n[j])
			mu[j] = pr.g[j] + pr.ln_p + pt.ln_n[j] - pt.ln_total + ln_phi[j]
			sum += n[j]
		}
		mut nmu := 0.0
		for j in 0 .. s {
			nmu += n[j] * mu[j]
		}
		for k in 0 .. m {
			mut bk := 0.0
			mut cmu := 0.0
			for j in 0 .. s {
				bk += pr.c[k][j] * n[j]
				cmu += pr.c[k][j] * n[j] * mu[j]
			}
			for i in 0 .. m {
				mut x := 0.0
				for j in 0 .. s {
					x += pr.c[k][j] * pr.c[i][j] * n[j]
				}
				a[k * w + i] = x
			}
			a[k * w + m] = bk
			a[m * w + k] = bk
			rhs[k] = pr.b[k] - bk + cmu
		}
		a[m * w + m] = sum - total
		rhs[m] = total - sum + nmu
		lu := numeric.lu_factor(a, w) or {
			return fail(.not_converged, 'the constraints are singular where the species are: ${err.msg()}')
		}
		x := lu.solve(rhs) or { return fail(.not_converged, err.msg()) }
		dln_total := x[m]
		for j in 0 .. s {
			mut cpi := 0.0
			for i in 0 .. m {
				cpi += pr.c[i][j] * x[i]
			}
			dln[j] = -mu[j] + cpi + dln_total
		}
		pt.pi = x[..m].clone()
		// CEA's control factor: no major species' ln n moves by more than 2,
		// nor ln n by more than 0.4, and no trace species (below 1e-8 of the
		// total) rises past 1e-4 of it in one step.
		mut big := 5 * math.abs(dln_total)
		mut lambda2 := 1.0
		for j in 0 .. s {
			share := pt.ln_n[j] - pt.ln_total
			if share > trace_share {
				big = math.max(big, math.abs(dln[j]))
			} else if dln[j] >= 0 && dln[j] - dln_total != 0 {
				lambda2 = math.min(lambda2, math.abs((-share - 9.2103404) / (dln[j] - dln_total)))
			}
		}
		lambda := math.min(1.0, math.min(if big > 0 { 2 / big } else { 1.0 }, lambda2))
		for j in 0 .. s {
			pt.ln_n[j] += lambda * dln[j]
		}
		pt.ln_total += lambda * dln_total
		if dln.all(math.abs(it) <= rand_tol) {
			return
		}
	}
	return fail(.not_converged, 'the minimisation took ${max_steps} steps without converging')
}

// trace_share is ln 1e-8: a species below that share of the total is a
// trace one, for the control factor.
const trace_share = -18.420681

// numeric_lu factors a small matrix, or the identity where it's singular,
// for a stand-in that a later check would then refuse.
fn numeric_lu(a []f64, n int) numeric.LU {
	return numeric.lu_factor(a, n) or {
		mut id := []f64{len: n * n}
		for i in 0 .. n {
			id[i * n + i] = 1
		}
		numeric.lu_factor(id, n) or { panic(err) }
	}
}
