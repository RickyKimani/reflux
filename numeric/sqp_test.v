module numeric

import math
import os
import x.json2

// Hs is one of Hock and Schittkowski's problems (criterion 7), written again
// here from the book, as tools/gen_m8_hock_schittkowski.py writes it: f, its
// equalities then its inequalities c ≥ 0, and their gradients by hand,
// checked against differences in test_the_gradients_by_hand.
struct Hs {
	number int
mut:
	calls int
}

fn (h Hs) f(x []f64) f64 {
	return match h.number {
		1, 2, 15 {
			100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1.0 - x[0]) * (1.0 - x[0])
		}
		3 {
			x[1] + 1e-5 * (x[1] - x[0]) * (x[1] - x[0])
		}
		4 {
			math.pow(x[0] + 1.0, 3) / 3.0 + x[1]
		}
		5 {
			math.sin(x[0] + x[1]) + (x[0] - x[1]) * (x[0] - x[1]) - 1.5 * x[0] + 2.5 * x[1] + 1.0
		}
		6 {
			(1.0 - x[0]) * (1.0 - x[0])
		}
		7 {
			math.log(1.0 + x[0] * x[0]) - x[1]
		}
		9 {
			math.sin(math.pi * x[0] / 12.0) * math.cos(math.pi * x[1] / 16.0)
		}
		10 {
			x[0] - x[1]
		}
		11 {
			(x[0] - 5.0) * (x[0] - 5.0) + x[1] * x[1] - 25.0
		}
		12 {
			0.5 * x[0] * x[0] + x[1] * x[1] - x[0] * x[1] - 7.0 * x[0] - 7.0 * x[1]
		}
		14 {
			(x[0] - 2.0) * (x[0] - 2.0) + (x[1] - 1.0) * (x[1] - 1.0)
		}
		21 {
			x[0] * x[0] / 100.0 + x[1] * x[1] - 100.0
		}
		28 {
			(x[0] + x[1]) * (x[0] + x[1]) + (x[1] + x[2]) * (x[1] + x[2])
		}
		35 {
			9.0 - 8.0 * x[0] - 6.0 * x[1] - 4.0 * x[2] + 2.0 * x[0] * x[0] + 2.0 * x[1] * x[1] +
				x[2] * x[2] + 2.0 * x[0] * x[1] + 2.0 * x[0] * x[2]
		}
		39 {
			-x[0]
		}
		43 {
			x[0] * x[0] + x[1] * x[1] + 2.0 * x[2] * x[2] + x[3] * x[3] - 5.0 * x[0] - 5.0 * x[1] - 21.0 * x[2] +
				7.0 * x[3]
		}
		48 {
			(x[0] - 1.0) * (x[0] - 1.0) + (x[1] - x[2]) * (x[1] - x[2]) + (x[3] - x[4]) * (x[3] - x[4])
		}
		71 {
			x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]
		}
		100 {
			(x[0] - 10.0) * (x[0] - 10.0) + 5.0 * (x[1] - 12.0) * (x[1] - 12.0) + math.pow(x[2], 4) +
				3.0 * (x[3] - 11.0) * (x[3] - 11.0) + 10.0 * math.pow(x[4], 6) + 7.0 * x[5] * x[5] +
				math.pow(x[6], 4) - 4.0 * x[5] * x[6] - 10.0 * x[5] - 8.0 * x[6]
		}
		else {
			panic('no problem ${h.number}')
		}
	}
}

fn (h Hs) c(x []f64) []f64 {
	return match h.number {
		6 {
			[10.0 * (x[1] - x[0] * x[0])]
		}
		7 {
			[(1.0 + x[0] * x[0]) * (1.0 + x[0] * x[0]) + x[1] * x[1] - 4.0]
		}
		9 {
			[4.0 * x[0] - 3.0 * x[1]]
		}
		10 {
			[-3.0 * x[0] * x[0] + 2.0 * x[0] * x[1] - x[1] * x[1] + 1.0]
		}
		11 {
			[-x[0] * x[0] + x[1]]
		}
		12 {
			[25.0 - 4.0 * x[0] * x[0] - x[1] * x[1]]
		}
		14 {
			[x[0] - 2.0 * x[1] + 1.0, -x[0] * x[0] / 4.0 - x[1] * x[1] + 1.0]
		}
		15 {
			[x[0] * x[1] - 1.0, x[0] + x[1] * x[1]]
		}
		21 {
			[10.0 * x[0] - x[1] - 10.0]
		}
		28 {
			[x[0] + 2.0 * x[1] + 3.0 * x[2] - 1.0]
		}
		35 {
			[3.0 - x[0] - x[1] - 2.0 * x[2]]
		}
		39 {
			[x[1] - x[0] * x[0] * x[0] - x[2] * x[2], x[0] * x[0] - x[1] - x[3] * x[3]]
		}
		43 {
			[
				8.0 - x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3] - x[0] + x[1] - x[2] + x[3],
				10.0 - x[0] * x[0] - 2.0 * x[1] * x[1] - x[2] * x[2] - 2.0 * x[3] * x[3] + x[0] + x[3],
				5.0 - 2.0 * x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - 2.0 * x[0] + x[1] + x[3],
			]
		}
		48 {
			[x[0] + x[1] + x[2] + x[3] + x[4] - 5.0, x[2] - 2.0 * (x[3] + x[4]) + 3.0]
		}
		71 {
			[x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 40.0,
				x[0] * x[1] * x[2] * x[3] - 25.0]
		}
		100 {
			[
				127.0 - 2.0 * x[0] * x[0] - 3.0 * math.pow(x[1], 4) - x[2] - 4.0 * x[3] * x[3] - 5.0 * x[4],
				282.0 - 7.0 * x[0] - 3.0 * x[1] - 10.0 * x[2] * x[2] - x[3] + x[4],
				196.0 - 23.0 * x[0] - x[1] * x[1] - 6.0 * x[5] * x[5] + 8.0 * x[6],
				-4.0 * x[0] * x[0] - x[1] * x[1] + 3.0 * x[0] * x[1] - 2.0 * x[2] * x[2] - 5.0 * x[5] +
					11.0 * x[6],
			]
		}
		else {
			[]f64{}
		}
	}
}

// grad is ∇f and each ∇cᵢ.
fn (h Hs) grad(x []f64) ([]f64, [][]f64) {
	match h.number {
		1, 2, 15 {
			g := [-400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]),
				200.0 * (x[1] - x[0] * x[0])]
			if h.number == 15 {
				return g, [[x[1], x[0]], [1.0, 2.0 * x[1]]]
			}
			return g, [][]f64{}
		}
		3 {
			return [-2e-5 * (x[1] - x[0]), 1.0 + 2e-5 * (x[1] - x[0])], [][]f64{}
		}
		4 {
			return [(x[0] + 1.0) * (x[0] + 1.0), 1.0], [][]f64{}
		}
		5 {
			s := math.cos(x[0] + x[1])
			return [s + 2.0 * (x[0] - x[1]) - 1.5, s - 2.0 * (x[0] - x[1]) + 2.5], [][]f64{}
		}
		6 {
			return [-2.0 * (1.0 - x[0]), 0.0], [[-20.0 * x[0], 10.0]]
		}
		7 {
			q := 1.0 + x[0] * x[0]
			return [2.0 * x[0] / q, -1.0], [[4.0 * x[0] * q, 2.0 * x[1]]]
		}
		9 {
			a, b := math.pi * x[0] / 12.0, math.pi * x[1] / 16.0
			return [math.pi / 12.0 * math.cos(a) * math.cos(b),
				-math.pi / 16.0 * math.sin(a) * math.sin(b)], [
				[4.0, -3.0],
			]
		}
		10 {
			return [1.0, -1.0], [[-6.0 * x[0] + 2.0 * x[1], 2.0 * x[0] - 2.0 * x[1]]]
		}
		11 {
			return [2.0 * (x[0] - 5.0), 2.0 * x[1]], [[-2.0 * x[0], 1.0]]
		}
		12 {
			return [x[0] - x[1] - 7.0, 2.0 * x[1] - x[0] - 7.0], [
				[-8.0 * x[0], -2.0 * x[1]],
			]
		}
		14 {
			return [2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)], [[1.0, -2.0], [-x[0] / 2.0, -2.0 * x[1]]]
		}
		21 {
			return [x[0] / 50.0, 2.0 * x[1]], [[10.0, -1.0]]
		}
		28 {
			a, b := 2.0 * (x[0] + x[1]), 2.0 * (x[1] + x[2])
			return [a, a + b, b], [[1.0, 2.0, 3.0]]
		}
		35 {
			return [-8.0 + 4.0 * x[0] + 2.0 * x[1] + 2.0 * x[2], -6.0 + 4.0 * x[1] + 2.0 * x[0],
				-4.0 + 2.0 * x[2] + 2.0 * x[0]], [[-1.0, -1.0, -2.0]]
		}
		39 {
			return [-1.0, 0.0, 0.0, 0.0], [
				[-3.0 * x[0] * x[0], 1.0, -2.0 * x[2], 0.0],
				[2.0 * x[0], -1.0, 0.0, -2.0 * x[3]],
			]
		}
		43 {
			return [2.0 * x[0] - 5.0, 2.0 * x[1] - 5.0, 4.0 * x[2] - 21.0, 2.0 * x[3] + 7.0], [
				[-2.0 * x[0] - 1.0, -2.0 * x[1] + 1.0, -2.0 * x[2] - 1.0, -2.0 * x[3] + 1.0],
				[-2.0 * x[0] + 1.0, -4.0 * x[1], -2.0 * x[2], -4.0 * x[3] + 1.0],
				[-4.0 * x[0] - 2.0, -2.0 * x[1] + 1.0, -2.0 * x[2], 1.0],
			]
		}
		48 {
			a, b := 2.0 * (x[1] - x[2]), 2.0 * (x[3] - x[4])
			return [2.0 * (x[0] - 1.0), a, -a, b, -b], [[1.0, 1, 1, 1, 1], [0.0, 0, 1, -2, -2]]
		}
		71 {
			s := x[0] + x[1] + x[2]
			return [x[3] * (s + x[0]), x[0] * x[3], x[0] * x[3] + 1.0, x[0] * s], [
				[2.0 * x[0], 2.0 * x[1], 2.0 * x[2], 2.0 * x[3]],
				[x[1] * x[2] * x[3], x[0] * x[2] * x[3], x[0] * x[1] * x[3], x[0] * x[1] * x[2]],
			]
		}
		100 {
			return [2.0 * (x[0] - 10.0), 10.0 * (x[1] - 12.0), 4.0 * math.pow(x[2], 3),
				6.0 * (x[3] - 11.0), 60.0 * math.pow(x[4], 5), 14.0 * x[5] - 4.0 * x[6] - 10.0,
				4.0 * math.pow(x[6], 3) - 4.0 * x[5] - 8.0], [
				[-4.0 * x[0], -12.0 * math.pow(x[1], 3), -1.0, -8.0 * x[3], -5.0, 0.0, 0.0],
				[-7.0, -3.0, -20.0 * x[2], -1.0, 1.0, 0.0, 0.0],
				[-23.0, -2.0 * x[1], 0.0, 0.0, 0.0, -12.0 * x[5], 8.0],
				[-8.0 * x[0] + 3.0 * x[1], -2.0 * x[1] + 3.0 * x[0], -4.0 * x[2], 0.0, 0.0, -5.0,
					11.0],
			]
		}
		else {
			panic('no problem ${h.number}')
		}
	}
}

fn (mut h Hs) values(x []f64, mut c []f64) !f64 {
	h.calls++
	for i, v in h.c(x) {
		c[i] = v
	}
	return h.f(x)
}

fn (mut h Hs) gradients(x []f64, mut g []f64, mut a []f64) ! {
	gf, ga := h.grad(x)
	n := x.len
	for j in 0 .. n {
		g[j] = gf[j]
	}
	for i, row in ga {
		for j in 0 .. n {
			a[i * n + j] = row[j]
		}
	}
}

// HsCase is a problem from the fixture: its sizes, bounds, start, and the
// optimum mpmath found again at 40 digits, with its multipliers.
struct HsCase {
	number       int
	n            int
	equalities   int
	inequalities int
	lower        []f64
	upper        []f64
	x0           []f64
	x            []f64
	f            f64
	lam          []f64
	mu_lower     []f64
	mu_upper     []f64
}

fn floats(v json2.Any, none_is f64) []f64 {
	return v.as_array().map(if it is json2.Null { none_is } else { it.f64() })
}

fn hs_cases() ![]HsCase {
	fixture := json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata',
		'm8_hock_schittkowski.json'))!)!.as_map()
	mut out := []HsCase{}
	for p in fixture['problems']!.as_array() {
		m := p.as_map()
		out << HsCase{
			number:       m['number']!.int()
			n:            m['n']!.int()
			equalities:   m['equalities']!.int()
			inequalities: m['inequalities']!.int()
			lower:        floats(m['lower']!, math.inf(-1))
			upper:        floats(m['upper']!, math.inf(1))
			x0:           floats(m['x0']!, 0.0)
			x:            floats(m['x']!, 0.0)
			f:            m['f_value']!.f64()
			lam:          floats(m['multipliers']!, 0.0)
			mu_lower:     floats(m['lower_multipliers']!, 0.0)
			mu_upper:     floats(m['upper_multipliers']!, 0.0)
		}
	}
	return out
}

// The gradients written by hand against central differences, at each
// start and optimum: a slip in them would send the SQP wrong for the
// test's reason, not its own. Differences at a step of 1e-6 are good to
// about 1e-9 of the scale here.
fn test_the_gradients_by_hand() {
	for p in hs_cases()! {
		h := Hs{
			number: p.number
		}
		for x in [p.x0, p.x] {
			g, a := h.grad(x)
			assert a.len == p.equalities + p.inequalities
			for j in 0 .. p.n {
				step := 1e-6 * math.max(1.0, math.abs(x[j]))
				mut up := x.clone()
				mut down := x.clone()
				up[j] += step
				down[j] -= step
				df := (h.f(up) - h.f(down)) / (up[j] - down[j])
				assert math.abs(df - g[j]) <= 1e-6 * math.max(1.0, math.abs(g[j])), 'HS${p.number} ∂f/∂x${j + 1}: ${g[j]}, by differences ${df}'
				cu, cd := h.c(up), h.c(down)
				for i in 0 .. a.len {
					dc := (cu[i] - cd[i]) / (up[j] - down[j])
					assert math.abs(dc - a[i][j]) <= 1e-6 * math.max(1.0, math.abs(a[i][j])), 'HS${p.number} ∂c${i + 1}/∂x${j + 1}'
				}
			}
		}
	}
}

// The iterations each problem took, measured, which can't grow.
const hs_iterations = {
	1:   24
	2:   21
	3:   4
	4:   2
	5:   7
	6:   12
	7:   10
	9:   8
	10:  13
	11:  8
	12:  10
	14:  6
	15:  3
	21:  2
	28:  5
	35:  12
	39:  15
	43:  11
	48:  9
	71:  8
	100: 16
}

// Every problem from its published start reaches its published optimum
// (criterion 7), each answer's first-order conditions checked here from
// the problem itself, not from sqp's bookkeeping:
// - c met: equalities to 1e-8, inequalities to within −1e-8 (tol);
// - ∇f − Σ λᵢ·∇cᵢ − μ_lower + μ_upper zero to 1e-7 of max(1, ‖∇f‖∞), ten
//   times sqp's own tol for the rounding in a sum of a dozen terms;
// - each inequality's and bound's multiplier at least zero, and zero where
//   it isn't active.
// Then x, f and the multipliers against mpmath's.
fn test_hock_and_schittkowskis_problems() {
	mut total := 0
	for p in hs_cases()! {
		mut h := Hs{
			number: p.number
		}
		s := sqp(mut h, p.x0,
			lower:        p.lower
			upper:        p.upper
			equalities:   p.equalities
			inequalities: p.inequalities
		) or {
			assert false, 'HS${p.number}: ${err.msg()}'
			return
		}
		at := 'HS${p.number}'
		g, a := h.grad(s.x)
		c := h.c(s.x)
		gscale := math.max(1.0, max_abs(g))
		for i, ci in c {
			if i < p.equalities {
				assert math.abs(ci) <= 1e-8, '${at}: c${i + 1} = ${ci}'
			} else {
				assert ci >= -1e-8, '${at}: c${i + 1} = ${ci}'
				assert s.multipliers[i] >= 0.0, '${at}: λ${i + 1} = ${s.multipliers[i]}'
				assert math.abs(s.multipliers[i] * ci) <= 1e-8 * gscale, '${at}: λ${i + 1}·c${i + 1}'
			}
		}
		for j in 0 .. p.n {
			assert s.x[j] >= p.lower[j] && s.x[j] <= p.upper[j], '${at}: x${j + 1} out of bounds'
			assert s.lower_multipliers[j] >= 0.0 && s.upper_multipliers[j] >= 0.0
			if s.lower_multipliers[j] > 0.0 {
				assert s.x[j] - p.lower[j] <= 1e-8 * math.max(1.0, math.abs(s.x[j])), '${at}: x${j + 1}'
			}
			if s.upper_multipliers[j] > 0.0 {
				assert p.upper[j] - s.x[j] <= 1e-8 * math.max(1.0, math.abs(s.x[j])), '${at}: x${j + 1}'
			}
			mut r := g[j] - s.lower_multipliers[j] + s.upper_multipliers[j]
			for i in 0 .. a.len {
				r -= s.multipliers[i] * a[i][j]
			}
			assert math.abs(r) <= 1e-7 * gscale, '${at}: ∂L/∂x${j + 1} = ${r}'
		}
		// Against mpmath. f is quadratic near its minimum, so it's good to
		// about tol²·scale; x to tol over the curvature, slower where the
		// curvature is small, as HS3's x₁ (1e-5).
		assert math.abs(s.f - p.f) <= 1e-9 * math.max(1.0, math.abs(p.f)), '${at}: f ${s.f}, mpmath ${p.f}'
		// HS9's optima are a family, x* = (12k − 3, 16k − 4) for every whole
		// k, each with f* = −1/2 and the same multiplier; the fixture's is
		// k = 0, and the one reached is held to its own k. From (0, 0) it's
		// k = −3: after the first step, f's curvature along the constraint
		// is near zero, so the first update's scaling makes B small and the
		// next step long, into another valley.
		mut want := p.x.clone()
		if p.number == 9 {
			k := math.round((s.x[0] + 3.0) / 12.0)
			want = [12.0 * k - 3.0, 16.0 * k - 4.0]
		}
		for j in 0 .. p.n {
			slack := if p.number == 3 { 1e-3 } else { 1e-6 }
			assert math.abs(s.x[j] - want[j]) <= slack * math.max(1.0, math.abs(want[j])), '${at}: x${j + 1} ${s.x[j]}, mpmath ${want[j]}'
		}
		for i in 0 .. a.len {
			assert math.abs(s.multipliers[i] - p.lam[i]) <= 1e-5 * math.max(1.0, math.abs(p.lam[i])), '${at}: λ${i + 1} ${s.multipliers[i]}, mpmath ${p.lam[i]}'
		}
		for j in 0 .. p.n {
			assert math.abs(s.lower_multipliers[j] - p.mu_lower[j]) <= 1e-5 * math.max(1.0,
				p.mu_lower[j])
			assert math.abs(s.upper_multipliers[j] - p.mu_upper[j]) <= 1e-5 * math.max(1.0,
				p.mu_upper[j])
		}
		assert s.evaluations == h.calls
		assert s.history.len == s.iterations && s.history.last().step == 0.0
		most := hs_iterations[p.number] or { 0 }
		println('${at}: ${s.iterations} iterations, ${s.evaluations} evaluations')
		assert s.iterations <= most, '${at}: ${s.iterations} iterations, recorded ${most}'
		total += s.iterations
	}
	println('${total} iterations in all')
}

// Fns is a problem as a struct of function values, for the refusals.
type ValuesFn = fn (x []f64, mut c []f64) !f64

type GradientsFn = fn (x []f64, mut g []f64, mut a []f64) !

struct Fns {
	v ValuesFn    @[required]
	g GradientsFn @[required]
}

fn (mut p Fns) values(x []f64, mut c []f64) !f64 {
	return p.v(x, mut c)
}

fn (mut p Fns) gradients(x []f64, mut g []f64, mut a []f64) ! {
	p.g(x, mut g, mut a)!
}

fn refused(mut p Constrained, x0 []f64, kind SqpFailure, says string, params SqpParams) {
	if s := sqp(mut p, x0, params) {
		assert false, 'solved at ${s.x}'
	} else {
		assert err is OptimizeError, err.msg()
		e := err as OptimizeError
		assert e.kind == kind && e.reason.contains(says), e.msg()
	}
}

// The problems sqp refuses (criterion 7), each with its reason:
// - x₁ + x₂ ≥ 3 and x₁ + x₂ ≤ 1, which no point meets, and a disc and a
//   half-plane that don't meet, whose linearisations do at first;
// - Rosenbrock's function in 3 iterations;
// - f = x₁, unbounded below, in 50;
// - bounds crossed, a start that isn't finite, a tol of 0;
// - a problem that can't be evaluated at its start, and one whose
//   gradients fail.
fn test_refusals() {
	mut linear := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			c[0] = x[0] + x[1] - 3.0
			c[1] = 1.0 - x[0] - x[1]
			return x[0] * x[0] + x[1] * x[1]
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {
			g[0], g[1] = 2.0 * x[0], 2.0 * x[1]
			a[0], a[1], a[2], a[3] = 1.0, 1.0, -1.0, -1.0
		}
	}
	refused(mut linear, [0.0, 0], .infeasible,
		"no step lowers the constraints' violation, 2.0: the linearised constraints have no point in common, and the problem may have no feasible point near x",
		inequalities: 2
	)
	mut apart := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			c[0] = 1.0 - x[0] * x[0] - x[1] * x[1]
			c[1] = x[0] - 2.0
			return x[1]
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {
			g[0], g[1] = 0.0, 1.0
			a[0], a[1], a[2], a[3] = -2.0 * x[0], -2.0 * x[1], 1.0, 0.0
		}
	}
	refused(mut apart, [0.5, 0.5], .infeasible, 'no feasible point near x', inequalities: 2)
	mut rosen := Hs{
		number: 1
	}
	refused(mut rosen, [-2.0, 1], .too_many_iterations, 'in 3 iterations', max_iter: 3)
	if _ := sqp(mut rosen, [-2.0, 1], max_iter: 3) {
		assert false
	} else {
		e := err as OptimizeError
		assert e.history.len == 3 && e.history[0].x == [-2.0, 1] && e.x != e.history[0].x
	}
	mut down := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			return x[0]
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {
			g[0] = 1.0
		}
	}
	refused(mut down, [0.0], .too_many_iterations, 'in 50 iterations', max_iter: 50)
	refused(mut rosen, [0.0, 0], .bad_input, 'bounds', lower: [1.0, 0], upper: [0.0, 1])
	refused(mut rosen, [math.nan(), 0], .bad_input, 'not finite')
	refused(mut rosen, [0.0, 0], .bad_input, 'tol', tol: 0.0)
	refused(mut rosen, [0.0, 0], .bad_input, 'one a variable', lower: [0.0])
	mut broken := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			return error('no value at x')
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {}
	}
	refused(mut broken, [1.0], .evaluation, 'at the start: no value at x')
	mut no_gradient := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			return x[0] * x[0]
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {
			return error('no gradient')
		}
	}
	refused(mut no_gradient, [1.0], .evaluation, 'no gradient')
}

// A problem undefined beyond part of its bounds' box: x − ln x is defined
// for x > 0 only, and its minimum is at 1. The first step from 8 overshoots
// below zero, where it fails, and the line search shortens it.
fn test_a_step_outside_the_domain() {
	mut p := Fns{
		v: fn (x []f64, mut c []f64) !f64 {
			if x[0] <= 0.0 {
				return error('ln of ${x[0]}')
			}
			return x[0] - math.log(x[0])
		}
		g: fn (x []f64, mut g []f64, mut a []f64) ! {
			g[0] = 1.0 - 1.0 / x[0]
		}
	}
	s := sqp(mut p, [8.0])!
	assert math.abs(s.x[0] - 1.0) <= 1e-7 && math.abs(s.f - 1.0) <= 1e-14
	assert s.history.any(it.step < 1.0)
}
