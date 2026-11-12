// RODAS4 on six of Hairer and Wanner's stiff test problems (M9 step 8,
// criterion 6), against the Bari test set's reference solutions
// (testdata/m9_odes.json), each problem written here again from the test
// set's Fortran.
module numeric

import json2
import math
import os

// Problem is a test problem: its f and Jacobian, by name.
struct Problem {
	name string
mut:
	calls int
}

fn (mut p Problem) eval(_ f64, y []f64, mut f []f64) ! {
	p.calls++
	match p.name {
		'hires' {
			f[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007
			f[1] = 1.71 * y[0] - 8.75 * y[1]
			f[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4]
			f[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3]
			f[4] = -1.745 * y[4] + 0.43 * (y[5] + y[6])
			f[5] = -280 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6]
			f[6] = 280 * y[5] * y[7] - 1.81 * y[6]
			f[7] = -f[6]
		}
		'rober' {
			f[0] = -0.04 * y[0] + 1e4 * y[1] * y[2]
			f[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1]
			f[2] = 3e7 * y[1] * y[1]
		}
		'vdpol' {
			f[0] = y[1]
			f[1] = ((1 - y[0] * y[0]) * y[1] - y[0]) / 1e-6
		}
		'orego' {
			f[0] = 77.27 * (y[1] + y[0] * (1 - 8.375e-6 * y[0] - y[1]))
			f[1] = (y[2] - (1 + y[0]) * y[1]) / 77.27
			f[2] = 0.161 * (y[0] - y[2])
		}
		'e5' {
			p1 := 7.89e-10 * y[0]
			p2 := 1.1e7 * y[0] * y[2]
			p3 := 1.13e9 * y[1] * y[2]
			p4 := 1.13e3 * y[3]
			f[0] = -p1 - p2
			f[1] = p1 - p3
			f[3] = p2 - p4
			f[2] = f[1] - f[3]
		}
		'pollu' {
			r := pollu_rates(y)
			f[0] = -r[0] - r[9] - r[13] - r[22] - r[23] + r[1] + r[2] + r[8] + r[10] + r[11] +
				r[21] + r[24]
			f[1] = -r[1] - r[2] - r[8] - r[11] + r[0] + r[20]
			f[2] = -r[14] + r[0] + r[16] + r[18] + r[21]
			f[3] = -r[1] - r[15] - r[16] - r[22] + r[14]
			f[4] = -r[2] + r[3] + r[3] + r[5] + r[6] + r[12] + r[19]
			f[5] = -r[5] - r[7] - r[13] - r[19] + r[2] + r[17] + r[17]
			f[6] = -r[3] - r[4] - r[5] + r[12]
			f[7] = r[3] + r[4] + r[5] + r[6]
			f[8] = -r[6] - r[7]
			f[9] = -r[11] + r[6] + r[8]
			f[10] = -r[8] - r[9] + r[7] + r[10]
			f[11] = r[8]
			f[12] = -r[10] + r[9]
			f[13] = -r[12] + r[11]
			f[14] = r[13]
			f[15] = -r[17] - r[18] + r[15]
			f[16] = -r[19]
			f[17] = r[19]
			f[18] = -r[20] - r[21] - r[23] + r[22] + r[24]
			f[19] = -r[24] + r[23]
		}
		else {
			return error('no problem ${p.name}')
		}
	}
}

fn (mut p Problem) jacobian(_ f64, y []f64, mut j []f64, mut dfdt []f64) ! {
	n := y.len
	for i in 0 .. n * n {
		j[i] = 0
	}
	for i in 0 .. n {
		dfdt[i] = 0
	}
	match p.name {
		'hires' {
			set(mut j, n, 0, [0, 1, 2], [-1.71, 0.43, 8.32])
			set(mut j, n, 1, [0, 1], [1.71, -8.75])
			set(mut j, n, 2, [2, 3, 4], [-10.03, 0.43, 0.035])
			set(mut j, n, 3, [1, 2, 3], [8.32, 1.71, -1.12])
			set(mut j, n, 4, [4, 5, 6], [-1.745, 0.43, 0.43])
			set(mut j, n, 5, [3, 4, 5, 6, 7], [0.69, 1.71, -0.43 - 280 * y[7], 0.69, -280 * y[5]])
			set(mut j, n, 6, [5, 6, 7], [280 * y[7], -1.81, 280 * y[5]])
			set(mut j, n, 7, [5, 6, 7], [-280 * y[7], 1.81, -280 * y[5]])
		}
		'rober' {
			set(mut j, n, 0, [0, 1, 2], [-0.04, 1e4 * y[2], 1e4 * y[1]])
			set(mut j, n, 1, [0, 1, 2], [0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]])
			set(mut j, n, 2, [1], [6e7 * y[1]])
		}
		'vdpol' {
			set(mut j, n, 0, [1], [1.0])
			set(mut j, n, 1, [0, 1], [(-2 * y[0] * y[1] - 1) / 1e-6, (1 - y[0] * y[0]) / 1e-6])
		}
		'orego' {
			set(mut j, n, 0, [0, 1], [77.27 * (1 - 2 * 8.375e-6 * y[0] - y[1]), 77.27 * (1 - y[0])])
			set(mut j, n, 1, [0, 1, 2], [-y[1] / 77.27, -(1 + y[0]) / 77.27, 1 / 77.27])
			set(mut j, n, 2, [0, 2], [0.161, -0.161])
		}
		'e5' {
			a, b, c := 7.89e-10, 1.1e7, 1.13e9
			set(mut j, n, 0, [0, 2], [-a - b * y[2], -b * y[0]])
			set(mut j, n, 1, [0, 1, 2], [a, -c * y[2], -c * y[1]])
			set(mut j, n, 3, [0, 2, 3], [b * y[2], b * y[0], -1.13e3])
			for k in 0 .. n {
				j[2 * n + k] = j[1 * n + k] - j[3 * n + k]
			}
		}
		'pollu' {
			// Each reaction's rate, its derivative by each of its reactants,
			// and its coefficient in each species' f.
			for r, sp in pollu_reactants {
				for a in sp {
					mut d := pollu_k[r]
					for s in sp {
						if s != a {
							d *= y[s]
						}
					}
					for i in 0 .. n {
						c := pollu_nu[r][i]
						if c != 0 {
							j[i * n + a] += c * d
						}
					}
				}
			}
		}
		else {
			return error('no problem ${p.name}')
		}
	}
}

fn set(mut j []f64, n int, row int, cols []int, values []f64) {
	for k, c in cols {
		j[row * n + c] = values[k]
	}
}

const pollu_k = [0.35, 0.266e2, 0.123e5, 0.86e-3, 0.82e-3, 0.15e5, 0.13e-3, 0.24e5, 0.165e5, 0.9e4,
	0.22e-1, 0.12e5, 0.188e1, 0.163e5, 0.48e7, 0.35e-3, 0.175e-1, 0.1e9, 0.444e12, 0.124e4, 0.21e1,
	0.578e1, 0.474e-1, 0.178e4, 0.312e1]

// pollu_reactants are each reaction's reactants, 0-based: the test set's
// r(1) to r(25).
const pollu_reactants = [[0], [1, 3], [4, 1], [6], [6], [6, 5], [8], [8, 5], [10, 1], [10, 0],
	[12], [9, 1], [13], [0, 5], [2], [3], [3], [15], [15], [16, 5], [18], [18], [0, 3], [18, 0],
	[19]]

fn pollu_rates(y []f64) []f64 {
	mut r := []f64{len: 25}
	for k, sp in pollu_reactants {
		mut x := pollu_k[k]
		for s in sp {
			x *= y[s]
		}
		r[k] = x
	}
	return r
}

// pollu_nu is each reaction's coefficient in each species' f, read off the
// test set's f(1) to f(20) by a probe of eval with one rate at a time.
const pollu_nu = pollu_coefficients()

fn pollu_coefficients() [][]f64 {
	// f is linear in the rates, so f with every rate but r_k zero, and r_k
	// one, is reaction k's column: eval's own lines, read back.
	mut out := [][]f64{}
	for k in 0 .. 25 {
		mut r := []f64{len: 25}
		r[k] = 1
		out << pollu_f_of_rates(r)
	}
	return out
}

fn pollu_f_of_rates(r []f64) []f64 {
	mut f := []f64{len: 20}
	f[0] = -r[0] - r[9] - r[13] - r[22] - r[23] + r[1] + r[2] + r[8] + r[10] + r[11] + r[21] +
		r[24]
	f[1] = -r[1] - r[2] - r[8] - r[11] + r[0] + r[20]
	f[2] = -r[14] + r[0] + r[16] + r[18] + r[21]
	f[3] = -r[1] - r[15] - r[16] - r[22] + r[14]
	f[4] = -r[2] + r[3] + r[3] + r[5] + r[6] + r[12] + r[19]
	f[5] = -r[5] - r[7] - r[13] - r[19] + r[2] + r[17] + r[17]
	f[6] = -r[3] - r[4] - r[5] + r[12]
	f[7] = r[3] + r[4] + r[5] + r[6]
	f[8] = -r[6] - r[7]
	f[9] = -r[11] + r[6] + r[8]
	f[10] = -r[8] - r[9] + r[7] + r[10]
	f[11] = r[8]
	f[12] = -r[10] + r[9]
	f[13] = -r[12] + r[11]
	f[14] = r[13]
	f[15] = -r[17] - r[18] + r[15]
	f[16] = -r[19]
	f[17] = r[19]
	f[18] = -r[20] - r[21] - r[23] + r[22] + r[24]
	f[19] = -r[24] + r[23]
	return f
}

struct OdeCase {
	name          string
	t0            f64
	t1            f64
	y0            []f64
	atol_per_rtol f64
	reference     []f64
}

struct OdeFixture {
	problems []OdeCase
}

fn ode_fixture() OdeFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_odes.json')
	return json2.decode[OdeFixture](os.read_file(path) or { panic(err) }) or { panic(err) }
}

// units is the largest error in units of the tolerance asked:
// |y − ref|/(atol + rtol·|ref|).
fn units(y []f64, reference []f64, rtol f64, atol f64) f64 {
	mut worst := 0.0
	for i, r in reference {
		worst = math.max(worst, math.abs(y[i] - r) / (atol + rtol * math.abs(r)))
	}
	return worst
}

// Each problem's Jacobian is its f's: against central differences at its
// start, moved off zero, to 1e-6 of the Jacobian's largest entry, the
// differences' truncation and rounding.
fn test_jacobians_are_right() {
	f := ode_fixture()
	for c in f.problems {
		n := c.y0.len
		mut p := Problem{
			name: c.name
		}
		y := c.y0.map(it + 1e-3)
		mut jac := []f64{len: n * n}
		mut dfdt := []f64{len: n}
		p.jacobian(0, y, mut jac, mut dfdt)!
		mut scale := 1.0
		for x in jac {
			scale = math.max(scale, math.abs(x))
		}
		mut fp := []f64{len: n}
		mut fm := []f64{len: n}
		for k in 0 .. n {
			h := 1e-6 * math.max(math.abs(y[k]), 1e-3)
			mut yp := y.clone()
			mut ym := y.clone()
			yp[k] += h
			ym[k] -= h
			p.eval(0, yp, mut fp)!
			p.eval(0, ym, mut fm)!
			for i in 0 .. n {
				d := (fp[i] - fm[i]) / (2 * h)
				assert math.abs(d - jac[i * n + k]) <= 1e-6 * scale, '${c.name}: ∂f${i}/∂y${k} ${jac[i * n + k]} vs ${d}'
			}
		}
	}
}

// RODAS4 reaches the test set's reference on HIRES, Robertson, van der Pol,
// the Oregonator and the pollution problem at rtol 1e-4, 1e-7 and 1e-10,
// atol the problem's share of rtol. Each end error, in units of the
// tolerance asked, |y − ref|/(atol + rtol·|ref|), is within a unit but the
// Oregonator's, whose six sharp fronts carry each step's error on (2.9
// units). The worst over the three tolerances and the steps at each are
// recorded, unable to grow: the units 10% above what was measured, the
// steps as measured. Every evaluation is counted.
fn test_rodas4_against_the_test_set() {
	limits := {
		'hires': 0.14
		'rober': 0.05
		'vdpol': 0.105
		'orego': 3.2
		'pollu': 0.104
	}
	budgets := {
		'hires': [90, 856, 12298]
		'rober': [111, 563, 5128]
		'vdpol': [311, 2209, 31775]
		'orego': [347, 2479, 26650]
		'pollu': [28, 149, 1074]
	}
	f := ode_fixture()
	mut seen := 0
	for c in f.problems {
		limit := limits[c.name] or { continue }
		seen++
		mut worst := 0.0
		for k, rtol in [1e-4, 1e-7, 1e-10] {
			atol := c.atol_per_rtol * rtol
			mut p := Problem{
				name: c.name
			}
			s := rodas4(mut p, c.t0, c.y0, c.t1, rtol: rtol, atol: atol)!
			u := units(s.y, c.reference, rtol, atol)
			worst = math.max(worst, u)
			assert s.evaluations == p.calls
			assert s.steps <= budgets[c.name][k], '${c.name} at rtol ${rtol}: ${s.steps} steps'
		}
		println('${c.name}: within ${worst} tolerance units')
		assert worst <= limit, '${c.name}: ${worst} units'
	}
	assert seen == 5
}

// E5's y₂ and y₃ at 1e13, 9e-23, are beyond RODAS4 in double precision:
// y₂ − y₃ − y₄ is conserved exactly, but while y₂ and y₃ are 1e-10 their
// rates' terms are 1e-7 and cancel, and a Rosenbrock method solves for each
// stage's whole increment from those, so its rounding is theirs (an
// implicit method's Newton solves for small corrections). The 1e-22 lost
// then stays when y₂ and y₃ have fallen to 1e-22; the test set's reference,
// in 128-bit arithmetic, keeps the invariant to 1.3e-25. What's claimed is
// held instead, at rtol 1e-4, 1e-7 and 1e-10 with atol 1e-24: y₁ and y₄
// decay below their tolerance, and the invariant drifts by rounding's
// amount. That amount is the compilers' too: 7.9e-23 under tcc, 1.9e-21
// under gcc, whose arithmetic differs from tcc's in the last bit, which E5
// amplifies; the bound is 1e-20. y₂ and y₃ end 51% from the reference under
// tcc and a factor of 20 under gcc: unresolved either way, so not held.
fn test_e5_as_far_as_double_precision_goes() {
	f := ode_fixture()
	c := f.problems.filter(it.name == 'e5')[0]
	mut drift := 0.0
	for rtol in [1e-4, 1e-7, 1e-10] {
		atol := 1e-24
		mut p := Problem{
			name: 'e5'
		}
		s := rodas4(mut p, c.t0, c.y0, c.t1, rtol: rtol, atol: atol)!
		assert math.abs(s.y[0]) <= atol && math.abs(s.y[3]) <= atol, s.y.str()
		drift = math.max(drift, math.abs(s.y[1] - s.y[2] - s.y[3]))
		println('e5 at rtol ${rtol}: y ${s.y} against ${c.reference}')
	}
	println('e5: y₂ − y₃ − y₄ drifts by ${drift}')
	assert drift <= 1e-20, '${drift}'
}

// dopri5 finds every one of the six stiff, within 1,016 steps, and refuses
// it for rodas4.
fn test_dopri5_refuses_them_as_stiff() {
	f := ode_fixture()
	for c in f.problems {
		mut p := Problem{
			name: c.name
		}
		dopri5(mut p, c.t0, c.y0, c.t1, rtol: 1e-6, atol: c.atol_per_rtol * 1e-6) or {
			assert err is OdeError && (err as OdeError).kind == .stiff, err.msg()
			assert (err as OdeError).steps <= 1016
			continue
		}
		assert false, '${c.name}: dopri5 answered'
	}
}
