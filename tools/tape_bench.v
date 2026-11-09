// tape_bench times the expression tape against the same rates written in V
// (M9 step 3, criterion 3):
//
//     v -prod -cc gcc run tools/tape_bench.v
//
// Xu and Froment's three steam reforming rates, written as a user writes
// them, with definitions and units, compiled by expr into one tape; and
// written out by hand in SI, with their gradients. It checks the two agree,
// then times both at 1,000 reforming states, the best of nine rounds in
// which they take turns, so that each sees the same clock speeds.
module main

import expr
import math
import quantity
import rand
import time

const n_states = 1000
const calls_per_round = 200 * n_states
const rounds = 9
const r_gas = 8.31446261815324

fn context() expr.Context {
	species := ['CH4', 'H2O', 'H2', 'CO', 'CO2']
	mut keys := map[string]int{}
	for i, s in species {
		keys[s] = i
	}
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, species, keys},
		]
	}
}

fn tape() !expr.Tape {
	ctx := context()
	defs := expr.define(ctx, [
		expr.Source{'k1', '4.22e15 kmol*bar^0.5/(kg*h) * exp(-240.1 kJ/mol / (R*T))'},
		expr.Source{'k2', '1.96e6 kmol/(kg*h*bar) * exp(-67.13 kJ/mol / (R*T))'},
		expr.Source{'k3', '1.02e15 kmol*bar^0.5/(kg*h) * exp(-243.9 kJ/mol / (R*T))'},
		expr.Source{'K_CO', '8.23e-5 bar^-1 * exp(70.65 kJ/mol / (R*T))'},
		expr.Source{'K_H2', '6.12e-9 bar^-1 * exp(82.9 kJ/mol / (R*T))'},
		expr.Source{'K_CH4', '6.65e-4 bar^-1 * exp(38.28 kJ/mol / (R*T))'},
		expr.Source{'K_H2O', '1.77e5 * exp(-88.68 kJ/mol / (R*T))'},
		expr.Source{'K1', 'exp(30.42 - 27106 K / T) * 1 bar^2'},
		expr.Source{'K2', 'exp(-3.798 + 4160 K / T)'},
		expr.Source{'K3', 'exp(34.218 - 31266 K / T) * 1 bar^2'},
		expr.Source{'DEN', '1 + K_CO*p[CO] + K_H2*p[H2] + K_CH4*p[CH4] + K_H2O*p[H2O]/p[H2]'},
	])!
	mut rates := []expr.Expr{}
	for text in [
		'k1 / p[H2]^2.5 * (p[CH4]*p[H2O] - p[H2]^3*p[CO]/K1) / DEN^2',
		'k2 / p[H2] * (p[CO]*p[H2O] - p[H2]*p[CO2]/K2) / DEN^2',
		'k3 / p[H2]^3.5 * (p[CH4]*p[H2O]^2 - p[H2]^4*p[CO2]/K3) / DEN^2',
	] {
		rates << expr.parse(ctx, &defs, text)!
	}
	return expr.compile(ctx, &defs, rates)
}

// The constants in SI: rates in mol/(kg·s), pressures in Pa.
const per_hour = 1e3 / 3600.0
const a1 = 4.22e15 * per_hour * math.sqrt(1e5)
const a2 = 1.96e6 * per_hour / 1e5
const a3 = 1.02e15 * per_hour * math.sqrt(1e5)

// hand is the three rates and their gradients, written out: x is T, P and
// the partial pressures of CH4, H2O, H2, CO and CO2; g[k*7 + j] is
// d r[k]/d x[j].
@[direct_array_access]
fn hand(x []f64, mut r []f64, mut g []f64, grad bool) {
	temp, ch4, h2o, h2, co, co2 := x[0], x[2], x[3], x[4], x[5], x[6]
	rt := r_gas * temp
	dt := 1.0 / (rt * temp)
	t2 := 1.0 / (temp * temp)
	kco := 8.23e-10 * math.exp(70.65e3 / rt)
	kh2 := 6.12e-14 * math.exp(82.9e3 / rt)
	kch4 := 6.65e-9 * math.exp(38.28e3 / rt)
	kh2o := 1.77e5 * math.exp(-88.68e3 / rt)
	den := 1.0 + kco * co + kh2 * h2 + kch4 * ch4 + kh2o * h2o / h2
	den2 := den * den
	h2_2 := h2 * h2
	h2_3 := h2_2 * h2
	h2_4 := h2_2 * h2_2
	keq1 := math.exp(30.42 - 27106.0 / temp) * 1e10
	keq2 := math.exp(-3.798 + 4160.0 / temp)
	keq3 := math.exp(34.218 - 31266.0 / temp) * 1e10
	q1 := a1 * math.exp(-240.1e3 / rt) / math.pow(h2, 2.5) / den2
	q2 := a2 * math.exp(-67.13e3 / rt) / h2 / den2
	q3 := a3 * math.exp(-243.9e3 / rt) / math.pow(h2, 3.5) / den2
	r1 := q1 * (ch4 * h2o - h2_3 * co / keq1)
	r2 := q2 * (co * h2o - h2 * co2 / keq2)
	r3 := q3 * (ch4 * h2o * h2o - h2_4 * co2 / keq3)
	r[0], r[1], r[2] = r1, r2, r3
	if !grad {
		return
	}
	c := -2.0 / den
	d_ch4 := c * kch4
	d_h2o := c * kh2o / h2
	d_h2 := c * (kh2 - kh2o * h2o / h2_2)
	d_co := c * kco
	d_t := c * dt * (-70.65e3 * kco * co - 82.9e3 * kh2 * h2 - 38.28e3 * kch4 * ch4 + 88.68e3 * kh2o * h2o / h2)
	g[0] = q1 * h2_3 * co / keq1 * 27106.0 * t2 + r1 * (240.1e3 * dt) + r1 * d_t
	g[1] = 0.0
	g[2] = q1 * h2o + r1 * d_ch4
	g[3] = q1 * ch4 + r1 * d_h2o
	g[4] = -q1 * 3.0 * h2_2 * co / keq1 - 2.5 * r1 / h2 + r1 * d_h2
	g[5] = -q1 * h2_3 / keq1 + r1 * d_co
	g[6] = 0.0
	g[7] = -q2 * h2 * co2 / keq2 * 4160.0 * t2 + r2 * (67.13e3 * dt) + r2 * d_t
	g[8] = 0.0
	g[9] = r2 * d_ch4
	g[10] = q2 * co + r2 * d_h2o
	g[11] = -q2 * co2 / keq2 - r2 / h2 + r2 * d_h2
	g[12] = q2 * h2o + r2 * d_co
	g[13] = -q2 * h2 / keq2
	g[14] = q3 * h2_4 * co2 / keq3 * 31266.0 * t2 + r3 * (243.9e3 * dt) + r3 * d_t
	g[15] = 0.0
	g[16] = q3 * h2o * h2o + r3 * d_ch4
	g[17] = q3 * 2.0 * ch4 * h2o + r3 * d_h2o
	g[18] = -q3 * 4.0 * h2_3 * co2 / keq3 - 3.5 * r3 / h2 + r3 * d_h2
	g[19] = r3 * d_co
	g[20] = -q3 * h2_4 / keq3
}

// states are reforming conditions: 0.5 to 40 bar, shares of each reactant
// and product, and 800 to 1300 K.
fn states() []f64 {
	rand.seed([u32(2026), 1110])
	mut out := []f64{cap: n_states * 7}
	for _ in 0 .. n_states {
		p := (0.5 + 39.5 * rand.f64()) * 1e5
		y := [0.05 + rand.f64(), 0.05 + 2 * rand.f64(), 0.01 + rand.f64(), 1e-4 + 0.3 * rand.f64(),
			1e-4 + 0.2 * rand.f64()]
		s := y[0] + y[1] + y[2] + y[3] + y[4]
		out << 800.0 + 500.0 * rand.f64()
		out << p
		for j in 0 .. 5 {
			out << p * y[j] / s
		}
	}
	return out
}

struct Timing {
mut:
	ns  f64 = 1e300
	sum f64
}

fn main() {
	t := tape()!
	xs := states()
	mut buf := t.buffer()
	mut g := []f64{len: 3 * t.inputs}
	mut r := []f64{len: 3}
	mut h := []f64{len: 21}
	mut x := []f64{len: 7}
	// The two agree: rates to 1e-12, and each gradient row to 1e-10 of its
	// largest term, |∂r/∂x_j · x_j|.
	mut worst_r, mut worst_g := 0.0, 0.0
	for s in 0 .. n_states {
		for j in 0 .. 7 {
			x[j] = xs[s * 7 + j]
		}
		t.eval(x, mut buf)
		t.gradient(mut buf, mut g)
		hand(x, mut r, mut h, true)
		for k in 0 .. 3 {
			worst_r = math.max(worst_r, math.abs(t.output(&buf, k) - r[k]) / math.abs(r[k]))
			mut row := 0.0
			for j in 0 .. 7 {
				row = math.max(row, math.abs(h[k * 7 + j] * x[j]))
			}
			for j in 0 .. 7 {
				worst_g = math.max(worst_g, math.abs(g[k * 7 + j] - h[k * 7 + j]) * x[j] / row)
			}
		}
	}
	println('The tape: ${t.instrs.len} operations for three rates (${t.merged} merged, ${t.folded} folded, ${t.fused} fused)')
	println('Against the rates written out: rates within ${worst_r:.1e}, gradients within ${worst_g:.1e}')
	if worst_r > 1e-12 || worst_g > 1e-10 {
		eprintln('they disagree')
		exit(1)
	}
	for grad in [false, true] {
		mut tape_t, mut hand_t := Timing{}, Timing{}
		for _ in 0 .. rounds {
			t0 := time.sys_mono_now()
			for c := 0; c < calls_per_round; c += n_states {
				for s in 0 .. n_states {
					for j in 0 .. 7 {
						x[j] = xs[s * 7 + j]
					}
					t.eval(x, mut buf)
					if grad {
						t.gradient(mut buf, mut g)
					}
					tape_t.sum += t.output(&buf, 0) + g[0]
				}
			}
			tape_t.ns = math.min(tape_t.ns, f64(time.sys_mono_now() - t0) / calls_per_round)
			t1 := time.sys_mono_now()
			for c := 0; c < calls_per_round; c += n_states {
				for s in 0 .. n_states {
					for j in 0 .. 7 {
						x[j] = xs[s * 7 + j]
					}
					hand(x, mut r, mut h, grad)
					hand_t.sum += r[0] + h[0]
				}
			}
			hand_t.ns = math.min(hand_t.ns, f64(time.sys_mono_now() - t1) / calls_per_round)
		}
		what := if grad { 'rates and gradients' } else { 'rates' }
		println('${what}: the tape ${tape_t.ns:.1f} ns, by hand ${hand_t.ns:.1f} ns, ${tape_t.ns / hand_t.ns:.2f} times (checksums ${tape_t.sum:.6e}, ${hand_t.sum:.6e})')
	}
}
