// M9 step 0 spike: is a tape as fast as rates written out in V?
//
//     v -prod -cc gcc run spike/tape
//
// Xu and Froment's three methane steam reforming rates (AIChE J. 35, 88,
// 1989), with their one shared adsorption term, and the rates' gradients
// with respect to the five partial pressures and T, computed at the same
// states:
//
// - by a tape, as M9's D5 plans: the expressions built once into a flat
//   array of operations on slots, constants folded and repeats merged as
//   it's built, evaluated into a buffer the caller owns, and differentiated
//   by one reverse sweep for all three rates at once. Three variants: every
//   node evaluated every time; the constants put in the buffer once; and
//   fused operations, with a constant operand held in the operation (x + c,
//   c·x, c/x) and a·exp(b/T) as one operation, which is what Arrhenius and
//   van 't Hoff terms are;
// - written out by hand, the way the structured forms would be evaluated if
//   the tape were too slow (the plan's fallback).
//
// It checks that they agree, the gradients against central differences,
// and that none allocates, then times them. Not part of the simulator.
module main

import math
import rand
import time

const r_gas = 8.31446261815324 // J/(mol·K), the 2019 SI's exact value
const n_vars = 6 // p_CH4, p_H2O, p_H2, p_CO, p_CO2 (bar), then T (K)
const n_rates = 3
const n_states = 1000
const calls_per_round = 200 * n_states
const rounds = 9
const names = ['p_CH4', 'p_H2O', 'p_H2', 'p_CO', 'p_CO2', 'T']!

enum Op as u8 {
	variable // the caller's variable number a
	constant // c
	add      // a + b
	sub      // a − b
	mul      // a·b
	div      // a/b
	exp      // exp(a)
	powi     // a^c, c a small integer
	pow      // a^c
	addc     // a + c
	mulc     // c·a
	cdiv     // c/a
	aexp     // c·exp(d/a): Arrhenius and van 't Hoff, a being T
}

struct Node {
	op Op
	a  int
	b  int
	c  f64
	d  f64
}

// Tape is a flat array of operations on slots, one slot per node, in an
// order where every operand comes before its use. After finish, the
// constants come first, then the variables, then the operations from first.
struct Tape {
mut:
	nodes   []Node
	seen    map[string]int // each node's key, to merge repeats
	asked   int            // nodes asked for
	merged  int
	folded  int
	n_const int
	first   int
}

fn (mut t Tape) node(op Op, a int, b int, c f64, d f64) int {
	t.asked++
	mut x, mut y := a, b
	if (op == .add || op == .mul) && x > y {
		x, y = y, x
	}
	if op != .variable && op != .constant && t.nodes[x].op == .constant
		&& (y < 0 || t.nodes[y].op == .constant) {
		t.folded++
		return t.constant(fold(op, t.nodes[x].c, if y >= 0 { t.nodes[y].c } else { 0.0 }, c, d))
	}
	key := '${int(op)} ${x} ${y} ${math.f64_bits(c)} ${math.f64_bits(d)}'
	if key in t.seen {
		t.merged++
		return t.seen[key]
	}
	t.nodes << Node{op, x, y, c, d}
	t.seen[key] = t.nodes.len - 1
	return t.nodes.len - 1
}

fn fold(op Op, x f64, y f64, c f64, d f64) f64 {
	return match op {
		.add { x + y }
		.sub { x - y }
		.mul { x * y }
		.div { x / y }
		.exp { math.exp(x) }
		.powi { powi(x, int(c)) }
		.pow { math.pow(x, c) }
		.addc { x + c }
		.mulc { c * x }
		.cdiv { c / x }
		.aexp { c * math.exp(d / x) }
		else { math.nan() }
	}
}

fn (mut t Tape) variable(i int) int {
	return t.node(.variable, i, -1, 0, 0)
}

fn (mut t Tape) constant(c f64) int {
	return t.node(.constant, -1, -1, c, 0)
}

fn (mut t Tape) add(a int, b int) int {
	return t.node(.add, a, b, 0, 0)
}

fn (mut t Tape) sub(a int, b int) int {
	return t.node(.sub, a, b, 0, 0)
}

fn (mut t Tape) mul(a int, b int) int {
	return t.node(.mul, a, b, 0, 0)
}

fn (mut t Tape) div(a int, b int) int {
	return t.node(.div, a, b, 0, 0)
}

fn (mut t Tape) exp(a int) int {
	return t.node(.exp, a, -1, 0, 0)
}

fn (mut t Tape) powi(a int, n int) int {
	return t.node(.powi, a, -1, f64(n), 0)
}

fn (mut t Tape) pow(a int, c f64) int {
	return t.node(.pow, a, -1, c, 0)
}

// finish puts the constants first and the variables next, so that an
// evaluation can start after them with the constants already in its
// buffer, and gives the outputs' new slots.
fn (mut t Tape) finish(outs []int) []int {
	mut order := []int{}
	for i, n in t.nodes {
		if n.op == .constant {
			order << i
		}
	}
	t.n_const = order.len
	for i, n in t.nodes {
		if n.op == .variable {
			order << i
		}
	}
	t.first = order.len
	for i, n in t.nodes {
		if n.op != .constant && n.op != .variable {
			order << i
		}
	}
	mut slot := []int{len: t.nodes.len}
	for j, i in order {
		slot[i] = j
	}
	mut nodes := []Node{cap: t.nodes.len}
	for i in order {
		n := t.nodes[i]
		a := if n.op == .variable || n.a < 0 { n.a } else { slot[n.a] }
		nodes << Node{n.op, a, if n.b < 0 { -1 } else { slot[n.b] }, n.c, n.d}
	}
	t.nodes = nodes
	t.seen.clear()
	return outs.map(slot[it])
}

// prepare puts the constants into an evaluation's buffer.
fn (t &Tape) prepare(mut v []f64) {
	for i in 0 .. t.n_const {
		v[i] = t.nodes[i].c
	}
}

@[inline]
fn powi(x f64, n int) f64 {
	return match n {
		0 { 1.0 }
		1 { x }
		2 { x * x }
		3 { x * x * x }
		4 { (x * x) * (x * x) }
		else { math.pow(x, f64(n)) }
	}
}

// eval writes the slots from `from` on into v, from the variables x: from 0
// to evaluate every node, or from n_const if prepare has filled the
// constants.
@[direct_array_access]
fn (t &Tape) eval(x []f64, mut v []f64, from int) {
	for i := from; i < t.nodes.len; i++ {
		n := t.nodes[i]
		v[i] = match n.op {
			.variable { x[n.a] }
			.constant { n.c }
			.add { v[n.a] + v[n.b] }
			.sub { v[n.a] - v[n.b] }
			.mul { v[n.a] * v[n.b] }
			.div { v[n.a] / v[n.b] }
			.exp { math.exp(v[n.a]) }
			.powi { powi(v[n.a], int(n.c)) }
			.pow { math.pow(v[n.a], n.c) }
			.addc { v[n.a] + n.c }
			.mulc { n.c * v[n.a] }
			.cdiv { n.c / v[n.a] }
			.aexp { n.c * math.exp(n.d / v[n.a]) }
		}
	}
}

// gradient sweeps the tape backwards once from all the outputs: adj holds m
// adjoints per slot, and g[k*n_vars + j] gets d out[k]/d x[j]. v is eval's.
@[direct_array_access]
fn (t &Tape) gradient(v []f64, outs []int, mut adj []f64, mut g []f64) {
	m := outs.len
	for i := t.n_const * m; i < adj.len; i++ {
		adj[i] = 0.0
	}
	for k, o in outs {
		adj[o * m + k] = 1.0
	}
	for i := t.nodes.len - 1; i >= t.first; i-- {
		n := t.nodes[i]
		w := i * m
		a := n.a * m
		match n.op {
			.add {
				b := n.b * m
				for k in 0 .. m {
					adj[a + k] += adj[w + k]
					adj[b + k] += adj[w + k]
				}
			}
			.sub {
				b := n.b * m
				for k in 0 .. m {
					adj[a + k] += adj[w + k]
					adj[b + k] -= adj[w + k]
				}
			}
			.mul {
				b := n.b * m
				va, vb := v[n.a], v[n.b]
				for k in 0 .. m {
					adj[a + k] += adj[w + k] * vb
					adj[b + k] += adj[w + k] * va
				}
			}
			.div {
				b := n.b * m
				q := 1.0 / v[n.b]
				dq := v[i] * q
				for k in 0 .. m {
					adj[a + k] += adj[w + k] * q
					adj[b + k] -= adj[w + k] * dq
				}
			}
			else {
				// One operand: d is d v[i]/d v[a].
				x := v[n.a]
				d := match n.op {
					.exp { v[i] }
					.powi { n.c * powi(x, int(n.c) - 1) }
					.pow {
						if x != 0.0 { n.c * v[i] / x } else { n.c * math.pow(x, n.c - 1) }
					}
					.addc { 1.0 }
					.mulc { n.c }
					.cdiv { -v[i] / x }
					.aexp { -v[i] * n.d / (x * x) }
					else { 0.0 }
				}
				for k in 0 .. m {
					adj[a + k] += adj[w + k] * d
				}
			}
		}
	}
	for i in t.n_const .. t.first {
		for k in 0 .. m {
			g[k * n_vars + t.nodes[i].a] = adj[i * m + k]
		}
	}
}

// Builder makes the rates' terms either from the plain operations, as a
// parsed expression gives them, or fused.
struct Builder {
	fused bool
mut:
	t    Tape
	temp int
	rt   int
}

// arrhenius is a·exp(−e/RT): e is an activation energy, or a heat of
// adsorption for an adsorption constant.
fn (mut b Builder) arrhenius(a f64, e f64) int {
	if b.fused {
		return b.t.node(.aexp, b.temp, -1, a, -e / r_gas)
	}
	return b.t.mul(b.t.constant(a), b.t.exp(b.t.div(b.t.constant(-e), b.rt)))
}

// vant_hoff is exp(a + b/T).
fn (mut b Builder) vant_hoff(a f64, c f64) int {
	if b.fused {
		return b.t.node(.aexp, b.temp, -1, math.exp(a), c)
	}
	return b.t.exp(b.t.add(b.t.constant(a), b.t.div(b.t.constant(c), b.temp)))
}

fn (mut b Builder) one_plus(x int) int {
	if b.fused {
		return b.t.node(.addc, x, -1, 1.0, 0)
	}
	return b.t.add(b.t.constant(1), x)
}

// build makes Xu and Froment's three rates, each as it would be read from
// its own expression: every rate builds the adsorption term again, and the
// tape merges the copies. Rates in kmol/(kg_cat·h), pressures in bar.
fn build(fused bool) (Tape, []int) {
	mut b := Builder{
		fused: fused
	}
	ch4, h2o, h2, co, co2 := b.t.variable(0), b.t.variable(1), b.t.variable(2), b.t.variable(3), b.t.variable(4)
	b.temp = b.t.variable(5)
	mut outs := []int{}
	for rate in 0 .. n_rates {
		if !fused {
			b.rt = b.t.mul(b.t.constant(r_gas), b.temp)
		}
		kco := b.arrhenius(8.23e-5, -70.65e3)
		kh2 := b.arrhenius(6.12e-9, -82.9e3)
		kch4 := b.arrhenius(6.65e-4, -38.28e3)
		kh2o := b.arrhenius(1.77e5, 88.68e3)
		den := b.one_plus(b.t.add(b.t.add(b.t.add(b.t.mul(kco, co), b.t.mul(kh2, h2)), b.t.mul(kch4,
			ch4)), b.t.div(b.t.mul(kh2o, h2o), h2)))
		den2 := b.t.powi(den, 2)
		r := match rate {
			0 {
				// CH4 + H2O = CO + 3 H2
				k := b.arrhenius(4.22e15, 240.1e3)
				keq := b.vant_hoff(30.42, -27106.0)
				drive := b.t.sub(b.t.mul(ch4, h2o), b.t.div(b.t.mul(b.t.powi(h2, 3), co), keq))
				b.t.div(b.t.mul(b.t.div(k, b.t.pow(h2, 2.5)), drive), den2)
			}
			1 {
				// CO + H2O = CO2 + H2
				k := b.arrhenius(1.96e6, 67.13e3)
				keq := b.vant_hoff(-3.798, 4160.0)
				drive := b.t.sub(b.t.mul(co, h2o), b.t.div(b.t.mul(h2, co2), keq))
				b.t.div(b.t.mul(b.t.div(k, h2), drive), den2)
			}
			else {
				// CH4 + 2 H2O = CO2 + 4 H2
				k := b.arrhenius(1.02e15, 243.9e3)
				keq := b.vant_hoff(34.218, -31266.0)
				drive := b.t.sub(b.t.mul(ch4, b.t.powi(h2o, 2)), b.t.div(b.t.mul(b.t.powi(h2,
					4), co2), keq))
				b.t.div(b.t.mul(b.t.div(k, b.t.pow(h2, 3.5)), drive), den2)
			}
		}
		outs << r
	}
	slots := b.t.finish(outs)
	return b.t, slots
}

// hand is the same three rates and their gradients, written out: r[k], and
// g[k*n_vars + j] = d r[k]/d x[j].
@[direct_array_access]
fn hand(x []f64, mut r []f64, mut g []f64) {
	ch4, h2o, h2, co, co2, temp := x[0], x[1], x[2], x[3], x[4], x[5]
	rt := r_gas * temp
	dt := 1.0 / (rt * temp) // d(−e/RT)/dT = e/(RT²)
	t2 := 1.0 / (temp * temp)
	kco := 8.23e-5 * math.exp(70.65e3 / rt)
	kh2 := 6.12e-9 * math.exp(82.9e3 / rt)
	kch4 := 6.65e-4 * math.exp(38.28e3 / rt)
	kh2o := 1.77e5 * math.exp(-88.68e3 / rt)
	den := 1.0 + kco * co + kh2 * h2 + kch4 * ch4 + kh2o * h2o / h2
	den2 := den * den
	// The adsorption term's derivatives, as d ln(1/den²).
	c := -2.0 / den
	d_ch4 := c * kch4
	d_h2o := c * kh2o / h2
	d_h2 := c * (kh2 - kh2o * h2o / (h2 * h2))
	d_co := c * kco
	d_t := c * dt * (-70.65e3 * kco * co - 82.9e3 * kh2 * h2 - 38.28e3 * kch4 * ch4 + 88.68e3 * kh2o * h2o / h2)
	h2_2 := h2 * h2
	h2_3 := h2_2 * h2
	h2_4 := h2_2 * h2_2
	// Steam reforming.
	k1 := 4.22e15 * math.exp(-240.1e3 / rt)
	keq1 := math.exp(30.42 + -27106.0 / temp)
	q1 := k1 / math.pow(h2, 2.5) / den2
	r1 := q1 * (ch4 * h2o - h2_3 * co / keq1)
	r[0] = r1
	g[0] = q1 * h2o + r1 * d_ch4
	g[1] = q1 * ch4 + r1 * d_h2o
	g[2] = -q1 * 3.0 * h2_2 * co / keq1 - 2.5 * r1 / h2 + r1 * d_h2
	g[3] = -q1 * h2_3 / keq1 + r1 * d_co
	g[4] = 0.0
	g[5] = q1 * h2_3 * co / keq1 * 27106.0 * t2 + r1 * (240.1e3 * dt) + r1 * d_t
	// The water-gas shift.
	k2 := 1.96e6 * math.exp(-67.13e3 / rt)
	keq2 := math.exp(-3.798 + 4160.0 / temp)
	q2 := k2 / h2 / den2
	r2 := q2 * (co * h2o - h2 * co2 / keq2)
	r[1] = r2
	g[6] = r2 * d_ch4
	g[7] = q2 * co + r2 * d_h2o
	g[8] = -q2 * co2 / keq2 - r2 / h2 + r2 * d_h2
	g[9] = q2 * h2o + r2 * d_co
	g[10] = -q2 * h2 / keq2
	g[11] = -q2 * h2 * co2 / keq2 * 4160.0 * t2 + r2 * (67.13e3 * dt) + r2 * d_t
	// The overall reaction.
	k3 := 1.02e15 * math.exp(-243.9e3 / rt)
	keq3 := math.exp(34.218 + -31266.0 / temp)
	q3 := k3 / math.pow(h2, 3.5) / den2
	r3 := q3 * (ch4 * h2o * h2o - h2_4 * co2 / keq3)
	r[2] = r3
	g[12] = q3 * h2o * h2o + r3 * d_ch4
	g[13] = q3 * 2.0 * ch4 * h2o + r3 * d_h2o
	g[14] = -q3 * 4.0 * h2_3 * co2 / keq3 - 3.5 * r3 / h2 + r3 * d_h2
	g[15] = r3 * d_co
	g[16] = -q3 * h2_4 / keq3
	g[17] = q3 * h2_4 * co2 / keq3 * 31266.0 * t2 + r3 * (243.9e3 * dt) + r3 * d_t
}

// hand_rates is hand's rates alone, for timing evaluation without gradients.
@[direct_array_access]
fn hand_rates(x []f64, mut r []f64) {
	ch4, h2o, h2, co, co2, temp := x[0], x[1], x[2], x[3], x[4], x[5]
	rt := r_gas * temp
	kco := 8.23e-5 * math.exp(70.65e3 / rt)
	kh2 := 6.12e-9 * math.exp(82.9e3 / rt)
	kch4 := 6.65e-4 * math.exp(38.28e3 / rt)
	kh2o := 1.77e5 * math.exp(-88.68e3 / rt)
	den := 1.0 + kco * co + kh2 * h2 + kch4 * ch4 + kh2o * h2o / h2
	den2 := den * den
	h2_3 := h2 * h2 * h2
	r[0] = 4.22e15 * math.exp(-240.1e3 / rt) / math.pow(h2, 2.5) / den2 * (ch4 * h2o - h2_3 * co / math.exp(30.42 +
		-27106.0 / temp))
	r[1] = 1.96e6 * math.exp(-67.13e3 / rt) / h2 / den2 * (co * h2o - h2 * co2 / math.exp(-3.798 + 4160.0 / temp))
	r[2] = 1.02e15 * math.exp(-243.9e3 / rt) / math.pow(h2, 3.5) / den2 * (ch4 * h2o * h2o - h2_3 * h2 * co2 / math.exp(34.218 +
		-31266.0 / temp))
}

// states are reforming conditions: 0.5 to 40 bar, shares of each reactant
// and product, and 800 to 1300 K.
fn states() []f64 {
	rand.seed([u32(2026), 1109])
	mut out := []f64{cap: n_states * n_vars}
	for _ in 0 .. n_states {
		p := 0.5 + 39.5 * rand.f64()
		y := [0.05 + rand.f64(), 0.05 + 2 * rand.f64(), 0.01 + rand.f64(), 1e-4 + 0.3 * rand.f64(),
			1e-4 + 0.2 * rand.f64()]
		s := y[0] + y[1] + y[2] + y[3] + y[4]
		for j in 0 .. 5 {
			out << p * y[j] / s
		}
		out << 800.0 + 500.0 * rand.f64()
	}
	return out
}

// scale is a gradient row's size: the largest |d r/d x_j · x_j|, so that
// an error is measured against the terms the gradient is made of.
fn scale(g []f64, x []f64, k int) f64 {
	mut s := 0.0
	for j in 0 .. n_vars {
		s = math.max(s, math.abs(g[k * n_vars + j] * x[j]))
	}
	return s
}

// Work is an evaluation's buffers, owned by the caller.
struct Work {
mut:
	x   []f64 = []f64{len: n_vars}
	v   []f64
	adj []f64
	r   []f64 = []f64{len: n_rates}
	g   []f64 = []f64{len: n_rates * n_vars}
}

fn Work.new(t &Tape) Work {
	mut w := Work{
		v:   []f64{len: t.nodes.len}
		adj: []f64{len: t.nodes.len * n_rates}
	}
	t.prepare(mut w.v)
	return w
}

fn (mut w Work) set(xs []f64, s int) {
	for j in 0 .. n_vars {
		w.x[j] = xs[s * n_vars + j]
	}
}

fn (mut w Work) run(t &Tape, outs []int, from int, grad bool) {
	t.eval(w.x, mut w.v, from)
	for k, o in outs {
		w.r[k] = w.v[o]
	}
	if grad {
		t.gradient(w.v, outs, mut w.adj, mut w.g)
	}
}

struct Timing {
mut:
	ns    f64 = 1e300
	bytes u64
	sum   f64
}

// round_tape times one round of calls of a tape, adding to its Timing: the
// best round's time per call, the bytes allocated in all rounds, and a
// checksum that keeps the work from being optimised away.
fn round_tape(t &Tape, outs []int, from int, grad bool, xs []f64, mut w Work, mut tm Timing) {
	before := gc_heap_usage().total_bytes
	t0 := time.sys_mono_now()
	for c := 0; c < calls_per_round; c += n_states {
		for s in 0 .. n_states {
			w.set(xs, s)
			w.run(t, outs, from, grad)
			tm.sum += w.r[0] + w.g[5]
		}
	}
	tm.ns = math.min(tm.ns, f64(time.sys_mono_now() - t0) / calls_per_round)
	tm.bytes += u64(gc_heap_usage().total_bytes - before)
}

fn round_hand(grad bool, xs []f64, mut w Work, mut tm Timing) {
	before := gc_heap_usage().total_bytes
	t0 := time.sys_mono_now()
	for c := 0; c < calls_per_round; c += n_states {
		for s in 0 .. n_states {
			w.set(xs, s)
			if grad {
				hand(w.x, mut w.r, mut w.g)
			} else {
				hand_rates(w.x, mut w.r)
			}
			tm.sum += w.r[0] + w.g[5]
		}
	}
	tm.ns = math.min(tm.ns, f64(time.sys_mono_now() - t0) / calls_per_round)
	tm.bytes += u64(gc_heap_usage().total_bytes - before)
}

// check holds a tape against hand: the largest relative difference in the
// rates, and in the gradients against the largest term in their row.
fn check(t &Tape, outs []int, xs []f64) (f64, f64) {
	mut w := Work.new(t)
	mut h := Work{}
	mut worst_r, mut worst_g := 0.0, 0.0
	for s in 0 .. n_states {
		w.set(xs, s)
		w.run(t, outs, t.n_const, true)
		h.set(xs, s)
		hand(h.x, mut h.r, mut h.g)
		for k in 0 .. n_rates {
			worst_r = math.max(worst_r, math.abs(w.r[k] - h.r[k]) / math.abs(h.r[k]))
			sc := scale(h.g, h.x, k)
			for j in 0 .. n_vars {
				worst_g = math.max(worst_g, math.abs(w.g[k * n_vars + j] - h.g[k * n_vars + j]) * h.x[j] / sc)
			}
		}
	}
	return worst_r, worst_g
}

// central holds hand's gradients against central differences of its rates,
// h = 1e-5·x: the largest error against the largest term in its row.
fn central(xs []f64) (f64, string) {
	mut w, mut p, mut m := Work{}, Work{}, Work{}
	mut worst := 0.0
	mut at := ''
	for s in 0 .. n_states {
		w.set(xs, s)
		hand(w.x, mut w.r, mut w.g)
		for j in 0 .. n_vars {
			p.set(xs, s)
			m.set(xs, s)
			p.x[j] += 1e-5 * w.x[j]
			m.x[j] -= 1e-5 * w.x[j]
			hand_rates(p.x, mut p.r)
			hand_rates(m.x, mut m.r)
			for k in 0 .. n_rates {
				fd := (p.r[k] - m.r[k]) / (p.x[j] - m.x[j])
				e := math.abs(fd - w.g[k * n_vars + j]) * w.x[j] / scale(w.g, w.x, k)
				if e > worst {
					worst = e
					at = 'rate ${k + 1}, d/d${names[j]}'
				}
			}
		}
	}
	return worst, at
}

fn main() {
	xs := states()
	sw := time.new_stopwatch()
	plain, plain_outs := build(false)
	built := sw.elapsed()
	fused, fused_outs := build(true)
	println('The plain tape: ${plain.nodes.len} nodes (${plain.n_const} constants, ${plain.first - plain.n_const} variables, ${plain.nodes.len - plain.first} operations) of ${plain.asked} asked for: ${plain.merged} merged as repeats, ${plain.folded} folded as constants; built in ${built.microseconds()} µs')
	println('The fused tape: ${fused.nodes.len} nodes (${fused.n_const} constants, ${fused.nodes.len - fused.first} operations) of ${fused.asked} asked for, ${fused.merged} merged')
	fd, fd_at := central(xs)
	println('At ${n_states} states, the hand-written gradients against central differences (h = 1e-5·x): within ${fd:.1e} of the largest term in their row, the worst ${fd_at}')
	pr, pg := check(&plain, plain_outs, xs)
	println('the plain tape against the hand-written: rates within ${pr:.1e} (relative), gradients within ${pg:.1e}')
	fr, fg := check(&fused, fused_outs, xs)
	println('the fused tape against the hand-written: rates within ${fr:.1e} (relative), gradients within ${fg:.1e}')
	// The variants take turns in each round, so that each sees the same
	// clock speeds; each keeps its best round.
	for grad in [false, true] {
		what := if grad { 'rates and gradients' } else { 'rates' }
		mut hd, mut every, mut once, mut fu := Timing{}, Timing{}, Timing{}, Timing{}
		mut wh, mut wp, mut wf := Work{}, Work.new(&plain), Work.new(&fused)
		for _ in 0 .. rounds {
			round_hand(grad, xs, mut wh, mut hd)
			round_tape(&plain, plain_outs, 0, grad, xs, mut wp, mut every)
			round_tape(&plain, plain_outs, plain.n_const, grad, xs, mut wp, mut once)
			round_tape(&fused, fused_outs, fused.n_const, grad, xs, mut wf, mut fu)
		}
		println('${what}, ns a call: by hand ${hd.ns:.1f}; the tape, every node ${every.ns:.1f} (${every.ns / hd.ns:.2f}x), constants once ${once.ns:.1f} (${once.ns / hd.ns:.2f}x), fused ${fu.ns:.1f} (${fu.ns / hd.ns:.2f}x)')
		println('  bytes allocated in ${rounds * calls_per_round} calls each: ${hd.bytes}, ${every.bytes}, ${once.bytes}, ${fu.bytes} (checksums ${hd.sum:.6e}, ${every.sum:.6e}, ${once.sum:.6e}, ${fu.sum:.6e})')
	}
}
