module expr

import math

// A tape is expressions compiled for speed (D5, as step 0 amended it): a
// flat array of operations on slots, one slot per operation, every operand
// before its use. Compiling folds constants, merges repeats (a definition
// used in many places is one slot), fuses a constant operand into its
// operation, and drops what no output uses. Evaluating walks the tape into
// a buffer the caller owns, allocating nothing; one reverse sweep then
// gives the gradient of every output with respect to every input at once.
//
// Every fusion keeps the tree's own arithmetic, so a tape gives the bits
// eval gives: x − c is x + (−c), exactly; c·exp(d/x) is the tree's
// A·exp(−E/(R·T)) with x the shared slot of R·T; exp(c + d/x) is its
// exp(30.42 − 27106 K/T).

// Op is one operation of a tape.
pub enum Op as u8 {
	input    // the caller's input a
	constant // c
	add      // a + b
	sub      // a − b
	mul      // a·b
	div      // a/b
	pow      // a^b, b varying
	neg      // −a
	exp      // exp(a)
	ln       // ln(a)
	log10    // log10(a)
	sqrt     // √a
	abs      // |a|
	powc     // a^c, c fixed
	addc     // a + c
	csub     // c − a
	mulc     // c·a
	divc     // a/c
	cdiv     // c/a
	aexp     // c·exp(d/a): an Arrhenius or van 't Hoff term, a being R·T or T
	expl     // exp(c + d/a): an equilibrium constant's correlation
}

// Instr is one operation, its operands' slots a and b, and its constants.
pub struct Instr {
pub:
	op Op
	a  int = -1
	b  int = -1
	c  f64
	d  f64
}

// Tape is compiled expressions: the constants come first, from slot 0, then
// the inputs used, from n_const, then the operations, from first. inputs is
// how many inputs an evaluation takes (the context's slots), and outputs
// each expression's slot, in the order given. merged, folded and fused
// count what compiling did.
pub struct Tape {
pub:
	instrs  []Instr
	n_const int
	first   int
	inputs  int
	outputs []int
	merged  int
	folded  int
	fused   int
}

// Buffer is an evaluation's workspace, owned by the caller: each slot's
// value, and its adjoints, one per output. Tape.buffer makes one, with the
// constants in place.
pub struct Buffer {
pub mut:
	v   []f64
	adj []f64
}

struct Builder {
	ctx  Context
	defs Definitions
mut:
	instrs  []Instr
	seen    map[string]int
	defs_at map[int]int
	merged  int
	folded  int
	fused   int
}

// compile compiles expressions, each bound and checked in ctx with defs,
// into one tape, so that what they share is computed once.
pub fn compile(ctx Context, defs &Definitions, exprs []Expr) Tape {
	mut b := Builder{
		ctx:  ctx
		defs: *defs
	}
	mut outs := []int{}
	for e in exprs {
		outs << b.tree(e)
	}
	return b.finish(ctx.slots(), outs)
}

// tree compiles one expression, its nodes in order, and gives its root's
// slot.
fn (mut b Builder) tree(e Expr) int {
	mut at := []int{len: e.nodes.len, init: -1}
	for i, n in e.nodes {
		at[i] = b.node(n, at)
	}
	return at[e.root]
}

fn (mut b Builder) node(n Node, at []int) int {
	return match n.kind {
		.num {
			b.constant(n.value)
		}
		.var, .member {
			b.emit(.input, n.index, -1, 0, 0)
		}
		.constant {
			b.constant(constants[n.index].value)
		}
		.def {
			// Once each: merging would make a definition one slot anyway, but
			// a chain of definitions each using the next twice would compile
			// 2ⁿ times without this.
			if slot := b.defs_at[n.index] {
				slot
			} else {
				s := b.tree(b.defs.exprs[n.index])
				b.defs_at[n.index] = s
				s
			}
		}
		.neg {
			b.emit(.neg, at[n.a], -1, 0, 0)
		}
		.add {
			b.sum(at[n.a], at[n.b])
		}
		.sub {
			b.difference(at[n.a], at[n.b])
		}
		.mul {
			b.product(at[n.a], at[n.b])
		}
		.div {
			b.quotient(at[n.a], at[n.b])
		}
		.pow {
			if b.is_constant(at[n.b]) {
				b.emit(.powc, at[n.a], -1, b.instrs[at[n.b]].c, 0)
			} else {
				b.emit(.pow, at[n.a], at[n.b], 0, 0)
			}
		}
		.call {
			b.call(functions[n.index], at[n.a])
		}
		.reading {
			// read_on's (si − post)/scale − pre, an operation skipped only
			// where it changes no bit: x − 0, x/1.
			mut x := at[n.a]
			if n.post != 0.0 {
				x = b.difference(x, b.constant(n.post))
			}
			if n.value != 1.0 {
				x = b.quotient(x, b.constant(n.value))
			}
			if n.pre != 0.0 {
				x = b.difference(x, b.constant(n.pre))
			}
			x
		}
		.scaled {
			b.product(at[n.a], b.constant(n.value))
		}
	}
}

fn (b &Builder) is_constant(slot int) bool {
	return b.instrs[slot].op == .constant
}

fn (mut b Builder) constant(c f64) int {
	return b.emit(.constant, -1, -1, c, 0)
}

// sum is a + b, a constant on either side fused as a + c.
fn (mut b Builder) sum(x int, y int) int {
	if b.is_constant(y) && !b.is_constant(x) {
		b.fused++
		return b.emit(.addc, x, -1, b.instrs[y].c, 0)
	}
	if b.is_constant(x) && !b.is_constant(y) {
		b.fused++
		return b.emit(.addc, y, -1, b.instrs[x].c, 0)
	}
	return b.emit(.add, x, y, 0, 0)
}

// difference is a − b: x − c as x + (−c), which IEEE arithmetic makes the
// same bits, and c − x as c − x.
fn (mut b Builder) difference(x int, y int) int {
	if b.is_constant(y) && !b.is_constant(x) {
		b.fused++
		return b.emit(.addc, x, -1, -b.instrs[y].c, 0)
	}
	if b.is_constant(x) && !b.is_constant(y) {
		b.fused++
		return b.emit(.csub, y, -1, b.instrs[x].c, 0)
	}
	return b.emit(.sub, x, y, 0, 0)
}

// product is a·b. A constant times exp(d/x) is c·exp(d/x), one operation;
// any other constant factor is fused as c·x.
fn (mut b Builder) product(x int, y int) int {
	for pair in [[x, y], [y, x]] {
		k, other := pair[0], pair[1]
		if b.is_constant(k) && !b.is_constant(other) {
			o := b.instrs[other]
			b.fused++
			if o.op == .aexp && o.c == 1.0 {
				return b.emit(.aexp, o.a, -1, b.instrs[k].c, o.d)
			}
			return b.emit(.mulc, other, -1, b.instrs[k].c, 0)
		}
	}
	return b.emit(.mul, x, y, 0, 0)
}

// quotient is a/b, a constant on either side fused.
fn (mut b Builder) quotient(x int, y int) int {
	if b.is_constant(y) && !b.is_constant(x) {
		b.fused++
		return b.emit(.divc, x, -1, b.instrs[y].c, 0)
	}
	if b.is_constant(x) && !b.is_constant(y) {
		b.fused++
		return b.emit(.cdiv, y, -1, b.instrs[x].c, 0)
	}
	return b.emit(.div, x, y, 0, 0)
}

// call is a function. exp(d/x) is 1·exp(d/x), which a constant factor joins
// later; exp(c + d/x) and exp(c − d/x) are exp(c + d/x), with −d for the
// second, since c − d/x is c + (−d)/x to the bit.
fn (mut b Builder) call(name string, x int) int {
	if name == 'exp' && !b.is_constant(x) {
		arg := b.instrs[x]
		if arg.op == .cdiv {
			b.fused++
			return b.emit(.aexp, arg.a, -1, 1.0, arg.c)
		}
		if arg.op == .addc && b.instrs[arg.a].op == .cdiv {
			inner := b.instrs[arg.a]
			b.fused++
			return b.emit(.expl, inner.a, -1, arg.c, inner.c)
		}
		if arg.op == .csub && b.instrs[arg.a].op == .cdiv {
			inner := b.instrs[arg.a]
			b.fused++
			return b.emit(.expl, inner.a, -1, arg.c, -inner.c)
		}
	}
	op := match name {
		'exp' { Op.exp }
		'ln' { Op.ln }
		'log10' { Op.log10 }
		'sqrt' { Op.sqrt }
		else { Op.abs }
	}
	return b.emit(op, x, -1, 0, 0)
}

// emit adds an operation, or finds it already there, or folds it into a
// constant if its operands are all constants.
fn (mut b Builder) emit(op Op, a int, bb int, c f64, d f64) int {
	mut x, mut y := a, bb
	if (op == .add || op == .mul) && x > y {
		x, y = y, x
	}
	if op != .input && op != .constant && b.instrs[x].op == .constant
		&& (y < 0 || b.instrs[y].op == .constant) {
		b.folded++
		return b.constant(fold_op(op, b.instrs[x].c, if y >= 0 { b.instrs[y].c } else { 0.0 }, c))
	}
	key := '${int(op)} ${x} ${y} ${math.f64_bits(c)} ${math.f64_bits(d)}'
	if slot := b.seen[key] {
		if op != .constant && op != .input {
			b.merged++
		}
		return slot
	}
	b.instrs << Instr{op, x, y, c, d}
	b.seen[key] = b.instrs.len - 1
	return b.instrs.len - 1
}

// fold_op is an operation of constants, with the arithmetic eval uses. A
// fused operation is never folded: it's made only with a varying operand.
fn fold_op(op Op, x f64, y f64, c f64) f64 {
	return match op {
		.add { x + y }
		.sub { x - y }
		.mul { x * y }
		.div { x / y }
		.pow { raise(x, y) }
		.neg { -x }
		.exp { math.exp(x) }
		.ln { math.log(x) }
		.log10 { math.log10(x) }
		.sqrt { math.sqrt(x) }
		.abs { math.abs(x) }
		.powc { raise(x, c) }
		else { math.nan() }
	}
}

// finish keeps what the outputs use, puts the constants first and the
// inputs next, and gives the tape.
fn (mut b Builder) finish(inputs int, outs []int) Tape {
	mut used := []bool{len: b.instrs.len}
	for o in outs {
		used[o] = true
	}
	for i := b.instrs.len - 1; i >= 0; i-- {
		if !used[i] {
			continue
		}
		n := b.instrs[i]
		if n.op != .input && n.op != .constant {
			used[n.a] = true
			if n.b >= 0 {
				used[n.b] = true
			}
		}
	}
	mut order := []int{}
	for i, n in b.instrs {
		if used[i] && n.op == .constant {
			order << i
		}
	}
	n_const := order.len
	for i, n in b.instrs {
		if used[i] && n.op == .input {
			order << i
		}
	}
	first := order.len
	for i, n in b.instrs {
		if used[i] && n.op != .constant && n.op != .input {
			order << i
		}
	}
	mut slot := []int{len: b.instrs.len, init: -1}
	for j, i in order {
		slot[i] = j
	}
	mut instrs := []Instr{cap: order.len}
	for i in order {
		n := b.instrs[i]
		if n.op == .input || n.op == .constant {
			instrs << n
		} else {
			instrs << Instr{n.op, slot[n.a], if n.b >= 0 { slot[n.b] } else { -1 }, n.c, n.d}
		}
	}
	return Tape{
		instrs:  instrs
		n_const: n_const
		first:   first
		inputs:  inputs
		outputs: outs.map(slot[it])
		merged:  b.merged
		folded:  b.folded
		fused:   b.fused
	}
}

// buffer makes an evaluation's workspace, with the constants in place.
pub fn (t &Tape) buffer() Buffer {
	mut buf := Buffer{
		v:   []f64{len: t.instrs.len}
		adj: []f64{len: t.instrs.len * t.outputs.len}
	}
	for i in 0 .. t.n_const {
		buf.v[i] = t.instrs[i].c
	}
	return buf
}

// eval evaluates the tape at the inputs x, in SI, into buf. output then
// reads each expression's value.
@[direct_array_access]
pub fn (t &Tape) eval(x []f64, mut buf Buffer) {
	for i := t.n_const; i < t.instrs.len; i++ {
		n := t.instrs[i]
		buf.v[i] = match n.op {
			.input { x[n.a] }
			.constant { n.c }
			.add { buf.v[n.a] + buf.v[n.b] }
			.sub { buf.v[n.a] - buf.v[n.b] }
			.mul { buf.v[n.a] * buf.v[n.b] }
			.div { buf.v[n.a] / buf.v[n.b] }
			.pow { raise(buf.v[n.a], buf.v[n.b]) }
			.neg { -buf.v[n.a] }
			.exp { math.exp(buf.v[n.a]) }
			.ln { math.log(buf.v[n.a]) }
			.log10 { math.log10(buf.v[n.a]) }
			.sqrt { math.sqrt(buf.v[n.a]) }
			.abs { math.abs(buf.v[n.a]) }
			.powc { raise(buf.v[n.a], n.c) }
			.addc { buf.v[n.a] + n.c }
			.csub { n.c - buf.v[n.a] }
			.mulc { n.c * buf.v[n.a] }
			.divc { buf.v[n.a] / n.c }
			.cdiv { n.c / buf.v[n.a] }
			.aexp { n.c * math.exp(n.d / buf.v[n.a]) }
			.expl { math.exp(n.c + n.d / buf.v[n.a]) }
		}
	}
}

// output is expression k's value after eval.
pub fn (t &Tape) output(buf &Buffer, k int) f64 {
	return buf.v[t.outputs[k]]
}

// gradient sweeps the tape back once, after eval, and writes into g every
// output's derivative with respect to every input: g[k*inputs + j] is
// d output k / d x[j]. An input no output uses has a derivative of 0.
@[direct_array_access]
pub fn (t &Tape) gradient(mut buf Buffer, mut g []f64) {
	m := t.outputs.len
	for i in 0 .. g.len {
		g[i] = 0.0
	}
	for i := t.n_const * m; i < buf.adj.len; i++ {
		buf.adj[i] = 0.0
	}
	for k, o in t.outputs {
		buf.adj[o * m + k] = 1.0
	}
	for i := t.instrs.len - 1; i >= t.first; i-- {
		n := t.instrs[i]
		w := i * m
		a := n.a * m
		v := buf.v[i]
		va := buf.v[n.a]
		match n.op {
			.add, .sub, .mul, .div, .pow {
				b := n.b * m
				vb := buf.v[n.b]
				// The derivatives with respect to a and b.
				mut da, mut db := 1.0, 1.0
				match n.op {
					.sub {
						db = -1.0
					}
					.mul {
						da, db = vb, va
					}
					.div {
						da = 1.0 / vb
						db = -v / vb
					}
					.pow {
						da = if va != 0.0 { vb * v / va } else { vb * raise(va, vb - 1.0) }
						db = if v != 0.0 { v * math.log(va) } else { 0.0 }
					}
					else {}
				}
				for k in 0 .. m {
					adj := buf.adj[w + k]
					buf.adj[a + k] += adj * da
					buf.adj[b + k] += adj * db
				}
			}
			else {
				d := match n.op {
					.neg { -1.0 }
					.exp { v }
					.ln { 1.0 / va }
					.log10 { 1.0 / (va * math.ln10) }
					.sqrt { 0.5 / v }
					.abs {
						if va > 0.0 {
							1.0
						} else if va < 0.0 {
							-1.0
						} else {
							0.0
						}
					}
					.powc {
						if va != 0.0 { n.c * v / va } else { n.c * raise(va, n.c - 1.0) }
					}
					.addc { 1.0 }
					.csub { -1.0 }
					.mulc { n.c }
					.divc { 1.0 / n.c }
					.cdiv { -v / va }
					.aexp, .expl { -v * n.d / (va * va) }
					else { 0.0 }
				}
				for k in 0 .. m {
					buf.adj[a + k] += buf.adj[w + k] * d
				}
			}
		}
	}
	for i in t.n_const .. t.first {
		j := t.instrs[i].a
		for k in 0 .. m {
			g[k * t.inputs + j] = buf.adj[i * m + k]
		}
	}
}
