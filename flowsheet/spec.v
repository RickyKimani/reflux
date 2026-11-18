module flowsheet

import expr
import math
import numeric
import quantity
import thermo
import toml.ast

// Spec is a design specification (M8's D5): vary an input, within its
// bounds, until an equation over the flowsheet holds. "Vary the cooler's
// outlet T between 200 and 260 K until the drum's vapour is 60 mol/s."
pub struct Spec {
pub:
	name  string
	input Input
	lo    f64 // the least the input may take, SI
	hi    f64 // the most
	until expr.Equation
	line  int // its table's line in the file
}

// SpecResult is a spec as solved: the value it varied its input to, its
// equation's two sides there, and how far apart they are as a fraction of
// their scale (spec_scale), within the solver's tolerance.
pub struct SpecResult {
pub:
	name     string
	key      string
	value    f64
	lhs      f64
	rhs      f64
	residual f64
}

// read_specs reads [spec.NAME] (D5, D16): vary, a key the file gives
// (D4); between, its least and most values, with their units; and until,
// an equation over the flowsheet's variables (D3).
fn read_specs(mut specs Section, l Loaded) ![]Spec {
	mut out := []Spec{}
	for name in specs.table.keys() {
		mut t := specs.child(name)!
		key := t.text('vary')!
		input := l.input(key) or { return at_key(err, t.line('vary'), t.full('vary')) }
		v := t.require('between')!
		if v !is []ast.Value || (v as []ast.Value).len != 2 {
			return t.fail('between', 'it needs two values, the least and the most ${key} may take, as ["200 K", "260 K"]',
				[])
		}
		pair := v as []ast.Value
		lo := bound(t, input, pair[0])!
		hi := bound(t, input, pair[1])!
		if !(lo < hi) {
			return t.fail('between', 'its least value must be below its most', [])
		}
		if !(input.value >= lo && input.value <= hi) {
			return t.fail('between', "the file's ${key}, ${input.text}, is outside it, and a spec starts from the file's value",
				[])
		}
		text := t.text('until')!
		q := expr.equation(l.variables.context, &no_definitions, text) or {
			if err is expr.ExprError {
				e := err as expr.ExprError
				return file_error(t.line('until'), t.full('until'), 'column ${e.column}: ${e.reason}',
					e.suggestions)
			}
			return t.fail('until', err.msg(), [])
		}
		if q.relation != .eq {
			return t.fail('until', "a spec's until is an equation, with =: <= and >= are an optimisation's constraints",
				[])
		}
		t.finish()!
		for other in out {
			if other.input.key == key {
				return t.fail('vary', 'spec ${other.name} varies ${key} already', [])
			}
		}
		out << Spec{
			name:  name
			input: input
			lo:    lo
			hi:    hi
			until: q
			line:  specs.line_of_root('spec.${name}')
		}
	}
	return out
}

// bound is one of a spec's bounds, read as its input's maker reads it.
fn bound(t Section, input Input, v ast.Value) !f64 {
	if v is ast.Number {
		if input.dim.is_dimensionless() {
			return v.f64()
		}
		return t.fail('between', 'a bound of ${input.key} needs its unit, in quotes, as "${v.text} ${si_text(input.dim)}"',
			[])
	}
	if v is ast.Quoted {
		x := input.parse(v.text) or {
			return file_error(v.pos.line_nr, t.full('between'), err.msg(), []string{})
		}
		return x
	}
	return t.fail('between', 'it needs two values, numbers or quantities in quotes', [])
}

// parse reads a value of this input from text, as its maker would: in one
// of M4's units for its dimension, a percentage for a fraction, or any
// unit of its dimension.
pub fn (i Input) parse(text string) !f64 {
	if named := i.named {
		return quantity.parse(text, named) or {
			if err is quantity.QuantityError {
				return error(err.reason)
			}
			return err
		}
	}
	if i.dim.is_dimensionless() {
		return quantity.parse(text, .fraction) or {
			if err is quantity.QuantityError {
				return error(err.reason)
			}
			return err
		}
	}
	e := expr.parse(expr.Context{}, &no_definitions, text) or { return error(reason_of(err)) }
	e.must_be(i.dim, i.key) or { return error(reason_of(err)) }
	return e.eval(&no_definitions, []f64{})
}

// at_key gives an error a key and line in the file.
fn at_key(err IError, line int, key string) IError {
	if err is FlowsheetError {
		return FlowsheetError{
			kind:        .bad_file
			reason:      err.reason
			line:        line
			key:         key
			suggestions: err.suggestions
		}
	}
	return file_error(line, key, err.msg(), [])
}

// links are the specs' dependences, for planning (D5): from each unit that
// makes something a spec's equation reads to the unit it varies, or that
// the feed it varies enters.
fn (l Loaded) links() ![]Link {
	mut out := []Link{}
	for k, spec in l.specs {
		mut from := []int{}
		for e in [spec.until.lhs, spec.until.rhs] {
			for n in e.nodes {
				if n.kind != .var && n.kind != .member {
					continue
				}
				slot := l.variables.slots[n.index]
				u := if slot.measure == .energy { slot.unit } else { l.sheet.maker_of(slot.stream) }
				if u >= 0 && u !in from {
					from << u
				}
			}
		}
		if from.len == 0 {
			return file_error(spec.line, 'spec.${spec.name}', "it reads only feeds, which nothing it varies can change: name what a unit makes, as a product's flow",
				[])
		}
		mut to := []int{}
		if spec.input.table == 'unit' {
			to << spec.input.index
		} else {
			for i, p in l.sheet.units {
				if spec.input.name in p.inlets {
					to << i
				}
			}
		}
		for a in from {
			for b in to {
				out << Link{
					from: a
					to:   b
					spec: k
					name: spec.name
					line: spec.line
				}
			}
		}
	}
	return out
}

// maker_of is the unit a stream leaves, or -1 for a feed.
fn (f Flowsheet) maker_of(stream string) int {
	for i, p in f.units {
		if stream in p.outlets {
			return i
		}
	}
	return -1
}

// with_specs solves a block that has specs (D5). One spec and no tears:
// Brent's method within its bounds, each evaluation a pass. Otherwise
// Broyden on the tears and the specs' values together, each evaluation a
// pass; and if that fails where there are tears, once more nested, the
// specs' values around the tears converged.
fn (mut s Solver) with_specs(b Block) ! {
	s.block = b
	s.rows = [][]f64{}
	s.spec_rows = [][]f64{}
	s.inner_started = false
	s.flow, s.t_in = s.entering(b)
	s.warm_values()
	if b.tears.len == 0 {
		if b.specs.len == 1 {
			s.spec_brent(b.specs[0]) or { return s.failed(err) }
		} else {
			s.spec_broyden() or { return s.failed(err) }
		}
		s.keep_history()
		return
	}
	s.together() or {
		if err is FlowsheetError && err.kind !in [.no_convergence, .spec_not_met] {
			return s.failed(err)
		}
		l := s.file or { return err }
		for k in b.specs {
			s.values[k] = l.specs[k].input.value
		}
		s.warm_values()
		s.inner_started = false
		if b.specs.len == 1 {
			s.spec_brent(b.specs[0]) or { return s.failed(err) }
		} else {
			s.spec_broyden() or { return s.failed(err) }
		}
	}
	s.keep_history()
	s.add_gap()
}

// entering is the flow into a block from outside it, and its hottest T,
// for the scale of a spec's residual.
fn (s Solver) entering(b Block) (f64, f64) {
	mut flow, mut t := 0.0, 0.0
	mut made := []string{}
	for i in b.units {
		made << s.sheet.units[i].outlets
	}
	for i in b.units {
		for name in s.sheet.units[i].inlets {
			if name in made {
				continue
			}
			if st := s.streams[name] {
				flow += st.total()
				t = math.max(t, st.t)
			}
		}
	}
	return flow, t
}

// apply remakes the units and feeds the block's specs vary, at their
// values now (D4).
fn (mut s Solver) apply() ! {
	l := s.file or { return }
	mut inputs := []Input{}
	mut values := []f64{}
	for k in s.block.specs {
		inputs << l.specs[k].input
		values << s.values[k]
	}
	made := l.with(inputs, values)!
	for input in inputs {
		if input.table == 'unit' {
			s.sheet.units[input.index] = made.sheet.units[input.index]
		} else {
			s.sheet.feeds[input.index] = made.sheet.feeds[input.index]
			s.streams[input.name] = made.sheet.feeds[input.index].stream
		}
	}
}

// inner calculates the block at its specs' values now: a pass, or, nested,
// its tears converged from where they were left.
fn (mut s Solver) inner() ! {
	if s.block.tears.len == 0 {
		s.pass(s.block)!
		return
	}
	if !s.inner_started {
		s.begin(s.block)!
		s.inner_started = true
	}
	if s.params.method == .broyden {
		s.by_broyden()!
	} else {
		s.by_substitution()!
	}
}

// spec_scale is what a spec's residual is a fraction of: the larger of its
// two sides, or, where both are near zero, the block's own scale for that
// dimension, a duty's F·RT, a flow's F, a pure number's 1.
fn (s Solver) spec_scale(d quantity.Dim, lhs f64, rhs f64, flow f64, t f64) f64 {
	mut floor := 0.0
	if d.eq(quantity.Dimension.power.dim()) {
		floor = flow * thermo.r_gas * t
	} else if d.eq(quantity.Dimension.molar_flow.dim()) {
		floor = flow
	} else if d.eq(quantity.Dimension.molar_enthalpy.dim()) {
		floor = thermo.r_gas * t
	} else if d.is_dimensionless() {
		floor = 1.0
	}
	scale := math.max(floor, math.max(math.abs(lhs), math.abs(rhs)))
	return if scale > 0.0 { scale } else { 1.0 }
}

// residual is a spec's two sides now, and their difference as a fraction
// of their scale.
fn (mut s Solver) residual(k int) !(f64, f64, f64) {
	l := s.file or { return error('no specs') }
	spec := l.specs[k]
	x, why := l.variables.fill(s.streams, s.runs)
	lhs := spec.until.lhs.value(&no_definitions, x, why) or { return s.no_value(spec, err) }
	rhs := spec.until.rhs.value(&no_definitions, x, why) or { return s.no_value(spec, err) }
	s.floors[k] = [s.flow, s.t_in]
	scale := s.spec_scale(spec.until.lhs.dim, lhs, rhs, s.flow, s.t_in)
	return lhs, rhs, (lhs - rhs) / scale
}

fn (s Solver) no_value(spec Spec, err IError) IError {
	return FlowsheetError{
		kind:   .spec_not_met
		reason: 'at ${spec.input.key} = ${shown(s.values[s.index_of(spec)], spec.input.dim)}, its equation has no value: ${reason_of(err)}'
		line:   spec.line
		key:    'spec.${spec.name}'
	}
}

fn (s Solver) index_of(spec Spec) int {
	l := s.file or { return 0 }
	for k, other in l.specs {
		if other.name == spec.name {
			return k
		}
	}
	return 0
}

// try is the block calculated with spec k at v: its residual, after the
// pass or the converged tears.
fn (mut s Solver) try(k int, v f64) !f64 {
	s.values[k] = v
	s.apply()!
	s.inner()!
	lhs, rhs, f := s.residual(k)!
	s.samples << Sample{v, lhs, rhs, f}
	s.spec_rows << [f]
	return f
}

// Sample is a spec's equation at one value of its input, for saying what
// its bounds reach.
struct Sample {
	value    f64
	lhs      f64
	rhs      f64
	residual f64
}

// SpecRoot is a spec's residual as a function of its value, for brent.
struct SpecRoot {
mut:
	s &Solver
	k int
}

// at is the residual at v. A value sampled already isn't calculated again:
// brent starts from the bracket's two ends, which the sampling has.
fn (r SpecRoot) at(v f64) !f64 {
	for x in r.s.samples {
		if x.value == v {
			return x.residual
		}
	}
	return r.s.try(r.k, v)
}

// spec_brent meets one spec by Brent's method within its bounds (D5): the
// residual at each bound, and if they have one sign, at eight points
// between, for a bracket. A target that no point brackets is refused, with
// the range the points reached.
fn (mut s Solver) spec_brent(k int) ! {
	l := s.file or { return }
	spec := l.specs[k]
	s.samples = []Sample{}
	root := SpecRoot{
		s: &s
		k: k
	}
	// From a neighbour's value, a bracket is looked for near it first.
	if s.start != none {
		if s.near(root, spec)! {
			return
		}
	}
	mut a, mut fa := spec.lo, root.at(spec.lo)!
	mut b, mut fb := spec.hi, root.at(spec.hi)!
	if bracket := s.bracket() {
		a, b = bracket[0], bracket[1]
		fa, fb = root.at(a)!, root.at(b)!
	} else if fa != 0.0 && fb != 0.0 && (fa > 0.0) == (fb > 0.0) {
		mut found := false
		mut last, mut f_last := spec.lo, fa
		for i in 1 .. 10 {
			v := if i < 9 { spec.lo + (spec.hi - spec.lo) * f64(i) / 9.0 } else { spec.hi }
			f := if i < 9 { root.at(v)! } else { fb }
			if f == 0.0 || (f > 0.0) != (f_last > 0.0) {
				a, fa, b, fb = last, f_last, v, f
				found = true
				break
			}
			last, f_last = v, f
		}
		if !found {
			return s.out_of_reach(spec)
		}
	}
	r := numeric.brent(root, a, b, ftol: s.params.tol, max_iter: s.params.max_iterations) or {
		return FlowsheetError{
			kind:   .spec_not_met
			reason: 'Brent within [${shown(a, spec.input.dim)}, ${shown(b, spec.input.dim)}]: ${err.msg()}'
			line:   spec.line
			key:    'spec.${spec.name}'
		}
	}
	// The state must be the root's, which a cached value isn't.
	if s.values[k] != r.x {
		s.try(k, r.x)!
	}
	_ = fa
	_ = fb
}

// warm_values starts the block's specs from the last point's values, or
// the secant through the last two, where the solve starts from a neighbour.
fn (mut s Solver) warm_values() {
	l := s.file or { return }
	for k in s.block.specs {
		spec := l.specs[k]
		if v := s.warm_value(k, spec.lo, spec.hi) {
			s.values[k] = v
		}
	}
}

// near meets a spec from its value now, a neighbour's: met there already,
// or bracketed by secant steps from it, each at most four times the last,
// the first as large as the last two points' change, then Brent. It says
// whether it met the spec; if it found no bracket within its bounds, the
// cold search follows, with what it sampled.
fn (mut s Solver) near(root SpecRoot, spec Spec) !bool {
	k := root.k
	v0 := s.values[k]
	f0 := root.at(v0)!
	if math.abs(f0) <= s.params.tol {
		return true
	}
	mut step := 1e-3 * (spec.hi - spec.lo)
	if start := s.start {
		if before := start.before {
			if a := before.specs[k] {
				if b := start.last.specs[k] {
					if a.value != b.value {
						step = math.abs(b.value - a.value)
					}
				}
			}
		}
	}
	mut x0, mut g0 := v0, f0
	mut x1 := math.min(spec.hi, v0 + step)
	if x1 == v0 {
		x1 = math.max(spec.lo, v0 - step)
	}
	for _ in 0 .. 12 {
		root.at(x1)!
		if bracket := s.bracket() {
			s.brent_in(root, spec, bracket[0], bracket[1])!
			return true
		}
		g1 := s.samples.last().residual
		mut next := x1 - g1 * (x1 - x0) / (g1 - g0)
		reach := 4.0 * math.abs(x1 - x0)
		if !math.is_finite(next) {
			next = x1 + 2.0 * (x1 - x0)
		}
		next = math.max(x1 - reach, math.min(x1 + reach, next))
		next = math.max(spec.lo, math.min(spec.hi, next))
		if next == x1 {
			return false
		}
		x0, g0, x1 = x1, g1, next
	}
	return false
}

// bracket is the first pair of neighbouring values sampled so far whose
// residuals have opposite signs, or one that is zero.
fn (s Solver) bracket() ?[]f64 {
	mut samples := s.samples.clone()
	samples.sort(a.value < b.value)
	for i in 1 .. samples.len {
		a, b := samples[i - 1], samples[i]
		if a.residual == 0.0 || b.residual == 0.0 || (a.residual > 0.0) != (b.residual > 0.0) {
			return [a.value, b.value]
		}
	}
	return none
}

// brent_in meets a spec by Brent's method in a bracket, and leaves the
// block at its root.
fn (mut s Solver) brent_in(root SpecRoot, spec Spec, a f64, b f64) ! {
	r := numeric.brent(root, a, b, ftol: s.params.tol, max_iter: s.params.max_iterations) or {
		return FlowsheetError{
			kind:   .spec_not_met
			reason: 'Brent within [${shown(a, spec.input.dim)}, ${shown(b, spec.input.dim)}]: ${err.msg()}'
			line:   spec.line
			key:    'spec.${spec.name}'
		}
	}
	// The state must be the root's, which a cached value isn't.
	if s.values[root.k] != r.x {
		s.try(root.k, r.x)!
	}
}

// out_of_reach refuses a target its bounds don't bracket, with what they
// reach.
fn (s Solver) out_of_reach(spec Spec) IError {
	mut samples := s.samples.clone()
	samples.sort(a.value < b.value)
	lo, hi := shown(spec.lo, spec.input.dim), shown(spec.hi, spec.input.dim)
	rhs_fixed := samples.all(it.rhs == samples[0].rhs)
	d := spec.until.lhs.dim
	mut least, mut most := math.inf(1), math.inf(-1)
	for x in samples {
		y := if rhs_fixed { x.lhs } else { x.lhs - x.rhs }
		least = math.min(least, y)
		most = math.max(most, y)
	}
	what := if rhs_fixed {
		'${spec.until.lhs.plain()} goes from ${shown(least, d)} to ${shown(most, d)}, never ${shown(samples[0].rhs, d)}'
	} else {
		'its left side less its right goes from ${shown(least, d)} to ${shown(most, d)}, never 0'
	}
	return FlowsheetError{
		kind:   .spec_not_met
		reason: 'its target is out of reach: with ${spec.input.key} from ${lo} to ${hi}, ${what} (at ${samples.len} values)'
		line:   spec.line
		key:    'spec.${spec.name}'
	}
}

// shown is a value for a message, to six figures, with its SI unit.
fn shown(x f64, d quantity.Dim) string {
	mut v := x
	if x != 0.0 && math.is_finite(x) {
		m := math.pow(10.0, 5.0 - math.floor(math.log10(math.abs(x))))
		v = math.round(x * m) / m
	}
	unit := if d.is_dimensionless() { '' } else { ' ${si_text(d)}' }
	number := if v == math.trunc(v) && math.abs(v) < 1e15 { i64(v).str() } else { v.str() }
	return number + unit
}

// Specs is a block's specs' values as Broyden's unknowns, each evaluation
// the block calculated there: a pass, or its tears converged.
struct Specs {
mut:
	s    &Solver
	base []f64
	m    []f64 // the values are base + m·u, m = −J⁻¹ so that ∂F/∂u ≈ −I at first
}

fn (mut o Specs) eval(u []f64, mut f []f64) ! {
	o.s.set_values(u, o.base, o.m)
	o.s.apply()!
	o.s.inner()!
	for j, k in o.s.block.specs {
		_, _, r := o.s.residual(k)!
		f[j] = r
	}
	o.s.spec_rows << f.clone()
}

// set_values sets the block's specs' values for Broyden's unknowns u:
// base + m·u, each held within its bounds. It says whether any was held.
fn (mut s Solver) set_values(u []f64, base []f64, m []f64) bool {
	l := s.file or { return false }
	n := base.len
	mut held := false
	for j, k in s.block.specs {
		mut v := base[j]
		for i in 0 .. n {
			v += m[j * n + i] * u[i]
		}
		spec := l.specs[k]
		s.values[k] = math.max(spec.lo, math.min(spec.hi, v))
		held = held || s.values[k] != v
	}
	return held
}

// spec_jacobian is how the block's specs' residuals move with their values,
// by a step of 1e-4 of each one's bounds' width from base, and from it m =
// −J⁻¹, which makes Broyden's first step Newton's on the specs. Each step is
// from the tears' start where they're given, the pass's Jacobian, or else
// from the block as the last step left it. Specs whose targets don't move
// independently with what they vary can't be met from there.
fn (mut s Solver) spec_jacobian(base []f64, f0 []f64, start []Stream) ![]f64 {
	l := s.file or { return []f64{} }
	specs := s.block.specs
	n := specs.len
	mut jac := []f64{len: n * n}
	for j, k in specs {
		spec := l.specs[k]
		mut h := 1e-4 * (spec.hi - spec.lo)
		if base[j] + h > spec.hi {
			h = -h
		}
		for t, name in s.block.tears {
			if t < start.len {
				s.streams[name] = start[t]
			}
		}
		s.values[k] = base[j] + h
		s.apply()!
		if start.len > 0 {
			s.pass(s.block)!
		} else {
			s.inner()!
		}
		for i, other in specs {
			_, _, r := s.residual(other)!
			jac[i * n + j] = (r - f0[i]) / h
		}
		s.values[k] = base[j]
	}
	first := l.specs[specs[0]]
	stuck := FlowsheetError{
		kind:   .spec_not_met
		reason: if n == 1 {
			'${first.until.plain()} does not move with ${first.input.key} at ${shown(base[0], first.input.dim)}'
		} else {
			"the targets of specs ${specs.map(l.specs[it].name).join(', ')} don't move independently with what they vary, so they can't be met together from here"
		}
		line:   first.line
		key:    'spec.${first.name}'
	}
	lu := numeric.lu_factor(jac, n) or { return stuck }
	mut m := []f64{len: n * n}
	for j in 0 .. n {
		mut e := []f64{len: n}
		e[j] = 1.0
		col := lu.solve(e) or { return stuck }
		for i in 0 .. n {
			m[i * n + j] = -col[i]
			if !math.is_finite(col[i]) {
				return stuck
			}
		}
	}
	return m
}

// spec_broyden meets a block's specs by Broyden's method on their values,
// each evaluation the block at those values (D5).
fn (mut s Solver) spec_broyden() ! {
	l := s.file or { return }
	specs := s.block.specs
	base := specs.map(s.values[it])
	s.apply()!
	s.inner()!
	mut f0 := []f64{}
	for k in specs {
		_, _, r := s.residual(k)!
		f0 << r
	}
	m := s.spec_jacobian(base, f0, []Stream{})!
	s.apply()!
	s.inner()!
	mut o := Specs{
		s:    &s
		base: base
		m:    m
	}
	numeric.broyden(mut o, []f64{len: specs.len},
		ftol:     s.params.tol
		max_iter: s.params.max_iterations
	) or {
		spec := l.specs[specs[0]]
		return FlowsheetError{
			kind:   .spec_not_met
			reason: 'specs ${specs.map(l.specs[it].name).join(', ')}: ${err.msg()}'
			line:   spec.line
			key:    'spec.${spec.name}'
		}
	}
}

// together meets a block's specs with its tears, by Broyden's method on the
// tears' variables and the specs' values at once (D5): each evaluation is
// one pass. The specs' first Jacobian is a pass's, from the tears' start.
fn (mut s Solver) together() ! {
	b := s.block
	s.begin(b)!
	start := b.tears.map(s.streams[it])
	s.apply()!
	s.pass(b)!
	mut f0 := []f64{}
	for k in b.specs {
		_, _, r := s.residual(k)!
		f0 << r
	}
	base := b.specs.map(s.values[it])
	m := s.spec_jacobian(base, f0, start)!
	for t, name in b.tears {
		s.streams[name] = start[t]
	}
	s.apply()!
	s.pass(b)!
	x := s.tear_values()
	n := s.nc + 2
	s.scale = []f64{len: x.len}
	for t in 0 .. b.tears.len {
		for i in 0 .. s.nc {
			s.scale[t * n + i] = s.flow
		}
		first := s.streams[b.tears[t]]
		s.scale[t * n + s.nc] = first.p
		s.scale[t * n + s.nc + 1] = thermo.r_gas * first.t
	}
	mut y := []f64{len: x.len + b.specs.len}
	for i, v in x {
		y[i] = v / s.scale[i]
	}
	s.base = base
	s.m = m
	s.fresh = true
	s.joint = true
	s.held = 0
	defer {
		s.joint = false
	}
	numeric.broyden(mut s, y,
		ftol:     s.params.tol
		max_iter: s.params.max_iterations - 2 - b.specs.len
	) or {
		if err is numeric.BroydenError {
			return s.not_converged(err.reason)
		}
		return err
	}
}

// spec_results checks each spec in a solution from scratch, as every
// balance is checked (D8): its two sides from the solution's streams and
// units, apart by no more than the tolerance of their scale, and its value
// within its bounds.
fn (s Solver) spec_results(sol Solution) ![]SpecResult {
	l := s.file or { return []SpecResult{} }
	x, why := l.variables.values(sol)
	mut out := []SpecResult{}
	for k, spec in l.specs {
		lhs := spec.until.lhs.value(&no_definitions, x, why) or { return s.no_value(spec, err) }
		rhs := spec.until.rhs.value(&no_definitions, x, why) or { return s.no_value(spec, err) }
		floor := s.floors[k] or { [0.0, 0.0] }
		scale := s.spec_scale(spec.until.lhs.dim, lhs, rhs, floor[0], floor[1])
		r := (lhs - rhs) / scale
		v := s.values[k]
		if !(math.abs(r) <= s.params.tol) || !(v >= spec.lo && v <= spec.hi) {
			return FlowsheetError{
				kind:    .spec_not_met
				reason:  'at ${spec.input.key} = ${shown(v, spec.input.dim)}, ${spec.until.plain()} is off by ${r} of its scale'
				line:    spec.line
				key:     'spec.${spec.name}'
				history: s.history
			}
		}
		out << SpecResult{
			name:     spec.name
			key:      spec.input.key
			value:    v
			lhs:      lhs
			rhs:      rhs
			residual: r
		}
	}
	return out
}
