// The solver (exit criteria 3 and 4): loops with exact answers worked by hand,
// the reference flowsheet by each method, and the same answer whichever stream
// is torn. Every solved flowsheet is checked from scratch by check: each unit's
// component, element and energy balances to 1e-10 on the inlets it used, with H
// recomputed from T, P and each phase's composition, and the whole to the tear
// tolerance. Failures are typed, and carry the history.
import components
import flash
import flowsheet { Flowsheet, FlowsheetError, Outcome, Solution, Stream, Unit }
import math
import thermo

const db = components.load_default() or { panic(err) }
const names = ['methane', 'ethane', 'propane', 'n-butane', 'n-pentane']
const feed_z = [0.70, 0.10, 0.08, 0.07, 0.05]
// Atoms of carbon and hydrogen in each, for the element balances.
const carbon = [1.0, 2.0, 3.0, 4.0, 5.0]
const hydrogen = [4.0, 6.0, 8.0, 10.0, 12.0]

fn fluid() !thermo.Fluid {
	return db.fluid(.pr, names)!
}

fn flows(total f64) []f64 {
	return feed_z.map(it * total)
}

// molar_h is a stream's molar H, each phase's volume root solved again.
fn molar_h(fl thermo.Fluid, s Stream) !f64 {
	mut h := 0.0
	for phase in s.result.phases {
		if phase.fraction == 0.0 {
			continue
		}
		root := if phase.label == .vapor { thermo.Root.vapor } else { thermo.Root.liquid }
		h += phase.fraction * fl.props(thermo.solve_tp(fl.model, s.t, s.p, phase.state.n, root)!).h
	}
	return h
}

// Tally adds up flows, elements and enthalpy from scratch.
struct Tally {
mut:
	comp   []f64
	carbon f64
	hydro  f64
	energy f64
	flow   f64
}

fn (mut t Tally) add(fl thermo.Fluid, s Stream, sign f64) ! {
	for i, n in s.flows {
		t.comp[i] += sign * n
		t.carbon += sign * n * carbon[i]
		t.hydro += sign * n * hydrogen[i]
	}
	t.energy += sign * s.total() * molar_h(fl, s)!
	t.flow += s.total()
}

// check checks a solution from scratch (criterion 3), and returns it.
fn check(fl thermo.Fluid, sheet Flowsheet, sol Solution, tol f64, at string) !Solution {
	for r in sol.units {
		mut t := Tally{
			comp: []f64{len: names.len}
		}
		mut hot := 0.0
		for s in r.inlets {
			t.add(fl, s, 1.0)!
			hot = math.max(hot, s.t)
		}
		for s in r.outlets {
			t.add(fl, s, -1.0)!
		}
		for e in r.energy {
			t.energy += e.rate
		}
		if t.flow == 0.0 {
			continue
		}
		for i, c in t.comp {
			assert math.abs(c) <= 1e-10 * t.flow, '${at}: ${r.name}: component ${i} is off by ${c}'
		}
		assert math.abs(t.carbon) <= 1e-9 * t.flow, '${at}: ${r.name}: carbon'
		assert math.abs(t.hydro) <= 1e-9 * t.flow, '${at}: ${r.name}: hydrogen'
		assert math.abs(t.energy) <= 1e-10 * t.flow * thermo.r_gas * hot, '${at}: ${r.name}: energy is off by ${t.energy}'
	}
	// The whole, feeds against products, to the tear tolerance: a tear's H
	// carries |H| with any wrong flow, and H is about 40 RT here.
	mut all := Tally{
		comp: []f64{len: names.len}
	}
	mut feed := 0.0
	for x in sheet.feeds {
		all.add(fl, x.stream, 1.0)!
		feed += x.stream.total()
	}
	for name in sol.products {
		all.add(fl, sol.streams[name], -1.0)!
	}
	mut loop := feed
	for b in sol.plan.blocks {
		for name in b.tears {
			loop += sol.streams[name].total()
		}
	}
	for r in sol.units {
		for e in r.energy {
			all.energy += e.rate
		}
	}
	for i, c in all.comp {
		assert math.abs(c) <= 10.0 * tol * loop, '${at}: component ${i} overall is off by ${c}'
	}
	assert math.abs(all.carbon) <= 100.0 * tol * loop, '${at}: carbon overall'
	assert math.abs(all.hydro) <= 100.0 * tol * loop, '${at}: hydrogen overall'
	assert math.abs(all.energy) <= 100.0 * tol * loop * 2e3 * 40.0, '${at}: energy overall is off by ${all.energy}'
	return sol
}

// solved solves a flowsheet and checks it.
fn solved(fl thermo.Fluid, f Flowsheet, params flowsheet.SolveParams) !Solution {
	return check(fl, f, f.solve(fl, params)!, params.tol, '${params.method}')!
}

// passes is the iterations of a solution's one block with tears.
fn passes(s Solution) int {
	return s.history[0].residuals.len
}

// recycle is a feed into a mixer, an adiabatic heater and a splitter that sends
// a fraction back.
fn recycle(fl thermo.Fluid, back f64) !Flowsheet {
	mut f := Flowsheet{}
	f.add_feed('feed', flowsheet.feed_tp(fl, flows(10.0), 300.0, 5e5)!)!
	f.add_unit(flowsheet.new_mixer('mixer', 2)!, {
		'in1': 'feed'
		'in2': 'back'
		'out': 'mixed'
	})!
	f.add_unit(flowsheet.new_heater('heater', mode: .duty, value: 0.0)!, {
		'in':  'mixed'
		'out': 'warm'
	})!
	f.add_unit(flowsheet.new_splitter('split', [back, 1.0 - back])!, {
		'in':   'warm'
		'out1': 'back'
		'out2': 'product'
	})!
	return f
}

fn test_a_recycle_has_its_exact_answer() {
	// With a fraction φ sent back, the loop carries F/(1 − φ), and the product
	// is the feed: here, 25 mol/s round the loop for 10 in, and 15 back. The
	// state doesn't change, since nothing heats it. A single gas at one P.
	fl := fluid()!
	f := recycle(fl, 0.6)!
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		s := solved(fl, f, method: m)!
		// The loop's flow goes to 1e-9 of itself, in what moved in the last pass,
		// and direct substitution's error is 1.5 times that at a gain of 0.6.
		assert math.abs(s.streams['mixed'].total() - 25.0) <= 1e-7, '${m}'
		assert math.abs(s.streams['back'].total() - 15.0) <= 1e-7, '${m}'
		assert math.abs(s.streams['product'].total() - 10.0) <= 1e-7, '${m}'
		assert math.abs(s.streams['product'].t - 300.0) <= 1e-7, '${m}'
		assert s.products == ['product']
		assert s.plan.blocks.len == 1 && s.plan.blocks[0].tears == ['back']
	}
}

fn test_each_method_converges_as_it_should() {
	// Direct substitution on g(R) = φ(F + R) closes the error by φ a pass. So
	// its residuals fall by 0.6 each, and it takes ln(tol)/ln φ passes: 41 here,
	// from an empty recycle.
	fl := fluid()!
	f := recycle(fl, 0.6)!
	d := solved(fl, f, method: .direct)!
	rows := d.history[0].residuals
	assert d.history[0].tears == ['back']
	for k in 20 .. rows.len {
		ratio := rows[k][0] / rows[k - 1][0]
		assert math.abs(ratio - 0.6) <= 0.01, 'pass ${k}: ${ratio}'
	}
	assert rows.len >= 38 && rows.len <= 44, '${rows.len}'
	assert rows.last()[0] <= 1e-9
	assert rows[rows.len - 2][0] > 1e-9

	// The loop is linear, so Wegstein's first step lands on the answer, and the
	// next pass shows it.
	w := solved(fl, f, method: .wegstein)!
	assert passes(w) <= 4, '${passes(w)}'
	// Broyden's first step is direct, and then it's Newton's on a linear loop.
	b := solved(fl, f, method: .broyden)!
	assert passes(b) <= 9, '${passes(b)}'
}

fn test_q_is_bounded() {
	// A gain near 1 is where Wegstein's unbounded step would go furthest. At a
	// fraction of 0.95 the answer is 200 mol/s round the loop, and its q is −19,
	// which would land there at once. Held to −5, each acceleration leaves
	// −5 + 6·0.95 = 0.7 of the error and the direct step 0.95 of that, so the
	// residual falls by 0.665 every two passes, from where the bound takes hold.
	fl := fluid()!
	f := recycle(fl, 0.95)!
	s := solved(fl, f, method: .wegstein, max_iterations: 200)!
	assert math.abs(s.streams['mixed'].total() - 200.0) <= 1e-5
	rows := s.history[0].residuals
	k := rows.len - 3
	ratio := rows[k][0] / rows[k - 2][0]
	assert math.abs(ratio - 0.665) <= 0.01, '${ratio}'
	assert rows.len <= 88, '${rows.len}'
}

// reference is M4's reference flowsheet: feed, mixer, cooler, flash drum,
// splitter, a pump back to the mixer, and a valve to a second drum.
fn reference(fl thermo.Fluid) !Flowsheet {
	mut f := Flowsheet{}
	f.add_feed('feed', flowsheet.feed_tp(fl, flows(100.0), 300.0, 30e5)!)!
	f.add_unit(flowsheet.new_mixer('mixer', 2)!, {
		'in1': 'feed'
		'in2': 'recycle'
		'out': 'mixed'
	})!
	f.add_unit(flowsheet.new_heater('cooler', mode: .outlet_t, value: 230.0)!, {
		'in':  'mixed'
		'out': 'cold'
	})!
	f.add_unit(flowsheet.new_flash_drum('drum', mode: .duty, value: 0.0)!, {
		'in1':    'cold'
		'vapour': 'gas'
		'liquid': 'liquid'
	})!
	f.add_unit(flowsheet.new_splitter('split', [0.5, 0.5])!, {
		'in':   'liquid'
		'out1': 'back'
		'out2': 'on'
	})!
	f.add_unit(flowsheet.new_pump('pump', mode: .outlet_p, value: 31e5, efficiency: 0.75)!, {
		'in':  'back'
		'out': 'recycle'
	})!
	f.add_unit(flowsheet.new_valve('valve', mode: .outlet_p, value: 5e5)!, {
		'in':  'on'
		'out': 'let down'
	})!
	f.add_unit(flowsheet.new_flash_drum('second', mode: .duty, value: 0.0)!, {
		'in1':    'let down'
		'vapour': 'flash gas'
		'liquid': 'product'
	})!
	return f
}

// same asserts that two solutions agree on every stream.
fn same(a Solution, b Solution, tol f64, at string) {
	for name, x in a.streams {
		y := b.streams[name]
		assert math.abs(x.t - y.t) <= tol * x.t, '${at}: ${name}: T ${x.t} vs ${y.t}'
		assert math.abs(x.p - y.p) <= tol * x.p, '${at}: ${name}: P'
		assert math.abs(x.beta - y.beta) <= tol * 10.0, '${at}: ${name}: β ${x.beta} vs ${y.beta}'
		for i, n in x.flows {
			assert math.abs(n - y.flows[i]) <= tol * math.max(x.total(), 1.0), '${at}: ${name}: flow ${i}: ${n} vs ${y.flows[i]}'
		}
	}
}

fn test_every_method_gets_the_reference_flowsheet_to_one_answer() {
	// Criterion 4. The iterations are numbers that can't grow.
	fl := fluid()!
	f := reference(fl)!
	d := solved(fl, f, method: .direct)!
	w := solved(fl, f, method: .wegstein)!
	b := solved(fl, f, method: .broyden)!
	println('reference flowsheet: direct ${passes(d)}, Wegstein ${passes(w)}, Broyden ${passes(b)}; balance ${w.balance.component} ${w.balance.energy}, worst unit ${w.worst.component} ${w.worst.energy}')
	assert passes(d) <= 29 && passes(w) <= 4 && passes(b) <= 9, '${passes(d)} ${passes(w)} ${passes(b)}'
	assert w.plan.blocks.map(it.tears) == [['recycle'], [], []]
	// Within the tear tolerance. Direct substitution stops with up to g/(1 − g)
	// of its last move still to go, which at this loop's gain of 0.5 is tol of
	// the tear's flow, and more of a smaller stream's: 1.6e-9 measured. Wegstein
	// and Broyden land closer, and agree within tol: 6.3e-12 measured.
	same(d, w, 1e-8, 'direct and Wegstein')
	same(d, b, 1e-8, 'direct and Broyden')
	same(w, b, 1e-9, 'Wegstein and Broyden')
	// Every residual history ends within tol, and no sooner.
	for s in [d, w, b] {
		rows := s.history[0].residuals
		assert rows.last()[0] <= 1e-9
		assert rows[rows.len - 2][0] > 1e-9
	}
	// Each unit reports what it took, on the units it ought to.
	work := w.units.filter(it.name == 'pump')[0].energy
	assert work[0].name == 'work' && work[0].rate > 0.0
	assert w.units.filter(it.name == 'cooler')[0].energy[0].rate < 0.0
	assert w.products == ['gas', 'flash gas', 'product']
	assert w.worst.component <= 1e-10 && w.worst.energy <= 1e-10
}

fn test_a_loop_that_is_truly_nonlinear_gets_one_answer_too() {
	// The reference loop is kinder than it looks. Its recycle has the drum
	// liquid's composition at the drum's T and P, so adding it back leaves the
	// equilibrium as it was: the loop is affine, with a gain of exactly 0.5, and
	// Wegstein's one secant lands on the answer. With the cooler on a duty, T
	// moves with the recycle, and the phases with it.
	fl := fluid()!
	mut f := reference(fl)!
	base := f.solve(fl)!
	q := base.units.filter(it.name == 'cooler')[0].energy[0].rate
	f.units = f.units.map(if it.unit.name == 'cooler' { replace_cooler(it, q) } else { it })
	d := solved(fl, f, method: .direct, max_iterations: 200)!
	w := solved(fl, f, method: .wegstein)!
	b := solved(fl, f, method: .broyden)!
	println('with the duty fixed: direct ${passes(d)}, Wegstein ${passes(w)}, Broyden ${passes(b)}')
	assert passes(w) <= 14 && passes(b) <= 16, '${passes(w)} ${passes(b)}'
	assert passes(w) < passes(d) && passes(b) < passes(d)
	// As for the reference flowsheet, whose gain this loop's 29 direct passes
	// share: 1.5e-9 and 5.3e-11 measured.
	same(d, w, 1e-8, 'direct and Wegstein')
	same(d, b, 1e-8, 'direct and Broyden')
	same(w, b, 1e-9, 'Wegstein and Broyden')
}

// replace_cooler is the cooler on a duty, on the same ports.
fn replace_cooler(p flowsheet.Placed, q f64) flowsheet.Placed {
	return flowsheet.Placed{
		unit:    flowsheet.new_heater('cooler', mode: .duty, value: q) or { panic(err) }
		inlets:  p.inlets
		outlets: p.outlets
	}
}

fn test_the_answer_does_not_depend_on_the_tear() {
	// Whichever stream of the loop is torn, the plant is the same. Only 'back'
	// and 'liquid' need an estimate: a gas started there would reach the pump.
	fl := fluid()!
	f := reference(fl)!
	base := solved(fl, f, method: .wegstein, tol: 1e-10)!
	for name in ['mixed', 'cold', 'recycle'] {
		s := solved(fl, f, method: .wegstein, tol: 1e-10, tears: [name])!
		assert s.plan.blocks[0].tears == [name]
		same(base, s, 1e-8, 'torn at ${name}')
	}
	for name in ['liquid', 'back'] {
		f.solve(fl, tears: [name]) or {
			assert err is FlowsheetError, err.msg()
			assert (err as FlowsheetError).unit == 'pump', err.msg()
			assert err.msg().contains('takes liquid'), err.msg()
			s := solved(fl, f,
				method:    .wegstein
				tol:       1e-10
				estimates: {
					name: base.streams[name]
				}
			)!
			assert s.plan.blocks[0].tears == [name]
			assert passes(s) <= 12, '${name}: ${passes(s)}'
			same(base, s, 1e-8, 'torn at ${name}')
			continue
		}
		assert false, 'a gas reached the pump at ${name}'
	}
}

// two_loops is a mixer with two recycles, from two splitters in turn.
fn two_loops(fl thermo.Fluid, a f64, b f64) !Flowsheet {
	mut f := Flowsheet{}
	f.add_feed('feed', flowsheet.feed_tp(fl, flows(10.0), 300.0, 5e5)!)!
	f.add_unit(flowsheet.new_mixer('mixer', 3)!, {
		'in1': 'feed'
		'in2': 'first back'
		'in3': 'second back'
		'out': 'mixed'
	})!
	f.add_unit(flowsheet.new_splitter('first', [a, 1.0 - a])!, {
		'in':   'mixed'
		'out1': 'first back'
		'out2': 'on'
	})!
	f.add_unit(flowsheet.new_splitter('second', [b, 1.0 - b])!, {
		'in':   'on'
		'out1': 'second back'
		'out2': 'product'
	})!
	return f
}

fn test_interlocking_loops_have_their_exact_answer() {
	// The loops share 'mixed', and tearing it breaks both. The mixer carries
	// L = F + aL + b(1 − a)L, so L = F/(1 − a − b(1 − a)): 10/(1 − 0.5 − 0.15)
	// = 28.571… with a = 0.5 and b = 0.3.
	fl := fluid()!
	f := two_loops(fl, 0.5, 0.3)!
	want := 10.0 / (1.0 - 0.5 - 0.3 * 0.5)
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		s := solved(fl, f, method: m)!
		assert s.plan.blocks[0].tears.len == 1
		assert math.abs(s.streams['mixed'].total() - want) <= 1e-6, '${m}'
		assert math.abs(s.streams['product'].total() - 10.0) <= 1e-6, '${m}'
		// Both recycles torn: two sets of variables, converged together.
		t := solved(fl, f, method: m, tears: ['first back', 'second back'])!
		assert t.plan.blocks[0].tears == ['first back', 'second back']
		assert t.history[0].residuals[0].len == 2
		assert math.abs(t.streams['mixed'].total() - want) <= 1e-6, '${m} with two tears'
	}
}

fn test_loops_in_series_solve_in_order() {
	// Two loops, one after the other: two blocks, each with its own history,
	// and the second fed by the first's product.
	fl := fluid()!
	mut f := recycle(fl, 0.5)!
	f.add_unit(flowsheet.new_mixer('mixer2', 2)!, {
		'in1': 'product'
		'in2': 'back2'
		'out': 'mixed2'
	})!
	f.add_unit(flowsheet.new_splitter('split2', [0.25, 0.75])!, {
		'in':   'mixed2'
		'out1': 'back2'
		'out2': 'end'
	})!
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		s := solved(fl, f, method: m)!
		assert s.history.len == 2 && s.history[1].tears == ['back2']
		// Each block's history is its own: the second loop's gain is a quarter, so
		// it's no slower.
		assert s.history[1].residuals.len <= s.history[0].residuals.len, '${m}'
		assert math.abs(s.streams['mixed'].total() - 20.0) <= 1e-6
		assert math.abs(s.streams['mixed2'].total() - 40.0 / 3.0) <= 1e-6
		assert math.abs(s.streams['end'].total() - 10.0) <= 1e-6
	}
}

fn test_a_flowsheet_with_no_loop_passes_once() {
	fl := fluid()!
	mut f := Flowsheet{}
	feed := flowsheet.feed_tp(fl, flows(10.0), 300.0, 30e5)!
	f.add_feed('feed', feed)!
	f.add_unit(flowsheet.new_valve('valve', mode: .outlet_p, value: 5e5)!, {
		'in':  'feed'
		'out': 'product'
	})!
	s := solved(fl, f, flowsheet.SolveParams{})!
	assert s.history.len == 0
	assert s.streams['product'].p == 5e5
	assert s.products == ['product']
	assert s.balance.component <= 1e-12 && s.balance.energy <= 1e-10
}

fn test_an_empty_recycle_converges_at_once() {
	// Nothing goes back, so the first pass is the answer, and the second shows
	// it: the recycle's state settles, though it has no flow.
	fl := fluid()!
	f := recycle(fl, 0.0)!
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		s := solved(fl, f, method: m)!
		assert s.streams['back'].total() == 0.0, '${m}'
		assert math.abs(s.streams['product'].total() - 10.0) <= 1e-12, '${m}'
		assert passes(s) <= 3, '${m}: ${passes(s)}'
	}
}

fn test_a_loop_with_no_steady_state_is_a_typed_failure() {
	// Everything goes back, so nothing leaves and the loop only fills. There's
	// no answer, and no method may claim one. The history comes with the error.
	fl := fluid()!
	f := recycle(fl, 1.0)!
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		f.solve(fl, method: m, max_iterations: 20) or {
			assert err is FlowsheetError, err.msg()
			e := err as FlowsheetError
			assert e.kind == .no_convergence, '${m}: ${err.msg()}'
			assert e.history.len == 1 && e.history[0].tears == ['back']
			assert e.history[0].residuals.len >= 2 && e.history[0].residuals.len <= 20, '${m}'
			continue
		}
		assert false, '${m}: an answer with no steady state'
	}
}

fn test_running_out_of_passes_is_an_error_with_the_history() {
	fl := fluid()!
	f := reference(fl)!
	f.solve(fl, method: .direct, max_iterations: 3) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .no_convergence, err.msg()
		assert e.history[0].residuals.len == 3
		assert e.history[0].residuals[2][0] > 1e-9
		assert err.msg().contains('tears recycle'), err.msg()
		return
	}
	assert false, 'an answer in 3 passes'
}

// Leak loses a tenth of what goes through it.
struct Leak {
	name string
}

fn (l Leak) inlet_ports() []string {
	return ['in']
}

fn (l Leak) outlet_ports() []string {
	return ['out']
}

fn (l Leak) calculate(_ thermo.Fluid, inlets []Stream) !Outcome {
	return Outcome{
		outlets: [Stream{
			...inlets[0]
			flows: inlets[0].flows.map(it * 0.9)
		}]
	}
}

// Warm adds a thousandth to its stream's H and says nothing of it: the moles
// balance, and the energy doesn't.
struct Warm {
	name string
}

fn (w Warm) inlet_ports() []string {
	return ['in']
}

fn (w Warm) outlet_ports() []string {
	return ['out']
}

fn (w Warm) calculate(_ thermo.Fluid, inlets []Stream) !Outcome {
	return Outcome{
		outlets: [Stream{
			...inlets[0]
			h: inlets[0].h + 1e-3 * math.abs(inlets[0].h)
		}]
	}
}

// Turn makes a tenth of its first component into its second, with the same
// total flow and H: the energy balances, and the components don't.
struct Turn {
	name string
}

fn (t Turn) inlet_ports() []string {
	return ['in']
}

fn (t Turn) outlet_ports() []string {
	return ['out']
}

fn (t Turn) calculate(_ thermo.Fluid, inlets []Stream) !Outcome {
	mut flows := inlets[0].flows.clone()
	flows[1] += 0.1 * flows[0]
	flows[0] *= 0.9
	return Outcome{
		outlets: [Stream{
			...inlets[0]
			flows: flows
		}]
	}
}

// Greedy is a unit that gives more outlets than it has ports.
struct Greedy {
	name string
}

fn (g Greedy) inlet_ports() []string {
	return ['in']
}

fn (g Greedy) outlet_ports() []string {
	return ['out']
}

fn (g Greedy) calculate(_ thermo.Fluid, inlets []Stream) !Outcome {
	return Outcome{
		outlets: [inlets[0], inlets[0]]
	}
}

fn one_unit(fl thermo.Fluid, u Unit) !Flowsheet {
	mut f := Flowsheet{}
	f.add_feed('feed', flowsheet.feed_tp(fl, flows(10.0), 300.0, 5e5)!)!
	f.add_unit(u, {
		'in':  'feed'
		'out': 'product'
	})!
	return f
}

fn test_a_unit_that_loses_mass_is_caught() {
	// No Solution is returned unless every balance closes (D8). This is what the
	// check is for: a unit from outside this module with a bug in it.
	fl := fluid()!
	f := one_unit(fl, Leak{'leak'})!
	f.solve(fl) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .unbalanced && e.unit == 'leak', err.msg()
		g := one_unit(fl, Greedy{'greedy'})!
		g.solve(fl) or {
			assert err is FlowsheetError, err.msg()
			f2 := err as FlowsheetError
			assert f2.kind == .bad_specification && f2.unit == 'greedy', err.msg()
			assert err.msg().contains('2 outlets for 1 ports'), err.msg()
			return
		}
		assert false, 'two outlets from one port'
		return
	}
	assert false, 'a unit that loses a tenth'
}

fn test_each_balance_is_checked_on_its_own() {
	// One unit that makes energy and one that makes one compound of another:
	// each is caught by its own balance, the other closing.
	fl := fluid()!
	for u in [Unit(Warm{'warm'}), Turn{'turn'}] {
		f := one_unit(fl, u)!
		f.solve(fl) or {
			e := err as FlowsheetError
			assert e.kind == .unbalanced && e.unit == u.name, err.msg()
			continue
		}
		assert false, '${u.name} was let through'
	}
}

fn test_a_failure_in_a_loop_keeps_its_cause_and_history() {
	// pt_test's water, n-decane and methane at 350 K and 5 bar has three phases.
	// The cooler in this loop reaches it on the first pass.
	fl := db.fluid(.pr, ['water', 'n-decane', 'methane'],
		kij:       false
		overrides: [
			kij_of('7732-18-5', '124-18-5', 0.5),
			kij_of('7732-18-5', '74-82-8', 0.5),
			kij_of('124-18-5', '74-82-8', 0.04),
		]
	)!
	mut f := Flowsheet{}
	f.add_feed('feed', flowsheet.feed_tp(fl, [5.0, 3.0, 2.0], 500.0, 5e5)!)!
	f.add_unit(flowsheet.new_mixer('mixer', 2)!, {
		'in1': 'feed'
		'in2': 'back'
		'out': 'mixed'
	})!
	f.add_unit(flowsheet.new_heater('cooler', mode: .outlet_t, value: 350.0)!, {
		'in':  'mixed'
		'out': 'cold'
	})!
	f.add_unit(flowsheet.new_splitter('split', [0.5, 0.5])!, {
		'in':   'cold'
		'out1': 'back'
		'out2': 'product'
	})!
	f.solve(fl) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .flash_failed && e.unit == 'cooler', err.msg()
		assert (e.cause or { flash.Failure.bad_input }) == .more_phases, err.msg()
		assert e.history.len == 1 && e.history[0].tears == ['back'], err.msg()
		return
	}
	assert false, 'a two-phase answer for a three-phase state'
}

fn kij_of(cas1 string, cas2 string, k f64) components.Kij {
	return components.Kij{
		cas1: cas1
		cas2: cas2
		k:    k
	}
}

fn test_a_bad_call_is_refused() {
	fl := fluid()!
	f := recycle(fl, 0.5)!
	for params in [
		flowsheet.SolveParams{
			tol: 0.0
		},
		flowsheet.SolveParams{
			tol: math.nan()
		},
		flowsheet.SolveParams{
			max_iterations: 0
		},
	] {
		f.solve(fl, params) or {
			assert err.msg().contains('tol and max_iterations'), err.msg()
			continue
		}
		assert false, 'a bad call was accepted'
	}
	// An estimate tears its stream, and that has to be in a loop.
	first := f.solve(fl)!
	f.solve(fl,
		estimates: {
			'product': first.streams['product']
		}
	) or {
		assert err.msg().contains('tear product is not a stream between two units'), err.msg()
		// A flowsheet that isn't right is refused before anything is calculated.
		bad := Flowsheet{}
		bad.solve(fl) or {
			assert err.msg().contains('needs a unit'), err.msg()
			return
		}
		assert false, 'an empty flowsheet solved'
		return
	}
	assert false, 'an estimate for a stream in no loop'
}

fn test_an_estimate_is_where_a_tear_starts() {
	// The answer as its own estimate converges in the fewest passes there are:
	// one, to show it, within what it was converged to.
	fl := fluid()!
	f := reference(fl)!
	s := f.solve(fl, tol: 1e-11)!
	for m in [flowsheet.Method.direct, .wegstein, .broyden] {
		again := solved(fl, f,
			method:    m
			tol:       1e-9
			estimates: {
				'recycle': s.streams['recycle']
			}
		)!
		assert passes(again) == 1, '${m}: ${passes(again)}'
	}
	// A pressure 1% out is the first pass's whole change, as a fraction of P.
	off := Stream{
		...s.streams['recycle']
		p: 31.31e5
	}
	wrong := solved(fl, f,
		tol:       1e-9
		estimates: {
			'recycle': off
		}
	)!
	first := wrong.history[0].residuals[0][0]
	assert math.abs(first - 0.31 / 31.0) <= 1e-8, '${first}'
}

// Gate sends back what's left of its inlet after a fixed amount, and the rest
// out. It has a kink at zero, and above that a slope of 0.9, so an estimate far
// too high falls slowly, and Wegstein's step, q = −5, goes straight past zero.
struct Gate {
	name string
	kept f64
}

fn (g Gate) inlet_ports() []string {
	return ['in']
}

fn (g Gate) outlet_ports() []string {
	return ['back', 'out']
}

fn (g Gate) calculate(_ thermo.Fluid, inlets []Stream) !Outcome {
	s := inlets[0]
	total := s.total()
	back := math.max(0.0, 0.9 * total - g.kept)
	return Outcome{
		outlets: [
			Stream{
				...s
				flows: s.flows.map(it * back / total)
			},
			Stream{
				...s
				flows: s.flows.map(it * (total - back) / total)
			},
		]
	}
}

fn test_a_step_past_zero_is_held_at_zero() {
	// One mol/s in, and a gate keeping 2, so nothing goes back at the answer. From
	// an estimate of 10 the recycle falls to 7.9 and 6.01, and the slope of 0.9
	// has q at −5: −5·7.9 + 6·6.01 is below zero. A flow can't be, so that's
	// held at 0, which here is the answer, with the tear empty.
	fl := fluid()!
	mut f := Flowsheet{}
	feed := flowsheet.feed_tp(fl, flows(1.0), 300.0, 5e5)!
	f.add_feed('feed', feed)!
	f.add_unit(flowsheet.new_mixer('mixer', 2)!, {
		'in1': 'feed'
		'in2': 'back'
		'out': 'mixed'
	})!
	f.add_unit(Gate{'gate', 2.0}, {
		'in':   'mixed'
		'back': 'back'
		'out':  'product'
	})!
	high := flowsheet.feed_tp(fl, flows(10.0), 300.0, 5e5)!
	s := solved(fl, f,
		method:    .wegstein
		estimates: {
			'back': high
		}
	)!
	assert s.streams['back'].total() == 0.0
	assert math.abs(s.streams['product'].total() - 1.0) <= 1e-12
}
