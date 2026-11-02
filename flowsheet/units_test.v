// Streams and units, each alone, by exit criterion 2's checks from scratch:
// every outlet's phases an equilibrium by M3's criterion 3, every unit's
// balances closed to 1e-10 with H recomputed from T, P and the compositions,
// and every outlet at its specification. The comparison with DWSIM is step 8.
import components
import flash
import flowsheet { FlowsheetError, Outcome, Stream, Unit }
import math
import thermo

const db = components.load_default() or { panic(err) }

// The reference flowsheet's feed.
const gases = ['methane', 'ethane', 'propane', 'n-butane', 'n-pentane']
const feed_z = [0.70, 0.10, 0.08, 0.07, 0.05]

fn fluid() !thermo.Fluid {
	return db.fluid(.pr, gases)!
}

fn flows(total f64) []f64 {
	return feed_z.map(it * total)
}

// Molar is a stream's H, S and V per mole.
struct Molar {
mut:
	h f64
	s f64
	v f64
}

// root_of is the volume root a phase sits on. `.stable` won't do: at a pure
// component's saturation point the two roots have the same G, so rounding
// would choose between them.
fn root_of(phase flash.Phase) thermo.Root {
	return if phase.label == .vapor { thermo.Root.vapor } else { thermo.Root.liquid }
}

// from_scratch sums H, S and V over the phases, each phase's volume root
// solved again from T, P and its composition.
fn from_scratch(fl thermo.Fluid, s Stream) !Molar {
	mut m := Molar{}
	for phase in s.result.phases {
		if phase.fraction == 0.0 {
			continue
		}
		st := thermo.solve_tp(fl.model, s.t, s.p, phase.state.n, root_of(phase))!
		pr := fl.props(st)
		m.h += phase.fraction * pr.h
		m.s += phase.fraction * pr.s
		m.v += phase.fraction * st.v
	}
	return m
}

// check_stream checks a stream from scratch: its flows agree with its
// composition, its H, S and V are its phases', and it's an equilibrium by
// M3's criterion 3: fugacities equal to 1e-10, each component's mass balance
// to 1e-12, Σx = Σy = 1, and no phase unstable.
fn check_stream(fl thermo.Fluid, s Stream, at string) ! {
	r := s.result
	assert r.t == s.t && r.p == s.p && r.beta == s.beta, at
	total := s.total()
	if total > 0.0 {
		for i, n in s.flows {
			assert math.abs(n / total - r.z[i]) <= 4e-16 + 1e-15 * r.z[i], '${at}: composition ${i}'
		}
	}
	m := from_scratch(fl, s)!
	// The volume is solved again to 1e-12 or so: H and S follow it.
	assert math.abs(m.h - s.h) <= 1e-9 * thermo.r_gas * s.t, '${at}: H ${s.h} vs ${m.h}'
	assert math.abs(m.s - s.s) <= 1e-9 * thermo.r_gas, '${at}: S ${s.s} vs ${m.s}'
	assert math.abs(m.v - s.v) <= 1e-10 * s.v, '${at}: V ${s.v} vs ${m.v}'
	if r.phases.len == 1 {
		st := flash.stability(fl.model, s.t, s.p, r.z)!
		assert st.stable, '${at}: a single phase that is unstable'
		return
	}
	assert r.phases.len == 2, at
	y, x := r.phases[0].state.n, r.phases[1].state.n
	sy := thermo.solve_tp(fl.model, s.t, s.p, y, root_of(r.phases[0]))!
	sx := thermo.solve_tp(fl.model, s.t, s.p, x, root_of(r.phases[1]))!
	mut sum_x, mut sum_y := 0.0, 0.0
	for i, zi in r.z {
		sum_x += x[i]
		sum_y += y[i]
		assert math.abs(r.phases[0].fraction * y[i] + r.phases[1].fraction * x[i] - zi) <= 1e-12, '${at}: mass balance of component ${i}'
		if zi > 0.0 {
			res := math.abs(math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
			assert res <= 1e-10, '${at}: |ln fV − ln fL| of component ${i} is ${res}'
		}
	}
	assert math.abs(sum_x - 1.0) <= 1e-14 && math.abs(sum_y - 1.0) <= 1e-14, at
	for pair in [[x, y], [y, x]] {
		st := flash.stability(fl.model, s.t, s.p, pair[0], trials: [pair[1]], tol: 1e-8)!
		assert st.stable, '${at}: a phase of the split is unstable'
	}
}

// check_unit checks a unit's outcome from scratch: each outlet by
// check_stream, and its component and energy balances to 1e-10 (criterion
// 2), of the flow and of F·RT.
fn check_unit(fl thermo.Fluid, inlets []Stream, o Outcome, at string) ! {
	mut comp := []f64{len: inlets[0].flows.len}
	mut energy := 0.0
	mut total := 0.0
	mut t := 0.0
	for s in inlets {
		energy += s.total() * from_scratch(fl, s)!.h
		total += s.total()
		t = math.max(t, s.t)
		for i, n in s.flows {
			comp[i] += n
		}
	}
	for e in o.energy {
		energy += e.rate
	}
	for s in o.outlets {
		check_stream(fl, s, at)!
		energy -= s.total() * from_scratch(fl, s)!.h
		for i, n in s.flows {
			comp[i] -= n
		}
	}
	for i, c in comp {
		assert math.abs(c) <= 1e-10 * total, '${at}: the balance of component ${i} is off by ${c} mol/s'
	}
	assert math.abs(energy) <= 1e-10 * total * thermo.r_gas * t, '${at}: the energy balance is off by ${energy} W'
}

// run calculates a unit and checks its outcome.
fn run(fl thermo.Fluid, u Unit, inlets []Stream) !Outcome {
	o := u.calculate(fl, inlets)!
	assert o.outlets.len == u.outlet_ports().len, u.name
	check_unit(fl, inlets, o, u.name)!
	return o
}

// refused asserts that f fails with a FlowsheetError of this kind and unit,
// its message mentioning want.
fn refused(f fn () !, kind flowsheet.Failure, unit string, want string) {
	f() or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == kind && e.unit == unit, err.msg()
		assert err.msg().contains(want), err.msg()
		return
	}
	assert false, 'accepted: ${want}'
}

fn test_a_feed_is_the_same_state_by_each_specification() {
	// D2's four ways of giving a feed. From a two-phase state at (T, P), each
	// other specification comes back to it: flash tolerances, so 1e-8.
	fl := fluid()!
	a := flowsheet.feed_tp(fl, flows(100.0), 230.0, 30e5)!
	assert a.beta > 0.05 && a.beta < 0.95, 'not two-phase: ${a.beta}'
	check_stream(fl, a, 'PT')!
	b := flowsheet.feed_ph(fl, flows(100.0), 30e5, a.h)!
	check_stream(fl, b, 'PH')!
	assert math.abs(b.t - a.t) <= 1e-8 && math.abs(b.beta - a.beta) <= 1e-8
	c := flowsheet.feed_tvf(fl, flows(100.0), 230.0, a.beta)!
	check_stream(fl, c, 'T, β')!
	assert math.abs(c.p - a.p) <= 1e-8 * a.p
	d := flowsheet.feed_pvf(fl, flows(100.0), 30e5, a.beta)!
	check_stream(fl, d, 'P, β')!
	assert math.abs(d.t - a.t) <= 1e-8
}

fn test_a_feed_by_h_finds_a_two_phase_state_the_ideal_gas_puts_near_0_k() {
	// At these H, about 215 K and two-phase, the feed as an ideal gas would be
	// at 11 to 38 K. The PH flash once started there and failed on PR's
	// liquid–liquid splits, with no estimate to save it.
	fl := fluid()!
	for k in 0 .. 7 {
		h := -97000.0 + 250.0 * k
		s := flowsheet.feed_ph(fl, flows(100.0), 30e5, h)!
		check_stream(fl, s, 'H = ${h}')!
		assert s.t > 195.0 && s.t < 225.0, '${s.t}'
		assert s.beta > 0.3 && s.beta < 0.7, '${s.beta}'
		assert math.abs(from_scratch(fl, s)!.h - h) <= 1e-10 * thermo.r_gas * s.t
	}
}

fn test_a_stream_totals_and_keeps_its_own_flows() {
	// The flows belong to the stream: changing the array it was built from
	// doesn't change it.
	fl := fluid()!
	mut f := flows(100.0)
	s := flowsheet.feed_tp(fl, f, 300.0, 30e5)!
	f[0] = 1e6
	assert math.abs(s.total() - 100.0) <= 1e-12
	assert s.flows[0] == 70.0
}

fn test_a_bad_feed_is_refused_for_its_reason() {
	fl := fluid()!
	bad := [
		[1.0, 2.0],
		[-1.0, 1.0, 1.0, 1.0, 1.0],
		[math.nan(), 1.0, 1.0, 1.0, 1.0],
		[math.inf(1), 1.0, 1.0, 1.0, 1.0],
		[0.0, 0.0, 0.0, 0.0, 0.0],
		[1e308, 1e308, 1e308, 1e308, 1e308],
	]
	wants := ['got 2 flows for 5 components', 'component 0', 'component 0', 'component 0', 'above zero',
		'finite total flow']
	for i, f in bad {
		refused(fn [fl, f] () ! {
			_ := flowsheet.feed_tp(fl, f, 300.0, 30e5)!
		}, .bad_specification, '', wants[i])
	}
	// A bad state is the flash's to refuse, and its reason comes through.
	refused(fn [fl] () ! {
		_ := flowsheet.feed_tp(fl, flows(1.0), -5.0, 30e5)!
	}, .flash_failed, '', 'flash: ')
	refused(fn [fl] () ! {
		_ := flowsheet.feed_ph(fl, flows(1.0), -5.0, 0.0)!
	}, .flash_failed, '', 'flash: ')
	refused(fn [fl] () ! {
		_ := flowsheet.feed_tvf(fl, flows(1.0), 300.0, 1.5)!
	}, .flash_failed, '', 'flash: ')
	refused(fn [fl] () ! {
		_ := flowsheet.feed_pvf(fl, flows(1.0), -5.0, 0.5)!
	}, .flash_failed, '', 'flash: ')
}

fn test_units_have_ports_and_hold_no_state() {
	// Through the interface: ports by name, and the same answer each time.
	fl := fluid()!
	a := flowsheet.feed_tp(fl, flows(10.0), 300.0, 30e5)!
	b := flowsheet.feed_tp(fl, flows(20.0), 250.0, 30e5)!
	units := [
		Unit(flowsheet.new_mixer('M', 3)!),
		flowsheet.new_splitter('S', [0.5, 0.25, 0.25])!,
		flowsheet.new_valve('V', mode: .delta_p, value: 1e5)!,
		flowsheet.new_heater('H', mode: .delta_t, value: 5.0)!,
		flowsheet.new_compressor('C', mode: .delta_p, value: 1e5, efficiency: 0.8)!,
		flowsheet.new_expander('E', mode: .delta_p, value: 1e5, efficiency: 0.8)!,
		flowsheet.new_flash_drum('D', mode: .duty, value: 0.0)!,
		flowsheet.new_component_separator('X',
			fractions: [1.0, 0.5, 0.5, 0.5, 0.0]
			t:         300.0
			p:         20e5
		)!,
	]
	assert units.map(it.name) == ['M', 'S', 'V', 'H', 'C', 'E', 'D', 'X']
	assert units[0].inlet_ports() == ['in1', 'in2', 'in3']
	assert units[0].outlet_ports() == ['out']
	assert units[1].inlet_ports() == ['in']
	assert units[1].outlet_ports() == ['out1', 'out2', 'out3']
	assert units[2].inlet_ports() == ['in'] && units[2].outlet_ports() == ['out']
	assert units[3].inlet_ports() == ['in'] && units[3].outlet_ports() == ['out']
	assert units[4].inlet_ports() == ['in'] && units[4].outlet_ports() == ['out']
	assert units[5].inlet_ports() == ['in'] && units[5].outlet_ports() == ['out']
	assert units[6].inlet_ports() == ['in1'] && units[6].outlet_ports() == ['vapour', 'liquid']
	assert units[7].inlet_ports() == ['in'] && units[7].outlet_ports() == ['out1', 'out2']
	three := flowsheet.new_flash_drum('D', inlets: 3, mode: .duty, value: 0.0)!
	assert three.inlet_ports() == ['in1', 'in2', 'in3']
	for u in units {
		inlets := if u.inlet_ports().len == 3 { [a, b, a] } else { [a] }
		first := run(fl, u, inlets)!
		second := u.calculate(fl, inlets)!
		for i, s in first.outlets {
			assert s.t == second.outlets[i].t && s.h == second.outlets[i].h, u.name
		}
		// The wrong number of inlets, or the wrong fluid's.
		refused(fn [fl, u, a] () ! {
			_ := u.calculate(fl, [a, a, a, a])!
		}, .bad_specification, u.name, 'got 4 inlets')
		other := db.fluid(.pr, ['methane', 'ethane'])!
		refused(fn [other, u, inlets] () ! {
			_ := u.calculate(other, inlets)!
		}, .bad_specification, u.name, 'for a fluid of 2')
	}
}

fn test_mixing_one_state_keeps_it() {
	// Two parts of one two-phase state: the outlet is that state. H is held
	// to 1e-10·RT, and C_p in the dome is about 200 J/(mol·K).
	fl := fluid()!
	a := flowsheet.feed_tp(fl, flows(30.0), 230.0, 30e5)!
	b := flowsheet.feed_tp(fl, flows(70.0), 230.0, 30e5)!
	m := flowsheet.new_mixer('M', 2)!
	s := run(fl, m, [a, b])!.outlets[0]
	assert math.abs(s.t - 230.0) <= 1e-8, '${s.t}'
	assert math.abs(s.beta - a.beta) <= 1e-8, '${s.beta} vs ${a.beta}'
	assert s.p == 30e5
	assert math.abs(s.total() - 100.0) <= 1e-12
}

fn test_a_mixer_balances_unlike_inlets_at_the_lowest_pressure() {
	fl := fluid()!
	// A gas into a pentane-rich liquid, which is two-phase mixed.
	gas := flowsheet.feed_tp(fl, [60.0, 5.0, 1.0, 0.5, 0.1], 320.0, 40e5)!
	liquid := flowsheet.feed_tp(fl, [1.0, 1.0, 5.0, 10.0, 25.0], 300.0, 35e5)!
	assert gas.beta == 1.0 && liquid.beta == 0.0, '${gas.beta} ${liquid.beta}'
	s := run(fl, flowsheet.new_mixer('M', 2)!, [gas, liquid])!.outlets[0]
	assert s.p == 35e5
	assert s.beta > 0.0 && s.beta < 1.0, '${s.beta}'
	assert s.t > 290.0 && s.t < 330.0, '${s.t}'
	// The inlets in the other order are the same mixture.
	swapped := run(fl, flowsheet.new_mixer('M', 2)!, [liquid, gas])!.outlets[0]
	assert swapped.p == 35e5
	assert math.abs(swapped.t - s.t) <= 1e-8, '${swapped.t} vs ${s.t}'

	// Three inlets, with the outlet pressure given below all of them.
	third := flowsheet.feed_tp(fl, flows(7.0), 400.0, 50e5)!
	t := run(fl, flowsheet.new_mixer('M3', 3, outlet_p: 20e5)!, [gas, liquid, third])!.outlets[0]
	assert t.p == 20e5

	// Components that one inlet lacks, and one that all do.
	light := flowsheet.feed_tp(fl, [1.0, 0.0, 0.0, 0.0, 0.0], 300.0, 30e5)!
	heavy := flowsheet.feed_tp(fl, [0.0, 0.0, 1.0, 0.0, 2.0], 300.0, 30e5)!
	u := run(fl, flowsheet.new_mixer('M2', 2)!, [light, heavy])!.outlets[0]
	assert u.flows == [1.0, 0.0, 1.0, 0.0, 2.0]
}

fn test_mixing_and_heating_commute() {
	// H is a state function: heat two streams and mix them, or mix them and
	// add both duties. The two outlets agree within the flashes' tolerance.
	fl := fluid()!
	a := flowsheet.feed_tp(fl, flows(40.0), 280.0, 30e5)!
	b := flowsheet.feed_tp(fl, [10.0, 5.0, 5.0, 5.0, 5.0], 260.0, 30e5)!
	qa := -60e3
	qb := -90e3
	ha := run(fl, flowsheet.new_heater('A', mode: .duty, value: qa)!, [a])!.outlets[0]
	hb := run(fl, flowsheet.new_heater('B', mode: .duty, value: qb)!, [b])!.outlets[0]
	m := flowsheet.new_mixer('M', 2)!
	first := run(fl, m, [ha, hb])!.outlets[0]
	mixed := run(fl, m, [a, b])!.outlets[0]
	second := run(fl, flowsheet.new_heater('C', mode: .duty, value: qa + qb)!, [mixed])!.outlets[0]
	assert math.abs(first.t - second.t) <= 1e-8, '${first.t} vs ${second.t}'
	assert math.abs(first.beta - second.beta) <= 1e-8, '${first.beta} vs ${second.beta}'
}

fn test_a_mixer_ignores_an_empty_inlet() {
	// A closed branch at 2 bar neither adds to the flow nor sets the pressure.
	fl := fluid()!
	low := flowsheet.feed_tp(fl, flows(10.0), 300.0, 2e5)!
	split := flowsheet.new_splitter('S', [1.0, 0.0])!
	empty := run(fl, split, [low])!.outlets[1]
	assert empty.total() == 0.0
	a := flowsheet.feed_tp(fl, flows(50.0), 280.0, 30e5)!
	m := flowsheet.new_mixer('M', 2)!
	s := run(fl, m, [a, empty])!.outlets[0]
	assert s.p == 30e5
	assert math.abs(s.t - a.t) <= 1e-8
	assert s.flows == a.flows
	refused(fn [fl, m, empty] () ! {
		_ := m.calculate(fl, [empty, empty])!
	}, .bad_specification, 'M', 'no inlet has any flow')
}

fn test_a_mixer_is_refused_for_its_reason() {
	fl := fluid()!
	refused(fn () ! {
		_ := flowsheet.new_mixer('M', 1)!
	}, .bad_specification, 'M', 'at least 2 inlets')
	refused(fn () ! {
		_ := flowsheet.new_mixer('M', 2, outlet_p: -1.0)!
	}, .bad_specification, 'M', 'outlet pressure')
	refused(fn () ! {
		_ := flowsheet.new_mixer('M', 2, outlet_p: math.inf(1))!
	}, .bad_specification, 'M', 'outlet pressure')
	// A mixer has no way to raise a pressure.
	a := flowsheet.feed_tp(fl, flows(1.0), 300.0, 30e5)!
	b := flowsheet.feed_tp(fl, flows(1.0), 300.0, 20e5)!
	m := flowsheet.new_mixer('M', 2, outlet_p: 25e5)!
	refused(fn [fl, m, a, b] () ! {
		_ := m.calculate(fl, [a, b])!
	}, .bad_specification, 'M', 'above the lowest inlet')
}

fn test_a_splitter_keeps_the_state_and_divides_the_flow() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(100.0), 230.0, 30e5)!
	o := run(fl, flowsheet.new_splitter('S', [0.5, 0.2, 0.3])!, [inlet])!
	for k, f in [0.5, 0.2, 0.3] {
		s := o.outlets[k]
		assert s.t == inlet.t && s.p == inlet.p && s.beta == inlet.beta
		assert s.h == inlet.h && s.s == inlet.s && s.v == inlet.v
		for i, n in s.flows {
			assert math.abs(n - f * inlet.flows[i]) <= 1e-13 * inlet.flows[i]
		}
	}
	assert o.energy.len == 0
	// Mixed back, the outlets are the inlet.
	back := run(fl, flowsheet.new_mixer('M', 3)!, o.outlets)!.outlets[0]
	assert math.abs(back.t - inlet.t) <= 1e-8
	assert math.abs(back.beta - inlet.beta) <= 1e-8
}

fn test_a_splitter_is_refused_for_its_reason() {
	cases := [[1.0], [1.2, -0.2], [0.5, math.nan()], [0.5, 0.4], [0.5, 0.5 + 2e-9]]
	wants := ['at least 2 outlets', 'a split fraction of -0.2', 'a split fraction of nan', 'sum to 0.9',
		'sum to 1.000000002']
	for i, c in cases {
		refused(fn [c] () ! {
			_ := flowsheet.new_splitter('S', c)!
		}, .bad_specification, 'S', wants[i])
	}
	// Rounding in the fractions isn't a mistake, and is taken out.
	s := flowsheet.new_splitter('S', [0.1, 0.2, 0.7 + 1e-12])!
	assert math.abs(s.fractions[0] + s.fractions[1] + s.fractions[2] - 1.0) <= 2.3e-16
}

fn test_a_valve_is_joule_thomson() {
	// A small expansion of a gas falls by μ_JT·ΔP, taken at the midpoint so
	// that the difference is third order. From Fluid.props, not the flash.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 400.0, 30e5)!
	assert inlet.beta == 1.0
	drop := 0.2e5
	s := run(fl, flowsheet.new_valve('V', mode: .delta_p, value: drop)!, [inlet])!.outlets[0]
	mid := flowsheet.feed_tp(fl, flows(10.0), 0.5 * (inlet.t + s.t), inlet.p - 0.5 * drop)!
	mu := fl.props(mid.result.phases[0].state).mu_jt
	assert mu > 0.0
	want := mu * -drop
	got := s.t - inlet.t
	assert math.abs(got - want) <= 1e-6 * math.abs(want), '${got} vs ${want}'
}

fn test_a_valve_keeps_h_and_raises_s() {
	// A bubble-point liquid flashes across the valve: H the same, T lower, and
	// the entropy up, since throttling is irreversible.
	fl := fluid()!
	inlet := flowsheet.feed_pvf(fl, flows(50.0), 30e5, 0.0)!
	v := flowsheet.new_valve('V', mode: .outlet_p, value: 5e5)!
	o := run(fl, v, [inlet])!
	s := o.outlets[0]
	assert o.energy.len == 0
	assert s.p == 5e5
	assert s.beta > 0.05 && s.beta < 0.95, '${s.beta}'
	assert s.t < inlet.t
	assert s.s > inlet.s
	assert math.abs(from_scratch(fl, s)!.h - from_scratch(fl, inlet)!.h) <= 1e-10 * thermo.r_gas * s.t
	// The same drop stated the other way.
	same := run(fl, flowsheet.new_valve('V', mode: .delta_p, value: 25e5)!, [inlet])!.outlets[0]
	assert same.p == s.p && same.t == s.t && same.beta == s.beta
	// No drop is no change.
	none_ := run(fl, flowsheet.new_valve('V', mode: .delta_p, value: 0.0)!, [inlet])!.outlets[0]
	assert math.abs(none_.t - inlet.t) <= 1e-8
	assert math.abs(none_.beta - inlet.beta) <= 1e-8
}

fn test_an_empty_stream_passes_through_a_valve_and_a_heater() {
	// It keeps its composition, so it has a state, and carries no energy.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 300.0, 30e5)!
	empty := run(fl, flowsheet.new_splitter('S', [0.0, 1.0])!, [inlet])!.outlets[0]
	s := run(fl, flowsheet.new_valve('V', mode: .outlet_p, value: 10e5)!, [empty])!.outlets[0]
	assert s.total() == 0.0 && s.p == 10e5
	o := run(fl, flowsheet.new_heater('H', mode: .outlet_t, value: 250.0)!, [empty])!
	assert o.outlets[0].t == 250.0 && o.energy[0].rate == 0.0
	duty := flowsheet.new_heater('H', mode: .duty, value: 1e3)!
	refused(fn [fl, duty, empty] () ! {
		_ := duty.calculate(fl, [empty])!
	}, .bad_specification, 'H', 'a duty on an empty stream')
}

fn test_a_valve_is_refused_for_its_reason() {
	fl := fluid()!
	refused(fn () ! {
		_ := flowsheet.new_valve('V', mode: .outlet_p, value: 0.0)!
	}, .bad_specification, 'V', 'outlet_p')
	refused(fn () ! {
		_ := flowsheet.new_valve('V', mode: .delta_p, value: -1.0)!
	}, .bad_specification, 'V', 'delta_p')
	refused(fn () ! {
		_ := flowsheet.new_valve('V', mode: .delta_p, value: math.nan())!
	}, .bad_specification, 'V', 'delta_p')
	refused(fn () ! {
		_ := flowsheet.new_valve('V', mode: .outlet_p, value: math.inf(1))!
	}, .bad_specification, 'V', 'outlet_p')
	inlet := flowsheet.feed_tp(fl, flows(1.0), 300.0, 30e5)!
	up := flowsheet.new_valve('V', mode: .outlet_p, value: 31e5)!
	refused(fn [fl, up, inlet] () ! {
		_ := up.calculate(fl, [inlet])!
	}, .bad_specification, 'V', "can't raise")
	big := flowsheet.new_valve('V', mode: .delta_p, value: 30e5)!
	refused(fn [fl, big, inlet] () ! {
		_ := big.calculate(fl, [inlet])!
	}, .bad_specification, 'V', 'leaves no pressure')
}

fn test_a_heater_reaches_each_specification() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(100.0), 300.0, 30e5)!

	// Outlet T, with the duty the balance gives.
	a := run(fl, flowsheet.new_heater('H', mode: .outlet_t, value: 230.0)!, [inlet])!
	s := a.outlets[0]
	q := a.energy[0].rate
	assert a.energy[0].name == 'duty'
	assert s.t == 230.0 && s.p == 30e5
	assert q < 0.0

	// The same as a temperature rise.
	b := run(fl, flowsheet.new_heater('H', mode: .delta_t, value: -70.0)!, [inlet])!
	assert b.outlets[0].t == 230.0 && b.energy[0].rate == q

	// That duty comes back to that state. H is held to 1e-10·RT, so the duty
	// to 1e-10 of F·RT.
	c := run(fl, flowsheet.new_heater('H', mode: .duty, value: q)!, [inlet])!
	assert math.abs(c.outlets[0].t - 230.0) <= 1e-8, '${c.outlets[0].t}'
	assert math.abs(c.outlets[0].beta - s.beta) <= 1e-8
	assert math.abs(c.energy[0].rate - q) <= 1e-10 * 100.0 * thermo.r_gas * 230.0

	// And that outlet's β gives that T.
	d := run(fl, flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: s.beta)!, [inlet])!
	assert math.abs(d.outlets[0].t - 230.0) <= 1e-8, '${d.outlets[0].t}'
	assert math.abs(d.energy[0].rate - q) <= 1e-7 * math.abs(q)

	// β at its ends: the bubble and dew points.
	e := run(fl, flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: 1.0)!, [inlet])!
	assert e.outlets[0].beta == 1.0
	g := run(fl, flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: 0.0)!, [inlet])!
	assert g.outlets[0].beta == 0.0
	assert g.outlets[0].t < e.outlets[0].t

	// A zero duty does nothing but drop the pressure, at constant H.
	f := run(fl, flowsheet.new_heater('H', mode: .duty, value: 0.0, pressure_drop: 5e5)!, [inlet])!
	assert f.outlets[0].p == 25e5
	assert f.outlets[0].t < inlet.t
	assert math.abs(f.energy[0].rate) <= 1e-10 * 100.0 * thermo.r_gas * 300.0
}

fn test_a_heater_drops_the_pressure_in_every_mode() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(100.0), 300.0, 30e5)!
	specs := [
		flowsheet.HeaterParams{ mode: .outlet_t, value: 250.0, pressure_drop: 5e5 },
		flowsheet.HeaterParams{ mode: .delta_t, value: -50.0, pressure_drop: 5e5 },
		flowsheet.HeaterParams{ mode: .duty, value: -200e3, pressure_drop: 5e5 },
		flowsheet.HeaterParams{ mode: .outlet_vapour_fraction, value: 0.5, pressure_drop: 5e5 },
	]
	for spec in specs {
		h := flowsheet.new_heater('H', spec)!
		s := run(fl, h, [inlet])!.outlets[0]
		assert s.p == 25e5, '${spec.mode}'
	}
}

fn test_a_heater_is_cp_dt() {
	// A small rise in one phase is F·Cp·ΔT with Cp at the midpoint, so third
	// order in ΔT. Cp is from Fluid.props, not from the flash.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 400.0, 30e5)!
	assert inlet.beta == 1.0
	q := run(fl, flowsheet.new_heater('H', mode: .delta_t, value: 1.0)!, [inlet])!.energy[0].rate
	mid := flowsheet.feed_tp(fl, flows(10.0), 400.5, 30e5)!
	cp := fl.props(mid.result.phases[0].state).cp
	want := 10.0 * cp * 1.0
	assert math.abs(q - want) <= 1e-6 * want, '${q} vs ${want}'
}

fn test_a_pure_component_boils_by_clapeyron() {
	// Propane at 20 bar, from saturated liquid to saturated vapour. The duty
	// is the heat of vaporisation, which Clapeyron's equation gives from the
	// saturation curve alone: ΔH = T·Δv·dP/dT, with dT/dP by central difference
	// over ±0.1%.
	fl := db.fluid(.pr, ['propane'])!
	p := 20e5
	liquid := flowsheet.feed_pvf(fl, [10.0], p, 0.0)!
	v := run(fl, flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: 1.0)!, [liquid])!
	q := v.energy[0].rate
	vapour := v.outlets[0]
	assert math.abs(vapour.t - liquid.t) <= 1e-9, '${vapour.t} vs ${liquid.t}'
	up := flowsheet.feed_pvf(fl, [10.0], p * 1.001, 0.0)!
	down := flowsheet.feed_pvf(fl, [10.0], p * 0.999, 0.0)!
	dt_dp := (up.t - down.t) / (0.002 * p)
	want := 10.0 * liquid.t * (vapour.v - liquid.v) / dt_dp
	assert math.abs(q - want) <= 1e-6 * want, '${q} vs ${want}'

	// Half that duty sits inside the jump at Tsat: half the moles boiled.
	half := run(fl, flowsheet.new_heater('H', mode: .duty, value: 0.5 * q)!, [liquid])!.outlets[0]
	assert math.abs(half.t - liquid.t) <= 1e-9, '${half.t} vs ${liquid.t}'
	assert math.abs(half.beta - 0.5) <= 1e-9, '${half.beta}'
}

fn test_a_heater_is_refused_for_its_reason() {
	fl := fluid()!
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .outlet_t, value: 0.0)!
	}, .bad_specification, 'H', 'outlet_t')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: 1.5)!
	}, .bad_specification, 'H', 'outlet_vapour_fraction')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .outlet_vapour_fraction, value: -0.5)!
	}, .bad_specification, 'H', 'outlet_vapour_fraction')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .duty, value: math.nan())!
	}, .bad_specification, 'H', 'duty')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .delta_t, value: math.inf(1))!
	}, .bad_specification, 'H', 'delta_t')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .duty, value: 0.0, pressure_drop: -1.0)!
	}, .bad_specification, 'H', 'pressure drop')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .duty, value: 0.0, pressure_drop: math.nan())!
	}, .bad_specification, 'H', 'pressure drop')
	refused(fn () ! {
		_ := flowsheet.new_heater('H', mode: .duty, value: 0.0, pressure_drop: math.inf(1))!
	}, .bad_specification, 'H', 'pressure drop')
	inlet := flowsheet.feed_tp(fl, flows(1.0), 300.0, 30e5)!
	big := flowsheet.new_heater('H', mode: .duty, value: 0.0, pressure_drop: 30e5)!
	refused(fn [fl, big, inlet] () ! {
		_ := big.calculate(fl, [inlet])!
	}, .bad_specification, 'H', 'leaves no pressure')
	// Past absolute zero is the flash's to refuse.
	cold := flowsheet.new_heater('H', mode: .delta_t, value: -400.0)!
	refused(fn [fl, cold, inlet] () ! {
		_ := cold.calculate(fl, [inlet])!
	}, .flash_failed, 'H', 'flash: ')
}

fn kij_of(cas1 string, cas2 string, k f64) components.Kij {
	return components.Kij{
		cas1: cas1
		cas2: cas2
		k:    k
	}
}

fn test_a_flash_failure_names_the_unit_and_keeps_its_kind() {
	// pt_test's water, n-decane and methane at 350 K and 5 bar has three
	// phases, which M4 doesn't handle (M6): the unit says so, by name, and
	// the flash's own kind comes through.
	fl := db.fluid(.pr, ['water', 'n-decane', 'methane'],
		kij:       false
		overrides: [
			kij_of('7732-18-5', '124-18-5', 0.5),
			kij_of('7732-18-5', '74-82-8', 0.5),
			kij_of('124-18-5', '74-82-8', 0.04),
		]
	)!
	inlet := flowsheet.feed_tp(fl, [5.0, 3.0, 2.0], 500.0, 5e5)!
	assert inlet.beta == 1.0
	// An enthalpy no temperature reaches fails in the mixer's and the valve's
	// PH flash, and each names itself.
	unreachable := Stream{
		...inlet
		h: 1e30
	}
	m := flowsheet.new_mixer('mixer', 2)!
	refused(fn [fl, m, inlet, unreachable] () ! {
		_ := m.calculate(fl, [inlet, unreachable])!
	}, .flash_failed, 'mixer', 'flash: ')
	v := flowsheet.new_valve('valve', mode: .delta_p, value: 1e5)!
	refused(fn [fl, v, unreachable] () ! {
		_ := v.calculate(fl, [unreachable])!
	}, .flash_failed, 'valve', 'flash: ')
	lost := Stream{
		...inlet
		s: 1e30
	}
	pump := flowsheet.new_pump('pump', mode: .delta_p, value: 1e5, efficiency: 0.8)!
	liquid := flowsheet.feed_tp(fl, [5.0, 0.0, 0.0], 300.0, 5e5)!
	refused(fn [fl, pump, liquid] () ! {
		_ := pump.calculate(fl, [Stream{
			...liquid
			h: 1e30
		}])!
	}, .flash_failed, 'pump', 'flash: ')
	c := flowsheet.new_compressor('compressor', mode: .delta_p, value: 1e5, efficiency: 0.8)!
	refused(fn [fl, c, lost] () ! {
		_ := c.calculate(fl, [lost])!
	}, .flash_failed, 'compressor', 'flash: ')
	ex := flowsheet.new_expander('expander', mode: .delta_p, value: 1e5, efficiency: 0.8)!
	refused(fn [fl, ex, lost] () ! {
		_ := ex.calculate(fl, [lost])!
	}, .flash_failed, 'expander', 'flash: ')
	d := flowsheet.new_flash_drum('drum', mode: .outlet_t, value: 350.0)!
	refused(fn [fl, d, inlet] () ! {
		_ := d.calculate(fl, [inlet])!
	}, .flash_failed, 'drum', 'flash: ')
	x := flowsheet.new_component_separator('separator',
		fractions: [1.0, 1.0, 1.0]
		t:         350.0
		p:         5e5
	)!
	refused(fn [fl, x, inlet] () ! {
		_ := x.calculate(fl, [inlet])!
	}, .flash_failed, 'separator', 'flash: ')
	h := flowsheet.new_heater('cooler', mode: .outlet_t, value: 350.0)!
	if _ := h.calculate(fl, [inlet]) {
		assert false, 'a two-phase answer for a three-phase state'
	} else {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .flash_failed && e.unit == 'cooler', err.msg()
		assert (e.cause or { flash.Failure.bad_input }) == .more_phases, err.msg()
		assert err.msg().starts_with('flowsheet: cooler: flash: '), err.msg()
	}
}

fn test_a_pump_is_v_dp_over_eta() {
	// D4's model, then the thermodynamics under it. Ideally dH = v dP at
	// constant S. D4 takes v at the inlet, so with η = 1 what's left over is
	// T·ΔS = v_in·ΔP − ∫v dP, or ½(v_in − v_out)·ΔP to third order. That
	// compares S from the Fluid with v from the volume root. The third-order
	// term is 3.2e-3 of it at this ΔP.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(50.0), 150.0, 30e5)!
	assert inlet.beta == 0.0
	dp := 10e5
	o := run(fl, flowsheet.new_pump('P', mode: .delta_p, value: dp, efficiency: 0.75)!, [inlet])!
	s := o.outlets[0]
	assert o.energy[0].name == 'work'
	want := 50.0 * inlet.v * dp / 0.75
	assert math.abs(o.energy[0].rate - want) <= 1e-10 * 50.0 * thermo.r_gas * 150.0
	assert s.p == 40e5
	assert s.t > inlet.t, 'a pump warms what it pumps'
	assert s.s > inlet.s

	// η = 1.
	i := run(fl, flowsheet.new_pump('P', mode: .delta_p, value: dp, efficiency: 1.0)!, [inlet])!.outlets[0]
	lost := 0.5 * (i.t + inlet.t) * (i.s - inlet.s)
	left := 0.5 * (inlet.v - i.v) * dp
	assert left > 0.0
	assert math.abs(lost - left) <= 5e-3 * left, '${lost} vs ${left}'

	// The same pressure stated the other way.
	same := run(fl, flowsheet.new_pump('P', mode: .outlet_p, value: 40e5, efficiency: 0.75)!, [inlet])!
	assert same.outlets[0].t == s.t && same.energy[0].rate == o.energy[0].rate
	// No rise is no work.
	none_ := run(fl, flowsheet.new_pump('P', mode: .delta_p, value: 0.0, efficiency: 0.75)!, [inlet])!
	assert math.abs(none_.energy[0].rate) <= 1e-10 * 50.0 * thermo.r_gas * 150.0
}

fn test_a_pump_takes_liquid_only() {
	fl := fluid()!
	p := flowsheet.new_pump('P', mode: .delta_p, value: 1e5, efficiency: 0.8)!
	assert Unit(p).name == 'P'
	assert p.inlet_ports() == ['in'] && p.outlet_ports() == ['out']
	gas := flowsheet.feed_tp(fl, flows(1.0), 400.0, 30e5)!
	refused(fn [fl, p, gas] () ! {
		_ := p.calculate(fl, [gas])!
	}, .bad_specification, 'P', 'takes liquid')
	wet := flowsheet.feed_tp(fl, flows(1.0), 230.0, 30e5)!
	refused(fn [fl, p, wet] () ! {
		_ := p.calculate(fl, [wet])!
	}, .bad_specification, 'P', 'takes liquid')
	// A bubble-point liquid is liquid, and so is what a drum sends out.
	bubble := flowsheet.feed_pvf(fl, flows(10.0), 30e5, 0.0)!
	assert run(fl, p, [bubble])!.outlets[0].p == 31e5
	drum := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0)!
	liquid := run(fl, drum, [wet])!.outlets[1]
	assert run(fl, p, [liquid])!.outlets[0].p == 31e5
	// An empty inlet holds nothing to refuse, whatever its state, and takes no
	// work.
	empty := run(fl, flowsheet.new_splitter('S', [0.0, 1.0])!, [gas])!.outlets[0]
	assert run(fl, p, [empty])!.energy[0].rate == 0.0
	refused(fn [fl, p, bubble] () ! {
		_ := p.calculate(fl, [bubble, bubble])!
	}, .bad_specification, 'P', 'got 2 inlets')
	// A pump doesn't lower pressure.
	low := flowsheet.new_pump('P', mode: .outlet_p, value: 1e5, efficiency: 0.8)!
	refused(fn [fl, low, bubble] () ! {
		_ := low.calculate(fl, [bubble])!
	}, .bad_specification, 'P', "can't lower")
}

fn test_a_compressor_is_isentropic_over_eta() {
	// Criterion 2: at η = 1 the outlet has the inlet's S. Below that, ΔH is
	// ΔH_s/η, and the outlet is hotter than the isentropic one. PS holds S to
	// 1e-10·R and PH holds H to 1e-10·RT.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 350.0, 5e5)!
	assert inlet.beta == 1.0
	best := run(fl, flowsheet.new_compressor('C', mode: .outlet_p, value: 20e5, efficiency: 1.0)!,
		[inlet])!
	s := best.outlets[0]
	assert s.p == 20e5
	assert best.energy[0].name == 'work' && best.energy[0].rate > 0.0
	assert math.abs(from_scratch(fl, s)!.s - from_scratch(fl, inlet)!.s) <= 1e-9 * thermo.r_gas
	real := run(fl, flowsheet.new_compressor('C', mode: .delta_p, value: 15e5, efficiency: 0.7)!,
		[inlet])!
	r := real.outlets[0]
	assert r.p == 20e5
	assert math.abs(real.energy[0].rate - best.energy[0].rate / 0.7) <= 1e-9 * 10.0 * thermo.r_gas * r.t
	assert r.t > s.t && r.s > inlet.s
}

fn test_isentropic_work_is_the_integral_of_v_dp() {
	// dH = v dP at constant S, so ΔH_s = ∫v dP along the isentrope. H comes
	// from the ideal gas and F's T-derivative, and v from the volume root, so
	// this checks one against the other. Simpson in ln P over 8 panels, each
	// point an isentropic compression from the inlet.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 350.0, 5e5)!
	n := 8
	h := math.log(4.0) / n
	mut sum := 0.0
	mut work := 0.0
	for k in 0 .. n + 1 {
		p := 5e5 * math.exp(k * h)
		o := run(fl, flowsheet.new_compressor('C', mode: .outlet_p, value: p, efficiency: 1.0)!,
			[inlet])!
		weight := if k == 0 || k == n {
			1.0
		} else if k % 2 == 1 {
			4.0
		} else {
			2.0
		}
		sum += weight * o.outlets[0].v * p
		work = o.energy[0].rate / 10.0
	}
	integral := sum * h / 3.0
	assert math.abs(integral - work) <= 1e-6 * work, '${integral} vs ${work}'
}

fn test_an_expander_gives_up_eta_of_the_isentropic_drop() {
	// η·ΔH_s, so it takes less work than the isentropic expansion and leaves
	// warmer. This one's outlet is two-phase, as a turboexpander's can be.
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 340.0, 60e5)!
	assert inlet.beta == 1.0, '${inlet.beta}'
	best := run(fl, flowsheet.new_expander('E', mode: .outlet_p, value: 10e5, efficiency: 1.0)!,
		[inlet])!
	s := best.outlets[0]
	assert s.p == 10e5
	assert s.beta > 0.05 && s.beta < 0.99, '${s.beta}'
	assert best.energy[0].name == 'work' && best.energy[0].rate < 0.0
	assert math.abs(from_scratch(fl, s)!.s - from_scratch(fl, inlet)!.s) <= 1e-9 * thermo.r_gas
	real := run(fl, flowsheet.new_expander('E', mode: .delta_p, value: 50e5, efficiency: 0.8)!,
		[inlet])!
	r := real.outlets[0]
	assert r.p == 10e5
	assert math.abs(real.energy[0].rate - 0.8 * best.energy[0].rate) <= 1e-9 * 10.0 * thermo.r_gas * r.t
	assert r.t > s.t && r.s > inlet.s
	// An expander drops more than a valve: it takes the work out.
	v := run(fl, flowsheet.new_valve('V', mode: .outlet_p, value: 10e5)!, [inlet])!.outlets[0]
	assert v.t > r.t
}

fn test_compressing_and_expanding_isentropically_is_reversible() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(10.0), 350.0, 5e5)!
	up := run(fl, flowsheet.new_compressor('C', mode: .delta_p, value: 15e5, efficiency: 1.0)!,
		[inlet])!
	down := run(fl, flowsheet.new_expander('E', mode: .delta_p, value: 15e5, efficiency: 1.0)!,
		up.outlets)!
	assert math.abs(down.outlets[0].t - inlet.t) <= 1e-7, '${down.outlets[0].t}'
	assert math.abs(up.energy[0].rate + down.energy[0].rate) <= 1e-9 * 10.0 * thermo.r_gas * 350.0
	// An empty stream passes through both, keeping its composition, and takes
	// no work.
	empty := run(fl, flowsheet.new_splitter('S', [0.0, 1.0])!, [inlet])!.outlets[0]
	c := run(fl, flowsheet.new_compressor('C', mode: .delta_p, value: 15e5, efficiency: 0.8)!,
		[empty])!
	assert c.outlets[0].total() == 0.0 && c.outlets[0].p == 20e5 && c.energy[0].rate == 0.0
	e := run(fl, flowsheet.new_expander('E', mode: .delta_p, value: 1e5, efficiency: 0.8)!, [
		empty,
	])!
	assert e.outlets[0].total() == 0.0 && e.outlets[0].p == 4e5 && e.energy[0].rate == 0.0
}

fn test_a_machine_is_refused_for_its_reason() {
	fl := fluid()!
	bad := [
		flowsheet.MachineParams{
			mode:       .outlet_p
			value:      0.0
			efficiency: 0.8
		},
		flowsheet.MachineParams{
			mode:       .delta_p
			value:      -1.0
			efficiency: 0.8
		},
		flowsheet.MachineParams{
			mode:       .delta_p
			value:      math.nan()
			efficiency: 0.8
		},
		flowsheet.MachineParams{
			mode:       .outlet_p
			value:      math.inf(1)
			efficiency: 0.8
		},
		flowsheet.MachineParams{
			mode:       .delta_p
			value:      1e5
			efficiency: 0.0
		},
		flowsheet.MachineParams{
			mode:       .delta_p
			value:      1e5
			efficiency: 1.5
		},
		flowsheet.MachineParams{
			mode:       .delta_p
			value:      1e5
			efficiency: math.nan()
		},
	]
	wants := ['outlet_p', 'delta_p', 'delta_p', 'outlet_p', 'an efficiency of', 'an efficiency of',
		'an efficiency of']
	for i, spec in bad {
		refused(fn [spec] () ! {
			_ := flowsheet.new_pump('M', spec)!
		}, .bad_specification, 'M', wants[i])
		refused(fn [spec] () ! {
			_ := flowsheet.new_compressor('M', spec)!
		}, .bad_specification, 'M', wants[i])
		refused(fn [spec] () ! {
			_ := flowsheet.new_expander('M', spec)!
		}, .bad_specification, 'M', wants[i])
	}
	inlet := flowsheet.feed_tp(fl, flows(1.0), 350.0, 5e5)!
	down := flowsheet.new_compressor('C', mode: .outlet_p, value: 1e5, efficiency: 0.8)!
	refused(fn [fl, down, inlet] () ! {
		_ := down.calculate(fl, [inlet])!
	}, .bad_specification, 'C', "can't lower")
	up := flowsheet.new_expander('E', mode: .outlet_p, value: 10e5, efficiency: 0.8)!
	refused(fn [fl, up, inlet] () ! {
		_ := up.calculate(fl, [inlet])!
	}, .bad_specification, 'E', "can't raise")
	big := flowsheet.new_expander('E', mode: .delta_p, value: 5e5, efficiency: 0.8)!
	refused(fn [fl, big, inlet] () ! {
		_ := big.calculate(fl, [inlet])!
	}, .bad_specification, 'E', 'leaves no pressure')
}

fn test_a_drum_sends_each_phase_out_saturated() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(100.0), 230.0, 30e5)!
	o := run(fl, flowsheet.new_flash_drum('D', mode: .outlet_t, value: 230.0)!, [inlet])!
	v, l := o.outlets[0], o.outlets[1]
	assert v.t == 230.0 && l.t == 230.0 && v.p == 30e5 && l.p == 30e5
	assert v.beta == 1.0 && l.beta == 0.0
	assert math.abs(v.total() - 100.0 * inlet.beta) <= 1e-9
	// Each outlet is the split: the same flash, so no heat.
	assert math.abs(o.energy[0].rate) <= 1e-10 * 100.0 * thermo.r_gas * 230.0
	assert o.energy[0].name == 'duty'
	// Each is on its phase boundary, shown by flashes that know nothing of the
	// drum: the liquid's bubble point and the vapour's dew point at 230 K are
	// the drum's pressure.
	bubble := flowsheet.feed_tvf(fl, l.flows, 230.0, 0.0)!
	dew := flowsheet.feed_tvf(fl, v.flows, 230.0, 1.0)!
	assert math.abs(bubble.p - 30e5) <= 1e-8 * 30e5, '${bubble.p}'
	assert math.abs(dew.p - 30e5) <= 1e-8 * 30e5, '${dew.p}'
}

fn test_a_drum_with_two_inlets_balances() {
	// DWSIM's vessel gets this wrong, flashing at the mean T or creating 62.7
	// kW (docs/spike-dwsim-flowsheet.md). Here the drum is the mixer it
	// contains: the same state as a mixer and then a one-inlet drum.
	fl := fluid()!
	gas := flowsheet.feed_tp(fl, [60.0, 5.0, 1.0, 0.5, 0.1], 320.0, 40e5)!
	liquid := flowsheet.feed_tp(fl, [1.0, 1.0, 5.0, 10.0, 25.0], 300.0, 35e5)!
	o := run(fl, flowsheet.new_flash_drum('D', inlets: 2, mode: .duty, value: 0.0)!, [
		gas,
		liquid,
	])!
	v := o.outlets[0]
	assert v.p == 35e5
	assert math.abs(o.energy[0].rate) <= 1e-10 * 108.6 * thermo.r_gas * 320.0
	m := run(fl, flowsheet.new_mixer('M', 2)!, [gas, liquid])!.outlets[0]
	assert math.abs(v.t - m.t) <= 1e-8, '${v.t} vs ${m.t}'
	assert math.abs(v.total() - m.beta * m.total()) <= 1e-7
	d := run(fl, flowsheet.new_flash_drum('D', mode: .duty, value: 0.0)!, [m])!
	assert math.abs(d.outlets[0].t - v.t) <= 1e-8
	// A drum with a duty is a heater and then a drum.
	q := -40e3
	hot := run(fl, flowsheet.new_flash_drum('D', inlets: 2, mode: .duty, value: q)!, [
		gas,
		liquid,
	])!
	assert math.abs(hot.energy[0].rate - q) <= 1e-10 * 108.6 * thermo.r_gas * 320.0
	cooled := run(fl, flowsheet.new_heater('H', mode: .duty, value: q)!, [m])!.outlets[0]
	assert math.abs(hot.outlets[0].t - cooled.t) <= 1e-8
	assert math.abs(hot.outlets[0].total() - cooled.beta * cooled.total()) <= 1e-7
	// At a given pressure, below all the inlets.
	low := run(fl, flowsheet.new_flash_drum('D',
		inlets:   2
		mode:     .outlet_t
		value:    280.0
		outlet_p: 20e5
	)!, [gas, liquid])!
	assert low.outlets[0].p == 20e5 && low.outlets[1].p == 20e5
	// With a duty there, it's a heater with that pressure drop.
	lowq := run(fl, flowsheet.new_flash_drum('D',
		inlets:   2
		mode:     .duty
		value:    q
		outlet_p: 20e5
	)!, [gas, liquid])!
	assert lowq.outlets[0].p == 20e5
	dropped := run(fl, flowsheet.new_heater('H', mode: .duty, value: q, pressure_drop: 15e5)!,
		[m])!.outlets[0]
	assert math.abs(lowq.outlets[0].t - dropped.t) <= 1e-8
	// An empty inlet is left out, as in a mixer.
	far := flowsheet.feed_tp(fl, flows(1.0), 300.0, 2e5)!
	empty := run(fl, flowsheet.new_splitter('S', [1.0, 0.0])!, [far])!.outlets[1]
	e := run(fl, flowsheet.new_flash_drum('D', inlets: 2, mode: .outlet_t, value: 280.0)!, [
		gas,
		empty,
	])!
	assert e.outlets[0].p == 40e5
}

fn test_a_drum_with_one_phase_has_an_empty_outlet() {
	fl := fluid()!
	gas := flowsheet.feed_tp(fl, flows(10.0), 400.0, 30e5)!
	drum := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0)!
	a := run(fl, drum, [gas])!
	assert a.outlets[0].total() > 9.999 && a.outlets[1].total() == 0.0
	cold := flowsheet.feed_tp(fl, flows(10.0), 150.0, 30e5)!
	b := run(fl, drum, [cold])!
	assert b.outlets[0].total() == 0.0 && b.outlets[1].total() > 9.999
	// The empty outlet keeps a state, so it can go on through the flowsheet.
	v := flowsheet.new_valve('V', mode: .delta_p, value: 1e5)!
	assert run(fl, v, [a.outlets[1]])!.outlets[0].p == 29e5
}

fn test_a_drum_boils_a_pure_component() {
	// Propane at 20 bar, half boiled: β from the energy balance at Tsat. The
	// vapour and the liquid are pure, and differ only in their root.
	fl := db.fluid(.pr, ['propane'])!
	liquid := flowsheet.feed_pvf(fl, [10.0], 20e5, 0.0)!
	vapour := flowsheet.feed_pvf(fl, [10.0], 20e5, 1.0)!
	q := 10.0 * (vapour.h - liquid.h) / 2.0
	o := run(fl, flowsheet.new_flash_drum('D', mode: .duty, value: q)!, [liquid])!
	assert math.abs(o.outlets[0].total() - 5.0) <= 1e-8, '${o.outlets[0].total()}'
	assert math.abs(o.outlets[1].total() - 5.0) <= 1e-8
	assert math.abs(o.outlets[0].h - vapour.h) <= 1e-8 * thermo.r_gas * liquid.t
	assert math.abs(o.outlets[1].h - liquid.h) <= 1e-8 * thermo.r_gas * liquid.t
}

fn test_a_drum_is_refused_for_its_reason() {
	fl := fluid()!
	refused(fn () ! {
		_ := flowsheet.new_flash_drum('D', inlets: 0, mode: .duty, value: 0.0)!
	}, .bad_specification, 'D', 'needs an inlet')
	refused(fn () ! {
		_ := flowsheet.new_flash_drum('D', mode: .outlet_t, value: 0.0)!
	}, .bad_specification, 'D', 'outlet_t')
	refused(fn () ! {
		_ := flowsheet.new_flash_drum('D', mode: .duty, value: math.nan())!
	}, .bad_specification, 'D', 'duty')
	refused(fn () ! {
		_ := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0, outlet_p: -1.0)!
	}, .bad_specification, 'D', 'outlet pressure')
	refused(fn () ! {
		_ := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0, outlet_p: math.inf(1))!
	}, .bad_specification, 'D', 'outlet pressure')
	inlet := flowsheet.feed_tp(fl, flows(1.0), 300.0, 30e5)!
	high := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0, outlet_p: 31e5)!
	refused(fn [fl, high, inlet] () ! {
		_ := high.calculate(fl, [inlet])!
	}, .bad_specification, 'D', 'above the lowest inlet')
	empty := run(fl, flowsheet.new_splitter('S', [1.0, 0.0])!, [inlet])!.outlets[1]
	drum := flowsheet.new_flash_drum('D', mode: .duty, value: 0.0)!
	refused(fn [fl, drum, empty] () ! {
		_ := drum.calculate(fl, [empty])!
	}, .bad_specification, 'D', 'no inlet has any flow')
}

fn test_a_component_separator_splits_by_component() {
	fl := fluid()!
	inlet := flowsheet.feed_tp(fl, flows(100.0), 300.0, 30e5)!
	fractions := [1.0, 0.9, 0.5, 0.1, 0.0]
	o := run(fl, flowsheet.new_component_separator('X',
		fractions: fractions
		t:         250.0
		p:         10e5
	)!, [inlet])!
	a, b := o.outlets[0], o.outlets[1]
	for i, f in fractions {
		assert math.abs(a.flows[i] - f * inlet.flows[i]) <= 1e-13 * inlet.flows[i]
	}
	assert a.t == 250.0 && b.t == 250.0 && a.p == 10e5 && b.p == 10e5
	assert o.energy[0].name == 'duty'
	// The duty is the balance's, here worked out by hand.
	want := a.total() * a.h + b.total() * b.h - inlet.total() * inlet.h
	assert o.energy[0].rate == want

	// A component held back entirely leaves an outlet empty, with the inlet's
	// composition, so it has a state.
	all := run(fl, flowsheet.new_component_separator('X',
		fractions: [1.0, 1.0, 1.0, 1.0, 1.0]
		t:         300.0
		p:         10e5
	)!, [inlet])!
	assert all.outlets[1].total() == 0.0 && all.outlets[1].t == 300.0
	assert all.outlets[0].flows == inlet.flows

	// It keeps its own fractions: changing the array afterwards changes nothing.
	mut given := [1.0, 1.0, 1.0, 1.0, 1.0]
	kept := flowsheet.new_component_separator('X', fractions: given, t: 300.0, p: 10e5)!
	given[0] = 0.0
	assert kept.fractions[0] == 1.0

	// It's a black box, so its outlets can be above the inlet's pressure.
	up := run(fl, flowsheet.new_component_separator('X',
		fractions: fractions
		t:         300.0
		p:         40e5
	)!, [inlet])!
	assert up.outlets[0].p == 40e5
}

fn test_a_component_separator_is_refused_for_its_reason() {
	fl := fluid()!
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [1.5, 0.0], t: 300.0, p: 1e5)!
	}, .bad_specification, 'X', 'a fraction of 1.5')
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [-0.5, 0.0], t: 300.0, p: 1e5)!
	}, .bad_specification, 'X', 'a fraction of -0.5')
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [math.nan()], t: 300.0, p: 1e5)!
	}, .bad_specification, 'X', 'a fraction of nan')
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [0.5], t: 0.0, p: 1e5)!
	}, .bad_specification, 'X', 'an outlet at')
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [0.5], t: 300.0, p: -1.0)!
	}, .bad_specification, 'X', 'an outlet at')
	refused(fn () ! {
		_ := flowsheet.new_component_separator('X', fractions: [0.5], t: math.inf(1), p: 1e5)!
	}, .bad_specification, 'X', 'an outlet at')
	short := flowsheet.new_component_separator('X', fractions: [0.5, 0.5], t: 300.0, p: 1e5)!
	inlet := flowsheet.feed_tp(fl, flows(1.0), 300.0, 30e5)!
	refused(fn [fl, short, inlet] () ! {
		_ := short.calculate(fl, [inlet])!
	}, .bad_specification, 'X', 'got 2 fractions for 5 components')
}
