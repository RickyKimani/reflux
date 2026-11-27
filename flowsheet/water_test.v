// Water by IAPWS-95 in a file (M5a step 7): `equation = "IAPWS-95"` in a
// format 4 [model], its refusals at their lines, and the Rankine cycle of
// criterion 8 against its hand calculation (testdata/m5_rankine.json, from
// the iapws package's states by tools/gen_m5_rankine.py).
import components
import flowsheet { FlowsheetError, Loaded }
import math
import os
import water.if97
import json2

const db = components.load_default() or { panic(err) }

fn rankine_text() string {
	return os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'rankine.toml')) or {
		panic(err)
	}
}

fn load(text string) !Loaded {
	return flowsheet.load_text(text, db: db)!
}

// with is the Rankine file with a text replaced, which must be there once.
fn with(a string, b string) string {
	t := rankine_text()
	assert t.count(a) == 1, a
	return t.replace(a, b)
}

// lined is the 1-based line of a text's first occurrence in the Rankine file.
fn lined(text string) int {
	return rankine_text().all_before(text).count('\n') + 1
}

fn refused(text string, line int, key string, want string, suggestion string) {
	load(text) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.line == line, 'line ${e.line}, not ${line}: ${err.msg()}'
		assert e.key == key, 'key ${e.key}, not ${key}: ${err.msg()}'
		assert e.reason.contains(want), err.msg()
		if suggestion != '' {
			assert suggestion in e.suggestions, err.msg()
		}
		return
	}
	assert false, 'accepted: ${want}'
}

fn number(m map[string]json2.Any, key string) f64 {
	v := m[key] or { panic('no ${key}') }
	return v.f64()
}

fn object(m map[string]json2.Any, key string) map[string]json2.Any {
	v := m[key] or { panic('no ${key}') }
	return v.as_map()
}

fn power(sol flowsheet.Solution, unit string) f64 {
	for u in sol.units {
		if u.name == unit {
			return u.energy[0].rate
		}
	}
	panic('no unit ${unit}')
}

// The cycle by hand: per kg on IAPWS's reference, the pump's v1·ΔP/η, the
// boiler's and condenser's changes of h, and the turbine's η times its
// isentropic drop, each state from the iapws package's IAPWS-95 (a wet one
// by the lever rule on its saturated phases). reflux's powers are molar,
// on its own reference, so each is compared as a change per kg, at the
// file's 50 mol/s of IAPWS's molar mass.
//
// Measured: every power within 1.9e-13 of the hand calculation's,
// relative (1.8e-10 kJ/kg), the feedwater's T within 3.8e-12 K, the
// exhaust's T within 1.1e-12 K and its x within 2.3e-15. The bounds are
// the flashes' tolerance: each outlet's PH or PS flash meets its energy
// balance to 1e-10·RT, which is 3.6e-8 kJ/kg at the boiler's 773 K,
// 3.5e-9 K of the feedwater's T and 6e-12 of the exhaust's x; two such
// flashes go into each power.
fn test_the_rankine_cycle_against_its_hand_calculation() {
	l := load(rankine_text())!
	sol := l.solve()!
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5_rankine.json'))!
	want := (json2.decode[json2.Any](text)!).as_map()
	flow := number(want, 'flow_kg_s')
	per_kg := object(want, 'per_kg')
	states := object(want, 'states')
	mut worst := 0.0 // kJ/kg
	for unit, key in {
		'pump':      'pump'
		'boiler':    'boiler'
		'turbine':   'turbine'
		'condenser': 'condenser'
	} {
		got := power(sol, unit) / flow * 1e-3 // kJ/kg
		w := number(per_kg, key)
		// The turbine's and pump's work are what the cycle gets and gives:
		// the expander's energy stream is the work done on it, so negative.
		sign := if unit == 'turbine' { -1.0 } else { 1.0 }
		d := math.abs(sign * got - w)
		worst = math.max(worst, d)
		assert d <= 1e-7, '${unit}: ${sign * got} kJ/kg against ${w}'
	}
	println('worst power: ${worst:.2e} kJ/kg')
	feedwater := sol.streams['feedwater'] or { panic('no feedwater') }
	f := object(states, 'feedwater')
	assert math.abs(feedwater.t - number(f, 'T')) <= 1e-8, 'feedwater: ${feedwater.t} K'
	// The exhaust is wet, at Tsat(10 kPa): the generator's Newton on the
	// package's Psat(T) and reflux's VF flash, each to its rounding.
	exhaust := sol.streams['exhaust'] or { panic('no exhaust') }
	e := object(states, 'exhaust')
	assert math.abs(exhaust.t - number(e, 'T')) <= 1e-9, 'exhaust: ${exhaust.t} K'
	assert math.abs(exhaust.beta - number(e, 'x')) <= 2e-11, 'exhaust: x ${exhaust.beta}'
	println('feedwater T ${feedwater.t - number(f, 'T'):.2e}, exhaust T ${exhaust.t - number(e, 'T'):.2e}, x ${exhaust.beta - number(e, 'x'):.2e}')
	// The cycle's efficiency, (w_turbine − w_pump)/q_boiler.
	eff := (-power(sol, 'turbine') - power(sol, 'pump')) / power(sol, 'boiler')
	assert math.abs(eff - number(want, 'efficiency')) <= 1e-10, 'efficiency ${eff}'
	// The condenser brings the exhaust back to the condensate: the cycle
	// closes on the feed's state.
	ret := sol.streams['return'] or { panic('no return') }
	feed := sol.streams['condensate'] or { panic('no condensate') }
	assert math.abs(ret.h - feed.h) <= 1e-9 * math.abs(feed.h), 'return: ${ret.h} J/mol against ${feed.h}'
	// The report names the model.
	rep := l.report(sol)
	m := rep.model or { panic('no model') }
	assert m.equation == 'IAPWS-95' && m.pairs.len == 0
}

// What IAPWS-95 can't be is refused at its line: a file before format 4, a
// cubic's keys, another compound or a mixture.
fn test_iapws95_in_a_file_is_refused_at_its_line() {
	model := lined('equation =')
	names := lined('names =')
	refused(with('format = 4', 'format = 3'), model, 'model.equation', 'IAPWS-95 needs format 4',
		'')
	refused(with('equation = "IAPWS-95"', 'equation = "IAPWS-95"\nalpha = "Twu"'), model + 1,
		'model.alpha', "alpha is a cubic's", '')
	refused(with('equation = "IAPWS-95"', 'equation = "IAPWS-95"\nkij = true'), model + 1,
		'model.kij', "kij is a cubic's", '')
	refused(with('equation = "IAPWS-95"', 'equation = "IAPWS95"'), model, 'model.equation',
		'IAPWS95', 'IAPWS-95')
	refused(with('names = ["water"]', 'names = ["methanol"]').replace('water = 1', 'methanol = 1'),
		names, 'components.names', "IAPWS-95 is water's alone", '')
	refused(with('names = ["water"]', 'names = ["water", "methanol"]').replace('water = 1',
		'water = 0.5, methanol = 0.5'), names, 'components.names', "IAPWS-95 is water's alone",
		'')
	refused(with('names = ["water"]', 'names = ["watr"]').replace('water = 1', 'watr = 1'),
		names, 'components.names', 'no compound watr', 'Water (7732-18-5)')
}

// DWSIM 9's same cycle (testdata/m5_dwsim.json, tools/gen_m5_dwsim.py), by
// state and per kg, its CoolProp package against reflux's solution and its
// Steam Tables package against IF97 worked through the same unit models
// (the pump's v₁ΔP/η, the turbine's η times its isentropic drop). DWSIM's
// water is 18.01528 g/mol, ChemSep's, so its 50 mol/s are 0.900764 kg/s,
// 6.7e-7 above reflux's by IAPWS's molar mass: each power is compared per
// kg of each one's own flow.
//
// Measured, and why:
// - CoolProp: the pump within 2.2e-12 kJ/kg, the turbine 3.2e-11, the
//   condenser 1.0e-10, the exhaust's T 1.8e-12 K and x 5.6e-16; the boiler
//   3.7e-7 kJ/kg low, because DWSIM's PH flash at the pump's outlet stops
//   8.9e-8 K high and the stream carries h at that T (cp·8.9e-8 K), which
//   the boiler then has less to add; the condenser's outlet 2.8e-8 K low.
// - Steam Tables: the pump exact; the boiler 1.2e-6 kJ/kg high, its pump
//   outlet's flash 3.0e-7 K low; the turbine and condenser 8.3e-8 kJ/kg,
//   its exhaust's x 3.5e-11 high (times h_vap, 2392 kJ/kg); the condenser's
//   outlet 1.9e-5 K high, its flash's tolerance there.
fn test_the_rankine_cycle_against_dwsim() {
	l := load(rankine_text())!
	sol := l.solve()!
	flow := 50.0 * 0.018015268 // kg/s
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5_dwsim.json'))!
	d := (json2.decode[json2.Any](text)!).as_map()
	for pkg in ['coolprop', 'steam_tables'] {
		rk := object(object(d, pkg), 'rankine')
		units := object(rk, 'units')
		streams := object(rk, 'streams')
		mf := number(object(streams, 'condensate'), 'mass_flow')
		assert mf == 0.900764 // DWSIM's own 50 mol/s
		mut want := map[string]f64{}
		mut tf, mut tx, mut xx := 0.0, 0.0, 0.0
		if pkg == 'coolprop' {
			for unit in ['pump', 'boiler', 'turbine', 'condenser'] {
				want[unit] = math.abs(power(sol, unit)) / flow * 1e-3
			}
			f := sol.streams['feedwater'] or { panic('no feedwater') }
			e := sol.streams['exhaust'] or { panic('no exhaust') }
			tf, tx, xx = f.t, e.t, e.beta
		} else {
			// IF97 through the same models, per kg on IAPWS's reference.
			s1 := if97.state_pt(1e4, 318.15)!
			w := s1.v * (8e6 - 1e4) / 0.75
			s3 := if97.state_pt(8e6, 773.15)!
			s4s := if97.state_ps(1e4, s3.s)!
			h4 := s3.h - 0.85 * (s3.h - s4s.h)
			want['pump'] = w * 1e-3
			want['boiler'] = (s3.h - s1.h - w) * 1e-3
			want['turbine'] = (s3.h - h4) * 1e-3
			want['condenser'] = (h4 - s1.h) * 1e-3
			e := if97.state_ph(1e4, h4)!
			tf, tx, xx = if97.state_ph(8e6, s1.h + w)!.t, e.t, e.x
		}
		// kJ/kg, by unit; then K, K, x and K for the streams
		bounds := if pkg == 'coolprop' {
			[1e-10, 1e-6, 1e-9, 1e-9, 1e-6, 1e-10, 1e-12, 1e-6]
		} else {
			[1e-10, 5e-6, 1e-6, 1e-6, 2e-6, 1e-10, 1e-9, 1e-4]
		}
		for k, unit in ['pump', 'boiler', 'turbine', 'condenser'] {
			got := number(units, unit) / mf // DeltaQ, kW, per kg/s
			assert math.abs(got - want[unit]) <= bounds[k], '${pkg} ${unit}: ${got} kJ/kg against ${want[unit]}'
		}
		feed, ex := object(streams, 'feedwater'), object(streams, 'exhaust')
		ret := object(streams, 'return')
		assert math.abs(number(feed, 'T') - tf) <= bounds[4], '${pkg} feedwater: ${number(feed, 'T')} K'
		assert math.abs(number(ex, 'T') - tx) <= bounds[5], '${pkg} exhaust: ${number(ex, 'T')} K'
		assert math.abs(number(ex, 'beta') - xx) <= bounds[6], '${pkg} exhaust: x ${number(ex, 'beta')}'
		assert math.abs(number(ret, 'T') - 318.15) <= bounds[7], '${pkg} return: ${number(ret, 'T')} K'
	}
}
