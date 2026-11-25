// Sensitivities (M8 step 8, D10, criterion 6): each derivative against its
// closed form, within its own error's estimate, with the estimate's ratio to
// the true error recorded so that a loose one shows.
module study

import flowsheet
import json2
import math
import os
import reactors
import thermo

fn load_text(text string) !flowsheet.Loaded {
	return flowsheet.load_text(text, makers: reactors.makers())
}

const pfr_file = 'format = 3

[components]
names = ["n-butane", "isobutane"]

[model]
equation = "PR"

[feed.f]
T = "400 K"
P = "1 bar"
flow = "1 mol/s"
mole_fractions = { n-butane = 1, isobutane = 0 }

[reaction.iso]
equation = "n-butane = isobutane"
law = "power_law"
k = { A = "70 mol/(m^3*s*Pa)", Ea = "50 kJ/mol" }
orders = { n-butane = 1 }
K = "none"
per = "volume"

[unit.tube]
type = "pfr"
in = "f"
out = "product"
volume = "1 m^3"
rtol = 1e-12

[study.grid]
vary = [{ key = "feed.f.T", values = ["380 K", "400 K", "420 K"] }, { key = "unit.tube.volume", values = ["0.5 m^3", "1 m^3", "2 m^3"] }]
sensitivities = { outputs = ["X"], inputs = ["unit.tube.volume", "feed.f.T"] }
' +
	'report = { X = "(f.n[\'n-butane\'] - product.n[\'n-butane\']) / f.n[\'n-butane\']" }\n'

// honest is a derivative's check: within its estimate of its error, the
// estimate itself within 1e-5 of the derivative's scale, max(|f′|, |f/v|),
// so that one grown loose shows, and the estimate's ratio to the true error.
// The scale isn't f′ alone: where f′ is tiny beside f/v, as a conversion's
// near 1, the solves' noise over the step is large beside f′ itself.
fn honest(d Derivative, exact f64, f f64, v f64, at string) f64 {
	wrong := math.abs(d.value - exact)
	assert d.reason == '' && math.is_finite(d.value), '${at}: ${d.reason}'
	assert wrong <= d.error, '${at}: off by ${wrong}, estimated ${d.error}'
	scale := math.max(math.abs(exact), math.abs(f / v))
	assert d.error <= 1e-5 * scale, '${at}: estimated ${d.error}'
	return d.error / math.max(wrong, 1e-300)
}

// An isothermal first-order PFR (criterion 6): X = 1 − exp(−kPV/F), its
// derivatives by V and by T against mpmath's closed forms at nine points.
// The PFR integrates to 1e-12, so its X is good to about that, and each
// derivative to about 1e-8 of its value, which is what D10 asks.
fn test_a_pfrs_derivatives_against_their_closed_forms() {
	fixture := json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata',
		'm8_sensitivities.json'))!)!.as_map()
	cases := fixture['cases']!.as_array()
	l := load_text(pfr_file)!
	r := run(l, l.studies[0])
	assert r.points.len == 9
	mut least, mut most := math.inf(1), 0.0
	for k, p in r.points {
		c := cases[k].as_map()
		at := 'T ${c['T']!.f64()} K, V ${c['V']!.f64()} m³'
		assert p.inputs == [c['T']!.f64(), c['V']!.f64()]
		assert math.abs(p.outputs[0] - c['X']!.f64()) <= 1e-10, '${at}: X ${p.outputs[0]}'
		assert p.derivatives.len == 2
		for i, name in ['dX/dV', 'dX/dT'] {
			exact := c[name]!.f64()
			ratio := honest(p.derivatives[i], exact, p.outputs[0], p.inputs[1 - i], '${at} ${name}')
			least = math.min(least, ratio)
			most = math.max(most, ratio)
			// Good to 1e-6 of itself, which an optimiser's gradient needs.
			assert math.abs(p.derivatives[i].value - exact) <= 1e-6 * math.abs(exact), '${at} ${name}'
		}
	}
	println('the PFR: each estimate is ${least:.3g} to ${most:.3g} times its true error')
	// Measured: 273 at least under tcc, 85 under gcc, whose rounding in the
	// PFR's steps differs. The extrapolation removes the truncation the
	// estimate's first part measures, so a derivative that skipped it would
	// sit at about a third of its estimate.
	assert least >= 20.0
}

// cp is a stream's heat capacity at its T and P, W/K for its flows, by
// thermo's analytic second derivatives: no difference in it.
fn cp(fl thermo.Fluid, st flowsheet.Stream) f64 {
	eos := fl.helmholtz() or { panic(err) }
	state := thermo.solve_tp(eos, st.t, st.p, st.flows, .vapor) or { panic(err) }
	return fl.props(state).cp
}

const heater_file = 'format = 3

[components]
names = ["methane", "propane", "n-pentane"]

[model]
equation = "PR"

[feed.f]
T = "350 K"
P = "10 bar"
flow = "100 mol/s"
mole_fractions = { methane = 0.6, propane = 0.3, n-pentane = 0.1 }

[unit.heat]
type = "heater"
in = "f"
out = "out"
outlet_T = "400 K"

[study.one]
vary = [{ key = "unit.heat.outlet_T", values = ["400 K", "450 K"] }]
report = { duty = "heat.duty" }
sensitivities = { outputs = ["duty"], inputs = ["unit.heat.outlet_T", "feed.f.T"] }
'

// A heater's duty, F·(h(T_out) − h(T_in)), moves with its outlet's T by
// F·Cp there and with its feed's by −F·Cp at the feed: each against thermo's
// analytic Cp, all of it one phase of vapour.
fn test_a_heaters_derivatives_against_cp() {
	l := load_text(heater_file)!
	r := run(l, l.studies[0])
	feed := l.sheet.feeds[0].stream
	assert feed.beta == 1.0
	mut least, mut most := math.inf(1), 0.0
	for p in r.points {
		out := flowsheet.feed_tp(l.fluid, feed.flows, p.inputs[0], feed.p)!
		assert out.beta == 1.0
		for i, exact in [cp(l.fluid, out), -cp(l.fluid, feed)] {
			ratio := honest(p.derivatives[i], exact, p.outputs[0], [p.inputs[0], feed.t][i],
				'T_out ${p.inputs[0]}, input ${i}')
			least = math.min(least, ratio)
			most = math.max(most, ratio)
			assert math.abs(p.derivatives[i].value - exact) <= 1e-6 * math.abs(exact)
		}
	}
	println('the heater: each estimate is ${least:.3g} to ${most:.3g} times its true error')
	// Measured: 45 until 0.6.0, and 10.6 once R moved by 1.84e-11. The true
	// error here is the solves' rounding, 3e-12 to 8e-12 of the derivative,
	// which any change in the last digits moves; the estimate, 5e-11 to
	// 5e-10 of it, holds it with this margin.
	assert least >= 5.0
}

// The same feed's derivative with a spec meeting the outlet's T: each
// moved feed meets its spec again at 1e-12, so the duty's derivative is
// −F·Cp at the feed still.
fn test_a_derivative_with_a_spec_met() {
	text := heater_file.replace('outlet_T = "400 K"', 'duty = "0 W"').all_before('[study.one]') +
		'[spec.hot]\nvary = "unit.heat.duty"\nbetween = ["-1 MW", "2 MW"]\nuntil = "out.T = 400 K"\n\n[study.feed]\nvary = [{ key = "feed.f.T", values = ["350 K", "360 K"] }]\nreport = { duty = "heat.duty" }\nsensitivities = { outputs = ["duty"], inputs = ["feed.f.T"] }\n'
	l := load_text(text)!
	r := run(l, l.studies[0])
	for p in r.points {
		feed := flowsheet.feed_tp(l.fluid, l.sheet.feeds[0].stream.flows, p.inputs[0],
			10e5)!
		ratio := honest(p.derivatives[0], -cp(l.fluid, feed), p.outputs[0], p.inputs[0],
			'T_in ${p.inputs[0]}')
		println('the spec: the estimate is ${ratio:.3g} times its true error')
		// Measured: 3.7 at either point, the solves' noise dominating.
		assert ratio >= 2.0
	}
}

// A derivative with no value says why: an input at 0 gives its step no
// scale.
fn test_a_derivative_with_no_value() {
	text := heater_file.replace('outlet_T = "400 K"', 'duty = "0 W"').replace('vary = [{ key = "unit.heat.outlet_T", values = ["400 K", "450 K"] }]',
		'vary = [{ key = "feed.f.T", values = ["350 K"] }]').replace('inputs = ["unit.heat.outlet_T", "feed.f.T"]',
		'inputs = ["unit.heat.duty"]')
	l := load_text(text)!
	r := run(l, l.studies[0])
	d := r.points[0].derivatives[0]
	assert math.is_nan(d.value) && math.is_nan(d.error)
	assert d.reason == 'unit.heat.duty is 0 here, which gives its step no scale', d.reason
}

fn refused_sensitivity(extra string, want string) {
	text := heater_file.all_before('sensitivities = ') + extra
	load_text(text) or {
		assert err.msg().contains(want), err.msg()
		return
	}
	assert false, 'it loaded: ${want}'
}

fn test_sensitivities_refused() {
	refused_sensitivity('sensitivities = { outputs = ["dut"], inputs = ["feed.f.T"] }\n', 'no output dut in report')
	refused_sensitivity('sensitivities = { outputs = ["duty", "duty"], inputs = ["feed.f.T"] }\n',
		'duty is listed twice')
	refused_sensitivity('sensitivities = { outputs = ["duty"], inputs = ["feed.f.T", "feed.f.T"] }\n',
		'feed.f.T is listed twice')
	refused_sensitivity('sensitivities = { outputs = ["duty"], inputs = ["feed.f.Q"] }\n', 'feed.f has no Q in the file')
	refused_sensitivity('sensitivities = { outputs = ["duty"], inputs = [] }\n', 'sensitivities need outputs and inputs')
	refused_sensitivity('sensitivities = { outputs = ["duty"] }\n', 'inputs is missing')
	refused_sensitivity('sensitivities = { outputs = ["duty"], inputs = ["feed.f.T"], step = 1 }\n',
		'no key step here')
	text := heater_file.replace('outlet_T = "400 K"', 'duty = "0 W"').all_before('[study.one]') +
		'[spec.hot]\nvary = "unit.heat.duty"\nbetween = ["-1 MW", "2 MW"]\nuntil = "out.T = 400 K"\n\n[study.feed]\nvary = [{ key = "feed.f.T", values = ["350 K"] }]\nreport = { duty = "heat.duty" }\nsensitivities = { outputs = ["duty"], inputs = ["unit.heat.duty"] }\n'
	load_text(text) or {
		assert err.msg().contains("spec hot varies unit.heat.duty, to meet its target, so it can't be a sensitivity's input"), err.msg()
		return
	}
	assert false
}

// cold_gas is the reference flowsheet's gas at a cooler T, solved cold at
// 1e-13, which no warm start's bias reaches.
fn cold_gas(l flowsheet.Loaded, input flowsheet.Input, t f64) f64 {
	at := l.with([input], [t]) or { panic(err) }
	tight := flowsheet.Loaded{
		...at
		params: flowsheet.SolveParams{
			...at.params
			tol: 1e-13
		}
	}
	sol := tight.solve() or { panic(err) }
	return sol.streams['gas'].total()
}

// A derivative through a recycle: the reference flowsheet's gas against its
// cooler's T, by direct substitution, against Richardson's over cold solves
// at 1e-13 and a step ten times as large. Each warm solve stops where its
// tears move by the tolerance, which leaves an error in proportion to how
// far it started, so to the step: a bias of about tol/1e-4 of the
// derivative, which neither part of the estimate sees. At 1e-12 that's
// within the estimate; at the file's 1e-9 it wouldn't be.
fn test_a_derivative_through_a_recycle() {
	text := os.read_file(os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'reference.toml'))!.replace('format = 1',
		'format = 3').replace('method = "wegstein"', 'method = "direct"') +
		'\n[study.s]\nvary = [{ key = "unit.cooler.outlet_T", values = ["230 K"] }]\nreport = { gas = "gas.F" }\nsensitivities = { outputs = ["gas"], inputs = ["unit.cooler.outlet_T"] }\n'
	l := load_text(text)!
	d := run(l, l.studies[0]).points[0].derivatives[0]
	assert d.reason == '' && d.value > 0
	input := l.input('unit.cooler.outlet_T')!
	h := 1e-3 * 230.0
	wide := (cold_gas(l, input, 230.0 + h) - cold_gas(l, input, 230.0 - h)) / (2 * h)
	narrow := (cold_gas(l, input, 230.0 + h / 2) - cold_gas(l, input, 230.0 - h / 2)) / h
	reference := narrow + (narrow - wide) / 3
	println('through the recycle: ${d.value} against ${reference}, estimated within ${d.error}')
	assert math.abs(d.value - reference) <= d.error, '${d.value} ${reference} ${d.error}'
	// Measured: estimated within 3.6e-8, which can't grow.
	assert d.error <= 5e-8, '${d.error}'
}

// One-sided differences (D13), as the optimiser takes them where a bound is
// within the step: the PFR's dX/dV and dX/dT at T = 400 K, V = 1 m³,
// against mpmath's closed forms, from above and from below, each within its
// estimate and the estimate within 1e-5 of the derivative's scale.
fn test_one_sided_differences() {
	fixture := json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata',
		'm8_sensitivities.json'))!)!.as_map()
	c := fixture['cases']!.as_array()[4].as_map()
	assert c['T']!.f64() == 400.0 && c['V']!.f64() == 1.0
	base := load_text(pfr_file)!
	l := flowsheet.Loaded{
		...base
		params: flowsheet.SolveParams{
			...base.params
			tol: sensitivity_tol
		}
	}
	sol := l.solve()!
	x := "(f.n['n-butane'] - product.n['n-butane']) / f.n['n-butane']"
	e := l.expression(x)!
	f := l.value(e, sol)!
	for key, name in {
		'unit.tube.volume': 'dX/dV'
		'feed.f.T':         'dX/dT'
	} {
		input := l.input(key)!
		exact := c[name]!.f64()
		for side in [Side.above, .below] {
			d := differences(l, sol, [e], ['X'], input, 1e-4 * input.value, side)
			ratio := honest(d[0], exact, f, input.value, '${name} ${side}')
			println('${name} from ${side}: the estimate is ${ratio:.3g} times its true error')
			assert ratio >= 2.0
		}
	}
}
