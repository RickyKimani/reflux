// M8 against DWSIM 9 (step 13, D18, criteria 8 and 11): a spec against its
// Adjust, an optimisation against its optimiser, and criterion 8's series
// reactions against its PFR and its optimiser, each difference explained.
// tools/gen_m8_dwsim.py writes DWSIM's answers.
module study

import flowsheet
import json2
import math
import os
import reactors

fn dwsim() map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm8_dwsim.json')) or {
		panic(err)
	}) or { panic(err) }.as_map()
}

fn flowsheet_file(name string, extra string) !flowsheet.Loaded {
	path := os.join_path(@VMODROOT, 'testdata', 'flowsheets', name)
	text := os.read_file(path)!.replace('format = 1', 'format = 3') + extra
	return flowsheet.load_text(text, makers: reactors.makers(), dir: os.dir(path))
}

// tightened is a file solved to 1e-12, so that its own error is far below what's
// compared.
fn tightened(l flowsheet.Loaded) flowsheet.Loaded {
	return flowsheet.Loaded{
		...l
		params: flowsheet.SolveParams{
			...l.params
			tol: 1e-12
		}
	}
}

// The reference flowsheet's spec (D18): the cooler's outlet T between 200
// and 260 K until the first drum's vapour is 60 mol/s, on DWSIM's constants
// and kij. DWSIM's answer, by its Adjust's Brent at 1e-10 (7 solves) and by
// SciPy's brentq on its solves, is 219.7546705409 K. reflux's is 2.8 mK
// above it, which the two flowsheets' states explain: at DWSIM's answer
// reflux's vapour is 1.6e-3 mol/s short of 60, the difference M4 found
// between them (its pump's liquid density, through the recycle), and over
// the vapour's slope there that's the 2.8 mK, to 1e-6 K.
fn test_the_spec_against_dwsims_adjust() {
	d := dwsim()['spec']!.as_map()
	t_d, slope := d['T']!.f64(), d['slope']!.f64()
	plain := tightened(flowsheet_file('reference_dwsim.toml', '')!)
	gas := plain.expression('gas.F')!
	at := fn [plain, gas] (t f64) f64 {
		l := plain.with([plain.input('unit.cooler.outlet_T') or { panic(err) }], [t]) or { panic(err) }
		return l.value(gas, l.solve() or { panic(err) }) or { panic(err) }
	}
	// By state: the vapour at 230 K, within M4's 2.2e-5.
	v230 := at(230.0)
	assert math.abs(v230 - d['vapour_at_230']!.f64()) <= 2.5e-5 * v230, '${v230}'
	specced := tightened(flowsheet_file('reference_dwsim.toml', '\n[spec.vapour]\nvary = "unit.cooler.outlet_T"\nbetween = ["200 K", "260 K"]\nuntil = "gas.F = 60 mol/s"\n')!)
	t_r := specced.solve()!.specs[0].value
	short := 60.0 - at(t_d)
	predicted := short / slope
	println("the spec: reflux ${t_r} K, DWSIM ${t_d} K, ${t_r - t_d} K apart; ${short} mol/s short at DWSIM's, which predicts ${predicted} K")
	assert math.abs((t_r - t_d) - predicted) <= 1e-6, '${t_r - t_d} against ${predicted}'
}

// The two-stage compression's least work (D18): the middle pressure
// between 1.5 and 12 bar, on DWSIM's constants and kij. DWSIM's optimum by
// SciPy on its solves, and by its own optimiser with bounded Brent (18
// evaluations), is 4.05502501 bar and 106171.3175 W.
//
// The two works differ by 1.5e-5 to 3.2e-5, each stage's by 1.4e-5 to
// 5.4e-5 and its outlet by 0.1 to 4.6 mK, scattered with the pressure: M3
// found DWSIM's ideal-gas H up to 1e-5 off an exact quadrature of the same
// Cp, and its S alike, and a work of about 5 kJ/mol between two such
// states carries both. Through the two works' slopes that moves the
// optimum: reflux's own slope at DWSIM's optimum, over the work's
// curvature, predicts reflux's optimum 1.5 Pa above DWSIM's, 3.6e-6 of the
// pressure, to 0.3 Pa, the two optimisers' tolerances.
fn test_the_compression_against_dwsims_optimiser() {
	d := dwsim()['compression']!.as_map()
	l := tightened(flowsheet_file('compression_dwsim.toml', '')!)
	work := l.expression('c1.work + c2.work')!
	at := fn [l, work] (p f64) f64 {
		x := l.with([l.input('unit.c1.outlet_P') or { panic(err) }], [p]) or { panic(err) }
		return x.value(work, x.solve() or { panic(err) }) or { panic(err) }
	}
	for key, w_d in d['work_at']!.as_map() {
		p := key.f64()
		w := at(p)
		assert math.abs(w - w_d.f64()) <= 3.5e-5 * w_d.f64(), '${p} Pa: ${w} W, DWSIM ${w_d.f64()}'
	}
	inside := optimize(l, l.optimizations[0])
	assert inside.ok, inside.error
	p_d, w_d := d['P']!.f64(), d['work']!.f64()
	slope := (at(p_d + 100.0) - at(p_d - 100.0)) / 200.0
	predicted := -slope / d['curvature']!.f64()
	println("the compression: reflux ${inside.inputs[0]} Pa, ${inside.objective} W; DWSIM ${p_d} Pa, ${w_d} W; reflux's slope at DWSIM's optimum predicts ${predicted} Pa")
	assert math.abs((inside.inputs[0] - p_d) - predicted) <= 0.3, '${inside.inputs[0] - p_d} Pa against ${predicted}'
	assert math.abs(inside.objective - w_d) <= 2e-5 * w_d
	// On a bound, 6 bar: reflux there exactly, DWSIM's bounded Brent
	// within rounding of it in 55 evaluations; the bound's multiplier is
	// the work's slope there, DWSIM's 0.0112548 W/Pa by differences, within
	// 3e-4, the slopes' difference at 6 bar.
	bound := optimize(l, l.optimizations[1])
	assert bound.ok && bound.inputs[0] == 6e5
	forms := d['forms']!.as_map()
	assert math.abs(forms['bound']!.as_map()['P']!.f64() - 6e5) <= 1e-3
	at_6 := d['slope_at_6_bar']!.f64()
	assert math.abs(bound.lower[0] - at_6) <= 3e-4 * at_6, '${bound.lower[0]} against ${at_6}'
}

// Criterion 8's series reactions in DWSIM's PFR, typed per kg on 1 kg/m³,
// against the closed form's V* = 0.654998 m³, where reflux's optimum is
// within 1e-5 (optimize_test.v).
//
// DWSIM's isopentane is the closed form's to 5e-11 at most volumes, and off
// by up to 1.3e-4 at others: near the peak, 0.64 m³ is 1.3e-4 low and SciPy's
// optimum 4.5e-5 high. Its PFR integrates each section at a relative
// tolerance of 1e-3 that nothing sets (M9 step 0), and where that
// tolerance lets a step through, the outflow jumps. On a peak this flat,
// B″ = −1.11 mol/(s·m⁶), an error e can move an apparent maximum by up to
// √(2e/|B″|), 0.015 m³ for 1.3e-4: SciPy's search on DWSIM's solves ends
// 0.0095 m³ off (1.5%), DWSIM's own bounded Brent 0.0029 m³ (0.44%) in 22
// evaluations.
fn test_the_series_reactions_against_dwsim() {
	d := dwsim()['series']!.as_map()
	s := series(400.0)
	mut worst := 0.0
	mut exact := 0
	mut volumes := d['outflow_at']!.as_map().clone()
	for key, b in d['near_peak']!.as_map() {
		volumes[key] = b
	}
	for key, b in volumes {
		v := key.f64()
		e := math.abs(b.f64() - s.b(v))
		worst = math.max(worst, e)
		if e <= 1e-10 {
			exact++
		}
	}
	// Measured: 1.32e-4 at worst, and 7 of 12 volumes within 1e-10.
	assert worst <= 1.4e-4 && worst > 1e-5, '${worst}'
	assert exact >= 7 && volumes.len == 12, '${exact}'
	v_star := math.log(s.k2 / s.k1) / (s.k2 - s.k1)
	reach := math.sqrt(2.0 * worst / math.abs(s.b2(v_star)))
	form := d['form']!.as_map()
	println('the series: DWSIM ${d['V']!.f64()} m³ by SciPy, ${form['V']!.f64()} by bounded Brent in ${form['evaluations']!.int()} evaluations; exact ${v_star}; its noise reaches ${reach} m³')
	for v in [d['V']!.f64(), form['V']!.f64()] {
		assert math.abs(v - v_star) <= reach && math.abs(v - v_star) > 1e-3, '${v}'
	}
}

// series is the closed form of optimize_test.v, at T.
struct SeriesAt {
	k1 f64
	k2 f64
}

fn series(t f64) SeriesAt {
	r := 1.380649e-23 * 6.02214076e23
	return SeriesAt{
		k1: 70.0 * math.exp(-50000.0 / (r * t)) * 1e5
		k2: 0.09 * math.exp(-30000.0 / (r * t)) * 1e5
	}
}

fn (s SeriesAt) b(v f64) f64 {
	return s.k1 / (s.k2 - s.k1) * (math.exp(-s.k1 * v) - math.exp(-s.k2 * v))
}

// b2 is B's second derivative by V.
fn (s SeriesAt) b2(v f64) f64 {
	return s.k1 / (s.k2 - s.k1) * (s.k1 * s.k1 * math.exp(-s.k1 * v) - s.k2 * s.k2 * math.exp(-s.k2 * v))
}
