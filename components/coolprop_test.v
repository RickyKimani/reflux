module components

import json2
import math
import os

// M2 exit criterion 6: ChemSep's data against CoolProp's reference
// equations, for the 85 CoolProp fluids ChemSep also has.
// tools/gen_m2_fixtures.py recorded ChemSep's deviation from CoolProp in
// testdata/m2_coolprop.json, evaluating ChemSep's correlations with its own
// formulas. Here the deviations are computed again through reflux's
// database, unit conversions and tdep. They must match, and they can't have
// grown, so a change to the data or its conversion that makes things worse
// fails. Large recorded deviations are ChemSep's, and step 8's report lists
// them.

struct CoolPropSource {
	generated string
	generator string
	python    string
	coolprop  string
	chemsep   string
}

struct Ref {
	coolprop  f64
	deviation f64
}

struct RefPoint {
	t         f64
	coolprop  f64
	deviation f64
}

// An absent curve in the JSON decodes as an empty one.
struct RefCurve {
	key    string
	points []RefPoint
}

struct RefFluid {
	coolprop  string
	cas       string
	name      string
	constants map[string]Ref
	cp0       RefCurve
	psat      RefCurve
}

struct CoolPropFixture {
	source CoolPropSource
	note   string
	fluids []RefFluid
}

fn coolprop_fixture() !CoolPropFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm2_coolprop.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m2_fixtures.py`.') }
	return json2.decode[CoolPropFixture](text)!
}

// same checks a deviation against the recorded one: equal, to rounding, and
// not larger.
fn same(what string, got f64, recorded f64) {
	tol := 1e-10 * math.max(1.0, math.abs(recorded))
	assert math.abs(got - recorded) <= tol, '${what}: deviation ${got}, recorded ${recorded}'
	assert math.abs(got) <= math.abs(recorded) + tol, '${what}: deviation grew from ${recorded} to ${got}'
}

fn constant(c Constants, key string) ?f64 {
	return match key {
		'tc' { c.tc }
		'pc' { c.pc }
		'omega' { c.omega }
		'mw' { c.mw }
		else { none }
	}
}

fn test_chemsep_against_coolprop() {
	fx := coolprop_fixture()!
	db := load_default()!
	assert fx.fluids.len == 85
	mut constants, mut cp0, mut psat := 0, 0, 0
	for f in fx.fluids {
		r := db.find(f.cas)!
		assert r.name == f.name, f.cas
		for key, ref in f.constants {
			v := constant(r.constants, key) or { panic('${f.cas} has no ${key}') }
			dev := if key == 'omega' { v - ref.coolprop } else { (v - ref.coolprop) / ref.coolprop }
			same('${f.coolprop} ${key}', dev, ref.deviation)
			constants++
		}
		if f.cp0.points.len > 0 {
			ig := r.ideal_gas()!
			for p in f.cp0.points {
				same('${f.coolprop} Cp⁰ at ${p.t} K', (ig.cp.cp(p.t) - p.coolprop) / p.coolprop,
					p.deviation)
				cp0++
			}
		}
		if f.psat.points.len > 0 {
			vp := r.property('vapor_pressure')!
			// reflux falls back to the same correlation the generator chose.
			assert vp.key == f.psat.key, '${f.coolprop}: ${vp.key}, the fixture used ${f.psat.key}'
			for p in f.psat.points {
				same('${f.coolprop} Psat at ${p.t} K', (vp.at(p.t).value - p.coolprop) / p.coolprop,
					p.deviation)
				psat++
			}
		}
	}
	// Every fluid has its four constants; all but R-1234yf, which ChemSep
	// gives no ideal-gas Cp, have a Cp⁰ curve; all have a vapour pressure.
	assert constants == 85 * 4
	assert cp0 == 84 * 4
	assert psat == 85 * 4
}

fn test_typical_agreement() {
	// What the comparison says about ChemSep as a whole: the median
	// deviations, as recorded at generation. A regression in the data would
	// move these; the outliers are ChemSep's, listed by step 8's report.
	fx := coolprop_fixture()!
	mut by := map[string][]f64{}
	for f in fx.fluids {
		for key, ref in f.constants {
			by[key] << math.abs(ref.deviation)
		}
		for p in f.cp0.points {
			by['cp0'] << math.abs(p.deviation)
		}
		for p in f.psat.points {
			by['psat'] << math.abs(p.deviation)
		}
	}
	median := fn [by] (key string) f64 {
		mut v := by[key].clone()
		v.sort()
		return v[v.len / 2]
	}
	assert median('tc') < 1e-3 // 0.03 %
	assert median('pc') < 1e-2 // 0.4 %
	assert median('omega') < 5e-3 // absolute
	assert median('mw') < 1e-6
	assert median('cp0') < 1e-2 // 0.3 %
	assert median('psat') < 2e-2 // 0.8 %
}
