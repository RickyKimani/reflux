// The gas viscosity (M9 step 9, D10 as amended in step 0): each compound's
// from its component file, mixed by Wilke's rule, against chemicals'
// (testdata/m9_viscosity.json).
module reactors

import components
import json2
import math
import os

struct ViscosityCase {
	name    string
	cas     []string
	y       []f64
	t       f64
	pure    []f64
	mixture f64
}

struct ViscosityFixture {
	cases []ViscosityCase
}

// Every pure gas and mixture is chemicals' to 1e-13: equation 102's power
// and Wilke's square roots, in values whose rounding is a few parts in
// 1e16 (measured: 2.2e-16, one unit of it).
fn test_against_chemicals() {
	db := components.load_default()!
	path := os.join_path(@VMODROOT, 'testdata', 'm9_viscosity.json')
	f := json2.decode[ViscosityFixture](os.read_file(path)!)!
	assert f.cases.len == 37
	mut worst := 0.0
	for c in f.cases {
		records := c.cas.map(db.find(it) or { panic(err) })
		w := new_wilke(records)!
		for i, want in c.pure {
			got := w.pure[i].at(c.t).value
			worst = math.max(worst, math.abs(got / want - 1))
			assert math.abs(got / want - 1) <= 1e-13, '${c.name}, compound ${i}: ${got} vs ${want}'
		}
		got := w.viscosity(c.t, 1e5, c.y)!
		worst = math.max(worst, math.abs(got / c.mixture - 1))
		assert math.abs(got / c.mixture - 1) <= 1e-13, '${c.name}: ${got} vs ${c.mixture}'
	}
	println('viscosities: within ${worst} of chemicals')
}

fn test_refused() {
	db := components.load_default()!
	w := new_wilke([db.find('methane')!, db.find('hydrogen')!])!
	w.viscosity(500, 1e5, [1.0]) or {
		assert err.msg() == '1 mole fractions for 2 compounds'
		return
	}
	assert false
}
