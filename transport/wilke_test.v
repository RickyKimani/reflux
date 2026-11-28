// The gas viscosity at low pressure (M9 step 9, D10 as amended in step 0,
// moved to transport by M5b's D2): each compound's from its component
// file, mixed by Wilke's rule, against chemicals'
// (testdata/m9_viscosity.json).
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
		w := components.wilke(records)!
		for i, want in c.pure {
			got := w.pure[i].at(c.t)
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
	w := components.wilke([db.find('methane')!, db.find('hydrogen')!])!
	mut small := components.wilke([db.find('methane')!])!.room()
	if _ := w.at(500, mut small) {
		assert false, 'room for one viscosity of two'
	} else {
		assert err.msg() == 'room for 1 viscosities for 2 compounds'
	}
	w.viscosity(500, 1e5, [1.0]) or {
		assert err.msg() == '1 mole fractions for 2 compounds'
		return
	}
	assert false
}

// What Wilke's rule keeps at a T and mixes later (M8's D11) gives the whole
// viscosity, to the bit, at every case of the fixture: mixed at the case's
// own composition and at others from the one set of values, as an
// isothermal reactor does along its length. A T it already holds isn't
// computed again, which a value overwritten by hand shows; another T
// replaces it; and a T it refuses leaves it holding nothing.
fn test_kept_and_mixed_is_the_whole() {
	db := components.load_default()!
	path := os.join_path(@VMODROOT, 'testdata', 'm9_viscosity.json')
	f := json2.decode[ViscosityFixture](os.read_file(path)!)!
	mut mixtures := 0
	for c in f.cases {
		records := c.cas.map(db.find(it) or { panic(err) })
		w := components.wilke(records)!
		n := c.y.len
		mut k := w.room()
		assert !k.kept
		w.at(c.t, mut k)!
		assert k.kept && k.t == c.t
		for i, want in c.pure {
			assert math.abs(k.mu[i] / want - 1) <= 1e-13
		}
		mut reversed := c.y.clone()
		reversed.reverse_in_place()
		even := []f64{len: n, init: 1.0 / n}
		mut first := []f64{len: n}
		first[0] = 1.0
		for y in [c.y, reversed, even, first] {
			whole := w.viscosity(c.t, 1e5, y)!
			assert math.f64_bits(w.mixed(k, y)) == math.f64_bits(whole), '${c.name} at ${y}'
			mixtures++
		}
		assert w.mixed(k, first) == k.mu[0]
		// Kept, replaced and dropped.
		right := k.mu[0]
		k.mu[0] = 2 * right
		w.at(c.t, mut k)!
		assert k.mu[0] == 2 * right, 'computed again at an unmoved T'
		w.at(c.t + 1, mut k)!
		assert k.t == c.t + 1 && k.mu[0] != 2 * right && k.mu[0] != right
		assert math.f64_bits(w.mixed(k, c.y)) == math.f64_bits(w.viscosity(c.t + 1, 1e5, c.y)!)
		if _ := w.at(-1, mut k) {
			assert false, '${c.name}: a viscosity at −1 K'
		}
		assert !k.kept
		w.at(c.t, mut k)!
		assert k.kept && k.mu[0] == right
	}
	assert mixtures == 4 * 37
}
