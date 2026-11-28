// M5b step 3, criteria 2 and 3: a liquid's transport (D6, D7, D8) against
// chemicals and mpmath (testdata/m5b_transport.json, from
// tools/gen_m5b_transport.py): Lucas's and Missenard's compressed-liquid
// factors, the Grunberg–Nissan rule, Li's rule and Winterfeld, Scriven and
// Davis'; each factor 1 at its own limit.
import components
import json2
import math
import os
import transport

fn fixture() !map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5b_transport.json'))!)!.as_map()
}

fn floats(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

// Lucas's factor is chemicals' Lucas on the same vapour pressure, within
// 1e-14 (measured: 0, the same bits), and the vapour pressure tdep's of chemicals'
// within 1e-13 (DIPPR 101's exponent, as pure_test.v); at or below Pˢᵃᵗ it's
// 1 to the bit.
fn test_lucas_liquid_against_chemicals() {
	db := components.load_default()!
	f := fixture()!
	mut worst := 0.0
	for c in f['lucas_liquid']!.as_array() {
		r := c.as_array()
		p := db.find(r[0].str())!.transport_pure()!
		t, pr, psat := r[1].f64(), r[2].f64(), r[3].f64()
		assert math.abs(p.value(.vapour_pressure, t)!.value / psat - 1) <= 1e-13
		got := transport.lucas_liquid(t, pr, p.tc, p.pc, p.omega, psat)
		want := r[4].f64()
		e := math.abs(got / want - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, '${p.name} at ${t} K and ${pr} Pa: ${got} vs ${want}'
		assert transport.lucas_liquid(t, psat, p.tc, p.pc, p.omega, psat) == 1.0
		assert transport.lucas_liquid(t, psat / 2, p.tc, p.pc, p.omega, psat) == 1.0
	}
	println("Lucas's liquid factors within ${worst} of chemicals")
}

// Missenard's factor is chemicals' spline of degree 1 on his table, held
// at its edges, within 1e-14 (measured: 0, the same bits), at Tr and Pr inside the
// table, on its lines and beyond each edge; at zero pressure, 1.
fn test_missenard_against_chemicals() {
	f := fixture()!
	mut worst := 0.0
	for c in f['missenard']!.as_array() {
		r := floats(c)
		got := transport.missenard(r[0], r[1], r[2], r[3])
		e := math.abs(got / r[4] - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, 'Tr ${r[0] / r[2]}, Pr ${r[1] / r[3]}: ${got} vs ${r[4]}'
		if r[1] == 0 {
			assert got == 1.0
		}
	}
	println("Missenard's factors within ${worst} of chemicals")
}

// The three mixing rules on given numbers, against chemicals' DIPPR9I and
// Winterfeld_Scriven_Davis and mpmath's Grunberg–Nissan at 50 digits, within
// 1e-14: sums of positive terms, and one exponential (measured: 5.6e-16).
fn test_mixing_rules() {
	f := fixture()!
	mut worst := 0.0
	for c in f['li']!.as_array() {
		m := c.as_map()
		got := transport.li(floats(m['x']!), floats(m['v']!), floats(m['k']!))
		e := math.abs(got / m['value']!.f64() - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, 'Li: ${got} vs ${m['value']!.f64()}'
	}
	for c in f['winterfeld']!.as_array() {
		m := c.as_map()
		got := transport.winterfeld(floats(m['x']!), floats(m['v']!), floats(m['sigma']!))
		e := math.abs(got / m['value']!.f64() - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, 'Winterfeld–Scriven–Davis: ${got} vs ${m['value']!.f64()}'
	}
	for c in f['grunberg_nissan']!.as_array() {
		m := c.as_map()
		mut g := []f64{}
		for row in m['g']!.as_array() {
			g << floats(row)
		}
		got := transport.grunberg_nissan(floats(m['x']!), floats(m['mu']!), g)
		e := math.abs(got / m['value']!.f64() - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, 'Grunberg–Nissan: ${got} vs ${m['value']!.f64()}'
	}
	println('mixing rules within ${worst} of chemicals and mpmath')
}
