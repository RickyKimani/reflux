// M5b step 2, criteria 2 and 3: a gas's transport (D4, D5) against
// chemicals and mpmath (testdata/m5b_transport.json, from
// tools/gen_m5b_transport.py): Mason and Saxena's conductivity on Wilke's
// factors, Lucas's high-pressure factor and its polar and quantum terms,
// Stiel and Thodos' residual on Yorizane's pseudo-critical constants; each
// correction 1 or 0 at zero pressure, and each refusal.
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

fn pures(db components.Database, cas []json2.Any) ![]transport.Pure {
	mut out := []transport.Pure{}
	for c in cas {
		out << db.find(c.str())!.transport_pure()!
	}
	return out
}

// Mason and Saxena's conductivity, Wassiljewa's form on Wilke's Φᵢⱼ of the
// viscosities, is mpmath's at 50 digits on chemicals' pure values, within
// 1e-14: the pure values are tdep's and chemicals' to a few ulps, and the
// sums of positive terms add nothing (measured: 2.2e-16). One compound
// alone is its own conductivity, to the bit.
fn test_mason_saxena_against_mpmath() {
	db := components.load_default()!
	f := fixture()!
	mut worst := 0.0
	for c in f['mason_saxena']!.as_array() {
		m := c.as_map()
		cs := pures(db, m['cas']!.as_array())!
		y := floats(m['y']!)
		t := m['t']!.f64()
		mut corr := []transport.Correlated{}
		mut lambda := []f64{}
		for p in cs {
			corr << p.sources[int(transport.Quantity.vapour_viscosity)].c
			lambda << p.value(.vapour_conductivity, t)!.value
		}
		w := transport.new_wilke(corr, cs.map(it.mw))!
		mut k := w.room()
		w.at(t, mut k)!
		got := w.conductivity(k, y, lambda)
		want := m['k']!.f64()
		e := math.abs(got / want - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, '${m['name']!.str()}: ${got} vs ${want}'
		if cs.len == 1 {
			assert got == lambda[0]
		}
	}
	println('Mason–Saxena within ${worst} of mpmath')
}

// Lucas's factor is mpmath's at 50 digits from the published equations,
// within 1e-13: powers of Pr to 5.5 and Tr to −37.7 (measured: 4.4e-16).
// A's bounds are compared as doubles, as here: the reforming gas's CO₂ is
// exactly 0.05, on its bound, where mpmath's decimal 0.05 would take it in
// and move the factor by 3e-4.
fn test_lucas_factor_against_mpmath() {
	db := components.load_default()!
	f := fixture()!
	mut worst := 0.0
	mut refused := 0
	for c in f['lucas_factor']!.as_array() {
		m := c.as_map()
		cs := pures(db, m['cas']!.as_array())!
		y := floats(m['y']!)
		t, p := m['t']!.f64(), m['p']!.f64()
		name := m['name']!.str()
		if m['factor']! is json2.Null {
			if got := transport.lucas_factor(cs, y, t, p) {
				assert false, "${name}: ${got} where the form doesn't hold"
			}
			refused++
			continue
		}
		got := transport.lucas_factor(cs, y, t, p)!
		want := m['factor']!.f64()
		e := math.abs(got / want - 1)
		worst = math.max(worst, e)
		assert e <= 1e-13, '${name}: ${got} vs ${want}'
		// At zero pressure, 1 to the bit (D5): Y is 1 exactly, and Fp and
		// Fq are written to be 1 at it.
		assert transport.lucas_factor(cs, y, t, 0)! == 1.0, name
	}
	assert refused == 0
	println('Lucas factors within ${worst} of mpmath')
}

// Lucas's own low-pressure viscosity, η° = Z₁/ξ, is chemicals' Lucas_gas,
// which holds Fp° and Fq°, within 1e-14: polar compounds (water, methanol,
// ammonia) on each branch of Fp°, hydrogen and helium with Fq° on each side
// of Tr = 12 (measured: 4.4e-16).
fn test_lucas_polar_and_quantum_terms_against_chemicals() {
	db := components.load_default()!
	f := fixture()!
	mut worst := 0.0
	for c in f['lucas_dilute']!.as_array() {
		r := c.as_array()
		p := db.find(r[0].str())!.transport_pure()!
		t := r[1].f64()
		got := transport.lucas_dilute(p, t)!
		want := r[2].f64()
		e := math.abs(got / want - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, '${p.name} at ${t} K: ${got} vs ${want}'
	}
	println("Lucas's dilute viscosities within ${worst} of chemicals")
}

// Yorizane's pseudo-critical constants are mpmath's, and Stiel and Thodos'
// residual chemicals' on them, within 1e-14: cube roots and sums of
// positive terms, then a power and an exponential (measured: 1.3e-15).
// Past ρr = 2.8 it's refused where the fixture has none.
fn test_stiel_thodos_against_chemicals() {
	db := components.load_default()!
	f := fixture()!
	mut worst := 0.0
	mut refused := 0
	for c in f['stiel_thodos']!.as_array() {
		m := c.as_map()
		cs := pures(db, m['cas']!.as_array())!
		y := floats(m['y']!)
		v := m['v']!.f64()
		name := m['name']!.str()
		ps := transport.yorizane(cs, y)!
		want := m['pseudo']!.as_map()
		gots := [ps.tc, ps.pc, ps.vc, ps.zc, ps.mw]
		for i, key in ['tc', 'pc', 'vc', 'zc', 'mw'] {
			got := gots[i]
			e := math.abs(got / want[key]!.f64() - 1)
			worst = math.max(worst, e)
			assert e <= 1e-14, '${name}: ${key} ${got} vs ${want[key]!.f64()}'
		}
		if m['residual']! is json2.Null {
			if got := transport.stiel_thodos(ps, v) {
				assert false, '${name}: ${got} past ρr = 2.8'
			}
			refused++
			continue
		}
		got := transport.stiel_thodos(ps, v)!
		r := m['residual']!.f64()
		e := math.abs(got / r - 1)
		worst = math.max(worst, e)
		assert e <= 1e-14, '${name}: ${got} vs ${r}'
		// At zero density, none.
		assert transport.stiel_thodos(ps, math.inf(1))! == 0.0
	}
	println('Stiel–Thodos within ${worst} of chemicals; ${refused} refused')
}

// What Lucas's factor and Stiel and Thodos' residual refuse, each saying
// why: a gas below its pseudo-critical temperature above its pseudo-critical
// pressure; Tr past 40 or Pr past 100; a polar compound's Fp° without its
// dipole moment, or with Zc at 0.292 or more; ρr past 2.8.
fn test_refusals() {
	db := components.load_default()!
	propane := [db.find('propane')!.transport_pure()!]
	transport.lucas_factor(propane, [1.0], 300, 5e6) or {
		assert err.msg().contains('is for a vapour below its saturation pressure'), err.msg()
	}
	if f := transport.lucas_factor(propane, [1.0], 300, 5e6) {
		assert false, 'a factor of ${f}'
	}
	h2 := [db.find('hydrogen')!.transport_pure()!]
	if f := transport.lucas_factor(h2, [1.0], 40 * h2[0].tc, 1e5) {
		assert false, 'a factor of ${f} at Tr 40'
	} else {
		assert err.msg().contains('is for Tr below 40 and Pr to 100'), err.msg()
	}
	if f := transport.lucas_factor(h2, [1.0], 300, 101 * h2[0].pc) {
		assert false, 'a factor of ${f} at Pr 101'
	}
	nodipole := db.find('100-18-5')!.transport_pure()!
	if f := transport.lucas_factor([nodipole], [1.0], 800, 1e5) {
		assert false, 'a factor of ${f}'
	} else {
		assert err.msg().contains("needs P-diisopropylbenzene's dipole moment, which its component file hasn't"), err.msg()
	}
	w := db.find('water')!.transport_pure()!
	odd := transport.Pure{
		...w
		name: 'a polar compound'
		zc:   0.3
	}
	if f := transport.lucas_factor([odd], [1.0], 800, 1e5) {
		assert false, 'a factor of ${f}'
	} else {
		assert err.msg().contains("needs Zc below 0.292, and a polar compound's is 0.3"), err.msg()
	}
	ps := transport.yorizane([w], [1.0])!
	if r := transport.stiel_thodos(ps, ps.vc / 2.9) {
		assert false, 'a residual of ${r}'
	} else {
		assert err.msg().contains("is for ρr up to 2.8, and the gas's is"), err.msg()
	}
	assert transport.stiel_thodos(ps, ps.vc / 2.8)! > 0
}
