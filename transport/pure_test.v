// M5b step 1, criterion 1: each compound's transport values are its
// component file's correlations, in SI, against chemicals' DIPPR forms on
// every one of them (testdata/m5b_transport.json, from
// tools/gen_m5b_transport.py); and D3's refusals.
import components
import json2
import math
import os
import transport

fn fixture() !map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5b_transport.json'))!)!.as_map()
}

fn quantity_of(key string) transport.Quantity {
	for q in transport.quantities {
		if q.key() == key {
			return q
		}
	}
	panic('no quantity ${key}')
}

// Every point of every transport correlation, 11,058 of them, is
// chemicals', within 1e-15 of the size of its terms, which the fixture
// gives: tdep and chemicals write each form differently, equation 16's
// A + e^(…) and the polynomials cancel near a zero, and an exponential
// form carries its exponent's rounding into the value, DIPPR 101's terms
// reaching 500 for a liquid's viscosity (measured: 2.8e-16 of the terms;
// p-cresol's viscosity at 563.72 K is 1.1e-13 of itself from chemicals').
// A point chemicals finds no positive value at is refused. The one compound without Tc, Pc, ω or M, which makes no Pure,
// has no correlations to check.
fn test_every_correlation_against_chemicals() {
	db := components.load_default()!
	f := fixture()!
	keys := f['keys']!.as_array().map(it.str())
	mut pures := map[string]transport.Pure{}
	mut unmade := map[string]bool{}
	mut worst, mut checked, mut refused := 0.0, 0, 0
	for row in f['pure']!.as_array() {
		r := row.as_array()
		cas := r[0].str()
		q := quantity_of(keys[r[1].int()])
		t := r[2].f64()
		if cas !in pures && cas !in unmade {
			pures[cas] = db.find(cas)!.transport_pure() or {
				assert err.msg().contains('so no equation of state can use it'), err.msg()
				unmade[cas] = true
				continue
			}
		}
		if cas in unmade {
			continue
		}
		p := pures[cas] or { panic('no ${cas}') }
		want := if r[3] is json2.Null { math.nan() } else { r[3].f64() }
		if !(want > 0) {
			if got := p.value(q, t) {
				assert false, '${cas} ${q.key()} at ${t} K: ${got.value} where chemicals has ${want}'
			}
			refused++
			continue
		}
		scale := if r.len > 4 { r[4].f64() } else { math.abs(want) }
		got := p.value(q, t)!
		assert got.in_range, '${cas} ${q.key()} at ${t} K'
		e := math.abs(got.value - want) / scale
		worst = math.max(worst, e)
		assert e <= 1e-15, '${cas} ${q.key()} at ${t} K: ${got.value} vs ${want}'
		checked++
	}
	assert unmade.len == 0, '${unmade.keys()}'
	assert checked + refused > 11000
	println('${checked} points within ${worst} of their terms of chemicals, ${refused} refused with it')
}

// D3: a saturated liquid's quantity is none at or above the compound's Tc,
// whatever its correlation would give (water's liquid density runs to
// 647.29 K, past its Tc of 647.14 K); a value that isn't positive is none,
// saying so (190 surface tensions go negative by Tc, water's among them, by
// 646.14 K); a quantity the file hasn't, or in a form tdep can't evaluate,
// is none with the file's lack; a vapour's quantity runs past Tc.
fn test_refusals() {
	db := components.load_default()!
	w := db.find('water')!.transport_pure()!
	tc := w.tc
	for q in [transport.Quantity.liquid_viscosity, .liquid_conductivity, .surface_tension,
		.liquid_density, .vapour_pressure] {
		if r := w.value(q, tc) {
			assert false, '${q.key()} at Tc: ${r.value}'
		} else {
			assert err.msg().contains("is the saturated liquid's, and 647.14 K isn't below its Tc"), err.msg()
		}
		assert w.value(q, tc - 10)!.value > 0
	}
	assert w.value(.vapour_viscosity, 1000)!.value > 0
	assert w.value(.vapour_conductivity, 1000)!.value > 0
	// Outside its fitted range a value is the correlation's, so marked.
	r := w.value(.vapour_viscosity, 2000)!
	assert !r.in_range && r.value > 0
	if s := w.value(.surface_tension, tc - 1) {
		assert false, 'a surface tension of ${s.value} 1 K below Tc'
	} else {
		assert err.msg().contains("which isn't one"), err.msg()
	}
	// P-diisopropylbenzene's surface tension is negative at 666.21 K, below
	// its Tc (data-quality.md).
	dipb := db.find('100-18-5')!.transport_pure()!
	if s := dipb.value(.surface_tension, 666.21) {
		assert false, 'a surface tension of ${s.value}'
	} else {
		assert err.msg().contains("which isn't one"), err.msg()
	}
	// Squalane's liquid viscosity is equation 118.
	sq := db.find('111-01-3')!.transport_pure()!
	if s := sq.value(.liquid_viscosity, 400) {
		assert false, 'a viscosity of ${s.value}'
	} else {
		assert err.msg().contains("is equation 118, which reflux can't evaluate"), err.msg()
	}
}

// The constants transport takes, in SI: Vc and Zc from the file, the
// dipole moment in debye, Lucas's Q for the quantum gases alone; a constant
// the file hasn't is NaN.
fn test_constants() {
	db := components.load_default()!
	w := db.find('water')!.transport_pure()!
	assert w.cas == '7732-18-5' && w.mw == 0.0180153 && w.vc == 5.595e-05 && w.zc == 0.229
	assert math.abs(w.dipole - 6.17e-30 / 3.33564e-30) <= 1e-15
	assert w.quantum == 0
	assert db.find('hydrogen')!.transport_pure()!.quantum == 0.76
	assert db.find('helium-4')!.transport_pure()!.quantum == 1.38
	assert math.is_nan(db.find('100-18-5')!.transport_pure()!.dipole)
	assert w.sources.len == transport.quantities.len
	for i, q in transport.quantities {
		assert w.sources[i].has && w.sources[i].c.key == q.key()
	}
}

// A number in a message is to 6 significant figures, the double nearest
// that decimal, so it prints as the decimal; 0, NaN and ±∞ as they are.
fn test_a_messages_number_is_shown_to_6_figures() {
	for x, want in {
		0.13192658013779612: 0.131927
		320.00000000000006:  320.0
		318.15000000000003:  318.15
		1.2345649999e-5:     1.23456e-5
		1.2345650001e-5:     1.23457e-5
		-98765.4321:         -98765.4
		6.02214076e23:       6.02214e23
		1e-300:              1e-300
		0.0:                 0.0
		999999.5:            1e6
	} {
		got := transport.shown(x)
		assert got == want, '${x}: ${got}, not ${want}'
		assert got.str() == want.str()
	}
	assert math.is_nan(transport.shown(math.nan()))
	assert transport.shown(math.inf(1)) == math.inf(1)
}

// How many compounds have each quantity reflux can use, and all six: numbers
// that can change only with the data (D3, step 1's coverage). Of the 638
// liquid viscosities in the files, 15 have impossible ranges and one is
// equation 118 (data/components/CHANGES.md); 594 compounds have all six in
// their files, 588 usable.
fn test_coverage() {
	db := components.load_default()!
	mut each := []int{len: transport.quantities.len}
	mut all, mut made := 0, 0
	for e in db.list() {
		p := db.find(e.cas)!.transport_pure() or { continue }
		made++
		mut every := true
		for i, s in p.sources {
			if s.has {
				each[i]++
			} else if i < 6 {
				every = false
			}
		}
		if every {
			all++
		}
	}
	println('${made} compounds; with each: ${each}; with all six: ${all}')
	assert made == 641
	assert each == [622, 606, 599, 598, 637, 626, 641]
	assert all == 588
}
