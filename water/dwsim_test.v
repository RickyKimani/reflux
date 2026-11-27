// Internal test (module water): DWSIM 9's two packages for water by state
// (M5a step 9, D14), testdata/m5_dwsim.json from tools/gen_m5_dwsim.py.
// Its "CoolProp" package is IAPWS-95 by the CoolProp library, held to
// reflux's IAPWS-95; its "Steam Tables (IAPWS-IF97)" to reflux's IF97. Every
// difference is either within the package's own iteration or a failure of
// DWSIM's, shown to be one by its own numbers, and each failure is listed:
// a list that can't grow.
module water

import flash
import json2
import math
import os
import thermo
import water.if97

fn dwsim() map[string]json2.Any {
	text := os.read_file(os.join_path(os.dir(@FILE), '..', 'testdata', 'm5_dwsim.json')) or {
		panic(err)
	}
	return (json2.decode[json2.Any](text) or { panic(err) }).as_map()
}

fn field(m map[string]json2.Any, k string) f64 {
	v := m[k] or { return math.nan() }
	return v.f64()
}

fn rows(d map[string]json2.Any, pkg string, key string) []map[string]json2.Any {
	p := d[pkg] or { panic('no ${pkg}') }
	list := p.as_map()[key] or { panic('no ${pkg}.${key}') }
	return list.as_array().map(it.as_map())
}

// h95 and s95 are IAPWS-95's per kg at (T, P) on IAPWS's reference, kJ/kg
// and kJ/(kg K), as DWSIM reports them; a wet state's by its phases.
fn hs95(m thermo.PhaseModel, t f64, p f64) (f64, f64) {
	r := flash.pt(m, t, p, [1.0]) or { return math.nan(), math.nan() }
	mut h, mut s := 0.0, 0.0
	for ph in r.phases {
		pr := props95(t, 1.0 / ph.state.v)
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
	}
	return h / molar_mass * 1e-3, s / molar_mass * 1e-3
}

fn key_of(t f64, p f64) string {
	return '${t} K ${p} Pa'
}

// The CoolProp package against reflux's IAPWS-95 at all 98 states.
//
// Measured: h within 7.5e-7 kJ/kg (3.0e-10 of h, at 800 K and 100 MPa) and s
// within 1.2e-9 kJ/(kg K), five states above 1e-8 kJ/kg, all gas or
// supercritical fluid at 10 to 100 MPa, where DWSIM's copy of the CoolProp
// library stops its density solve (CoolProp 8.0.0's own densities agree
// with reflux's to 3.3e-12, testdata/m5_iapws95.json).
//
// Its PH flash back lands within 1e-4 K but at four states near the critical
// point and in the dense supercritical fluid, where it returns, with no
// error, a T 0.23 K to 1.5 K off, whose h is 5.9 to 90 kJ/kg (0.24% to
// 4.0%) from the h it was given: DWSIM's single-compound flash, which the
// package uses for a PH flash, not CoolProp's own.
fn test_coolprop_package_is_iapws95() {
	d := dwsim()
	m := thermo.by_roots(iapws95())
	mut misses := []string{}
	mut n := 0
	for row in rows(d, 'coolprop', 'states') {
		assert 'error' !in row, '${row}'
		t, p := field(row, 'T'), field(row, 'P')
		h, s := field(row, 'h'), field(row, 's')
		hr, sr := hs95(m, t, p)
		assert math.abs(h - hr) <= 1e-6, '${key_of(t, p)}: h ${h} against ${hr}'
		assert math.abs(s - sr) <= 2e-9, '${key_of(t, p)}: s ${s} against ${sr}'
		back := field(row, 'ph_T')
		if !(math.abs(back - t) <= 1e-4) {
			at, _ := hs95(m, back, p)
			// a wrong answer, not a loose one: its T's h is over 0.1% off
			assert math.abs(at - h) > 1e-3 * h, '${key_of(t, p)}: PH back at ${back} K'
			misses << key_of(t, p)
		}
		n++
	}
	assert n == 98
	assert misses == [key_of(660.0, 2.5e7), key_of(800.0, 1e8), key_of(647.0, 2.2e7),
		key_of(648.0, 2.21e7)], '${misses}'
}

// The Steam Tables package against reflux's IF97 at all 98 states.
//
// Measured: within 2.1e-11 kJ/kg in h and 3.2e-14 in s wherever it answers
// right: regions 1 and 2, and region 3 at 640 K (25 and 50 MPa), 648 K and
// 700 K. Its failures, every one silent:
// - region 3 at 650 K, 25 to 100 MPa: h = s = −1, a sentinel, returned as
//   a value;
// - region 3 at 640 K and 100 MPa: h 539 kJ/kg low;
// - region 5: 5.75 kJ/kg low at 1500 K and 0.5 MPa, 3875 at 30 MPa, and
//   −46,354 kJ/kg at 2000 K and 30 MPa.
// Its PH flash back is within 1e-3 K (4.3e-4 measured) but at every state of
// 25 MPa and above up to 650 K, and in region 5: its PH flash runs Brent on
// H − h(T) between 273.15 K and 2000 K (SteamTables.vb, DW_CalcEquilibrium),
// on its own h, which is −1 in region 3 and wrong in region 5, so the
// bracket doesn't hold one root, and it returns 1729.42 K (25 to 30 MPa),
// 2000 K (50 MPa) or about 1525 K (100 MPa), whatever h it's given.
fn test_steam_tables_package_is_if97_where_it_answers() {
	d := dwsim()
	mut wrong := []string{}
	mut sentinel := []string{}
	mut ph_misses := []string{}
	mut n := 0
	for row in rows(d, 'steam_tables', 'states') {
		assert 'error' !in row, '${row}'
		t, p := field(row, 'T'), field(row, 'P')
		h, s := field(row, 'h'), field(row, 's')
		st := if97.state_pt(p, t)!
		if h == -1.0 && s == -1.0 {
			sentinel << key_of(t, p)
		} else if !(math.abs(h - st.h * 1e-3) <= 1e-10 && math.abs(s - st.s * 1e-3) <= 1e-13) {
			assert math.abs(h - st.h * 1e-3) > 1.0, '${key_of(t, p)}: h ${h} against ${st.h * 1e-3}'
			wrong << '${key_of(t, p)} region ${st.region}'
		}
		back := field(row, 'ph_T')
		if !(math.abs(back - t) <= 1e-3) {
			// its T's h, where IF97 has one, is far from the h it was given
			if at := if97.state_pt(p, back) {
				assert math.abs(at.h * 1e-3 - h) > 10.0, '${key_of(t, p)}: PH back at ${back} K'
			}
			ph_misses << key_of(t, p)
		}
		n++
	}
	assert n == 98
	assert sentinel == [key_of(650.0, 2.5e7), key_of(650.0, 5e7), key_of(650.0, 1e8),
		key_of(650.0, 2.55837018e7)], '${sentinel}'
	assert wrong == ['${key_of(640.0, 1e8)} region 3', '${key_of(1500.0, 5e5)} region 5',
		'${key_of(1500.0, 3e7)} region 5', '${key_of(2000.0, 3e7)} region 5'], '${wrong}'
	mut want := []string{}
	for t in [280.0, 300.0, 350.0, 400.0, 500.0, 600.0, 620.0, 640.0, 650.0] {
		for p in [2.5e7, 5e7, 1e8] {
			want << key_of(t, p)
		}
	}
	want << [key_of(650.0, 2.55837018e7), key_of(1500.0, 5e5), key_of(1500.0, 3e7),
		key_of(2000.0, 3e7)]
	assert ph_misses == want, '${ph_misses}'
}

// Inside the dome and on it: each package's saturated phases' h and its Ts
// at 0.1 and 1 MPa, its PH flash at a quality back to that quality, and its
// bubble temperature at five pressures, against reflux's IAPWS-95 and IF97.
//
// Measured: CoolProp's Ts within 2.0e-12 K, its saturated h within
// 1.3e-10 kJ/kg and its bubble T within 6.5e-11 K (at 20 MPa, its own
// iteration); Steam Tables' Ts and bubble T exactly IF97's Eq. (31), its h
// within 9.1e-13 kJ/kg; each PH flash's x within 1.1e-16 of the quality.
// (println gives the worsts.)
fn test_each_package_in_the_dome() {
	d := dwsim()
	m := thermo.by_roots(iapws95())
	mut worst := map[string]f64{}
	for pkg in ['coolprop', 'steam_tables'] {
		for row in rows(d, pkg, 'dome') {
			assert 'error' !in row, '${row}'
			p, x := field(row, 'P'), field(row, 'x')
			mut ts, mut hl, mut hv := 0.0, 0.0, 0.0
			if pkg == 'coolprop' {
				r := flash.pvf(m, p, 0.0, [1.0])!
				ts = r.t
				for ph in r.phases {
					h := props95(r.t, 1.0 / ph.state.v).h / molar_mass * 1e-3
					if ph.label == .vapor {
						hv = h
					} else {
						hl = h
					}
				}
			} else {
				l, v := if97.state_px(p, 0.0)!, if97.state_px(p, 1.0)!
				ts, hl, hv = l.t, l.h * 1e-3, v.h * 1e-3
			}
			worst['${pkg} T'] = math.max(worst['${pkg} T'], math.abs(field(row, 'T') - ts))
			worst['${pkg} h'] = math.max(worst['${pkg} h'], math.max(math.abs(field(row, 'h_liquid') - hl),
				math.abs(field(row, 'h_vapour') - hv)))
			worst['${pkg} x'] = math.max(worst['${pkg} x'], math.abs(field(row, 'beta') - x))
			assert math.abs(field(row, 'T') - ts) <= 1e-10, '${pkg} at ${p} Pa: Ts ${field(row, 'T')} against ${ts}'
			assert math.abs(field(row, 'h_liquid') - hl) <= 1e-9, '${pkg} at ${p} Pa: h′'
			assert math.abs(field(row, 'h_vapour') - hv) <= 1e-9, '${pkg} at ${p} Pa: h″'
			assert math.abs(field(row, 'beta') - x) <= 1e-12, '${pkg} at ${p} Pa, x ${x}'
		}
		for row in rows(d, pkg, 'bubble') {
			p := field(row, 'P')
			ts := if pkg == 'coolprop' { flash.pvf(m, p, 0.0, [1.0])!.t } else { if97.tsat(p)! }
			worst['${pkg} bubble'] = math.max(worst['${pkg} bubble'], math.abs(field(row, 'T') - ts))
			assert math.abs(field(row, 'T') - ts) <= 1e-9, '${pkg} bubble at ${p} Pa: ${field(row, 'T')} against ${ts}'
		}
	}
	println(worst)
}
