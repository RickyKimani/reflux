// Internal test (module if97): IF97 against the iapws package's
// (testdata/m5_if97.json, tools/gen_m5_if97.py) over grids in every region
// (M5a's criterion 5): the basic equations, the backward equations as
// published and the boundaries to rounding, and the steam-table functions
// within the package's own iteration (SciPy's newton and fsolve at their
// default tolerances). Each bound is the measured worst with a margin, a
// number that can't grow.
module if97

import json2
import math
import os

fn fixture_states() !json2.Any {
	path := os.join_path(os.dir(@FILE), '..', '..', 'testdata', 'm5_if97.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m5_if97.py`.')
	}
	root := json2.decode[json2.Any](text)!
	return root.as_map()['states'] or { return error('${path}: no states') }
}

fn rows_of(a json2.Any, key string) [][]f64 {
	list := a.as_map()[key] or { return [][]f64{} }
	return list.as_array().map(fn (r json2.Any) []f64 {
		return r.as_array().map(fn (x json2.Any) f64 {
			return if x is json2.Null { math.nan() } else { x.f64() }
		})
	})
}

fn rel(a f64, b f64) f64 {
	return math.abs(a / b - 1.0)
}

// The package's row: region, T, P (MPa), v, h (kJ/kg), s, cp, cv
// (kJ/(kg K)), w, x; reflux's state in the same units.
fn ours_row(st State) []f64 {
	return [f64(st.region), st.t, st.p / 1e6, st.v, st.h / 1e3, st.s / 1e3, st.cp / 1e3, st.cv / 1e3,
		st.w, st.x]
}

// diff is the largest relative difference of T, v, cp, cv and w between
// two rows, and of h and s against the larger of their size and 1 kJ/kg
// (or kJ/(kg K)), so that the triple point's h ≈ 0 doesn't divide.
fn diff(a []f64, b []f64) f64 {
	mut w := 0.0
	for k in [1, 3, 6, 7, 8] {
		if math.is_nan(a[k]) && math.is_nan(b[k]) {
			continue
		}
		w = math.max(w, rel(a[k], b[k]))
	}
	for k in [4, 5] {
		w = math.max(w, math.abs(a[k] - b[k]) / math.max(math.abs(b[k]), 1.0))
	}
	return w
}

fn worst_of(m map[int]f64, k int, x f64) map[int]f64 {
	mut out := m.clone()
	out[k] = math.max(m[k] or { 0.0 }, x)
	return out
}

// The basic equations, from (p, T), region by region. Regions 1, 2 and 5:
// the same sums, ours by repeated squaring and the package's by NumPy's
// powers, measured 7.5e-13, 3.1e-15 and 1.1e-15. Region 3: both Newton on
// p(ρ, T) = p, measured 2.6e-12.
fn test_states_from_p_and_t() {
	f := fixture_states()!
	mut worst := map[int]f64{}
	mut n := 0
	for full in rows_of(f, 'pt') {
		row := full[2..]
		st := state_pt(full[0] * 1e6, full[1])!
		o := ours_row(st)
		assert int(o[0]) == int(row[0]), '${full[0]} MPa ${full[1]} K: region ${o[0]}, the package ${row[0]}'
		worst = worst_of(worst, int(row[0]), diff(o, row))
		n++
	}
	assert n == 831
	bound := {
		1: 3e-12
		2: 2e-14
		3: 1e-11
		5: 1e-14
	}
	for reg, w in worst {
		assert w <= bound[reg], 'region ${reg}: ${w}'
	}
}

// From (p, h) and (p, s) at the grid's states, against the package's
// answers. Ours are Newton on the basic equation to rounding; the
// package's SciPy's newton (1.48e-8 in T) and, in region 3, fsolve (1.49e-8
// relative). Measured: regions 1, 2 and 5 within 2.2e-12, region 3 within
// 6.7e-11. The package refuses 8 of region 3's 33 states from (p, h) and 9
// from (p, s): its fsolve fails there.
fn test_states_from_p_and_h_or_s() {
	f := fixture_states()!
	pt := rows_of(f, 'pt')
	for key in ['ph', 'ps'] {
		mut worst := map[int]f64{}
		for row in rows_of(f, key) {
			src := pt[int(row[0])]
			theirs := row[1..]
			st := if key == 'ph' {
				state_ph(src[0] * 1e6, src[6] * 1e3)!
			} else {
				state_ps(src[0] * 1e6, src[7] * 1e3)!
			}
			o := ours_row(st)
			// every state comes back in the region of its (p, T), as the
			// package's does
			assert int(o[0]) == int(src[2]), '${key} ${src[0]} MPa ${src[1]} K: region ${o[0]}'
			assert int(o[0]) == int(theirs[0]), '${key} ${src[0]} MPa ${src[1]} K: the package ${theirs[0]}'
			worst = worst_of(worst, int(o[0]), diff(o, theirs))
		}
		for reg, w in worst {
			assert w <= if reg == 3 { 2e-10 } else { 1e-11 }, '${key} region ${reg}: ${w}'
		}
	}
}

// Saturated states from (T, x) and (p, x), x of 0, 0.35 and 1. Up to
// 623.15 K both are regions 1 and 2 at region 4's (p, T): measured within
// 1e-12. Above it the phases are region 3's:
// - the package's (T, x) of 0 and 1 saturates by Maxwell's criterion on
//   Eq. (28), not at Eq. (30)'s p as ours does (D9): they differ by up to
//   5.0e-6 in p, within Table 44's 0.0026% between the two, and 1.7e-6 in
//   v, h and s. Its mixtures there sit between ends at two pressures, so
//   aren't compared;
// - its mixtures from (p, x) use SR5-05's backward densities unfinished: at
//   22 MPa and x = 0.35, v is 1% off the mixture of its own saturated
//   phases, which agree with ours within 1e-10. So a mixture of ours is held
//   to the package's two ends, mixed.
fn test_saturated_states() {
	f := fixture_states()!
	for key in ['tx', 'px'] {
		rows := rows_of(f, key)
		mut low, mut high_p, mut high := 0.0, 0.0, 0.0
		for k, row in rows {
			st := if key == 'tx' {
				state_tx(row[1], row[9])!
			} else {
				state_px(row[2] * 1e6, row[9])!
			}
			o := ours_row(st)
			mut want := row.clone()
			if row[1] > t13 && row[9] == 0.35 {
				if key == 'tx' {
					continue
				}
				// the package's ends are the rows either side
				a, b := rows[k - 1], rows[k + 1]
				for j in [3, 4, 5] {
					want[j] = a[j] + 0.35 * (b[j] - a[j])
				}
			}
			dp := math.max(rel(o[1], want[1]), rel(o[2], want[2]))
			mut w := 0.0
			for j in [3, 4, 5] {
				w = math.max(w, math.abs(o[j] - want[j]) / math.max(math.abs(want[j]), 1.0))
			}
			if row[1] <= t13 {
				low = math.max(low, math.max(dp, w))
			} else {
				high_p = math.max(high_p, dp)
				high = math.max(high, w)
			}
		}
		assert low <= 1e-11, '${key} up to 623.15 K: ${low}'
		if key == 'tx' {
			assert high_p <= 2.6e-5, 'tx above 623.15 K, p: ${high_p}'
			assert high <= 5e-6, 'tx above 623.15 K, v, h and s: ${high}'
		} else {
			assert math.max(high_p, high) <= 1e-10, 'px above 623.15 K: ${high_p} ${high}'
		}
	}
}

// The backward equations as published, at each grid state of their region,
// and region 3's v(p, T) on a fine grid near the critical point: the same
// sums, measured within 4.6e-15, and 2.3e-12 for v(p, T) at 22.1 MPa and
// 647.2 K, where subregion 3v's terms cancel.
fn test_backward_equations() {
	f := fixture_states()!
	back := f.as_map()['backward'] or { panic('no backward') }
	for key in ['T1_ph', 'T1_ps', 'T2_ph', 'T2_ps', 'T3_ph', 'T3_ps', 'v3_ph', 'v3_ps', 'v3_pt'] {
		rows := rows_of(back, key)
		assert rows.len > 0
		mut worst := 0.0
		for row in rows {
			p := row[0] * 1e6
			got := match key {
				'T1_ph' { t1_ph(p, row[1] * 1e3) }
				'T1_ps' { t1_ps(p, row[1] * 1e3) }
				'T2_ph' { t2_ph(p, row[1] * 1e3) }
				'T2_ps' { t2_ps(p, row[1] * 1e3) }
				'T3_ph' { t3_ph(p, row[1] * 1e3) }
				'T3_ps' { t3_ps(p, row[1] * 1e3) }
				'v3_ph' { v3_ph(p, row[1] * 1e3) }
				'v3_ps' { v3_ps(p, row[1] * 1e3) }
				else { v3_backward(p, row[1]) }
			}
			worst = math.max(worst, rel(got, row[2]))
		}
		assert worst <= if key == 'v3_pt' { 5e-12 } else { 2e-14 }, '${key}: ${worst}'
	}
}

// Region 4 and the boundaries: measured within 7.6e-15.
fn test_boundaries() {
	f := fixture_states()!
	mut worst := map[string]f64{}
	for row in rows_of(f, 'psat') {
		worst['psat'] = math.max(worst['psat'] or { 0.0 }, rel(psat97(row[0]) / 1e6, row[1]))
	}
	for row in rows_of(f, 'tsat') {
		worst['tsat'] = math.max(worst['tsat'] or { 0.0 }, rel(tsat97(row[0] * 1e6), row[1]))
	}
	for row in rows_of(f, 'p23') {
		worst['p23'] = math.max(worst['p23'] or { 0.0 }, rel(p23(row[0]) / 1e6, row[1]))
	}
	for row in rows_of(f, 'p3sat_h') {
		worst['p3sat_h'] = math.max(worst['p3sat_h'] or { 0.0 }, rel(p3sat_h(row[0] * 1e3) / 1e6,
			row[1]))
	}
	for row in rows_of(f, 'p3sat_s') {
		worst['p3sat_s'] = math.max(worst['p3sat_s'] or { 0.0 }, rel(p3sat_s(row[0] * 1e3) / 1e6,
			row[1]))
	}
	for row in rows_of(f, 'h3ab') {
		worst['h3ab'] = math.max(worst['h3ab'] or { 0.0 }, rel(h3ab(row[0] * 1e6) / 1e3, row[1]))
	}
	assert worst.len == 6
	for k, w in worst {
		assert w <= 2e-14, '${k}: ${w}'
	}
}
