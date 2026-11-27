// IF97's steam-table functions (M5a's criterion 5): each backward answer
// finished on its basic equation, so that a state given by (p, h) or (p, s)
// is the state at (p, T) to rounding; the dome by region 4 and the
// saturated phases; each refusal with its reason.
module if97

import math

// grid is states from (p, T) in every region: T every 10 K (20 K above
// 1073.15 K), and p from 1 Pa to 100 MPa, each kept in IF97's range.
fn grid() []State {
	mut out := []State{}
	ps := [1.0, 100.0, 1e3, 1e4, 1e5, 5e5, 1e6, 3e6, 5e6, 1e7, 1.5e7, 1.7e7, 2e7, 2.2e7, 2.3e7,
		2.5e7, 3e7, 4e7, 5e7, 7e7, 1e8]
	mut t := t_min
	for t <= t_max {
		for p in ps {
			st := state_pt(p, t) or { continue }
			out << st
		}
		t += if t < t5 { 10.0 } else { 20.0 }
	}
	return out
}

// Every state of the grid back from (p, h) and from (p, s): the same
// region, and T and v to rounding. A Newton finish stops where the residual's
// rounding moves T: measured worst 3.6e-15 of T in regions 1, 2 and 5, and
// 4.4e-14 in region 3, whose Newton in (ρ, T) is less well conditioned near
// the critical point; v within 3.9e-14 and 7.6e-14.
fn test_back_from_h_and_s() {
	states := grid()
	assert states.len == 2841
	mut worst := map[int][]f64{}
	for st in states {
		b := state_ph(st.p, st.h)!
		c := state_ps(st.p, st.s)!
		assert b.region == st.region && c.region == st.region, '${st.p} Pa ${st.t} K: ${b.region} ${c.region}'
		assert b.p == st.p && c.p == st.p
		e := [math.abs(b.t / st.t - 1.0), math.abs(c.t / st.t - 1.0), math.abs(b.v / st.v - 1.0),
			math.abs(c.v / st.v - 1.0)]
		mut w := worst[st.region] or { [0.0, 0.0, 0.0, 0.0] }
		for k in 0 .. 4 {
			w[k] = math.max(w[k], e[k])
		}
		worst[st.region] = w
	}
	for reg, w in worst {
		t_bound := if reg == 3 { 1e-13 } else { 1e-14 }
		v_bound := if reg == 3 { 2e-13 } else { 1e-13 }
		assert w[0] <= t_bound && w[1] <= t_bound, 'region ${reg} T: ${w}'
		assert w[2] <= v_bound && w[3] <= v_bound, 'region ${reg} v: ${w}'
	}
}

// The dome: a state of quality x at p, back from its h and its s, gives x
// and Ts back; from 611.213 Pa to 22 MPa, both sides of 623.15 K.
fn test_the_dome() {
	ps := [611.3, 1e3, 1e4, 1e5, 1e6, 5e6, 1e7, 1.6e7, 1.66e7, 1.7e7, 1.9e7, 2.1e7, 2.2e7, 2.205e7]
	mut worst_x := 0.0
	for p in ps {
		ts := tsat97(p)
		for x in [0.05, 0.25, 0.5, 0.75, 0.95] {
			st := state_px(p, x)!
			assert st.region == 4 && st.t == ts && st.p == p
			assert math.is_nan(st.cp) && math.is_nan(st.w)
			for b in [state_ph(p, st.h)!, state_ps(p, st.s)!] {
				assert b.region == 4 && b.t == ts, '${p} Pa x ${x}: region ${b.region}'
				worst_x = math.max(worst_x, math.abs(b.x - x))
			}
		}
	}
	// x from h is (h − h')/(h'' − h'), its rounding a few ε·h/(h'' − h'), and
	// h'' − h' is 4% of h at 22.05 MPa: measured 4.2e-15.
	assert worst_x <= 2e-14, '${worst_x}'
}

// The saturated phases, by (T, x) of 0 and 1: each a single phase's
// state, and in region 3 each density on Eq. (28) at Eq. (30)'s p, the
// liquid's above ρc and the vapour's below, to 1 µK below Tc.
fn test_saturated_phases() {
	for t in [273.15, 300.0, 400.0, 500.0, 600.0, 623.15, 623.16, 630.0, 640.0, 645.0, 647.0, 647.09,
		647.095, 647.0959, 647.09599, 647.095999] {
		l := state_tx(t, 0.0)!
		v := state_tx(t, 1.0)!
		assert l.region == 4 && v.region == 4 && l.x == 0.0 && v.x == 1.0
		assert l.p == psat97(t) && v.p == l.p && l.t == t
		assert l.rho() > v.rho() && math.is_finite(l.cp) && math.is_finite(v.w)
		if t > t13 {
			assert l.rho() > rhoc && v.rho() < rhoc, '${t} K: ${l.rho()} ${v.rho()}'
			for st in [l, v] {
				assert math.abs(region3(st.rho(), t).p / st.p - 1.0) <= 1e-9
			}
		}
	}
	// At Tc both are the critical point, Eq. (28)'s exactly.
	c := state_tx(tc, 0.5)!
	assert c.rho() == rhoc && c.p == pc && c.t == tc
	d := state_px(pc, 0.0)!
	assert d.rho() == rhoc && d.t == tc
}

fn kind(e IError) Failure {
	return (e as If97Error).kind
}

// Each function's failure kind at a state, or .no_answer for an answer
// where a refusal was due.
fn pt_fails(p f64, t f64) Failure {
	_ := state_pt(p, t) or { return kind(err) }
	return .no_answer
}

fn ph_fails(p f64, h f64) Failure {
	_ := state_ph(p, h) or { return kind(err) }
	return .no_answer
}

fn ps_fails(p f64, s f64) Failure {
	_ := state_ps(p, s) or { return kind(err) }
	return .no_answer
}

fn tx_fails(t f64, x f64) Failure {
	_ := state_tx(t, x) or { return kind(err) }
	return .no_answer
}

fn px_fails(p f64, x f64) Failure {
	_ := state_px(p, x) or { return kind(err) }
	return .no_answer
}

// Each state outside the formulation's range is refused as such, saying
// what the range is.
fn test_refusals() {
	kinds := [pt_fails(100.000001e6, 300.0), pt_fails(50.0001e6, 1500.0), pt_fails(1e5, 273.14),
		pt_fails(1e5, 2273.16), pt_fails(0.0, 400.0), pt_fails(-1.0, 400.0),
		pt_fails(math.nan(), 400.0), ph_fails(1e5, -100e3), ph_fails(1e5, 1e8), ph_fails(101e6, 1e6),
		ps_fails(1e5, -1e3), ps_fails(1e5, 50e3), tx_fails(273.0, 0.5), tx_fails(647.1, 0.5),
		tx_fails(400.0, 1.01), tx_fails(400.0, -0.01), px_fails(600.0, 0.5), px_fails(22.065e6, 0.5),
		px_fails(1e5, math.nan())]
	for k, f in kinds {
		assert f == .out_of_range, 'case ${k}: ${f}'
	}
	mut msg := ''
	_ := state_pt(100.000001e6, 300.0) or {
		msg = err.msg()
		State{}
	}
	assert msg.contains('273.15 K to 1073.15 K up to 100 MPa, and to 2273.15 K up to 50 MPa')
	_ := psat(273.0) or { 0.0 }
	_ := tsat(22.1e6) or { 0.0 }
	assert (psat(273.0) or { -1.0 }) == -1.0 && (tsat(22.1e6) or { -1.0 }) == -1.0
}

// The range's edges are in it: 273.15 K and 2273.15 K from their own h and
// s, though each is a hair outside by rounding.
fn test_the_edges() {
	for p in [1e3, 1e5, 1e7, 1e8] {
		a := state_pt(p, t_min)!
		assert math.abs(state_ph(p, a.h)!.t - t_min) < 1e-9
		assert math.abs(state_ps(p, a.s)!.t - t_min) < 1e-9
	}
	for p in [1e3, 1e5, 1e7, 5e7] {
		a := state_pt(p, t_max)!
		assert math.abs(state_ph(p, a.h)!.t - t_max) < 1e-9
		assert math.abs(state_ps(p, a.s)!.t - t_max) < 1e-9
	}
}

// rho3 refuses a pressure with no root on its side. At 647 K the vapour's
// branch of Eq. (28) ends at its spinodal 555 Pa above Eq. (30)'s p, and
// the liquid's 646 Pa below it: 1 kPa past each, the search closes on the
// spinodal with p off by 445 Pa and 354 Pa (2e-5 of p), and refuses.
fn test_no_root_on_a_side() {
	ps := psat97(647.0)
	assert (rho3(ps + 1e3, 647.0, .vapour) or { -1.0 }) == -1.0
	assert (rho3(ps - 1e3, 647.0, .liquid) or { -1.0 }) == -1.0
	// and either side of those, each has its root
	assert rho3(ps + 500.0, 647.0, .vapour)! < rhoc
	assert rho3(ps - 600.0, 647.0, .liquid)! > rhoc
}

// The critical point from (p, T): Eq. (28) reproduces it exactly, and rho3
// gives ρc, where ∂p/∂ρ is 0 and Newton has no step.
fn test_the_critical_point() {
	c := state_pt(pc, tc)!
	assert c.region == 3 && c.rho() == rhoc && c.p == pc && c.t == tc
}

// A vapour of subregion 3k below Tc: between 20.5 and about 21.6 MPa, T3jk
// is below Tc, so 3k holds vapours there too, on the vapour's side of ρc.
fn test_a_vapour_in_3k() {
	p, t := 21e6, 645.0
	assert subregion3(p, t, .by_t) == `k` && t < tc && t > tsat97(p)
	st := state_pt(p, t)!
	assert st.region == 3 && st.rho() < rhoc
	assert math.abs(region3(st.rho(), t).p / p - 1.0) <= 1e-12
}

// A value just outside the range's edge by its rounding is taken; one past
// that room is refused.
fn test_a_hair_outside_the_edge() {
	for p in [1e3, 1e5, 1e7, 1e8] {
		a := state_pt(p, t_min)!
		inside := a.h - 0.5e-12 * r * t_min
		assert math.abs(state_ph(p, inside)!.t - t_min) < 1e-9
		assert ph_fails(p, a.h - 2e-12 * r * t_min) == .out_of_range
		assert math.abs(state_ps(p, a.s - 0.5e-12 * r)!.t - t_min) < 1e-9
	}
}

// Between region 2's h at 1073.15 K and region 5's, where region 5's is the
// higher (by 15 to 90 J/kg at 1 kPa, 0.1 MPa, 30 MPa and 50 MPa, Table 43's
// inconsistency), a state is region 5's a little below 1073.15 K.
fn test_the_gap_at_1073_k() {
	mut n := 0
	for p in [1e3, 1e5, 3e7, 5e7] {
		a := region2(p, t5)
		b := region5(p, t5)
		assert b.h > a.h
		st := state_ph(p, 0.5 * (a.h + b.h))!
		assert st.region == 5 && st.t < t5 && st.t > t5 - 1.0
		n++
	}
	assert n == 4
}
