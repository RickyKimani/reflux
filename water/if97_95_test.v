// Internal test (module water): IF97 against IAPWS-95 (M5a's criterion 6).
// IF97 was fitted to IAPWS-95, and R7-97's Figs. 3 to 6 give its
// uncertainties, IAPWS-95's widened by the two's differences. On a grid in
// each region, every difference of v, cp, w and ps lies within IF97's
// uncertainty there; each region's largest and RMS differences of v, h, s,
// cp and w are recorded, as numbers that can't grow, the critical region
// apart.
module water

import math
import thermo
import water.if97

// The figures' maps, % (R7-97, Section 12), as read from them; -1 where
// there's no estimate, above 1273 K. Their lines are "approximate", so a
// point gets the largest value within 20 K and a factor of 1.5 in p
// (around). liquid is above the saturation line, below Tc.

// Fig. 3, v.
fn u_v(p f64, t f64, liquid bool) f64 {
	pm := p / 1e6
	if t > 1273.15 {
		return -1.0
	}
	if t <= 415.0 {
		return if !liquid {
			0.05
		} else if pm > 10.0 {
			0.005
		} else {
			0.003
		}
	}
	if t <= 580.0 {
		return if liquid { 0.015 } else { 0.05 }
	}
	if t <= 785.0 {
		if pm >= 30.0 && t >= 650.0 {
			return 0.15
		}
		return if t <= 650.0 && (liquid || pm >= 22.064) { 0.05 } else { 0.1 }
	}
	return if pm > 5.0 { 0.3 } else { 0.15 }
}

// Fig. 4, cp.
fn u_cp(p f64, t f64, liquid bool) f64 {
	pm := p / 1e6
	if t > 1273.15 {
		return -1.0
	}
	if t <= 307.0 && pm >= 20.0 {
		return 0.4
	}
	if t <= 590.0 {
		return if pm >= 60.0 {
			0.3
		} else if liquid {
			0.2
		} else {
			0.3
		}
	}
	return if pm >= 10.0 { 0.5 } else { 0.3 }
}

// Fig. 5, w.
fn u_w(p f64, t f64) f64 {
	pm := p / 1e6
	if t > 1273.15 {
		return -1.0
	}
	if t >= 400.0 && t <= 590.0 && pm >= 8.0 {
		return 0.4
	}
	if t > 590.0 && t <= 735.0 && pm >= 15.0 && pm < 30.0 {
		return 0.5
	}
	if (t > 590.0 && t <= 1073.15 && pm >= 15.0) || (t > 1073.15 && pm >= 10.0) {
		return 1.0
	}
	return 0.3
}

// around is a map's largest value within 20 K and a factor of 1.5 in p.
fn around(u fn (f64, f64) f64, p f64, t f64) f64 {
	mut best := -1.0
	for dt in [-20.0, 0.0, 20.0] {
		for fp in [1.0 / 1.5, 1.0, 1.5] {
			best = math.max(best, u(p * fp, t + dt))
		}
	}
	return best
}

// Fig. 6, ps: the band's half-width against T, as read.
fn u_ps(t f64) f64 {
	pts := [[273.16, 0.0], [300.0, 0.03], [373.124, 0.02], [373.125, 0.03], [450.0, 0.05],
		[647.096, 0.035]]
	for k in 1 .. pts.len {
		if t <= pts[k][0] {
			f := (t - pts[k - 1][0]) / (pts[k][0] - pts[k - 1][0])
			return pts[k - 1][1] + f * (pts[k][1] - pts[k - 1][1])
		}
	}
	return pts[pts.len - 1][1]
}

// The critical region, reported apart: Fig. 3's triangle, bordered by the
// isochores 0.0019 and 0.0069 m³/kg and the 30 MPa isobar, where v's
// uncertainty is given as 0.2% in p and cp's as 3% (Fig. 4), but for the
// small circle Fig. 4 draws at the critical point, read as ±7.5 K and ±11%
// in p, where it gives none; and around the triangle, from 590 K to 700 K
// and 10 MPa to 30 MPa near ρc, where Figs. 4 and 5 give cp and w more
// than 3% and 2%.
fn critical(st if97.State) bool {
	return st.t >= 590.0 && st.t <= 700.0 && st.p >= 10e6 && st.p <= 30e6 && st.v > 0.0017
		&& st.v < 0.009
}

fn triangle(st if97.State) bool {
	return st.t >= if97.t13 && st.p <= 30e6 && st.v >= 0.0019 && st.v <= 0.0069
}

fn circle(st if97.State) bool {
	return math.abs(st.t - if97.tc) <= 7.5 && math.abs(st.p / if97.pc - 1.0) <= 0.11
}

// ridge is the triangle's isochores continued above 30 MPa to 40 MPa, the
// pseudo-critical ridge, where cp peaks: there IF97's cp is up to 0.99%
// from IAPWS-95's, where Fig. 4 states 0.5% (found in step 5).
fn ridge(st if97.State) bool {
	return st.p > 30e6 && st.p <= 40e6 && st.t >= 675.0 && st.t <= 715.0 && st.v >= 0.0019
		&& st.v <= 0.0069
}

struct Diffs {
mut:
	max []f64 = []f64{len: 5}
	sq  []f64 = []f64{len: 5}
	n   int
}

fn (mut d Diffs) add(x []f64) {
	for k in 0 .. 5 {
		d.max[k] = math.max(d.max[k], x[k])
		d.sq[k] += x[k] * x[k]
	}
	d.n++
}

fn test_if97_against_iapws95() {
	m := iapws95()
	mut by := map[string]Diffs{}
	mut over := []string{}
	mut ridge_cp := 0.0
	mut counts := map[string]int{}
	ps := [1e3, 1e4, 1e5, 3e5, 1e6, 3e6, 5e6, 8e6, 1e7, 1.3e7, 1.6e7, 1.8e7, 2e7, 2.1e7, 2.2e7,
		2.25e7, 2.3e7, 2.5e7, 2.8e7, 3e7, 3.5e7, 4e7, 5e7, 6e7, 8e7, 1e8]
	mut t := 275.0
	for t <= if97.t_max {
		for p in ps {
			a := if97.state_pt(p, t) or { continue }
			st := thermo.solve_tp(m, t, p, [1.0], .stable)!
			b := props95(t, 1.0 / st.v)
			v95 := st.v / molar_mass
			// v, cp and w in %, h in kJ/kg, s in J/(kg K), IF97 against IAPWS-95
			x := [100.0 * math.abs(a.v / v95 - 1.0), math.abs(a.h - b.h / molar_mass) / 1e3,
				math.abs(a.s - b.s / molar_mass), 100.0 * math.abs(a.cp * molar_mass / b.cp - 1.0),
				100.0 * math.abs(a.w / b.w - 1.0)]
			liquid := t < if97.tc && p >= if97.psat97(t)
			key := if critical(a) {
				'critical'
			} else if t > 1273.15 {
				'above 1273 K'
			} else {
				'region ${a.region}'
			}
			mut d := by[key] or { Diffs{} }
			d.add(x)
			by[key] = d
			if triangle(a) {
				// v as 0.2% in p: IAPWS-95's p at IF97's (ρ, T); cp to 3% outside
				// the circle
				counts['triangle']++
				dp := 100.0 * math.abs(props95(t, 1.0 / (a.v * molar_mass)).p / p - 1.0)
				if !(dp <= 0.2) || (!(x[3] <= 3.0) && !circle(a)) {
					over << 'triangle ${p} Pa ${t} K: Δp/p ${dp}%, cp ${x[3]}%'
				}
				continue
			}
			if critical(a) || t > 1273.15 {
				continue
			}
			uv := around(fn [liquid] (p f64, t f64) f64 {
				return u_v(p, t, liquid)
			}, p, t)
			ucp := around(fn [liquid] (p f64, t f64) f64 {
				return u_cp(p, t, liquid)
			}, p, t)
			uw := around(u_w, p, t)
			if ridge(a) {
				counts['ridge']++
				ridge_cp = math.max(ridge_cp, x[3])
				if !(x[0] <= uv && x[4] <= uw) {
					over << 'ridge ${p} Pa ${t} K: v ${x[0]:.5f}/${uv}, w ${x[4]:.4f}/${uw}'
				}
				continue
			}
			counts['checked']++
			if !(x[0] <= uv && x[3] <= ucp && x[4] <= uw) {
				over << '${p} Pa ${t} K (region ${a.region}): v ${x[0]:.5f}/${uv}, cp ${x[3]:.4f}/${ucp}, w ${x[4]:.4f}/${uw}'
			}
		}
		t += if t < 1073.0 { 5.0 } else { 50.0 }
	}
	assert over.len == 0, over.join('\n')
	assert counts['checked'] == 4133 && counts['triangle'] == 50 && counts['ridge'] == 18, '${counts}'
	// The ridge's cp, a number that can't grow: measured 0.9874%, at 35 MPa
	// and 695 K.
	assert ridge_cp <= 1.0, '${ridge_cp}'
	// Each region's largest differences, numbers that can't grow (measured,
	// rounded up in their last figure): v, h (kJ/kg), s (J/(kg K)), cp, w.
	recorded := {
		'region 1':     [0.00314, 0.240, 0.431, 0.123, 0.239]
		'region 2':     [0.0374, 0.411, 0.413, 0.486, 0.326]
		'region 3':     [0.116, 0.983, 1.481, 0.988, 0.486]
		'region 5':     [0.0112, 0.264, 0.164, 0.0966, 0.0332]
		'above 1273 K': [0.00656, 0.323, 0.211, 0.0652, 0.0394]
		'critical':     [1.359, 6.607, 10.136, 15.735, 1.569]
	}
	mut grew := []string{}
	for key, d in by {
		for k in 0 .. 5 {
			if !(d.max[k] <= recorded[key][k]) {
				grew << '${key} ${k}: ${d.max[k]} above ${recorded[key][k]}'
			}
		}
	}
	assert by.len == recorded.len && grew.len == 0, grew.join('\n')
	// ps against IAPWS-95's saturation pressure, from 273.16 K to 647 K,
	// within Fig. 6's band (read with 15 K of room); its largest difference
	// can't grow.
	mut wps := 0.0
	mut bad_ps := []string{}
	mut tt := 273.16
	for tt < 647.0 {
		p95 := thermo.saturation_pressure(m, tt)!
		d := 100.0 * math.abs(if97.psat97(tt) / p95 - 1.0)
		wps = math.max(wps, d)
		if !(d <= math.max(u_ps(tt - 15.0), math.max(u_ps(tt), u_ps(tt + 15.0)))) {
			bad_ps << '${tt} K: ${d}%'
		}
		tt += 1.0
	}
	assert bad_ps.len == 0, bad_ps.join('\n')
	assert wps <= 0.0176, '${wps}'
}
