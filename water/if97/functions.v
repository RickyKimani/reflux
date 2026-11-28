module if97

import math

// The steam-table functions: a state from (p, T), (p, h), (p, s), (T, x)
// or (p, x). Each chooses its region by R7-97's Fig. 1 (and SR3-03's Fig. 2
// in region 3), refuses a state outside the formulation's range saying
// which, and finishes every backward equation's answer on the region's
// basic equation by Newton (D9), so that h(p, T(p, h)) = h to rounding.

const eps = 2.220446049250313e-16

// settled says a Newton step has reached x's rounding: below 4 ulps, or
// small and no longer halving, where the residual's own rounding sets the
// step.
fn settled(step f64, prev f64, x f64) bool {
	a := math.abs(step)
	return a <= 4.0 * eps * math.abs(x) || (a <= 1e-10 * math.abs(x) && a >= 0.5 * math.abs(prev))
}

fn range_text(p f64, what string) string {
	return 'IF97 covers 273.15 K to 1073.15 K up to 100 MPa, and to 2273.15 K up to 50 MPa: not ${what} at ${p} Pa'
}

// state_pt is the state at (p, T).
pub fn state_pt(p f64, t f64) !State {
	if !(t >= t_min && t <= t_max) || !(p > 0.0 && p <= p_max) || (t > t5 && p > p5_max) {
		return out_of_range(range_text(p, '${t} K'))
	}
	if t > t5 {
		return region5(p, t)
	}
	if t <= t13 {
		return if p >= psat97(t) { region1(p, t) } else { region2(p, t) }
	}
	if t <= t23_max && p > p23(t) {
		return at_p(region3(rho3(p, t, .by_t)!, t), p)
	}
	return region2(p, t)
}

// at_p is a region-3 state found by Newton at p, carrying p as asked rather
// than as its density gives it back, which differs by rounding.
fn at_p(st State, p f64) State {
	return State{
		...st
		p: p
	}
}

// rho3 is region 3's density at (p, T): SR5-05's v(p, T) as the start,
// then Newton on p(ρ, T) = p, bracketed. Below the critical pressure the
// isotherm has a loop whose spinodals straddle ρc (as IAPWS-95's do, D7),
// so a liquid keeps to ρ > ρc and a vapour to ρ < ρc; where a step lands in
// the loop, where p falls as ρ rises, the bracket moves past it and the
// search bisects. Near Tc the two roots at Eq. (30)'s p are within a few
// pascals of the loop's ends (0.5 Pa at 647.095 K), where Newton from the
// start alone jumps across the loop.
//
// The answer's p(ρ, T) is within 1e-9 of p, or there's none: to rounding
// everywhere but within about 10 µK of Tc, where Eq. (30)'s p is above the
// vapour's spinodal on Eq. (28) by up to 8e-4 Pa (4e-11 of p), the two
// equations' inconsistency there (Table 44 permits 0.0026%), and the
// search ends on the spinodal.
fn rho3(p f64, t f64, side Side) !f64 {
	if p == pc && t == tc {
		return rhoc // Eq. (28) reproduces the critical point exactly (R7-97)
	}
	code := subregion3(p, t, side)
	sub := p < pc && t < tc
	on_vapour := sub && code in [`r`, `t`, `x`, `z`, `k`]
	on_liquid := sub && !on_vapour
	mut lo := if on_liquid { rhoc } else { 0.0 }
	mut hi := if on_vapour { rhoc } else { math.inf(1) }
	mut rho := math.max(lo, math.min(hi, 1.0 / v3_sub(sub_of(code), p, t)))
	tau := tc / t
	mut prev := math.inf(1)
	for _ in 0 .. 200 {
		delta := rho / rhoc
		f := phi3(delta, tau)
		resid := rho * r * t * delta * f.fd - p
		slope := r * t * (2.0 * delta * f.fd + delta * delta * f.fdd)
		if slope > 0.0 {
			if resid < 0.0 {
				lo = rho
			} else {
				hi = rho
			}
		} else if on_vapour {
			hi = rho // in the loop: the vapour's branch is below
		} else {
			lo = rho // the liquid's above
		}
		if hi - lo <= 4.0 * eps * rho && math.abs(resid) <= 1e-9 * p {
			return rho // the bracket has closed to rounding
		}
		next := if slope > 0.0 { rho - resid / slope } else { math.nan() }
		if !(next >= lo && next <= hi) {
			// a bisection, or a step out when the bracket is open above
			rho = if math.is_inf(hi, 1) { lo * 1.1 } else { 0.5 * (lo + hi) }
			prev = math.inf(1)
			continue
		}
		step := rho - next
		rho = next
		if settled(step, prev, rho) && math.abs(resid) <= 1e-9 * p {
			return rho
		}
		prev = step
	}
	return no_answer('IF97 region 3 at ${p} Pa, ${t} K: no density from subregion 3${code.ascii_str()}, last ${rho} kg/m3')
}

// at_t is a region's state at (p, T): 1, 2 or 5.
fn at_t(region int, p f64, t f64) State {
	return match region {
		1 { region1(p, t) }
		2 { region2(p, t) }
		else { region5(p, t) }
	}
}

// finish_t is Newton on T at fixed p for h (by_s false) or s, from t0, in
// a region where it rises with T, bisecting in [lo, hi] when a step would
// leave the bracket. The bracket has 1e-9 of room at its ends, so that an
// answer on one (273.15 K, say) is reached by Newton, not by halving.
fn finish_t(region int, p f64, t0 f64, lo f64, hi f64, target f64, by_s bool) !State {
	mut a, mut b := lo * (1.0 - 1e-9), hi * (1.0 + 1e-9)
	mut t := math.max(lo, math.min(hi, t0))
	mut prev := math.inf(1)
	for _ in 0 .. 100 {
		st := at_t(region, p, t)
		g := if by_s { st.s } else { st.h }
		slope := if by_s { st.cp / t } else { st.cp }
		f := g - target
		if f > 0.0 {
			b = t
		} else if f < 0.0 {
			a = t
		} else {
			return st
		}
		next := t - f / slope
		if !(next >= a && next <= b) {
			// A step out of a bracket already closed to rounding: h's own
			// rounding, about 6e-9 J/kg in region 1 at 273 K and 25 MPa, is
			// more than the step left, so Newton jitters across it, and
			// halving a closed bracket never settles (M5b bundle B found it:
			// h 1.3e-9 J/kg below h(273.16 K) at 25 MPa had no T).
			if b - a <= 4.0 * eps * t {
				return at_t(region, p, t)
			}
			// a bisection: never taken for rounding's stall
			t = 0.5 * (a + b)
			prev = math.inf(1)
			continue
		}
		step := t - next
		t = next
		if settled(step, prev, t) {
			return at_t(region, p, t)
		}
		prev = step
	}
	what := if by_s { 's = ${target} J/(kg K)' } else { 'h = ${target} J/kg' }
	return no_answer('IF97 region ${region} at ${p} Pa, ${what}: T not found, last ${t} K')
}

// finish3 is Newton in (ρ, T) on region 3's p and h (or s), from a backward
// equation's (ρ0, T0).
fn finish3(p f64, rho0 f64, t0 f64, target f64, by_s bool) !State {
	mut rho, mut t := rho0, t0
	mut prev_r, mut prev_t := math.inf(1), math.inf(1)
	for _ in 0 .. 60 {
		delta := rho / rhoc
		tau := tc / t
		f := phi3(delta, tau)
		rt := r * t
		// p, and its derivatives in ρ and T (Table 31)
		pp := rho * rt * delta * f.fd
		p_r := rt * (2.0 * delta * f.fd + delta * delta * f.fdd)
		p_t := rho * r * (delta * f.fd - delta * tau * f.fdt)
		mut g, mut g_r, mut g_t := 0.0, 0.0, 0.0
		if by_s {
			g = r * (tau * f.ft - f.f)
			g_r = r / rho * (delta * tau * f.fdt - delta * f.fd)
			g_t = -r * tau * tau * f.ftt / t
		} else {
			g = rt * (tau * f.ft + delta * f.fd)
			g_r = rt / rho * (delta * tau * f.fdt + delta * f.fd + delta * delta * f.fdd)
			g_t = r * (delta * f.fd - tau * tau * f.ftt - delta * tau * f.fdt)
		}
		e1 := pp - p
		e2 := g - target
		det := p_r * g_t - p_t * g_r
		mut dr := (e1 * g_t - e2 * p_t) / det
		mut dt := (p_r * e2 - g_r * e1) / det
		if !math.is_finite(dr) || !math.is_finite(dt) {
			break
		}
		dr = math.max(-0.1 * rho, math.min(0.1 * rho, dr))
		dt = math.max(-0.01 * t, math.min(0.01 * t, dt))
		rho -= dr
		t -= dt
		if settled(dr, prev_r, rho) && settled(dt, prev_t, t) {
			return at_p(region3(rho, t), p)
		}
		prev_r, prev_t = dr, dt
	}
	what := if by_s { 's = ${target} J/(kg K)' } else { 'h = ${target} J/kg' }
	return no_answer('IF97 region 3 at ${p} Pa, ${what}: no state, last ${rho} kg/m3 and ${t} K')
}

// saturated is the two saturated phases at (p, T) on the saturation line:
// regions 1 and 2 up to 623.15 K, region 3's two densities at p above it,
// and its critical point at Tc (R7-97: Eq. (28) reproduces it exactly).
fn saturated(p f64, t f64) !(State, State) {
	if t <= t13 {
		return region1(p, t), region2(p, t)
	}
	if t >= tc || p >= pc {
		c := region3(rhoc, tc)
		return c, c
	}
	return at_p(region3(rho3(p, t, .liquid)!, t), p), at_p(region3(rho3(p, t, .vapour)!, t), p)
}

// mixture is the two-phase state of quality x between saturated phases. A
// saturated phase (x of 0 or 1) keeps its cp, cv and w.
fn mixture(liq State, vap State, p f64, x f64) State {
	if x == 0.0 {
		return State{
			...liq
			region: 4
			p:      p
			x:      0.0
		}
	}
	if x == 1.0 {
		return State{
			...vap
			region: 4
			p:      p
			x:      1.0
		}
	}
	nan := math.nan()
	return State{
		region: 4
		t:      liq.t
		p:      p
		v:      liq.v + x * (vap.v - liq.v)
		h:      liq.h + x * (vap.h - liq.h)
		s:      liq.s + x * (vap.s - liq.s)
		u:      liq.u + x * (vap.u - liq.u)
		cp:     nan
		cv:     nan
		w:      nan
		drhodp: nan
		dvdt:   nan
		x:      x
	}
}

// state_tx is the saturated state at T with quality x, 273.15 K to Tc.
pub fn state_tx(t f64, x f64) !State {
	if !(t >= t_min && t <= tc) || !(x >= 0.0 && x <= 1.0) {
		return out_of_range('IF97 has saturated states from 273.15 K to ${tc} K with 0 <= x <= 1: not ${t} K, x = ${x}')
	}
	p := if t >= tc { pc } else { psat97(t) }
	liq, vap := saturated(p, t)!
	return mixture(liq, vap, p, x)
}

// state_px is the saturated state at p with quality x, psat(273.15 K) to pc.
pub fn state_px(p f64, x f64) !State {
	if !(p >= psat97(t_min) && p <= pc) || !(x >= 0.0 && x <= 1.0) {
		return out_of_range('IF97 has saturated states from ${psat97(t_min)} Pa to ${pc} Pa with 0 <= x <= 1: not ${p} Pa, x = ${x}')
	}
	t := if p >= pc { tc } else { tsat97(p) }
	liq, vap := saturated(p, t)!
	return mixture(liq, vap, p, x)
}

// state_ph is the state at (p, h).
pub fn state_ph(p f64, h f64) !State {
	return by_p(p, h, false)
}

// state_ps is the state at (p, s).
pub fn state_ps(p f64, s f64) !State {
	return by_p(p, s, true)
}

fn of(st State, by_s bool) f64 {
	return if by_s { st.s } else { st.h }
}

// by_p is (p, h) or (p, s): the region by comparing g with its values on
// the regions' boundaries at p, then the backward equation's answer
// finished on the basic equation.
fn by_p(p f64, g f64, by_s bool) !State {
	what := if by_s { 's = ${g} J/(kg K)' } else { 'h = ${g} J/kg' }
	if !(p > 0.0 && p <= p_max) {
		return out_of_range(range_text(p, what))
	}
	p_t0 := psat97(t_min)
	low := if p >= p_t0 { region1(p, t_min) } else { region2(p, t_min) }
	high := if p <= p5_max { region5(p, t_max) } else { region2(p, t5) }
	// A g on the range's edge, within its rounding there (1e-12 of R·T for
	// h, of R for s), is in the range: 273.15 K or 2273.15 K given back.
	if g < of(low, by_s) - edge(low, by_s) || g > of(high, by_s) + edge(high, by_s) {
		return out_of_range(range_text(p, what))
	}
	if p <= psat97(t13) {
		if p >= p_t0 {
			ts := tsat97(p)
			liq, vap := region1(p, ts), region2(p, ts)
			if g <= of(liq, by_s) {
				t0 := if by_s { t1_ps(p, g) } else { t1_ph(p, g) }
				return finish_t(1, p, t0, t_min, ts, g, by_s)
			}
			if g < of(vap, by_s) {
				return mixture(liq, vap, p, (g - of(liq, by_s)) / (of(vap, by_s) - of(liq, by_s)))
			}
			return vapour(p, g, ts, by_s)
		}
		return vapour(p, g, t_min, by_s)
	}
	if g <= of(region1(p, t13), by_s) {
		t0 := if by_s { t1_ps(p, g) } else { t1_ph(p, g) }
		return finish_t(1, p, t0, t_min, t13, g, by_s)
	}
	t_23 := t23(p)
	if g >= of(region2(p, t_23), by_s) {
		return vapour(p, g, t_23, by_s)
	}
	if p < pc {
		ts := tsat97(p)
		liq, vap := saturated(p, ts)!
		if g > of(liq, by_s) && g < of(vap, by_s) {
			return mixture(liq, vap, p, (g - of(liq, by_s)) / (of(vap, by_s) - of(liq, by_s)))
		}
	}
	t0 := if by_s { t3_ps(p, g) } else { t3_ph(p, g) }
	v0 := if by_s { v3_ps(p, g) } else { v3_ph(p, g) }
	return finish3(p, 1.0 / v0, t0, g, by_s)
}

fn edge(st State, by_s bool) f64 {
	return 1e-12 * if by_s { r } else { r * st.t }
}

// vapour is regions 2 and 5 at (p, h) or (p, s), split at 1073.15 K by
// region 2's value; region 2 from lo, its lower edge at p (Ts or T23).
// Region 5's search may go to 1 K below 1073.15 K: its h there is within
// Table 43's inconsistency (0.1 kJ/kg) of region 2's, so a g just above
// region 2's can lie below region 5's.
fn vapour(p f64, g f64, lo f64, by_s bool) !State {
	top := region2(p, t5)
	// Above 50 MPa there's no region 5: a g within the edge's rounding above
	// region 2's top (by_p lets it in) is region 2's at 1073.15 K, not region
	// 5's, which took it to 1072.9 K at 100 MPa (M5b bundle B found it).
	if g <= of(top, by_s) || p > p5_max {
		t0 := if by_s { t2_ps(p, g) } else { t2_ph(p, g) }
		return finish_t(2, p, t0, lo, t5, g, by_s)
	}
	// Region 5 has no backward equation: start where g lies between its
	// ends, linearly.
	bottom := region5(p, t5)
	end := region5(p, t_max)
	t0 := t5 + (t_max - t5) * (g - of(bottom, by_s)) / (of(end, by_s) - of(bottom, by_s))
	return finish_t(5, p, t0, t5 - 1.0, t_max, g, by_s)
}
