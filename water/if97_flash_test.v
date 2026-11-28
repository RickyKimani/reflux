// Internal test (module water): IF97 as a phase model (M5b steps 5 and 6,
// D12 to D14, criterion 7). Its phases are IF97's own states, by
// water.if97's functions, to rounding, on reflux's reference by one pair of
// constants; ln φ's derivatives are its H and V's; the flash's paths for
// one compound take region 4's line; and every flash agrees with the
// steam-table functions over IF97's range.
module water

import flash
import math
import thermo
import water.if97

const dhf_water = -241814.0 // J/mol, water's record (ChemSep)

const s0_water = 188.724 // J/(mol·K)

fn fluid97() thermo.Fluid {
	return if97_fluid(dhf_water, s0_water) or { panic(err) }
}

// offsets97 is the pair of constants between the fluid's h and s per mole
// and IF97's per kg times M, taken at one region-1 state (D14, as M5a's D6).
fn offsets97(fl thermo.Fluid) (f64, f64) {
	t, p := 500.0, 10e6
	pr := fl.props(fl.model.phase(t, p, [1.0], .stable) or { panic(err) })
	st := if97.state_pt(p, t) or { panic(err) }
	return pr.h - molar_mass * st.h, pr.s - molar_mass * st.s
}

// The grid: every region, both sides of B23 and of the saturation line,
// the critical region, and IF97's edges.
fn grid97() [][]f64 {
	mut out := [][]f64{}
	for t in [273.16, 300.0, 373.15, 450.0, 550.0, 623.15, 630.0, 640.0, 647.0, 650.0, 700.0, 800.0,
		863.0, 1000.0, 1073.15, 1200.0, 1500.0, 2000.0, 2273.15] {
		for p in [1e3, 1e5, 1e6, 1e7, 2e7, 2.5e7, 5e7, 1e8] {
			if t > if97.t5 && p > if97.p5_max {
				continue
			}
			out << [t, p]
		}
	}
	return out
}

// Each phase at (T, P) is IF97's state there, through the fluid: v to the
// bit (it's M·v), h and s per mole M·h and M·s plus the same two constants
// at every state, within 2e-8 J/mol and 2e-11 J/(mol·K) (the ideal gas's
// polynomial integrated from 298.15 K, its T⁶ terms reaching 7e5 J/mol at
// 2273 K, against IF97's own γ°; measured below), and cp, cv and w within
// 1e-12 of IF97's. The label is the region's (D12).
fn test_phases_are_if97_states() {
	fl := fluid97()
	m := if97_phases()
	h0, s0 := offsets97(fl)
	mut wh, mut ws, mut wc := 0.0, 0.0, 0.0
	mut count := 0
	for tp in grid97() {
		t, p := tp[0], tp[1]
		want := if97.state_pt(p, t) or { continue }
		st := m.phase(t, p, [1.0], .stable)!
		assert st.v == molar_mass * want.v, '${t} K, ${p} Pa'
		// Z is PV/nRT with reflux's R, as every model's (a planted mistake with
		// IF97's R survived until this).
		assert st.z == p * st.v / (thermo.r_gas * t)
		pr := fl.props(st)
		eh := math.abs(pr.h - molar_mass * want.h - h0)
		es := math.abs(pr.s - molar_mass * want.s - s0)
		wh = math.max(wh, eh)
		ws = math.max(ws, es)
		assert eh <= 2e-8, '${t} K, ${p} Pa: h off by ${eh} J/mol'
		assert es <= 2e-11, '${t} K, ${p} Pa: s off by ${es} J/(mol K)'
		for pair in [[pr.cp, molar_mass * want.cp], [pr.cv, molar_mass * want.cv], [pr.w, want.w]] {
			e := math.abs(pair[0] / pair[1] - 1)
			wc = math.max(wc, e)
			assert e <= 1e-12, '${t} K, ${p} Pa: ${pair[0]} vs ${pair[1]}'
		}
		label := m.phase_label(st)
		want_label := match want.region {
			1 { thermo.PhaseLabel.liquid }
			3 {
				if want.rho() > if97.rhoc { thermo.PhaseLabel.liquid } else { .vapor }
			}
			else { .vapor }
		}
		assert label == want_label, '${t} K, ${p} Pa: region ${want.region}, ${label}'
		count++
	}
	assert count > 120
	println('${count} IF97 states: h within ${wh} J/mol, s ${ws} J/(mol K), cp, cv, w ${wc}')
}

// M1's identities on its phases: ln φ's T- and P-derivatives against
// central differences of ln φ itself (Gibbs–Helmholtz, and Vʳ/RT), within
// 1e-6 of the derivative plus what the differences' rounding allows, ln φ's
// terms' 1e-14 over twice the step (at 450 K and 1 kPa ∂ln φ/∂P is 5.6e-8
// per Pa, and a difference over 0.02 Pa holds 1e-13 of it); n·∂ln φ/∂n is
// 0 for one compound.
fn test_ln_phi_derivatives() {
	m := if97_phases()
	mut worst := 0.0
	for tp in grid97() {
		t, p := tp[0], tp[1]
		if t <= 273.2 || t >= 2273.0 || (t >= 630 && t < 660 && p > 1.5e7 && p < 3e7) {
			continue // a difference needs room, and none in the critical region
		}
		st := m.phase(t, p, [1.0], .stable) or { continue }
		d := m.ln_phi_derivs(st)!
		assert d.dn == [0.0]
		ht := 1e-4 * t
		hp := 1e-5 * p
		a := m.phase(t + ht, p, [1.0], .stable) or { continue }
		b := m.phase(t - ht, p, [1.0], .stable) or { continue }
		c := m.phase(t, p + hp, [1.0], .stable) or { continue }
		e := m.phase(t, p - hp, [1.0], .stable) or { continue }
		region := m.state_of_phase(st).region
		if [a, b, c, e].any(m.state_of_phase(it).region != region) {
			continue // a difference across a region's boundary
		}
		dt := (a.ln_phi[0] - b.ln_phi[0]) / (2 * ht)
		dp := (c.ln_phi[0] - e.ln_phi[0]) / (2 * hp)
		et := math.abs(dt - d.dt[0]) / (1e-6 * math.abs(d.dt[0]) + 1e-13 / ht)
		ep := math.abs(dp - d.dp[0]) / (1e-6 * math.abs(d.dp[0]) + 1e-13 / hp)
		worst = math.max(worst, math.max(et, ep))
		assert et <= 1, '${t} K, ${p} Pa: dlnφ/dT ${d.dt[0]} vs ${dt}'
		assert ep <= 1, '${t} K, ${p} Pa: dlnφ/dP ${d.dp[0]} vs ${dp}'
	}
	println('ln φ derivatives within ${worst} of their bound from differences')
}

// The volume's derivatives a phase's residual gives, (∂V/∂T)_P from each
// region's own derivatives (if97.State.dvdt) and (∂P/∂V)_T, against central
// differences of the phase's own volume, within 1e-6 plus the differences'
// rounding (10 ulps of V over twice the step), each kept in its region.
fn test_volume_derivatives() {
	m := if97_phases()
	mut worst := 0.0
	for tp in grid97() {
		t, p := tp[0], tp[1]
		if t <= 273.2 || t >= 2273.0 || (t >= 630 && t < 660 && p > 1.5e7 && p < 3e7) {
			continue
		}
		st := m.phase(t, p, [1.0], .stable) or { continue }
		res := m.phase_residual(st)
		ht := 1e-4 * t
		hp := 1e-5 * p
		a := m.phase(t + ht, p, [1.0], .stable) or { continue }
		b := m.phase(t - ht, p, [1.0], .stable) or { continue }
		c := m.phase(t, p + hp, [1.0], .stable) or { continue }
		e := m.phase(t, p - hp, [1.0], .stable) or { continue }
		region := m.state_of_phase(st).region
		if [a, b, c, e].any(m.state_of_phase(it).region != region) {
			continue
		}
		dvdt := (a.v - b.v) / (2 * ht)
		dvdp := (c.v - e.v) / (2 * hp)
		et := math.abs(dvdt - res.dv_dt) / (1e-6 * math.abs(res.dv_dt) + 10 * 2.2e-16 * st.v / ht)
		ep := math.abs(1.0 / dvdp - res.dp_dv) / (1e-6 * math.abs(res.dp_dv) + 10 * 2.2e-16 * st.v / hp * res.dp_dv * res.dp_dv)
		worst = math.max(worst, math.max(et, ep))
		assert et <= 1, '${t} K, ${p} Pa: (∂V/∂T)_P ${res.dv_dt} vs ${dvdt}'
		assert ep <= 1, '${t} K, ${p} Pa: (∂P/∂V)_T ${res.dp_dv} vs ${1.0 / dvdp}'
	}
	println('volume derivatives within ${worst} of their bound from differences')
}

// D13: the vapour-fraction flashes on IF97 water take region 4's line, T
// and P each other's Pˢᵃᵗ and Tˢᵃᵗ exactly, the two phases IF97's saturated
// states there, their ln φ apart by IF97's inconsistency alone, within the
// Prague value; at and above the critical point, no line.
fn test_saturation_is_region_4() {
	fl := fluid97()
	h0, _ := offsets97(fl)
	mut gap := 0.0
	for t in [273.16, 300.0, 350.0, 450.0, 550.0, 600.0, 623.15, 630.0, 640.0, 646.0, 647.0] {
		for beta in [0.0, 0.5, 1.0] {
			r := flash.tvf(fl.model, t, beta, [1.0])!
			assert r.path == .own_line
			assert r.t == t && r.p == if97.psat97(t), '${t} K'
			liq := r.phases[1].state
			vap := r.phases[0].state
			gap = math.max(gap, math.abs(liq.ln_phi[0] - vap.ln_phi[0]))
			wl := if97.state_tx(t, 0.0)!
			wv := if97.state_tx(t, 1.0)!
			assert math.abs(fl.props(liq).h - molar_mass * wl.h - h0) <= 2e-8
			assert math.abs(fl.props(vap).h - molar_mass * wv.h - h0) <= 2e-8
			back := flash.pvf(fl.model, r.p, beta, [1.0])!
			assert back.p == r.p && back.t == if97.tsat97(r.p)
			assert math.abs(back.t - t) <= 1e-12 * t
		}
	}
	assert gap <= 5e-4
	println('on region 4, the two phases ln φ within ${gap}')
	for t in [if97.tc, 650.0] {
		if r := flash.tvf(fl.model, t, 0.0, [1.0]) {
			assert false, 'a saturation at ${t} K: ${r.p} Pa'
		} else {
			assert err.msg().contains("IF97's saturation line runs from 273.15 K to its critical point"), err.msg()
		}
	}
}

// on_boundary says (T, P) is on B13 or at 1073.15 K, where IF97's regions
// meet within Table 43's inconsistency.
fn on_boundary(t f64, p f64) bool {
	return (t == if97.t13 && p > if97.p23(t)) || t == if97.t5
}

// Criterion 7: every flash on IF97 water agrees with the steam-table
// functions over its range: a PT flash is the region map's state exactly
// (one compound's PT is its stable phase); PH and PS back from it give
// IF97's own T within 3e-11·T (the flash's energy balance at 1e-10·RT, a
// gas's 2e-8 K at 863 K; IF97's finished backward equations to rounding),
// as M5a held IAPWS-95's; inside the dome a PH
// or PS gives region 4's T and IF97's quality within 1e-10.
//
// But on two region boundaries one h (or s) is had at two T, each an IF97
// state with the value asked, within R7-97's Table 43's inconsistency
// (0.2 kJ/kg in h): on B13, 623.15 K above 16.53 MPa, regions 1 and 3 up
// to 2.1 mK apart, where the flash's search can take region 3's; and at
// 1073.15 K, regions 2 and 5 up to 38 mK apart (110 J/kg at 10 MPa), where
// IF97's own search may take region 5's up to 1 K below 1073.15 K. From a
// state on either, IF97's own backward functions take either side by
// rounding. There the flash's state is held to the value asked, the two
// within Table 43's 0.2 kJ/kg over steam's cp of 2 kJ/(kg·K) or more,
// 0.1 K, and the states it happens at are counted. The PH and PS flashes'
// iterations are recorded, a number that can't grow: 2,118.
fn test_every_flash_against_the_steam_tables() {
	fl := fluid97()
	m97 := if97_phases()
	h0, s0 := offsets97(fl)
	mut wt, mut wb := 0.0, 0.0
	mut count, mut across, mut iterations := 0, 0, 0
	for tp in grid97() {
		t, p := tp[0], tp[1]
		want := if97.state_pt(p, t) or { continue }
		r := flash.pt(fl.model, t, p, [1.0])!
		assert r.phases.len == 1 && r.phases[0].state.v == molar_mass * want.v
		pr := fl.props(r.phases[0].state)
		for g in ['h', 's'] {
			a := if g == 'h' {
				flash.ph(fl, p, pr.h, [1.0]) or { panic('PH at ${t} K, ${p} Pa: ${err}') }
			} else {
				flash.ps(fl, p, pr.s, [1.0]) or { panic('PS at ${t} K, ${p} Pa: ${err}') }
			}
			st := if g == 'h' {
				if97.state_ph(p, (pr.h - h0) / molar_mass)!
			} else {
				if97.state_ps(p, (pr.s - s0) / molar_mass)!
			}
			iterations += a.iterations
			e := math.abs(a.t - st.t)
			if on_boundary(t, p) && e > 3e-11 * t {
				got := fl.props(a.phases[0].state)
				off := if g == 'h' { got.h - pr.h } else { (got.s - pr.s) * t }
				assert math.abs(off) <= 1e-9 * thermo.r_gas * t, '${g} at ${t} K, ${p} Pa: off by ${off}'
				assert m97.state_of_phase(a.phases[0].state).region != st.region || math.abs(a.t - t) > 3e-11 * t
				assert e <= 0.1, '${g} at ${t} K, ${p} Pa: ${a.t} vs IF97 ${st.t}'
				across++
				wb = math.max(wb, e)
				continue
			}
			wt = math.max(wt, e)
			assert e <= 3e-11 * t, '${g} back at ${t} K, ${p} Pa: ${a.t} vs IF97 ${st.t}'
			if on_boundary(t, p) {
				continue // IF97's own may land across it too
			}
			assert math.abs(a.t - t) <= 3e-11 * t
		}
		count++
	}
	for p in [1e3, 1e5, 1e6, 1e7, 2e7] {
		for x in [0.1, 0.5, 0.9] {
			want := if97.state_px(p, x)!
			h := molar_mass * want.h + h0
			a := flash.ph(fl, p, h, [1.0])!
			assert a.t == if97.tsat97(p), '${p} Pa, x ${x}'
			assert math.abs(a.beta - x) <= 1e-10, '${p} Pa: β ${a.beta} vs x ${x}'
		}
	}
	println('${count} states: PH and PS back within ${wt} K of IF97; ${across} across a boundary, within ${wb} K; ${iterations} iterations')
	assert across <= 10
	assert iterations <= 2118, '${iterations} iterations, more than recorded'
}

// What IF97 has no state for is refused with its range, and a PH search
// that steps past 2273.15 K comes back into it (thermo.NoState, D12); the
// TV and UV flashes refuse IF97, which has no F(T, V); the liquid and vapour
// sides are refused where their equations have no state.
fn test_range_and_refusals() {
	fl := fluid97()
	m := if97_phases()
	for tp in [[270.0, 1e5], [2300.0, 1e5], [1200.0, 6e7], [500.0, 1.1e8]] {
		if st := m.phase(tp[0], tp[1], [1.0], .stable) {
			assert false, 'a state at ${tp}: ${st.v}'
		} else {
			assert err.msg().contains('IF97 covers 273.15 K to 1073.15 K up to 100 MPa'), err.msg()
		}
		assert m.no_state(tp[0], tp[1])
	}
	// Steam at 1 bar and 2200 K: the PH search from Tsat steps by 1.5 to
	// 2829 K, past IF97's range, and back.
	h := fl.props(m.phase(2200, 1e5, [1.0], .stable)!).h
	r := flash.ph(fl, 1e5, h, [1.0])!
	assert math.abs(r.t - 2200) <= 3e-11 * 2200
	v := m.phase(500, 1e6, [1.0], .stable)!.v
	if a := flash.tv(fl.helmholtz() or {
		assert err.msg().contains("an equation of state's Helmholtz energy"), err.msg()
		return
	}, 500, v, [1.0])
	{
		assert false, 'a TV flash: ${a.p}'
	}
}

// IF97 water's transport (D10): R12-08's and R15-11's industrial forms on
// the phase's IF97 state, and R1-76's σ for a liquid, to the bit.
fn test_transport_if97() {
	m := if97_phases()
	tm := TransportIf97{}
	for tp in [[300.0, 1e5], [450.0, 1e6], [640.0, 2.1e7], [700.0, 3e7], [1500.0, 1e6]] {
		st := m.phase(tp[0], tp[1], [1.0], .stable)!
		label := m.phase_label(st)
		got := tm.phase(label, tp[0], tp[1], [1.0], st.v)
		want := if97.state_pt(tp[1], tp[0])!
		assert got.viscosity.value == viscosity_industrial(tp[0], molar_mass / st.v)!
		assert got.conductivity.value == conductivity_industrial(m.state_at(tp[0], tp[1], molar_mass / st.v))!.lambda
		assert math.abs(got.conductivity.value / conductivity_industrial(want)!.lambda - 1) <= 1e-12
		if label == .vapor {
			assert got.surface_tension.why == 'a gas has no surface tension'
		} else {
			assert got.surface_tension.value == surface_tension(tp[0])!
		}
	}
}

// Skewed is IF97 with its Tˢᵃᵗ moved by 1e-9 of itself, or its gap bound
// tightened: what own_line must refuse, though IF97 itself never trips
// either check (a planted mistake loosening the first survived until this).
struct Skewed {
	If97Phases
	shift f64
	gap   f64
}

fn (s Skewed) saturation_t(p f64) !f64 {
	return s.If97Phases.saturation_t(p)! * (1 + s.shift)
}

fn (s Skewed) saturation_gap() f64 {
	return s.gap
}

// own_line checks the model's line by its own inverse, to 1e-12, and the
// phases' ln φ by the model's own bound, refusing either, saying why.
fn test_own_line_refuses_an_inconsistent_line() {
	off := Skewed{
		If97Phases: if97_phases()
		shift:      1e-9
		gap:        5e-4
	}
	if r := flash.tvf(off, 400, 0.0, [1.0]) {
		assert false, 'a saturation at ${r.p} Pa'
	} else {
		assert err.msg().contains("the model's saturation line gives back"), err.msg()
	}
	tight := Skewed{
		If97Phases: if97_phases()
		gap:        1e-6
	}
	if r := flash.tvf(tight, 350, 0.0, [1.0]) {
		assert false, 'a saturation at ${r.p} Pa'
	} else {
		assert err.msg().contains("on the model's own saturation line its phases' ln φ are"), err.msg()
	}
	assert flash.tvf(Skewed{
		If97Phases: if97_phases()
		gap:        5e-4
	}, 400, 0.0, [1.0])!.path == .own_line
}
