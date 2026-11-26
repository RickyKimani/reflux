// Internal test (module water): IAPWS-95 water through thermo's volume roots
// and saturation, and through every flash on main (M5a's step 2, criterion
// 3), each answer checked from scratch: the volume roots against a scan of
// the isotherm's two branches, saturation by its equal fugacities and both
// roots on their branches, and the states against the iapws package's and
// CoolProp's (testdata/m5_iapws95.json, tools/gen_m5_iapws95.py).
module water

import flash
import json2
import math
import os
import thermo

struct StateRef {
	ok    bool
	rho   f64 // mol/m³
	h     f64 // J/mol, on its own reference
	s     f64
	phase string // the iapws package's
	error string // CoolProp's, where it has no state
}

struct StateRow {
	t        f64 @[json: 'T']
	p        f64 @[json: 'P']
	iapws    StateRef
	coolprop StateRef
}

struct SatRef {
	ok    bool
	p     f64
	rho_l f64
	rho_v f64
	h_l   f64
	h_v   f64
	s_l   f64
	s_v   f64
	error string
}

struct SatRow {
	t        f64 @[json: 'T']
	iapws    SatRef
	coolprop SatRef
}

struct Fixture {
	states     []StateRow
	saturation []SatRow
}

fn fixture() !Fixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm5_iapws95.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m5_iapws95.py`.')
	}
	return json2.decode[Fixture](text)!
}

fn water_fluid() thermo.Fluid {
	return fluid(-241814.0, 188.724) or { panic(err) }
}

// Isotherm is a scan of P(ρ) at one T: 60,000 densities, log-spaced from
// 1e-9 ρc to the model's δ_max, and where each physical branch ends.
struct Isotherm {
	t   f64
	f   AtT // the model at T, one mole: P = RT(1/V − F_V)
	rho []f64
	p   []f64
	ev  int // the vapour's branch is rho[..ev+1], P rising from ρ → 0
	sl  int // the liquid's is rho[sl..], P rising to δ_max
}

fn (iso Isotherm) pressure(rho f64) f64 {
	return thermo.r_gas * iso.t * (rho - iso.f.at(1.0 / rho).f_v)
}

fn scan(t f64) Isotherm {
	n := 60000
	lo := math.log(1e-9 * critical_rho)
	hi := math.log(delta_max * critical_rho)
	mut iso := Isotherm{
		t: t
		f: AtT{
			tau:  tau_parts(critical_t / t)
			t:    t
			ntot: 1.0
		}
	}
	rho := []f64{len: n, init: math.exp(lo + (hi - lo) * f64(index) / f64(n - 1))}
	p := rho.map(iso.pressure(it))
	mut ev := 0
	for ev + 1 < n && p[ev + 1] > p[ev] {
		ev++
	}
	mut sl := n - 1
	for sl > 0 && p[sl - 1] < p[sl] {
		sl--
	}
	iso = Isotherm{
		...iso
		rho: rho
		p:   p
		ev:  ev
		sl:  sl
	}
	return iso
}

// bisect refines a root of P(ρ) = p between two scanned densities.
fn (iso Isotherm) bisect(p f64, a f64, b f64) f64 {
	mut lo, mut hi := a, b
	flo := iso.pressure(lo) - p
	for _ in 0 .. 200 {
		mid := 0.5 * (lo + hi)
		fm := iso.pressure(mid) - p
		if (fm < 0.0) == (flo < 0.0) {
			lo = mid
		} else {
			hi = mid
		}
		if hi - lo <= 4.0 * f64_eps * hi {
			break
		}
	}
	return 0.5 * (lo + hi)
}

const f64_eps = 2.220446049250313e-16

// branch_root is the root of P(ρ) = p on one branch, scanned from i to j,
// or 0 where the branch doesn't reach p.
fn (iso Isotherm) branch_root(p f64, i int, j int) f64 {
	if !(iso.p[i] <= p && p <= iso.p[j]) {
		return 0.0
	}
	for k in i .. j {
		if iso.p[k] <= p && p <= iso.p[k + 1] {
			return iso.bisect(p, iso.rho[k], iso.rho[k + 1])
		}
	}
	return 0.0
}

// gibbs is Gʳ/RT of one mole at (T, ρ) and p: F + Z − 1 − ln Z.
fn gibbs(t f64, rho f64, p f64) f64 {
	v := 1.0 / rho
	z := p * v / (thermo.r_gas * t)
	return iapws95().residual(t, v, [1.0]).f + z - 1.0 - math.log(z)
}

// stable_root is the branches' answer: of their roots, the one of lower
// Gibbs energy, or 0 where neither branch reaches p.
fn (iso Isotherm) stable_root(p f64) f64 {
	rv := iso.branch_root(p, 0, iso.ev)
	rl := iso.branch_root(p, iso.sl, iso.rho.len - 1)
	if rv == 0.0 {
		return rl
	}
	if rl == 0.0 || rl == rv {
		return rv
	}
	return if gibbs(iso.t, rl, p) < gibbs(iso.t, rv, p) { rl } else { rv }
}

// Table 8: the saturation states at 275, 450 and 625 K.
struct Row8 {
	t     f64
	p     string // MPa
	rho_l string // kg/m³
	rho_v string
	h_l   string // kJ/kg
	h_v   string
	s_l   string // kJ/(kg·K)
	s_v   string
}

const table8 = [
	Row8{275.0, '0.698451167e-3', '0.999887406e3', '0.550664919e-2', '0.775972202e1', '0.250428995e4', '0.283094670e-1', '0.910660121e1'},
	Row8{450.0, '0.932203564', '0.890341250e3', '0.481200360e1', '0.749161585e3', '0.277441078e4', '0.210865845e1', '0.660921221e1'},
	Row8{625.0, '0.169082693e2', '0.567090385e3', '0.118290280e3', '0.168626976e4', '0.255071625e4', '0.380194683e1', '0.518506121e1'},
]

// half_unit and printed_ok are iapws95_test.v's: each test file is its own
// program, built with the module's sources but not with the other tests.
fn half_unit(printed string) f64 {
	mant := printed.all_before('e')
	exp := if printed.contains('e') { printed.all_after('e').int() } else { 0 }
	return 0.5 * math.pow(10.0, f64(exp - mant.all_after('.').len))
}

fn printed_ok(got f64, printed string) bool {
	return math.abs(got - printed.f64()) <= half_unit(printed)
}

fn per_kg_k(x f64) f64 {
	return x / molar_mass * 1e-3
}

fn test_table_8() {
	// The saturation line through the flash's bubble point (tvf at β = 0,
	// Newton on ln φ to 1e-12), h and s through Fluid.props on IAPWS's
	// reference by D6's constants (reference_offsets). thermo's own
	// saturation_pressure stops at 1e-10 in ln φ, 1.1e-10 in P at 625 K,
	// which puts ρ'' there at 118.2902805018 kg/m³, a hair past the half
	// unit of the table's 118.290280 (the iapws package and CoolProp give
	// …804506 and …804512).
	fl := water_fluid()
	m := fl.model
	dh, ds := reference_offsets(fl)
	for row in table8 {
		r := flash.tvf(m, row.t, 0.0, [1.0])!
		ps := r.p
		v := r.phases[0].state
		l := r.phases[1].state
		pl := fl.props(l)
		pv := fl.props(v)
		assert printed_ok(ps * 1e-6, row.p), 'p at ${row.t} K: ${ps * 1e-6} against ${row.p}'
		assert printed_ok(molar_mass / l.v, row.rho_l), "ρ' at ${row.t} K: ${molar_mass / l.v}"
		assert printed_ok(molar_mass / v.v, row.rho_v), "ρ'' at ${row.t} K: ${molar_mass / v.v}"
		assert printed_ok(per_kg_k(pl.h - dh), row.h_l), "h' at ${row.t} K: ${per_kg_k(pl.h - dh)}"
		assert printed_ok(per_kg_k(pv.h - dh), row.h_v), "h'' at ${row.t} K: ${per_kg_k(pv.h - dh)}"
		assert printed_ok(per_kg_k(pl.s - ds), row.s_l), "s' at ${row.t} K: ${per_kg_k(pl.s - ds)}"
		assert printed_ok(per_kg_k(pv.s - ds), row.s_v), "s'' at ${row.t} K: ${per_kg_k(pv.s - ds)}"
	}
}

fn test_saturation_up_to_a_tenth_of_a_millikelvin_below_tc() {
	// D7: the bracketed start and the branches' split. Each answer is checked
	// from scratch: its liquid and vapour roots are the scan's two branches'
	// roots, and their fugacities are equal; and it's held to the iapws
	// package's and CoolProp's. Measured: Psat within 5.7e-9 of both, the
	// worst near Tc, where each source stops at its own tolerance.
	fx := fixture()!
	m := iapws95()
	mut worst_iapws := 0.0
	mut worst_cool := 0.0
	for row in fx.saturation {
		ps := thermo.saturation_pressure(m, row.t) or {
			assert false, 'Psat at ${row.t} K: ${err}'
			continue
		}
		l := thermo.solve_tp(m, row.t, ps, [1.0], .liquid)!
		v := thermo.solve_tp(m, row.t, ps, [1.0], .vapor)!
		assert math.abs(l.ln_phi[0] - v.ln_phi[0]) <= 1e-9, 'ln φ at ${row.t} K'
		iso := scan(row.t)
		rl := iso.branch_root(ps, iso.sl, iso.rho.len - 1)
		rv := iso.branch_root(ps, 0, iso.ev)
		assert math.abs(1.0 / l.v / rl - 1.0) < 1e-9, "ρ' at ${row.t} K: ${1.0 / l.v} against the scan's ${rl}"
		assert math.abs(1.0 / v.v / rv - 1.0) < 1e-9, "ρ'' at ${row.t} K: ${1.0 / v.v} against the scan's ${rv}"
		worst_iapws = math.max(worst_iapws, math.abs(ps / row.iapws.p - 1.0))
		assert row.coolprop.ok, 'CoolProp has no saturation at ${row.t} K: ${row.coolprop.error}'
		worst_cool = math.max(worst_cool, math.abs(ps / row.coolprop.p - 1.0))
	}
	// And at 647.0959 K, 0.1 mK below Tc, which the fixture doesn't reach.
	thermo.saturation_pressure(m, 647.0959)!
	println('Psat against iapws ${worst_iapws:.2e}, CoolProp ${worst_cool:.2e}')
	assert worst_iapws < 1e-8
	assert worst_cool < 1e-8
}

fn test_acentric_factor_is_the_model_s_own() {
	ps := thermo.saturation_pressure(iapws95(), 0.7 * critical_t)!
	omega := -1.0 - math.log10(ps / critical_p)
	assert math.abs(omega - acentric) < 5e-5, 'ω = ${omega}'
}

fn test_the_volume_roots_are_the_branches() {
	// Criterion 3's states: thermo.solve_tp's stable root against the scan's
	// branches, and the iapws package's and CoolProp's densities. Measured:
	// the scan's to 1.2e-14, CoolProp's to 3.3e-12 and the package's to
	// 6.4e-13, 51 states skipped (below).
	fx := fixture()!
	m := iapws95()
	mut t_last := 0.0
	mut iso := Isotherm{}
	mut worst_scan := 0.0
	mut worst_iapws := 0.0
	mut worst_cool := 0.0
	mut refused := 0
	mut skipped := 0
	for row in fx.states {
		if row.t != t_last {
			iso = scan(row.t)
			t_last = row.t
		}
		want := iso.stable_root(row.p)
		st := thermo.solve_tp(m, row.t, row.p, [1.0], .stable) or {
			assert want == 0.0, '${row.t} K, ${row.p} Pa: ${err}, the scan has ${want}'
			refused++
			continue
		}
		rho := 1.0 / st.v
		assert want != 0.0, '${row.t} K, ${row.p} Pa: a root at ${rho}, none in the scan'
		worst_scan = math.max(worst_scan, math.abs(rho / want - 1.0))
		// CoolProp refuses a state at its Tmin, 273.16 K, below the triple
		// point's pressure; the rest it has.
		if !row.coolprop.ok {
			skipped++
			continue
		}
		worst_cool = math.max(worst_cool, math.abs(rho / row.coolprop.rho - 1.0))
		// The iapws package's own density solver goes wrong at low pressure
		// near Tc: from 640 to 647 K and 1 Pa to a few hundred, it returns a
		// density near the saturated vapour's (9,840 to 15,930 mol/m³ against
		// 2e-4 to 0.1), where CoolProp agrees with this to 1e-14. It's held
		// only where it agrees with CoolProp.
		if row.iapws.ok && math.abs(row.iapws.rho / row.coolprop.rho - 1.0) < 1e-6 {
			worst_iapws = math.max(worst_iapws, math.abs(rho / row.iapws.rho - 1.0))
		} else {
			skipped++
		}
	}
	println('roots: scan ${worst_scan:.2e}, iapws ${worst_iapws:.2e} (${skipped} skipped), CoolProp ${worst_cool:.2e}; ${refused} refused')
	assert worst_scan < 1e-12
	assert worst_cool < 1e-10
	assert worst_iapws < 1e-8
	assert skipped < 60
}

fn test_d7_refusals_and_the_compressed_liquid() {
	m := iapws95()
	// Above P(T, vmin), 1.38 GPa at 300 K, there's no root, not one at vmin.
	if st := thermo.solve_tp(m, 300.0, 2e9, [1.0], .stable) {
		assert false, 'a root at ${molar_mass / st.v} kg/m³ above P(T, vmin)'
	}
	// Liquid water at 300 K and 42 MPa, whose ideal-gas start lies in the
	// dome: step 0's solvers, with no split, took it to a spurious root at
	// 322 kg/m³.
	st := thermo.solve_tp(m, 300.0, 42.169650342858225e6, [1.0], .stable)!
	assert math.abs(molar_mass / st.v - 1014.6453001059319) < 1e-6, '${molar_mass / st.v} kg/m³'
}

// mean is a flash result's molar h, s, u and v, its phases' by fraction.
fn mean(fl thermo.Fluid, r flash.FlashResult) (f64, f64, f64, f64) {
	mut h, mut s, mut u, mut v := 0.0, 0.0, 0.0, 0.0
	for ph in r.phases {
		pr := fl.props(ph.state)
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
		u += ph.fraction * pr.u
		v += ph.fraction * ph.state.v
	}
	return h, s, u, v
}

// by_kind is the flash a kind names, at T or P and h, s, u or v. (An `or`
// can't follow a match whose arms return Results, so each arm forwards with
// `!` and the caller's `or` takes the function's.)
fn by_kind(kind string, fl thermo.Fluid, t f64, p f64, h f64, s f64, u f64, v f64) !flash.FlashResult {
	return match kind {
		'ph' { flash.ph(fl, p, h, [1.0])! }
		'ps' { flash.ps(fl, p, s, [1.0])! }
		'tv' { flash.tv(fl.model, t, v, [1.0])! }
		else { flash.uv(fl, u, v, [1.0])! }
	}
}

fn test_every_flash_from_each_state() {
	// Criterion 3: at each of the fixture's states, the PT flash, and PH,
	// PS, TV and UV back to it from its h, s, u and v. Each must land on the
	// PT flash's T and P, the stable root's state. Iterations are summed,
	// numbers that can't grow. Measured: every sixth state, 620 flashes, T to
	// 2.4e-11 and P to 3.1e-8 (TV's P, by a liquid's bulk modulus), 3,343
	// iterations; all of them, 3,700 flashes and 20,094 iterations.
	// Every flash costs about 0.06 s under tcc (step 10 is its speed), so the
	// routine run takes every sixth state; REFLUX_WATER_SWEEP=all takes all.
	fx := fixture()!
	fl := water_fluid()
	m := fl.model
	every := if os.getenv('REFLUX_WATER_SWEEP') == 'all' { 1 } else { 6 }
	mut worst_t := 0.0
	mut worst_p := 0.0
	mut iterations := 0
	mut flashes := 0
	for k, row in fx.states {
		if k % every != 0 {
			continue
		}
		r := flash.pt(m, row.t, row.p, [1.0]) or {
			// Above P(T, vmin), or a state the formulation has no root for:
			// solve_tp refuses it too (test_the_volume_roots_are_the_branches).
			thermo.solve_tp(m, row.t, row.p, [1.0], .stable) or { continue }
			assert false, 'PT at ${row.t} K, ${row.p} Pa: ${err}'
			continue
		}
		assert r.phases.len == 1, 'PT at ${row.t} K, ${row.p} Pa: ${r.phases.len} phases'
		h, s, u, v := mean(fl, r)
		for kind in ['ph', 'ps', 'tv', 'uv'] {
			back := by_kind(kind, fl, row.t, row.p, h, s, u, v) or {
				assert false, '${kind} back to ${row.t} K, ${row.p} Pa: ${err}'
				continue
			}
			worst_t = math.max(worst_t, math.abs(back.t / row.t - 1.0))
			worst_p = math.max(worst_p, math.abs(back.p / row.p - 1.0))
			iterations += back.iterations
			flashes++
		}
	}
	println('${flashes} flashes back: T ${worst_t:.2e}, P ${worst_p:.2e}; ${iterations} iterations')
	assert worst_t < 1e-10
	assert worst_p < 1e-7
	assert iterations <= if every == 1 { 20094 } else { 3343 }, '${iterations} iterations'
}

fn test_inside_the_dome() {
	// Two phases of water at T, by each flash that can reach them: TV and UV
	// at a quality's volume (where A is minimised over the phases' volumes
	// with no volume roots, the risk being the dome's spurious loops), PH
	// and PS at Psat, and the vapour-fraction flashes. Each must give Psat
	// (or Tsat), and the quality. Measured: T to 7.3e-13, P to 5.7e-12 and
	// β to 6.8e-11, each against the bubble point's own saturation.
	fx := fixture()!
	fl := water_fluid()
	m := fl.model
	mut worst_t := 0.0
	mut worst_p := 0.0
	mut worst_beta := 0.0
	for row in fx.saturation {
		if row.t > 646.0 {
			continue // the dome narrows to nothing; Psat itself is tested above
		}
		sat := flash.tvf(m, row.t, 0.0, [1.0])!
		ps := sat.p
		vap := sat.phases[0].state
		l := sat.phases[1].state
		pl := fl.props(l)
		pv := fl.props(vap)
		for x in [1e-4, 0.3, 0.9] {
			v := x * vap.v + (1.0 - x) * l.v
			h := x * pv.h + (1.0 - x) * pl.h
			s := x * pv.s + (1.0 - x) * pl.s
			u := x * pv.u + (1.0 - x) * pl.u
			for kind in ['tv', 'uv', 'ph', 'ps'] {
				r := by_kind(kind, fl, row.t, ps, h, s, u, v) or {
					assert false, '${kind} at ${row.t} K, quality ${x}: ${err}'
					continue
				}
				assert r.phases.len == 2, '${kind} at ${row.t} K, quality ${x}: ${r.phases.len} phases'
				worst_t = math.max(worst_t, math.abs(r.t / row.t - 1.0))
				worst_p = math.max(worst_p, math.abs(r.p / ps - 1.0))
				worst_beta = math.max(worst_beta, math.abs(r.beta - x))
			}
		}
		dew := flash.pvf(m, ps, 1.0, [1.0])!
		worst_t = math.max(worst_t, math.abs(dew.t / row.t - 1.0))
	}
	println('inside the dome: T ${worst_t:.2e}, P ${worst_p:.2e}, β ${worst_beta:.2e}')
	assert worst_t < 1e-10
	assert worst_p < 1e-9
	assert worst_beta < 1e-9
}

fn test_labels_are_the_branches() {
	// Below Tc a single phase is labelled by its branch (D7, as amended in
	// step 2): the phase identification parameter calls liquid water vapour
	// near its density maximum, Π = −7 at 300 K and 1 bar, −130 at 280 K.
	// Above Tc there's no split, and Π decides as for any model.
	m := iapws95()
	cases := [
		[280.0, 1e5, 0.0],
		[300.0, 1e5, 0.0],
		[273.16, 5.62e8, 0.0],
		[640.0, 2.5e7, 0.0],
		[400.0, 1e5, 1.0],
		[640.0, 1e7, 1.0],
	]
	for c in cases {
		st := thermo.solve_tp(m, c[0], c[1], [1.0], .stable)!
		want := if c[2] == 0.0 { thermo.PhaseLabel.liquid } else { .vapor }
		assert thermo.phase_label(m, st) == want, '${c[0]} K, ${c[1]} Pa: ${thermo.phase_label(m, st)}'
		r := flash.pt(m, c[0], c[1], [1.0])!
		assert r.phases[0].label == want
	}
}

fn test_no_water_below_251_k() {
	// Below 251.165 K, the melting curve's lowest point, there's no state
	// (step 2): not the liquid's extrapolation, whose H falls as T rises, nor
	// a vapour above the sublimation curve.
	m := iapws95()
	for p in [10.0, 1e5, 5.62e8] {
		if st := thermo.solve_tp(m, 240.0, p, [1.0], .stable) {
			assert false, 'a state at 240 K and ${p} Pa: ${molar_mass / st.v} kg/m³'
		}
	}
	// And the searches that step past it come back inside: PH to liquid at
	// 273.16 K and 562 MPa stepped to 215.7 K, where H was above its value
	// at 273.16 K, and UV at quality 0.9 found a vapour at 233.8 K with the
	// right U and V.
	fl := water_fluid()
	p := 5.62341325190349e8
	r := flash.pt(m, 273.16, p, [1.0])!
	h, _, _, _ := mean(fl, r)
	back := flash.ph(fl, p, h, [1.0])!
	assert math.abs(back.t / 273.16 - 1.0) < 1e-10, 'PH back: ${back.t} K'
}
