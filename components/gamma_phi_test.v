// The γ–φ liquid (M6 step 5, criteria 3 and 7) against thermo 0.6.1's
// GibbsExcessLiquid on the same basis, at the seeded states of
// testdata/m6_gamma_phi.json (tools/gen_m6_gamma_phi.py): ln φ and its
// derivatives, and H, S and Cp less the ideal gas's; and Henry's law for
// gases in water and in a mixed solvent.
import activity
import components
import flash
import json2
import math
import numeric
import os
import tdep
import thermo

struct Answer {
	ln_phi []f64
	dt     []f64
	dp     []f64
	dn     []f64
	h      f64
	s      f64
	cp     f64
	v      ?f64
}

struct PhiSat {
	thermo  [][]f64 // thermo's ln φˢᵃᵗ, at PR's own Pˢᵃᵗ, and its two T-derivatives
	at_psat [][]f64 @[json: 'at_psat'] // PR's ln φ at the correlation's Pˢᵃᵗ, likewise
}

struct LiquidCase {
	model   string
	basis   string
	cas     []string
	t       f64 @[json: 'T']
	p       f64 @[json: 'P']
	x       []f64
	thermo  Answer
	phi_sat ?PhiSat @[json: 'phi_sat']
}

// moved is thermo's answer with its φˢᵃᵗ, taken at PR's own vapour pressure,
// replaced by PR's φ at the correlation's, which reflux takes: only that
// makes a pure liquid's fugacity its vapour's at Pˢᵃᵗ. Each ln φᵢ and its
// T-derivatives move by the difference, and H, S and Cp as they're made of
// them, on thermo's R.
fn moved(c LiquidCase, r f64) Answer {
	ps := c.phi_sat or { return c.thermo }
	a := c.thermo
	mut lnphi, mut dt := a.ln_phi.clone(), a.dt.clone()
	mut s0, mut s1, mut s2 := 0.0, 0.0, 0.0
	for i, x in c.x {
		d0 := ps.at_psat[i][0] - ps.thermo[i][0]
		d1 := ps.at_psat[i][1] - ps.thermo[i][1]
		d2 := ps.at_psat[i][2] - ps.thermo[i][2]
		lnphi[i] += d0
		dt[i] += d1
		s0 += x * d0
		s1 += x * d1
		s2 += x * d2
	}
	t := c.t
	dh := -r * t * t * s1
	return Answer{
		...a
		ln_phi: lnphi
		dt:     dt
		h:      a.h + dh
		s:      a.s + (dh - r * t * s0) / t
		cp:     a.cp - r * t * (2.0 * s1 + t * s2)
	}
}

struct HenryCase {
	cas    []string
	t      f64 @[json: 'T']
	p      f64 @[json: 'P']
	x      []f64
	thermo Answer
}

struct Source {
	r_thermo f64 @[json: 'R_thermo']
}

struct HenryConstant {
	solute  int
	solvent int
	abcde   []f64
	tmin    f64
	tmax    f64
}

struct HenryState {
	t      f64 @[json: 'T']
	p      f64 @[json: 'P']
	x      []f64
	thermo Answer
}

struct IapwsGas {
	name    string
	cas     string
	abc     []f64
	tmin    f64
	tmax    f64
	t       []f64 @[json: 'T']
	ln_kh   []f64
	ln_psat []f64
}

struct IapwsFixture {
	gases []IapwsGas
}

struct HenryMixed {
	cas       []string
	constants []HenryConstant
	states    []HenryState
}

struct FlashAnswer {
	t  f64 @[json: 'T']
	p  f64 @[json: 'P']
	vf f64 @[json: 'VF']
	x  []f64
	y  []f64
}

struct FlashCase {
	model  string
	cas    []string
	kind   string
	z      []f64
	spec   []f64
	thermo FlashAnswer
}

struct Azeotrope {
	x1 f64
	t  f64 @[json: 'T']
}

struct AzeotropeCase {
	model      string
	cas        []string
	p          f64 @[json: 'P']
	azeotropes []Azeotrope
}

struct Fixture {
	source      Source
	liquid      []LiquidCase
	henry       []HenryCase
	flash       []FlashCase
	azeotropes  []AzeotropeCase
	henry_mixed HenryMixed @[json: 'henry_mixed']
}

fn fixture() !Fixture {
	return json2.decode[Fixture](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_gamma_phi.json'))!)!
}

fn liquid_of(name string) components.LiquidModel {
	return match name {
		'NRTL' { .nrtl }
		'UNIQUAC' { .uniquac }
		else { .unifac }
	}
}

// gap is |a − b| over max(|b|, floor).
fn gap(a f64, b f64, floor f64) f64 {
	return math.abs(a - b) / math.max(math.abs(b), floor)
}

struct Worst {
mut:
	m map[string]f64
}

fn (mut w Worst) note(key string, v f64) {
	if v > (w.m[key] or { 0.0 }) {
		w.m[key] = v
	}
}

// compare notes how far the phase at (t, p, x) is from thermo's answer.
// thermo's energies are on its R, N_A·k_B to every digit, and reflux's on
// thermo.r_gas, 1.84e-11 below: scaled.
fn compare(mut worst Worst, tag string, fl thermo.Fluid, t f64, p f64, x []f64, want Answer, r_thermo f64) ! {
	m := fl.model
	st := m.phase(t, p, x, .liquid)!
	assert m.phase_label(st) == .liquid
	d := m.ln_phi_derivs(st)!
	res := m.phase_residual(st)
	nc := x.len
	for i in 0 .. nc {
		worst.note('${tag} ln φ', gap(st.ln_phi[i], want.ln_phi[i], 1.0))
		worst.note('${tag} ∂T', gap(d.dt[i], want.dt[i], 1e-3))
		worst.note('${tag} ∂P', gap(d.dp[i], want.dp[i], 1e-9))
		for j in 0 .. nc {
			worst.note('${tag} dn', math.abs(d.dn[i * nc + j] - want.dn[i * nc + j]) / math.max(math.abs(want.dn[i * nc + j]),
				1e-3))
		}
	}
	scale := thermo.r_gas / r_thermo
	worst.note('${tag} H', gap(res.h, want.h * scale, 10.0))
	worst.note('${tag} S', gap(res.s, want.s * scale, 1.0))
	worst.note('${tag} Cp', gap(res.cp, want.cp * scale, 1.0))
	if v := want.v {
		worst.note('${tag} V', gap(st.v, v, 1e-6))
	}
}

// cp_against_h is the gap between the liquid's Cp and Richardson's
// extrapolation of central differences of its own H in T, at fixed P and x:
// a check of the exact second derivatives from what the first ones give.
fn cp_against_h(fl thermo.Fluid, t f64, p f64, x []f64) !f64 {
	m := fl.model
	cp := m.phase_residual(m.phase(t, p, x, .liquid)!).cp
	h := 1e-3 * t
	mut hs := []f64{len: 4}
	for k, dt in [h, -h, h / 2.0, -h / 2.0] {
		hs[k] = m.phase_residual(m.phase(t + dt, p, x, .liquid)!).h
	}
	fd := (4.0 * (hs[2] - hs[3]) / h - (hs[0] - hs[1]) / (2.0 * h)) / 3.0
	return gap(cp, fd, 1.0)
}

fn test_the_liquid_against_thermo() {
	f := fixture()!
	db := components.load_default()!
	mut worst := Worst{}
	mut counted := map[string]int{}
	for c in f.liquid {
		s := db.gamma_phi(c.cas,
			liquid:   liquid_of(c.model)
			poynting: c.basis != 'Psat'
			vapour:   if c.basis == 'Poynting&PhiSat' { ?thermo.CubicKind(.pr) } else { none }
			kij:      false
		)!
		compare(mut worst, c.basis, s.fluid, c.t, c.p, c.x, moved(c, f.source.r_thermo), f.source.r_thermo)!
		worst.note('${c.basis} Cp against ∂H/∂T', cp_against_h(s.fluid, c.t, c.p, c.x)!)
		counted['${c.model} ${c.basis}']++
	}
	println('states: ${counted}')
	mut keys := worst.m.keys()
	keys.sort()
	for k in keys {
		println('  ${k}: ${worst.m[k]:.2e}')
	}
	// Measured, under tcc. Without Poynting's factor, a few ulps of ln φ's
	// terms, and of dn's (∂ᵢⱼ − Σₖ xₖ∂ᵢₖ). With it, V(P − Pˢᵃᵗ)/RT is on
	// thermo's R, 1.84e-11 from reflux's, which moves ln φ and ∂P by that
	// share of the term: 4.7e-13 and 1.6e-12. Cp is the closed form's, held
	// also to differences of reflux's own H (3.2e-12 to 1.7e-11, their
	// truncation and rounding); with φˢᵃᵗ, both sides take its second
	// derivative by differences, which leaves 8e-9.
	limits := {
		'Psat ln φ':                        1e-13
		'Psat ∂T':                          1e-13
		'Psat dn':                          1e-12
		'Psat H':                           1e-13
		'Psat S':                           1e-13
		'Psat Cp':                          1e-12
		'Psat V':                           1e-14
		'Psat Cp against ∂H/∂T':            1e-10
		'Poynting ln φ':                    5e-12
		'Poynting ∂T':                      5e-12
		'Poynting ∂P':                      1e-11
		'Poynting dn':                      1e-12
		'Poynting H':                       5e-12
		'Poynting S':                       5e-12
		'Poynting Cp':                      5e-11
		'Poynting V':                       1e-14
		'Poynting Cp against ∂H/∂T':        1e-10
		'Poynting&PhiSat ln φ':             5e-12
		'Poynting&PhiSat ∂T':               5e-12
		'Poynting&PhiSat ∂P':               1e-11
		'Poynting&PhiSat dn':               1e-12
		'Poynting&PhiSat H':                5e-12
		'Poynting&PhiSat S':                5e-12
		'Poynting&PhiSat Cp':               1e-7
		'Poynting&PhiSat V':                1e-14
		'Poynting&PhiSat Cp against ∂H/∂T': 1e-9
	}
	for k, v in worst.m {
		assert v <= limits[k] or { 0.0 }, '${k}: ${v}'
	}
	assert worst.m.len == limits.len, '${worst.m.keys()}'
}

fn test_henrys_law_against_thermo() {
	f := fixture()!
	db := components.load_default()!
	mut worst := Worst{}
	for c in f.henry {
		s := db.gamma_phi(c.cas, henry: c.cas[1..], poynting: false)!
		assert s.solutes == c.cas[1..], '${s.solutes}'
		compare(mut worst, 'Henry', s.fluid, c.t, c.p, c.x, c.thermo, f.source.r_thermo)!
	}
	mut keys := worst.m.keys()
	keys.sort()
	for k in keys {
		println('  ${k}: ${worst.m[k]:.2e}')
	}
	// Measured: 2.5e-14 at most. Water alone is ideal, and each gas's ln H
	// is a sum of four terms.
	for k, v in worst.m {
		assert v <= 1e-12, '${k}: ${v}'
	}
}

// Henry's law in a mixed solvent (criterion 7): ethylene and propane in
// water and ethylene glycol, mixed as thermo's 'solvents_with_parameters' mode
// mixes them, Σ xⱼ ln Hᵢⱼ / Σ xⱼ. ChemSep's bank has water alone as a
// solvent, so each gas's constant in the glycol is the fixture's, synthetic,
// and water's is the bank's, the fixture's too. No bank lists water with the
// glycol, so the solvents are an ideal solution (ideal: true), as thermo's
// is here; its dn is by differences of its ln φ.
fn test_henrys_law_in_a_mixed_solvent() {
	f := fixture()!
	db := components.load_default()!
	hm := f.henry_mixed
	s := db.gamma_phi(hm.cas, henry: hm.cas[2..], ideal: true, poynting: false)!
	assert s.solutes == hm.cas[2..], '${s.solutes}'
	mut solutes := []thermo.Solute{}
	for sol in s.model.solutes {
		mut henrys := sol.henrys.clone()
		assert henrys.len == 1 && henrys[0].solvent == 0
		for c in hm.constants.filter(it.solute == sol.i) {
			k := c.abcde
			h := thermo.Henry{
				solvent: c.solvent
				a:       k[0]
				b:       k[1]
				c:       k[2]
				d:       k[3]
				e:       k[4]
				tmin:    c.tmin
				tmax:    c.tmax
			}
			if c.solvent == 0 {
				assert h == henrys[0], "the bank's constant in water, not the fixture's"
			} else {
				henrys << h
			}
		}
		solutes << thermo.Solute{sol.i, henrys}
	}
	model := thermo.new_gamma_phi(thermo.GammaPhi{
		...s.model
		solutes: solutes
	})!
	fl := thermo.new_phase_fluid(model, s.fluid.ideal)!
	mut worst := Worst{}
	for c in hm.states {
		compare(mut worst, 'mixed', fl, c.t, c.p, c.x, c.thermo, f.source.r_thermo)!
		// dn apart: a solvent's column, thermo's step 1e-3 mol, and a gas's,
		// its step half the gas's amount, where ln φ is flat (each depends on
		// the solvents' proportions alone) and reflux's is exactly 0.
		st := fl.model.phase(c.t, c.p, c.x, .liquid)!
		dn := fl.model.ln_phi_derivs(st)!.dn
		for i in 0 .. 4 {
			for j in 0 .. 4 {
				want := c.thermo.dn[i * 4 + j]
				if j < 2 {
					worst.note('mixed dn, a solvent column', gap(dn[i * 4 + j], want, 1e-3))
				} else {
					assert dn[i * 4 + j] == 0.0
					worst.note('mixed dn, a gas column, thermo absolute', math.abs(want))
				}
			}
		}
	}
	mut keys := worst.m.keys()
	keys.sort()
	for k in keys {
		println('  ${k}: ${worst.m[k]:.2e}')
	}
	// Measured: a few ulps of each ln φ's terms and of H, S and Cp, and ∂P
	// is −1/P on both sides, to the bit (so never noted). In dn,
	// a solvent's column to Richardson's truncation and rounding; a gas's,
	// thermo's differences' rounding alone, an ln φ near 8 over a step down
	// to 5e-6 mol, which against compare's floor of 1e-3 is 1.3e-7.
	limits := {
		'mixed ln φ':                              1e-13
		'mixed ∂T':                                1e-12
		'mixed dn':                                1e-6
		'mixed dn, a solvent column':              1e-9
		'mixed dn, a gas column, thermo absolute': 1e-9
		'mixed H':                                 1e-12
		'mixed S':                                 1e-12
		'mixed Cp':                                1e-11
	}
	for k, v in worst.m {
		assert v <= limits[k] or { 0.0 }, '${k}: ${v}'
	}
	assert worst.m.len == limits.len, '${worst.m.keys()}'
}

// A pure liquid's fugacity is its vapour's at Pˢᵃᵗ, by the ideal gas or by
// Peng–Robinson (φˢᵃᵗ at the correlation's Pˢᵃᵗ, which is why it's taken
// there); and its heat of vaporisation from H is Clapeyron's, T·ΔV·dPˢᵃᵗ/dT,
// on its own Pˢᵃᵗ and both phases' volumes (criterion 3).
fn test_a_pure_liquid_is_consistent() {
	db := components.load_default()!
	mut worst_f, mut worst_c := 0.0, 0.0
	for name in ['methanol', 'water', 'acetone'] {
		for vapour in [?thermo.CubicKind(none), ?thermo.CubicKind(.pr)] {
			s := db.gamma_phi([name], vapour: vapour)!
			m := s.fluid.model
			rec := s.records[0]
			vp := rec.property('vapor_pressure')!
			for t in [300.0, 330.0, 360.0] {
				ps := vp.at(t).value
				liq := m.phase(t, ps, [1.0], .liquid)!
				vap := m.phase(t, ps, [1.0], .vapor)!
				worst_f = math.max(worst_f, math.abs(liq.ln_phi[0] - vap.ln_phi[0]))
				dh := m.phase_residual(vap).h - m.phase_residual(liq).h
				clapeyron := t * (vap.v - liq.v) * vp.scale * vp.c.deriv(t)
				worst_c = math.max(worst_c, math.abs(dh / clapeyron - 1.0))
			}
		}
	}
	println('pure: ln f, the liquid against the vapour, at Pˢᵃᵗ ${worst_f:.1e}; ΔH against Clapeyron ${worst_c:.1e}')
	// Measured, 8.9e-16 and 4.4e-16: ln φ's few terms, and ΔH's, each a few
	// ulps of a quantity some 40 kJ/mol.
	assert worst_f <= 1e-13 && worst_c <= 1e-12
}

// A mixture's H less its pure liquids', each at its own Pˢᵃᵗ and Poynting
// factor, is the activity model's Hᴱ = −RT²∂(Gᴱ/RT)/∂T (criterion 3).
fn test_the_excess_enthalpy_is_the_models() {
	db := components.load_default()!
	cas := ['67-64-1', '67-56-1', '7732-18-5']
	x := [0.2, 0.3, 0.5]
	mut worst := 0.0
	for t in [300.0, 330.0] {
		s := db.gamma_phi(cas)!
		m := s.fluid.model
		mix := m.phase_residual(m.phase(t, 1e5, x, .liquid)!).h
		mut pure := 0.0
		for i, c in cas {
			one := db.gamma_phi([c])!.fluid.model
			pure += x[i] * one.phase_residual(one.phase(t, 1e5, [1.0], .liquid)!).h
		}
		g := s.model.liquid
		mut w := g.work()
		activity.eval(g, t, x, .t, mut w)!
		he := activity.excess(t, x, &w).h
		worst = math.max(worst, math.abs(mix - pure - he) / math.abs(he))
	}
	println('Hᴱ: ${worst:.1e}')
	// Measured, 2.2e-15: Hᴱ is a few hundred J/mol, the difference of H's
	// some 40 kJ/mol each.
	assert worst <= 1e-11
}

// checked holds a split to what it claims, from scratch: equal fugacities,
// the mass balance, and every phase stable, by a tangent-plane distance
// scanned over compositions from the model's phases alone (criterion 5).
// It returns the lowest TPD found.
fn checked(m thermo.PhaseModel, r flash.FlashResult, at string) !f64 {
	nc := r.z.len
	mut lowest := 0.0
	for ph in r.phases {
		x := ph.state.n
		sx := m.phase(r.t, r.p, x, .stable)!
		assert math.abs(sx.v - ph.state.v) <= 1e-12 * ph.state.v, '${at}: a phase not on its phase of lower G'
		// A scan of the simplex, at 1/40 in a binary or ternary.
		steps := 40
		for a in 1 .. steps {
			for b in 0 .. (if nc == 3 { steps - a } else { 1 }) {
				mut w := []f64{len: nc}
				if nc == 2 {
					w = [f64(a) / steps, 1.0 - f64(a) / steps]
				} else {
					w[0] = f64(a) / steps
					w[1] = f64(b) / steps
					w[2] = 1.0 - w[0] - w[1]
					if w[2] <= 0.0 {
						continue
					}
				}
				sw := m.phase(r.t, r.p, w, .stable)!
				mut tpd := 0.0
				for i in 0 .. nc {
					if w[i] > 0.0 {
						tpd += w[i] * (math.log(w[i]) + sw.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
					}
				}
				lowest = math.min(lowest, tpd)
			}
		}
	}
	if r.phases.len == 2 {
		y, x := r.phases[0].state.n, r.phases[1].state.n
		for i in 0 .. nc {
			assert math.abs(r.beta * y[i] + (1.0 - r.beta) * x[i] - r.z[i]) <= 1e-12, '${at}: the balance of ${i}'
			f := math.abs(math.log(y[i]) + r.phases[0].state.ln_phi[i] - math.log(x[i]) - r.phases[1].state.ln_phi[i])
			assert f <= 1e-10, '${at}: |ln fV − ln fL| of ${i} is ${f}'
		}
	}
	return lowest
}

// The flashes on a γ–φ fluid against thermo's FlashVL, on the 'Poynting'
// basis with the ideal gas (criterion 5): PT, bubble and dew T and P, each
// answer checked from scratch; and PH and PS taken back to each PT
// flash's T.
fn test_the_flashes_against_flashvl() {
	f := fixture()!
	db := components.load_default()!
	mut fluids := map[string]thermo.Fluid{}
	mut worst := Worst{}
	mut counted := map[string]int{}
	mut lowest := 0.0
	for c in f.flash {
		key := '${c.model} ${c.cas}'
		if key !in fluids {
			fluids[key] = db.gamma_phi(c.cas, liquid: liquid_of(c.model))!.fluid
		}
		fl := fluids[key] or { continue }
		m := fl.model
		at := '${key} ${c.kind} ${c.spec}'
		r := match c.kind {
			'pt' { flash.pt(m, c.spec[0], c.spec[1], c.z)! }
			'bubble_t' { flash.bubble_t(m, c.spec[0], c.z)! }
			'dew_t' { flash.dew_t(m, c.spec[0], c.z)! }
			'bubble_p' { flash.bubble_p(m, c.spec[0], c.z)! }
			else { flash.dew_p(m, c.spec[0], c.z)! }
		}
		lowest = math.min(lowest, checked(m, r, at)!)
		want := c.thermo
		counted['${c.kind} ${r.phases.len}']++
		if c.kind == 'pt' {
			if want.vf == 0.0 || want.vf == 1.0 {
				assert r.phases.len == 1 && r.beta == want.vf, '${at}: β ${r.beta}, thermo ${want.vf}'
			} else {
				assert r.phases.len == 2, '${at}: one phase, thermo VF ${want.vf}'
				worst.note('PT β', math.abs(r.beta - want.vf))
				for i in 0 .. c.z.len {
					worst.note('PT x, y', math.abs(r.phases[1].state.n[i] - want.x[i]))
					worst.note('PT x, y', math.abs(r.phases[0].state.n[i] - want.y[i]))
				}
			}
			// PH and PS from the answer's own H and S, back to its T.
			mut h, mut s := 0.0, 0.0
			for ph in r.phases {
				pr := fl.props(ph.state)
				h += ph.fraction * pr.h
				s += ph.fraction * pr.s
			}
			for kind, back in {
				'PH': flash.ph(fl, c.spec[1], h, c.z)!
				'PS': flash.ps(fl, c.spec[1], s, c.z)!
			} {
				worst.note('${kind} T', math.abs(back.t / c.spec[0] - 1.0))
				assert back.phases.len == r.phases.len, '${at}: ${kind} gave ${back.phases.len} phases'
				worst.note('${kind} β', math.abs(back.beta - r.beta))
			}
			continue
		}
		// a bubble or dew point: T or P, and the incipient phase
		if c.kind.ends_with('_t') {
			worst.note('boundary T', math.abs(r.t / want.t - 1.0))
		} else {
			worst.note('boundary P', math.abs(r.p / want.p - 1.0))
		}
		incipient := if c.kind.starts_with('bubble') {
			r.phases[0].state.n
		} else {
			r.phases[1].state.n
		}
		other := if c.kind.starts_with('bubble') { want.y } else { want.x }
		for i in 0 .. c.z.len {
			worst.note('boundary composition', math.abs(incipient[i] - other[i]))
		}
	}
	println('flashes: ${counted}; lowest TPD in the scans ${lowest:.1e}')
	mut keys := worst.m.keys()
	keys.sort()
	for k in keys {
		println('  ${k}: ${worst.m[k]:.2e}')
	}
	// Measured, under tcc: β 6.7e-11, where a split lies within 0.01 K of an
	// azeotrope and β moves as 1/(K − 1); every composition 1e-12; bubble
	// and dew points 8.6e-13 in P, 5.7e-14 in T; PH and PS to their T within
	// 1.1e-11, etol's 1e-10 of RT in H over Cp. Criterion 5 asks 1e-8.
	limits := {
		'PT β':                 1e-9
		'PT x, y':              1e-10
		'boundary T':           1e-12
		'boundary P':           1e-11
		'boundary composition': 1e-11
		'PH T':                 1e-10
		'PH β':                 1e-9
		'PS T':                 1e-10
		'PS β':                 1e-9
	}
	for k, v in worst.m {
		assert v <= limits[k] or { 0.0 }, '${k}: ${v}'
	}
	assert worst.m.len == limits.len, '${worst.m.keys()}'
	// No composition in any scan lies below a phase's tangent plane.
	assert lowest >= -1e-10, 'a TPD of ${lowest}'
	assert (counted['pt 2'] or { 0 }) >= 20 && (counted['pt 1'] or { 0 }) >= 20
}

// Each azeotrope the model has at 1 atm, found where a scan of bubble
// points changes the sign of y₁ − x₁, at thermo's x and T (criterion 5).
fn test_each_azeotrope_is_found() {
	f := fixture()!
	db := components.load_default()!
	mut worst_x, mut worst_t := 0.0, 0.0
	mut found := 0
	for c in f.azeotropes {
		m := db.gamma_phi(c.cas, liquid: liquid_of(c.model))!.fluid.model
		mut ours := []flash.FlashResult{}
		mut last := 0.0
		for k in 1 .. 50 {
			x1 := f64(k) / 50.0
			b := flash.bubble_t(m, c.p, [x1, 1.0 - x1])!
			d := b.phases[0].state.n[0] - x1
			if k > 1 && (d > 0.0) != (last > 0.0) {
				ours << flash.azeotrope_t(m, c.p, [x1, 1.0 - x1])!
			}
			last = d
		}
		assert ours.len == c.azeotropes.len, '${c.model} ${c.cas}: ${ours.len} azeotropes, thermo ${c.azeotropes.len}'
		for k, a in ours {
			assert a.iterations <= 3, 'Newton took ${a.iterations} iterations, recorded 3'
			worst_x = math.max(worst_x, math.abs(a.z[0] - c.azeotropes[k].x1))
			worst_t = math.max(worst_t, math.abs(a.t / c.azeotropes[k].t - 1.0))
			found++
		}
	}
	println('azeotropes: ${found}, x to ${worst_x:.1e}, T to ${worst_t:.1e}')
	// Measured, 2.0e-14 and 7.8e-16: both are Newton's to tol, and thermo's
	// Brent's to 1e-14 in x.
	assert found == 2 && worst_x <= 1e-12 && worst_t <= 1e-13
}

// azeotrope_t from the azeotrope itself, where its bubble point collapses
// onto the feed, and azeotrope_p at its T, give back its P and x. And a
// heterogeneous azeotrope, water and 1-butanol's, whose one liquid would
// split, is refused as one that needs three phases (step 7's): from the
// butanol-rich side, every start reaches y = x at 365.86 K, where the
// liquid's tm is −0.031; inside the gap, no one-liquid bubble point starts
// it.
fn test_azeotropes_at_a_given_t_and_heterogeneous_ones() {
	f := fixture()!
	db := components.load_default()!
	for c in f.azeotropes {
		m := db.gamma_phi(c.cas, liquid: liquid_of(c.model))!.fluid.model
		for a in c.azeotropes {
			x := [a.x1, 1.0 - a.x1]
			at_p := flash.azeotrope_t(m, c.p, x)!
			back := flash.azeotrope_p(m, at_p.t, [0.5, 0.5])!
			// Newton's iterations, recorded, which can't grow: 1 and 2 from the
			// azeotrope itself, 5 and 4 from x = 0.5 (a wrong Jacobian column
			// still converges, more slowly).
			assert at_p.iterations <= 2 && back.iterations <= 5, '${at_p.iterations}, ${back.iterations}'
			assert math.abs(back.p / c.p - 1.0) <= 1e-12, '${c.cas}: P ${back.p}'
			assert math.abs(back.z[0] - at_p.z[0]) <= 1e-12, '${c.cas}: x ${back.z}'
			assert back.phases[0].state.n == back.phases[1].state.n
		}
	}
	m := db.gamma_phi(['water', '1-butanol'])!.fluid.model
	// One component is its own boiling point, no azeotrope.
	if r := flash.azeotrope_t(m, 101325.0, [1.0, 0.0]) {
		assert false, 'an azeotrope of water alone at ${r.t} K'
	} else {
		assert err.msg().contains('an azeotrope needs two components or more, and x0 has 1'), err.msg()
	}
	if r := flash.azeotrope_t(m, 101325.0, [0.3, 0.7]) {
		assert false, 'a homogeneous azeotrope of water and 1-butanol at ${r.t} K, x ${r.z}'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .more_phases, err.msg()
		assert err.msg().contains('a heterogeneous azeotrope'), err.msg()
	}
}

// A correlation in other units, by its scale: each vapour pressure in kPa
// (A less ln 1000, scale 1000) and each density in mol/m³ (A times 1000,
// scale 1) gives the same liquid. ChemSep's are all in Pa and kmol/m³, so
// only this sees whether the vapour pressure's scale, 1 there, is applied.
fn test_a_correlation_in_other_units() {
	db := components.load_default()!
	s := db.gamma_phi(['methanol', 'water'])!
	g := s.model
	mut psat, mut density := g.psat.clone(), g.density.clone()
	for i in g.solvents {
		c := g.psat[i].c
		assert c.eq == 101 && g.psat[i].scale == 1.0
		psat[i] = thermo.Scaled{tdep.new(101, [c.a - math.log(1000.0), c.b, c.c, c.d, c.e], c.tmin,
			c.tmax, c.tc)!, 1000.0}
		r := g.density[i].c
		assert g.density[i].scale == 1000.0
		density[i] = thermo.Scaled{tdep.new(r.eq, [1000.0 * r.a, r.b, r.c, r.d, r.e], r.tmin, r.tmax,
			r.tc)!, 1.0}
	}
	other := thermo.new_phase_fluid(thermo.new_gamma_phi(thermo.GammaPhi{
		...g
		psat:    psat
		density: density
	})!, s.fluid.ideal)!
	for t in [300.0, 340.0] {
		a := s.fluid.model.phase(t, 2e5, [0.3, 0.7], .liquid)!
		b := other.model.phase(t, 2e5, [0.3, 0.7], .liquid)!
		assert math.abs(b.v / a.v - 1.0) <= 1e-14, '${t} K: V'
		for i in 0 .. 2 {
			assert math.abs(b.ln_phi[i] - a.ln_phi[i]) <= 1e-13, '${t} K: ln φ'
		}
		ra, rb := s.fluid.model.phase_residual(a), other.model.phase_residual(b)
		assert gap(rb.h, ra.h, 1.0) <= 1e-12 && gap(rb.cp, ra.cp, 1.0) <= 1e-12, '${t} K: H, Cp'
	}
}

// n·∂ln φᵢ/∂nⱼ against Richardson's extrapolation of differences of ln φ
// itself, in every nⱼ at one mole (a check of the claim that thermo can't
// make: it has no composition derivatives):
// - methanol and water by NRTL with CO₂ and hydrogen dissolved, so a
//   solvent's dn is the activity model's over the solvents' own fractions,
//   divided by their share, and a solute's is its mixed H's;
// - the mixed solvent's fluid, its glycol absent, whose column, one-sided,
//   is a solute's H moving towards the glycol's.
fn test_dn_against_differences() {
	f := fixture()!
	db := components.load_default()!
	m := db.gamma_phi(['methanol', 'water', 'carbon dioxide', 'hydrogen'],
		henry: ['carbon dioxide']
	)!.fluid.model
	mut worst := 0.0
	for t in [300.0, 330.0] {
		worst = math.max(worst, dn_gap(m, t, 3e6, [0.3, 0.6, 0.06, 0.04])!)
	}
	hm := f.henry_mixed
	s := db.gamma_phi(hm.cas, henry: hm.cas[2..], ideal: true, poynting: false)!
	mut solutes := []thermo.Solute{}
	for sol in s.model.solutes {
		mut henrys := sol.henrys.clone()
		for c in hm.constants.filter(it.solute == sol.i && it.solvent == 1) {
			k := c.abcde
			henrys << thermo.Henry{
				solvent: c.solvent
				a:       k[0]
				b:       k[1]
				c:       k[2]
				d:       k[3]
				e:       k[4]
				tmin:    c.tmin
				tmax:    c.tmax
			}
		}
		solutes << thermo.Solute{sol.i, henrys}
	}
	mixed := thermo.new_gamma_phi(thermo.GammaPhi{
		...s.model
		solutes: solutes
	})!
	one_sided := dn_gap(mixed, 320.0, 3e6, [0.9, 0.0, 0.06, 0.04])!
	println('dn against differences: ${worst:.1e}, an absent solvent one-sided ${one_sided:.1e}')
	// Measured: central differences to Richardson's truncation, O(h⁴), and
	// ln φ's rounding over h; one-sided to O(h²).
	assert worst <= 1e-9, '${worst}'
	assert one_sided <= 1e-6, '${one_sided}'
}

// dn_gap is the worst gap between n·∂ln φᵢ/∂nⱼ and Richardson's
// extrapolation of differences of ln φ in nⱼ, step h = 1e-4 mol and h/2:
// central, or forward where nⱼ is 0.
fn dn_gap(m thermo.PhaseModel, t f64, p f64, x []f64) !f64 {
	nc := x.len
	dn := m.ln_phi_derivs(m.phase(t, p, x, .liquid)!)!.dn
	base := m.phase(t, p, x, .liquid)!.ln_phi
	mut worst := 0.0
	for j in 0 .. nc {
		h := 1e-4
		mut d := [][]f64{}
		for step in [h, h / 2.0] {
			mut up := x.clone()
			up[j] += step
			fu := m.phase(t, p, up, .liquid)!.ln_phi
			if x[j] > step {
				mut down := x.clone()
				down[j] -= step
				fd := m.phase(t, p, down, .liquid)!.ln_phi
				d << []f64{len: nc, init: (fu[index] - fd[index]) / (2.0 * step)}
			} else {
				d << []f64{len: nc, init: (fu[index] - base[index]) / step}
			}
		}
		for i in 0 .. nc {
			r := if x[j] > h { (4.0 * d[1][i] - d[0][i]) / 3.0 } else { 2.0 * d[1][i] - d[0][i] }
			worst = math.max(worst, math.abs(dn[i * nc + j] - r) / math.max(math.abs(r), 1e-3))
		}
	}
	return worst
}

// A phase of 2.5 mol is 2.5 times a phase of one: V, H, S and Cp scale,
// Z, ln φ and its derivatives don't, liquid and vapour, with a solute. And
// each phase's (∂V/∂T)_P, for the Joule–Thomson coefficient, is its own
// volume's slope.
fn test_a_phase_is_extensive() {
	db := components.load_default()!
	m := db.gamma_phi(['methanol', 'water', 'carbon dioxide'],
		henry:  ['carbon dioxide']
		vapour: .pr
	)!.fluid.model
	x := [0.3, 0.69, 0.01]
	n := x.map(2.5 * it)
	for root in [thermo.Root.liquid, .vapor] {
		a := m.phase(330.0, 2e5, x, root)!
		b := m.phase(330.0, 2e5, n, root)!
		assert math.abs(b.v / a.v - 2.5) <= 1e-14, '${root}: V'
		assert math.abs(b.z / a.z - 1.0) <= 1e-14, '${root}: Z ${a.z} ${b.z}'
		for i in 0 .. 3 {
			assert math.abs(b.ln_phi[i] - a.ln_phi[i]) <= 1e-13, '${root}: ln φ'
		}
		da, db_ := m.ln_phi_derivs(a)!, m.ln_phi_derivs(b)!
		for i in 0 .. 3 {
			assert math.abs(db_.dt[i] - da.dt[i]) <= 1e-13 * math.max(math.abs(da.dt[i]), 1e-3)
		}
		ra, rb := m.phase_residual(a), m.phase_residual(b)
		assert gap(rb.h, 2.5 * ra.h, 1.0) <= 1e-13 && gap(rb.cp, 2.5 * ra.cp, 1.0) <= 1e-13, '${root}: H, Cp'
		assert gap(rb.dv_dt, 2.5 * ra.dv_dt, 1e-12) <= 1e-13, '${root}: dV/dT'
		// (∂V/∂T)_P against Richardson's central differences of V itself.
		h := 0.01
		mut v := []f64{}
		for dt in [h, -h, h / 2.0, -h / 2.0] {
			v << m.phase(330.0 + dt, 2e5, x, root)!.v
		}
		fd := (4.0 * (v[2] - v[3]) / h - (v[0] - v[1]) / (2.0 * h)) / 3.0
		assert gap(ra.dv_dt, fd, 1e-12) <= 1e-9, '${root}: dV/dT ${ra.dv_dt} against ${fd}'
	}
}

// With the ideal gas for its vapour, the vapour's Joule–Thomson coefficient
// is 0 and its speed of sound √(γRT/M), γ = Cp/(Cp − R), from what the seam
// gives (its (∂V/∂T)_P and (∂P/∂V)_T).
fn test_the_ideal_vapour_is_the_ideal_gas() {
	db := components.load_default()!
	fl := db.gamma_phi(['methanol', 'water'])!.fluid
	x := [0.4, 0.6]
	st := fl.model.phase(380.0, 1e5, x, .vapor)!
	pr := fl.props(st)
	mw := 0.4 * fl.model.components()[0].mw + 0.6 * fl.model.components()[1].mw
	r := thermo.r_gas
	w := math.sqrt(pr.cp / (pr.cp - r) * r * 380.0 / mw)
	assert math.abs(pr.w / w - 1.0) <= 1e-14, 'w ${pr.w} against ${w}'
	// μ_JT = (T(∂V/∂T)_P − V)/Cp, whose bracket cancels to V's ulps.
	assert math.abs(pr.mu_jt) <= 1e-15 * st.v / pr.cp, 'μ_JT ${pr.mu_jt}'
	assert math.abs(pr.cp - pr.cv - r) <= 1e-12 * r, 'Cp − Cv'
}

// new_gamma_phi refuses what the builder never makes: a solvent's
// correlation without exact derivatives, and a solute listed twice.
fn test_new_gamma_phi_refuses_bad_parts() {
	db := components.load_default()!
	g := db.gamma_phi(['water', 'nitrogen'])!.model
	mut psat := g.psat.clone()
	psat[0] = thermo.Scaled{tdep.new(100, [1.0, 0.0, 0.0, 0.0, 0.0], 273.15, 373.15, 0.0)!, 1.0}
	if _ := thermo.new_gamma_phi(thermo.GammaPhi{
		...g
		psat: psat
	})
	{
		assert false, 'a vapour pressure of form 100'
	} else {
		assert err.msg().contains('form 100, which has no exact derivatives'), err.msg()
	}
	if _ := thermo.new_gamma_phi(thermo.GammaPhi{
		...g
		solutes: [g.solutes[0], g.solutes[0]]
	})
	{
		assert false, 'a solute listed twice'
	} else {
		assert err.msg().contains('or listed twice'), err.msg()
	}
}

// The builder's guards on a Henry's-law bank that ChemSep's can't reach,
// every row of which is a gas in water by equation 121, on another data
// root's: a gas's constant in itself, which mustn't make it a solute (with
// it, water would be one, and nothing a solvent); a gas's constant in
// another solute, skipped; and another equation, refused by its line
// (ammonia's: a gas IAPWS has takes IAPWS's constant in water).
fn test_a_henry_bank_of_another_data_root() {
	db := components.load_default()!
	dir := os.join_path(os.temp_dir(), 'reflux_henry_${os.getpid()}')
	os.mkdir_all(dir)!
	defer {
		os.rmdir_all(dir) or {}
	}
	os.cp(os.join_path(db.vendor, 'nrtl.ipd'), os.join_path(dir, 'nrtl.ipd'))!
	base := os.read_file(os.join_path(db.vendor, 'henrys_law.ipd'))!.trim_right('\r\n')
	other := components.Database{
		...db
		vendor: dir
	}
	bank := os.join_path(dir, 'henrys_law.ipd')
	os.write_file(bank, base + '\n 7732-18-5    7732-18-5   1 121      273.15   373.15   20.0  0  0  0  1  Water\n')!
	s := other.gamma_phi(['water', 'nitrogen'])!
	assert s.solutes == ['7727-37-9'], '${s.solutes}'
	os.write_file(bank, base + '\n 7727-37-9    1333-74-0   1 121      273.15   373.15   20.0  0  0  0  1  Nitrogen\n')!
	s2 := other.gamma_phi(['water', 'nitrogen', 'hydrogen'])!
	assert s2.solutes == ['7727-37-9', '1333-74-0'], '${s2.solutes}'
	assert s2.model.solutes[0].henrys.map(it.solvent) == [0]
	os.write_file(bank, base + '\n 7664-41-7    7732-18-5   1 120      273.15   373.15   20.0  0  0  0  1  Ammonia\n')!
	if _ := other.gamma_phi(['water', 'ammonia'], henry: ['ammonia']) {
		assert false, 'equation 120'
	} else {
		assert err.msg().contains('is equation 120, not 121'), err.msg()
	}
}

// The vapour by PR is the cubic fluid's own vapour root, bit for bit, with
// its bank's kᵢⱼ or without, as asked; and methanol–water's kᵢⱼ moves it.
fn test_the_vapour_is_the_cubics() {
	db := components.load_default()!
	x := [0.4, 0.6]
	mut seen := []f64{}
	for kij in [true, false] {
		g := db.gamma_phi(['methanol', 'water'], vapour: .pr, kij: kij)!.fluid.model
		eos := db.fluid(.pr, ['methanol', 'water'], kij: kij)!.helmholtz()!
		a := g.phase(400.0, 2e5, x, .vapor)!
		b := thermo.solve_tp(eos, 400.0, 2e5, x, .vapor)!
		assert a.v == b.v && a.ln_phi == b.ln_phi, 'kij ${kij}'
		seen << a.ln_phi[0]
	}
	assert seen[0] != seen[1], "kᵢⱼ didn't move the vapour"
}

// IAPWS's Henry's constants for gases in water (G7-04) against the iapws
// package's (testdata/m6_henry_iapws.json, tools/gen_m6_henry.py): the
// table's A, B, C and ranges as the package has them, ln kH at 25
// temperatures across each gas's range, and its T-derivatives against
// Richardson's differences of itself. And a gas in water takes it, saying
// so: hydrogen's ln φ in a liquid of water is ln kH − ln P.
fn test_iapws_henrys_constants() {
	f := json2.decode[IapwsFixture](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_henry_iapws.json'))!)!
	assert f.gases.len == 14
	mut worst, mut worst_d := 0.0, 0.0
	for g in f.gases {
		h := components.iapws_henry_in_water(g.cas, 0) or {
			assert false, '${g.name}: not in the table'
			continue
		}
		assert [h.a, h.b, h.c] == g.abc && h.tmin == g.tmin && h.tmax == g.tmax, g.name
		for k, t in g.t {
			worst = math.max(worst, math.abs(h.jet(numeric.along(t)).v / g.ln_kh[k] - 1.0))
		}
		for t in [g.tmin + 5.0, 0.5 * (g.tmin + g.tmax), g.tmax - 5.0] {
			j := h.jet(numeric.along(t))
			step := 0.05
			mut v := []f64{}
			mut d := []f64{}
			for dt in [step, -step, step / 2.0, -step / 2.0] {
				q := h.jet(numeric.along(t + dt))
				v << q.v
				d << q.d
			}
			dv := (4.0 * (v[2] - v[3]) / step - (v[0] - v[1]) / (2.0 * step)) / 3.0
			dd := (4.0 * (d[2] - d[3]) / step - (d[0] - d[1]) / (2.0 * step)) / 3.0
			worst_d = math.max(worst_d, math.max(math.abs(j.d - dv) / math.abs(j.d), math.abs(j.dd - dd) / math.max(math.abs(j.dd),
				1e-6)))
		}
	}
	println("IAPWS's ln kH: ${worst:.1e}; its derivatives against differences ${worst_d:.1e}")
	// Measured: the package and reflux take τ^b and e^τ each their own way,
	// a few ulps of ln kH; the differences' truncation and rounding.
	assert worst <= 1e-14, '${worst}'
	assert worst_d <= 1e-7, '${worst_d}'
	db := components.load_default()!
	s := db.gamma_phi(['water', 'hydrogen'])!
	assert s.notes.any(it == "Hydrogen in Water: IAPWS's Henry's constant (G7-04), in place of henrys_law.ipd's"), '${s.notes}'
	// Sulfur hexafluoride has no row in the bank, so IAPWS's replaces none.
	sf6 := db.gamma_phi(['water', 'sulfur hexafluoride'], henry: ['sulfur hexafluoride'])!
	assert sf6.notes.any(it == "Sulfur hexafluoride in Water: IAPWS's Henry's constant (G7-04)"), '${sf6.notes}'
	assert sf6.model.solutes[0].henrys.map(it.form) == [thermo.HenryForm.iapws]
	h2 := f.gases.filter(it.name == 'H2')[0]
	for k in [0, 12, 24] {
		st := s.fluid.model.phase(h2.t[k], 3e7, [0.999, 0.001], .liquid)!
		assert math.abs(st.ln_phi[1] + math.log(3e7) - h2.ln_kh[k]) <= 1e-13 * h2.ln_kh[k], '${h2.t[k]} K'
	}
}

// A workspace for one thread (thermo.Scoped) keeps what depends on T alone
// between calls, and gives the bits of a model that makes its own each
// time, at one T, at another, and back.
fn test_a_workspace_gives_the_same_bits() {
	db := components.load_default()!
	plain := db.gamma_phi(['67-64-1', '67-56-1', '7732-18-5'], vapour: .pr)!.fluid.model
	kept := thermo.for_one_thread(plain)
	x := [0.2, 0.3, 0.5]
	for t in [320.0, 345.0, 320.0, 320.0, 360.0] {
		for root in [thermo.Root.liquid, .vapor, .stable] {
			a := plain.phase(t, 2e5, x, root)!
			b := kept.phase(t, 2e5, x, root)!
			assert a.v == b.v && a.ln_phi == b.ln_phi, '${t} K ${root}'
			da, db_ := plain.ln_phi_derivs(a)!, kept.ln_phi_derivs(b)!
			assert da.dt == db_.dt && da.dp == db_.dp && da.dn == db_.dn, '${t} K ${root}'
			ra, rb := plain.phase_residual(a), kept.phase_residual(b)
			assert ra.h == rb.h && ra.s == rb.s && ra.cp == rb.cp, '${t} K ${root}'
		}
	}
	// A refused T leaves nothing kept: at 250 K acetone's and methanol's
	// parts are filled before water's range refuses, and the next call at
	// the T before must not take them for its own.
	if _ := kept.phase(250.0, 2e5, x, .liquid) {
		assert false, 'a liquid below the range of water'
	}
	a := plain.phase(360.0, 2e5, x, .liquid)!
	b := kept.phase(360.0, 2e5, x, .liquid)!
	assert a.ln_phi == b.ln_phi, 'after a refusal'
}

// What the γ–φ fluid refuses, it refuses by name, and the flash passes the
// reason on: a T outside a correlation's range (never the vapour in the
// liquid's place, D4), a liquid of solutes alone, a pair no bank lists, a
// solute with no constant in the solvents present; and a model's reason
// reaches the caller through numeric.newton.
fn test_the_refusals() {
	db := components.load_default()!
	m := db.gamma_phi(['methanol', 'water'])!.fluid.model
	for root in [thermo.Root.liquid, .stable] {
		if _ := m.phase(250.0, 1e5, [0.5, 0.5], root) {
			assert false, 'a ${root} phase below the range of the vapour pressure of water'
		} else {
			assert err.msg().contains('fitted from 263.15 K to 647.29 K, not at 250.0 K'), err.msg()
		}
	}
	assert m.phase(250.0, 1e5, [0.5, 0.5], .vapor)!.z == 1.0
	if _ := flash.pt(m, 250.0, 1e5, [0.5, 0.5]) {
		assert false, 'a flash below the range'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .no_volume_root, err.msg()
		assert err.msg().contains('not at 250.0 K'), err.msg()
	}
	// The bubble point's Newton can't start below the range: its reason,
	// not "the system cannot be evaluated" alone.
	if _ := flash.bubble_p(m, 250.0, [0.5, 0.5]) {
		assert false, 'a bubble point below the range'
	} else {
		assert err.msg().contains('not at 250.0 K'), err.msg()
	}
	// A solute outside its constant's range: nitrogen's in water, IAPWS's,
	// ends at 636.46 K.
	n2 := db.gamma_phi(['water', 'nitrogen'])!.fluid.model
	if _ := n2.phase(640.0, 5e6, [0.999, 0.001], .liquid) {
		assert false, "nitrogen beyond its Henry's-law range"
	} else {
		assert err.msg().contains("Nitrogen's Henry's-law constant in Water is fitted from 278.12 K to 636.46 K, not at 640.0 K"), err.msg()
	}
	h2 := db.gamma_phi(['water', 'hydrogen'])!.fluid.model
	if _ := h2.phase(330.0, 1e5, [0.0, 1.0], .liquid) {
		assert false, 'a liquid of hydrogen alone'
	} else {
		assert err.msg().contains("Henry's-law solutes alone"), err.msg()
	}
	if _ := db.gamma_phi(['methanol', 'water', 'carbon dioxide']) {
		assert false, 'NRTL without the pairs'
	} else {
		assert err.msg().contains("doesn't list Methanol with Carbon dioxide"), err.msg()
	}
	// Hydrogen has no constant in methanol, so with it no T has a liquid.
	if _ := db.gamma_phi(['methanol', 'hydrogen']) {
		assert false, 'hydrogen as a solvent of methanol'
	} else {
		assert err.msg().starts_with('no temperature has a liquid of Hydrogen and Methanol: Hydrogen'), err.msg()
		assert err.msg().contains("henrys_law.ipd has no constant for Hydrogen in the others, so it can't be a Henry's-law solute"), err.msg()
	}
	if _ := db.gamma_phi(['water', 'nitrogen'], henry: ['hydrogen']) {
		assert false, 'a solute named that is not a compound'
	} else {
		assert err.msg().contains("Hydrogen, named as a Henry's-law solute, isn't one of the compounds"), err.msg()
	}
	if _ := db.gamma_phi(['methanol', 'hydrogen'], henry: ['hydrogen']) {
		assert false, 'a solute with no constant in methanol'
	} else {
		assert err.msg().contains('has no constant for it'), err.msg()
	}
}

// PH's bracket steps back from a T the model refuses, rather than fail:
// methanol and water's liquid at 266 K, 3 K inside water's range, from
// the ideal gas's start well below it.
fn test_ph_steps_back_inside_the_range() {
	db := components.load_default()!
	fl := db.gamma_phi(['methanol', 'water'])!.fluid
	z := [0.4, 0.6]
	r := flash.pt(fl.model, 266.0, 101325.0, z)!
	h := fl.props(r.phases[0].state).h
	back := flash.ph(fl, 101325.0, h, z)!
	assert math.abs(back.t / 266.0 - 1.0) <= 1e-10, '${back.t}'
}
