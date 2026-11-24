// Against DWSIM 9.0.5's activity packages (M6 step 12, D17, criterion 11),
// by state, every difference explained (testdata/m6_dwsim.json, from
// tools/gen_m6_dwsim.py, which reads DWSIM's numbers and works its gamma
// out again from its own code):
// - NRTL, UNIFAC and Dortmund's are reflux's, to its Poynting term's R;
// - DWSIM's Poynting term has |P − Psat|, so above Psat it raises φ;
// - its UNIQUAC for methanol–water takes two of its four τ's from a second,
//   inconsistent record of its own;
// - its Wilson takes R = 1.9872 and volumes from the liquid's density at
//   298.15 K, not ChemSep's Wilson volumes;
// - its UNIFAC-LL is its VLE UNIFAC with each compound's group fraction
//   normalised by the Q of its Dortmund group numbers;
// - its flashes are its model's answers, and its three-phase flash misses a
//   split the model has.
import activity
import components
import flash
import json2
import math
import os
import thermo

const db = components.load_default() or { panic(err) }
const cas = ['67-56-1', '7732-18-5']
const r_si = 8.314462618

fn fixture() !map[string]json2.Any {
	return json2.decode[json2.Any](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_dwsim.json'))!)!.as_map()
}

fn floats(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

// At is a state, T and x, with a value there, a γ.
struct At {
	x []f64
	t f64
	v []f64
}

// gammas is one package's γ by state, the named value of each.
fn gammas(f map[string]json2.Any, package string, key string) []At {
	mut out := []At{}
	for g in f['packages']!.as_map()[package]!.as_map()['gamma']!.as_array() {
		m := g.as_map()
		out << At{floats(m['x']!), m['T']!.f64(), floats(m[key]!)}
	}
	return out
}

// worst is the largest |ln γ − ln want| over a package's states, of a
// model of reflux's.
fn worst(m activity.Model, states []At) !f64 {
	mut w := m.work()
	mut out := 0.0
	for s in states {
		activity.eval(m, s.t, s.x, activity.ln_gamma_only, mut w)!
		for i in 0 .. s.x.len {
			out = math.max(out, math.abs(w.ln_gamma[i] - math.log(s.v[i])))
		}
	}
	return out
}

// worst_of is the largest |ln a − ln b| between two of a package's values.
fn worst_of(a []At, b []At) f64 {
	mut out := 0.0
	for k, s in a {
		for i in 0 .. s.v.len {
			out = math.max(out, math.abs(math.log(s.v[i]) - math.log(b[k].v[i])))
		}
	}
	return out
}

// NRTL, UNIFAC and Dortmund's UNIFAC are reflux's, state by state. What's
// left, 1.7e-9 in ln γ, is DWSIM's Poynting term on its own R, against the
// SI's taken out here.
fn test_nrtl_and_unifacs_are_reflux_s() {
	f := fixture()!
	for package, m in {
		'NRTL':                       activity.Model(db.nrtl(cas)!.model)
		'UNIFAC':                     db.unifac(.original, cas)!
		'Modified UNIFAC (Dortmund)': db.unifac(.dortmund, cas)!
	} {
		w := worst(m, gammas(f, package, 'gamma'))!
		println('${package}: ${w:.1e}')
		assert w <= 5e-9, '${package}: ${w}'
	}
}

// DWSIM's Poynting term is exp(V_L·|P − Psat|/RT): where Psat is above P,
// methanol's at 360 K and 1 atm, the textbook's signed form leaves exactly
// twice the term between them, and elsewhere nothing.
fn test_dwsim_s_poynting_term_loses_its_sign() {
	f := fixture()!
	p := f['P']!.f64()
	comps := f['methanol_water']!.as_map()['compounds']!.as_array()
	pk := f['packages']!.as_map()['NRTL']!.as_map()
	mut psat := map[string][]f64{}
	mut rho := map[string][]f64{}
	for r in pk['psat']!.as_array() {
		psat['${r.as_map()['T']!.f64()}'] = floats(r.as_map()['psat']!)
	}
	for r in pk['rho']!.as_array() {
		rho['${r.as_map()['T']!.f64()}'] = floats(r.as_map()['rho']!)
	}
	signed := gammas(f, 'NRTL', 'gamma_signed')
	plain := gammas(f, 'NRTL', 'gamma')
	mut above := 0
	for k, s in plain {
		key := '${s.t}'
		for i in 0 .. s.x.len {
			ps := psat[key][i]
			v := comps[i].as_map()['mw']!.f64() / rho[key][i]
			gap := math.log(signed[k].v[i]) - math.log(s.v[i])
			if ps > p {
				assert math.abs(gap - 2.0 * v * (ps - p) / (r_si * s.t)) <= 1e-12, '${s.t} ${i}: ${gap}'
				above++
			} else {
				assert gap == 0.0, '${s.t} ${i}: ${gap}'
			}
		}
	}
	assert above == 5 // methanol at 360 K, at each x
}

// DWSIM's UNIQUAC is its code on the two records it holds for
// methanol–water: ChemSep's, −337.1298 and 549.2958 cal/mol, and a second
// of its own, −6.359 and −1.192, from which two of the four τ's come. On
// ChemSep's alone it's DECHEMA's form, reflux's; the second moves ln γ by
// up to 0.36.
fn test_dwsim_s_uniquac_takes_a_second_record() {
	f := fixture()!
	rec := f['packages']!.as_map()['UNIQUAC']!.as_map()['records']!.as_map()
	assert rec['Methanol|Water']!.as_map()['A12']!.f64() == -337.1298
	assert rec['Water|Methanol']!.as_map()['A12']!.f64() == -6.359072
	dwsim := gammas(f, 'UNIQUAC', 'gamma')
	formula := gammas(f, 'UNIQUAC', 'formula')
	one := gammas(f, 'UNIQUAC', 'formula_one_record')
	assert worst_of(dwsim, formula) <= 5e-9
	m := db.uniquac(cas)!.model
	assert worst(m, one)! <= 1e-12
	gap := worst(m, dwsim)!
	println('UNIQUAC: reflux against DWSIM ${gap:.3f}, explained to ${worst_of(dwsim, formula):.1e}')
	assert gap > 0.3
}

// DWSIM's Wilson is Λᵢⱼ = (Vⱼ/Vᵢ)·exp(−Aᵢⱼ/(1.9872·T)), A ChemSep's, V the
// liquid's molar volume from its density at 298.15 K. Reflux's Wilson on
// those, written as expressions, is DWSIM's; on ChemSep's Wilson volumes and
// DECHEMA's 1.98721, its own, it's 1.5e-3 away.
fn test_dwsim_s_wilson_volumes_and_r() {
	f := fixture()!
	pk := f['packages']!.as_map()['Wilson']!.as_map()
	bips := pk['bips']!.as_map()
	v := floats(pk['volumes_298']!)
	a12, a21 := bips['Methanol|Water']!.f64(), bips['Water|Methanol']!.f64()
	assert a12 == 205.3 && a21 == 482.16 // wilson.ipd's, cal/mol
	dw := activity.new_wilson(['methanol', 'water'], [
		activity.WilsonPair{
			i:         0
			j:         1
			lambda_ij: activity.plain('${v[1] / v[0]} * exp(-${a12} K / (1.9872*T))', 'Λ')!
			lambda_ji: activity.plain('${v[0] / v[1]} * exp(-${a21} K / (1.9872*T))', 'Λ')!
		},
	])!
	dwsim := gammas(f, 'Wilson', 'gamma')
	assert worst_of(dwsim, gammas(f, 'Wilson', 'formula')) <= 5e-9
	assert worst(dw, dwsim)! <= 5e-9
	own := worst(db.wilson(cas)!.model, dwsim)!
	println('Wilson: reflux against DWSIM ${own:.2e}')
	assert own > 1e-3 && own < 2e-3
}

// DWSIM's UNIFAC-LL is its VLE UNIFAC, methanol's CH3OH group and water's
// H2O with the VLE table's a, each group fraction ν·Q/q divided not by its
// own Q but by the Q its Dortmund group number has in the VLE subgroup
// table: OH's 1.2 for methanol (Dortmund 15), CH3OH's 1.432 for water (16).
// With fractions of 1 it's reflux's original UNIFAC. Magnussen's LLE table,
// reflux's UNIFAC-LLE, takes methanol as CH3 and OH, 0.24 away.
fn test_dwsim_s_unifac_ll_is_its_vle_unifac() {
	f := fixture()!
	g := f['packages']!.as_map()['UNIFAC-LL']!.as_map()['groups']!.as_map()['compounds']!.as_array()
	assert g.map(it.as_map()['e']!.f64()) == [1.432 / 1.2, 1.4 / 1.432]
	dwsim := gammas(f, 'UNIFAC-LL', 'gamma')
	assert worst_of(dwsim, gammas(f, 'UNIFAC-LL', 'formula')) <= 5e-9
	assert worst(db.unifac(.original, cas)!, gammas(f, 'UNIFAC-LL', 'formula_e1'))! <= 1e-12
	lle := worst(db.unifac(.lle, cas)!, dwsim)!
	println('UNIFAC-LL: reflux UNIFAC-LLE against DWSIM ${lle:.3f}')
	assert lle > 0.2
}

// Each of DWSIM's NRTL flashes is its model's answer: at its own x and y,
// yᵢP = xᵢγᵢPᵢˢᵃᵗ·exp(Vᵢ|P − Pᵢˢᵃᵗ|/RT), γ reflux's NRTL (which is DWSIM's,
// above), Psat and V DWSIM's. Its UNIQUAC calls each feed one liquid, as its
// two-record γ makes it: the bubble pressure it gives is below 1 atm.
// Reflux's own flashes differ by DWSIM's two departures, the Poynting sign
// and ChemSep's vapour pressures to more figures than ChemSep LITE's.
fn test_dwsim_s_flashes_are_its_model_s() {
	f := fixture()!
	p := f['P']!.f64()
	m := db.nrtl(cas)!.model
	mut w := m.work()
	fl := db.gamma_phi(cas)!.fluid
	mut worst_res, mut worst_beta := 0.0, 0.0
	for e in f['flashes']!.as_array() {
		r := e.as_map()
		t := r['T']!.f64()
		z := floats(r['z']!)
		psat, rho := floats(r['psat']!), floats(r['rho']!)
		mw := [0.03204186, 0.01801528]
		poy := [0, 1].map(math.exp(mw[it] / rho[it] * math.abs(p - psat[it]) / (r_si * t)))
		if r['package']!.str() == 'UNIQUAC' {
			gz := floats(r['gamma_z']!)
			bubble := z[0] * gz[0] * psat[0] * poy[0] + z[1] * gz[1] * psat[1] * poy[1]
			assert r['beta']!.f64() == 0.0 && bubble < p, '${t}: ${bubble}'
			continue
		}
		x, y := floats(r['x1']!), floats(r['y']!)
		activity.eval(m, t, x, activity.ln_gamma_only, mut w)!
		for i in 0 .. 2 {
			res := math.log(y[i] * p) - math.log(x[i] * psat[i] * poy[i]) - w.ln_gamma[i]
			worst_res = math.max(worst_res, math.abs(res))
		}
		own := flash.pt(fl.model, t, p, z)!
		worst_beta = math.max(worst_beta, math.abs(own.beta - r['beta']!.f64()))
	}
	println('DWSIM NRTL flashes: its own equations to ${worst_res:.1e}; reflux β differs by ${worst_beta:.1e}')
	// Measured: 8.7e-7, its nested loops' K-values, for all the 1e-10 asked
	// of them.
	assert worst_res <= 2e-6, '${worst_res}'
	assert worst_beta > 1e-5 && worst_beta < 1e-2
}

// DWSIM's three-phase flash of water–1-butanol by NRTL, at 1 atm. Water 0.5
// is outside the split ChemSep's NRTL has (0.6008 to 0.9944 water at
// 298.15 K), so one liquid is right there, as reflux finds, and as step 0
// found DWSIM giving and took for a miss. At 0.8 water, inside it, DWSIM's
// default (NestedLoops) finds reflux's two liquids; InsideOut adds a vapour
// below the heteroazeotrope's T, where the liquids' bubble pressure is below
// 1 atm (0.79 atm at 360 K); GibbsMinimization misses the split at 298.15 and
// 330 K. Its three-phase stability test changes none of it.
fn test_dwsim_s_three_phase_flash_by_approach() {
	f := fixture()!
	fl := db.gamma_phi(['7732-18-5', '71-36-3'])!.fluid
	p := 101325.0
	for t in [298.15, 330.0, 360.0] {
		one := flash.pt(fl.model, t, p, [0.5, 0.5])!
		assert one.phases.len == 1, '${t}'
	}
	mut worst := 0.0
	mut runs := 0
	for run in f['three_phase']!.as_map()['runs']!.as_array() {
		r := run.as_map()
		t := r['T']!.f64()
		z := floats(r['z']!)
		approach := r['approach']!.str()
		l2 := r['liquid2']!.f64()
		x1 := floats(r['x1']!)
		x2 := if l2 > 0.0 { floats(r['x2']!) } else { x1 }
		split := l2 > 1e-6 && math.abs(x1[0] - x2[0]) > 0.01
		runs++
		if z[0] == 0.5 {
			assert !split, '${approach} ${t}: ${r}'
			continue
		}
		own := flash.pt(fl.model, t, p, z)!
		assert own.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2], '${t}'
		match approach {
			'NestedLoops' {
				assert split && r['beta']!.f64() == 0.0, '${t}'
				lo, hi := if x1[0] < x2[0] { x1[0], x2[0] } else { x2[0], x1[0] }
				worst = math.max(worst, math.abs(lo - own.phases[0].state.n[0]))
				worst = math.max(worst, math.abs(hi - own.phases[1].state.n[0]))
			}
			'InsideOut' {
				assert split && r['beta']!.f64() > 1e-4, '${t}'
				b := flash.bubble_p(fl.model, t, own.phases[0].state.n)!
				assert b.p < p, '${t}: ${b.p}'
			}
			else {
				if t < 340.0 {
					assert !split, '${approach} ${t}'
				}
			}
		}
	}
	assert runs == 36
	println("DWSIM's NestedLoops split against reflux's: ${worst:.1e}")
	// Measured: 1.5e-6, its loops' convergence and its vapour pressures'
	// more figures.
	assert worst <= 1e-5, '${worst}'
}

// A pure liquid's enthalpy less its ideal gas's at 330 K and 1 atm, by each
// of DWSIM's modes, against reflux's from its fugacities (D5): Lee–Kesler's
// corresponding states; Ideal and Excess, minus its ΔHvap correlation to a
// few J/mol; ExpData and LiqCp_Excess, its liquid-Cp data. Recorded.
fn test_dwsim_s_enthalpy_modes() {
	f := fixture()!
	fl := db.gamma_phi(cas)!.fluid
	mut own := []f64{}
	for i in 0 .. 2 {
		n := [0, 1].map(if it == i { 1.0 } else { 0.0 })
		st := fl.model.phase(330.0, 101325.0, n, .liquid)!
		own << fl.model.phase_residual(st).h
	}
	for row in f['enthalpy']!.as_array() {
		r := row.as_map()
		mode := r['mode']!.str()
		h := floats(r['liquid_less_gas']!)
		hvap := floats(r['hvap']!)
		println('${mode}: ${h.map(math.round(it))} J/mol; reflux ${own.map(math.round(it))}; ΔHvap correlation ${hvap.map(math.round(it))}')
		if mode in ['Ideal', 'Excess'] {
			for i in 0 .. 2 {
				assert math.abs(h[i] + hvap[i]) < 5.0, '${mode} ${i}'
			}
		}
	}
	// Reflux's is Clapeyron's on its vapour pressure, about 37.5 and 42.8
	// kJ/mol, between the modes.
	assert own[0] < -37000.0 && own[0] > -38000.0 && own[1] < -42500.0 && own[1] > -43500.0
}
