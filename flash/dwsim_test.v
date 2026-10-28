// Checks reflux against DWSIM 9's Peng–Robinson (exit criterion 6), using
// testdata/m3_dwsim.json, written by hand on Windows by
// tools/gen_m3_dwsim.py. Each case has DWSIM's answer, the exact answer
// for the same PR with DWSIM's constants and kᵢⱼ (thermo, and mpmath for
// the ideal-gas Cp), and the PR with DWSIM's rounded constants
// (Ωa 0.45724, Ωb 0.0778, R 8.314), which explains the difference.
//
// reflux must match the exact answer as tightly as it matches thermo
// elsewhere. Its distance from DWSIM is then the exact answer's, which the
// fixture records, so it can't grow unnoticed.
import flash
import json2
import math
import os
import tdep
import thermo { Component }

struct DwComponent {
	name    string
	tc      f64
	pc      f64
	omega   f64
	mw      f64
	cp      []f64
	psat_eq int
	psat    []f64
}

struct DwSystem {
	name       string
	components []DwComponent
	kij        []f64
	z          []f64
}

// DwState is a split as DWSIM reports one: β, y and x, with β 1 or 0 and
// zeros for the absent phase when there's one phase.
struct DwState {
	t        f64
	p        f64
	beta     f64
	y        []f64
	x        []f64
	h        f64
	s        f64
	residual f64
}

struct BoundaryAnswer {
	t         f64
	p         f64
	incipient []f64
	residual  f64
	error     string
}

struct ByState {
	state   DwState
	spec    f64
	landing DwState
}

struct DwCase {
	system         int
	kind           string
	t              f64
	p              f64
	vf             f64
	dwsim          DwState
	exact          DwState
	rounded        DwState
	s_mixing_error f64
}

struct DwByStateCase {
	system int
	kind   string
	t      f64
	p      f64
	vf     f64
	dwsim  ByState
	exact  ByState
}

struct DwBoundaryCase {
	system  int
	kind    string
	t       f64
	p       f64
	dwsim   BoundaryAnswer
	exact   BoundaryAnswer
	rounded BoundaryAnswer
}

struct DwFixture {
	t_ref   f64
	p_ref_s f64
	systems []DwSystem
}

fn read_text() !string {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_dwsim.json')
	return os.read_file(path) or {
		error('${path}: ${err}. Regenerate on Windows with DWSIM 9: `uv run --python 3.13 tools/gen_m3_dwsim.py`.')
	}
}

// read_cases decodes the cases of one kind. The kinds have different shapes,
// so each decodes the whole list as its own type and keeps its own.
fn read_cases[T](kinds []string) ![]T {
	all := json2.decode[DwCases[T]](read_text()!)!
	return all.cases.filter(it.kind in kinds)
}

struct DwCases[T] {
	cases []T
}

fn cubic_of(s DwSystem) !thermo.Cubic {
	comps := s.components.map(Component{
		name:  it.name
		tc:    it.tc
		pc:    it.pc
		omega: it.omega
		mw:    it.mw
	})
	return thermo.new_cubic(.pr, comps, kij: s.kij)!
}

// fluid_of adds DWSIM's ideal-gas Cp, ChemSep's eq. 16 in J/(kmol·K), with
// no formation data: H is then relative to the ideal gas at 298.15 K, as
// DWSIM's is.
fn fluid_of(s DwSystem, m thermo.Cubic) !thermo.Fluid {
	return thermo.new_fluid(m, s.components.map(thermo.IdealGas{
		cp: thermo.Correlated{
			c:     tdep.Correlation{
				eq:   16
				a:    it.cp[0]
				b:    it.cp[1]
				c:    it.cp[2]
				d:    it.cp[3]
				e:    it.cp[4]
				tmax: math.inf(1)
			}
			scale: 1e-3
		}
	}))!
}

// Molar is one mole's H and S, with S on DWSIM's reference pressure,
// 101325 Pa (reflux's is 1 bar).
struct Molar {
	h f64
	s f64
}

fn molar_of(fl thermo.Fluid, p_ref_s f64, r flash.FlashResult) Molar {
	mut h, mut s := 0.0, 0.0
	for ph in r.phases {
		pr := fl.props(ph.state)
		h += ph.fraction * pr.h
		s += ph.fraction * (pr.s - thermo.r_gas * math.log(thermo.p_ref / p_ref_s))
	}
	return Molar{h, s}
}

// split_of puts reflux's result in DWSIM's form.
fn split_of(r flash.FlashResult) DwState {
	nc := r.z.len
	if r.phases.len == 2 {
		return DwState{
			beta: r.beta
			y:    r.phases[0].state.n
			x:    r.phases[1].state.n
		}
	}
	zero := []f64{len: nc}
	n := r.phases[0].state.n
	return if r.beta == 1.0 {
		DwState{
			beta: 1.0
			y:    n
			x:    zero
		}
	} else {
		DwState{
			beta: 0.0
			y:    zero
			x:    n
		}
	}
}

// split_dev is the largest difference in β, y or x.
fn split_dev(a DwState, b DwState) f64 {
	mut d := math.abs(a.beta - b.beta)
	for i in 0 .. a.y.len {
		d = math.max(d, math.max(math.abs(a.y[i] - b.y[i]), math.abs(a.x[i] - b.x[i])))
	}
	return d
}

fn test_pt_matches_the_exact_pr_and_dwsim() {
	// Against the exact PR: the split to 1e-8, as against thermo (exit
	// criterion 2; measured 1.9e-13), and H and S to 1e-9, relative to
	// |H| + RT and |S| + R (measured 1.9e-11 and 2.6e-11). The H and S
	// include integrating ChemSep's eq. 16 Cp, which tdep does numerically.
	//
	// Against DWSIM: the split within 1e-4 (criterion 6; DWSIM's worst is
	// 3e-5, all but 1e-6 of it the rounded constants), and reflux no further
	// from DWSIM than the exact answer is, give or take 1e-8, in the split, H
	// and S.
	fx := json2.decode[DwFixture](read_text()!)!
	cases := read_cases[DwCase](['pt'])!
	mut worst := map[string]f64{}
	for c in cases {
		at := '${fx.systems[c.system].name} at ${c.t} K, ${c.p} Pa'
		sys := fx.systems[c.system]
		m := cubic_of(sys)!
		fl := fluid_of(sys, m)!
		r := flash.pt(m, c.t, c.p, sys.z)!
		got := split_of(r)
		mol := molar_of(fl, fx.p_ref_s, r)
		hs, ss := math.abs(c.exact.h) + thermo.r_gas * c.t, math.abs(c.exact.s) + thermo.r_gas
		dev := {
			'split': split_dev(got, c.exact)
			'H':     math.abs(mol.h - c.exact.h) / hs
			'S':     math.abs(mol.s - c.exact.s) / ss
		}
		for k, v in dev {
			worst[k] = math.max(worst[k], v)
		}
		assert dev['split'] <= 1e-8, '${at}: split ${dev['split']}'
		assert dev['H'] <= 1e-9 && dev['S'] <= 1e-9, '${at}: ${dev}'
		assert split_dev(got, c.dwsim) <= 1e-4, at
		assert split_dev(got, c.dwsim) <= split_dev(c.exact, c.dwsim) + 1e-8, at
		assert math.abs(mol.h - c.dwsim.h) <= math.abs(c.exact.h - c.dwsim.h) + 1e-8 * hs, at
		assert math.abs(mol.s - c.dwsim.s) <= math.abs(c.exact.s - c.dwsim.s) + 1e-8 * ss, at
	}
	println('against the exact PR: ${worst}')
	assert cases.len == 10
}

fn test_dwsims_differences_are_explained() {
	// No reflux here: this checks the fixture's account of DWSIM.
	// - Its split is its rounded-constant PR's, to its convergence floor:
	//   1e-6 (DWSIM's own fugacity residual, at most 1.4e-6).
	// - Its H is the exact one to 2e-4 of |H| + RT: the rounded constants
	//   move departures by 1e-4, and its ideal-gas integral is good to 1e-5.
	// - Its S is the exact one plus its mixing error, −R Σ zᵢ ln zᵢ (M/Mᵢ − 1),
	//   to 0.02 J/(mol·K): R = 8.314 in −R ln(P/P°) and the rounded
	//   departures make the rest.
	fx := json2.decode[DwFixture](read_text()!)!
	for c in read_cases[DwCase](['pt'])! {
		at := '${fx.systems[c.system].name} at ${c.t} K, ${c.p} Pa'
		assert split_dev(c.dwsim, c.rounded) <= 1e-6, at
		assert c.dwsim.residual <= 1.5e-6, at
		assert math.abs(c.dwsim.h - c.exact.h) <= 2e-4 * (math.abs(c.exact.h) + thermo.r_gas * c.t), at
		assert math.abs(c.dwsim.s - c.exact.s - c.s_mixing_error) <= 0.02, at
	}
}

fn test_ph_and_ps_by_state() {
	// Each program's own H or S at the same state, and where DWSIM's flash
	// from it lands. reflux's H and S at the state match the exact ones to
	// 1e-9, so step 7's flashes start from the right number; DWSIM lands
	// within 1e-7 K of its own state, and of the exact one within criterion
	// 6's 0.01 K, except on the pure component's saturation curve (below).
	fx := json2.decode[DwFixture](read_text()!)!
	cases := read_cases[DwByStateCase](['ph', 'ps'])!
	for c in cases {
		sys := fx.systems[c.system]
		at := '${c.kind} ${sys.name} at ${c.t} K, ${c.p} Pa, vf ${c.vf}'
		m := cubic_of(sys)!
		fl := fluid_of(sys, m)!
		st := c.exact.state
		mut mol := Molar{}
		if sys.z.len == 1 {
			// The saturated state: reflux's two roots at the exact Tsat.
			liq := thermo.solve_tp(m, st.t, c.p, [1.0], .liquid)!
			vap := thermo.solve_tp(m, st.t, c.p, [1.0], .vapor)!
			assert math.abs(liq.ln_phi[0] - vap.ln_phi[0]) <= 1e-10, at
			l, g := fl.props(liq), fl.props(vap)
			shift := thermo.r_gas * math.log(thermo.p_ref / fx.p_ref_s)
			mol = Molar{(1.0 - c.vf) * l.h + c.vf * g.h, (1.0 - c.vf) * l.s + c.vf * g.s - shift}
		} else {
			mol = molar_of(fl, fx.p_ref_s, flash.pt(m, st.t, c.p, sys.z)!)
		}
		want := c.exact.spec
		scale := if c.kind == 'ph' { thermo.r_gas * st.t } else { thermo.r_gas }
		got := if c.kind == 'ph' { mol.h } else { mol.s }
		assert math.abs(got - want) <= 1e-9 * (math.abs(want) + scale), at
		land := c.dwsim.landing
		assert math.abs(land.t - c.dwsim.state.t) <= 1e-7, '${at}: DWSIM lands off its own state'
		if sys.z.len > 1 {
			assert math.abs(land.t - st.t) <= 0.01, at
			assert math.abs(land.beta - st.beta) <= 1e-4, at
		}
	}
	assert cases.len == 9
}

fn test_dwsim_saturates_a_pure_component_on_its_correlation() {
	// For a pure component, DWSIM's PR package doesn't use the EOS's
	// saturation: its Tsat is where its vapour-pressure correlation (ChemSep
	// eq. 101) gives P. So it sits 0.118 K below PR's Tsat for n-butane at
	// 5 bar, which criterion 6's 0.01 K can't absorb, and its own PR
	// fugacities don't balance there. Recorded, and it can't grow.
	fx := json2.decode[DwFixture](read_text()!)!
	for c in read_cases[DwByStateCase](['ph', 'ps'])! {
		sys := fx.systems[c.system]
		if sys.z.len > 1 {
			continue
		}
		k := sys.components[0]
		pv := tdep.Correlation{
			eq:   k.psat_eq
			a:    k.psat[0]
			b:    k.psat[1]
			c:    k.psat[2]
			d:    k.psat[3]
			e:    k.psat[4]
			tmax: math.inf(1)
		}
		t_dw := c.dwsim.state.t
		assert math.abs(pv.eval(t_dw) - c.p) <= 1e-6 * c.p, 'DWSIM Tsat ${t_dw} K is on its correlation'
		gap := c.exact.state.t - t_dw
		assert gap > 0.1 && gap <= 0.12, 'PR Tsat minus DWSIM Tsat: ${gap} K'
	}
}

fn test_bubble_and_dew_points() {
	// The exact bubble and dew points satisfy reflux's equilibrium condition
	// (criterion 3's 1e-10). DWSIM's are within 5e-5 K and 1e-6 in P of its
	// rounded-constant PR's (measured 3.2e-5 K and 8.2e-7, its convergence
	// floor), and of the exact PR's within criterion 6's
	// 0.01 K and 0.01% in P, but one: methane/n-butane's dew pressure at
	// 300 K, 1.8e-4 away, all of it the rounded constants. Recorded, with a
	// bound that can't grow.
	fx := json2.decode[DwFixture](read_text()!)!
	cases := read_cases[DwBoundaryCase](['bubble', 'dew'])!
	mut over := 0
	for c in cases {
		sys := fx.systems[c.system]
		at := '${c.kind} ${sys.name} at ${c.t} K, ${c.p} Pa'
		m := cubic_of(sys)!
		e := c.exact
		feed_root, new_root := if c.kind == 'bubble' {
			thermo.Root.liquid, thermo.Root.vapor
		} else {
			thermo.Root.vapor, thermo.Root.liquid
		}
		feed := thermo.solve_tp(m, e.t, e.p, sys.z, feed_root)!
		new := thermo.solve_tp(m, e.t, e.p, e.incipient, new_root)!
		for i, w in e.incipient {
			res := math.abs(math.log(sys.z[i]) + feed.ln_phi[i] - math.log(w) - new.ln_phi[i])
			assert res <= 1e-10, '${at}: component ${i}: ${res}'
		}
		d, r := c.dwsim, c.rounded
		assert d.error == '', at
		assert math.abs(d.t - r.t) <= 5e-5 && math.abs(d.p - r.p) <= 1e-6 * r.p, at
		dt, dp := math.abs(d.t - e.t), math.abs(d.p - e.p) / e.p
		assert dt <= 0.01, '${at}: ${dt} K'
		if dp > 1e-4 {
			assert dp <= 1.8e-4, '${at}: P ${dp}'
			over++
		}
	}
	assert cases.len == 6 && over == 1
}

fn test_reflux_bubble_and_dew_points() {
	// reflux's own bubble and dew points, from the same specification: the
	// exact PR's to 1e-10 relative in T or P, and so no further from DWSIM's
	// than the exact answer is, give or take 1e-8 (exit criterion 6, as
	// amended: within 0.01 K and 0.01% in P, but for the dew pressure that
	// DWSIM's rounded constants move by 1.8e-4).
	fx := json2.decode[DwFixture](read_text()!)!
	for c in read_cases[DwBoundaryCase](['bubble', 'dew'])! {
		sys := fx.systems[c.system]
		at := '${c.kind} ${sys.name} at ${c.t} K, ${c.p} Pa'
		m := cubic_of(sys)!
		beta := if c.kind == 'bubble' { 0.0 } else { 1.0 }
		e, d := c.exact, c.dwsim
		if c.t > 0.0 {
			r := flash.tvf(m, c.t, beta, sys.z)!
			assert math.abs(r.p - e.p) <= 1e-10 * e.p, '${at}: ${r.p} vs ${e.p}'
			assert math.abs(r.p - d.p) <= math.abs(e.p - d.p) + 1e-8 * e.p, at
		} else {
			r := flash.pvf(m, c.p, beta, sys.z)!
			assert math.abs(r.t - e.t) <= 1e-10 * e.t, '${at}: ${r.t} vs ${e.t}'
			assert math.abs(r.t - d.t) <= math.abs(e.t - d.t) + 1e-8 * e.t, at
			assert math.abs(r.t - d.t) <= 0.01, at
		}
	}
}

fn test_reflux_ph_and_ps_by_state() {
	// reflux's own PH and PS flashes from its own H or S at each state (the
	// exact one, to 1e-9: test_ph_and_ps_by_state) land back on the state:
	// T to 1e-8 relative and β to 1e-8. So for a mixture reflux sits no
	// further from DWSIM's landing than the exact state does, give or take
	// 1e-6 K: within criterion 6's 0.01 K. For pure n-butane DWSIM lands on
	// its own saturation temperature, 0.118 K away, as explained above.
	fx := json2.decode[DwFixture](read_text()!)!
	for c in read_cases[DwByStateCase](['ph', 'ps'])! {
		sys := fx.systems[c.system]
		at := '${c.kind} ${sys.name} at ${c.t} K, ${c.p} Pa, vf ${c.vf}'
		m := cubic_of(sys)!
		fl := fluid_of(sys, m)!
		st := c.exact.state
		r := if c.kind == 'ph' {
			flash.ph(fl, c.p, c.exact.spec, sys.z)!
		} else {
			flash.ps(fl, c.p, c.exact.spec + thermo.r_gas * math.log(thermo.p_ref / fx.p_ref_s), sys.z)!
		}
		assert math.abs(r.t - st.t) <= 1e-8 * st.t, '${at}: ${r.t} vs ${st.t}'
		assert math.abs(r.beta - st.beta) <= 1e-8, '${at}: β ${r.beta} vs ${st.beta}'
		land := c.dwsim.landing
		if sys.z.len > 1 {
			assert math.abs(r.t - land.t) <= math.abs(st.t - land.t) + 1e-6, at
			assert math.abs(r.t - land.t) <= 0.01, at
		} else {
			gap := r.t - land.t
			assert gap > 0.1 && gap <= 0.12, '${at}: ${gap} K'
		}
	}
}
