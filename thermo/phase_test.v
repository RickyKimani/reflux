// Checks phase identification and the Soave α guard against
// testdata/m3_pip.json, written by tools/gen_m3_fixtures.py (run by hand
// with uv) from Caleb Bell's `thermo` 0.6.1: the same PR and SRK, the same
// constants and kᵢⱼ.
import json2
import math
import os
import thermo { Component, CubicKind, Root }

struct Source {
	generated string
	generator string
	python    string
	thermo    string
}

struct FixtureComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
}

// RootValues is one volume root of a state. The Soave-zero states also
// carry what M1's fixtures check.
struct RootValues {
	root   string
	v      f64
	pip    f64
	label  string
	z      f64
	ln_phi []f64
	h_dep  f64
	s_dep  f64
	cp_dep f64
}

struct PipCase {
	eos        string
	components []int
	kij        []f64
	n          []f64
	t          f64
	p          f64
	t0         f64
	roots      []RootValues
}

struct PipFixture {
	source     Source
	note       string
	components []FixtureComponent
	cases      []PipCase
	soave_zero []PipCase
}

fn read_pip_fixture() !PipFixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_pip.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py`.') }
	return json2.decode[PipFixture](text)!
}

fn model_for(f PipFixture, c PipCase) !thermo.Cubic {
	comps := c.components.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if c.eos == 'pr' { CubicKind.pr } else { CubicKind.srk }
	return thermo.new_cubic(kind, comps, kij: c.kij)!
}

fn root_of(name string) Root {
	return if name == 'liquid' { Root.liquid } else { Root.vapor }
}

fn agree(got f64, want f64, tol f64, floor f64) bool {
	return math.abs(got - want) <= tol * math.max(floor, math.abs(want))
}

fn test_pip_and_labels_match_thermo() {
	// 1,327 roots of PR and SRK states: pure components and mixtures, liquid,
	// vapour and supercritical, either side of the saturation curve and the
	// critical point, and liquids squeezed to within 0.02% of the covolume.
	// The volume is M1's 1e-9. Π uses finite differences for its second
	// derivatives (D6), so 1e-7; it agrees to about 5e-9.
	f := read_pip_fixture()!
	mut roots := 0
	for c in f.cases {
		m := model_for(f, c)!
		for r in c.roots {
			at := '${c.eos} ${c.components} ${c.n} at ${c.t} K, ${c.p} Pa, ${r.root} root'
			st := thermo.solve_tp(m, c.t, c.p, c.n, root_of(r.root))!
			assert agree(st.molar_volume(), r.v, 1e-9, 0.0), '${at}: V ${st.molar_volume()} vs ${r.v}'
			p := thermo.pip(m, c.t, st.v, c.n)
			assert agree(p, r.pip, 1e-7, 0.0), '${at}: Π ${p} vs ${r.pip}'
			label := thermo.phase_label(m, st)
			assert label == if r.label == 'liquid' { thermo.PhaseLabel.liquid } else { .vapor }, '${at}: ${label} vs ${r.label}'
			roots++
		}
	}
	assert roots == 1327
}

fn test_soave_zero_matches_thermo() {
	// Pure CO₂ from 0.9 to 1.5 times T₀ = Tc(1 + 1/m)², where Soave's α is
	// zero: 1,775 K with PR, 1,491 K with SRK. For a pure component thermo's
	// α has no division in it, so it's a valid reference there, and reflux
	// must match it to M1's 1e-9 right through the zero.
	f := read_pip_fixture()!
	assert f.soave_zero.len == 48
	for c in f.soave_zero {
		m := model_for(f, c)!
		rt := thermo.r_gas * c.t
		for r in c.roots {
			at := '${c.eos} CO₂ at ${c.t / c.t0} T₀, ${c.p} Pa, ${r.root} root'
			st := thermo.solve_tp(m, c.t, c.p, c.n, root_of(r.root))!
			rp := thermo.residual_props(m, st)
			assert agree(st.z, r.z, 1e-9, 0.0), '${at}: Z ${st.z} vs ${r.z}'
			assert agree(st.ln_phi[0], r.ln_phi[0], 1e-9, 1.0), '${at}: ln φ'
			assert agree(rp.h, r.h_dep, 1e-9, rt), '${at}: Hʳ ${rp.h} vs ${r.h_dep}'
			assert agree(rp.s, r.s_dep, 1e-9, thermo.r_gas), '${at}: Sʳ ${rp.s} vs ${r.s_dep}'
			assert agree(rp.cp, r.cp_dep, 1e-9, thermo.r_gas), '${at}: Cpʳ ${rp.cp} vs ${r.cp_dep}'
			assert agree(thermo.pip(m, c.t, st.v, c.n), r.pip, 1e-7, 0.0), '${at}: Π'
		}
	}
}

fn test_a_mixture_stays_smooth_through_the_soave_zero() {
	// m = 0.5 puts the zero exactly at Tr = 9, where the old √(a_c·α), with
	// derivatives that divided by it, gave 0/0. Through it, F and its
	// T-derivatives must be finite and consistent with each other: the
	// central differences of F and of F_T match F_T and F_TT. Taking
	// √α = |k| instead would make F_T jump here.
	methane := Component{
		name:  'methane'
		tc:    190.6
		pc:    45.99e5
		omega: 0.012
		mw:    0.016043
	}
	heavy := Component{
		name:  'heavy'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	m := thermo.new_cubic(.pr, [methane, heavy],
		kij:    [0.0, 0.02, 0.02, 0.0]
		alphas: [thermo.soave_pr(methane.omega), thermo.Alpha(thermo.Soave{
			m: 0.5
		})]
	)!
	t0 := 9.0 * heavy.tc
	n := [0.5, 0.5]
	v := 1e-3
	r := m.residual2(t0, v, n)
	for x in [r.f, r.f_v, r.f_t, r.f_tt, r.f_tv] {
		assert math.is_finite(x)
	}
	h := 1e-4 * t0
	up := m.residual2(t0 + h, v, n)
	dn := m.residual2(t0 - h, v, n)
	assert agree((up.f - dn.f) / (2.0 * h), r.f_t, 1e-7, 0.0)
	assert agree((up.f_t - dn.f_t) / (2.0 * h), r.f_tt, 1e-7, 0.0)
	st := thermo.solve_tp(m, t0, 1e6, n, .vapor)!
	rp := thermo.residual_props(m, st)
	assert math.is_finite(rp.h) && math.is_finite(rp.cp)
}

fn test_the_volume_solver_reaches_a_critical_point() {
	// At (Tc, Pc) the cubic in Z is (Z − Zc)³: a triple root, where Newton
	// slows to linear convergence and rounding limits V to about ε^(1/3).
	// Each branch must still find its root, the liquid's from below Vc and
	// the vapour's from above, within 1e-4. They must differ: solve_tp falls
	// back to the other branch when one fails, and that would return the
	// same V twice. Zc is (1 − Ωb)/3 for PR and 1/3 for SRK.
	f := read_pip_fixture()!
	zc_pr := (1.0 - 0.07779607390388844) / 3.0
	for kind in [CubicKind.pr, .srk] {
		zc := if kind == .pr { zc_pr } else { 1.0 / 3.0 }
		for fc in f.components {
			c := Component{
				name:  fc.name
				tc:    fc.tc
				pc:    fc.pc
				omega: fc.omega
				mw:    fc.mw
			}
			m := thermo.new_cubic(kind, [c])!
			vc := zc * thermo.r_gas * c.tc / c.pc
			liquid := thermo.solve_tp(m, c.tc, c.pc, [1.0], .liquid)!
			vapor := thermo.solve_tp(m, c.tc, c.pc, [1.0], .vapor)!
			assert agree(liquid.v, vc, 1e-4, 0.0), '${kind} ${c.name}: liquid ${liquid.v} vs Vc ${vc}'
			assert agree(vapor.v, vc, 1e-4, 0.0), '${kind} ${c.name}: vapour ${vapor.v} vs Vc ${vc}'
			assert liquid.v < vc && vc < vapor.v, '${kind} ${c.name}: ${liquid.v} and ${vapor.v} around Vc ${vc}'
		}
	}
}

fn test_the_volume_solver_works_around_a_critical_point() {
	// Within 1e-6 of Tc and 1e-8 of Pc, either side of both, dP/dV at the
	// root is tiny and V is poorly determined. Before the rounding floor
	// (PDpdv.noise) the solver failed on 276 of these 528 requests.
	f := read_pip_fixture()!
	mut solved := 0
	for kind in [CubicKind.pr, .srk] {
		for fc in f.components[1..3] {
			c := Component{
				name:  fc.name
				tc:    fc.tc
				pc:    fc.pc
				omega: fc.omega
				mw:    fc.mw
			}
			m := thermo.new_cubic(kind, [c])!
			for dt in [0.0, 1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, -1e-12, -1e-10, -1e-9, -1e-8] {
				for dp in [0.0, 1e-12, 1e-10, 1e-8, -1e-10, -1e-8] {
					for root in [Root.liquid, .vapor] {
						t, p := c.tc * (1.0 + dt), c.pc * (1.0 + dp)
						st := thermo.solve_tp(m, t, p, [1.0], root) or {
							assert false, '${kind} ${c.name} at Tc(1 + ${dt}), Pc(1 + ${dp}), ${root}: ${err}'
							continue
						}
						assert agree(thermo.pressure(m, t, st.v, [1.0]), p, 1e-12, 0.0)
						solved++
					}
				}
			}
		}
	}
	assert solved == 528
}

fn test_phase_label_on_textbook_states() {
	// n-butane as a cold liquid, a warm vapour, and a
	// supercritical fluid that is gas-like at low density and liquid-like
	// at high density.
	butane := Component{
		name:  'n-butane'
		tc:    425.1
		pc:    37.96e5
		omega: 0.2
		mw:    0.058123
	}
	m := thermo.new_cubic(.pr, [butane])!
	cases := [
		[250.0, 1e5, 0.0], // liquid
		[400.0, 1e5, 1.0], // vapour
		[600.0, 1e6, 1.0], // supercritical, dilute
		[430.0, 3e7, 0.0], // supercritical, dense
	]
	for c in cases {
		st := thermo.solve_tp(m, c[0], c[1], [1.0], .stable)!
		want := if c[2] == 0.0 { thermo.PhaseLabel.liquid } else { .vapor }
		assert thermo.phase_label(m, st) == want, '${c[0]} K, ${c[1]} Pa'
	}
}
