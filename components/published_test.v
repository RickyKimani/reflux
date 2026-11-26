// Published data (M6 step 10, criterion 8): the activity models on ChemSep's
// parameters against data sets from the ThermoML archive
// (testdata/m6_published.json, from tools/gen_m6_published.py). Each VLE set
// is first given the point test for Gibbs–Duhem consistency, here by
// fit.point_test and there by SciPy, and one that fails is recorded so and
// not used. Then each model's deviations from the rest, recorded, which
// can't grow (below, `recorded`).
import components
import fit
import flash
import json2
import math
import os
import thermo

struct VlePoint {
	x1 string
	t  string @[json: 'T']
	y1 string
}

struct PyPointTest {
	terms        int
	points       int
	coefficients []f64
	mean_dy      f64
	max_dy       f64
	consistent   bool
}

struct VleSet {
	system     []string
	cas        []string
	doi        string
	p_kpa      string @[json: 'P_kPa']
	points     []VlePoint
	point_test PyPointTest
}

struct LlePoint {
	doi     string
	t       string @[json: 'T']
	p_kpa   string @[json: 'P_kPa']
	organic map[string]string // water's or 1-butanol's fraction
	aqueous map[string]string
}

struct Vlle {
	doi     string
	p_kpa   string @[json: 'P_kPa']
	t       string @[json: 'T']
	organic map[string]string
	aqueous map[string]string
	vapour  map[string]string
}

struct TieLine {
	organic map[string]string
	aqueous map[string]string
}

struct TernaryLle {
	system    []string
	doi       string
	t         string @[json: 'T']
	p_kpa     string @[json: 'P_kPa']
	tie_lines []TieLine
}

struct Published {
	vle         []VleSet
	lle         []LlePoint
	vlle        []Vlle
	ternary_lle []TernaryLle @[json: 'ternary_lle']
}

const db = components.load_default() or { panic(err) }

fn published() !Published {
	return json2.decode[Published](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_published.json'))!)!
}

// The deviations, recorded: each the measured value rounded up to two
// figures, so that a change that makes a model worse fails here. Keyed by
// what, the model and the statistic.
const recorded = {
	'methanol-water nrtl mean |ΔT|':                0.38
	'methanol-water nrtl max |ΔT|':                 1.8
	'methanol-water nrtl mean |Δy|':                0.011
	'methanol-water nrtl max |Δy|':                 0.041
	'methanol-water uniquac mean |ΔT|':             0.41
	'methanol-water uniquac max |ΔT|':              2.0
	'methanol-water uniquac mean |Δy|':             0.011
	'methanol-water uniquac max |Δy|':              0.047
	'methanol-water wilson mean |ΔT|':              0.49
	'methanol-water wilson max |ΔT|':               2.4
	'methanol-water wilson mean |Δy|':              0.013
	'methanol-water wilson max |Δy|':               0.056
	'methanol-water unifac mean |ΔT|':              0.42
	'methanol-water unifac max |ΔT|':               1.8
	'methanol-water unifac mean |Δy|':              0.012
	'methanol-water unifac max |Δy|':               0.043
	'methanol-water unifac_dortmund mean |ΔT|':     0.41
	'methanol-water unifac_dortmund max |ΔT|':      1.9
	'methanol-water unifac_dortmund mean |Δy|':     0.011
	'methanol-water unifac_dortmund max |Δy|':      0.045
	'ethanol-water nrtl mean |ΔT|':                 0.21
	'ethanol-water nrtl max |ΔT|':                  0.78
	'ethanol-water nrtl mean |Δy|':                 0.0044
	'ethanol-water nrtl max |Δy|':                  0.026
	'ethanol-water uniquac mean |ΔT|':              0.18
	'ethanol-water uniquac max |ΔT|':               0.82
	'ethanol-water uniquac mean |Δy|':              0.0044
	'ethanol-water uniquac max |Δy|':               0.024
	'ethanol-water wilson mean |ΔT|':               0.17
	'ethanol-water wilson max |ΔT|':                0.67
	'ethanol-water wilson mean |Δy|':               0.0051
	'ethanol-water wilson max |Δy|':                0.021
	'ethanol-water unifac mean |ΔT|':               0.28
	'ethanol-water unifac max |ΔT|':                1.3
	'ethanol-water unifac mean |Δy|':               0.0071
	'ethanol-water unifac max |Δy|':                0.033
	'ethanol-water unifac_dortmund mean |ΔT|':      0.26
	'ethanol-water unifac_dortmund max |ΔT|':       0.55
	'ethanol-water unifac_dortmund mean |Δy|':      0.0044
	'ethanol-water unifac_dortmund max |Δy|':       0.019
	'water-1-butanol LLE nrtl max |Δx|':            0.094
	'water-1-butanol LLE uniquac max |Δx|':         0.12
	'water-1-butanol LLE unifac max |Δx|':          0.016
	'water-1-butanol LLE unifac_lle max |Δx|':      0.032
	'water-1-butanol LLE unifac_dortmund max |Δx|': 0.15
	'heteroazeotrope nrtl |ΔT|':                    0.34
	'heteroazeotrope nrtl organic |Δx|':            0.05
	'heteroazeotrope nrtl aqueous |Δx|':            0.0031
	'heteroazeotrope nrtl vapour |Δy|':             0.013
	'heteroazeotrope uniquac |ΔT|':                 0.19
	'heteroazeotrope uniquac organic |Δx|':         0.18
	'heteroazeotrope uniquac aqueous |Δx|':         0.0012
	'heteroazeotrope uniquac vapour |Δy|':          0.0065
	'heteroazeotrope unifac |ΔT|':                  0.39
	'heteroazeotrope unifac organic |Δx|':          0.089
	'heteroazeotrope unifac aqueous |Δx|':          0.011
	'heteroazeotrope unifac vapour |Δy|':           0.008
	'heteroazeotrope unifac_lle |ΔT|':              0.81
	'heteroazeotrope unifac_lle organic |Δx|':      0.074
	'heteroazeotrope unifac_lle aqueous |Δx|':      0.0096
	'heteroazeotrope unifac_lle vapour |Δy|':       0.021
	'heteroazeotrope unifac_dortmund |ΔT|':         0.11
	'heteroazeotrope unifac_dortmund organic |Δx|': 0.18
	'heteroazeotrope unifac_dortmund aqueous |Δx|': 0.0084
	'heteroazeotrope unifac_dortmund vapour |Δy|':  6.9e-05
	'ternary LLE unifac splits':                    0.0
	'ternary LLE unifac mean |Δx|':                 0.013
	'ternary LLE unifac max |Δx|':                  0.065
	'ternary LLE unifac_lle splits':                0.0
	'ternary LLE unifac_lle mean |Δx|':             0.0074
	'ternary LLE unifac_lle max |Δx|':              0.033
	'ternary LLE unifac_dortmund splits':           0.0
	'ternary LLE unifac_dortmund mean |Δx|':        0.04
	'ternary LLE unifac_dortmund max |Δx|':         0.15
}

// within asserts a deviation is no larger than recorded, and prints it.
fn within(key string, value f64) {
	println('${key}: ${value:.3e}')
	if bound := recorded[key] {
		assert value <= bound, '${key}: ${value} against ${bound} recorded'
	}
}

fn psat(name string, t f64) f64 {
	p := db.find(name) or { panic(err) }.property('vapor_pressure') or { panic(err) }
	return p.c.eval(t) * p.scale
}

// fit.point_test against SciPy's, set by set: the same Legendre coefficients
// to 1e-8 of their scale and the same mean and largest |Δy| to 1e-9, and the
// same verdict. Both meet their gradient tests on one objective, reflux's
// at 1e-10, which leaves θ that close to the minimum.
fn test_the_point_test_against_scipy() {
	data := published()!
	mut consistent := 0
	mut worst_coefficient, mut worst_dy := 0.0, 0.0
	for s in data.vle {
		inside := s.points.filter(it.x1.f64() > 0.0 && it.x1.f64() < 1.0)
		p := s.p_kpa.f64() * 1e3
		pt := fit.point_test(fit.PointData{
			x1:    inside.map(it.x1.f64())
			y1:    inside.map(it.y1.f64())
			p:     []f64{len: inside.len, init: p}
			psat1: inside.map(psat(s.cas[0], it.t.f64()))
			psat2: inside.map(psat(s.cas[1], it.t.f64()))
		})!
		want := s.point_test
		assert want.points == inside.len && want.terms == fit.point_test_terms
		mut scale := 0.0
		for a in want.coefficients {
			scale = math.max(scale, math.abs(a))
		}
		for k, a in want.coefficients {
			worst_coefficient = math.max(worst_coefficient, math.abs(pt.coefficients[k] - a) / scale)
		}
		worst_dy = math.max(worst_dy, math.abs(pt.mean_dy - want.mean_dy))
		worst_dy = math.max(worst_dy, math.abs(pt.max_dy - want.max_dy))
		assert pt.consistent == want.consistent, s.doi
		if pt.consistent {
			consistent++
		}
	}
	println('point test against SciPy: coefficients ${worst_coefficient:.1e}, |Δy| ${worst_dy:.1e}')
	// Measured: 6.1e-9 and 3.9e-10, the gradient test's 1e-10 through the
	// Jacobian's condition, up to 1.3e3.
	assert worst_coefficient <= 3e-8 && worst_dy <= 1e-9, '${worst_coefficient} ${worst_dy}'
	// Measured: 9 of the 11 consistent; the 2 that aren't have mean |Δy| of
	// 0.011 and 0.022. One more, ethanol–water's azeotropic region alone,
	// has two minima, 0.0025 and 0.050, and is consistent by the lower.
	assert consistent == 9
}

// The VLE models against the consistent sets at 101.3 kPa: each point's
// bubble T and vapour by the model at the measured P and x, the vapour the
// ideal gas and the Poynting correction on (D4's defaults).
fn test_vle_against_the_data() {
	data := published()!
	for system in ['methanol', 'ethanol'] {
		sets := data.vle.filter(it.system[0] == system && it.point_test.consistent)
		assert sets.len == if system == 'methanol' { 4 } else { 5 }
		for model in [components.LiquidModel.nrtl, .uniquac, .wilson, .unifac, .unifac_dortmund] {
			fl := db.gamma_phi(sets[0].cas, liquid: model)!.fluid
			mut sum_t, mut worst_t, mut sum_y, mut worst_y := 0.0, 0.0, 0.0, 0.0
			mut n := 0
			for s in sets {
				p := s.p_kpa.f64() * 1e3
				for pt in s.points {
					x1 := pt.x1.f64()
					b := flash.bubble_t(fl.model, p, [x1, 1.0 - x1])!
					dt := math.abs(b.t - pt.t.f64())
					vapour := b.phases.filter(it.label == .vapor)[0]
					dy := math.abs(vapour.state.n[0] - pt.y1.f64())
					sum_t += dt
					sum_y += dy
					worst_t = math.max(worst_t, dt)
					worst_y = math.max(worst_y, dy)
					n++
				}
			}
			within('${system}-water ${model} mean |ΔT|', sum_t / f64(n))
			within('${system}-water ${model} max |ΔT|', worst_t)
			within('${system}-water ${model} mean |Δy|', sum_y / f64(n))
			within('${system}-water ${model} max |Δy|', worst_y)
		}
	}
}

// Water and 1-butanol's mutual solubilities: a PT flash of the two
// measured liquids' mean, at the measured T and P, against each liquid's
// water fraction. Every model splits it.
fn test_lle_against_the_data() {
	data := published()!
	names := ['water', '1-butanol']
	for model in [components.LiquidModel.nrtl, .uniquac, .unifac, .unifac_lle, .unifac_dortmund] {
		fl := db.gamma_phi(names, liquid: model)!.fluid
		mut worst := 0.0
		for pt in data.lle {
			organic := water_of(pt.organic)
			aqueous := water_of(pt.aqueous)
			z := (organic + aqueous) / 2.0
			r := flash.pt(fl.model, pt.t.f64(), pt.p_kpa.f64() * 1e3, [z, 1.0 - z])!
			assert r.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2], '${model} ${pt.t}'
			worst = math.max(worst, math.abs(r.phases[0].state.n[0] - organic))
			worst = math.max(worst, math.abs(r.phases[1].state.n[0] - aqueous))
		}
		within('water-1-butanol LLE ${model} max |Δx|', worst)
	}
}

// water_of is a liquid's water fraction, given as water's or 1-butanol's.
fn water_of(m map[string]string) f64 {
	if w := m['water'] {
		return w.f64()
	}
	return 1.0 - m['1-butanol'].f64()
}

// The heteroazeotrope at 101.3 kPa, where a vapour first appears over the
// two liquids: bisection on T for a feed between them, two liquids below
// and a vapour above, to 1e-8 K; the liquids from below and the vapour
// from above.
fn test_the_heteroazeotrope_against_the_data() {
	data := published()!
	names := ['water', '1-butanol']
	v := data.vlle[0]
	p := v.p_kpa.f64() * 1e3
	for model in [components.LiquidModel.nrtl, .uniquac, .unifac, .unifac_lle, .unifac_dortmund] {
		fl := db.gamma_phi(names, liquid: model)!.fluid
		z := [0.8, 0.2]
		mut lo, mut hi := 330.0, 380.0
		assert !flash.pt(fl.model, lo, p, z)!.phases.any(it.label == .vapor)
		assert flash.pt(fl.model, hi, p, z)!.phases.any(it.label == .vapor)
		for hi - lo > 1e-8 {
			mid := (lo + hi) / 2.0
			if flash.pt(fl.model, mid, p, z)!.phases.any(it.label == .vapor) {
				hi = mid
			} else {
				lo = mid
			}
		}
		below := flash.pt(fl.model, lo, p, z)!
		above := flash.pt(fl.model, hi, p, z)!
		assert below.phases.map(it.label) == [thermo.PhaseLabel.liquid, .liquid2], '${model}'
		vapour := above.phases.filter(it.label == .vapor)[0]
		within('heteroazeotrope ${model} |ΔT|', math.abs(lo - v.t.f64()))
		within('heteroazeotrope ${model} organic |Δx|', math.abs(below.phases[0].state.n[1] - v.organic['1-butanol'].f64()))
		within('heteroazeotrope ${model} aqueous |Δx|', math.abs(below.phases[1].state.n[1] - v.aqueous['1-butanol'].f64()))
		within('heteroazeotrope ${model} vapour |Δy|', math.abs(vapour.state.n[1] - v.vapour['1-butanol'].f64()))
	}
}

// Water, 1-butanol and p-xylene's tie lines at 313.15 K: a PT flash of each
// tie line's mean against both liquids' water and 1-butanol fractions. The
// banks pair neither alcohol nor water with p-xylene, so NRTL and UNIQUAC
// are refused, saying so, and the UNIFACs are compared.
fn test_a_ternary_lle_against_the_data() {
	data := published()!
	tern := data.ternary_lle[0]
	names := ['water', '1-butanol', 'p-xylene']
	for model in [components.LiquidModel.nrtl, .uniquac] {
		if _ := db.gamma_phi(names, liquid: model) {
			assert false, '${model} has every pair'
		} else {
			assert err.msg().contains("doesn't list"), err.msg()
		}
	}
	t := tern.t.f64()
	p := tern.p_kpa.f64() * 1e3
	for model in [components.LiquidModel.unifac, .unifac_lle, .unifac_dortmund] {
		fl := db.gamma_phi(names, liquid: model)!.fluid
		mut worst, mut sum := 0.0, 0.0
		mut n, mut splits := 0, 0
		for line in tern.tie_lines {
			o := [line.organic['water'].f64(), line.organic['1-butanol'].f64()]
			a := [line.aqueous['water'].f64(), line.aqueous['1-butanol'].f64()]
			z := [(o[0] + a[0]) / 2.0, (o[1] + a[1]) / 2.0, 1.0 - (o[0] + a[0] + o[1] + a[1]) / 2.0]
			r := flash.pt(fl.model, t, p, z)!
			if r.phases.len != 2 {
				continue
			}
			splits++
			// The organic liquid is the lighter, with the p-xylene.
			for i in 0 .. 2 {
				for d in [r.phases[0].state.n[i] - o[i], r.phases[1].state.n[i] - a[i]] {
					worst = math.max(worst, math.abs(d))
					sum += math.abs(d)
					n++
				}
			}
		}
		within('ternary LLE ${model} splits', 11.0 - f64(splits))
		within('ternary LLE ${model} mean |Δx|', sum / f64(n))
		within('ternary LLE ${model} max |Δx|', worst)
	}
}
