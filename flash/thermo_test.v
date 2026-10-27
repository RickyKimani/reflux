// Checks reflux against the flash references in testdata/m3_thermo.json,
// written by tools/gen_m3_fixtures.py from Caleb Bell's `thermo` 0.6.1 with
// the same PR or SRK, constants, kᵢⱼ and ideal-gas Cp (a quadratic, so
// total H and S compare as exactly as the departures do).
//
// The PT states are what the PH, PS, TV, UV and vapour-fraction flashes of
// steps 6 to 8 will be given and must recover. Here the PT flash reproduces
// them, and every saturation state and bubble or dew point is checked to be
// one with reflux's own model.
import flash
import json2
import math
import os
import thermo { Component, CubicKind }

struct ThermoComponent {
	name  string
	tc    f64
	pc    f64
	omega f64
	mw    f64
	cp    []f64
}

struct ThermoPhase {
	fraction f64
	z        []f64
	v        f64
	h        f64
	s        f64
	label    string
}

// Landing is where thermo's own PH, PS or TV flash went from a PT state's
// H, S or V, or the error it raised.
struct Landing {
	t      f64
	p      f64
	phases int
	error  string
}

struct Landings {
	ph Landing
	ps Landing
	tv Landing
}

struct PtCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	t          f64
	p          f64
	phases     []ThermoPhase
	v          f64
	h          f64
	s          f64
	u          f64
	residual   f64
	thermo     Landings
}

struct SaturationCase {
	components []int
	kij        []f64
	spec       string
	vf         f64
	t          f64
	p          f64
	phases     []ThermoPhase
	v          f64
	h          f64
	s          f64
	u          f64
	residual   f64
}

struct BoundaryCase {
	eos        string
	components []int
	kij        []f64
	z          []f64
	spec       string
	vf         f64
	t          f64
	p          f64
	incipient  []f64
	residual   f64
	v          f64
	h          f64
	s          f64
	error      string
	rejected   string
}

struct ThermoFixture {
	t_ref      f64
	p_ref_s    f64
	components []ThermoComponent
	pt         []PtCase
	saturation []SaturationCase
	boundaries []BoundaryCase
}

fn read_thermo() !ThermoFixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_thermo.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py thermo`.')
	}
	return json2.decode[ThermoFixture](text)!
}

fn cubic_of(f ThermoFixture, eos string, idx []int, kij []f64) !thermo.Cubic {
	comps := idx.map(Component{
		name:  f.components[it].name
		tc:    f.components[it].tc
		pc:    f.components[it].pc
		omega: f.components[it].omega
		mw:    f.components[it].mw
	})
	kind := if eos == 'srk' { CubicKind.srk } else { CubicKind.pr }
	return thermo.new_cubic(kind, comps, kij: kij)!
}

// fluid_of adds the fixture's ideal-gas Cp, cp[0] + cp[1]·T + cp[2]·T², with
// no formation data: H is then relative to the ideal gas at 298.15 K, as
// thermo's is.
fn fluid_of(f ThermoFixture, m thermo.Cubic, idx []int) !thermo.Fluid {
	return thermo.new_fluid(m, idx.map(thermo.IdealGas{
		cp: thermo.Polynomial{
			terms: [
				thermo.PowerTerm{f.components[it].cp[0], 0},
				thermo.PowerTerm{f.components[it].cp[1], 1},
				thermo.PowerTerm{f.components[it].cp[2], 2},
			]
		}
	}))!
}

// Molar is one mole's V, H, S and U, with S on thermo's reference pressure.
struct Molar {
	v f64
	h f64
	s f64
	u f64
}

// molar_of returns a phase's molar V, H, S and U. thermo's S is relative to
// the ideal gas at 101325 Pa, reflux's at P° = 1 bar, so S moves by
// R ln(1 bar / 101325 Pa) between them.
fn molar_of(fl thermo.Fluid, fx ThermoFixture, st thermo.State) Molar {
	pr := fl.props(st)
	return Molar{
		v: st.molar_volume()
		h: pr.h
		s: pr.s - thermo.r_gas * math.log(thermo.p_ref / fx.p_ref_s)
		u: pr.u
	}
}

// Dev is a deviation and where it happened, so that the worst can be shown.
struct Dev {
mut:
	worst f64
	at    string
}

fn (mut d Dev) see(dev f64, at string) {
	if dev > d.worst {
		d.worst = dev
		d.at = at
	}
}

fn test_the_reference_states_match() {
	f := read_thermo()!
	assert f.t_ref == thermo.t_ref
	assert f.pt.len == 83 && f.saturation.len == 36 && f.boundaries.len == 42
}

fn test_pt_matches_thermo() {
	// β and the compositions to 1e-8 (exit criterion 2), molar V, H, S and U
	// to 1e-9.
	f := read_thermo()!
	mut d_split, mut d_v, mut d_h, mut d_s := Dev{}, Dev{}, Dev{}, Dev{}
	mut splits := 0
	for c in f.pt {
		at := '${c.eos} ${c.components} z ${c.z} at ${c.t} K, ${c.p} Pa'
		m := cubic_of(f, c.eos, c.components, c.kij)!
		fl := fluid_of(f, m, c.components)!
		r := flash.pt(m, c.t, c.p, c.z)!
		assert r.phases.len == c.phases.len, '${at}: ${r.phases.len} phases'
		if r.phases.len == 1 {
			assert r.phases[0].label.str() == c.phases[0].label, at
		} else {
			splits++
		}
		mut v, mut h, mut s, mut u := 0.0, 0.0, 0.0, 0.0
		for k, ph in r.phases {
			want := c.phases[k]
			d_split.see(math.abs(ph.fraction - want.fraction), at)
			for i, x in ph.state.n {
				d_split.see(math.abs(x - want.z[i]), at)
			}
			mol := molar_of(fl, f, ph.state)
			d_v.see(math.abs(mol.v - want.v) / want.v, at)
			d_h.see(math.abs(mol.h - want.h) / (math.abs(want.h) + thermo.r_gas * c.t), at)
			d_s.see(math.abs(mol.s - want.s) / (math.abs(want.s) + thermo.r_gas), at)
			v += ph.fraction * mol.v
			h += ph.fraction * mol.h
			s += ph.fraction * mol.s
			u += ph.fraction * mol.u
		}
		d_v.see(math.abs(v - c.v) / c.v, at)
		d_h.see(math.abs(h - c.h) / (math.abs(c.h) + thermo.r_gas * c.t), at)
		d_h.see(math.abs(u - c.u) / (math.abs(c.u) + thermo.r_gas * c.t), at)
		d_s.see(math.abs(s - c.s) / (math.abs(c.s) + thermo.r_gas), at)
		// thermo's own PH, PS and TV flashes, where they ran, land back on
		// this state: so the PT state is the answer steps 6 to 8 must find.
		for l in [c.thermo.ph, c.thermo.ps, c.thermo.tv] {
			if l.error == '' {
				assert math.abs(l.t - c.t) <= 1e-9 * c.t && math.abs(l.p - c.p) <= 1e-9 * c.p, at
				assert l.phases == c.phases.len, at
			}
		}
	}
	// Measured: 2e-12 in the split; 1.9e-11 in V, 8.1e-11 in H and 9.8e-11
	// in S, each relative to |X| plus RT (H, U) or R (S).
	println('split ${d_split.worst:.2g} (${d_split.at}); V ${d_v.worst:.2g} (${d_v.at}); H ${d_h.worst:.2g} (${d_h.at}); S ${d_s.worst:.2g} (${d_s.at})')
	assert d_split.worst <= 1e-8, d_split.at
	assert d_v.worst <= 1e-9, d_v.at
	assert d_h.worst <= 1e-9, d_h.at
	assert d_s.worst <= 1e-9, d_s.at
	assert splits == 39
}

fn test_saturation_matches_thermo() {
	// A pure component's saturation states, which no PT flash reaches: the
	// two roots have equal fugacity in reflux's model too (exit criterion
	// 3's 1e-10; measured 5.6e-13), saturation_pressure finds the same P to
	// 1e-9 (it stops at a relative change of 1e-10, measured 3.8e-10 for CO₂
	// at Tr = 0.99), and V, H and S at the vapour fraction match to 1e-9.
	f := read_thermo()!
	mut d_res, mut d_psat, mut d_v, mut d_h, mut d_s := Dev{}, Dev{}, Dev{}, Dev{}, Dev{}
	for c in f.saturation {
		at := '${c.components} ${c.spec} given, vf ${c.vf}: ${c.t} K, ${c.p} Pa'
		m := cubic_of(f, 'pr', c.components, c.kij)!
		fl := fluid_of(f, m, c.components)!
		liq := thermo.solve_tp(m, c.t, c.p, [1.0], .liquid)!
		vap := thermo.solve_tp(m, c.t, c.p, [1.0], .vapor)!
		assert liq.v < vap.v, '${at}: two roots'
		d_res.see(math.abs(liq.ln_phi[0] - vap.ln_phi[0]), at)
		psat := thermo.saturation_pressure(m, c.t)!
		d_psat.see(math.abs(psat - c.p) / c.p, at)
		l, g := molar_of(fl, f, liq), molar_of(fl, f, vap)
		d_v.see(math.abs((1.0 - c.vf) * l.v + c.vf * g.v - c.v) / c.v, at)
		d_h.see(math.abs((1.0 - c.vf) * l.h + c.vf * g.h - c.h) / (math.abs(c.h) + thermo.r_gas * c.t), at)
		d_s.see(math.abs((1.0 - c.vf) * l.s + c.vf * g.s - c.s) / (math.abs(c.s) + thermo.r_gas), at)
	}
	println('ln φL − ln φV ${d_res.worst:.2g} (${d_res.at}); Psat ${d_psat.worst:.2g} (${d_psat.at}); V ${d_v.worst:.2g}; H ${d_h.worst:.2g}; S ${d_s.worst:.2g}')
	assert d_res.worst <= 1e-10, d_res.at
	assert d_psat.worst <= 1e-9, d_psat.at
	assert d_v.worst <= 1e-9 && d_h.worst <= 1e-9 && d_s.worst <= 1e-9
}

// incipient_fraction is the share of the feed in the phase nearer the
// incipient composition w: 0 for a single phase.
fn incipient_fraction(r flash.FlashResult, w []f64) f64 {
	if r.phases.len == 1 {
		return 0.0
	}
	mut best, mut frac := math.inf(1), 0.0
	for ph in r.phases {
		mut d := 0.0
		for i, x in ph.state.n {
			d += math.abs(x - w[i])
		}
		if d < best {
			best, frac = d, ph.fraction
		}
	}
	return frac
}

fn test_bubble_and_dew_points_are_boundaries() {
	// thermo's bubble and dew points, polished on its own ln φ, satisfy
	// reflux's equilibrium condition (exit criterion 3's 1e-10; measured
	// 1.3e-13), the feed's V, H and S match to 1e-9, and reflux's PT flash
	// sees a phase boundary there: 1e-7 to one side a single phase, to the
	// other a split with a trace of the incipient phase (at most 2.5e-4).
	// The two thermo refused are the condensate's bubble points at a given
	// T, where FlashVL raises an UnboundLocalError.
	f := read_thermo()!
	mut d_res, mut d_v, mut d_h, mut d_s := Dev{}, Dev{}, Dev{}, Dev{}
	mut refused, mut most := 0, 0.0
	for c in f.boundaries {
		at := '${c.eos} ${c.components} z ${c.z}, ${c.spec} given, vf ${c.vf}: ${c.t} K, ${c.p} Pa'
		if c.error != '' || c.rejected != '' {
			refused++
			continue
		}
		m := cubic_of(f, c.eos, c.components, c.kij)!
		fl := fluid_of(f, m, c.components)!
		feed_root, new_root := if c.vf == 0.0 {
			thermo.Root.liquid, thermo.Root.vapor
		} else {
			thermo.Root.vapor, thermo.Root.liquid
		}
		feed := thermo.solve_tp(m, c.t, c.p, c.z, feed_root)!
		new := thermo.solve_tp(m, c.t, c.p, c.incipient, new_root)!
		mut sum := 0.0
		for i, w in c.incipient {
			sum += w
			d_res.see(math.abs(math.log(c.z[i]) + feed.ln_phi[i] - math.log(w) - new.ln_phi[i]), at)
		}
		assert math.abs(sum - 1.0) <= 1e-14, at
		mol := molar_of(fl, f, feed)
		d_v.see(math.abs(mol.v - c.v) / c.v, at)
		d_h.see(math.abs(mol.h - c.h) / (math.abs(c.h) + thermo.r_gas * c.t), at)
		d_s.see(math.abs(mol.s - c.s) / (math.abs(c.s) + thermo.r_gas), at)
		mut sides := []f64{}
		for d in [1e-7, -1e-7] {
			t, p := if c.spec == 'p' { c.t * (1.0 + d), c.p } else { c.t, c.p * (1.0 + d) }
			r := flash.pt(m, t, p, c.z)!
			sides << incipient_fraction(r, c.incipient)
		}
		sides.sort()
		assert sides[0] == 0.0 && sides[1] > 0.0 && sides[1] < 1e-3, '${at}: incipient fractions either side ${sides}'
		most = math.max(most, sides[1])
	}
	println('ln f ${d_res.worst:.2g} (${d_res.at}); V ${d_v.worst:.2g}; H ${d_h.worst:.2g}; S ${d_s.worst:.2g}; largest trace ${most:.2g}')
	assert d_res.worst <= 1e-10, d_res.at
	assert d_v.worst <= 1e-9 && d_h.worst <= 1e-9 && d_s.worst <= 1e-9
	assert refused == 2
}
