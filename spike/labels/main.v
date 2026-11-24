module main

// M6 step 7's measurement for D9: how M3's sweep's two-phase splits would be
// labelled by Π (with the four-covolume rule) rather than by density alone,
// as today, where the less dense phase is the vapour. Each split is one of:
// - same: the lighter phase vapour-like, the denser liquid-like;
// - both liquid-like: by D9 a liquid–liquid split, liquid and liquid2;
// - both vapour-like: by D9 the lighter is the vapour and the denser a liquid,
//   the labels as today;
// - inverted: the lighter liquid-like and the denser vapour-like.
//
//     v -prod -cc gcc run spike/labels
import components
import flash
import math
import sweep
import thermo

fn density(fl thermo.Fluid, st thermo.State) f64 {
	comps := fl.model.components()
	mut mass := 0.0
	for i, n in st.n {
		mass += n * comps[i].mw
	}
	return mass / st.v
}

fn main() {
	golden_labels()!
	db := components.load_default()!
	pool := sweep.candidates(db)!
	pts := sweep.points(pool, sweep.full_count)
	mut counts := map[string]int{}
	mut gaps := map[string][]f64{}
	mut examples := map[string][]string{}
	for pt in pts {
		fl := db.fluid(.pr, pt.cas) or { continue }
		r := flash_pt(fl, pt) or { continue }
		if r.phases.len != 2 {
			continue
		}
		eos := fl.helmholtz()!
		light, heavy := r.phases[0].state, r.phases[1].state
		ll, lh := thermo.phase_label(eos, light), thermo.phase_label(eos, heavy)
		kind := match true {
			ll == .vapor && lh == .liquid { 'same' }
			ll == .liquid && lh == .liquid { 'both liquid-like' }
			ll == .vapor && lh == .vapor { 'both vapour-like' }
			else { 'inverted' }
		}
		counts[kind]++
		dl, dh := density(fl, light), density(fl, heavy)
		gap := (dh - dl) / dh
		gaps[kind] << gap
		if kind != 'same' && examples[kind].len < 6 {
			examples[kind] << '${pt.index}: ${pt.names} ${math.round(pt.t * 10) / 10} K ${math.round(pt.p / 1e3) / 100} bar, ${math.round(dl)} and ${math.round(dh)} kg/m³, Π ${math.round(thermo.pip(eos, light.t, light.v, light.n) * 1000) / 1000} and ${math.round(thermo.pip(eos, heavy.t, heavy.v, heavy.n) * 1000) / 1000}'
		}
	}
	for k, n in counts {
		mut g := gaps[k].clone()
		g.sort()
		median := g[g.len / 2]
		println('${k}: ${n} splits; density gap (heavy − light)/heavy from ${math.round(g[0] * 1e4) / 1e4} to ${math.round(g[g.len - 1] * 1e4) / 1e4}, median ${math.round(median * 1e4) / 1e4}')
		for e in examples[k] {
			println('    ${e}')
		}
	}
}

fn flash_pt(fl thermo.Fluid, pt sweep.Point) !flash.FlashResult {
	return flash.pt(fl.model, pt.t, pt.p, pt.z)
}
