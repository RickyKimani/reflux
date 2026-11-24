module components

import activity
import tdep
import thermo

// A γ–φ fluid from the database (M6's D4 to D6): the activity model from
// ChemSep's banks or UNIFAC tables, each solvent's vapour pressure and
// liquid density from its component file, each solute's Henry's-law
// constants from henrys_law.ipd, and the vapour as the ideal gas or a
// cubic.

// LiquidModel names an activity model.
pub enum LiquidModel {
	nrtl
	uniquac
	wilson
	unifac
	unifac_lle
	unifac_dortmund
}

// GammaPhiOptions say how a γ–φ fluid is built.
@[params]
pub struct GammaPhiOptions {
pub:
	liquid LiquidModel = .nrtl
	// vapour is a cubic for the vapour, or none for the ideal gas (D4's
	// default).
	vapour ?thermo.CubicKind
	// kij takes the cubic's kᵢⱼ from its bank, as a cubic fluid does.
	kij      bool = true
	poynting bool = true
	// henry names the Henry's-law solutes, by name or CAS number (D6). A
	// compound with constants in a solvent present whose Tc is below every
	// solvent's vapour pressure's range is one too, whether named or not:
	// wherever the liquid can be had, it's above its critical point. In
	// water, a gas IAPWS covers takes IAPWS's constant, the bank's elsewhere.
	henry []string
	// ideal lets a pair the bank doesn't list have no energy of
	// interaction (D7 as amended).
	ideal bool
	// pairs are a binary model's pairs given by the caller, a file's
	// [[model.pairs]] (D12), in place of the bank's: each of two solvents,
	// by CAS number.
	pairs []PairOverride
}

// GammaPhiSetup is a γ–φ fluid and how it was built.
pub struct GammaPhiSetup {
pub:
	fluid   thermo.Fluid
	model   thermo.GammaPhi
	records []Record
	solutes []string  // by CAS number
	pairs   []PairUse // each solvent pair's source, for the binary models
	notes   []string  // a vapour pressure from its Antoine record; a gas's constant from IAPWS
}

// gamma_phi builds a γ–φ fluid of the compounds named.
pub fn (db Database) gamma_phi(names []string, opts GammaPhiOptions) !GammaPhiSetup {
	if names.len == 0 {
		return error('a fluid needs at least one compound')
	}
	mut records := []Record{}
	mut comps := []thermo.Component{}
	mut ideal := []thermo.IdealGas{}
	for name in names {
		r := db.find(name)!
		comps << r.component()!
		ideal << r.ideal_gas()!
		records << r
	}
	cas := records.map(it.ids.cas)
	nc := cas.len
	mut named := []int{}
	for h in opts.henry {
		r := db.find(h)!
		i := cas.index(r.ids.cas)
		if i < 0 {
			return error("${r.name}, named as a Henry's-law solute, isn't one of the compounds")
		}
		named << i
	}
	bank := load_henry(db.vendor)!
	// Each solvent's vapour pressure and density, for every compound not
	// named a solute; then which of the rest are solutes by D6's rule.
	mut psat := []thermo.Scaled{len: nc}
	mut density := []thermo.Scaled{len: nc}
	mut notes := []string{}
	mut lowest := 1e300
	mut usable := []bool{len: nc}
	for i, r in records {
		if i in named {
			continue
		}
		vp := r.property('vapor_pressure') or { continue }
		rho := r.property('liquid_density') or { continue }
		if vp.c.eq !in tdep.jet_forms || rho.c.eq !in tdep.jet_forms {
			continue
		}
		psat[i] = thermo.Scaled{vp.c, vp.scale}
		density[i] = thermo.Scaled{rho.c, rho.scale}
		usable[i] = true
		if vp.note != '' {
			notes << '${r.name}: ${vp.note}'
		}
	}
	// A compound with constants in another one present can be a solute;
	// one without can only be a solvent, and the lowest T any of those has a
	// vapour pressure at is where the liquid can first be had.
	mut solute := []bool{len: nc}
	for i in named {
		solute[i] = true
	}
	mut candidate := []bool{len: nc}
	for i, r in records {
		candidate[i] = bank.entries.any(it.gas == r.ids.cas && it.solvent in cas
			&& it.solvent != r.ids.cas) || (r.ids.cas in iapws_henry && water_cas in cas)
	}
	for i in 0 .. nc {
		if !solute[i] && !candidate[i] && usable[i] && psat[i].c.tmin < lowest {
			lowest = psat[i].c.tmin
		}
	}
	for i in 0 .. nc {
		if candidate[i] && comps[i].tc < lowest {
			solute[i] = true
		}
	}
	mut solvents := []int{}
	mut solutes := []thermo.Solute{}
	for i, r in records {
		if !solute[i] {
			if !usable[i] {
				return error("${r.name} has no vapour pressure and liquid density reflux can use for a γ–φ liquid (forms ${tdep.jet_forms}), and isn't a Henry's-law solute")
			}
			solvents << i
			continue
		}
		mut henrys := []thermo.Henry{}
		// In water, IAPWS's constant where it has one, not the bank's row
		// (D6 as amended after step 6).
		w := cas.index(water_cas)
		mut by_iapws := false
		if w >= 0 && !solute[w] {
			if h := iapws_henry_in_water(r.ids.cas, w) {
				henrys << h
				by_iapws = true
			}
		}
		if by_iapws {
			replaced := bank.entries.any(it.gas == r.ids.cas && it.solvent == water_cas)
			notes << "${r.name} in ${records[w].name}: IAPWS's Henry's constant (G7-04)" +
				if replaced { ", in place of ${bank.file}'s" } else { '' }
		}
		for e in bank.entries {
			j := cas.index(e.solvent)
			if e.gas != r.ids.cas || j < 0 || solute[j] || (by_iapws && j == w) {
				continue
			}
			if e.eq != 121 || e.coeffs.len != 5 {
				return error('${bank.file} line ${e.line}: ${r.name} in ${records[j].name} is equation ${e.eq}, not 121')
			}
			k := e.coeffs
			henrys << thermo.Henry{
				solvent: j
				form:    .chemsep_121
				a:       k[0]
				b:       k[1]
				c:       k[2]
				d:       k[3]
				e:       k[4]
				tmin:    e.tmin
				tmax:    e.tmax
			}
		}
		if henrys.len == 0 {
			return error("${r.name} is a Henry's-law solute, and ${bank.file} has no constant for it in any of the solvents")
		}
		solutes << thermo.Solute{i, henrys}
	}
	if solvents.len == 0 {
		return error("a γ–φ liquid needs a solvent, and every compound is a Henry's-law solute")
	}
	// The liquid can be had only where every solvent's correlations hold
	// (D4). Where no T does, one of them is above its critical point
	// wherever the others are liquid, with no Henry's-law constant to be a
	// solute by (criterion 7): refused here, by name, not at every phase.
	mut ends, mut begins := solvents[0], solvents[0]
	for i in solvents {
		if liquid_top(psat[i], density[i]) < liquid_top(psat[ends], density[ends]) {
			ends = i
		}
		if liquid_bottom(psat[i], density[i]) > liquid_bottom(psat[begins], density[begins]) {
			begins = i
		}
	}
	hi := liquid_top(psat[ends], density[ends])
	lo := liquid_bottom(psat[begins], density[begins])
	if hi < lo {
		a, b := records[ends].name, records[begins].name
		return error("no temperature has a liquid of ${a} and ${b}: ${a}'s vapour pressure and liquid density end at ${hi} K (its Tc is ${comps[ends].tc} K), and ${b}'s begin at ${lo} K; ${bank.file} has no constant for ${a} in the others, so it can't be a Henry's-law solute")
	}
	solvent_cas := solvents.map(cas[it])
	for o in opts.pairs {
		for c in [o.a, o.b] {
			i := cas.index(c)
			if i < 0 {
				return error("a pair names ${c}, which isn't one of the compounds")
			}
			if solute[i] {
				return error("a pair names ${records[i].name}, a Henry's-law solute, which is outside the activity model: a pair is of two solvents")
			}
		}
	}
	if opts.pairs.len > 0 && opts.liquid in [.unifac, .unifac_lle, .unifac_dortmund] {
		return error("UNIFAC's parameters are its groups', from its table: pairs are NRTL's, UNIQUAC's and Wilson's")
	}
	mut pairs := []PairUse{}
	liquid := match opts.liquid {
		.nrtl {
			b := db.nrtl(solvent_cas, ideal: opts.ideal, overrides: opts.pairs)!
			pairs = b.pairs.clone()
			activity.Model(b.model)
		}
		.uniquac {
			b := db.uniquac(solvent_cas, ideal: opts.ideal, overrides: opts.pairs)!
			pairs = b.pairs.clone()
			activity.Model(b.model)
		}
		.wilson {
			b := db.wilson(solvent_cas, ideal: opts.ideal, overrides: opts.pairs)!
			pairs = b.pairs.clone()
			activity.Model(b.model)
		}
		.unifac {
			activity.Model(db.unifac(.original, solvent_cas)!)
		}
		.unifac_lle {
			activity.Model(db.unifac(.lle, solvent_cas)!)
		}
		.unifac_dortmund {
			activity.Model(db.unifac(.dortmund, solvent_cas)!)
		}
	}
	mut vapour := ?thermo.Model(none)
	mut pures := []thermo.Model{}
	if kind := opts.vapour {
		eos := db.fluid(kind, names, kij: opts.kij)!.helmholtz()!
		if eos is thermo.Cubic {
			for i in 0 .. nc {
				pures << thermo.new_cubic(kind, [eos.comps[i]], alphas: [eos.alphas[i]])!
			}
		}
		vapour = eos
	}
	model := thermo.new_gamma_phi(thermo.GammaPhi{
		comps:    comps
		liquid:   liquid
		solvents: solvents
		psat:     psat
		density:  density
		vapour:   vapour
		pures:    pures
		poynting: opts.poynting
		solutes:  solutes
	})!
	return GammaPhiSetup{
		fluid:   thermo.new_phase_fluid(model, ideal)!
		model:   model
		records: records
		solutes: solutes.map(cas[it.i])
		pairs:   pairs
		notes:   notes
	}
}

// liquid_top is where a solvent's vapour pressure or liquid density first ends,
// and liquid_bottom where both have begun: its range as a liquid.
fn liquid_top(a thermo.Scaled, b thermo.Scaled) f64 {
	return if a.c.tmax < b.c.tmax { a.c.tmax } else { b.c.tmax }
}

fn liquid_bottom(a thermo.Scaled, b thermo.Scaled) f64 {
	return if a.c.tmin > b.c.tmin { a.c.tmin } else { b.c.tmin }
}

const water_cas = '7732-18-5'

// IapwsHenry is one gas's Henry's constant in water by IAPWS (G7-04).
struct IapwsHenry {
	a    f64
	b    f64
	c    f64
	tmin f64 // K, the range of the data it was fitted to
	tmax f64
}

// iapws_henry is IAPWS's Guideline on the Henry's Constant and Vapor–Liquid
// Distribution Constant for Gases in H2O and D2O at High Temperatures
// (G7-04, 2004; Fernández-Prini, Alvarez and Harvey, J. Phys. Chem. Ref.
// Data 32, 903, 2003), Table 2, for H2O, by each gas's CAS number. A gas in
// water takes it in place of henrys_law.ipd's row: there hydrogen's is 5.95
// times IAPWS's at 298 K and neon's 4.5e32 times, and the other eleven's
// drift to 10 to 55% off at their hottest (data-quality.md).
const iapws_henry = {
	'7440-59-7': IapwsHenry{-3.52839, 7.12983, 4.47770, 273.21, 553.18} // helium
	'7440-01-9': IapwsHenry{-3.18301, 5.31448, 5.43774, 273.20, 543.36} // neon
	'7440-37-1': IapwsHenry{-8.40954, 4.29587, 10.52779, 273.19, 568.36} // argon
	'7439-90-9': IapwsHenry{-8.97358, 3.61508, 11.29963, 273.19, 525.56} // krypton
	'7440-63-3': IapwsHenry{-14.21635, 4.00041, 15.60999, 273.22, 574.85} // xenon
	'1333-74-0': IapwsHenry{-4.73284, 6.08954, 6.06066, 273.15, 636.09} // hydrogen
	'7727-37-9': IapwsHenry{-9.67578, 4.72162, 11.70585, 278.12, 636.46} // nitrogen
	'7782-44-7': IapwsHenry{-9.44833, 4.43822, 11.42005, 274.15, 616.52} // oxygen
	'630-08-0':  IapwsHenry{-10.52862, 5.13259, 12.01421, 278.15, 588.67} // carbon monoxide
	'124-38-9':  IapwsHenry{-8.55445, 4.01195, 9.52345, 274.19, 642.66} // carbon dioxide
	'7783-06-4': IapwsHenry{-4.51499, 5.23538, 4.42126, 273.15, 533.09} // hydrogen sulfide
	'74-82-8':   IapwsHenry{-10.44708, 4.66491, 12.12986, 275.46, 633.11} // methane
	'74-84-0':   IapwsHenry{-19.67563, 4.51222, 20.62567, 275.44, 473.46} // ethane
	'2551-62-4': IapwsHenry{-16.56118, 2.15289, 20.35440, 283.14, 505.55} // sulfur hexafluoride
}

// iapws_henry_in_water is IAPWS's Henry's constant for a gas in water, by
// its CAS number, with water at index solvent, or none if IAPWS hasn't it:
// for the quality report and the tests.
pub fn iapws_henry_in_water(cas string, solvent int) ?thermo.Henry {
	c := iapws_henry[cas] or { return none }
	return thermo.Henry{
		solvent: solvent
		form:    .iapws
		a:       c.a
		b:       c.b
		c:       c.c
		tmin:    c.tmin
		tmax:    c.tmax
	}
}
