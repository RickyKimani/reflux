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
	// wherever the liquid can be had, it's above its critical point.
	henry []string
	// ideal lets a pair the bank doesn't list have no energy of
	// interaction (D7 as amended).
	ideal bool
}

// GammaPhiSetup is a γ–φ fluid and how it was built.
pub struct GammaPhiSetup {
pub:
	fluid   thermo.Fluid
	model   thermo.GammaPhi
	records []Record
	solutes []string  // by CAS number
	pairs   []PairUse // each solvent pair's source, for the binary models
	notes   []string  // a vapour pressure taken from its Antoine record, and why
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
			&& it.solvent != r.ids.cas)
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
		for e in bank.entries {
			j := cas.index(e.solvent)
			if e.gas != r.ids.cas || j < 0 || solute[j] {
				continue
			}
			if e.eq != 121 || e.coeffs.len != 5 {
				return error('${bank.file} line ${e.line}: ${r.name} in ${records[j].name} is equation ${e.eq}, not 121')
			}
			k := e.coeffs
			henrys << thermo.Henry{j, k[0], k[1], k[2], k[3], k[4], e.tmin, e.tmax}
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
	mut pairs := []PairUse{}
	liquid := match opts.liquid {
		.nrtl {
			b := db.nrtl(solvent_cas, ideal: opts.ideal)!
			pairs = b.pairs.clone()
			activity.Model(b.model)
		}
		.uniquac {
			b := db.uniquac(solvent_cas, ideal: opts.ideal)!
			pairs = b.pairs.clone()
			activity.Model(b.model)
		}
		.wilson {
			b := db.wilson(solvent_cas, ideal: opts.ideal)!
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
