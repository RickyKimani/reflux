// Standard reaction properties and K(T) against mpmath on ChemSep's own data
// (M9 step 4, D7), with thermo's data and the owner's DWSIM correlations
// beside them, from testdata/m9_reactions.json.
module kinetics

import components
import expr
import json2
import math
import os
import thermo

struct FixtureSpecies {
	cas string
	nu  f64
}

struct FixturePoint {
	t           f64
	dh          f64
	ds          f64
	dg          f64
	ln_k        f64
	dln_k_dt    f64
	thermo_ln_k f64
}

struct FixtureReaction {
	name     string
	equation string
	species  []FixtureSpecies
	points   []FixturePoint
	dgf_sum  f64
	gap      f64
}

struct FixtureCorrelation {
	reaction     string
	a            f64
	b            f64
	dn           f64
	temperatures []f64
	diffs        []f64
}

struct Fixture {
	reactions    []FixtureReaction
	correlations []FixtureCorrelation
}

const db = components.load_default() or { panic(err) }

fn load() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_reactions.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

// species_of is a compound as reactions know it, from its component file.
fn species_of(cas string) !Species {
	r := db.find(cas)!
	return Species{
		name:    r.name
		cas:     r.ids.cas
		aliases: r.ids.aliases
		formula: r.formula()!
		charge:  r.constants.charge or { 0 }
		gas:     r.ideal_gas()!
		dgf:     r.constants.dgf
	}
}

// fixture_reaction is a fixture's reaction read from its equation, over its
// own species.
fn fixture_reaction(f FixtureReaction) !(Reaction, []Species) {
	mut species := []Species{}
	for s in f.species {
		species << species_of(s.cas)!
	}
	r := parse_equation(f.name, f.equation, species)!
	return r, species
}

// magnitude is what a standard property's rounding is measured against: each
// term's magnitude summed, H°ᵢ's parts and T·S°ᵢ's, in J/mol.
fn magnitude(r Reaction, species []Species, t f64) f64 {
	mut s := 0.0
	for j, k in r.species {
		g := species[k].gas or { panic('no gas') }
		h := math.abs(g.dhf or { 0 }) + math.abs(g.cp.int_cp(thermo.t_ref, t))
		e := math.abs(g.s0 or { 0 }) + math.abs(g.cp.int_cp_over_t(thermo.t_ref, t))
		s += math.abs(r.nu[j]) * (h + t * e)
	}
	return s
}

// Every reaction reads from its equation to the fixture's species and
// coefficients, and its ΔH°, ΔS°, ΔG°, ln K and d ln K/dT are mpmath's to
// rounding. Each species' H° and S° add numbers of up to 1e6 J/mol, and the
// Cp integrals are good to 1e-13 of themselves (tdep's QAG tolerance), so
// each property is held to 1e-13 of the sum of its terms' magnitudes
// (largest measured: 1.3e-16).
fn test_standard_properties_are_mpmaths() {
	f := load()!
	assert f.reactions.len == 12
	mut worst := 0.0
	for fr in f.reactions {
		r, species := fixture_reaction(fr)!
		assert r.nu == fr.species.map(it.nu), fr.name
		assert r.species == []int{len: species.len, init: index}, fr.name
		for p in fr.points {
			st := r.standard(species, p.t)!
			sc := magnitude(r, species, p.t)
			rt := gas_constant * p.t
			checks := [
				[st.dh, p.dh, sc],
				[st.ds, p.ds, sc / p.t],
				[st.dg, p.dg, sc],
				[st.ln_k, p.ln_k, sc / rt],
				[st.dln_k_dt, p.dln_k_dt, sc / (rt * p.t)],
			]
			for i, c in checks {
				e := math.abs(c[0] - c[1]) / c[2]
				worst = math.max(worst, e)
				assert e <= 1e-13, '${fr.name} at ${p.t} K, property ${i}: ${c[0]} vs ${c[1]} (${e} of the scale)'
			}
		}
	}
	println('standard properties: worst ${worst} of the scale')
}

// Σν ΔG°f agrees with ΔG° at 298.15 K from ΔH°f and S° to within each
// file's rounding of its ΔG°f, ΔH°f and S°: the gaps are mpmath's to 1e-9
// J/mol, and every one is under 0.2 kJ/mol, which moves ln K by at most
// 0.07 at 298.15 K.
fn test_gibbs_energies_of_formation_agree() {
	f := load()!
	for fr in f.reactions {
		r, species := fixture_reaction(fr)!
		gap := r.formation_gap(species)!
		assert math.abs(gap - fr.gap) <= 1e-9 * (1 + math.abs(fr.dgf_sum)), '${fr.name}: ${gap} vs ${fr.gap}'
		assert math.abs(gap) < 200, '${fr.name}: ${gap}'
	}
}

// ln K from ChemSep's data against thermo's own data (chemicals' ΔH°f and
// S° at 1 bar, thermo's Cp), the largest distance over 298.15 to 1500 K for
// each reaction. These are the data's differences, recorded so that they
// can't grow: ChemSep's S° is at 1 atm for some compounds (H₂, H₂O, CO,
// CO₂, N₂, NH₃: 0.10 to 0.11 J/(mol·K) below thermo's) and at 1 bar for
// others (CH₄, CH₃OH), which moves reforming's ln K by about −0.04; thermo
// estimates one of ethylbenzene dehydrogenation's species by Joback's
// groups; and the butanes' formation data differ by 1.9 kJ/mol.
const thermo_distance = {
	'steam reforming':              0.065
	'water-gas shift':              0.016
	'overall reforming':            0.050
	'methanol from CO2':            0.129
	'methanol from CO':             0.141
	'ammonia':                      0.300
	'methane combustion':           0.083
	'ethylene hydration':           0.115
	'SO2 oxidation':                0.102
	'butane isomerisation':         0.769
	'propane dehydrogenation':      0.108
	'ethylbenzene dehydrogenation': 0.478
}

fn test_k_against_thermos_data() {
	f := load()!
	for fr in f.reactions {
		r, species := fixture_reaction(fr)!
		mut far := 0.0
		for p in fr.points {
			st := r.standard(species, p.t)!
			far = math.max(far, math.abs(st.ln_k - p.thermo_ln_k))
		}
		println("${fr.name}: ln K within ${far} of thermo's")
		assert far <= thermo_distance[fr.name], fr.name
	}
}

// The equilibrium constants the owner's DWSIM case types into its rates,
// against reflux's from ChemSep's data, as the fixture records them: steam
// reforming's and the shift's within 0.1 of ln K from 700 to 1300 K, but the
// overall reaction's isn't steam reforming's times the shift's, and is off
// by 4.45 at 700 K, a factor of 86, and by −1.21 at 1300 K.
fn test_the_owners_correlations() {
	f := load()!
	mut by_name := map[string]FixtureReaction{}
	for fr in f.reactions {
		by_name[fr.name] = fr
	}
	mut far := map[string]f64{}
	for c in f.correlations {
		r, species := fixture_reaction(by_name[c.reaction])!
		for i, t in c.temperatures {
			st := r.standard(species, t)!
			file := c.a + c.b / t + c.dn * math.log(1.01325)
			d := st.ln_k - file
			assert math.abs(d - c.diffs[i]) <= 1e-9, '${c.reaction} at ${t} K: ${d} vs ${c.diffs[i]}'
			far[c.reaction] = math.max(far[c.reaction], math.abs(d))
		}
	}
	assert far['steam reforming'] < 0.07
	assert far['water-gas shift'] < 0.1
	assert far['overall reforming'] > 4.4
	sr, wgs, overall := f.correlations[0], f.correlations[1], f.correlations[2]
	assert overall.a != sr.a + wgs.a && overall.b != sr.b + wgs.b
}

// K is a reaction's: a combination's ln K is its parts' combined, to
// rounding, as the overall reforming reaction's is steam reforming's plus
// the shift's.
fn test_k_combines_as_its_reactions_do() {
	f := load()!
	sr, species := fixture_reaction(f.reactions[0])!
	wgs, s2 := fixture_reaction(f.reactions[1])!
	overall, s3 := fixture_reaction(f.reactions[2])!
	for t in [298.15, 700.0, 1000.0, 1500.0] {
		a := sr.standard(species, t)!
		b := wgs.standard(s2, t)!
		c := overall.standard(s3, t)!
		assert math.abs(a.ln_k + b.ln_k - c.ln_k) <= 1e-13 * magnitude(overall, s3, t) / (gas_constant * t)
		assert math.abs(a.dh + b.dh - c.dh) <= 1e-13 * magnitude(overall, s3, t)
	}
}

fn test_missing_data_is_refused() {
	mut species := [species_of('74-82-8')!, species_of('7732-18-5')!, species_of('630-08-0')!,
		species_of('1333-74-0')!]
	species[0] = Species{
		...species[0]
		gas: none
	}
	r := parse_equation('sr', 'CH4 + H2O = CO + 3 H2', species)!
	if _ := r.standard(species, 1000) {
		assert false
	} else {
		assert err.msg() == 'sr: CH4 has no ideal-gas data'
	}
	ok := [species_of('74-82-8')!, species[1], species[2], species[3]]
	for t in [0.0, -1, math.nan(), math.inf(1)] {
		if _ := r.standard(ok, t) {
			assert false, '${t}'
		}
	}
	mut nodgf := ok.clone()
	nodgf[2] = Species{
		...ok[2]
		dgf: none
	}
	if _ := r.formation_gap(nodgf) {
		assert false
	} else {
		assert err.msg() == 'sr: CO has no ΔG°f'
	}
}

// K uses the gas constant rates do, expr's built-in R.
fn test_the_gas_constant_is_exprs() {
	c := expr.constants.filter(it.name == 'R')
	assert c.len == 1 && c[0].value == gas_constant
}
