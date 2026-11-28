// DWSIM 9.0.5's transport of mixtures (M5b step 10, criterion 9;
// testdata/m5b_dwsim.json from tools/gen_m5b_dwsim.py): its PR and NRTL
// packages by state, every number it reports worked out again from its own
// formulas, and every difference from reflux's put down to the data or to
// the method.
//
// - Its formulas on its own compound values give what it reports, to
//   rounding: every number within 1e-15 (4.7e-16 measured), each a sum of a
//   few terms.
// - Its formulas on reflux's component files' values ("on files") are its
//   numbers to 3e-6: its ChemSep correlations are the files', so the data
//   explain nothing of the difference (benzene–toluene's vapour
//   conductivity is 3e-6 off, its vapour-conductivity function's own).
// - So the difference between reflux's number and DWSIM's is the methods',
//   D4 to D9's against DWSIM's: a liquid's viscosity by the logarithmic
//   rule, not the mole-fraction mean, its Lucas factor 1 below Psat, not
//   at |P − Psat| (D6); a gas's by Wilke's rule times Lucas's factor, not
//   the mole-fraction mean plus Jossi–Stiel–Thodos at every density (D4,
//   D5); a liquid's conductivity by Li's rule on the liquids' volumes with
//   Missenard's factor, not on critical volumes (D7); a gas's by
//   Mason–Saxena plus Stiel–Thodos, not the mole-fraction mean (D4, D5);
//   the surface tension by Winterfeld, Scriven and Davis, not the
//   mole-fraction mean (D8); a compound above its Tc left out of a liquid
//   (D9), where DWSIM keeps its liquid viscosity extrapolated with its
//   Henry's constant (in water) for Psat, and drops it from the surface
//   tension in an order-dependent way. Each is recorded (`method`), a
//   number that can't grow, at DWSIM's own phase compositions and
//   densities.
// - A second liquid: DWSIM gives it the first's viscosity and conductivity
//   (its functions called without the phase), so there the difference is
//   that too.
import components
import json2
import math
import os
import thermo
import water

// The fixture's types, named apart from every other module's
// (docs/v-for-go.md).
struct DwPure {
	name  string
	cas   string
	tc    f64
	sigma ?f64
	henry ?f64
}

struct DwValues {
	mu    ?f64
	k     ?f64
	sigma ?f64
}

struct DwPhase {
	phase      string
	fraction   f64
	x          []f64
	density    f64
	molar_mass f64
	dwsim      DwValues
	reproduced DwValues
	on_files   DwValues
}

struct DwCase {
	label      string
	package    string
	compounds  []string
	cas        []string
	z          []f64
	t          f64 @[json: 'T']
	p          f64 @[json: 'P']
	pure_dwsim []DwPure
	phases     []DwPhase
}

struct DwWaterRow {
	t       f64 @[json: 'T']
	p       f64 @[json: 'P']
	phase   string
	density f64
	mu      f64
	k       f64
}

struct DwWater {
	steam_tables []DwWaterRow
	coolprop     []DwWaterRow
}

struct DwFile {
	dwsim string
	cases []DwCase
	water DwWater
}

const db = components.load_default() or { panic(err) }

fn dw() !DwFile {
	return json2.decode[DwFile](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5b_dwsim.json'))!)!
}

// The methods' differences, reflux's number over DWSIM's formulas' on the
// files' values, less 1, by case, phase and property: each recorded to four
// figures and held to 0.1% of itself, so that a change to either side's
// method, toward the other's or away, fails here until it's explained.
const method = {
	'benzene-toluene liquid, 1 atm liquid1 mu':      -0.0007246
	'benzene-toluene liquid, 1 atm liquid1 k':       0.001588
	'benzene-toluene liquid, 1 atm liquid1 sigma':   -0.0003658
	'benzene-toluene liquid, 100 bar liquid1 mu':    -0.0006176
	'benzene-toluene liquid, 100 bar liquid1 k':     0.02566
	'benzene-toluene liquid, 100 bar liquid1 sigma': -0.0003658
	'benzene-toluene, two phases vapour mu':         -0.02747
	'benzene-toluene, two phases vapour k':          0.01135
	'benzene-toluene, two phases liquid1 mu':        -0.0003346
	'benzene-toluene, two phases liquid1 k':         0.001637
	'benzene-toluene, two phases liquid1 sigma':     0.00115
	'methanol-water, two phases vapour mu':          -0.02964
	'methanol-water, two phases vapour k':           0.008133
	'methanol-water, two phases liquid1 mu':         -0.003392
	'methanol-water, two phases liquid1 k':          -0.03458
	'methanol-water, two phases liquid1 sigma':      -0.2687
	'methanol-water liquid liquid1 mu':              -0.02262
	'methanol-water liquid liquid1 k':               -0.01807
	'methanol-water liquid liquid1 sigma':           -0.2296
	'reforming gas, 30 bar vapour mu':               0.2273
	'reforming gas, 30 bar vapour k':                -0.1483
	'methane-propane, 100 bar vapour mu':            -0.02821
	'methane-propane, 100 bar vapour k':             0.3233
	'methane in hexane, 20 bar liquid1 mu':          0.04976
	'methane in hexane, 20 bar liquid1 k':           0.0289
	'methane in hexane, 20 bar liquid1 sigma':       0.0
	'hexane with methane, 20 bar liquid1 mu':        0.04976
	'hexane with methane, 20 bar liquid1 k':         0.0289
	'hexane with methane, 20 bar liquid1 sigma':     0.05263
	'water-1-butanol, two liquids liquid1 mu':       -0.1223
	'water-1-butanol, two liquids liquid1 k':        -0.008496
	'water-1-butanol, two liquids liquid1 sigma':    -0.3779
	'water-1-butanol, two liquids liquid2 mu':       -0.4153
	'water-1-butanol, two liquids liquid2 k':        1.81
	'water-1-butanol, two liquids liquid2 sigma':    -0.02007
}

fn value(v ?f64) f64 {
	return v or { math.nan() }
}

// Each number DWSIM reports is its formulas', and its formulas on the
// files' values are its numbers: the data explain none of the differences.
fn test_its_numbers_are_its_formulas_and_its_data_the_files() {
	f := dw()!
	assert f.dwsim.starts_with('DWSIM version 9.0.5')
	assert f.cases.len == 10
	mut phases := 0
	for c in f.cases {
		for ph in c.phases {
			at := '${c.label} ${ph.phase}'
			for r in [ph.reproduced.mu, ph.reproduced.k, ph.reproduced.sigma] {
				if d := r {
					assert math.abs(d) <= 1e-15, '${at}: reproduced to ${d}'
				}
			}
			// A second liquid's viscosity and conductivity are the first's.
			if ph.phase == 'liquid2' {
				first := c.phases.filter(it.phase == 'liquid1')[0]
				assert value(ph.dwsim.mu) == value(first.dwsim.mu), at
				assert value(ph.dwsim.k) == value(first.dwsim.k), at
				assert value(ph.dwsim.sigma) != value(first.dwsim.sigma), at
				phases++
				continue
			}
			for j, pair in [[ph.on_files.mu, ph.dwsim.mu], [ph.on_files.k, ph.dwsim.k],
				[ph.on_files.sigma, ph.dwsim.sigma]] {
				if j == 2 && ph.phase == 'vapour' {
					continue // a vapour has no surface tension
				}
				if got := pair[0] {
					want := value(pair[1])
					assert math.abs(got / want - 1.0) <= 1e-5, '${at}: ${got} on the files, ${want} DWSIM'
				}
			}
			phases++
		}
	}
	assert phases == 13
}

// reflux's transport at DWSIM's phase states: its composition, T and P, and
// its density for the molar volume, so that only the methods differ.
fn reflux_at(ph DwPhase) thermo.PhaseLabel {
	return match ph.phase {
		'vapour' { thermo.PhaseLabel.vapor }
		'liquid2' { .liquid2 }
		else { .liquid }
	}
}

// Every difference from reflux's is the methods' (above): each recorded.
fn test_every_difference_is_the_methods() {
	f := dw()!
	mut measured := []string{}
	for c in f.cases {
		mix := components.mixture(c.cas.map(db.find(it)!))!
		for ph in c.phases {
			label := reflux_at(ph)
			v := ph.molar_mass / ph.density / 1000.0
			tr := mix.phase(label, c.t, c.p, ph.x, v)
			// A second liquid's DWSIM numbers are the first's: against its
			// own composition by DWSIM's formulas, which is on_files' here
			// for σ, its μ and λ are compared with reflux's directly.
			mut against := ph.on_files
			if ph.phase == 'liquid2' {
				against = DwValues{
					mu:    ph.dwsim.mu
					k:     ph.dwsim.k
					sigma: ph.on_files.sigma
				}
			}
			names := ['mu', 'k', 'sigma']
			ours := [tr.viscosity, tr.conductivity, tr.surface_tension]
			theirs := [against.mu, against.k, against.sigma]
			for i, name in names {
				t := theirs[i] or { continue }
				if ph.phase == 'vapour' && name == 'sigma' {
					continue
				}
				key := '${c.label} ${ph.phase} ${name}'
				if !ours[i].ok() {
					measured << '${key}: none (${ours[i].why}); DWSIM ${t}'
					continue
				}
				d := ours[i].value / t - 1.0
				measured << '${key}: ${d}'
				bound := method[key] or {
					assert false, '${key}: nothing recorded (${d})'
					continue
				}
				assert math.abs(d - bound) <= 1e-3 * math.abs(bound) + 1e-7, '${key}: ${d} against ${bound}'
			}
		}
	}
	println(measured.join('\n'))
	assert measured.len == method.len
}

// Its water, recorded and explained, not compared as a method, each against
// M5a's R12-08 and R15-11 at its own T and density:
// - the Steam Tables package's viscosity and conductivity are IAPS's 1985
//   formulations (its source): within 0.31% and 2.1% of IAPWS's here;
// - its CoolProp package's are CoolProp's, not ChemSep's as step 0 read the
//   source (those are 0.7% to 20% off them): its vapour's viscosity is
//   R12-08 at the state, within 1e-10 (1.3e-11 measured, against the
//   `iapws` package); its liquid's conductivity at 300 K R15-11 at the
//   saturated liquid's density, to 1e-13 (by the `iapws` package), and its
//   liquid's viscosity near the saturated liquid's times a pressure factor
//   (1e-5 at 300 K, 1.6% at 500 K and 10 MPa, not explained); in all,
//   within 2.7% and 1.1% of IAPWS's at the state. Its liquid at 500 K and
//   10 MPa is the 2.7%: for the report to DWSIM.
fn test_its_water_against_iapws() {
	f := dw()!
	for pkg, rows in {
		'steam tables': f.water.steam_tables
		'coolprop':     f.water.coolprop
	} {
		mut worst_mu, mut worst_k := 0.0, 0.0
		for r in rows {
			mu := water.viscosity(r.t, r.density)!
			k := water.thermal_conductivity(r.t, r.density)!
			if pkg == 'coolprop' && r.phase == 'vapour' {
				assert math.abs(r.mu / mu - 1.0) <= 1e-10, '${r.t} K: ${r.mu} against ${mu}'
			}
			worst_mu = math.max(worst_mu, math.abs(r.mu / mu - 1.0))
			worst_k = math.max(worst_k, math.abs(r.k / k - 1.0))
		}
		println('${pkg}: viscosity within ${worst_mu}, conductivity within ${worst_k} of IAPWS')
		bounds := if pkg == 'coolprop' { [0.027, 0.011] } else { [0.0031, 0.021] }
		assert worst_mu <= bounds[0] && worst_k <= bounds[1], pkg
		assert rows.len == 5
	}
}
