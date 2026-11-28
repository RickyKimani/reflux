// Measured mixtures (M5b step 9, criterion 5): each phase's transport, by
// D4 to D9's methods on ChemSep's pure values, against data sets from the
// ThermoML archive (testdata/m5b_published.json, from
// tools/gen_m5b_published.py), each set's average and largest deviation and
// its points with no value recorded, which can't grow (`recorded`), with
// each alternative's beside it, printed; and ChemSep's water against
// IAPWS's.
import components
import flash
import json2
import math
import os
import thermo
import water

// The fixture's types, named so that no other module's type shares a name:
// json2 knows a type by its short name (docs/v-for-go.md).
struct MeasuredPoint {
	t     f64 @[json: 'T']
	p     f64 @[json: 'P']
	x     []f64
	value f64
	u     ?f64
}

struct MeasuredAlternative {
	method  string
	points  int
	refused int
	aad     ?f64
	bias    ?f64
	max     ?f64
}

struct MeasuredCompound {
	name string
	cas  string
}

struct MeasuredSet {
	key            string
	doi            string
	citation       string
	property       string
	phase          string
	compounds      []MeasuredCompound
	pressure_given bool
	points         []MeasuredPoint
	alternatives   []MeasuredAlternative
}

struct MeasuredFile {
	sets []MeasuredSet
}

const db = components.load_default() or { panic(err) }

fn measured() !MeasuredFile {
	return json2.decode[MeasuredFile](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm5b_published.json'))!)!
}

// The deviations, recorded: each as measured, rounded up to two figures,
// so that a change that makes a method worse fails here; and the points a
// set has no value at. Keyed by the set and the statistic.
const recorded = {
	'H2 + CO2 gas viscosity, 300 kPa aad':                             0.026
	'H2 + CO2 gas viscosity, 300 kPa max':                             0.068
	'H2 + CO2 gas viscosity, 300 kPa refused':                         0.0
	'CO2 + CH4 gas viscosity, to 32 MPa aad':                          0.063
	'CO2 + CH4 gas viscosity, to 32 MPa max':                          0.21
	'CO2 + CH4 gas viscosity, to 32 MPa refused':                      7.0
	'CH4 + C3H8 gas viscosity, to 31 MPa aad':                         0.028
	'CH4 + C3H8 gas viscosity, to 31 MPa max':                         0.072
	'CH4 + C3H8 gas viscosity, to 31 MPa refused':                     0.0
	'CH4 + C4H10 gas viscosity, to 30 MPa aad':                        0.023
	'CH4 + C4H10 gas viscosity, to 30 MPa max':                        0.079
	'CH4 + C4H10 gas viscosity, to 30 MPa refused':                    0.0
	'CO2 gas conductivity, to 15 MPa aad':                             0.033
	'CO2 gas conductivity, to 15 MPa max':                             0.23
	'CO2 gas conductivity, to 15 MPa refused':                         0.0
	'CH4 + CO2 gas conductivity, to 12 MPa aad':                       0.016
	'CH4 + CO2 gas conductivity, to 12 MPa max':                       0.051
	'CH4 + CO2 gas conductivity, to 12 MPa refused':                   0.0
	'heptane liquid viscosity, to 10 MPa aad':                         0.053
	'heptane liquid viscosity, to 10 MPa max':                         0.098
	'heptane liquid viscosity, to 10 MPa refused':                     0.0
	'octane liquid viscosity, to 10 MPa aad':                          0.043
	'octane liquid viscosity, to 10 MPa max':                          0.086
	'octane liquid viscosity, to 10 MPa refused':                      0.0
	'heptane + octane liquid viscosity, to 10 MPa aad':                0.06
	'heptane + octane liquid viscosity, to 10 MPa max':                0.11
	'heptane + octane liquid viscosity, to 10 MPa refused':            0.0
	'heptane liquid viscosity, to 513 K aad':                          0.043
	'heptane liquid viscosity, to 513 K max':                          0.076
	'heptane liquid viscosity, to 513 K refused':                      0.0
	'octane liquid viscosity, to 513 K aad':                           0.032
	'octane liquid viscosity, to 513 K max':                           0.055
	'octane liquid viscosity, to 513 K refused':                       0.0
	'heptane + octane liquid viscosity, to 513 K aad':                 0.048
	'heptane + octane liquid viscosity, to 513 K max':                 0.073
	'heptane + octane liquid viscosity, to 513 K refused':             0.0
	'ethanol liquid viscosity, to pressure aad':                       0.16
	'ethanol liquid viscosity, to pressure max':                       0.26
	'ethanol liquid viscosity, to pressure refused':                   0.0
	'ethanol + heptane liquid viscosity, 1 atm aad':                   0.15
	'ethanol + heptane liquid viscosity, 1 atm max':                   0.26
	'ethanol + heptane liquid viscosity, 1 atm refused':               0.0
	'ethanol + heptane liquid viscosity, to pressure aad':             0.11
	'ethanol + heptane liquid viscosity, to pressure max':             0.26
	'ethanol + heptane liquid viscosity, to pressure refused':         0.0
	'1-butanol + 2-methylheptane liquid viscosity aad':                0.2
	'1-butanol + 2-methylheptane liquid viscosity max':                0.44
	'1-butanol + 2-methylheptane liquid viscosity refused':            0.0
	'methanol + water liquid viscosity aad':                           0.41
	'methanol + water liquid viscosity max':                           0.61
	'methanol + water liquid viscosity refused':                       0.0
	'diethylene glycol liquid viscosity aad':                          0.032
	'diethylene glycol liquid viscosity max':                          0.05
	'diethylene glycol liquid viscosity refused':                      0.0
	'diethylene glycol + water liquid viscosity aad':                  0.46
	'diethylene glycol + water liquid viscosity max':                  0.74
	'diethylene glycol + water liquid viscosity refused':              0.0
	'diethylethanolamine + water liquid viscosity aad':                0.72
	'diethylethanolamine + water liquid viscosity max':                0.91
	'diethylethanolamine + water liquid viscosity refused':            0.0
	'2-butoxyethanol liquid conductivity, to 150 MPa aad':             0.081
	'2-butoxyethanol liquid conductivity, to 150 MPa max':             0.11
	'2-butoxyethanol liquid conductivity, to 150 MPa refused':         0.0
	'2-butoxyethanol + water liquid conductivity, to 150 MPa aad':     0.039
	'2-butoxyethanol + water liquid conductivity, to 150 MPa max':     0.082
	'2-butoxyethanol + water liquid conductivity, to 150 MPa refused': 0.0
	'CH4 + C2H6 surface tension aad':                                  0.13
	'CH4 + C2H6 surface tension max':                                  0.67
	'CH4 + C2H6 surface tension refused':                              28.0
	'O2 surface tension aad':                                          0.0055
	'O2 surface tension max':                                          0.013
	'O2 surface tension refused':                                      0.0
	'N2 surface tension aad':                                          0.023
	'N2 surface tension max':                                          0.074
	'N2 surface tension refused':                                      0.0
	'N2 + O2 surface tension aad':                                     0.073
	'N2 + O2 surface tension max':                                     0.44
	'N2 + O2 surface tension refused':                                 12.0
	'diethylethanolamine + water surface tension aad':                 0.33
	'diethylethanolamine + water surface tension max':                 0.69
	'diethylethanolamine + water surface tension refused':             0.0
	'water liquid conductivity to 550 K':                              0.02
	'water liquid conductivity to 630 K':                              0.018
	'water liquid surface tension to 550 K':                           0.0065
	'water liquid surface tension to 630 K':                           0.069
	'water liquid viscosity to 550 K':                                 0.021
	'water liquid viscosity to 630 K':                                 0.058
	'water vapour conductivity at 1 kPa':                              0.014
	'water vapour viscosity at 1 kPa':                                 0.018
}

// within asserts a deviation is no larger than recorded.
fn within(key string, value f64) {
	bound := recorded[key] or {
		assert false, '${key}: nothing recorded (${value})'
		return
	}
	assert value <= bound, '${key}: ${value} against ${bound} recorded'
}

// words is a reason up to its first number, so that reasons differing only
// in their numbers are one.
fn words(why string) string {
	for i, c in why {
		if c.is_digit() {
			return why[..i].trim_right(' (,=')
		}
	}
	return why
}

// Stats are a method's deviations over a set, as fractions of the measured
// values.
struct Stats {
	aad     f64
	bias    f64
	worst   f64 // the deviation of largest size, signed
	points  int
	refused int
	why     []string // each different reason for a point with no value
}

fn stats_of(devs []f64, refused int, why []string) Stats {
	mut sum, mut abs := 0.0, 0.0
	mut worst := 0.0
	for d in devs {
		sum += d
		abs += math.abs(d)
		if math.abs(d) > math.abs(worst) {
			worst = d
		}
	}
	n := f64(math.max(devs.len, 1))
	return Stats{abs / n, sum / n, worst, devs.len, refused, why}
}

// reflux_set is M5b's methods over a set: each point's phase at its label's
// volume root on reflux's Peng–Robinson (ChemSep's kᵢⱼ), as a flowsheet's
// stream would have it, then the mixture's transport there.
fn reflux_set(s MeasuredSet) !Stats {
	cas := s.compounds.map(it.cas)
	fluid := db.setup(.pr, cas)!.fluid
	mix := components.mixture(cas.map(db.find(it)!))!
	label := if s.phase == 'vapour' { thermo.PhaseLabel.vapor } else { thermo.PhaseLabel.liquid }
	root := if s.phase == 'vapour' { thermo.Root.vapor } else { thermo.Root.liquid }
	mut devs := []f64{}
	mut refused := 0
	mut why := []string{}
	for pt in s.points {
		st := fluid.model.phase(pt.t, pt.p, pt.x, root) or {
			refused++
			if words(err.msg()) !in why {
				why << words(err.msg())
			}
			continue
		}
		mut n := 0.0
		for xi in pt.x {
			n += xi
		}
		tr := mix.phase(label, pt.t, pt.p, pt.x, st.v / n)
		r := match s.property {
			'viscosity' { tr.viscosity }
			'thermal_conductivity' { tr.conductivity }
			else { tr.surface_tension }
		}
		if !r.ok() {
			refused++
			// A reason's numbers differ point to point: keep its words.
			if words(r.why) !in why {
				why << words(r.why)
			}
			continue
		}
		devs << r.value / pt.value - 1.0
	}
	return stats_of(devs, refused, why)
}

fn percent(x f64) string {
	return '${math.round(x * 1000.0) / 10.0}%'
}

// Criterion 5: every set's deviations, M5b's methods within what's
// recorded, each alternative's printed beside them.
fn test_each_set_against_the_data() {
	f := measured()!
	assert f.sets.len == 27
	mut points := 0
	for s in f.sets {
		got := reflux_set(s)!
		points += s.points.len
		assert got.points + got.refused == s.points.len
		println('${s.key} (${s.points.len} points; ${s.doi}):')
		println('    M5b: ${percent(got.aad)} average, ${percent(got.bias)} bias, ${percent(got.worst)} at worst, ${got.refused} with no value ${got.why}')
		for a in s.alternatives {
			if aad := a.aad {
				println('    ${a.method}: ${percent(aad)} average, ${percent(a.bias or { 0.0 })} bias, ${percent(a.max or {
					0.0
				})} at worst, ${a.refused} with no value')
			} else {
				println('    ${a.method}: no value at any point')
			}
		}
		within('${s.key} aad', got.aad)
		within('${s.key} max', math.abs(got.worst))
		within('${s.key} refused', got.refused)
	}
	assert points == 1662
}

// ChemSep's water against IAPWS's (D10's measurement, step 0's again in V):
// its liquid's viscosity, conductivity and surface tension on IAPWS-95's
// saturated liquid from 275 to 630 K, and its vapour's at 1 kPa from 300
// to 1073 K, against M5a's functions (R12-08, R15-11, R1-76, each held to
// IAPWS's tables), the largest deviation in each range recorded.
fn test_chemseps_water_against_iapws() {
	w := db.find('water')!.transport_pure()!
	iapws := db.iapws95('water')!.fluid
	mut worst := map[string]f64{}
	for k in 0 .. 72 {
		t := 275.0 + 5.0 * f64(k)
		sat := flash.tvf(iapws.model, t, 0.0, [1.0])!
		liq := sat.phases.filter(it.label == .liquid)[0].state
		rho := water.molar_mass / liq.v
		range := if t <= 550.0 { 'to 550 K' } else { 'to 630 K' }
		for name, pair in {
			'viscosity':       [w.value(.liquid_viscosity, t)!.value, water.viscosity(t, rho)!]
			'conductivity':    [w.value(.liquid_conductivity, t)!.value, water.thermal_conductivity(t,
				rho)!]
			'surface tension': [w.value(.surface_tension, t)!.value, water.surface_tension(t)!]
		} {
			key := 'liquid ${name} ${range}'
			worst[key] = math.max(worst[key], math.abs(pair[0] / pair[1] - 1.0))
		}
	}
	for k in 0 .. 78 {
		t := 300.0 + 10.0 * f64(k)
		st := iapws.model.phase(t, 1e3, [1.0], .vapor)!
		rho := water.molar_mass / st.v
		for name, pair in {
			'viscosity':    [w.value(.vapour_viscosity, t)!.value, water.viscosity(t, rho)!]
			'conductivity': [w.value(.vapour_conductivity, t)!.value, water.thermal_conductivity(t,
				rho)!]
		} {
			key := 'vapour ${name} at 1 kPa'
			worst[key] = math.max(worst[key], math.abs(pair[0] / pair[1] - 1.0))
		}
	}
	mut keys := worst.keys()
	keys.sort()
	for key in keys {
		println('ChemSep water, ${key}: ${percent(worst[key])} at worst')
		within('water ${key}', worst[key])
	}
	assert keys.len == 8
}
