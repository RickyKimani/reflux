// The reactors against DWSIM 9 (M9 step 10, D13, criteria 7 and 8), from
// testdata/m9_dwsim.json (tools/gen_m9_dwsim.py): each reactor alone from
// the same inlet, and the owner's reforming PFR against DWSIM's answers and
// the exact answer of DWSIM's own equations. Every difference from DWSIM is
// taken apart here, or bounded by what the fixture shows DWSIM's own answer
// to be worth.
module reactors

import components
import expr
import flowsheet
import json2
import kinetics
import math
import os
import thermo

const db = components.load_default() or { panic(err) }

struct DwOutlet {
	t          f64
	p          f64
	flows      []f64
	duty       f64
	liquid_nan bool
}

struct DwsimParts {
	departure     f64
	ideal_gas     f64
	reaction_heat f64
}

struct ThermoParts {
	departure_pr                   f64
	departure_pr78_dwsim_constants f64
	ideal_gas_reflux_cp            f64
}

struct ConversionPoint {
	t      f64
	p      f64
	feed   []f64
	outlet DwOutlet
	extent f64
	dwsim  DwsimParts
	thermo ThermoParts
}

struct Adiabatic {
	feed   []f64
	outlet DwOutlet
}

struct ConversionPart {
	cas       []string
	points    []ConversionPoint
	adiabatic Adiabatic
}

struct Answer {
	outlet DwOutlet
}

struct EquilibriumPoint {
	t     f64
	p     f64
	feed  []f64
	bar   []f64
	atm   []f64
	dflt  Answer @[json: 'default']
	tight Answer
}

struct EquilibriumPart {
	cas    []string
	points []EquilibriumPoint
}

struct GibbsAnswer {
	outlet            DwOutlet
	element_imbalance f64
	ln_q_over_k_atm   []f64
}

struct GibbsDwsim {
	ipopt GibbsAnswer
	other GibbsAnswer
}

struct GibbsPoint {
	t     f64
	p     f64
	feed  []f64
	dwsim GibbsDwsim
	bar   []f64
	atm   []f64
}

struct GibbsPart {
	cas    []string
	points []GibbsPoint
}

struct CstrDwsim {
	reacted f64
}

struct CstrPart {
	dwsim   CstrDwsim
	cas     []string
	feed    []f64
	t       f64
	p       f64
	volume  f64
	loading f64
	exact   []f64
}

struct PfrDwsim {
	dv    f64
	ergun bool
	flows []f64
	x     f64
	drop  f64
	duty  f64
}

struct Exact {
	flows []f64
	x     f64
	drop  f64
}

struct PfrPart {
	feed           []f64
	t              f64
	p              f64
	cas            []string
	volume         f64
	length         f64
	loading        f64
	void_fraction  f64 @[json: 'void']
	particle       f64
	dwsim          []PfrDwsim
	exact_isobaric Exact
	exact_ergun    Exact
}

struct IdealGasPart {
	cas    []string
	dwsim  []f64
	reflux []f64
}

struct Fixture {
	conversion  ConversionPart
	equilibrium EquilibriumPart
	gibbs       GibbsPart
	cstr        CstrPart
	pfr         PfrPart
	ideal_gas   IdealGasPart
}

fn load_fixture() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_dwsim.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

// system is the species and a Peng–Robinson fluid with reflux's own kᵢⱼ,
// over compounds by CAS number.
fn system(cas []string) ([]kinetics.Species, thermo.Fluid, []components.Record) {
	records := cas.map(db.find(it) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name)) or { panic(err) }
	return records.map(flowsheet.species_of(it)), fl, records
}

fn worst_relative(a []f64, b []f64) f64 {
	mut w := 0.0
	for i, x in a {
		w = math.max(w, math.abs(x - b[i]) / math.max(math.abs(b[i]), 1e-300))
	}
	return w
}

// The owner's Tutorial 2 (methanol synthesis in a conversion reactor, M9
// step 10), its conversion typed as the tutorial prints it, at every point
// of its two sweeps: the flows are DWSIM's to rounding, both computing the
// same conversion. The duties differ by up to 4.1e-5, and the difference
// is taken apart point by point: DWSIM's departures less reflux's (thermo's
// PR with exact constants, which reflux matches), and DWSIM's ideal-gas
// enthalpy less an exact quadrature of reflux's Cp. Measured: what's left
// is 1.6e-8 of the duty at most, rounding of the parts. Of DWSIM's
// departures, PR78's α for methanol is 6.0 W at the base case and DWSIM's
// rounded constants −0.15 W, both from thermo's PR78; 3.0 W is DWSIM's own,
// beyond thermo's PR78 with its constants, and not traced.
fn test_conversion_against_dwsim() {
	f := load_fixture()!
	c := f.conversion
	sp, fl, _ := system(c.cas)
	rx := kinetics.parse_equation('synthesis', 'CO + 2 H2 = CH3OH', sp)!
	ctx := kinetics.rate_context(sp, [rx])
	defs := expr.Definitions{}
	e := expr.parse(ctx, &defs, '0.85 - 0.0018*(T in degC)')!
	r := new_conversion_reactor('t2',
		species: sp
		steps:   [Conversion{
			reaction:   rx
			key:        0
			expression: e
		}]
		ctx:     ctx
		defs:    defs
	)!
	mut worst_flow, mut worst_duty, mut worst_left := 0.0, 0.0, 0.0
	for pt in c.points {
		inlet := flowsheet.feed_tp(fl, pt.feed, pt.t, pt.p)!
		o := r.calculate(fl, [inlet])!
		flows := o.outlets[0].flows
		duty := o.energy[0].rate
		worst_flow = math.max(worst_flow, worst_relative(flows, pt.outlet.flows))
		gap := duty - pt.outlet.duty
		worst_duty = math.max(worst_duty, math.abs(gap / pt.outlet.duty))
		explained := (pt.thermo.departure_pr - pt.dwsim.departure) +
			(pt.thermo.ideal_gas_reflux_cp - pt.dwsim.ideal_gas)
		left := math.abs(gap - explained) / math.abs(pt.outlet.duty)
		worst_left = math.max(worst_left, left)
		println('${pt.t} K, ${pt.p / 1e5} bar: duty ${duty:.3f} W, DWSIM ${pt.outlet.duty:.3f}; gap ${gap:.3f} W: departures ${pt.thermo.departure_pr - pt.dwsim.departure:.3f} (PR78 and its constants ${pt.thermo.departure_pr - pt.thermo.departure_pr78_dwsim_constants:.3f}), ideal gas ${pt.thermo.ideal_gas_reflux_cp - pt.dwsim.ideal_gas:.3f}; left ${gap - explained:.4f}; warnings ${o.warnings}')
		// The same heat of formation in both.
		assert math.abs(pt.dwsim.reaction_heat - (-90410.0)) <= 1e-9
		// Past equilibrium below 70 bar at 320 °C, and not at 70 bar.
		past := o.warnings.any(it.contains('past equilibrium'))
		assert past == (pt.p < 7e6 && pt.t == 593.15), '${pt.t} K, ${pt.p} Pa: ${o.warnings}'
	}
	println('Tutorial 2: flows within ${worst_flow:.1e}, duties within ${worst_duty:.1e}, ${worst_left:.1e} of the duty unexplained')
	assert worst_flow <= 1e-14
	assert worst_duty <= 5e-5
	assert worst_left <= 1e-7
}

// The constructor refuses a conversion expression it can't use, as the
// file's reader does before it: one with no context to read T from, one
// that isn't a pure number, and one that uses P.
fn test_a_conversion_expression_the_constructor_refuses() {
	sp, _, _ := system(['630-08-0', '1333-74-0', '67-56-1'])
	rx := kinetics.parse_equation('synthesis', 'CO + 2 H2 = CH3OH', sp)!
	ctx := kinetics.rate_context(sp, [rx])
	defs := expr.Definitions{}
	cases := [
		['0.85 - 0.0018*(T in degC)', 'no context', 'with no context to read T from'],
		['0.0018 * T', '', 'must be a pure number'],
		['0.003*(P in bar)', '', 'may depend on T alone, and this uses P'],
	]
	for c in cases {
		e := expr.parse(ctx, &defs, c[0]) or { panic(err) }
		new_conversion_reactor('t2',
			species: sp
			steps:   [Conversion{
				reaction:   rx
				key:        0
				expression: e
			}]
			ctx:     if c[1] == 'no context' { expr.Context{} } else { ctx }
			defs:    defs
		) or {
			assert err.msg().contains(c[2]), '${c[0]}: ${err.msg()}'
			continue
		}
		assert false, '${c[0]} accepted'
	}
}

// DWSIM's adiabatic conversion reactor evaluates its expression at the
// inlet's T: 27.4% at an outlet of 783.3 K, where the correlation gives
// −6.8%. reflux evaluates it at the outlet's (D10, as amended before step
// 10), where it finds a consistent outlet; given DWSIM's fixed 27.4%, it
// lands on DWSIM's outlet to its energy balance's difference, and warns
// that the outlet is far past equilibrium.
fn test_adiabatic_conversion_against_dwsim() {
	f := load_fixture()!
	c := f.conversion
	sp, fl, _ := system(c.cas)
	rx := kinetics.parse_equation('synthesis', 'CO + 2 H2 = CH3OH', sp)!
	ctx := kinetics.rate_context(sp, [rx])
	defs := expr.Definitions{}
	e := expr.parse(ctx, &defs, '0.85 - 0.0018*(T in degC)')!
	inlet := flowsheet.feed_tp(fl, c.adiabatic.feed, 593.15, 70e5)!
	at_outlet := new_conversion_reactor('t2',
		species: sp
		steps:   [Conversion{
			reaction:   rx
			key:        0
			expression: e
		}]
		ctx:     ctx
		defs:    defs
		thermal: .adiabatic
	)!
	o := at_outlet.calculate(fl, [inlet])!
	out := o.outlets[0]
	x := (inlet.flows[0] - out.flows[0]) / inlet.flows[0]
	assert math.abs(x - (0.85 - 0.0018 * (out.t - 273.15))) <= 1e-9, '${x} at ${out.t} K'
	assert out.t < c.adiabatic.outlet.t - 100
	fixed := new_conversion_reactor('t2',
		species: sp
		steps:   [Conversion{
			reaction:   rx
			key:        0
			conversion: 0.274
		}]
		thermal: .adiabatic
	)!
	g := fixed.calculate(fl, [inlet])!
	println('adiabatic: reflux at its outlet ${out.t:.4f} K, X ${x:.6f}; with 0.274 fixed ${g.outlets[0].t:.4f} K, DWSIM ${c.adiabatic.outlet.t:.4f} K; ${g.warnings}')
	assert math.abs(g.outlets[0].t - c.adiabatic.outlet.t) <= 0.05
	assert worst_relative(g.outlets[0].flows, c.adiabatic.outlet.flows) <= 1e-14
	assert g.warnings.any(it.contains('past equilibrium'))
}

// DWSIM's equilibrium reactor, at 1e-10, gives the equilibrium with each
// fugacity over 1 atm against K at 1 bar from the Gibbs energies: the
// fixture's atm solution, to 1e-5 in each flow (measured 2.9e-6, where
// DWSIM's loops stop at their tolerance of 1e-10). reflux's
// equilibrium reactor gives the fixture's bar solution, its own definition
// computed independently (thermo's PR, an exact quadrature of the
// component data), to 1e-9: so the whole difference, 0.19 to 0.55 points
// of conversion, is the standard state, 1.01325² on Q. At DWSIM's default
// tolerances (0.01 outside, 1e-3 inside) its answers move by up to 1.5e-2.
fn test_equilibrium_against_dwsim() {
	f := load_fixture()!
	eq := f.equilibrium
	sp, fl, _ := system(eq.cas)
	rx := kinetics.parse_equation('synthesis', 'CO + 2 H2 = CH3OH', sp)!
	r := new_equilibrium_reactor('eq', species: sp, reactions: [rx])!
	for pt in eq.points {
		inlet := flowsheet.feed_tp(fl, pt.feed, pt.t, pt.p)!
		o := r.calculate(fl, [inlet])!
		flows := o.outlets[0].flows
		ours := worst_relative(flows, pt.bar)
		theirs := worst_relative(pt.tight.outlet.flows, pt.atm)
		loose := worst_relative(pt.dflt.outlet.flows, pt.tight.outlet.flows)
		x := (pt.feed[0] - flows[0]) / pt.feed[0]
		xd := (pt.feed[0] - pt.tight.outlet.flows[0]) / pt.feed[0]
		println('${pt.p / 1e5} bar: reflux X ${x:.6f} (bar, within ${ours:.1e}), DWSIM ${xd:.6f} (atm, within ${theirs:.1e}; its default tolerances ${loose:.1e} away)')
		assert ours <= 1e-9
		assert theirs <= 1e-5
	}
}

// DWSIM's Gibbs reactor is good to about 1e-3. Its default solver, IPOPT,
// leaves carbon and oxygen unbalanced by 7.3e-4 at 1000 K and 1 atm (and
// 7.7e-5 at 900 K and 10 bar); its other method balances them to 3e-6, but
// the two disagree by 1e-3, and each is about 1e-3 from the equilibrium of
// fugacities over 1 atm (which explains DWSIM's equilibrium reactor
// exactly) and further from 1 bar's. reflux's Gibbs reactor gives the
// fixture's bar solution, independently computed, to 1e-9.
fn test_gibbs_against_dwsim() {
	f := load_fixture()!
	gp := f.gibbs
	sp, fl, _ := system(gp.cas)
	r := new_gibbs('gibbs', species: sp)!
	for pt in gp.points {
		inlet := flowsheet.feed_tp(fl, pt.feed, pt.t, pt.p)!
		o := r.calculate(fl, [inlet])!
		flows := o.outlets[0].flows
		ours := worst_relative(flows, pt.bar)
		ipopt := worst_relative(pt.dwsim.ipopt.outlet.flows, flows)
		other := worst_relative(pt.dwsim.other.outlet.flows, flows)
		println('${pt.t} K, ${pt.p / 1e5} bar: reflux within ${ours:.1e} of bar; DWSIM IPOPT ${ipopt:.1e} away (elements off ${pt.dwsim.ipopt.element_imbalance:.1e}, ln Q/K at 1 atm ${pt.dwsim.ipopt.ln_q_over_k_atm}), its other method ${other:.1e} (${pt.dwsim.other.element_imbalance:.1e}, ${pt.dwsim.other.ln_q_over_k_atm}); atm ${worst_relative(pt.atm, flows):.1e} away')
		assert ours <= 1e-9
		assert pt.dwsim.other.element_imbalance <= 1e-5
		assert ipopt <= 2e-2 && other <= 2e-2
	}
}

// reforming_kinetics is the owner's three rates as DWSIM has them typed
// (testdata/kinetics/dwsim_reforming.toml), and where each of its species
// is in a fixture's order of CAS numbers.
fn reforming_kinetics(cas []string) !(kinetics.Kinetics, []int) {
	kf := flowsheet.load_kinetics(os.join_path(@VMODROOT, 'testdata', 'kinetics', 'dwsim_reforming.toml'))!
	order := kf.kinetics.species.map(cas.index(it.cas))
	assert !order.contains(-1)
	return kf.kinetics, order
}

fn reorder(x []f64, order []int) []f64 {
	return order.map(x[it])
}

// DWSIM 9.0.5's CSTR reacts nothing in the vapour, catalytic or not, and
// says nothing (the fixture records it moving the feed by 9e-16 mol/s), so
// reflux's CSTR is held to the exact answer of DWSIM's own equations: the
// owner's rates as typed, isothermal at 1000 K and 1 atm, with 0.386 kg of
// catalyst in 1 m³, by mpmath at 30 digits. Newton's answer is within
// 1e-10 of each flow (measured below).
fn test_cstr_against_dwsims_equations() {
	f := load_fixture()!
	c := f.cstr
	assert c.dwsim.reacted <= 1e-12
	k, order := reforming_kinetics(c.cas)!
	records := k.species.map(db.find(it.cas) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name))!
	r := new_cstr('cstr',
		k:         k
		reactions: []int{len: k.reactions.len, init: index}
		volume:    c.volume
		loading:   c.loading
	)!
	inlet := flowsheet.feed_tp(fl, reorder(c.feed, order), c.t, c.p)!
	o := r.calculate(fl, [inlet])!
	w := worst_relative(o.outlets[0].flows, reorder(c.exact, order))
	println("CSTR: within ${w:.1e} of the exact answer of DWSIM's equations")
	assert w <= 1e-10
}

// DwsimViscosity is DWSIM's gas viscosity (M9 step 0): the mole-fraction
// mean of each compound's ChemSep correlation (eq. 102), corrected by Jossi,
// Stiel and Thodos with the mixture's Tc, Pc, Vc and M as mole-fraction
// means and its molar volume from the fluid's PR. The test gives it to the
// PFR, so that with Ergun the comparison measures the solver.
struct DwsimViscosity {
	model thermo.Model
	pure  []components.Property
	tc    []f64
	pc    []f64
	vc    []f64
	mw    []f64
}

fn new_dwsim_viscosity(fl thermo.Fluid, records []components.Record) !DwsimViscosity {
	return DwsimViscosity{
		model: fl.model
		pure:  records.map(it.property('vapor_viscosity')!)
		tc:    records.map(it.constants.tc or { 0.0 })
		pc:    records.map(it.constants.pc or { 0.0 })
		vc:    records.map(it.constants.vc or { 0.0 })
		mw:    records.map(it.constants.mw or { 0.0 })
	}
}

fn (d DwsimViscosity) viscosity(t f64, p f64, y []f64) !f64 {
	mut mean, mut tc, mut pc, mut vc, mut mm := 0.0, 0.0, 0.0, 0.0, 0.0
	for i, yi in y {
		mean += yi * d.pure[i].at(t).value
		tc += yi * d.tc[i]
		pc += yi * d.pc[i]
		vc += yi * d.vc[i]
		mm += yi * d.mw[i] * 1e3
	}
	st := thermo.solve_tp(d.model, t, p, y, .vapor)!
	xi := math.pow(tc / (mm * mm * mm * math.pow(pc / 101325, 4)), 1.0 / 6)
	rr := vc / st.v
	rhs := 1.023 + 0.23364 * rr + 0.58533 * rr * rr - 0.40758 * rr * rr * rr + 0.093324 * rr * rr * rr * rr
	return mean + (math.pow(rhs, 4) - 1) / xi * 1e-7
}

// The owner's reforming PFR (criterion 8, D13, step 0's targets), on DWSIM's
// setup: its three rates as typed, its catalyst and, with Ergun, DWSIM's
// viscosity. Against the exact answer of DWSIM's equations (SciPy's Radau
// at 1e-12, no sections): every outlet flow within 1e-8, isobaric and with
// Ergun, and the drop within 1e-3 Pa. Against DWSIM at dV = 0.001: each flow
// within 3e-4, methane's conversion within 5e-5, the drop within 1 Pa. At
// the file's dV = 0.05, DWSIM's sections put it further off: recorded.
fn test_the_owners_pfr_against_dwsim() {
	f := load_fixture()!
	pf := f.pfr
	k, order := reforming_kinetics(pf.cas)!
	records := k.species.map(db.find(it.cas) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name))!
	ch4 := k.species.map(it.cas).index('74-82-8')
	inlet := flowsheet.feed_tp(fl, reorder(pf.feed, order), pf.t, pf.p)!
	all := []int{len: k.reactions.len, init: index}
	conversion_of := fn [inlet, ch4] (flows []f64) f64 {
		return (inlet.flows[ch4] - flows[ch4]) / inlet.flows[ch4]
	}
	visc := new_dwsim_viscosity(fl, records)!
	wilke := new_wilke(records)!
	for ergun in [false, true] {
		r := new_pfr('pfr',
			k:         k
			reactions: all
			volume:    pf.volume
			length:    pf.length
			catalyst:  Catalyst{pf.loading, pf.void_fraction, pf.particle}
			drop:      if ergun { PressureDrop.ergun } else { PressureDrop.no_drop }
			viscosity: visc
			rtol:      1e-10
		)!
		o := r.calculate(fl, [inlet])!
		flows := o.outlets[0].flows
		drop := inlet.p - o.outlets[0].p
		x := conversion_of(flows)
		exact := if ergun { pf.exact_ergun } else { pf.exact_isobaric }
		w := worst_relative(flows, reorder(exact.flows, order))
		println('Ergun ${ergun}: X ${x:.10f}, drop ${drop:.6f} Pa, duty ${o.energy[0].rate:.3f} W; the exact answer within ${w:.1e}, its drop ${drop - exact.drop:.2e} Pa away')
		assert w <= 1e-8
		if ergun {
			assert math.abs(drop - exact.drop) <= 1e-3
		}
		for dw in pf.dwsim.filter(it.ergun == ergun) {
			wd := worst_relative(flows, reorder(dw.flows, order))
			println('  DWSIM at dV ${dw.dv}: flows ${wd:.2e} away, X ${x - dw.x:+.2e}, drop ${drop - dw.drop:+.4f} Pa, duty ${o.energy[0].rate - dw.duty:+.1f} W')
			if dw.dv == 0.001 {
				assert wd <= 3e-4
				assert math.abs(x - dw.x) <= 5e-5
				assert math.abs(drop - dw.drop) <= 1.0
				// The duty: reflux's enthalpies of DWSIM's own outlet give
				// DWSIM's duty but for each compound's ideal-gas enthalpy,
				// which DWSIM's copy of ChemSep has with water's Cp as eq 16
				// where reflux's file has eq 107 (72 J/mol at 1000 K).
				out_dw := flowsheet.feed_tp(fl, reorder(dw.flows, order), pf.t, inlet.p - dw.drop)!
				at_dw := out_dw.total() * out_dw.h - inlet.total() * inlet.h
				mut explained := 0.0
				for i, cas in f.ideal_gas.cas {
					made := dw.flows[i] - pf.feed[pf.cas.index(cas)]
					explained += made * (f.ideal_gas.reflux[i] - f.ideal_gas.dwsim[i])
				}
				left := at_dw - dw.duty - explained
				println("    its duty: reflux on its outlet ${at_dw - dw.duty:+.2f} W from it, the ideal gas ${explained:+.2f}, left ${left:+.3f} W; from its outlet to reflux's ${o.energy[0].rate - at_dw:+.2f} W")
				assert math.abs(left) <= 1e-5 * math.abs(dw.duty)
			}
		}
		if ergun {
			// reflux's own viscosity, Wilke's rule, beside it.
			own := new_pfr('pfr',
				k:         k
				reactions: all
				volume:    pf.volume
				length:    pf.length
				catalyst:  Catalyst{pf.loading, pf.void_fraction, pf.particle}
				drop:      .ergun
				viscosity: wilke
				rtol:      1e-10
			)!
			oo := own.calculate(fl, [inlet])!
			println("  with reflux's viscosity (Wilke's rule): X ${conversion_of(oo.outlets[0].flows):.10f}, drop ${inlet.p - oo.outlets[0].p:.4f} Pa")
		}
	}
}
