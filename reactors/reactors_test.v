// The reactors (M9 step 9, D10, criterion 7): closed forms for the PFR and
// the CSTR, the Gibbs and equilibrium reactors against mpmath's equilibria
// (testdata/m9_equilibrium.json), a conversion reactor's sequence, energy
// balances, and Ergun's drop.
module reactors

import components
import flowsheet
import json2
import kinetics
import math
import os
import thermo

const db = components.load_default() or { panic(err) }

// setup is species, a Peng–Robinson fluid with every kᵢⱼ zero, and the
// reactions written, over compounds by CAS number.
fn setup(cas []string, equations []string) ([]kinetics.Species, thermo.Fluid, []kinetics.Reaction) {
	records := cas.map(db.find(it) or { panic(err) })
	sp := records.map(flowsheet.species_of(it))
	fl := db.fluid(.pr, records.map(it.name), kij: false) or { panic(err) }
	mut rs := []kinetics.Reaction{}
	for i, e in equations {
		rs << kinetics.parse_equation('r${i + 1}', e, sp) or { panic(err) }
	}
	return sp, fl, rs
}

// first_order is a kinetics of one irreversible rate, k·p of a species, per
// volume.
fn first_order(sp []kinetics.Species, rs []kinetics.Reaction, k string, species string) kinetics.Kinetics {
	return kinetics.build(sp, rs, [], [
		kinetics.Rated{0, kinetics.PowerLaw{
			k:      kinetics.Fixed{k}
			orders: [kinetics.Power{species, 1}]
			per:    .volume
		}},
	], []) or { panic(err) }
}

fn feed(fl thermo.Fluid, flows []f64, t f64, p f64) flowsheet.Stream {
	return flowsheet.feed_tp(fl, flows, t, p) or { panic(err) }
}

// An isomerisation, n-butane to isobutane, first order in p: the moles
// don't change, so dF_A/dV = −kP·F_A/F, and F_A = F₀·exp(−kPV/F) in a PFR
// and F₀/(1 + kPV/F) in a CSTR, kPV/F = 2 here. The PFR's outlet is within
// 1e-9 of it at rtol 1e-10, both integrators (measured: 2.0e-11 and
// 3.2e-11), and the CSTR's within 1e-12, Newton's.
fn test_first_order_closed_forms() {
	sp, fl, rs := setup(['106-97-8', '75-28-5'], ['n-butane = isobutane'])
	k := first_order(sp, rs, '2e-5 mol/(m^3*s*Pa)', '106-97-8')
	inlet := feed(fl, [1.0, 0], 400, 1e5)
	for integ in [Integrator.rodas4, .dopri5] {
		r := new_pfr('pfr', k: k, reactions: [0], volume: 1, rtol: 1e-10, integrator: integ)!
		o := r.calculate(fl, [inlet])!
		f := o.outlets[0].flows
		println('${integ}: ${f[0] / math.exp(-2) - 1}')
		assert math.abs(f[0] / math.exp(-2) - 1) <= 1e-9, '${integ}: ${f[0]} vs ${math.exp(-2)}'
		assert math.abs(f[0] + f[1] - 1) <= 1e-12
	}
	c := new_cstr('cstr', k: k, reactions: [0], volume: 1)!
	o := c.calculate(fl, [inlet])!
	assert math.abs(o.outlets[0].flows[0] - 1.0 / 3) <= 1e-12, '${o.outlets[0].flows}'
}

// Ethane to ethylene and hydrogen, first order in p, from pure ethane: the
// moles grow, F = F₀(1 + X), and dX/dV = kP(1 − X)/(F₀(1 + X)), so
// −X − 2 ln(1 − X) = kPV/F₀, 1.5 here, solved for X by Newton.
fn test_moles_that_change() {
	sp, fl, rs := setup(['74-84-0', '74-85-1', '1333-74-0'], ['C2H6 = C2H4 + H2'])
	k := first_order(sp, rs, '1.5e-5 mol/(m^3*s*Pa)', 'C2H6')
	mut x := 0.5
	for _ in 0 .. 50 {
		x -= (-x - 2 * math.log(1 - x) - 1.5) / (-1 + 2 / (1 - x))
	}
	inlet := feed(fl, [1.0, 0, 0], 900, 1e5)
	for integ in [Integrator.rodas4, .dopri5] {
		r := new_pfr('pfr', k: k, reactions: [0], volume: 1, rtol: 1e-10, integrator: integ)!
		o := r.calculate(fl, [inlet])!
		got := 1 - o.outlets[0].flows[0]
		println('${integ}: X ${got} vs ${x}')
		assert math.abs(got / x - 1) <= 1e-9, '${integ}: ${got} vs ${x}'
	}
}

struct ECase {
	name      string
	species   []string
	model     string
	reactions []string
	inert     []int
	t         f64
	p         f64
	n0        []f64
	ln_y      []?f64
}

struct EFixture {
	cases []ECase
}

// The Gibbs reactor and the equilibrium reactor over chosen reactions land
// on mpmath's equilibria, Peng–Robinson's with every kᵢⱼ zero, to 1e-11 in
// ln y: the equilibrium module's answers, through a reactor's inlet and
// outlet. Each balances its elements and its energy, its duty the change in
// F·H.
fn test_against_the_equilibria() {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_equilibrium.json')
	f := json2.decode[EFixture](os.read_file(path)!)!
	mut done := 0
	for c in f.cases {
		if c.model != 'pr' || c.inert.len > 0 {
			continue
		}
		sp, fl, rs := setup(c.species, c.reactions)
		inlet := feed(fl, c.n0, c.t, c.p)
		unit := if rs.len == 0 {
			flowsheet.Unit(new_gibbs('gibbs', species: sp, t: c.t)!)
		} else {
			flowsheet.Unit(new_equilibrium_reactor('eq', species: sp, reactions: rs, t: c.t)!)
		}
		o := unit.calculate(fl, [inlet])!
		out := o.outlets[0]
		total := out.total()
		for j, want in c.ln_y {
			if w := want {
				assert math.abs(math.log(out.flows[j] / total) - w) <= 1e-11, '${c.name}, ${j}'
			}
		}
		check_balances(o, inlet)
		done++
	}
	assert done == 8
}

// check_balances holds a reactor's outcome to its elements, to 1e-12 of the
// flow, and its energy, to 1e-10 of F·RT.
fn check_balances(o flowsheet.Outcome, inlet flowsheet.Stream) {
	out := o.outlets[0]
	for row in o.conserved {
		mut x := 0.0
		for j, c in row {
			x += c * (inlet.flows[j] - out.flows[j])
		}
		assert math.abs(x) <= 1e-12 * inlet.total(), '${x}'
	}
	mut q := 0.0
	for e in o.energy {
		q += e.rate
	}
	gap := inlet.total() * inlet.h + q - out.total() * out.h
	assert math.abs(gap) <= 1e-10 * inlet.total() * thermo.r_gas * inlet.t, '${gap} W'
}

// Adiabatic, a Gibbs reactor's outlet carries the inlet's enthalpy, and is
// the equilibrium at the temperature that takes: reforming, endothermic,
// cools from 1100 K.
fn test_adiabatic_gibbs() {
	sp, fl, _ := setup(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'],
		[])
	inlet := feed(fl, [1.0, 3, 0.1, 0, 0], 1100, 5e5)
	g := new_gibbs('gibbs', species: sp, thermal: .adiabatic)!
	o := g.calculate(fl, [inlet])!
	out := o.outlets[0]
	assert out.t < 1100 && out.t > 600, '${out.t}'
	assert o.energy.len == 0
	check_balances(o, inlet)
	iso := new_gibbs('iso', species: sp, t: out.t)!
	same := iso.calculate(fl, [inlet])!
	for j in 0 .. 5 {
		assert math.abs(same.outlets[0].flows[j] - out.flows[j]) <= 1e-12 * inlet.total()
	}
}

// A conversion reactor runs its reactions in order, each on what the last
// left: half the methane reformed, then 80% of the CO shifted, by hand.
fn test_conversion_in_order() {
	sp, fl, rs := setup(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'],
		['CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2'])
	c := new_conversion_reactor('conv',
		species: sp
		steps:   [Conversion{ reaction: rs[0], key: 0, conversion: 0.5 },
			Conversion{ reaction: rs[1], key: 3, conversion: 0.8 }]
		t:       900
	)!
	inlet := feed(fl, [2.0, 5, 0, 0, 0], 900, 2e5)
	o := c.calculate(fl, [inlet])!
	// Reforming: CH4 2 → 1, H2O 5 → 4, CO 1, H2 3. Shift: 0.8 of CO 1.
	want := [1.0, 4 - 0.8, 3 + 0.8, 1 - 0.8, 0.8]
	for j, w in want {
		assert math.abs(o.outlets[0].flows[j] - w) <= 1e-14 * 8, '${j}: ${o.outlets[0].flows[j]}'
	}
	check_balances(o, inlet)
	// Asking for more than there is is refused.
	greedy := new_conversion_reactor('greedy',
		species: sp
		steps:   [Conversion{ reaction: rs[0], key: 0, conversion: 1.0 }]
	)!
	greedy.calculate(fl, [feed(fl, [2.0, 1, 0, 0, 0], 900, 2e5)]) or {
		assert err.msg().contains('r1 needs more Water than there is'), err.msg()
		return
	}
	assert false
}

// Adiabatic, the PFR's and the CSTR's outlets carry the inlet's enthalpy,
// the PFR's integrated T agreeing with its outlet's flash (which calculate
// checks), on ethane's dehydrogenation, endothermic, slow enough to stay a
// gas: its rate here doesn't fall with T.
fn test_adiabatic_energy() {
	sp, fl, rs := setup(['74-84-0', '74-85-1', '1333-74-0'], ['C2H6 = C2H4 + H2'])
	k := first_order(sp, rs, '1.5e-6 mol/(m^3*s*Pa)', 'C2H6')
	inlet := feed(fl, [1.0, 0, 0], 1100, 1e5)
	r := new_pfr('pfr', k: k, reactions: [0], volume: 1, thermal: .adiabatic, rtol: 1e-10)!
	o := r.calculate(fl, [inlet])!
	println('adiabatic: PFR to ${o.outlets[0].t} K')
	assert o.outlets[0].t < 1100
	check_balances(o, inlet)
	c := new_cstr('cstr', k: k, reactions: [0], volume: 1, thermal: .adiabatic)!
	oc := c.calculate(fl, [inlet])!
	assert oc.outlets[0].t < 1100
	check_balances(oc, inlet)
}

// Mu is a viscosity that doesn't move, for Ergun's closed form.
struct Mu {
	value f64
}

fn (m Mu) viscosity(_ f64, _ f64, _ []f64) !f64 {
	return m.value
}

// Ergun's drop with nothing reacting (a rate too slow to show) and μ fixed:
// both terms go as 1/P, so P·dP/dz = −C and P² = P₀² − 2CL. C from the
// inlet's Z, which moves by under 1e-4 over the drop: the outlet's P within
// 1e-3 of the drop of it.
fn test_ergun() {
	sp, fl, rs := setup(['106-97-8', '75-28-5'], ['n-butane = isobutane'])
	k := first_order(sp, rs, '1e-30 mol/(m^3*s*Pa)', '106-97-8')
	inlet := feed(fl, [1.0, 0], 500, 2e5)
	cat := Catalyst{
		loading:  1000
		void:     0.4
		diameter: 2e-3
	}
	r := new_pfr('bed',
		k:         k
		reactions: [0]
		volume:    0.1
		length:    2
		catalyst:  cat
		drop:      .ergun
		viscosity: Mu{1.5e-5}
		rtol:      1e-10
	)!
	o := r.calculate(fl, [inlet])!
	p := o.outlets[0].p
	area := 0.1 / 2
	z := inlet.p * inlet.v / (thermo.r_gas * inlet.t)
	mw := 0.0581222
	g := inlet.total() * mw / area
	e3 := 0.4 * 0.4 * 0.4
	cc := 150 * 1.5e-5 * 0.6 * 0.6 / (e3 * 4e-6) * inlet.total() * z * thermo.r_gas * inlet.t / area +
		1.75 * 0.6 / (e3 * 2e-3) * g * g * z * thermo.r_gas * inlet.t / mw
	want := math.sqrt(inlet.p * inlet.p - 2 * cc * 2)
	println('ergun: ${p} vs ${want}, a drop of ${inlet.p - p} Pa')
	assert math.abs(p - want) <= 1e-3 * (inlet.p - p), '${p} vs ${want}'
	assert o.profile != none
}

// xu_froment is Xu and Froment's kinetics, read from its file.
fn xu_froment() kinetics.Kinetics {
	path := os.join_path(@VMODROOT, 'testdata', 'kinetics', 'xu_froment.toml')
	f := flowsheet.load_kinetics(path, db: db) or { panic(err) }
	return f.kinetics
}

// The sources' derivatives, by each flow, T and P, are their central
// differences', to 1e-6 of each row's largest term: Xu and Froment's three
// rates at 900 K and 10 bar, per mass of catalyst, in Peng–Robinson's gas
// (their rates on partial pressures don't move with φ or Z).
fn test_source_derivatives() {
	k := xu_froment()
	records := k.species.map(db.find(it.cas) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name), kij: false)!
	mut src := new_source('src', k, [0, 1, 2], 1000, fl.model)!
	flows := [1.0, 3, 0.3, 0.2, 0.1]
	t, p := 900.0, 10e5
	n := 5
	st := src.state(t, p, flows)!
	mut s := []f64{len: n}
	src.at(st, mut s)!
	mut dflows := []f64{len: n * n}
	mut dt := []f64{len: n}
	mut dp := []f64{len: n}
	src.derivs(st, numeric_sum(flows), mut dflows, mut dt, mut dp)
	mut sp := []f64{len: n}
	mut sm := []f64{len: n}
	for col in 0 .. n + 2 {
		mut fp := flows.clone()
		mut fm := flows.clone()
		mut tp, mut tm, mut pp, mut pm := t, t, p, p
		h := if col < n {
			1e-6 * flows[col]
		} else if col == n {
			1e-6 * t
		} else {
			1e-6 * p
		}
		if col < n {
			fp[col] += h
			fm[col] -= h
		} else if col == n {
			tp += h
			tm -= h
		} else {
			pp += h
			pm -= h
		}
		stp := src.state(tp, pp, fp)!
		src.at(stp, mut sp)!
		stm := src.state(tm, pm, fm)!
		src.at(stm, mut sm)!
		for i in 0 .. n {
			fd := (sp[i] - sm[i]) / (2 * h)
			exact := if col < n {
				dflows[i * n + col]
			} else if col == n {
				dt[i]
			} else {
				dp[i]
			}
			mut row := 0.0
			for kk in 0 .. n {
				row = math.max(row, math.abs(dflows[i * n + kk]) * flows[kk])
			}
			row = math.max(row, math.max(math.abs(dt[i]) * t, math.abs(dp[i]) * p))
			scale := if col < n {
				flows[col]
			} else if col == n {
				t
			} else {
				p
			}
			assert math.abs(exact - fd) * scale <= 1e-6 * row, 'S${i} by column ${col}: ${exact} vs ${fd}'
		}
	}
}

fn numeric_sum(x []f64) f64 {
	mut s := 0.0
	for v in x {
		s += v
	}
	return s
}

// The CSTR's Newton runs on its system's own Jacobian, which is its central
// differences', to 1e-6 of each row's largest term: adiabatic, so with the
// energy balance's row and T's column, Xu and Froment's rates at 950 K.
fn test_cstr_jacobian() {
	k := xu_froment()
	records := k.species.map(db.find(it.cas) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name), kij: false)!
	c := new_cstr('cstr',
		k:         k
		reactions: [0, 1, 2]
		volume:    0.1
		loading:   100
		thermal:   .adiabatic
	)!
	inlet := feed(fl, [1.0, 3, 0.1, 0, 0], 1000, 5e5)
	sys := CstrSystem{
		c:       c
		fl:      fl
		inlet:   inlet
		p:       5e5
		t_fixed: 1000
		target:  inlet.total() * inlet.h
		scale:   inlet.total()
		src:     new_source('cstr', k, [0, 1, 2], 100, fl.model)!
	}
	x := [0.7, 2.6, 1.0, 0.25, 0.05, 950]
	m := x.len
	mut f := []f64{len: m}
	mut jac := []f64{len: m * m}
	sys.eval(x, mut f, mut jac)!
	mut f1 := []f64{len: m}
	mut f2 := []f64{len: m}
	mut j1 := []f64{len: m * m}
	for col in 0 .. m {
		h := 1e-6 * x[col]
		mut xp := x.clone()
		mut xm := x.clone()
		xp[col] += h
		xm[col] -= h
		sys.eval(xp, mut f1, mut j1)!
		sys.eval(xm, mut f2, mut j1)!
		for i in 0 .. m {
			mut row := 0.0
			for kk in 0 .. m {
				row = math.max(row, math.abs(jac[i * m + kk]) * x[kk])
			}
			fd := (f1[i] - f2[i]) / (2 * h)
			assert math.abs(jac[i * m + col] - fd) * x[col] <= 1e-6 * row, 'row ${i}, column ${col}: ${jac[i * m + col]} vs ${fd}'
		}
	}
}

// A given drop is a line in V, so the isomerisation's ln(F_A/F₀) is
// −(k/F)·(P_in·V − ΔP·V/2): −1.8 here. The PFR, whose ∂f/∂V then isn't
// zero, lands within 1e-9 of it at rtol 1e-10.
fn test_a_given_drop() {
	sp, fl, rs := setup(['106-97-8', '75-28-5'], ['n-butane = isobutane'])
	k := first_order(sp, rs, '2e-5 mol/(m^3*s*Pa)', '106-97-8')
	r := new_pfr('pfr',
		k:             k
		reactions:     [0]
		volume:        1
		drop:          .given
		pressure_drop: 2e4
		rtol:          1e-10
	)!
	o := r.calculate(fl, [feed(fl, [1.0, 0], 400, 1e5)])!
	assert math.abs(o.outlets[0].flows[0] / math.exp(-1.8) - 1) <= 1e-9, '${o.outlets[0].flows[0]}'
	assert o.outlets[0].p == 8e4
	profile := o.profile or { panic('no profile') }
	for pt in profile.points {
		assert math.abs(pt.p - (1e5 - 2e4 * pt.volume)) <= 1e-9
	}
}

// What a reactor refuses: an integration whose T and the energy balance's
// disagree (a rate that doesn't fall with T, cooling ethane past where it's
// a gas); an outlet that would split into two phases (methanol and water
// at 400 K and 50 bar); a specification that can't be; and an equilibrium
// over dependent reactions. An inert species passes through.
fn test_refused() {
	sp, fl, rs := setup(['74-84-0', '74-85-1', '1333-74-0'], ['C2H6 = C2H4 + H2'])
	k := first_order(sp, rs, '1.5e-5 mol/(m^3*s*Pa)', 'C2H6')
	r := new_pfr('pfr', k: k, reactions: [0], volume: 1, thermal: .adiabatic, rtol: 1e-10)!
	mut why := ''
	r.calculate(fl, [feed(fl, [1.0, 0, 0], 1100, 1e5)]) or { why = err.msg() }
	assert why.contains('by its energy balance and'), why
	ms, mfl, mrs := setup(['124-38-9', '1333-74-0', '67-56-1', '7732-18-5'], ['CO2 + 3 H2 = CH3OH + H2O'])
	c := new_conversion_reactor('conv',
		species: ms
		steps:   [Conversion{ reaction: mrs[0], key: 0, conversion: 0.5 }]
		t:       400
	)!
	why = ''
	c.calculate(mfl, [feed(mfl, [1.0, 3, 0, 0], 500, 50e5)]) or { why = err.msg() }
	assert why.contains('would split into two phases'), why
	gs, gfl, grs := setup(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'],
		['CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2', 'CH4 + 2 H2O = CO2 + 4 H2'])
	for bad in [
		fn [gs] () !flowsheet.Unit {
			return new_gibbs('g', species: gs, t: -5)
		},
		fn [gs] () !flowsheet.Unit {
			return new_gibbs('g', species: gs, thermal: .duty, duty: math.nan())
		},
		fn [gs] () !flowsheet.Unit {
			return new_gibbs('g', species: gs, drop: -1)
		},
		fn [gs] () !flowsheet.Unit {
			return new_gibbs('g', species: gs, inert: [7])
		},
		fn [gs, grs] () !flowsheet.Unit {
			return new_equilibrium_reactor('e', species: gs, reactions: grs)
		},
	] {
		bad() or { continue }
		assert false, 'accepted'
	}
	why = ''
	new_equilibrium_reactor('e', species: gs, reactions: grs) or { why = err.msg() }
	assert why.contains('must be independent, but r3 = r1 + r2'), why
	g := new_gibbs('g', species: gs, inert: [4], t: 900)!
	inlet := feed(gfl, [1.0, 2, 0.1, 0, 0.5], 900, 1e5)
	o := g.calculate(gfl, [inlet])!
	assert o.outlets[0].flows[4] == 0.5
}

// A reactor's conserved rows are its species' elements, in order of atomic
// number, and a row of its own for a component without a formula, which
// can't react: air, a mixture, here.
fn test_element_rows() {
	records := ['74-82-8', '7732-18-5', '132259-10-0'].map(db.find(it) or { panic(err) })
	rows := element_rows(records.map(flowsheet.species_of(it)))
	assert rows == [[4.0, 2, 0], [1.0, 0, 0], [0.0, 1, 0], [0.0, 0, 1]], '${rows}'
}

// The conversion reactor at its edges: a full conversion leaves none of its
// key, to rounding, its coefficient 3 here; a step that would take another
// reactant a millionth below zero is refused, past rounding; and a drop
// that leaves no pressure is refused.
fn test_conversion_edges() {
	sp, fl, rs := setup(['7727-37-9', '1333-74-0', '7664-41-7'], ['N2 + 3 H2 = 2 NH3'])
	c := new_conversion_reactor('conv',
		species: sp
		steps:   [Conversion{ reaction: rs[0], key: 1, conversion: 1.0 }]
	)!
	o := c.calculate(fl, [feed(fl, [1.0, 0.7, 0], 700, 1e6)])!
	assert math.abs(o.outlets[0].flows[1]) <= 1e-16, '${o.outlets[0].flows}'
	short := new_conversion_reactor('short',
		species: sp
		steps:   [Conversion{ reaction: rs[0], key: 1, conversion: 1.0 }]
	)!
	mut why := ''
	short.calculate(fl, [feed(fl, [0.1 * (1 - 1e-6), 0.3, 0], 700, 1e6)]) or { why = err.msg() }
	assert why.contains('r1 needs more Nitrogen than there is'), why
	steep := new_conversion_reactor('steep',
		species: sp
		steps:   [Conversion{ reaction: rs[0], key: 1, conversion: 0.5 }]
		drop:    2e6
	)!
	why = ''
	steep.calculate(fl, [feed(fl, [1.0, 3, 0], 700, 1e6)]) or { why = err.msg() }
	assert why.contains('a drop of 2e+06 Pa from 1e+06 Pa leaves no pressure'), why
}

// A flow that rounding takes below zero is read as zero in the rates.
fn test_a_flow_below_zero() {
	k := xu_froment()
	records := k.species.map(db.find(it.cas) or { panic(err) })
	fl := db.fluid(.pr, records.map(it.name), kij: false)!
	mut src := new_source('src', k, [0, 1, 2], 1000, fl.model)!
	st := src.state(900, 1e6, [1.0, 3, 0.3, -1e-20, 0.1])!
	assert st.y[3] == 0 && math.abs(st.y[0] - 1 / 4.4) <= 1e-16
}

// A PFR out of isothermal with a pressure drop, where the Joule–Thomson term
// enters its energy equation, and one with a duty spread over its volume:
// each outlet closes its energy, and the integrated T agrees with its flash
// (which calculate checks).
fn test_pfr_energy_with_drop_and_duty() {
	sp, fl, rs := setup(['74-84-0', '74-85-1', '1333-74-0'], ['C2H6 = C2H4 + H2'])
	k := first_order(sp, rs, '1.5e-6 mol/(m^3*s*Pa)', 'C2H6')
	inlet := feed(fl, [1.0, 0, 0], 1100, 5e5)
	a := new_pfr('a',
		k:             k
		reactions:     [0]
		volume:        0.5
		thermal:       .adiabatic
		drop:          .given
		pressure_drop: 2e5
		rtol:          1e-10
	)!
	oa := a.calculate(fl, [inlet])!
	check_balances(oa, inlet)
	assert oa.outlets[0].p == 3e5
	d := new_pfr('d', k: k, reactions: [0], volume: 0.5, thermal: .duty, duty: 5e4, rtol: 1e-10)!
	od := d.calculate(fl, [inlet])!
	check_balances(od, inlet)
	assert od.outlets[0].t > oa.outlets[0].t
}
