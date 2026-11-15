// Rate laws and sites (M9 step 5, D4, D6): Xu and Froment's and Vanden
// Bussche and Froment's kinetics from structured forms, against mpmath from
// the papers' own formulas (testdata/m9_kinetics.json); each form rendered
// as its paper prints it, and the same bits as the rate typed; the chain
// rule through K(T) and the state; and everything that's refused.
module kinetics

import chem
import components
import expr
import json2
import math
import os

fn reforming_species() []Species {
	mut out := []Species{}
	for name, f in {
		'Methane':         'CH4'
		'Water':           'H2O'
		'Hydrogen':        'H2'
		'Carbon monoxide': 'CO'
		'Carbon dioxide':  'CO2'
	} {
		out << Species{
			name:    name
			formula: chem.parse(f) or { panic(err) }
		}
	}
	return out
}

fn reforming_reactions(sp []Species) []Reaction {
	return [
		parse_equation('reforming', 'CH4 + H2O = CO + 3 H2', sp) or { panic(err) },
		parse_equation('shift', 'CO + H2O = CO2 + H2', sp) or { panic(err) },
		parse_equation('overall', 'CH4 + 2 H2O = CO2 + 4 H2', sp) or { panic(err) },
	]
}

// Xu and Froment's site: one adsorption term, every reaction's.
fn xu_froment_site() Site {
	return Site{
		name:   'Ni'
		symbol: 'DEN'
		terms:  [
			Monomial{'K_CO', VantHoff{'8.23e-5 bar^-1', '-70.65 kJ/mol', ''}, [
				Power{'CO', 1},
			]},
			Monomial{'K_H2', VantHoff{'6.12e-9 bar^-1', '-82.90 kJ/mol', ''}, [
				Power{'H2', 1},
			]},
			Monomial{'K_CH4', VantHoff{'6.65e-4 bar^-1', '-38.28 kJ/mol', ''}, [
				Power{'CH4', 1},
			]},
			Monomial{'K_H2O', VantHoff{'1.77e5', '88.68 kJ/mol', ''}, [
				Power{'H2O', 1},
				Power{'H2', -1},
			]},
		]
	}
}

// Xu and Froment's laws: A from their Table 6 and E from Table 5. The K's
// are the correlations later papers use with these kinetics, their paper
// naming K1, K2 and K3 with their units but giving no values.
fn xu_froment_laws(keq Reversibility) []Rated {
	corr := [
		Written{'exp(30.114 - 26830 K / T) * 1 bar^2'},
		Written{'exp(-4.036 + 4400 K / T)'},
		Written{'exp(30.114 - 26830 K / T) * exp(-4.036 + 4400 K / T) * 1 bar^2'},
	]
	ks := [
		Arrhenius{'4.225e15 kmol*bar^0.5/(kg*h)', '240.1 kJ/mol', ''},
		Arrhenius{'1.955e6 kmol/(kg*h*bar)', '67.13 kJ/mol', ''},
		Arrhenius{'1.020e15 kmol*bar^0.5/(kg*h)', '243.9 kJ/mol', ''},
	]
	orders := [
		[Power{'CH4', 1}, Power{'H2O', 1}, Power{'H2', -2.5}],
		[Power{'CO', 1}, Power{'H2O', 1}, Power{'H2', -1}],
		[Power{'CH4', 1}, Power{'H2O', 2}, Power{'H2', -3.5}],
	]
	mut out := []Rated{}
	for i in 0 .. 3 {
		out << Rated{i, Lhhw{
			k:           ks[i]
			orders:      orders[i]
			equilibrium: Equilibrium{
				kind:        keq
				correlation: corr[i]
			}
			form:        .difference
			site:        'Ni'
			power:       2
		}}
	}
	return out
}

fn xu_froment() Kinetics {
	sp := reforming_species()
	return build(sp, reforming_reactions(sp), [xu_froment_site()], xu_froment_laws(.correlation),
		[]) or { panic(err) }
}

fn methanol_species() []Species {
	mut out := []Species{}
	for name, f in {
		'Carbon dioxide':  'CO2'
		'Hydrogen':        'H2'
		'Methanol':        'CH4O'
		'Water':           'H2O'
		'Carbon monoxide': 'CO'
	} {
		out << Species{
			name:    name
			formula: chem.parse(f) or { panic(err) }
		}
	}
	return out
}

// Vanden Bussche and Froment's methanol kinetics: one site, its term cubed
// for methanol and to the first power for the reverse shift, with Graaf's
// equilibrium constants.
fn vanden_bussche_froment() Kinetics {
	sp := methanol_species()
	rs := [
		parse_equation('methanol', 'CO2 + 3 H2 = CH3OH + H2O', sp) or { panic(err) },
		parse_equation('RWGS', 'CO2 + H2 = CO + H2O', sp) or { panic(err) },
	]
	site := Site{
		name:  'Cu'
		terms: [
			Monomial{'K_a', Fixed{'3453.38'}, [Power{'H2O', 1}, Power{'H2', -1}]},
			Monomial{'K_b', VantHoff{'0.499 bar^-0.5', '-17197 J/mol', ''}, [
				Power{'H2', 0.5},
			]},
			Monomial{'K_c', VantHoff{'6.62e-11 bar^-1', '-124119 J/mol', ''}, [
				Power{'H2O', 1},
			]},
		]
	}
	laws := [
		Rated{0, Lhhw{
			k:           Arrhenius{'1.07 mol/(kg*s*bar^2)', '-36696 J/mol', ''}
			orders:      [Power{'CO2', 1}, Power{'H2', 1}]
			equilibrium: Equilibrium{
				kind:        .correlation
				correlation: Written{'10^(3066 K / T - 10.592) * 1 bar^-2'}
			}
			site:        'Cu'
			power:       3
		}},
		Rated{1, Lhhw{
			k:           Arrhenius{'1.22e10 mol/(kg*s*bar)', '94765 J/mol', ''}
			orders:      [Power{'CO2', 1}]
			equilibrium: Equilibrium{
				kind:        .correlation
				correlation: Written{'10^(2.029 - 2073 K / T)'}
			}
			site:        'Cu'
		}},
	]
	return build(sp, rs, [site], laws, []) or { panic(err) }
}

struct KCase {
	t      f64
	p      []f64
	rates  []f64
	drdt   []f64
	drdp   [][]f64
	scales []f64
}

struct KModel {
	name    string
	species []string
	cases   []KCase
}

struct KFixture {
	models []KModel
}

fn kinetics_fixture() !KFixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm9_kinetics.json')
	return json2.decode[KFixture](os.read_file(path)!)!
}

// at puts a fixture's state into a point: T, P, and the partial pressures,
// in the fixture's species order, which is the kinetics' own.
fn at(k &Kinetics, c KCase) Point {
	mut pt := k.point()
	pt.x[0] = c.t
	mut total := 0.0
	for i, p in c.p {
		pt.x[k.ctx.slot(0, i)] = p
		total += p
	}
	pt.x[1] = total
	return pt
}

// Each model's rates are mpmath's from the paper's own formulas, to 1e-13
// of each rate's largest term (the forward one where the backward nearly
// cancels it: the tape rounds a few ulps of each; 9.4e-15 at worst); and
// their derivatives with respect to T and every partial pressure are
// mpmath's to 1e-12 of the larger of the derivative and the scale over T or
// over the pressure (1.1e-14 at worst).
fn check_model(k &Kinetics, m KModel) {
	mut buf := k.buffer()
	mut out := []f64{len: k.rates.len}
	mut g := []f64{len: k.rates.len * k.tape.inputs}
	mut worst := 0.0
	mut worst_d := 0.0
	for c in m.cases {
		pt := at(k, c)
		k.rates(&pt, mut buf, mut out)
		k.tape.gradient(mut buf, mut g)
		for i, r in out {
			e := math.abs(r - c.rates[i]) / c.scales[i]
			worst = math.max(worst, e)
			assert e <= 1e-13, '${m.name} rate ${i} at ${c.t} K: ${r} vs ${c.rates[i]}'
			dt := g[i * k.tape.inputs]
			et := math.abs(dt - c.drdt[i]) / math.max(math.abs(c.drdt[i]), c.scales[i] / c.t)
			worst_d = math.max(worst_d, et)
			assert et <= 1e-12, '${m.name} dr${i}/dT at ${c.t} K: ${dt} vs ${c.drdt[i]}'
			for j in 0 .. c.p.len {
				dp := g[i * k.tape.inputs + k.ctx.slot(0, j)]
				want := c.drdp[j][i]
				ep := math.abs(dp - want) / math.max(math.abs(want), c.scales[i] / c.p[j])
				worst_d = math.max(worst_d, ep)
				assert ep <= 1e-12, '${m.name} dr${i}/dp${j} at ${c.t} K: ${dp} vs ${want}'
			}
		}
	}
	println('${m.name}: ${m.cases.len} states, rates within ${worst} of their scales, derivatives ${worst_d}')
}

fn test_xu_froment_against_mpmath() {
	f := kinetics_fixture()!
	k := xu_froment()
	assert f.models[0].species == ['CH4', 'H2O', 'H2', 'CO', 'CO2']
	check_model(&k, f.models[0])
}

fn test_vanden_bussche_froment_against_mpmath() {
	f := kinetics_fixture()!
	k := vanden_bussche_froment()
	assert f.models[1].species == ['CO2', 'H2', 'CH3OH', 'H2O', 'CO']
	check_model(&k, f.models[1])
}

// Each rate renders as its paper prints it: Xu and Froment's as a
// difference over DEN², Vanden Bussche and Froment's as an approach to
// equilibrium over θ³, with every constant a definition.
fn test_renders_as_published() {
	xf := xu_froment()
	assert xf.rates.map(it.unicode()) == [
		'k₁/p_H₂^2.5·(p_CH₄·p_H₂O − p_CO·p_H₂³/K₁)/DEN²',
		'k₂/p_H₂·(p_CO·p_H₂O − p_CO₂·p_H₂/K₂)/DEN²',
		'k₃/p_H₂^3.5·(p_CH₄·p_H₂O² − p_CO₂·p_H₂⁴/K₃)/DEN²',
	]
	den := xf.defs.exprs[xf.defs.find('DEN')]
	assert den.plain() == '1 + K_CO * p[CO] + K_H2 * p[H2] + K_CH4 * p[CH4] + K_H2O * p[H2O] / p[H2]'
	kco := xf.defs.exprs[xf.defs.find('K_CO')]
	assert kco.plain() == '8.23e-5 bar^-1 * exp(70.65 kJ/mol / (R * T))'
	vbf := vanden_bussche_froment()
	assert vbf.rates.map(it.plain()) == [
		'k1 * p[CO2] * p[H2] * (1 - p[CH3OH] * p[H2O] / (p[CO2] * p[H2]^3) / K1) / theta_Cu^3',
		'k2 * p[CO2] * (1 - p[CO] * p[H2O] / (p[CO2] * p[H2]) / K2) / theta_Cu',
	]
	theta := vbf.defs.exprs[vbf.defs.find('theta_Cu')]
	assert theta.plain() == '1 + K_a * p[H2O] / p[H2] + K_b * p[H2]^0.5 + K_c * p[H2O]'
}

// A structured rate and the same rate typed give the same bits (D4): Xu and
// Froment's three, typed as a user would, as step 3 did, over the same
// definitions.
fn test_structured_is_typed() {
	xf := xu_froment()
	sp := reforming_species()
	rs := reforming_reactions(sp)
	mut own := []expr.Source{}
	for i, name in xf.defs.names {
		own << expr.Source{name, xf.defs.exprs[i].text}
	}
	typed := [
		'k1 / p[H2]^2.5 * (p[CH4]*p[H2O] - p[H2]^3*p[CO]/K1) / DEN^2',
		'k2 / p[H2] * (p[CO]*p[H2O] - p[H2]*p[CO2]/K2) / DEN^2',
		'k3 / p[H2]^3.5 * (p[CH4]*p[H2O]^2 - p[H2]^4*p[CO2]/K3) / DEN^2',
	]
	mut laws := []Rated{}
	for i, t in typed {
		laws << Rated{i, Expression{
			text: t
		}}
	}
	k := build(sp, rs, [], laws, own)!
	f := kinetics_fixture()!
	mut b1 := xf.buffer()
	mut b2 := k.buffer()
	mut o1 := []f64{len: 3}
	mut o2 := []f64{len: 3}
	for c in f.models[0].cases {
		p1 := at(&xf, c)
		p2 := at(&k, c)
		xf.rates(&p1, mut b1, mut o1)
		k.rates(&p2, mut b2, mut o2)
		for i in 0 .. 3 {
			assert math.f64_bits(o1[i]) == math.f64_bits(o2[i]), 'rate ${i} at ${c.t} K'
		}
	}
}

const db = components.load_default() or { panic(err) }

fn species_from(cas []string) []Species {
	mut out := []Species{}
	for c in cas {
		r := db.find(c) or { panic(err) }
		out << Species{
			name:    r.name
			cas:     r.ids.cas
			formula: r.formula() or { panic(err) }
			gas:     r.ideal_gas() or { panic(err) }
			dgf:     r.constants.dgf
		}
	}
	return out
}

fn with_t(s State, t f64) State {
	return State{
		...s
		t: t
	}
}

fn with_p(s State, p f64) State {
	return State{
		...s
		p: p
	}
}

fn with_y(s State, j int, dy f64) State {
	mut y := s.y.clone()
	y[j] += dy
	return State{
		...s
		y: y
	}
}

fn rates_of(k &Kinetics, s State) []f64 {
	mut pt := k.point()
	k.fill(s, mut pt) or { panic(err) }
	mut buf := k.buffer()
	mut out := []f64{len: k.rates.len}
	k.rates(&pt, mut buf, mut out)
	return out
}

// With K from the Gibbs energies (D7), the rates use Keq[...] at 1 bar
// times p_std^Δn, and the jacobian with respect to T, P and the mole
// fractions, through every family and each K's d ln K/dT, is central
// differences', to 1e-6 of each row's largest term in the variable's units.
fn test_thermo_k_and_the_jacobian() {
	sp := species_from(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	rs := reforming_reactions(sp)
	k := build(sp, rs, [xu_froment_site()], xu_froment_laws(.thermo), [])!
	assert k.rates[0].plain().contains('/ (Keq[reforming] * p_std^2)')
	assert k.rates[1].plain().contains('/ Keq[shift]')
	assert k.k_used == [0, 1, 2]
	s := State{
		t: 900
		p: 15e5
		y: [0.2, 0.45, 0.2, 0.1, 0.05]
	}
	mut pt := k.point()
	k.fill(s, mut pt)!
	mut buf := k.buffer()
	mut out := []f64{len: 3}
	k.rates(&pt, mut buf, mut out)
	mut g := []f64{len: 3 * k.tape.inputs}
	w := 2 + sp.len
	mut d := []f64{len: 3 * w}
	k.jacobian(s, &pt, mut buf, mut g, mut d)
	units := [s.t, s.p, 1.0, 1.0, 1.0, 1.0, 1.0]
	for v in 0 .. w {
		h := [1e-3, 1.0, 1e-7][math.min(v, 2)]
		up := match v {
			0 { with_t(s, s.t + h) }
			1 { with_p(s, s.p + h) }
			else { with_y(s, v - 2, h) }
		}
		down := match v {
			0 { with_t(s, s.t - h) }
			1 { with_p(s, s.p - h) }
			else { with_y(s, v - 2, -h) }
		}
		ru := rates_of(&k, up)
		rd := rates_of(&k, down)
		for i in 0 .. 3 {
			fd := (ru[i] - rd[i]) / (2 * h)
			mut row := 0.0
			for q in 0 .. w {
				row = math.max(row, math.abs(d[i * w + q]) * units[q])
			}
			assert math.abs(d[i * w + v] - fd) <= 1e-6 * row / units[v], 'rate ${i}, variable ${v}: ${d[i * w + v]} vs ${fd}'
		}
	}
}

// Rates are evaluated into the caller's buffer, allocating nothing.
fn test_rates_allocate_nothing() {
	k := xu_froment()
	f := kinetics_fixture()!
	mut pt := at(&k, f.models[0].cases[0])
	mut buf := k.buffer()
	mut out := []f64{len: 3}
	before := gc_heap_usage().total_bytes
	mut sum := 0.0
	for i in 0 .. 100000 {
		pt.x[0] = 800.0 + f64(i % 300)
		k.rates(&pt, mut buf, mut out)
		sum += out[0]
	}
	assert gc_heap_usage().total_bytes - before == 0
	assert math.is_finite(sum)
}

fn single(sp []Species, eq string, law RateLaw, sites []Site) !Kinetics {
	r := parse_equation('r', eq, sp)!
	return build(sp, [r], sites, [Rated{0, law}], [])
}

// The other forms and bases, against the same rates worked by hand.
fn test_power_law_eley_rideal_and_bases() {
	sp := reforming_species()
	s := State{
		t: 800
		p: 2e5
		y: [0.25, 0.5, 0.1, 0.1, 0.05]
	}
	// A power law per volume on concentrations: r = k·c_CH4·c_H2O².
	pl := single(sp, 'CH4 + H2O = CO + 3 H2', PowerLaw{
		k:      Fixed{'2e-3 m^6/(mol^2*s)'}
		orders: [Power{'CH4', 1}, Power{'H2O', 2}]
		basis:  .concentration
		per:    .volume
	}, [])!
	assert pl.rates[0].plain() == 'k1 * c[CH4] * c[H2O]^2'
	c := s.p / (gas_constant * s.t)
	want := 2e-3 * (0.25 * c) * math.pow(0.5 * c, 2)
	assert math.abs(rates_of(&pl, s)[0] - want) <= 1e-15 * want
	// On mole fractions, with a centred Arrhenius constant: r = k(T)·y_CH4.
	fr := single(sp, 'CH4 + H2O = CO + 3 H2', PowerLaw{
		k:      Arrhenius{'0.5 mol/(kg*s)', '80 kJ/mol', '700 K'}
		orders: [Power{'CH4', 1}]
		basis:  .fraction
	}, [])!
	assert fr.rates[0].plain() == 'k1 * y[CH4]'
	kt := 0.5 * math.exp(-80e3 / gas_constant * (1 / 800.0 - 1 / 700.0))
	assert math.abs(rates_of(&fr, s)[0] - kt * 0.25) <= 1e-14 * kt
	// Eley–Rideal: CO adsorbed, H2O from the gas, irreversible.
	site := Site{
		name:  'S'
		terms: [
			Monomial{'', Fixed{'3 bar^-1'}, [Power{'CO', 1}]},
			Monomial{'', Fixed{'0.5 bar^-1'}, [Power{'H2', 1}]},
		]
	}
	er := single(sp, 'CO + H2O = CO2 + H2', EleyRideal{
		k:        Fixed{'1e-6 mol/(kg*s*Pa)'}
		adsorbed: 'CO'
		gas:      'H2O'
		site:     'S'
	}, [site])!
	assert er.rates[0].plain() == 'k1 * K_CO * p[CO] * p[H2O] / theta_S'
	pco, ph2o, ph2 := 0.1 * s.p, 0.5 * s.p, 0.1 * s.p
	ka, kb := 3e-5, 0.5e-5
	want_er := 1e-6 * ka * pco * ph2o / (1 + ka * pco + kb * ph2)
	assert math.abs(rates_of(&er, s)[0] - want_er) <= 1e-14 * want_er
	// A reversible power law on fugacities, from the Gibbs energies: zero
	// where Q is K, whatever the fugacity coefficients.
	gs := species_from(['630-08-0', '7732-18-5', '124-38-9', '1333-74-0'])
	shift := parse_equation('shift', 'CO + H2O = CO2 + H2', gs)!
	kk := build(gs, [shift], [], [
		Rated{0, PowerLaw{
			k:           Fixed{'1e-9 mol/(kg*s*Pa^2)'}
			orders:      [Power{'CO', 1}, Power{'H2O', 1}]
			equilibrium: Equilibrium{
				kind: .thermo
			}
			basis:       .fugacity
		}},
	], [])!
	assert kk.rates[0].plain() == 'k1 * f[CO] * f[H2O] * (1 - f[CO2] * f[H2] / (f[CO] * f[H2O]) / Keq[shift])'
	kt2 := math.exp(shift.standard(gs, 1000)!.ln_k)
	phi := [1.01, 0.97, 1.02, 1.05]
	// y_CO2·y_H2·φ_CO2·φ_H2 = K·y_CO·y_H2O·φ_CO·φ_H2O, with y_CO = y_H2O = 0.2.
	y2 := math.sqrt(kt2 * 0.2 * 0.2 * phi[0] * phi[1] / (phi[2] * phi[3]))
	at_eq := State{
		t:   1000
		p:   3e6
		y:   [0.2, 0.2, y2, y2]
		phi: phi
	}
	r_eq := rates_of(&kk, at_eq)[0]
	forward := 1e-9 * (0.2 * 3e6 * phi[0]) * (0.2 * 3e6 * phi[1])
	assert math.abs(r_eq) <= 1e-14 * forward, '${r_eq}'
}

// A reaction may override its site's term, as the owner's DWSIM case's
// three denominators do; the override is its own definition, and shares
// the site's constants where their texts agree.
fn test_an_override() {
	sp := reforming_species()
	rs := reforming_reactions(sp)
	site := xu_froment_site()
	mut laws := xu_froment_laws(.correlation)
	l1 := laws[1].law as Lhhw
	// The shift's denominator without the CH4 term, as in DWSIM's file.
	laws[1] = Rated{1, Lhhw{
		...l1
		terms: [site.terms[0], site.terms[1], site.terms[3]]
	}}
	k := build(sp, rs, [site], laws, [])!
	assert k.rates[1].plain().ends_with('/ DEN_2^2')
	own := k.defs.exprs[k.defs.find('DEN_2')]
	assert own.plain() == '1 + K_CO * p[CO] + K_H2 * p[H2] + K_H2O * p[H2O] / p[H2]'
	// The same name for another constant is refused.
	mut clash := site.terms.clone()
	clash[0] = Monomial{'K_CO', Fixed{'1 bar^-1'}, [Power{'CO', 1}]}
	laws[1] = Rated{1, Lhhw{
		...l1
		terms: clash
	}}
	if _ := build(sp, rs, [site], laws, []) {
		assert false
	} else {
		assert err.msg().contains('K_CO already names a term of site Ni')
	}
}

fn refused(sp []Species, eq string, law RateLaw, sites []Site) string {
	k := single(sp, eq, law, sites) or { return err.msg() }
	return 'built ${k.rates[0].plain()}'
}

fn test_refused() {
	sp := reforming_species()
	eq := 'CH4 + H2O = CO + 3 H2'
	k := Fixed{'1 mol/(kg*s*bar^2)'}
	// A rate per volume whose constant gives one per mass.
	assert refused(sp, eq, PowerLaw{
		k:      k
		orders: [Power{'CH4', 1}, Power{'H2O', 1}]
		per:    .volume
	}, []).contains('a rate per volume must be a')
	assert refused(sp, eq, PowerLaw{
		k:      k
		orders: [Power{'CH5', 1}]
	}, []).contains('no species here is CH5')
	assert refused(sp, eq, Lhhw{
		k:      k
		orders: [Power{'CH4', 1}]
		site:   'Pt'
	}, []) == 'r: no site Pt'
	// A fitted K in the wrong units: reforming's is in bar², not bar.
	assert refused(sp, eq, PowerLaw{
		k:           k
		orders:      [Power{'CH4', 1}, Power{'H2O', 1}]
		equilibrium: Equilibrium{
			kind:        .correlation
			correlation: Written{'3 bar'}
		}
	}, []).contains("its K, in its rate's activities, must be a")
	assert refused(sp, eq, EleyRideal{
		k:        k
		adsorbed: 'H2O'
		gas:      'CH4'
		site:     'S'
	}, [Site{'S', '', [Monomial{'', Fixed{'1 bar^-1'}, [Power{'CO', 1}]}]}]) == 'r: H2O has no term of its own on site S, K·a¹, for an Eley–Rideal rate'
	assert refused(sp, eq, Expression{
		text: 'k9 * p[CH4]'
	}, []).contains('k9')
	r := parse_equation('r', eq, sp)!
	if _ := build(sp, [r], [], [Rated{0, PowerLaw{
		k:      k
		orders: []
	}}, Rated{0, Expression{
		text: '1 mol/(kg*s)'
	}}], []) {
		assert false
	} else {
		assert err.msg() == 'r: two rate laws'
	}
	if _ := build(sp, [r], [], [Rated{3, Expression{
		text: '1 mol/(kg*s)'
	}}], []) {
		assert false
	} else {
		assert err.msg() == 'a rate law: no reaction 4'
	}
	if _ := build(sp, [r], [], [Rated{0, PowerLaw{
		k:      k
		orders: []
	}}], [expr.Source{'k1', '2'}]) {
		assert false
	} else {
		assert err.msg() == "r's rate constant: k1 already names a definition of your own: give this one a name of its own"
	}
}

// reforming_gas is reforming's species, with their gas data, and four
// reactions: reforming, the shift, the overall reaction and methanation.
fn reforming_gas() ([]Species, []Reaction) {
	sp := species_from(['74-82-8', '7732-18-5', '1333-74-0', '630-08-0', '124-38-9'])
	mut rs := reforming_reactions(sp)
	rs << parse_equation('methanation', 'CO + 3 H2 = CH4 + H2O', sp) or { panic(err) }
	return sp, rs
}

// One rate on each family, K from the Gibbs energies in each: the
// jacobian's chain rule through p, f (with φ ≠ 1), c (with Z ≠ 1), y and
// each K, against central differences, to 1e-6 of each row's largest term.
fn test_the_jacobian_on_every_family() {
	sp, rs := reforming_gas()
	thermo_k := Equilibrium{
		kind: .thermo
	}
	laws := [
		Rated{0, PowerLaw{
			k:           Fixed{'1e-12 mol/(kg*s*Pa^2)'}
			orders:      [Power{'CH4', 1}, Power{'H2O', 1}]
			equilibrium: thermo_k
			basis:       .pressure
		}},
		Rated{1, PowerLaw{
			k:           Fixed{'1e-11 mol/(kg*s*Pa^2)'}
			orders:      [Power{'CO', 1}, Power{'H2O', 1}]
			equilibrium: thermo_k
			basis:       .fugacity
		}},
		Rated{2, PowerLaw{
			k:           Fixed{'1e-3 m^3/(mol*s)'}
			orders:      [Power{'CH4', 1}, Power{'H2O', 1}]
			equilibrium: thermo_k
			basis:       .concentration
			per:         .volume
		}},
		Rated{3, PowerLaw{
			k:           Fixed{'2 mol/(kg*s)'}
			orders:      [Power{'CO', 1}, Power{'H2', 2}]
			equilibrium: thermo_k
			basis:       .fraction
		}},
	]
	k := build(sp, rs, [], laws, [])!
	assert k.k_used == [0, 1, 2, 3]
	s := State{
		t:   850
		p:   12e5
		y:   [0.15, 0.4, 0.3, 0.1, 0.05]
		phi: [0.98, 0.95, 1.03, 1.01, 0.97]
		z:   1.02
	}
	mut pt := k.point()
	k.fill(s, mut pt)!
	mut buf := k.buffer()
	mut out := []f64{len: 4}
	k.rates(&pt, mut buf, mut out)
	mut g := []f64{len: 4 * k.tape.inputs}
	w := 2 + sp.len
	mut d := []f64{len: 4 * w}
	k.jacobian(s, &pt, mut buf, mut g, mut d)
	units := [s.t, s.p, 1.0, 1.0, 1.0, 1.0, 1.0]
	for v in 0 .. w {
		h := [1e-3, 1.0, 1e-7][math.min(v, 2)]
		up := match v {
			0 { with_t(s, s.t + h) }
			1 { with_p(s, s.p + h) }
			else { with_y(s, v - 2, h) }
		}
		down := match v {
			0 { with_t(s, s.t - h) }
			1 { with_p(s, s.p - h) }
			else { with_y(s, v - 2, -h) }
		}
		ru := rates_of(&k, up)
		rd := rates_of(&k, down)
		for i in 0 .. 4 {
			fd := (ru[i] - rd[i]) / (2 * h)
			mut row := 0.0
			for q in 0 .. w {
				row = math.max(row, math.abs(d[i * w + q]) * units[q])
			}
			assert math.abs(d[i * w + v] - fd) <= 1e-6 * row / units[v], 'rate ${i}, variable ${v}: ${d[i * w + v]} vs ${fd}'
		}
	}
}

// At reforming's equilibrium, its rate is zero on every basis: K's standard
// state becomes each basis's own, p°^Δn, (p°/RT)^Δn or (p°/P)^Δn, with Δn
// = 2.
fn test_equilibrium_zeroes_every_basis() {
	sp, rs := reforming_gas()
	t, p := 900.0, 2e5
	kt := math.exp(rs[0].standard(sp, t)!.ln_k)
	// y_CO·y_H2³/(y_CH4·y_H2O)·(P/p°)² = K.
	y_ch4, y_h2o, y_co := 0.1, 0.3, 0.05
	y_h2 := math.cbrt(kt * math.pow(1e5 / p, 2) * y_ch4 * y_h2o / y_co)
	s := State{
		t: t
		p: p
		y: [y_ch4, y_h2o, y_h2, y_co, 0.0]
	}
	for basis in [Basis.pressure, .concentration, .fraction, .fugacity] {
		per := if basis == .concentration { Per.volume } else { Per.mass }
		unit := match basis {
			.pressure, .fugacity { 'mol/(kg*s*Pa^2)' }
			.concentration { 'm^3/(mol*s)' }
			.fraction { 'mol/(kg*s)' }
		}
		k := build(sp, [rs[0]], [], [
			Rated{0, PowerLaw{
				k:           Fixed{'1 ${unit}'}
				orders:      [Power{'CH4', 1}, Power{'H2O', 1}]
				equilibrium: Equilibrium{
					kind: .thermo
				}
				basis:       basis
				per:         per
			}},
		], [])!
		mut pt := k.point()
		k.fill(s, mut pt)!
		f := match basis {
			.pressure, .fugacity { y_ch4 * p * y_h2o * p }
			.concentration { y_ch4 * y_h2o * math.pow(p / (gas_constant * t), 2) }
			.fraction { y_ch4 * y_h2o }
		}
		r := rates_of(&k, s)[0]
		assert math.abs(r) <= 1e-13 * f, '${basis}: ${r} against a forward rate of ${f}'
	}
}
