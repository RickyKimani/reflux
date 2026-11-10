// The checks' edge cases (M9 step 6); the papers' and DWSIM's files are
// checked whole in flowsheet/reactions_test.v.
module kinetics

import chem
import components

const check_db = components.load_default() or { panic(err) }

fn gas_species(names ...string) []Species {
	mut out := []Species{}
	for n in names {
		r := check_db.find(n) or { panic(err) }
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

fn findings(sp []Species, eqs []string, laws []Rated, sites []Site, s ?State) []string {
	mut rs := []Reaction{}
	for i, e in eqs {
		rs << parse_equation('r${i + 1}', e, sp) or { panic(err) }
	}
	k := build(sp, rs, sites, laws, []) or { panic(err) }
	return k.check(s).map(it.str())
}

fn test_constants_at_or_below_zero() {
	sp := gas_species('methane', 'water', 'carbon monoxide', 'hydrogen')
	f := findings(sp, ['CH4 + H2O = CO + 3 H2'], [
		Rated{0, PowerLaw{
			k:      Arrhenius{'-2 mol/(kg*s*Pa^2)', '100 kJ/mol', ''}
			orders: [Power{'CH4', 1}, Power{'H2O', 1}]
		}},
	], [], none)
	assert f == ["warning: r1's rate constant is -2 mol/(kg*s*Pa^2): a rate or an adsorption constant is above zero"]
	// A written constant is the user's to keep above zero.
	g := findings(sp, ['CH4 + H2O = CO + 3 H2'], [
		Rated{0, PowerLaw{
			k:      Written{'-2 mol/(kg*s*Pa^2)'}
			orders: [Power{'CH4', 1}, Power{'H2O', 1}]
		}},
	], [], none)
	assert g == []
	// Zero is no more a rate constant than below it is.
	z := findings(sp, ['CH4 + H2O = CO + 3 H2'], [
		Rated{0, PowerLaw{
			k:      Fixed{'0 mol/(kg*s*Pa^2)'}
			orders: [Power{'CH4', 1}, Power{'H2O', 1}]
		}},
	], [], none)
	assert z == ["warning: r1's rate constant is 0 mol/(kg*s*Pa^2): a rate or an adsorption constant is above zero"]
}

// A reaction whose species' ΔG°f disagree with their ΔH°f and S° by over
// 1 kJ/mol is warned of; under it, nothing is said.
fn test_formation_gap() {
	mut sp := gas_species('methane', 'water', 'carbon monoxide', 'hydrogen')
	assert findings(sp, ['CH4 + H2O = CO + 3 H2'], [], [], none) == []
	sp[2] = Species{
		...sp[2]
		dgf: (sp[2].dgf or { 0 }) + 2000
	}
	f := findings(sp, ['CH4 + H2O = CO + 3 H2'], [], [], none)
	assert f.len == 1
	assert f[0].starts_with("warning: r1: its species' ΔG°f sum to 2.0121 kJ/mol from what their ΔH°f and S° give at 298.15 K")
	// A gap the other way is as wrong.
	sp[2] = Species{
		...sp[2]
		dgf: (sp[2].dgf or { 0 }) - 4000
	}
	g := findings(sp, ['CH4 + H2O = CO + 3 H2'], [], [], none)
	assert g.len == 1
	assert g[0].starts_with("warning: r1: its species' ΔG°f sum to -1.9879 kJ/mol")
}

// A fitted K is held to the Gibbs energies' only where the species have the
// data; otherwise that's noted.
fn test_a_fitted_k_without_data() {
	mut sp := gas_species('carbon monoxide', 'water', 'carbon dioxide', 'hydrogen')
	sp[0] = Species{
		...sp[0]
		gas: none
	}
	f := findings(sp, ['CO + H2O = CO2 + H2'], [
		Rated{0, PowerLaw{
			k:           Fixed{'1 mol/(kg*s*Pa^2)'}
			orders:      [Power{'CO', 1}, Power{'H2O', 1}]
			equilibrium: Equilibrium{
				kind:        .correlation
				correlation: Written{'exp(-4.036 + 4400 K / T)'}
			}
		}},
	], [], none)
	assert f == ["note: r1's fitted K can't be held to the Gibbs energies': r1: CO has no ideal-gas data"]
	// With the data, and no state: the range alone, warned of past a factor
	// of 2.
	sp2 := gas_species('carbon monoxide', 'water', 'carbon dioxide', 'hydrogen')
	g := findings(sp2, ['CO + H2O = CO2 + H2'], [
		Rated{0, PowerLaw{
			k:           Fixed{'1 mol/(kg*s*Pa^2)'}
			orders:      [Power{'CO', 1}, Power{'H2O', 1}]
			equilibrium: Equilibrium{
				kind:        .correlation
				correlation: Written{'2 * exp(-4.036 + 4400 K / T)'}
			}
		}},
	], [], none)
	assert g.len == 1
	assert g[0] == "warning: r1's fitted K, K1: from 300 to 1500 K its ln K is -0.12468 to +0.71993 from theirs"
}

// A typed constant in a rate expression itself is the reaction's; and an
// override with a term its site lacks says so.
fn test_look_alikes_and_overrides() {
	sp := [
		Species{
			name:    'A'
			formula: chem.parse('C2H4')!
		},
		Species{
			name:    'B'
			formula: chem.parse('CH2')!
		},
	]
	r := parse_equation('dimer', '2 B = A', sp)!
	k := build(sp, [r], [], [
		Rated{0, Expression{
			text: '1e-9 mol/(kg*s*Pa) * exp(-1000 J/mol / (8.31 J/(mol*K) * T)) * p[B]'
		}},
	], [])!
	f := k.check(State{
		t: 500
		p: 1e5
		y: [0.5, 0.5]
	})
	assert f.len == 1 && f[0].reaction == 0
	assert f[0].message.starts_with("dimer's rate has 8.31 J/(mol*K), 0.00053673 below R")
	assert f[0].message.ends_with("at 500 K it moves dimer's rate by -0.012917% from what R gives")
	site := Site{
		name:  'S'
		terms: [Monomial{'K_A', Fixed{'1 Pa^-1'}, [Power{'A', 1}]}]
	}
	l := Lhhw{
		k:      Fixed{'1 mol/(kg*s*Pa^2)'}
		orders: [Power{'B', 2}]
		site:   'S'
		terms:  [Monomial{'K_A', Fixed{'1 Pa^-1'}, [Power{'A', 1}]},
			Monomial{'K_B', Fixed{'2 Pa^-1'}, [Power{'B', 1}]}]
	}
	g := findings(sp, ['2 B = A'], [Rated{0, l}], [site], none)
	assert g == ["note: r1 has its own adsorption term in place of site S's, with K_B too: a site's term belongs to every reaction on it"]
}
