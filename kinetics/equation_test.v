module kinetics

import chem

fn sp(name string, formula string, aliases ...string) Species {
	return Species{
		name:    name
		formula: chem.parse(formula) or { panic(err) }
		aliases: aliases
	}
}

// A small fluid: reforming's species, with names as ChemSep gives them.
fn fluid() []Species {
	return [
		Species{
			...sp('Methane', 'CH4')
			cas: '74-82-8'
		},
		sp('Water', 'H2O'),
		sp('Carbon monoxide', 'CO'),
		sp('Hydrogen', 'H2'),
		sp('Carbon dioxide', 'CO2'),
		sp('Oxygen', 'O2'),
		sp('Nitrogen', 'N2'),
		sp('Methanol', 'CH4O'),
		sp('Ethanol', 'C2H6O'),
		sp('Dimethyl ether', 'C2H6O'),
		sp('1-butene', 'C4H8'),
		sp('n-butane', 'C4H10', 'butane'),
		Species{
			name: 'Air'
		},
	]
}

fn read(text string) string {
	r := parse_equation('r', text, fluid()) or { return err.msg() }
	return r.str() + ' | ' + r.species.map(fluid()[it].name).join(', ')
}

fn test_names_formulas_and_cas_numbers() {
	// A term that is a species' name is that species, even where a number
	// touching a name would read too: 13CO is carbon-13 monoxide, not 13 CO.
	labelled := [
		Species{
			name:    '13CO'
			formula: chem.parse('CO')!
		},
		sp('Carbon monoxide', 'CO', 'CO'),
		sp('Oxygen', 'O2'),
		sp('Carbon dioxide', 'CO2'),
	]
	r13 := parse_equation('r', '13CO + 0.5 O2 = CO2', labelled) or { panic(err) }
	assert r13.species == [0, 2, 3]
	assert r13.nu == [-1.0, -0.5, 1]
	two := parse_equation('r', '2CO2 = 2 Carbon monoxide + O2', labelled)!
	assert two.nu == [-2.0, 2, 1]
	assert read('CH4 + H2O = CO + 3 H2') == 'CH4 + H2O = CO + 3 H2 | Methane, Water, Carbon monoxide, Hydrogen'
	assert read('Methane + water = carbon monoxide + 3 Hydrogen') == 'Methane + water = carbon monoxide + 3 Hydrogen | Methane, Water, Carbon monoxide, Hydrogen'
	assert read('74-82-8 + HOH = CO + 3H2') == '74-82-8 + HOH = CO + 3 H2 | Methane, Water, Carbon monoxide, Hydrogen'
	assert read('CO2 + 3 H2 = CH3OH + H2O') == 'CO2 + 3 H2 = CH3OH + H2O | Carbon dioxide, Hydrogen, Methanol, Water'
	assert read('CO + 0.5 O2 = CO2') == 'CO + 0.5 O2 = CO2 | Carbon monoxide, Oxygen, Carbon dioxide'
	assert read('CO + 1/2 O2 = CO2') == 'CO + 0.5 O2 = CO2 | Carbon monoxide, Oxygen, Carbon dioxide'
	assert read('1-butene + H2 = butane') == '1-butene + H2 = butane | 1-butene, Hydrogen, n-butane'
	assert read("'1-butene' + H2 = 'n-butane'") == '1-butene + H2 = n-butane | 1-butene, Hydrogen, n-butane'
	assert read('2 CO + O2 = 2 CO2') == '2 CO + O2 = 2 CO2 | Carbon monoxide, Oxygen, Carbon dioxide'
	r := parse_equation('shift', 'CO + H2O = CO2 + H2', fluid())!
	assert r.nu == [-1.0, -1, 1, 1]
	assert r.species == [2, 1, 4, 3]
}

fn test_unbalanced_equations_say_how() {
	assert read('CH4 + H2O = CO + H2') == "r: CH4 + H2O = CO + H2 doesn't balance: H: 6 on the left, 2 on the right; balanced, it would be CH4 + H2O = CO + 3 H2"
	assert read('CH4 + O2 = CO2 + H2O') == "r: CH4 + O2 = CO2 + H2O doesn't balance: H: 4 on the left, 2 on the right; O: 2 on the left, 3 on the right; balanced, it would be CH4 + 2 O2 = CO2 + 2 H2O"
	assert read('CO + O2 = CO2') == "r: CO + O2 = CO2 doesn't balance: O: 3 on the left, 2 on the right; balanced, it would be 2 CO + O2 = 2 CO2"
	assert read('CH3OH = CO + H2') == "r: CH3OH = CO + H2 doesn't balance: H: 4 on the left, 2 on the right; balanced, it would be CH3OH = CO + 2 H2"
	// A species the elements leave no part for is left out, and one on the
	// wrong side is moved, each said.
	assert read('CH4 + O2 + N2 = CO2 + H2O') == "r: CH4 + O2 + N2 = CO2 + H2O doesn't balance: H: 4 on the left, 2 on the right; N: 2 on the left, 0 on the right; O: 2 on the left, 3 on the right; balanced, it would be CH4 + 2 O2 = CO2 + 2 H2O, leaving out N2"
	assert read('CO2 + H2 = CO + O2') == "r: CO2 + H2 = CO + O2 doesn't balance: H: 2 on the left, 0 on the right; O: 2 on the left, 3 on the right; balanced, it would be 2 CO2 = 2 CO + O2, leaving out H2"
	assert read('CO + H2O + CO2 = H2') == "r: CO + H2O + CO2 = H2 doesn't balance: C: 2 on the left, 0 on the right; O: 4 on the left, 0 on the right; balanced, it would be CO + H2O = CO2 + H2, with CO2 on the other side"
	// No coefficients, or many.
	assert read('CH4 = CO2 + H2O') == "r: CH4 = CO2 + H2O doesn't balance: H: 4 on the left, 2 on the right; O: 0 on the left, 3 on the right; no coefficients would balance it, as O is only on the right"
	assert read('CO = CO2') == "r: CO = CO2 doesn't balance: O: 1 on the left, 2 on the right; no coefficients would balance these species"
	assert read('CH4 + H2O = CO + CO2 + H2') == "r: CH4 + H2O = CO + CO2 + H2 doesn't balance: C: 1 on the left, 2 on the right; H: 6 on the left, 2 on the right; O: 1 on the left, 3 on the right; these species balance in 2 independent ways, so there's no one equation to suggest"
	// Fractions balance exactly, and are written as decimals where they are.
	assert read('CO + 0.4 O2 = CO2') == "r: CO + 0.4 O2 = CO2 doesn't balance: O: 1.8 on the left, 2 on the right; balanced, it would be 2 CO + O2 = 2 CO2"
	assert read('CO + 1/3 O2 = CO2') == "r: CO + 1/3 O2 = CO2 doesn't balance: O: 5/3 on the left, 2 on the right; balanced, it would be 2 CO + O2 = 2 CO2"
	assert read('1/3 CO + 1/6 O2 = 1/3 CO2') == '1/3 CO + 1/6 O2 = 1/3 CO2 | Carbon monoxide, Oxygen, Carbon dioxide'
}

fn test_refused() {
	assert read('CH4 + H2O -> CO + 3 H2') == 'r: column 11: an equation is written with =, not ->: whether a reaction runs both ways is for its rate law to say'
	assert read('CH4 + H2O ⇌ CO + 3 H2') == 'r: column 11: an equation is written with =, not ⇌: whether a reaction runs both ways is for its rate law to say'
	assert read('CH4 + H2O') == 'r: no = between the two sides'
	assert read('CH4 = CO = H2') == 'r: column 10: a second ='
	assert read('CH4 + = CO') == 'r: column 7: a term with nothing in it'
	assert read('CH4 + CH4 = C2H6') == 'r: column 7: CH4 appears twice on the left'
	assert read('CO + H2O + H2 = CO2 + 2 H2') == 'r: column 23: H2 is on both sides: write its net coefficient on one'
	assert read('C2H6O = C2H4 + H2O') == 'r: column 1: C2H6O is the formula of both Ethanol and Dimethyl ether: name one'
	assert read('CH4 + H2O = CO + 3 H2 + Air') == "r: Air has no formula, as a mixture's pseudo-compound, so it can't react"
	assert read('Methan + H2O = CO + 3 H2') == 'r: column 1: no species Methan (did you mean Methane, Methanol?)'
	assert read('C3H8 = C3H6 + H2') == 'r: column 1: no species here is C3H8 (did you mean C4H8, CH4?)'
	assert read('0 CO + O2 = CO2') == 'r: column 1: a coefficient of 0: leave the species out'
	assert read('1/0 CO + O2 = CO2') == 'r: column 1: a coefficient, 1/0, over zero'
	assert read('1.2.3 CO = CO') == "r: column 1: a coefficient, 1.2.3, that isn't a number: write 2, 0.5 or 1/2"
	assert read('3 = CO') == 'r: column 1: a coefficient, 3, with no species after it'
	assert read("'CO + O2 = CO2") == 'r: column 1: a quote never closed'
	assert read('12345678 CO = CO') == 'r: column 1: a coefficient, 12345678, out of range'
}

// Rng is a small seeded generator for the fuzzing: splitmix64.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) below(n int) int {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return int((z ^ (z >> 31)) % u64(n))
}

// 20,000 seeded equations of reforming's pieces and stray ones: each is
// refused with a ReactionError, with a column inside it where it has one,
// or reads to a balanced reaction that its own writing reads back to.
fn test_random_equations_never_panic() {
	pieces := ['CH4', 'H2O', 'CO', 'H2', 'CO2', 'O2', 'Methane', 'water', '74-82-8', '3H2', '2',
		'0.5', '1/2', '1/0', '0', ' + ', '+', ' = ', '=', '->', '⇌', "'", "'n-butane'", 'Air',
		'C2H6O', '1-butene', ' ', 'x', '(', '12345678']
	sp := fluid()
	mut rng := Rng{20261112}
	mut accepted := 0
	balanced := ['CH4 + H2O = CO + 3 H2', 'CO + H2O = CO2 + H2', 'Methane + 2 water = CO2 + 4 H2',
		'CO + 1/2 O2 = CO2', "'1-butene' + H2 = 'n-butane'"]
	for _ in 0 .. 20000 {
		// Half of them random pieces; half a balanced equation, as it is or
		// with a piece put in somewhere.
		mut t := ''
		if rng.below(2) == 0 {
			for _ in 0 .. 1 + rng.below(9) {
				t += pieces[rng.below(pieces.len)]
			}
		} else {
			t = balanced[rng.below(balanced.len)]
			if rng.below(3) > 0 {
				at := rng.below(t.len + 1)
				t = t[..at] + pieces[rng.below(pieces.len)] + t[at..]
			}
		}
		r := parse_equation('r', t, sp) or {
			assert err is ReactionError, '"${t}": ${err.msg()}'
			e := err as ReactionError
			assert e.column >= 0 && e.column <= t.runes().len + 1 && e.reason != '', '"${t}": ${err.msg()}'
			continue
		}
		back := parse_equation('r', r.str(), sp) or { panic('"${t}" wrote "${r}": ${err.msg()}') }
		assert back.species == r.species && back.exact == r.exact, '"${t}" wrote "${r}"'
		accepted++
	}
	assert accepted > 3000
}
