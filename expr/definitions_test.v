// Definitions (M9 step 2, D2): named expressions, in any order, using one
// another but never themselves, and never taking a name that means
// something already.
import expr
import quantity

fn ctx() expr.Context {
	species := ['CH4', 'H2O', 'H2', 'CO', 'CO2']
	mut keys := map[string]int{}
	for i, s in species {
		keys[s] = i
	}
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, species, keys},
		]
	}
}

// xu_froment is Xu and Froment's constants as definitions, listed out of
// order: the adsorption term before the constants it uses.
fn xu_froment() []expr.Source {
	return [
		expr.Source{'DEN', '1 + K_CO*p[CO] + K_H2*p[H2] + K_CH4*p[CH4] + K_H2O*p[H2O]/p[H2]'},
		expr.Source{'k1', '4.22e15 kmol*bar^0.5/(kg*h) * exp(-240.1 kJ/mol / (R*T))'},
		expr.Source{'K_CH4', '6.65e-4 bar^-1 * exp(38.28 kJ/mol / (R*T))'},
		expr.Source{'K_H2O', '1.77e5 * exp(-88.68 kJ/mol / (R*T))'},
		expr.Source{'K_H2', '6.12e-9 bar^-1 * exp(82.9 kJ/mol / (R*T))'},
		expr.Source{'K_CO', '8.23e-5 bar^-1 * exp(70.65 kJ/mol / (R*T))'},
		expr.Source{'Keq1', 'exp(30.42 - 27106 K / T) * 1 bar^2'},
	]
}

const rate = 'k1 / p[H2]^2.5 * (p[CH4]*p[H2O] - p[H2]^3*p[CO]/Keq1) / DEN^2'

fn test_definitions_are_used_by_name_in_any_order() {
	defs := expr.define(ctx(), xu_froment())!
	// Each comes after what it uses.
	position := fn [defs] (name string) int {
		return defs.order.index(defs.find(name))
	}
	for k in ['K_CO', 'K_H2', 'K_CH4', 'K_H2O'] {
		assert position(k) < position('DEN'), k
	}
	assert defs.exprs[defs.find('DEN')].dim.is_dimensionless()
	e := expr.parse(ctx(), &defs, rate)!
	e.must_be(quantity.unit('mol/(kg*s)')!.dim, 'a rate per mass of catalyst')!
	assert e.uses().map(defs.names[it]) == ['k1', 'Keq1', 'DEN']
	// The rate through its definitions is the rate written out, to the bit.
	// DEN first, since it holds the others, and K_H2O before K_H2, which is
	// in its name.
	mut whole := rate
	sources := xu_froment()
	for name in ['DEN', 'K_H2O', 'K_CH4', 'K_CO', 'K_H2', 'k1', 'Keq1'] {
		s := sources.filter(it.name == name)[0]
		whole = whole.replace(s.name, '(${s.text})')
	}
	inline := expr.parse(ctx(), &expr.Definitions{}, whole)!
	x := [900.0, 2e5, 4e4, 8e4, 2e4, 1e4, 5e3]
	assert e.eval(&defs, x) == inline.eval(&expr.Definitions{}, x)
	// A definition is written back by its name.
	assert e.plain() == rate.replace('*', ' * ').replace('/', ' / ').replace('-', ' - ').replace('  ', ' ')
	assert e.unicode().contains('k₁/p_H₂^2.5') && e.unicode().contains('/Keq₁')
	assert e.latex().starts_with('\\frac{\\frac{k_{1}}{p_{\\mathrm{H_{2}}}^{2.5}}')
	k_ch4 := expr.parse(ctx(), &defs, 'K_CH4 * p[CH4]')!
	assert k_ch4.unicode() == 'K_CH₄·p_CH₄'
	assert k_ch4.latex() == 'K_{\\mathrm{CH_{4}}} \\cdot p_{\\mathrm{CH_{4}}}'
}

// refused asserts that the definitions are refused, in definition def at a
// line and column, with reason mentioning want.
fn refused(sources []expr.Source, def string, line int, column int, want string) expr.ExprError {
	expr.define(ctx(), sources) or {
		assert err is expr.ExprError, err.msg()
		e := err as expr.ExprError
		assert e.definition == def, '${def}: ${err.msg()}'
		assert e.reason.contains(want), '${def}: ${err.msg()}'
		assert e.line == line && e.column == column, '${def}: line ${e.line}, column ${e.column}, not ${line}, ${column}: ${err.msg()}'
		return e
	}
	assert false, '${def} was accepted'
	return expr.ExprError{}
}

fn test_a_cycle_is_refused_with_its_chain() {
	// Reported in the first of the cycle found, at its use of the next.
	refused([expr.Source{'u', '2 * v + 1'}, expr.Source{'v', 'w * 2'}, expr.Source{'w', 'u'}], 'u', 1,
		5, 'u uses itself: u → v → w → u')
	refused([expr.Source{'x', '1 + x'}], 'x', 1, 5, 'x uses itself: x → x')
	// A chain that ends is fine, however long.
	mut chain := []expr.Source{}
	for i in 0 .. 50 {
		chain << expr.Source{'d${i}', if i == 49 { '1' } else { '2 * d${i + 1}' }}
	}
	defs := expr.define(ctx(), chain)!
	assert expr.parse(ctx(), &defs, 'd0')!.eval(&defs, [300.0, 1e5, 0, 0, 0, 0, 0]) == f64(u64(1) << 49)
}

fn test_a_definition_takes_no_name_already_meant() {
	refused([expr.Source{'R', '8.314 J/(mol*K)'}], 'R', 1, 1, "R is a built-in constant, the molar gas constant, N_A·k_B, so it can't name a definition")
	refused([expr.Source{'T', '300 K'}], 'T', 1, 1, 'T is a variable')
	refused([expr.Source{'p', '1 bar'}], 'p', 1, 1, 'p is a family of variables, as in p[CH4]')
	refused([expr.Source{'exp', '1'}], 'exp', 1, 1, 'exp is a function')
	refused([expr.Source{'2k', '1'}], '2k', 1, 1, "can't name a definition: a name is letters, digits and _")
	refused([expr.Source{'k', '1'}, expr.Source{'k', '2'}], 'k', 1, 1, 'k is defined twice')
}

fn test_an_error_in_a_definition_says_which() {
	refused([expr.Source{'k', 'T +'}], 'k', 1, 3, 'a + with nothing after it')
	e := refused([expr.Source{'K_CO', '1 bar^-1'}, expr.Source{'k', 'K_C0 * p[CO]'}], 'k', 1, 1,
		'no name K_C0 here')
	assert 'K_CO' in e.suggestions
	refused([expr.Source{'k', 'T + P'}], 'k', 1, 3, 'the two sides of + must have the same units')
	assert e.msg().starts_with('k: line 1, column 1: no name K_C0 here (did you mean K_CO')
}

fn test_a_fixed_definition_can_be_a_power() {
	defs := expr.define(ctx(), [expr.Source{'n', '2.5'}, expr.Source{'m', 'T / (1 K)'},
		expr.Source{'half', 'n / 5'}])!
	pressure := quantity.Dimension.pressure.dim()
	assert expr.parse(ctx(), &defs, 'P^n')!.dim.eq(pressure.pow(quantity.ratio(5, 2)!)!)
	assert expr.parse(ctx(), &defs, 'P^half')!.dim.eq(pressure.pow(quantity.ratio(1, 2)!)!)
	expr.parse(ctx(), &defs, 'P^m') or {
		assert err.msg().contains('can only be raised to a fixed power, and this power varies'), err.msg()
		return
	}
	assert false, 'a varying definition was a power'
}
