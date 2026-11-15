// Units in expressions (M9 step 2, criterion 2): every node's dimension
// derived and checked, numbers converted to SI as they're read, the place's
// dimension enforced, and typed numbers that look like a constant flagged.
import expr
import math
import quantity

fn context() expr.Context {
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

const no_defs = expr.Definitions{}

fn read(text string) !expr.Expr {
	return expr.parse(context(), &no_defs, text)
}

fn dim_of(text string) quantity.Dim {
	e := read(text) or { panic('${text}: ${err.msg()}') }
	return e.dim
}

fn test_dimensions_follow_the_algebra() {
	per_mass := quantity.unit('mol/(kg*s)')!.dim
	// Xu and Froment's steam reforming rate, with its constants' units: a
	// rate per mass of catalyst.
	k1 := '4.22e15 kmol*bar^0.5/(kg*h) * exp(-240.1 kJ/mol / (R*T))'
	den := '(1 + 6.65e-4 bar^-1 * exp(38.28 kJ/mol / (R*T)) * p[CH4] + 1.77e5 * exp(-88.68 kJ/mol / (R*T)) * p[H2O] / p[H2])'
	r1 := '${k1} / p[H2]^2.5 * (p[CH4] * p[H2O] - p[H2]^3 * p[CO] / (exp(30.42 - 27106 K / T) * 1 bar^2)) / ${den}^2'
	assert dim_of(r1).eq(per_mass), dim_of(r1).str()
	// Every function and operator.
	pressure := quantity.Dimension.pressure.dim()
	assert dim_of('sqrt(P * p[CH4])').eq(pressure)
	assert dim_of('abs(P - p[CH4])').eq(pressure)
	assert dim_of('-P').eq(pressure)
	assert dim_of('P^(1/2) * P^0.5').eq(pressure)
	assert dim_of('(P / p_std)^2.37').is_dimensionless()
	assert dim_of('exp(ln(P / p_std) * 2)').is_dimensionless()
	assert dim_of('R * T / P').eq(quantity.unit('m^3/mol')!.dim)
	assert dim_of('k_B * N_A / R').is_dimensionless()
	assert dim_of('F').eq(quantity.unit('C/mol')!.dim)
	assert dim_of('h * c / (k_B * T)').eq(quantity.unit('m')!.dim)
	assert dim_of('g_n').eq(quantity.unit('m/s^2')!.dim)
	assert dim_of('pi').is_dimensionless()
	// A power that's a constant expression: a fraction from 1/3, or from a
	// sum of numbers.
	assert dim_of('P^(1/3)').eq(pressure.pow(quantity.ratio(1, 3)!)!)
	assert dim_of('P^(2 - 0.5)').eq(pressure.pow(quantity.ratio(3, 2)!)!)
	// A fixed power keeps its sign.
	assert dim_of('P^-0.5').eq(pressure.pow(quantity.ratio(-1, 2)!)!)
}

fn test_every_constant_has_its_unit() {
	for c in expr.constants {
		if c.unit == '' {
			assert c.dim.is_dimensionless(), c.name
			continue
		}
		assert quantity.unit(c.unit)!.dim.eq(c.dim), '${c.name}: ${c.unit}'
		assert quantity.unit(c.unit)!.scale == 1.0, c.name
	}
}

fn test_numbers_are_converted_to_si() {
	// A number with a unit is its SI value, by quantity's scale.
	for text, si in {
		'240.1 kJ/mol':                240100.0
		'6.65e-4 bar^-1':              6.65e-4 * 1e-5
		'1.02e15 kmol*bar^0.5/(kg*h)': 1.02e15 * 1e3 * math.pow(1e5, 0.5) / 3600.0
		'2 mm':                        2e-3
		'0.386 kg/m3':                 0.386
		'3 cP':                        3e-3
		'1 atm':                       101325.0
	} {
		e := read(text)!
		assert math.abs(e.nodes[0].value - si) <= 4.0 * 2.220446049250313e-16 * math.abs(si), '${text}: ${e.nodes[0].value} vs ${si}'
	}
}

// refused asserts that text is refused at a column, with reason
// mentioning want.
fn refused(text string, column int, want string) {
	read(text) or {
		e := err as expr.ExprError
		assert e.reason.contains(want), '${text}: ${err.msg()}'
		assert e.line == 1 && e.column == column, '${text}: column ${e.column}, not ${column}: ${err.msg()}'
		return
	}
	assert false, '"${text}" was read'
}

fn test_a_wrong_dimension_is_refused_with_its_reason() {
	// DWSIM's text: 240100/8.314/T is in 1/K, not a pure number.
	refused('exp(-240100/8.314/T)', 1, 'the argument of exp must be a pure number, but it is in 1/K')
	refused('T + P', 3, 'the two sides of + must have the same units, but the left is a temperature, K and the right a pressure, kg/(m·s²)')
	refused('P - 1', 3, 'the two sides of − must have the same units')
	refused('ln(P)', 1, 'the argument of ln must be a pure number, but it is a pressure')
	refused('log10(2 m)', 1, 'the argument of log10 must be a pure number, but it is a length')
	refused('2^T', 3, 'a power must be a pure number, but this one is a temperature')
	refused('P^(T/(1 K))', 4, 'can only be raised to a fixed power, and this power varies')
	refused('P^(1/7654321)', 4, 'can only be raised to a fraction')
	refused('P^pi', 3, 'can only be raised to a fraction')
	// A sum's mismatch deep inside is found where it is.
	refused('exp(1 + (T + P) / T)', 12, 'the two sides of +')
}

fn test_a_place_asks_for_its_dimension() {
	e := read('4.22e15 kmol/(m3*h)')!
	e.must_be(quantity.unit('mol/(m^3*s)')!.dim, 'a rate per volume')!
	e.must_be(quantity.unit('mol/(kg*s)')!.dim, 'a rate per mass of catalyst') or {
		assert err.msg().contains('a rate per mass of catalyst must be a mol/(kg·s), but this is a mol/(m³·s)'), err.msg()
		return
	}
	assert false, 'a rate per volume passed as one per mass'
}

fn test_numbers_that_look_like_a_constant_are_flagged() {
	// 8.314 for R moves a rate with E/RT ≈ 29 by −0.16%; a transposed 8.341
	// by +9.6% (step 0's research).
	e := read('exp(-240.1 kJ/mol / (8.314 J/(mol*K) * T))')!
	found := e.look_alikes()
	assert found.len == 1 && found[0].constant == 'R'
	assert math.abs(found[0].relative - (8.314 - 8.31446261815324) / 8.31446261815324) <= 1e-15
	assert e.nodes[found[0].node].text == '8.314'
	transposed := read('exp(-240.1 kJ/mol / (8.341 J/(mol*K) * T))')!.look_alikes()
	assert transposed.len == 1 && transposed[0].relative > 0.003
	// R itself, a number far from any constant, a number of another
	// dimension, and 1 atm against the standard pressure, aren't flagged.
	assert read('exp(-240.1 kJ/mol / (R * T))')!.look_alikes().len == 0
	assert read('29.1 J/(mol*K) * T / (1 J/mol)')!.look_alikes().len == 0
	assert read('8.314 * P')!.look_alikes().len == 0
	assert read('1 atm')!.look_alikes().len == 0
	// Nor a pressure within 1% of 1 bar, nor a number that is the constant.
	assert read('1.005 bar')!.look_alikes().len == 0
	assert read('8.31446261815324 J/(mol*K) * T / (1 J/mol)')!.look_alikes().len == 0
	// π, typed short.
	pi := read('3.1416 * P')!.look_alikes()
	assert pi.len == 1 && pi[0].constant == 'pi'
	assert read('3.1 * P')!.look_alikes().len == 0
}
