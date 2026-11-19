// Writing an expression back (M9 step 2, criterion 1): plain text that
// reads back to the same tree, Unicode for a terminal and LaTeX for a
// notebook, each to the exact string, so a change to one is seen.
// sympy_test.v holds every rendering of 258 cases to SymPy's meaning.
import expr
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
			expr.Variable{'x', quantity.Dim{}},
		]
		families:  [
			expr.Family{'p', pressure, species, keys},
		]
	}
}

// renders asserts an expression's three renderings.
fn renders(text string, plain string, unicode string, latex string) {
	defs := expr.define(context(), [expr.Source{'theta', 'x'}, expr.Source{'k1', 'x'},
		expr.Source{'K_CO', '1 bar^-1'}, expr.Source{'k12', 'x'}]) or { panic(err.msg()) }
	e := expr.parse(context(), &defs, text) or { panic('${text}: ${err.msg()}') }
	assert e.plain() == plain, '${text}: plain ${e.plain()}'
	assert e.unicode() == unicode, '${text}: unicode ${e.unicode()}'
	assert e.latex() == latex, '${text}: latex ${e.latex()}'
}

fn test_operators_and_parentheses() {
	renders('x+x*x', 'x + x * x', 'x + x·x', 'x + x \\cdot x')
	renders('(x+x)*x', '(x + x) * x', '(x + x)·x', '\\left(x + x\\right) \\cdot x')
	renders('x-(x-x)', 'x - (x - x)', 'x − (x − x)', 'x - \\left(x - x\\right)')
	renders('x+(x+x)', 'x + (x + x)', 'x + (x + x)', 'x + \\left(x + x\\right)')
	renders('x*(x/x)', 'x * (x / x)', 'x·(x/x)', 'x \\cdot \\frac{x}{x}')
	renders('x/(x*x)', 'x / (x * x)', 'x/(x·x)', '\\frac{x}{x \\cdot x}')
	renders('-x^2', '-x^2', '−x²', '-x^{2}')
	renders('(-x)^2', '(-x)^2', '(−x)²', '\\left(-x\\right)^{2}')
	renders('x^-2', 'x^(-2)', 'x⁻²', 'x^{-2}')
	renders('x^2^3', 'x^2^3', 'x^2³', 'x^{2^{3}}')
	renders('(x^2)^3', '(x^2)^3', '(x²)³', '\\left(x^{2}\\right)^{3}')
	renders('x^0.5', 'x^0.5', 'x^0.5', 'x^{0.5}')
	renders('x + -x', 'x + (-x)', 'x + (−x)', 'x + \\left(-x\\right)')
	renders('- -x', '-(-x)', '−(−x)', '-\\left(-x\\right)')
	renders('-(x+x)', '-(x + x)', '−(x + x)', '-\\left(x + x\\right)')
}

fn test_functions() {
	renders('exp(-x)', 'exp(-x)', 'exp(−x)', '\\exp\\left(-x\\right)')
	renders('ln(x)+log10(x)', 'ln(x) + log10(x)', 'ln(x) + log₁₀(x)', '\\ln\\left(x\\right) + \\log_{10}\\left(x\\right)')
	renders('sqrt(x)', 'sqrt(x)', '√x', '\\sqrt{x}')
	renders('sqrt(x+1)', 'sqrt(x + 1)', '√(x + 1)', '\\sqrt{x + 1}')
	renders('abs(x-1)', 'abs(x - 1)', '|x − 1|', '\\left|x - 1\\right|')
	renders('exp(x)^2', 'exp(x)^2', 'exp(x)²', '\\left(\\exp\\left(x\\right)\\right)^{2}')
}

fn test_names_members_and_constants() {
	renders('p[CH4]/p[H2O]', 'p[CH4] / p[H2O]', 'p_CH₄/p_H₂O', '\\frac{p_{\\mathrm{CH_{4}}}}{p_{\\mathrm{H_{2}O}}}')
	renders('theta^2', 'theta^2', 'θ²', '\\theta^{2}')
	renders('k12', 'k12', 'k₁₂', 'k_{12}')
	renders('k1*K_CO*P', 'k1 * K_CO * P', 'k₁·K_CO·P', 'k_{1} \\cdot K_{\\mathrm{CO}} \\cdot P')
	renders('N_A*k_B/R', 'N_A * k_B / R', 'N_A·k_B/R', '\\frac{N_{A} \\cdot k_{B}}{R}')
	renders('pi*P/p_std', 'pi * P / p_std', 'π·P/p°', '\\frac{\\pi \\cdot P}{p^{\\circ}}')
	renders('g_n*h/c', 'g_n * h / c', 'gₙ·h/c', '\\frac{g_{n} \\cdot h}{c}')
}

fn test_numbers_and_units() {
	renders('4.22e15*x', '4.22e15 * x', '4.22·10¹⁵·x', '4.22 \\times 10^{15} \\cdot x')
	renders('1e-3*x', '1e-3 * x', '10⁻³·x', '10^{-3} \\cdot x')
	renders('x/2.5e-1', 'x / 2.5e-1', 'x/(2.5·10⁻¹)', '\\frac{x}{2.5 \\times 10^{-1}}')
	renders('(1e3)^2', '1e3^2', '(10³)²', '\\left(10^{3}\\right)^{2}')
	renders('6.65e-4 bar^-1 * P', '6.65e-4 bar^-1 * P', '(6.65·10⁻⁴ bar⁻¹)·P', '6.65 \\times 10^{-4}\\,\\mathrm{bar^{-1}} \\cdot P')
	renders('240.1 kJ/mol / (R*T)', '240.1 kJ/mol / (R * T)', '(240.1 kJ/mol)/(R·T)', '\\frac{240.1\\,\\mathrm{kJ/mol}}{R \\cdot T}')
	renders('2 [K] / T', '2 [K] / T', '(2 K)/T', '\\frac{2\\,\\mathrm{K}}{T}')
	renders('(1 bar)^2 / P^2', '(1 bar)^2 / P^2', '(1 bar)²/P²', '\\frac{\\left(1\\,\\mathrm{bar}\\right)^{2}}{P^{2}}')
	renders('1 kmol*bar^0.5/(kg*h) * x', '1 kmol*bar^0.5/(kg*h) * x', '(1 kmol·bar^0.5/(kg·h))·x', '1\\,\\mathrm{kmol\\cdot bar^{0.5}/(kg\\cdot h)} \\cdot x')
	renders('2 m3/s * x', '2 m3/s * x', '(2 m³/s)·x', '2\\,\\mathrm{m^{3}/s} \\cdot x')
}

fn test_xu_and_froment_render_as_published() {
	// The steam reforming rate, written out: k₁/p_H₂^2.5·(p_CH₄·p_H₂O −
	// p_H₂³·p_CO/K₁)/DEN², as Xu and Froment (1989) print it.
	defs := expr.define(context(), [
		expr.Source{'k1', '4.22e15 kmol*bar^0.5/(kg*h) * exp(-240.1 kJ/mol / (R*T))'},
		expr.Source{'K1', 'exp(30.42 - 27106 K / T) * 1 bar^2'},
		expr.Source{'DEN', '1 + 6.65e-4 bar^-1 * p[CH4]'},
	])!
	e := expr.parse(context(), &defs, 'k1 / p[H2]^2.5 * (p[CH4]*p[H2O] - p[H2]^3*p[CO]/K1) / DEN^2')!
	assert e.unicode() == 'k₁/p_H₂^2.5·(p_CH₄·p_H₂O − p_H₂³·p_CO/K₁)/DEN²'
	assert e.latex() == '\\frac{\\frac{k_{1}}{p_{\\mathrm{H_{2}}}^{2.5}} \\cdot \\left(p_{\\mathrm{CH_{4}}} \\cdot p_{\\mathrm{H_{2}O}} - \\frac{p_{\\mathrm{H_{2}}}^{3} \\cdot p_{\\mathrm{CO}}}{K_{1}}\\right)}{DEN^{2}}'
	k1 := defs.exprs[defs.find('k1')]
	// A unit never runs into the operator after it.
	assert k1.unicode() == '(4.22·10¹⁵ kmol·bar^0.5/(kg·h))·exp(−(240.1 kJ/mol)/(R·T))'
}
