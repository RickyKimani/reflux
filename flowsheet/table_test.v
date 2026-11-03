module flowsheet

// The stream table's numbers: 5 significant figures, rounded half away from
// zero whatever the C compiler, a decimal from a thousandth up to a million.

fn test_figures_are_five_significant() {
	for x, want in {
		300.0:         '300.00'
		30.0:          '30.000'
		0.7:           '0.70000'
		1.0:           '1.0000'
		-85.12345:     '-85.123'
		123456.0:      '123460'
		12345.6:       '12346'
		0.001:         '0.0010000'
		0.00099999:    '9.9999e-04'
		1e-17:         '1.0000e-17'
		6.0858898e-12: '6.0859e-12'
		1.5e7:         '1.5000e+07'
		-2.5e-5:       '-2.5000e-05'
		0.0:           '0'
		-0.0:          '0'
	} {
		assert figures(x) == want, '${x}: ${figures(x)}'
	}
}

fn test_a_rounding_that_carries_moves_the_point() {
	// Each rounds up into a new leading digit.
	assert figures(99999.7) == '100000'
	assert figures(999995.0) == '1.0000e+06'
	assert figures(9.99996) == '10.000'
	assert figures(0.000999996) == '0.0010000'
	// 1000 is 3 exactly, but a value an ulp below it is on the far side.
	assert figures(999.9999999999999) == '1000.0'
	assert figures(1000.0) == '1000.0'
}

fn test_a_power_below_one_divides() {
	// 846064.9999999999 · 0.1 rounds to 84606.5, a tie that would round up, but
	// the quotient by 10 is 84606.49999999999: the true value is below the tie.
	assert figures(846064.9999999999) == '846060'
}

fn test_a_tie_rounds_away_from_zero() {
	// 1.03125 and 2.5 are exact doubles, so each is a true tie. glibc's printf
	// rounds a tie to even, giving 1.0312, and the C library under tcc needn't
	// agree (docs/v-for-go.md, §14). math.round goes away from zero under both.
	assert figures(1.03125) == '1.0313'
	assert figures(-1.03125) == '-1.0313'
	assert scientific(2.5, 1) == '3e+00'
	assert scientific(0.125, 2) == '1.3e-01'
}

fn test_scientific_has_its_figures() {
	assert scientific(9.459766303621108e-14, 2) == '9.5e-14'
	assert scientific(1e-9, 2) == '1.0e-09'
	assert scientific(18.667412644357633, 3) == '1.87e+01'
	assert scientific(-0.0410062, 3) == '-4.10e-02'
	assert scientific(1e100, 3) == '1.00e+100'
	assert scientific(0.0, 3) == '0'
}
