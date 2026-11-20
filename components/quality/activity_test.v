// M6 step 1, criterion 2: each of ChemSep's UNIFAC tables against thermo's
// copy of the published one, every difference settled (docs/data-quality.md),
// and the rule that tells rounding from a difference.
module quality

import activity
import components

fn test_the_unifac_tables_against_thermo() {
	db := components.load_default()!
	th := read_thermo_tables()!
	mut got := map[string][]int{}
	for v, pair in {
		activity.Variant.original: [th.original, ThermoTable{}]
		.lle:                      [th.lle, ThermoTable{}]
		.dortmund:                 [th.dortmund, th.dortmund_2006]
	} {
		c := compare(v, components.load_unifac(db.vendor, v)!.table, pair[0], pair[1])
		got[v.name()] = [c.subgroups_same, c.subgroups_rounded, c.subgroups_diff.len, c.subgroups_own,
			c.params_same, c.params_rounded, c.params_older, c.params_diff.len, c.params_own]
		// Every difference left is settled: thermo's typo, or DDBST's page
		// printing ChemSep's value.
		for k in c.params_keys {
			assert v == .dortmund && k in ddbst_dortmund, '${v.name()}: ${k}'
		}
		for s in c.subgroups_diff {
			assert v == .lle && s.starts_with('30 FCH2O'), s
		}
	}
	println(got)
	assert got == {
		'UNIFAC':          [103, 7, 0, 23, 1049, 107, 0, 0, 208]
		'UNIFAC-LLE':      [56, 0, 1, 2, 512, 0, 0, 0, 4]
		'UNIFAC-Dortmund': [124, 0, 0, 3, 1042, 425, 3, 12, 46]
	}
}

fn test_rounding_is_told_from_a_difference() {
	assert sig_digits(287.7) == 4 && sig_digits(727.0) == 3 && sig_digits(-2.2499e-05) == 5
	assert sig_digits(0.0135) == 3
	assert rounds_to(8484.0, 8483.5) // rounded
	assert rounds_to(287.7, 287.77) // cut short
	assert rounds_to(0.0135, 0.013472)
	assert rounds_to(-0.4858, -0.4857) // the last digit
	assert !rounds_to(4228.0, 4288.0) // a digit wrong
	assert !rounds_to(54.47, 54.57)
	assert !rounds_to(-27.545, -27.45) // a digit too many
	// A c too small to move ψ: −7.4e-14 against −7.4e-16 changes the
	// exponent by 3e-11 at 373 K.
	assert all_round([717.485, -1.4852, -7.4379e-14], [717.485, -1.4852, -7.44e-16])
	assert !all_round([2622.04, -11.7013, 0.013472], [2622.04, -11.7013, 0.0145])
}
