module quality

import math
import os
import components
import tdep

// The sweep over the built-in data, once for the tests below.
const swept = sweep(components.load_default() or { panic(err) }) or { panic(err) }

fn test_the_committed_report_is_current() {
	// M2 exit criterion 5: failures go into a committed report rather than
	// failing the build, and the report can't silently grow. Any change in
	// what the sweep finds changes the report, and fails here until it's
	// regenerated, so the change shows in a diff.
	path := os.join_path(@VMODROOT, report_path)
	committed := os.read_file(path) or { '' }
	assert committed == generate()!, '${report_path} is out of date. Run `v run tools/data_quality.vsh`.'
}

fn test_what_the_sweep_found() {
	s := swept
	assert s.compounds == 642
	assert s.evaluated == 8932
	assert s.properties.len == 222
	// Nearly all at the critical point itself.
	assert s.near_tc == 211
	assert s.boiling_checked == 641 && s.boiling.len == 132
	assert s.no_ideal_gas.map(it.cas) == ['3967-55-3', '754-12-1']
	assert s.no_component.map(it.cas) == ['3967-55-3']
	assert s.fallbacks.len == 26
	assert s.past_tc.len == 99
	assert s.shared_inchikey.len == 1 && s.shared_inchikey[0].cas == '50-00-0'
	assert s.shared_inchikey[0].what.contains('Carbon monoxide (630-08-0)')
	// R-23 takes its vapour pressure from its second ChemSep copy, which
	// gives 1 atm at Tb; the first copy's gave 0.7% of it.
	assert !s.boiling.any(it.cas == '75-46-7')
}

fn test_pr_against_chemsep() {
	assert swept.pr.len == 20
	by := fn (cas string) PrPsat {
		for p in swept.pr {
			if p.cas == cas {
				return p
			}
		}
		panic(cas)
	}
	// PR with ChemSep's constants follows ChemSep's methane to 1% from
	// Tr 0.65 up, as PR should for a simple fluid.
	for d in by('74-82-8').deviations[1..] {
		assert math.abs(d) < 0.01
	}
	// Methanol at Tr 0.5 is PR-with-Soave's known weak spot.
	assert by('67-56-1').deviations[0] < -0.2
	// n-decane at Tr 0.95 is ChemSep's correlation going wrong, not PR.
	assert by('124-18-5').deviations[3] > 0.5
}

fn test_number_formatting_breaks_ties_the_same_way_everywhere() {
	// Exact ties, which C libraries round differently: half away from zero.
	assert g4(267.25) == '267.3'
	assert g4(-267.25) == '-267.3'
	assert f2(0.125) == '0.13'
	assert f3(1.0005) == '1.001' || f3(1.0005) == '1.000' // 1.0005 isn't exact in binary
	assert g5(101325.0) == '1.0132e+05' || g5(101325.0) == '1.0133e+05'
	assert g4(0.0) == '0'
	assert signed(f2(1.5)) == '+1.50' && signed(f2(-1.5)) == '-1.50'
}

fn prop(key string, eq int, coeffs []f64, tmin f64, tmax f64, tc f64) components.Property {
	return components.Property{
		key:   key
		c:     tdep.new(eq, coeffs, tmin, tmax, tc) or { panic(err) }
		scale: 1.0
	}
}

fn test_the_property_rules() {
	// A line through zero at 250 K, between the sample points at 240 and
	// 260 K, is first negative at 260 K.
	falling := prop('liquid_density', 2, [2.5, -0.01], 100.0, 300.0, 0.0)
	f := check_property('liquid_density', falling, 300.0) or { panic('not caught') }
	assert f.why.starts_with('negative at 260 K (T/Tc 0.867)'), f.why
	assert f.t == 260.0
	// The same line is fine for a second virial coefficient, which may be
	// negative.
	assert check_property('second_virial', falling, 300.0) == none
	// Eq. 106 reaches 0 at Tc: fine for a heat of vaporization, not for a
	// liquid density.
	hv := prop('heat_of_vaporization', 106, [3.0e7, 0.4], 100.0, 400.0, 400.0)
	assert check_property('heat_of_vaporization', hv, 400.0) == none
	zero := check_property('liquid_density', hv, 400.0) or { panic('not caught') }
	assert zero.why == 'zero at 400 K (T/Tc 1.000)'
	// Eq. 105 past its C is NaN.
	rho := prop('liquid_density', 105, [1.0, 0.26, 400.0, 0.28], 100.0, 410.0, 0.0)
	nan := check_property('liquid_density', rho, 0.0) or { panic('not caught') }
	assert nan.why == 'not finite at 410 K'
	// A range starting at 0 K skips T = 0, where forms in 1/T aren't defined.
	inv := prop('vapor_viscosity', 102, [1e-6, 0.5, 100.0, 0.0], 0.0, 1000.0, 0.0)
	assert check_property('vapor_viscosity', inv, 0.0) == none
	// A single-temperature range is checked once, there.
	one := prop('solid_density', 1, [-1.0], 117.97, 117.97, 0.0)
	assert (check_property('solid_density', one, 0.0) or { panic('not caught') }).t == 117.97
}
