module fit

// reflux has one R, quantity's (N_A·k_B): every module that names one names
// that double. Until 0.6.0 thermo's and the activity models' were
// 8.314462618, 1.84e-11 below the one expr, kinetics and the fits used, and
// a Gibbs reactor or a PFR mixed the two. fit imports most of the modules
// that name one, so the test is here.
import activity
import expr
import flowsheet
import kinetics
import quantity
import thermo

fn test_every_r_is_quantitys() {
	r := quantity.molar_gas_constant
	assert thermo.r_gas == r
	assert activity.r_gas == r
	assert kinetics.gas_constant == r
	assert gas_r == r
	assert flowsheet.molar_gas_constant() == r
	mut named := 0
	for c in expr.constants {
		if c.name == 'R' {
			assert c.value == r
			named++
		}
	}
	assert named == 1
}
