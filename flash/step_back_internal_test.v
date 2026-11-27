// Internal test (module flash): StepBack, which decides whether a search's
// step that failed left the model's range (M5a step 2). IAPWS-95 has no state
// below 251.165 K; a cubic has no Branches, so it never steps back.
module flash

import math
import thermo
import water

fn test_step_back_asks_the_model() {
	m := water.iapws95()
	z := [1.0]
	// A search in T, and one in ln T, each just either side of 251.165 K.
	in_t := StepBack{
		m: m
		z: z
	}
	in_ln_t := StepBack{
		m:    m
		z:    z
		ln_t: true
	}
	assert !in_t.outside(251.2)
	assert in_t.outside(251.1)
	assert !in_ln_t.outside(math.log(251.2))
	assert in_ln_t.outside(math.log(251.1))
	// Above Tc there's one branch and no limit; the zero value never steps
	// back, as TV's search in ln P uses it.
	assert !in_t.outside(700.0)
	assert !StepBack{}.outside(100.0)
	// A cubic has no Branches: never.
	c := thermo.new_cubic(.pr, [
		thermo.Component{
			name:  'water'
			tc:    647.14
			pc:    22.064e6
			omega: 0.344
			mw:    0.0180153
		},
	])!
	assert !StepBack{
		m: c
		z: z
	}.outside(1.0)
}
