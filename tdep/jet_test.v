// A correlation's jet (M6 step 5): its value eval's, its first derivative
// deriv's closed form, and its second the first's central difference, for
// every form a γ–φ liquid takes: ChemSep's own coefficients for methanol
// and water (water's density is a form 106), and a form 10 and a form 16.
import math
import numeric
import tdep

struct Case {
	what   string
	eq     int
	coeffs []f64
	tmin   f64
	tmax   f64
	tc     f64
}

const cases = [
	Case{'methanol, vapour pressure', 101, [73.4034, -6548.08, -7.40999, 5.72492e-06, 2.0], 175.47, 512.64, 0},
	Case{'water, vapour pressure', 101, [74.555, -7295.59, -7.44245, 4.2881e-06, 2.0], 263.15, 647.29, 0},
	Case{'methanol, liquid density', 105, [1.7918, 0.23929, 512.64, 0.21078], 175.47, 503.1, 0},
	Case{'water, liquid density', 106, [32.5162, -3.213, 7.92411, -7.3599, 2.70352], 253.1, 647.0, 647.14},
	Case{'an Antoine form', 10, [21.0, 3600.0, -40.0], 250.0, 450.0, 0},
	Case{'a ChemSep form 16', 16, [0.5, -4100.0, 22.0, -0.002, 1e-6], 250.0, 450.0, 0},
]

fn test_each_form_against_its_value_and_derivatives() {
	mut worst_v, mut worst_d, mut worst_dd := 0.0, 0.0, 0.0
	for c in cases {
		cor := tdep.new(c.eq, c.coeffs, c.tmin, c.tmax, c.tc)!
		for k in 0 .. 21 {
			t := c.tmin + (c.tmax - c.tmin) * (0.02 + 0.96 * f64(k) / 20.0)
			j := cor.jet(numeric.along(t))!
			v := cor.eval(t)
			d := cor.deriv(t)
			h := 1e-4 * t
			// Richardson's extrapolation of two central differences of the
			// closed-form first derivative.
			d1 := (cor.deriv(t + h) - cor.deriv(t - h)) / (2.0 * h)
			d2 := (cor.deriv(t + h / 2.0) - cor.deriv(t - h / 2.0)) / h
			dd := (4.0 * d2 - d1) / 3.0
			at := '${c.what} at ${t} K'
			worst_v = math.max(worst_v, math.abs(j.v - v) / math.abs(v))
			worst_d = math.max(worst_d, math.abs(j.d - d) / math.max(math.abs(d), 1e-300))
			worst_dd = math.max(worst_dd, math.abs(j.dd - dd) / math.max(math.abs(dd), math.abs(d) / t))
			assert math.abs(j.v - v) <= 1e-14 * math.abs(v), '${at}: ${j.v} vs ${v}'
		}
	}
	println('value ${worst_v:.1e}, first derivative ${worst_d:.1e}, second ${worst_dd:.1e}')
	// Measured: 3.3e-16, 9.3e-15 and 7.9e-10. The value and first
	// derivative are the closed forms' to a few ulps of terms that grow
	// near Tc (water's form 106 at 647 K); the second's reference
	// is a difference, good to its h⁴ truncation and ε/h² rounding.
	assert worst_d <= 5e-14 && worst_dd <= 1e-8
}

// Form 106 is 0 from Tc on, with derivatives that don't move; a form
// without a jet says so.
fn test_the_edges() {
	f106 := tdep.new(106, [5.2, 0.3, 0.1, -0.05, 0.02], 250.0, 600.0, 512.64)!
	j := f106.jet(numeric.along(520.0))!
	assert j.v == 0.0 && j.d == 0.0 && j.dd == 0.0
	f100 := tdep.new(100, [1.0, 2.0], 250.0, 600.0, 0)!
	if _ := f100.jet(numeric.along(300.0)) {
		assert false, 'a jet of form 100'
	} else {
		assert err.msg().contains('equation 100 has no jet'), err.msg()
	}
}
