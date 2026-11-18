// The intervals' coverage (M8 step 11, criterion 9): over many seeded
// synthetic data sets from known constants, how often each 95% interval
// holds the true value, which must be about as often as it claims.
module fit

import flowsheet
import math
import rand

const coverage_file = 'format = 3

[components]
names = ["ethane", "hydrogen", "methane"]

[reaction.hydrogenolysis]
equation = "C2H6 + H2 = 2 CH4"
law = "power_law"
k = { A = "1e4 mol/(kg*s*bar^1.5)", Ea = "75000 J/mol" }
orders = { C2H6 = 1, H2 = 0.5 }
K = "none"

[fit.rates]
data = "made in the test"
parameters = ["k1.A", "k1.Ea"]
weights = { "r[hydrogenolysis]" = "5 %" }
'

// The true constants, in SI: A = 2e4 mol/(kg·s·bar^1.5) = 2e4·(1e-5)^1.5
// mol/(kg·s·Pa^1.5), and Eₐ.
const true_a = 2e4 * math.pow(1e-5, 1.5)
const true_ea = 80000.0

// normal is a standard normal deviate, by Box and Muller.
fn normal() f64 {
	u := 1.0 - rand.f64()
	v := rand.f64()
	return math.sqrt(-2.0 * math.log(u)) * math.cos(2.0 * math.pi * v)
}

// data_set is 20 rates at five temperatures from 500 to 600 K, at random
// compositions and pressures, each with 5% Gaussian noise.
fn data_set() Data {
	mut rows := [][]f64{}
	for i in 0 .. 20 {
		t := 500.0 + 25.0 * f64(i % 5)
		p := 1e5 * (1.0 + 4.0 * rand.f64())
		raw := [0.1 + rand.f64(), 0.1 + rand.f64(), 0.1 + rand.f64()]
		total := raw[0] + raw[1] + raw[2]
		y := raw.map(it / total)
		r := true_a * math.exp(-true_ea / (gas_r * t)) * y[0] * p * math.sqrt(y[1] * p)
		rows << [t, p, y[0], y[1], y[2], r * (1.0 + 0.05 * normal())]
	}
	return Data{
		file:    'made in the test'
		columns: [Column{'T', 'K'}, Column{'P', 'Pa'}, Column{'y[ethane]', ''},
			Column{'y[hydrogen]', ''}, Column{'y[methane]', ''},
			Column{'r[hydrogenolysis]', 'mol/(kg*s)'}]
		rows:    rows
	}
}

// 400 data sets: each interval should hold its true value in 380 of them,
// binomially ± 4.4; the test allows ± 3.3 of those, 99.9% two-sided, 366 to
// 394. A, reported back from the centred fit by the delta method, is the
// less linear of the two.
fn test_the_intervals_cover_as_they_claim() {
	file := flowsheet.load_fits_text(coverage_file)!
	plan := file.fits[0]
	rand.seed([u32(2026), 1118])
	mut covered := [0, 0]
	truth := [true_a, true_ea]
	for _ in 0 .. 400 {
		f := fit_rates(file, plan, data_set())
		assert f.ok, f.error
		for q in 0 .. 2 {
			if f.lower[q] <= truth[q] && truth[q] <= f.upper[q] {
				covered[q]++
			}
		}
	}
	println('coverage of 400: A ${covered[0]}, Ea ${covered[1]}')
	for q in 0 .. 2 {
		assert covered[q] >= 366 && covered[q] <= 394, '${plan.parameters[q].name}: ${covered[q]} of 400'
	}
}
