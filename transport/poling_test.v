// The book's worked examples (M5b's criterion 3, and criterion 2's rules):
// Poling, Prausnitz and O'Connell, The Properties of Gases and Liquids, 5th
// ed. (McGraw-Hill, 2001), each example's inputs as it prints them, and
// each answer to its printed figures. The book stays out of the
// repository; the numbers here are its examples', with their pages.
import activity
import math
import tdep
import transport

// printed is half a unit in the last figure of a number printed with
// `figures` significant figures: how far a value may be from it and still
// print so.
fn printed(x f64, figures int) f64 {
	e := int(math.floor(math.log10(math.abs(x))))
	return 0.5 * math.pow10(e - figures + 1)
}

// near asserts got prints as want does, to its figures, with `slack` more
// where the book carries rounded intermediates (said at each use).
fn near(got f64, want f64, figures int, slack f64, at string) {
	d := math.abs(got - want)
	println('${at}: ${got} against ${want} (${d / printed(want, figures)} of half its last figure)')
	assert d <= printed(want, figures) + slack, '${at}: ${got} against ${want}'
}

// constant is a correlation that's the value v at every T, in units of
// `scale` SI: the book's pure values, given at the example's T alone.
fn constant(v f64, scale f64) transport.Correlated {
	return transport.Correlated{
		key:   'given'
		c:     tdep.Correlation{
			eq:   100
			a:    v
			tmin: 1
			tmax: 10000
		}
		scale: scale
	}
}

// Example 9-10 (p. 9.37): ammonia at 420 K and 300 bar by Lucas's method,
// D5's factor on a pure gas above its Tc, where taking Y against its value
// at zero pressure is Lucas's own Z₂/Z₁. Lucas's mixture rule gives a pure
// gas's Pc as R·Tc·Zc/Vc, so Vc is taken so that it's the book's Pc. The
// book prints η° = 147 µP (Z₁ = 0.7258), Z₂ = 3.466, Y = 4.775, Fp = 0.860
// and η = 602 µP; Y·Fp is 4.1065 from its rounded figures, which carry
// 0.0026 between them. It works its coefficients at Tr = 1.036 and Pr =
// 2.643, rounded (420/405.50 is 1.03576, 300/113.53 is 2.6425), and its c
// and d go as Tr^−37.7 and Tr^−7.6: so Y·Fp is held at those, and at the
// state itself (4.1096) only through η, which the book gives on ξ printed
// as 4.95e-3 and as 4.96e-3, 0.2% apart.
fn test_example_9_10_lucas_dense_gas() {
	nh3 := transport.Pure{
		name:   'ammonia'
		tc:     405.50
		pc:     113.53e5
		zc:     0.244
		vc:     8.31446261815324 * 405.50 * 0.244 / 113.53e5
		mw:     0.017031
		dipole: 1.47
	}
	m := transport.lucas_mixture([nh3], [1.0], 420.0)!
	near(m.fp0, 1.164, 4, 0, 'Fp°')
	near(m.pc, 113.53e5, 5, 1e-6, 'Pc')
	eta0 := transport.lucas_dilute(nh3, 420.0)!
	near(eta0 * 1e7, 147, 3, 0, 'η°, µP')
	at_book := transport.lucas_factor([nh3], [1.0], 1.036 * 405.50, 2.643 * 113.53e5)!
	near(at_book, 4.775 * 0.860, 4, 0.0026, "Y·Fp at the book's Tr and Pr")
	f := transport.lucas_factor([nh3], [1.0], 420.0, 300e5)!
	near(eta0 * f * 1e7, 602, 3, 0.002 * 602, 'η, µP')
}

// Example 9-7 (p. 9.24): ammonia with hydrogen, 67.7% ammonia, at 33 °C
// (306.2 K) and low pressure, Lucas's mixture rules (D5's): Tcm = 285.2 K,
// Pcm = 89.6 bar, Mm = 12.18, Fp°m = 1.107, Fq°m = 1.051 (A = 1, since
// M_H/M_L = 8.4 < 9), and η = 116.1 µP. The book's Fq°(H₂) takes Tr =
// 9.209 where its own table has 9.223, which prints the same 1.158; its
// Fp°m is truncated: its own (1.159)(0.677) + 0.323 is 1.10764; and its η
// is Z₁ rounded to 0.750 over ξ, where Z₁ is 0.7510 (0.15%).
fn test_example_9_7_lucas_mixture_rules() {
	nh3 := transport.Pure{
		name:   'ammonia'
		tc:     405.50
		pc:     113.5e5
		vc:     72.5e-6
		zc:     0.244
		mw:     0.017031
		dipole: 1.47
	}
	h2 := transport.Pure{
		name:    'hydrogen'
		tc:      33.2
		pc:      13.0e5
		vc:      64.3e-6
		zc:      0.306
		mw:      0.002016
		dipole:  0
		quantum: 0.76
	}
	m := transport.lucas_mixture([nh3, h2], [0.677, 0.323], 306.2)!
	near(m.tc, 285.2, 4, 0, 'Tcm')
	near(m.pc * 1e-5, 89.6, 3, 0, 'Pcm, bar')
	near(m.mw * 1e3, 12.18, 4, 0, 'Mm')
	assert m.fp0 >= 1.107 && m.fp0 < 1.108, 'Fp°m: ${m.fp0}, not 1.107 truncated'
	near(m.fq0, 1.051, 4, 0, 'Fq°m')
	near(m.dilute(306.2) * 1e7, 116.1, 4, 0.0015 * 116.1, 'η, µP')
}

// Example 9-5 (p. 9.20): methane with n-butane, 30.3% n-butane, at 293 K,
// by Wilke's rule (D4) on the book's pure viscosities: φ₁₂ = 2.268, φ₂₁ =
// 0.416, η = 92.26 µP. The book's η is its rounded φ₂₁'s (the rule on
// 2.268 and 0.416 gives 92.257); on φ₂₁ = 0.41625 unrounded it's 92.245.
fn test_example_9_5_wilke() {
	w := transport.new_wilke([constant(109.4, 1e-7), constant(72.74, 1e-7)], [
		0.016043,
		0.058123,
	])!
	mut k := w.room()
	w.at(293.0, mut k)!
	near(k.a[1] * k.a[1] / k.root[1], 2.268, 4, 0, 'φ12')
	near(k.a[2] * k.a[2] / k.root[2], 0.416, 3, 0, 'φ21')
	rounded := 0.697 * 109.4 / (0.697 + 0.303 * 2.268) + 0.303 * 72.74 / (0.303 + 0.697 * 0.416)
	near(rounded, 92.26, 4, 0, "η on the book's rounded φ, µP")
	near(w.mixed(k, [0.697, 0.303]) * 1e7, 92.245, 5, 0, 'η, µP')
}

// Example 10-5 (p. 10.32): benzene with argon, 25% benzene, at 373.8 K and
// 1 bar, by Wassiljewa's form with Mason and Saxena's Aᵢⱼ = φᵢⱼ (D4) on the
// book's pure viscosities and conductivities: A₁₂ = 0.459, A₂₁ = 2.630,
// λ = 1.84e-2 W/(m·K).
fn test_example_10_5_mason_saxena() {
	w := transport.new_wilke([constant(92.5, 1e-7), constant(271.0, 1e-7)], [0.078114, 0.039948])!
	mut k := w.room()
	w.at(373.8, mut k)!
	near(k.a[1] * k.a[1] / k.root[1], 0.459, 3, 0, 'A12')
	near(k.a[2] * k.a[2] / k.root[2], 2.630, 4, 0, 'A21')
	near(w.conductivity(k, [0.25, 0.75], [1.66e-2, 2.14e-2]), 1.84e-2, 3, 0, 'λ')
}

// Example 10-3 (p. 10.22): nitrous oxide at 105 °C and 138 bar by Stiel and
// Thodos (D5), at the book's V = 144 cm³/mol (its Z = 0.63, rounded):
// ρr = 0.674, Γ = 208, λ − λ° = 1.81e-2 W/(m·K), and λ = 4.15e-2 on its
// λ° = 2.34e-2.
fn test_example_10_3_stiel_thodos() {
	n2o := transport.Pseudo{
		tc: 309.6
		pc: 72.55e5
		vc: 97.0e-6
		zc: 0.273
		mw: 0.044013
	}
	r := transport.stiel_thodos(n2o, 144e-6)!
	near(r, 1.81e-2, 3, 0, 'λ − λ°')
	near(2.34e-2 + r, 4.15e-2, 3, 0, 'λ')
}

// Example 10-6 (p. 10.37): methane with carbon dioxide, 75.5% methane, at
// 370.8 K and 174.8 bar, Stiel and Thodos on Yorizane's rules (D5) at the
// measured V = 159 cm³/mol. The book prints Tcm = 215.4 K, Zcm = 0.286,
// Mm = 22.9 and Γ = 176, and has three slips:
// - it prints Vcm = 97.9 cm³/mol, where its own sum, and every later step,
//   is 97.5 (reflux's 97.48);
// - it prints Pcm = 52.4 bar, where its own (0.286)(8.314)(215.4)/97.5e-6
//   is 52.5 (reflux's 52.55, on its unrounded Tcm, 215.48, which the
//   book's Vc₁₂ and Vcm rounded to 96.3 and 97.5 move by 0.036%);
// - it takes Zcm = 0.285 in the residual, not the 0.286 it prints.
// So its λ − λ° = 1.51e-2 W/(m·K) is reflux's on its inputs as it took
// them, to Γ's rounding (176, not 175.71); reflux's own is 1.49e-2.
fn test_example_10_6_yorizane() {
	ch4 := transport.Pure{
		name:  'methane'
		tc:    190.56
		pc:    45.99e5
		vc:    98.6e-6
		zc:    0.286
		omega: 0.011
		mw:    0.016043
	}
	co2 := transport.Pure{
		name:  'carbon dioxide'
		tc:    304.12
		pc:    73.74e5
		vc:    94.07e-6
		zc:    0.274
		omega: 0.225
		mw:    0.04401
	}
	m := transport.yorizane([ch4, co2], [0.755, 0.245])!
	near(m.vc * 1e6, 97.5, 3, 0, 'Vcm, cm³/mol')
	near(m.tc, 215.4, 4, 0.00036 * 215.4, 'Tcm')
	near(m.zc, 0.286, 3, 0, 'Zcm')
	near(m.pc * 1e-5, 52.5, 3, 0.00036 * 52.5, 'Pcm, bar, as the book works it')
	near(m.mw * 1e3, 22.9, 3, 0, 'Mm')
	near(transport.stiel_thodos(m, 159e-6)!, 1.49e-2, 3, 0, "λ − λ°, reflux's")
	as_taken := transport.Pseudo{
		tc: 215.4
		pc: 52.4e5
		vc: 97.5e-6
		zc: 0.285
		mw: 22.9e-3
	}
	near(transport.stiel_thodos(as_taken, 159e-6)!, 1.51e-2, 3, 0.0017 * 1.51e-2, 'λ − λ°, the book as it took them')
}

// Example 9-15 (p. 9.55): methylcyclohexane at 300 K and 500 bar by
// Lucas's compressed-liquid factor (D6), its vapour pressure neglected as
// the book does: A = 0.9822, D = 0.1371, C = 0.0619, η/η_SL = 1.57, and η
// = 1.07 cP on η_SL = 0.68.
fn test_example_9_15_lucas_compressed_liquid() {
	f := transport.lucas_liquid(300.0, 500e5, 572.19, 34.71e5, 0.235, 0.0)
	near(f, 1.57, 3, 0, 'η/η_SL')
	near(0.68 * f, 1.07, 3, 0, 'η, cP')
}

// Example 10-13 (p. 10.54): toluene at 304 K and 6330 bar by Missenard's
// factor (D7), on λ = 0.129 W/(m·K) at 1 bar: Tr = 0.514, Pr = 154. The book
// reads Q = 0.0205 from its Table 10-7; between the table's rows at Tr 0.5
// and 0.6 (0.020 and 0.025 at Pr 100 and 200) it's 0.0207, which reflux
// takes, so its λ is 0.2197 against the book's 0.219 on 0.0205.
fn test_example_10_13_missenard() {
	f := transport.missenard(304.0, 6330e5, 591.75, 41.08e5)
	pr := 6330.0 / 41.08
	q := (f - 1.0) / math.pow(pr, 0.7)
	near(q, 0.0207, 3, 0, 'Q')
	near(q, 0.0205, 3, 0.0002, "Q against the book's reading")
	near(0.129 * f, 0.219, 3, 0.0007, 'λ')
}

// Example 9-21 (p. 9.79): acetic acid with acetone, 70% acetic acid, at
// 323 K, by Grunberg and Nissan's rule (D6) with Isdale's G: 0.184 at 298 K,
// and at 323 K by its Eq. 9-13.5, 0.258, as a file's G would write it;
// η = 0.588 cP on the book's 0.798 and 0.241 cP.
fn test_example_9_21_grunberg_nissan() {
	p := activity.plain('1 - (1 - 0.184) * (573 K - T) / 275 K', 'G')!
	g := p.e.value(&p.defs, [323.0], []string{})!
	near(g, 0.258, 3, 0, 'G at 323 K')
	mu := transport.grunberg_nissan([0.7, 0.3], [0.798e-3, 0.241e-3], [0.0, g, g, 0.0])
	near(mu * 1e3, 0.588, 3, 0, 'η, cP')
}
