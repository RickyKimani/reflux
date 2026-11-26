// eval answers or refuses (M6, corrected after step 3): an answer it returns
// is finite, in everything asked of it, and right; anything else is an
// error. The models here are binaries made by hand, so that each limit has a
// closed form that owes nothing to the formulas' code.
import activity
import math
import numeric

// least_exponent is the k at which exp(−k) leaves the normal doubles,
// 1022·ln 2: a smaller double has lost digits, and 0 has lost all of them.
const least_exponent = 708.3964185322641

const r1 = 0.9011
const q1 = 0.848
const r2 = 0.92
const q2 = 1.4

// binary is a model of two components at 1 K whose second parameter, for 2
// with 1, is exp(−k), and whose first, for 1 with 2, is exp(−j): NRTL's G
// (with α = 1, so that τ is k and j), UNIQUAC's τ, Wilson's Λ (its first
// and second the other way round, as its formula reads them), and UNIFAC's
// ψ between two groups, a component each.
fn binary(name string, k f64, j f64) !activity.Model {
	names := ['a', 'b']
	minus_j := activity.plain('exp((${-j} K)/T)', 'j')!
	minus_k := activity.plain('exp((${-k} K)/T)', 'k')!
	match name {
		'NRTL' {
			tau12 := activity.plain('(${j} K)/T', 'τ')!
			tau21 := activity.plain('(${k} K)/T', 'τ')!
			one := activity.plain('1', 'α')!
			return activity.Model(activity.new_nrtl(names, [
				activity.NrtlPair{0, 1, tau12, tau21, one},
			])!)
		}
		'UNIQUAC' {
			return activity.Model(activity.new_uniquac(names, [r1, r2], [q1, q2], [
				activity.UniquacPair{0, 1, minus_j, minus_k},
			])!)
		}
		'Wilson' {
			return activity.Model(activity.new_wilson(names, [
				activity.WilsonPair{0, 1, minus_k, minus_j},
			])!)
		}
		else {
			table := activity.Table{
				variant:      .original
				subgroups:    {
					1: activity.Subgroup{1, 'A', 1, r1, q1, 15.0}
					2: activity.Subgroup{2, 'B', 7, r2, q2, 18.0}
				}
				mains:        {
					1: 'A'
					7: 'B'
				}
				interactions: {
					activity.interaction_key(1, 7): activity.Interaction{j, 0.0, 0.0}
					activity.interaction_key(7, 1): activity.Interaction{k, 0.0, 0.0}
				}
			}
			return activity.Model(activity.new_unifac(table, names, [{
				1: 1
			}, {
				2: 1
			}])!)
		}
	}
}

// at_infinite_dilution is binary's ln γ₁ at x = (0, 1), by hand from each
// formula, in k and j themselves, so that no exp(−k) has to be held:
// - NRTL: τ₂₁ + τ₁₂G₁₂;
// - Wilson: 1 − ln Λ₁₂ − Λ₂₁;
// - UNIQUAC and UNIFAC: the combinatorial term, and q₁(1 − ln τ₂₁ − τ₁₂).
fn at_infinite_dilution(name string, k f64, j f64) f64 {
	if name == 'NRTL' {
		return k + j * math.exp(-j)
	}
	residual := 1.0 + k - math.exp(-j)
	if name == 'Wilson' {
		return residual
	}
	v, f := r1 / r2, q1 / q2
	return 1.0 - v + math.log(v) - 5.0 * q1 * (1.0 - v / f + math.log(v / f)) + q1 * residual
}

// its_derivatives is at_infinite_dilution's ∂/∂T and ∂²/∂T² at 1 K, where
// k and j are k/T and j/T: exactly, for each formula's closed form.
fn its_derivatives(name string, k f64, j f64) (f64, f64) {
	g := math.exp(-j)
	if name == 'NRTL' {
		// k/T + u·e^−u, u = j/T.
		return -k - (1.0 - j) * g * j, 2.0 * k + (j - 2.0) * g * j * j + (1.0 - j) * g * 2.0 * j
	}
	scale := if name == 'Wilson' { 1.0 } else { q1 }
	// 1 + k/T − e^(−j/T).
	return scale * (-k - g * j), scale * (2.0 * k - g * (j * j - 2.0 * j))
}

const model_names = ['NRTL', 'UNIQUAC', 'Wilson', 'UNIFAC']

// The whole range of a parameter, through both ends of the doubles: ln γ is
// returned exactly where the parameter is a normal double, and is then the
// closed form's; with derivatives, whatever is returned is finite.
fn test_an_answer_is_finite_and_right_or_refused() {
	mut ks := []f64{}
	for i in -80 .. 81 {
		ks << f64(i) * 10.0
	}
	for i in 1390 .. 1501 {
		ks << f64(i) * 0.5 + 0.3 // through the subnormals, off a whole multiple
		ks << -f64(i) * 0.5
	}
	mut wrong := []string{}
	mut answered, mut refused := 0, 0
	mut worst, mut worst_dt := 0.0, 0.0
	for name in model_names {
		for k in ks {
			for j in [-3.0, 0.0, 0.7, 5.0] {
				m := binary(name, k, j)!
				case := '${name}, k ${k}, j ${j}'
				mut w := m.work()
				mut ok := true
				activity.eval(m, 1.0, [0.0, 1.0], activity.ln_gamma_only, mut w) or { ok = false }
				if ok {
					answered++
					want := at_infinite_dilution(name, k, j)
					gap := math.abs(w.ln_gamma[0] - want) / math.max(math.abs(want), 1.0)
					if !math.is_finite(w.ln_gamma[0]) || !math.is_finite(w.ln_gamma[1])
						|| !(gap <= 1e-13) {
						wrong << '${case}: ${w.ln_gamma}, not ${want}'
					} else {
						worst = math.max(worst, gap)
					}
				} else {
					refused++
				}
				// With T-derivatives likewise: they're finite there, whatever
				// the square of a parameter would be.
				if k >= 0.0 && k < least_exponent - 1e-9 {
					mut wt := m.work()
					mut had := true
					activity.eval(m, 1.0, [0.0, 1.0], .t, mut wt) or {
						wrong << '${case}: T-derivatives refused, ${err}'
						had = false
					}
					if had {
						d, dd := its_derivatives(name, k, j)
						// In units of ε times the terms that cancel to each: of
						// s² to the first and s³ to the second, s = 1 + |k| + |j|,
						// where τ·G and G are differentiated and divided (NRTL's).
						size := 1.0 + math.abs(k) + math.abs(j)
						eps := math.pow(2.0, -52)
						gap := math.max(math.abs(wt.dt[0] - d) / (eps * size * size),
							math.abs(wt.dtt[0] - dd) / (eps * size * size * size))
						if !(gap <= 5.0) { // measured, 2.5 under tcc and gcc
							wrong << '${case}: ${wt.dt[0]}, ${wt.dtt[0]}, not ${d}, ${dd}'
						} else {
							worst_dt = math.max(worst_dt, gap)
						}
					}
				}
				// Above 0 the one limit is exp(−k) leaving the normal doubles.
				if k >= 0.0 && k < least_exponent - 1e-9 && !ok {
					wrong << '${case}: refused'
				}
				if k > least_exponent + 1e-9 && ok {
					wrong << '${case}: answered on a parameter below the normal doubles'
				}
				for x in [[0.0, 1.0], [0.4, 0.6], [1.0, 0.0]] {
					mut all := m.work()
					mut given := true
					activity.eval(m, 1.0, x, .t | .n, mut all) or { given = false }
					if !given {
						continue
					}
					for set in [all.ln_gamma, all.dt, all.dtt, all.dn] {
						if set.any(!math.is_finite(it)) {
							wrong << '${case}, x ${x}: returned ${set}'
						}
					}
				}
			}
		}
	}
	println('answered ${answered}, refused ${refused}; ln γ∞ against its closed form ${worst:.1e}, its T-derivatives ${worst_dt:.1f} units of the terms it cancels from')
	for line in wrong#[..6] {
		println('  ${line}')
	}
	// Measured, 7.1e-15: each a few terms of up to |k|, rounded, so a few tens
	// of ε of the answer. 1e-13 allows a compiler's own rounding and
	// is far below a subnormal parameter's lost digits.
	assert wrong.len == 0, '${wrong.len} wrong'
}

// What was found: NRTL's G = exp(900), an overflow, at any composition.
fn test_an_overflowing_g_is_refused_by_name() {
	m := activity.new_nrtl(['water', 'oil'], [
		activity.NrtlPair{0, 1, activity.plain('(-3000 K)/T', 'τ')!, activity.plain('0', 'τ')!, activity.plain('0.3', 'α')!},
	])!
	mut w := m.work()
	if _ := activity.eval(m, 1.0, [0.5, 0.5], activity.ln_gamma_only, mut w) {
		assert false, 'returned ${w.ln_gamma}'
	} else {
		assert err is activity.ActivityError, err.msg()
		assert err.msg().contains('G of water with oil'), err.msg()
		assert err.msg().contains('at 1.0 K'), err.msg()
	}
}

// A G whose value is 1 and whose second T-derivative overflows is refused
// as τ's own would be, whatever was asked.
fn test_a_parameters_derivatives_are_finite() {
	m := activity.new_nrtl(['a', 'b'], [
		activity.NrtlPair{0, 1, activity.plain('1e160*(T/(300 K) - 1)', 'τ')!, activity.plain('0', 'τ')!, activity.plain('0.3', 'α')!},
	])!
	mut w := m.work()
	if _ := activity.eval(m, 300.0, [0.5, 0.5], activity.ln_gamma_only, mut w) {
		assert false, 'returned ${w.ln_gamma}'
	} else {
		assert err.msg().contains('G of a with b'), err.msg()
	}
}

// A parameter that a formula takes the log of, or divides by, can't be
// negative: UNIQUAC's τ answered at some compositions and not at others.
fn test_a_negative_tau_is_refused() {
	m := activity.new_uniquac(['a', 'b'], [r1, r2], [q1, q2], [
		activity.UniquacPair{0, 1, activity.plain('-0.1', 'τ')!, activity.plain('1', 'τ')!},
	])!
	mut w := m.work()
	if _ := activity.eval(m, 300.0, [0.5, 0.5], activity.ln_gamma_only, mut w) {
		assert false, 'returned ${w.ln_gamma}'
	} else {
		assert err.msg().contains('τ of a with b'), err.msg()
	}
}

// Mole fractions that don't sum to 1, all zero among them, in every model;
// and those that do to rounding, which aren't refused.
fn test_mole_fractions_sum_to_one() {
	mut wrong := []string{}
	for name in model_names {
		m := binary(name, 0.8, 0.3)!
		mut w := m.work()
		for x in [[0.0, 0.0], [0.8, 1.2], [0.2, 0.3], [0.4, 0.6 + 1e-11], [0.4, 0.6 - 1e-11],
			[1e308, 1e308]] {
			mut said := ''
			activity.eval(m, 1.0, x, activity.ln_gamma_only, mut w) or { said = err.msg() }
			if !said.contains('sum to') {
				wrong << '${name}, x ${x}: "${said}", ${w.ln_gamma}'
			}
		}
		for d in [0.0, 1e-14, -1e-14, 5e-13, -5e-13] {
			activity.eval(m, 1.0, [0.4, 0.6 + d], .t | .n, mut w) or {
				wrong << '${name}: refused a sum ${d} from 1, ${err}'
			}
		}
	}
	assert wrong == [], wrong.join('; ')
}

// Broken is a model whose formula, or whose formula on jets, gives what
// isn't finite for component at, as eval must catch whatever model it is.
struct Broken {
	inner activity.Model
	which string
	at    int
	bad   f64
}

fn (b &Broken) size() int {
	return b.inner.size()
}

fn (b &Broken) name() string {
	return 'broken'
}

fn (b &Broken) work() activity.Work {
	return b.inner.work()
}

fn (b &Broken) params(t f64, mut w activity.Work) ! {
	b.inner.params(t, mut w)!
}

fn (b &Broken) lng(p []f64, x []f64, mut s []f64, mut out []f64) {
	b.inner.lng(p, x, mut s, mut out)
	if b.which == 'value' {
		out[b.at] = b.bad
	}
}

fn (b &Broken) lng_jet(p []numeric.Jet, x []numeric.Jet, mut s []numeric.Jet, mut out []numeric.Jet) {
	b.inner.lng_jet(p, x, mut s, mut out)
	if b.which == 'first' {
		out[b.at] = numeric.Jet{
			...out[b.at]
			d: b.bad
		}
	}
	if b.which == 'second' {
		out[b.at] = numeric.Jet{
			...out[b.at]
			dd: b.bad
		}
	}
}

// Each output is checked when it's asked for, and only then, for each
// component, NaN and either infinity; and a Work that was refused gives a
// fresh one's bits afterwards.
fn test_every_output_asked_for_is_checked() {
	inner := binary('NRTL', 0.8, 0.3)!
	x := [0.4, 0.6]
	mut fresh := inner.work()
	activity.eval(inner, 1.0, x, .t | .n, mut fresh)!
	for at in 0 .. 2 {
		for bad in [math.nan(), math.inf(1), math.inf(-1)] {
			// What's broken, then whether eval refuses it asked for ln γ
			// alone, for its T-derivatives, and for dn.
			for which, want in {
				'value':  [true, true, true]
				'first':  [false, true, true]
				'second': [false, true, false]
			} {
				m := activity.Model(Broken{inner, which, at, bad})
				for k, asked in [activity.ln_gamma_only, activity.Want.t, activity.Want.n] {
					mut w := m.work()
					mut refused := false
					activity.eval(m, 1.0, x, asked, mut w) or {
						assert err is activity.ActivityError, err.msg()
						assert err.msg().contains('component ${at}'), err.msg()
						refused = true
					}
					assert refused == want[k], '${which} ${bad} of ${at}, asked ${k}'
					activity.eval(inner, 1.0, x, .t | .n, mut w)!
					assert w.ln_gamma == fresh.ln_gamma && w.dt == fresh.dt && w.dn == fresh.dn
				}
			}
		}
	}
}
