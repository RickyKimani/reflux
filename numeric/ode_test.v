// The integrators (M9 step 8, D8, criterion 6): each solution's order, by
// one step at halving sizes; RODAS4's stability function; closed forms,
// stiff and not; a step that fails and is retried; the profile; the counts;
// and everything refused.
module numeric

import math

// Cosine is y' = y·cos t, y = y₀·exp(sin t): smooth, and with t in f, so
// rodas4's ∂f/∂t is used.
struct Cosine {
mut:
	calls int
}

fn (mut c Cosine) eval(t f64, y []f64, mut dy []f64) ! {
	c.calls++
	dy[0] = y[0] * math.cos(t)
}

fn (mut c Cosine) jacobian(t f64, y []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0] = math.cos(t)
	dfdt[0] = -y[0] * math.sin(t)
}

// Logistic is y' = y(1 − y), y = 1/(1 + (1/y₀ − 1)·e^(t₀−t)): nonlinear,
// with no t in f.
struct Logistic {}

fn (mut l Logistic) eval(_ f64, y []f64, mut dy []f64) ! {
	dy[0] = y[0] * (1 - y[0])
}

fn (mut l Logistic) jacobian(_ f64, y []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0] = 1 - 2 * y[0]
	dfdt[0] = 0
}

// Each solution's local error at h, h/2, h/4, … down to h/128, and the
// order its halving shows: log₂ of successive errors' ratio, p + 1 for a
// method of order p, at the smallest pair of steps whose errors are both
// above 1e-13, clear of rounding. The new point's, and the embedded
// solution's: the new point less the error estimate.
fn local_orders(mut ode StiffOde, t f64, y f64, exact fn (f64) f64, rodas bool) (f64, f64) {
	mut main_errs, mut embedded_errs := []f64{}, []f64{}
	for k in 0 .. 8 {
		h := 0.4 / math.pow(2, k)
		mut fresh, mut low := 0.0, 0.0
		if rodas {
			mut r := new_rodas4(1)
			r.prepare(mut ode, t, [y]) or { panic(err) }
			kind, _ := r.attempt(mut ode, t, [y], h)
			assert kind == .none
			fresh, low = r.ynew[0], r.ynew[0] - r.k6[0]
		} else {
			mut d := new_dopri5(1)
			ode.eval(t, [y], mut d.k1) or { panic(err) }
			assert d.attempt(mut ode, t, [y], h) == ''
			fresh, low = d.y1[0], d.y1[0] - d.est[0]
		}
		main_errs << math.abs(fresh - exact(t + h))
		embedded_errs << math.abs(low - exact(t + h))
	}
	return last_order(main_errs), last_order(embedded_errs)
}

fn last_order(errs []f64) f64 {
	mut p := math.nan()
	for k in 0 .. errs.len - 1 {
		if errs[k] > 1e-13 && errs[k + 1] > 1e-13 {
			p = math.log2(errs[k] / errs[k + 1])
		}
	}
	return p
}

// RODAS4's solution is of order 4 and its embedded one of order 3, so their
// local errors fall 2⁵ and 2⁴ times as the step halves; Dormand–Prince's
// are of orders 5 and 4, 2⁶ and 2⁵. On a problem with t in f and one
// without, each order is within 0.35 of its claim at the smallest steps
// clear of rounding, where h's next power still shows a little (measured:
// within 0.14 for RODAS4, 0.29 for Dormand–Prince's fifth-order solution),
// and a method an order lower would show a whole unit less.
fn test_orders() {
	t0 := 0.3
	cosine := fn [t0] (t f64) f64 {
		return math.exp(math.sin(t) - math.sin(t0))
	}
	logistic := fn [t0] (t f64) f64 {
		return 1 / (1 + 3 * math.exp(t0 - t))
	}
	for rodas in [true, false] {
		want_main := if rodas { 5.0 } else { 6.0 }
		mut c := Cosine{}
		mut s := Logistic{}
		for which in 0 .. 2 {
			main_p, embedded_p := if which == 0 {
				local_orders(mut c, t0, 1, cosine, rodas)
			} else {
				local_orders(mut s, t0, 0.25, logistic, rodas)
			}
			println('${if rodas { 'rodas4' } else { 'dopri5' }}, ${if which == 0 {
				'cosine'
			} else {
				'logistic'
			}}: ${main_p} and ${embedded_p}')
			assert math.abs(main_p - want_main) <= 0.35, '${main_p}'
			assert math.abs(embedded_p - (want_main - 1)) <= 0.35, '${embedded_p}'
		}
	}
}

// Decay is y' = −y.
struct Decay {}

fn (mut d Decay) eval(_ f64, y []f64, mut dy []f64) ! {
	dy[0] = -y[0]
}

// The first step is Hairer's HINIT, worked by hand for y' = −y from y = 1
// at rtol = atol = 1e-6: every sᵢ is 2e-6, so h₀ = 0.01·‖y‖/‖f‖ = 0.01, the
// Euler step's ‖f₁ − f₀‖/h₀ is 5e5, as ‖f₀‖ is, and h = (0.01/5e5)^(1/p),
// p the order, under 100·h₀ and hmax; backwards, the same with its sign.
fn test_first_step_is_hairers() {
	w := weights(0, [1.0], 10, OdeParams{ rtol: 1e-6, atol: 1e-6 })!
	mut d := Decay{}
	for order in [4, 5] {
		want := math.pow(2e-8, 1.0 / f64(order))
		forwards := hinit(mut d, 0, [1.0], [-1.0], 1, order, 10, w)
		backwards := hinit(mut d, 0, [1.0], [-1.0], -1, order, 10, w)
		assert math.abs(forwards / want - 1) <= 1e-12, '${forwards} vs ${want}'
		assert math.abs(backwards / -want - 1) <= 1e-12, '${backwards}'
	}
	assert hinit(mut d, 0, [1.0], [-1.0], 1, 5, 0.001, w) == 0.001
}

// Rotation is y' = A·y with A = [[a, −b], [b, a]], whose eigenvalues are
// a ± ib: one step multiplies |y| by RODAS4's |R(hλ)|.
struct Rotation {
	a f64
	b f64
}

fn (mut r Rotation) eval(_ f64, y []f64, mut dy []f64) ! {
	dy[0] = r.a * y[0] - r.b * y[1]
	dy[1] = r.b * y[0] + r.a * y[1]
}

fn (mut r Rotation) jacobian(_ f64, _ []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0], jac[1], jac[2], jac[3] = r.a, -r.b, r.b, r.a
	dfdt[0], dfdt[1] = 0, 0
}

// RODAS4 is A-stable, |R(z)| ≤ 1 over the left half-plane (here 7 radii
// from 0.01 to 1e12 at 12 angles from the imaginary axis to the negative
// real one), and L-stable: R(z) → 0 as z → −∞, being stiffly accurate. A
// step of hλ = −1e12 leaves 1e-11 of y at most.
fn test_rodas4_is_l_stable() {
	mut worst := 0.0
	for r in [0.01, 1, 10, 100, 1e4, 1e8, 1e12] {
		for k in 0 .. 12 {
			theta := math.pi / 2 + math.pi / 2 * f64(k + 1) / 12
			mut p := Rotation{r * math.cos(theta), r * math.sin(theta)}
			mut ro := new_rodas4(2)
			ro.prepare(mut p, 0, [1.0, 0]) or { panic(err) }
			kind, _ := ro.attempt(mut p, 0, [1.0, 0], 1)
			assert kind == .none
			g := math.hypot(ro.ynew[0], ro.ynew[1])
			worst = math.max(worst, g)
			assert g <= 1 + 1e-14, 'r ${r}, θ ${theta}: |R| = ${g}'
		}
	}
	mut far := Rotation{-1e12, 0}
	mut ro := new_rodas4(2)
	ro.prepare(mut far, 0, [1.0, 0]) or { panic(err) }
	_, _ := ro.attempt(mut far, 0, [1.0, 0], 1)
	assert math.abs(ro.ynew[0]) <= 1e-11, '${ro.ynew[0]}'
	println('rodas4: |R| at most ${worst} over the left half-plane, R(−1e12) = ${ro.ynew[0]}')
}

// Linear is y' = A·y with A = S·diag(−1, −1e3, −1e6)·S⁻¹: stiff, with
// y(t) = S·exp(Λt)·S⁻¹·y₀ in closed form.
struct Linear {}

const lin_s = [[1.0, 1, 0], [0.0, 1, 1], [1.0, 0, 1]]
// S⁻¹, exactly: S's inverse is a half of [[1, −1, 1], [1, 1, −1], [−1, 1, 1]].
const lin_si = [[0.5, -0.5, 0.5], [0.5, 0.5, -0.5], [-0.5, 0.5, 0.5]]
const lin_l = [-1.0, -1e3, -1e6]

fn lin_a(i int, j int) f64 {
	mut x := 0.0
	for k in 0 .. 3 {
		x += lin_s[i][k] * lin_l[k] * lin_si[k][j]
	}
	return x
}

fn (mut l Linear) eval(_ f64, y []f64, mut dy []f64) ! {
	for i in 0 .. 3 {
		dy[i] = lin_a(i, 0) * y[0] + lin_a(i, 1) * y[1] + lin_a(i, 2) * y[2]
	}
}

fn (mut l Linear) jacobian(_ f64, _ []f64, mut jac []f64, mut dfdt []f64) ! {
	for i in 0 .. 3 {
		for j in 0 .. 3 {
			jac[i * 3 + j] = lin_a(i, j)
		}
		dfdt[i] = 0
	}
}

fn lin_exact(t f64, y0 []f64) []f64 {
	mut out := []f64{len: 3}
	for i in 0 .. 3 {
		for k in 0 .. 3 {
			c := lin_si[k][0] * y0[0] + lin_si[k][1] * y0[1] + lin_si[k][2] * y0[2]
			out[i] += lin_s[i][k] * math.exp(lin_l[k] * t) * c
		}
	}
	return out
}

// Prothero and Robinson's problem, y' = λ(y − g(t)) + g'(t), y = g(t) =
// sin t from y(0) = 0, λ = −1e6: stiff, with t in f.
struct Prothero {}

fn (mut p Prothero) eval(t f64, y []f64, mut dy []f64) ! {
	dy[0] = -1e6 * (y[0] - math.sin(t)) + math.cos(t)
}

fn (mut p Prothero) jacobian(t f64, _ []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0] = -1e6
	dfdt[0] = 1e6 * math.cos(t) - math.sin(t)
}

// error_units is the largest error in units of the tolerance:
// |y − exact|/(atol + rtol·|exact|).
fn error_units(y []f64, exact []f64, rtol f64, atol f64) f64 {
	mut worst := 0.0
	for i, x in exact {
		worst = math.max(worst, math.abs(y[i] - x) / (atol + rtol * math.abs(x)))
	}
	return worst
}

// Stiff closed forms: the linear system from t = 0 to 10, and Prothero and
// Robinson's to t = 10, at rtol = atol of 1e-4, 1e-7 and 1e-10, each within
// a tolerance unit of the closed form (measured: 0.16 at worst), in steps
// recorded, unable to grow, and dopri5 refusing the linear one as stiff.
fn test_stiff_closed_forms() {
	y0 := [1.0, 2, 3]
	mut worst := 0.0
	linear_steps := [29, 128, 689]
	prothero_steps := [8, 21, 4420]
	for i, tol in [1e-4, 1e-7, 1e-10] {
		mut l := Linear{}
		a := rodas4(mut l, 0, y0, 10, rtol: tol, atol: tol)!
		ua := error_units(a.y, lin_exact(10, y0), tol, tol)
		mut p := Prothero{}
		b := rodas4(mut p, 0, [0.0], 10, rtol: tol, atol: tol)!
		ub := error_units(b.y, [math.sin(10)], tol, tol)
		println('rtol ${tol}: linear ${ua} units in ${a.steps} steps, Prothero–Robinson ${ub} in ${b.steps}')
		worst = math.max(worst, math.max(ua, ub))
		assert a.steps <= linear_steps[i] && b.steps <= prothero_steps[i]
	}
	assert worst <= 1, '${worst}'
	mut l := Linear{}
	dopri5(mut l, 0, y0, 10, rtol: 1e-6, atol: 1e-6) or {
		assert err is OdeError && (err as OdeError).kind == .stiff, err.msg()
		return
	}
	assert false, 'dopri5 integrated a stiff problem'
}

// Kepler is the two-body problem in the plane, eccentricity e: y = (x, y,
// x', y'), with y'' = −y/r³.
struct Kepler {}

fn (mut k Kepler) eval(_ f64, y []f64, mut dy []f64) ! {
	r3 := math.pow(y[0] * y[0] + y[1] * y[1], 1.5)
	dy[0], dy[1] = y[2], y[3]
	dy[2], dy[3] = -y[0] / r3, -y[1] / r3
}

// kepler_exact is the orbit at t from Kepler's equation E − e·sin E = t,
// solved by Newton to rounding: x = cos E − e, y = √(1 − e²)·sin E.
fn kepler_exact(e f64, t f64) []f64 {
	mut ea := t
	for _ in 0 .. 50 {
		ea -= (ea - e * math.sin(ea) - t) / (1 - e * math.cos(ea))
	}
	s := math.sqrt(1 - e * e)
	d := 1 - e * math.cos(ea)
	return [math.cos(ea) - e, s * math.sin(ea), -math.sin(ea) / d, s * math.cos(ea) / d]
}

// Arenstorf's orbit, the restricted three-body problem of Hairer, Nørsett
// and Wanner I, §II.0, which comes back to its start after its period.
struct Arenstorf {}

const ar_mu = 0.012277471

fn (mut a Arenstorf) eval(_ f64, y []f64, mut dy []f64) ! {
	mu := ar_mu
	mup := 1 - mu
	d1 := math.pow((y[0] + mu) * (y[0] + mu) + y[1] * y[1], 1.5)
	d2 := math.pow((y[0] - mup) * (y[0] - mup) + y[1] * y[1], 1.5)
	dy[0], dy[1] = y[2], y[3]
	dy[2] = y[0] + 2 * y[3] - mup * (y[0] + mu) / d1 - mu * (y[0] - mup) / d2
	dy[3] = y[1] - 2 * y[2] - mup * y[1] / d1 - mu * y[1] / d2
}

const ar_y0 = [0.994, 0, 0, -2.00158510637908252240537862224]
const ar_period = 17.0652165601579625588917206249

// Non-stiff closed forms for dopri5: Kepler's orbit with e = 0.5 over three
// periods, and Arenstorf's back to its start after one, at rtol = atol of
// 1e-4, 1e-7 and 1e-10. A step's error is held to the tolerance, not the
// orbit's: its phase drifts along it, and Arenstorf's close passes by the
// moon amplify every error, so their errors at the end are many tolerances
// (Kepler's 730 to 4,400 units, Arenstorf's 1,100 to 24,000; the errors
// themselves 0.44 to 7.3e-8 and 0.32 to 2.4e-6). Each falls as
// the tolerance tightens, and each is recorded, unable to grow: the largest
// |y − exact| at each tolerance, 10% above what was measured, and the
// steps as measured.
fn test_non_stiff_closed_forms() {
	e := 0.5
	k0 := kepler_exact(e, 0)
	t1 := 6 * math.pi
	tols := [1e-4, 1e-7, 1e-10]
	kepler_errors := [0.48, 1.4e-4, 8.1e-8]
	arenstorf_errors := [0.36, 1.6e-3, 2.7e-6]
	kepler_steps := [42, 137, 539]
	arenstorf_steps := [64, 216, 841]
	mut last_k, mut last_a := math.inf(1), math.inf(1)
	for i, tol in tols {
		mut k := Kepler{}
		a := dopri5(mut k, 0, k0, t1, rtol: tol, atol: tol)!
		ek := largest_difference(a.y, kepler_exact(e, t1))
		mut ar := Arenstorf{}
		b := dopri5(mut ar, 0, ar_y0, ar_period, rtol: tol, atol: tol)!
		ea := largest_difference(b.y, ar_y0)
		println('rtol ${tol}: Kepler ${ek} (${error_units(a.y, kepler_exact(e, t1), tol, tol)} units) in ${a.steps} steps, Arenstorf ${ea} (${error_units(b.y, ar_y0, tol, tol)} units) in ${b.steps}')
		assert ek <= kepler_errors[i] && ek < last_k, '${ek}'
		assert a.steps <= kepler_steps[i] && b.steps <= arenstorf_steps[i]
		assert ea <= arenstorf_errors[i] && ea < last_a, '${ea}'
		last_k, last_a = ek, ea
	}
}

fn largest_difference(a []f64, b []f64) f64 {
	mut worst := 0.0
	for i, x in a {
		worst = math.max(worst, math.abs(x - b[i]))
	}
	return worst
}

// Root is y' = −√y, y = (√y₀ − t/2)², whose f fails below zero: a step
// that would overshoot zero is rejected and halved, and the integration
// carries on to the closed form.
struct Root {
mut:
	failed int
}

fn (mut r Root) eval(_ f64, y []f64, mut dy []f64) ! {
	if y[0] < 0 {
		r.failed++
		return error('√ of ${y[0]}')
	}
	dy[0] = -math.sqrt(y[0])
}

fn (mut r Root) jacobian(_ f64, y []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0] = -0.5 / math.sqrt(y[0])
	dfdt[0] = 0
}

// The first step is the whole interval, 1.9, and a stage there falls below
// zero; each integrator halves it, rejecting 3 steps in all as measured
// (the failed one and two too long), and lands within 1e-9 of the closed
// form at rtol 1e-8 (measured: 4.3e-11 and 1.3e-10).
fn test_a_failing_evaluation_is_retried() {
	exact := math.pow(1 - 1.9 / 2, 2)
	for rodas in [true, false] {
		mut r := Root{}
		s := if rodas {
			rodas4(mut r, 0, [1.0], 1.9, rtol: 1e-8, atol: 1e-12, h0: 1.9)!
		} else {
			dopri5(mut r, 0, [1.0], 1.9, rtol: 1e-8, atol: 1e-12, h0: 1.9)!
		}
		assert r.failed > 0 && s.rejected > 0 && s.rejected <= 3, '${s.rejected} rejected'
		assert math.abs(s.y[0] - exact) <= 1e-9, '${s.y[0]} vs ${exact}'
	}
}

// Backwards in t, the profile, and the counts: every evaluation the
// problem saw, and for rodas4 a Jacobian at every accepted point but the
// last.
fn test_profile_direction_and_counts() {
	mut c := Cosine{}
	s := rodas4(mut c, 2, [1.0], 0, rtol: 1e-9, atol: 1e-12, record: true)!
	assert s.t == 0
	assert math.abs(s.y[0] - math.exp(math.sin(0) - math.sin(2))) <= 1e-8
	assert s.evaluations == c.calls
	assert s.jacobians == s.steps
	assert s.factorisations == s.steps + s.rejected
	assert s.ts.len == s.steps + 1 && s.ys.len == s.steps + 1
	assert s.ts[0] == 2 && s.ys[0] == [1.0] && s.ts.last() == 0 && s.ys.last() == s.y
	for k in 1 .. s.ts.len {
		assert s.ts[k] < s.ts[k - 1]
	}
	mut d := Cosine{}
	u := dopri5(mut d, 2, [1.0], 0, rtol: 1e-9, atol: 1e-12)!
	assert u.t == 0 && u.evaluations == d.calls && u.ts.len == 0
	assert math.abs(u.y[0] - math.exp(math.sin(0) - math.sin(2))) <= 1e-8
}

// hmax bounds every step, each integrator's, though the tolerance would
// allow far longer ones.
fn test_hmax_bounds_every_step() {
	for rodas in [true, false] {
		mut c := Cosine{}
		s := if rodas {
			rodas4(mut c, 0, [1.0], 1, rtol: 1e-3, atol: 1e-3, hmax: 0.05, record: true)!
		} else {
			dopri5(mut c, 0, [1.0], 1, rtol: 1e-3, atol: 1e-3, hmax: 0.05, record: true)!
		}
		assert s.steps >= 20
		for k in 1 .. s.ts.len {
			assert s.ts[k] - s.ts[k - 1] <= 0.05 * (1 + 1e-12), '${s.ts[k] - s.ts[k - 1]}'
		}
	}
}

// Bad is a problem whose Jacobian makes I/(hγ) − J singular at any step,
// rank one and too large for 1/(hγ) to show in it.
struct Bad {
	fail_at_start bool
}

fn (mut b Bad) eval(_ f64, y []f64, mut dy []f64) ! {
	if b.fail_at_start {
		return error('no f here')
	}
	dy[0], dy[1] = y[0], y[1]
}

fn (mut b Bad) jacobian(_ f64, _ []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0], jac[1], jac[2], jac[3] = 1e300, 1e300, 1e300, 1e300
	dfdt[0], dfdt[1] = 0, 0
}

// Unfinite is a problem not finite at its start: its f, or only its
// Jacobian.
struct Unfinite {
	jacobian_only bool
}

fn (mut u Unfinite) eval(_ f64, y []f64, mut dy []f64) ! {
	dy[0] = if u.jacobian_only { -y[0] } else { math.inf(1) }
}

fn (mut u Unfinite) jacobian(_ f64, _ []f64, mut jac []f64, mut dfdt []f64) ! {
	jac[0] = math.nan()
	dfdt[0] = 0
}

// Blowup is y' = y², which leaves for infinity at t = 1 from y(0) = 1.
struct Blowup {}

fn (mut b Blowup) eval(_ f64, y []f64, mut dy []f64) ! {
	dy[0] = y[0] * y[0]
}

fn kind_of(f fn () !OdeSolution) string {
	s := f() or {
		if err is OdeError {
			return '${err.kind}: ${err.msg()}'
		}
		return 'not an OdeError: ${err.msg()}'
	}
	return 'answered ${s.y}'
}

fn test_refused() {
	mut c := Cosine{}
	assert kind_of(fn [mut c] () !OdeSolution {
		return dopri5(mut c, 0, [1.0], 1, rtol: 0)
	}).starts_with('bad_input: ode: rtol 0.0 is within rounding')
	assert kind_of(fn [mut c] () !OdeSolution {
		return dopri5(mut c, 0, [1.0], 1, atol: -1)
	}).starts_with('bad_input: ode: atol -1.0 for y[0] is not above zero')
	assert kind_of(fn [mut c] () !OdeSolution {
		return rodas4(mut c, 0, [1.0], 1, atols: [1.0, 2])
	}).starts_with('bad_input: ode: 2 absolute tolerances for 1 equations')
	assert kind_of(fn [mut c] () !OdeSolution {
		return rodas4(mut c, 1, [1.0], 1)
	}).starts_with('bad_input: ode: no interval from 1.0 to 1.0')
	assert kind_of(fn [mut c] () !OdeSolution {
		return rodas4(mut c, 0, [math.nan()], 1)
	}).starts_with('bad_input: ode: y[0] is nan at the start')
	assert kind_of(fn [mut c] () !OdeSolution {
		return dopri5(mut c, 0, []f64{}, 1)
	}).starts_with('bad_input: ode: no equations')
	assert kind_of(fn [mut c] () !OdeSolution {
		return dopri5(mut c, 0, [1.0], 100, rtol: 1e-12, atol: 1e-12, max_steps: 10)
	}).starts_with('too_many_steps: ode: more than 10 steps at t = ')
	mut b := Bad{
		fail_at_start: true
	}
	assert kind_of(fn [mut b] () !OdeSolution {
		return rodas4(mut b, 0, [1.0, 1], 1)
	}).starts_with('bad_input: ode: no f here at the start')
	mut s := Bad{}
	assert kind_of(fn [mut s] () !OdeSolution {
		return rodas4(mut s, 0, [1.0, 1], 1)
	}).starts_with('singular: ode: I/(hγ) − J is singular, after halving the step 5 times at t = 0.0 (0 steps)')
	mut inf := Unfinite{}
	assert kind_of(fn [mut inf] () !OdeSolution {
		return dopri5(mut inf, 0, [1.0], 1)
	}) == "bad_input: ode: f isn't finite at the start at t = 0.0 (0 steps)"
	assert kind_of(fn [mut inf] () !OdeSolution {
		return rodas4(mut inf, 0, [1.0], 1)
	}) == "bad_input: ode: f isn't finite at the start at t = 0.0 (0 steps)"
	mut nan := Unfinite{
		jacobian_only: true
	}
	assert kind_of(fn [mut nan] () !OdeSolution {
		return rodas4(mut nan, 0, [1.0], 1)
	}) == "bad_input: ode: the Jacobian isn't finite at the start at t = 0.0 (0 steps)"
	mut up := Blowup{}
	got := kind_of(fn [mut up] () !OdeSolution {
		return dopri5(mut up, 0, [1.0], 2)
	})
	assert got.starts_with('step_too_small: ode: the step fell to '), got
}
