module flash

import math
import numeric
import thermo

// VolumeParams are the TV and UV flashes' options.
@[params]
pub struct VolumeParams {
pub:
	// The PT flashes' fugacity tolerance.
	tol f64 = 1e-12
	// Converged when |V − V_spec| ≤ vtol·V_spec,
	vtol f64 = 1e-10
	// and, for UV, when |U − U_spec| ≤ etol·RT.
	etol f64 = 1e-10
	// verify re-tests each phase of a split for stability, as the PT flash
	// does.
	verify    bool = true
	max_steps int  = 100
	// estimate is where to start: a pressure in Pa for tv, a temperature in
	// K for uv.
	estimate f64
}

// tv flashes a feed z at temperature T to its molar volume V, in m³ per
// mole of feed: it finds P and the split. A Model is enough, since the
// ideal gas's part of the Helmholtz energy doesn't depend on the split.
//
// At fixed T and V the equilibrium minimises the Helmholtz energy A, and
// with Q(P) = G − P·V_spec (Michelsen's state function, D5),
// dQ*/dP = V(P) − V_spec, where V(P) is the equilibrium volume, which
// falls as P rises. So:
//
// - if the feed as one phase at (T, V) is stable, its pressure from the
//   model is the answer;
// - otherwise ln P is stepped from P(T, V) (or params.estimate) by ln 2
//   until V(P) − V_spec changes sign, with a PT flash at each, and Brent's
//   method narrows that to 0.1%;
// - if the PT flash there splits, A itself is minimised over the vapour's
//   mole numbers and volume, the liquid taking the rest: no volume root is
//   solved, only F(T, V, n) and its derivatives;
// - the answer is a verified PT flash at the pressure found, held to
//   |V − V_spec| ≤ vtol·V_spec.
//
// A pure component inside its saturation dome has P = Psat(T) and the
// vapour fraction β = (V − Vᴸ)/(Vⱽ − Vᴸ).
pub fn tv(m thermo.Model, t f64, v f64, z []f64, params VolumeParams) !FlashResult {
	if !(t > 0.0) || !math.is_finite(t) {
		return fail(.bad_input, 'T = ${t} K is not a temperature')
	}
	s := volume_flash(m, t, v, z, params)!
	if s.present.len == 1 {
		return s.pure()!
	}
	return s.mixture()!
}

// VolumeFlash is one TV flash.
struct VolumeFlash {
	m       thermo.Model
	t       f64
	v       f64 // V_spec, m³/mol
	z       []f64
	present []int
	params  VolumeParams
}

fn volume_flash(m thermo.Model, t f64, v f64, z []f64, params VolumeParams) !VolumeFlash {
	if !(v > 0.0) || !math.is_finite(v) {
		return fail(.bad_input, 'V = ${v} m³/mol is not a volume')
	}
	if !(params.tol > 0.0) || !(params.vtol > 0.0) || !(params.etol > 0.0) || params.max_steps < 1 {
		return fail(.bad_input, 'tol, vtol, etol and max_steps must be positive')
	}
	if !(params.estimate >= 0.0) || !math.is_finite(params.estimate) {
		return fail(.bad_input, 'estimate = ${params.estimate} must be positive, or 0 for the default start')
	}
	check_feed(m, z)!
	feed := normalised(z)
	if !(v > m.min_volume(feed)) {
		return fail(.bad_input, "V = ${v} m³/mol is below the model's smallest volume for this feed, ${m.min_volume(feed)}")
	}
	return VolumeFlash{
		m:       m
		t:       t
		v:       v
		z:       feed
		present: []int{len: feed.len, init: index}.filter(feed[it] > 0.0)
		params:  params
	}
}

fn (s VolumeFlash) what() string {
	return 'the TV flash to V = ${s.v} m³/mol at ${s.t} K'
}

// residual is (V − V_spec)/V_spec for a flash result.
fn (s VolumeFlash) residual(r FlashResult) f64 {
	mut v := 0.0
	for ph in r.phases {
		v += ph.fraction * ph.state.molar_volume()
	}
	return (v - s.v) / s.v
}

fn (s VolumeFlash) flash_at(lnp f64) !FlashResult {
	p := math.exp(lnp)
	if !(p >= 1e-6 && p <= 1e12) {
		return fail(.bad_input, "${s.what()}: no pressure from 1e-6 Pa to 1e12 Pa gives it, so the specification is outside the model's range")
	}
	return pt(s.m, s.t, p, s.z, tol: s.params.tol, verify: s.params.verify)
}

fn (s VolumeFlash) mixture() !FlashResult {
	p0 := thermo.pressure(s.m, s.t, s.v, s.z)
	// The feed as one phase at its own pressure, if that's stable and the
	// PT flash there finds this volume root. If that PT flash splits
	// instead, its split is a good start for the Helmholtz minimisation.
	if p0 > 0.0 && math.is_finite(p0) && s.params.estimate == 0.0 {
		if r := s.flash_at(math.log(p0)) {
			if r.phases.len == 1 && math.abs(s.residual(r)) <= s.params.vtol {
				return FlashResult{
					...r
					iterations: 1
				}
			}
			if r.phases.len == 2 {
				if p, n := s.helmholtz(r) {
					if q := s.flash_at(math.log(p)) {
						if math.abs(s.residual(q)) <= s.params.vtol {
							return FlashResult{
								...q
								iterations: n + 2
							}
						}
					}
				}
			}
		}
	}
	start := if s.params.estimate > 0.0 {
		s.params.estimate
	} else if p0 > 0.0 && math.is_finite(p0) {
		p0
	} else {
		thermo.r_gas * s.t / s.v
	}
	f := fn [s] (lnp f64) !f64 {
		return s.residual(s.flash_at(lnp)!)
	}
	// V(P) falls as P rises, so the residual does too.
	lo, hi, mut steps := bracket_down(f, math.log(start), math.ln2, '${s.what()}')!
	coarse := numeric.brent(f, lo, hi, xtol: 1e-3, max_iter: s.params.max_steps)!
	steps += coarse.evaluations
	mut r := s.flash_at(coarse.x)!
	steps++
	if r.phases.len == 2 {
		if p, n := s.helmholtz(r) {
			steps += n
			if q := s.flash_at(math.log(p)) {
				r = q
				steps++
			}
		}
	}
	if math.abs(s.residual(r)) > s.params.vtol {
		fine := numeric.brent(f, coarse.a, coarse.b,
			ftol:     s.params.vtol
			max_iter: s.params.max_steps
		)!
		steps += fine.evaluations
		r = s.flash_at(fine.x)!
		steps++
	}
	res := s.residual(r)
	if math.abs(res) > s.params.vtol {
		return fail(.no_convergence, '${s.what()}: stopped at ${r.p} Pa with V off by ${res} relatively (vtol ${s.params.vtol})')
	}
	return FlashResult{
		...r
		iterations: steps
	}
}

// bracket_down steps x from x0 by dx until f, which falls as x rises,
// changes sign, and returns the bracket and the evaluations.
fn bracket_down(f numeric.ScalarFn, x0 f64, dx f64, what string) !(f64, f64, int) {
	mut x := x0
	f0 := f(x)!
	if f0 == 0.0 {
		return x, x, 1
	}
	// f too high: raise x.
	up := f0 > 0.0
	for k in 0 .. 80 {
		next := if up { x + dx } else { x - dx }
		fn_ := f(next)!
		if fn_ == 0.0 || (fn_ > 0.0) != (f0 > 0.0) {
			return if up { x, next, k + 2 } else { next, x, k + 2 }
		}
		x = next
	}
	return fail(.no_convergence, '${what}: no change of sign in 80 steps')
}

// helmholtz minimises A/RT of a split at (T, V) over the vapour's mole
// numbers nⱽ and volume Vⱽ, the liquid taking nᴸ = z − nⱽ and
// Vᴸ = V − Vⱽ, from the PT split r. With each phase's F(T, V, n) and the
// ideal gas's Σ nᵢ(ln(nᵢ/V) − 1) (its T-dependent part is the same for any
// split):
//
//     ∂(A/RT)/∂nᵢⱽ = (F_nᵢ + ln(nᵢ/V))ⱽ − (…)ᴸ = (μᵢⱽ − μᵢᴸ)/RT,
//     ∂(A/RT)/∂Vⱽ = (F_V − n/V)ⱽ − (…)ᴸ = (Pᴸ − Pⱽ)/RT,
//
// and the Hessian is the two phases' added: F_nn + δᵢⱼ/nᵢ, F_Vn − 1/V and
// F_VV + n/V². It returns the pressure of the minimum and the steps taken.
fn (s VolumeFlash) helmholtz(r FlashResult) !(f64, int) {
	y := r.phases[0].state.n
	mut x0 := s.present.map(r.beta * y[it])
	// The liquid, which hardly compresses, keeps its volume and the vapour
	// takes the rest of V, unless that leaves it too little.
	vv := s.v - (1.0 - r.beta) * r.phases[1].state.molar_volume()
	x0 << if vv > s.m.min_volume(s.split_moles(x0)) {
		vv
	} else {
		r.beta * r.phases[0].state.molar_volume()
	}
	objective := fn [s] (x []f64, mut g []f64, mut h []f64) !f64 {
		return s.helmholtz_eval(x, mut g, mut h)
	}
	mut g := []f64{len: x0.len}
	mut h := []f64{len: x0.len * x0.len}
	a0 := s.helmholtz_eval(x0, mut g, mut h)!
	sol := numeric.minimize(objective, x0,
		gtol:     s.params.tol
		fscale:   math.max(math.abs(a0), 1.0)
		max_iter: 50
	)!
	k := s.present.len
	nv := s.split_moles(sol.x)
	return thermo.pressure(s.m, s.t, sol.x[k], nv), sol.iterations
}

// split_moles is the vapour's mole numbers, over every component.
fn (s VolumeFlash) split_moles(x []f64) []f64 {
	mut nv := []f64{len: s.z.len}
	for j, i in s.present {
		nv[i] = x[j]
	}
	return nv
}

fn (s VolumeFlash) helmholtz_eval(x []f64, mut g []f64, mut h []f64) !f64 {
	nc := s.z.len
	k := s.present.len
	n := k + 1
	vv := x[k]
	vl := s.v - vv
	nv := s.split_moles(x)
	mut nl := []f64{len: nc}
	for j, i in s.present {
		if !(x[j] > 0.0) || !(x[j] < s.z[i]) {
			return error('a mole number left (0, zᵢ)')
		}
		nl[i] = s.z[i] - x[j]
	}
	if !(vv > s.m.min_volume(nv)) || !(vl > s.m.min_volume(nl)) {
		return error('a phase volume left its range')
	}
	rv := s.m.residual2(s.t, vv, nv)
	rl := s.m.residual2(s.t, vl, nl)
	totv, totl := numeric.sum(nv), numeric.sum(nl)
	mut a := rv.f + rl.f
	for j, i in s.present {
		a += nv[i] * (math.log(nv[i] / vv) - 1.0) + nl[i] * (math.log(nl[i] / vl) - 1.0)
		g[j] = rv.f_n[i] + math.log(nv[i] / vv) - rl.f_n[i] - math.log(nl[i] / vl)
		for l, q in s.present {
			mut hij := rv.f_nn[i * nc + q] + rl.f_nn[i * nc + q]
			if i == q {
				hij += 1.0 / nv[i] + 1.0 / nl[i]
			}
			h[j * n + l] = hij
		}
		hv := rv.f_vn[i] - 1.0 / vv + rl.f_vn[i] - 1.0 / vl
		h[j * n + k] = hv
		h[k * n + j] = hv
	}
	g[k] = rv.f_v - totv / vv - rl.f_v + totl / vl
	h[k * n + k] = rv.f_vv + totv / (vv * vv) + rl.f_vv + totl / (vl * vl)
	if !math.is_finite(a) {
		return error('A is not finite here')
	}
	return a
}

// pure is a single component's TV flash. Inside the saturation dome at T
// (from tvf), P is Psat and β follows from the volumes; outside it, or past
// the critical temperature, the one phase's pressure is the model's
// P(T, V).
fn (s VolumeFlash) pure() !FlashResult {
	m := s.m
	if sat := tvf(m, s.t, 0.5, s.z, tol: s.params.tol) {
		vap, liq := sat.phases[0].state, sat.phases[1].state
		vv, vl := vap.molar_volume(), liq.molar_volume()
		// Just outside the dome is a single phase at P(T, V), so no margin
		// is needed here, unlike in ph: that pressure is solved directly.
		if s.v >= vl && s.v <= vv {
			beta := (s.v - vl) / (vv - vl)
			return FlashResult{
				...sat
				beta:   beta
				phases: [Phase{.vapor, vap, beta}, Phase{.liquid, liq, 1.0 - beta}]
			}
		}
		return s.single(if s.v < vl { thermo.Root.liquid } else { .vapor })!
	}
	return s.single(.stable)!
}

// single is a pure component's one phase at (T, V): P from the model, and
// the root at that P that has this volume.
fn (s VolumeFlash) single(root thermo.Root) !FlashResult {
	p := thermo.pressure(s.m, s.t, s.v, s.z)
	if !(p > 0.0) || !math.is_finite(p) {
		return fail(.no_convergence, '${s.what()}: the one phase there has P = ${p} Pa')
	}
	st := thermo.solve_tp(s.m, s.t, p, s.z, root)!
	if math.abs(st.v - s.v) > s.params.vtol * s.v {
		return fail(.no_convergence, '${s.what()}: at ${p} Pa the root found has V = ${st.v}, not this one')
	}
	label := thermo.phase_label(s.m, st)
	return FlashResult{
		t:          s.t
		p:          p
		z:          s.z
		beta:       if label == .vapor { 1.0 } else { 0.0 }
		phases:     [Phase{label, st, 1.0}]
		path:       .single_phase
		iterations: 1
		stability:  stability(s.m, s.t, p, s.z)!
	}
}

// uv flashes a feed z to its molar internal energy U, in J per mole of
// feed on the Fluid's reference state, and its molar volume V: it finds T,
// P and the split. At fixed U and V the equilibrium maximises S, and with
// Q(T) = (A − U_spec)/T, dQ*/dT = −(U(T) − U_spec)/T², with U(T) the TV
// flash's U, which rises with T. So ln T is stepped out by ln 1.5 from the
// ideal gas's T for U (or params.estimate) until U − U_spec changes sign,
// with a TV flash at each, and Brent's method narrows it to 0.1%. A split
// there is finished by Newton on (nⱽ, Vⱽ, ln T) together, the TV flash's
// Helmholtz conditions plus the energy balance; otherwise Brent goes on.
// The answer is a verified TV flash at the final T, held to
// |U − U_spec| ≤ etol·RT. iterations counts the TV flashes, the expensive
// part: a Newton step costs a small fraction of one.
pub fn uv(fl thermo.Fluid, u f64, v f64, z []f64, params VolumeParams) !FlashResult {
	if !math.is_finite(u) {
		return fail(.bad_input, 'the specified U = ${u} is not a number')
	}
	// The temperature here only fixes the checks; each step sets its own.
	s := volume_flash(fl.model, 300.0, v, z, params)!
	e := EnergyFlash{
		fl:  fl
		u:   u
		tvs: s
	}
	return e.solve()!
}

// EnergyFlash is one UV flash.
struct EnergyFlash {
	fl  thermo.Fluid
	u   f64 // U_spec, J/mol
	tvs VolumeFlash
}

fn (e EnergyFlash) what() string {
	return 'the UV flash to U = ${e.u} J/mol, V = ${e.tvs.v} m³/mol'
}

// at is the TV flash at temperature T.
fn (e EnergyFlash) at(lnt f64) !FlashResult {
	t := math.exp(lnt)
	if !(t >= state_t_min && t <= state_t_max) {
		return fail(.bad_input, "${e.what()}: no temperature from ${state_t_min} K to ${state_t_max} K reaches it, so the specification is outside the model's range")
	}
	return tv(e.fl.model, t, e.tvs.v, e.tvs.z,
		tol:    e.tvs.params.tol
		vtol:   e.tvs.params.vtol
		verify: e.tvs.params.verify
	)
}

// residual is (U − U_spec)/RT for a flash result.
fn (e EnergyFlash) residual(r FlashResult) f64 {
	mut u := 0.0
	for ph in r.phases {
		u += ph.fraction * e.fl.props(ph.state).u
	}
	return (u - e.u) / (thermo.r_gas * r.t)
}

// ideal_t is the temperature at which the feed as an ideal gas has U_spec.
fn (e EnergyFlash) ideal_t() f64 {
	if e.tvs.params.estimate > 0.0 {
		return e.tvs.params.estimate
	}
	g := fn [e] (t f64) f64 {
		mut u := 0.0
		for i in e.tvs.present {
			ig := e.fl.ideal[i]
			u += e.tvs.z[i] * ((ig.dhf or { 0.0 }) + ig.cp.int_cp(thermo.t_ref, t) - thermo.r_gas * t)
		}
		return u - e.u
	}
	return numeric.bisect(g, state_t_min, state_t_max) or { 300.0 }
}

fn (e EnergyFlash) solve() !FlashResult {
	f := fn [e] (lnt f64) !f64 {
		return e.residual(e.at(lnt)!)
	}
	// U rises with T, so −residual falls.
	neg := fn [f] (lnt f64) !f64 {
		return -f(lnt)!
	}
	// The ideal gas's T can be far too low for a liquid's U, where no TV
	// flash at this V has an answer; then start from 300 K. (Starting with
	// Newton on the one phase's U(T) would be quicker, but it can run far:
	// from a split's U it found a one-phase answer at 5035 K, where the
	// test data's quadratic Cp has turned negative. The bracket finds the
	// nearest change of sign instead.)
	mut t0 := e.ideal_t()
	mut steps := 0
	if _ := e.at(math.log(t0)) {
	} else {
		t0 = 300.0
	}
	lo, hi, bsteps := bracket_down(neg, math.log(t0), math.log(1.5), e.what())!
	steps += bsteps
	coarse := numeric.brent(f, lo, hi, xtol: 1e-3, max_iter: e.tvs.params.max_steps)!
	steps += coarse.evaluations
	mut r := e.at(coarse.x)!
	steps++
	// Finish with Newton: on the split's Helmholtz conditions and the
	// energy balance together, or, for one phase, on U(T) alone at this V.
	mut t_fin := 0.0
	if r.phases.len == 2 {
		if t, _ := e.newton(r) {
			t_fin = t
		}
	} else {
		if t, _ := e.single_newton(r.t) {
			t_fin = t
		}
	}
	if t_fin > 0.0 {
		if q := e.at(math.log(t_fin)) {
			r = q
			steps++
		}
	}
	if math.abs(e.residual(r)) > e.tvs.params.etol {
		fine := numeric.brent(f, coarse.a, coarse.b,
			ftol:     e.tvs.params.etol
			max_iter: e.tvs.params.max_steps
		)!
		steps += fine.evaluations
		r = e.at(fine.x)!
		steps++
	}
	res := e.residual(r)
	if math.abs(res) > e.tvs.params.etol {
		return fail(.no_convergence, '${e.what()}: stopped at ${r.t} K with the energy balance off by ${res} (etol ${e.tvs.params.etol})')
	}
	return FlashResult{
		...r
		iterations: steps
	}
}

// single_newton solves U(T) = U_spec for the feed as one phase at V, by
// Newton in ln T from t0, from F alone: U = Σᵢ zᵢ(ΔH°f,ᵢ + ∫Cpᵢ dT − RT) −
// RT²·F_T, with dU/dT = Cv. It returns the temperature and the steps.
fn (e EnergyFlash) single_newton(t0 f64) !(f64, int) {
	system := fn [e] (x []f64, mut f []f64, mut jac []f64) ! {
		t := math.exp(x[0])
		if !(t >= state_t_min && t <= state_t_max) {
			return error('T left the range searched')
		}
		r := e.tvs.m.residual2(t, e.tvs.v, e.tvs.z)
		mut uig, mut cvig := 0.0, 0.0
		for i in e.tvs.present {
			ig := e.fl.ideal[i]
			uig += e.tvs.z[i] * ((ig.dhf or { 0.0 }) + ig.cp.int_cp(thermo.t_ref, t) - thermo.r_gas * t)
			cvig += e.tvs.z[i] * (ig.cp.cp(t) - thermo.r_gas)
		}
		rt := thermo.r_gas * t
		f[0] = (uig - rt * t * r.f_t - e.u) / rt
		cv := cvig - rt * (t * r.f_tt + 2.0 * r.f_t)
		jac[0] = cv / thermo.r_gas - f[0]
	}
	sol := numeric.newton(system, [math.log(t0)],
		ftol:     e.tvs.params.etol / 10.0
		max_iter: e.tvs.params.max_steps
	)!
	return math.exp(sol.x[0]), sol.iterations
}

// newton solves the UV conditions in (nⱽ, Vⱽ, ln T) together, from the
// TV split r: the Helmholtz gradient of VolumeFlash.helmholtz_eval, and
// r_U = (U − U_spec)/RT with U = Σᵢ zᵢ(ΔH°f,ᵢ + ∫Cpᵢ dT − RT) − RT²·Σ F_T over
// the phases. With cᵢ = T(F_Tnᵢⱽ − F_Tnᵢᴸ) and c_V = T(F_TVⱽ − F_TVᴸ):
//
//     ∂gᵢ/∂ln T = cᵢ,   ∂g_V/∂ln T = c_V,
//     ∂r_U/∂nᵢⱽ = −cᵢ,  ∂r_U/∂Vⱽ = −c_V,  ∂r_U/∂ln T = Cv/R − r_U,
//
// with Cv = Σᵢ zᵢ(Cpᵢ − R) − RT·Σ(T·F_TT + 2F_T). It returns the temperature
// and the steps taken.
fn (e EnergyFlash) newton(r FlashResult) !(f64, int) {
	y := r.phases[0].state.n
	mut x0 := e.tvs.present.map(r.beta * y[it])
	x0 << r.beta * r.phases[0].state.molar_volume()
	x0 << math.log(r.t)
	system := fn [e] (x []f64, mut f []f64, mut jac []f64) ! {
		e.eval(x, mut f, mut jac)!
	}
	sol := numeric.newton(system, x0, ftol: e.tvs.params.tol, max_iter: e.tvs.params.max_steps)!
	return math.exp(sol.x[x0.len - 1]), sol.iterations
}

fn (e EnergyFlash) eval(x []f64, mut f []f64, mut jac []f64) ! {
	k := e.tvs.present.len
	n := k + 2 // unknowns: nⱽ, Vⱽ, ln T
	nb := k + 1 // the Helmholtz block
	t := math.exp(x[k + 1])
	if !(t >= state_t_min && t <= state_t_max) {
		return error('T left the range searched')
	}
	s := VolumeFlash{
		...e.tvs
		t: t
	}
	mut g := []f64{len: nb}
	mut h := []f64{len: nb * nb}
	s.helmholtz_eval(x[..nb], mut g, mut h)!
	nc := s.z.len
	vv := x[k]
	vl := s.v - vv
	nv := s.split_moles(x[..nb])
	mut nl := []f64{len: nc}
	for i in s.present {
		nl[i] = s.z[i] - nv[i]
	}
	rv := s.m.residual2(t, vv, nv)
	rl := s.m.residual2(t, vl, nl)
	for row in 0 .. nb {
		f[row] = g[row]
		for col in 0 .. nb {
			jac[row * n + col] = h[row * nb + col]
		}
	}
	mut uig, mut cvig := 0.0, 0.0
	for j, i in s.present {
		c := t * (rv.f_tn[i] - rl.f_tn[i])
		jac[j * n + nb] = c
		jac[nb * n + j] = -c
		ig := e.fl.ideal[i]
		uig += s.z[i] * ((ig.dhf or { 0.0 }) + ig.cp.int_cp(thermo.t_ref, t) - thermo.r_gas * t)
		cvig += s.z[i] * (ig.cp.cp(t) - thermo.r_gas)
	}
	c_v := t * (rv.f_tv - rl.f_tv)
	jac[k * n + nb] = c_v
	jac[nb * n + k] = -c_v
	ft := rv.f_t + rl.f_t
	rt := thermo.r_gas * t
	u := uig - rt * t * ft
	f[nb] = (u - e.u) / rt
	heat := cvig - rt * (t * (rv.f_tt + rl.f_tt) + 2.0 * ft)
	jac[nb * n + nb] = heat / thermo.r_gas - f[nb]
}
