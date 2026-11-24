// The flash on a phase model with no F(T, V, n) at all (M6's D3). Raoult's
// law is the least γ–φ fluid there is: an ideal solution of incompressible
// liquids under an ideal-gas vapour, with ln Pˢᵃᵗ = a − b/T. Every answer
// has a closed form, so nothing here leans on a cubic, and a flash that
// still reached for one would fail or be wrong.
import flash
import math
import thermo { Component, Root }

const benzene = Component{
	name:  'benzene'
	tc:    562.05
	pc:    48.95e5
	omega: 0.21
	mw:    0.078112
}

const toluene = Component{
	name:  'toluene'
	tc:    591.75
	pc:    41.08e5
	omega: 0.264
	mw:    0.092138
}

const atm = 101325.0

// Each one's slope b, K, boiling point at 1 atm, liquid volume, m³/mol, and
// the ideal gas's constant Cp, ΔH°f and S°.
const slopes = [3695.0, 4000.0]
const boiling = [353.24, 383.78]
const liquid_volumes = [89.4e-6, 106.9e-6]
const heat_capacities = [82.4, 103.7]
const formation = [82.9e3, 50.2e3]
const entropies = [269.3, 320.7]

// Raoult is the model. cold, if set, is the T below which it has no liquid
// residual to give (NaN, as PhaseModel allows); hot, if set, the T above
// which it has no phase at all; flat, if set, means no derivatives.
struct Raoult {
	comps []Component
	a     []f64
	b     []f64
	vl    []f64
	cold  f64
	hot   f64
	flat  bool
}

@[params]
struct Broken {
	cold f64
	hot  f64
	flat bool
}

fn raoult(which []int, how Broken) Raoult {
	both := [benzene, toluene]
	return Raoult{
		comps: which.map(both[it])
		a:     which.map(math.log(atm) + slopes[it] / boiling[it])
		b:     which.map(slopes[it])
		vl:    which.map(liquid_volumes[it])
		cold:  how.cold
		hot:   how.hot
		flat:  how.flat
	}
}

fn fluid(which []int, how Broken) !thermo.Fluid {
	return thermo.new_phase_fluid(raoult(which, how), which.map(thermo.IdealGas{
		cp:  thermo.Polynomial{
			terms: [thermo.PowerTerm{heat_capacities[it], 0}]
		}
		dhf: formation[it]
		s0:  entropies[it]
	}))
}

fn (m Raoult) psat(i int, t f64) f64 {
	return math.exp(m.a[i] - m.b[i] / t)
}

fn (m Raoult) components() []Component {
	return m.comps
}

fn (m Raoult) phase(t f64, p f64, n []f64, root Root) !thermo.State {
	if m.hot > 0.0 && t > m.hot {
		return error("Raoult's vapour pressures end at ${m.hot} K")
	}
	mut total, mut g, mut v := 0.0, 0.0, 0.0
	mut lnphi := []f64{len: n.len}
	for i, ni in n {
		lnphi[i] = m.a[i] - m.b[i] / t - math.log(p)
		total += ni
		g += ni * lnphi[i]
		v += ni * m.vl[i]
	}
	liquid := root == .liquid || (root == .stable && g < 0.0)
	if !liquid {
		v = total * thermo.r_gas * t / p
		lnphi = []f64{len: n.len}
	}
	return thermo.State{
		t:      t
		p:      p
		n:      n.clone()
		v:      v
		z:      p * v / (total * thermo.r_gas * t)
		ln_phi: lnphi
	}
}

// is_liquid tells the model's two phases apart, below 20 bar.
fn is_liquid(st thermo.State) bool {
	return st.z < 0.5
}

fn (m Raoult) ln_phi_derivs(st thermo.State) !thermo.LnPhiDerivs {
	if m.flat {
		return error('Raoult has no derivatives here')
	}
	nc := st.n.len
	mut dt, mut dp := []f64{len: nc}, []f64{len: nc}
	if is_liquid(st) {
		for i in 0 .. nc {
			dt[i] = m.b[i] / (st.t * st.t)
			dp[i] = -1.0 / st.p
		}
	}
	return thermo.LnPhiDerivs{
		dt: dt
		dp: dp
		dn: []f64{len: nc * nc}
	}
}

fn (m Raoult) phase_label(st thermo.State) thermo.PhaseLabel {
	return if is_liquid(st) { .liquid } else { .vapor }
}

fn (m Raoult) phase_residual(st thermo.State) thermo.PhaseResidual {
	mut total, mut h, mut g := 0.0, 0.0, 0.0
	for i, ni in st.n {
		total += ni
		h -= ni * thermo.r_gas * m.b[i]
		g += ni * thermo.r_gas * st.t * st.ln_phi[i]
	}
	r := thermo.r_gas
	if !is_liquid(st) {
		return thermo.PhaseResidual{
			dv_dt: total * r / st.p
			dp_dv: -st.p / st.v
		}
	}
	if st.t < m.cold {
		none_ := math.nan()
		return thermo.PhaseResidual{none_, none_, none_, none_, none_, none_}
	}
	// An incompressible liquid's Cv is its Cp.
	return thermo.PhaseResidual{
		h:     h
		s:     (h - g) / st.t
		cv:    total * r
		dv_dt: 0.0
		dp_dv: math.inf(-1)
	}
}

// Split is Rachford–Rice's closed form for two components of constant K.
struct Split {
	beta f64
	x    []f64
	y    []f64
}

fn split_of(m Raoult, t f64, p f64, z []f64) Split {
	k := [m.psat(0, t) / p, m.psat(1, t) / p]
	c1, c2 := k[0] - 1.0, k[1] - 1.0
	beta := -(z[0] * c1 + z[1] * c2) / (c1 * c2)
	x := [z[0] / (1.0 + beta * c1), z[1] / (1.0 + beta * c2)]
	return Split{beta, x, [k[0] * x[0], k[1] * x[1]]}
}

// whole is the closed form's molar H and S of one mole at (T, P).
fn whole(m Raoult, t f64, p f64, z []f64) (f64, f64) {
	bubble := z[0] * m.psat(0, t) + z[1] * m.psat(1, t)
	dew := 1.0 / (z[0] / m.psat(0, t) + z[1] / m.psat(1, t))
	mut beta := if p >= bubble { 0.0 } else { 1.0 }
	mut x, mut y := z.clone(), z.clone()
	if p < bubble && p > dew {
		s := split_of(m, t, p, z)
		beta, x, y = s.beta, s.x.clone(), s.y.clone()
	}
	r := thermo.r_gas
	mut h, mut s := 0.0, 0.0
	for i in 0 .. 2 {
		hig := formation[i] + heat_capacities[i] * (t - thermo.t_ref)
		sig := entropies[i] + heat_capacities[i] * math.log(t / thermo.t_ref) - r * math.log(p / thermo.p_ref)
		h += beta * y[i] * hig + (1.0 - beta) * x[i] * (hig - r * m.b[i])
		if beta > 0.0 {
			s += beta * y[i] * (sig - r * math.log(y[i]))
		}
		if beta < 1.0 {
			// The liquid's residual S is −R(a − ln P), whatever T is.
			s += (1.0 - beta) * x[i] * (sig - r * math.log(x[i]) - r * (m.a[i] - math.log(p)))
		}
	}
	return h, s
}

fn test_pt_and_stability_by_raoults_law() {
	m := raoult([0, 1])
	seam := thermo.PhaseModel(m)
	mut splits, mut liquids, mut vapours := 0, 0, 0
	mut worst_beta, mut worst_x := 0.0, 0.0
	for k in 0 .. 26 {
		t := 345.0 + 2.0 * f64(k)
		for p in [0.8e5, atm, 1.3e5] {
			for z1 in [0.1, 0.3, 0.5, 0.7, 0.9] {
				z := [z1, 1.0 - z1]
				at := '${t} K, ${p} Pa, ${z}'
				bubble := z[0] * m.psat(0, t) + z[1] * m.psat(1, t)
				dew := 1.0 / (z[0] / m.psat(0, t) + z[1] / m.psat(1, t))
				r := flash.pt(seam, t, p, z)!
				stable := flash.stability(seam, t, p, z)!.stable
				assert stable == (p >= bubble || p <= dew), at
				if stable {
					assert r.phases.len == 1, at
					if p >= bubble {
						assert r.phases[0].label == .liquid && r.beta == 0.0, at
						liquids++
					} else {
						assert r.phases[0].label == .vapor && r.beta == 1.0, at
						vapours++
					}
					continue
				}
				splits++
				want := split_of(m, t, p, z)
				assert r.phases.len == 2 && r.phases[0].label == .vapor, at
				worst_beta = math.max(worst_beta, math.abs(r.beta - want.beta))
				for i in 0 .. 2 {
					worst_x = math.max(worst_x, math.abs(r.phases[0].state.n[i] - want.y[i]))
					worst_x = math.max(worst_x, math.abs(r.phases[1].state.n[i] - want.x[i]))
				}
			}
		}
	}
	println('${splits} splits, ${liquids} liquids, ${vapours} vapours: β ${worst_beta:.1e}, x and y ${worst_x:.1e}')
	assert splits >= 15 && liquids >= 15 && vapours >= 15
	// Measured, 1.6e-14 and 1.1e-15 under tcc and gcc. K doesn't depend on
	// composition, so substitution's first step is exact, and what's left
	// is Rachford–Rice's rounding over K − 1.
	assert worst_beta <= 1e-13 && worst_x <= 1e-14
}

fn test_bubble_and_dew_points_by_raoults_law() {
	m := raoult([0, 1])
	seam := thermo.PhaseModel(m)
	mut worst := 0.0
	for t in [345.0, 360.0, 375.0] {
		for z1 in [0.2, 0.5, 0.8] {
			z := [z1, 1.0 - z1]
			bubble := z[0] * m.psat(0, t) + z[1] * m.psat(1, t)
			dew := 1.0 / (z[0] / m.psat(0, t) + z[1] / m.psat(1, t))
			worst = math.max(worst, math.abs(flash.bubble_p(seam, t, z)!.p / bubble - 1.0))
			worst = math.max(worst, math.abs(flash.dew_p(seam, t, z)!.p / dew - 1.0))
			worst = math.max(worst, math.abs(flash.bubble_t(seam, bubble, z)!.t / t - 1.0))
			worst = math.max(worst, math.abs(flash.dew_t(seam, dew, z)!.t / t - 1.0))
			mid := 0.5 * (bubble + dew)
			beta := split_of(m, t, mid, z).beta
			worst = math.max(worst, math.abs(flash.tvf(seam, t, beta, z)!.p / mid - 1.0))
			worst = math.max(worst, math.abs(flash.pvf(seam, mid, beta, z)!.t / t - 1.0))
		}
	}
	for i in 0 .. 2 {
		one := thermo.PhaseModel(raoult([i]))
		worst = math.max(worst, math.abs(flash.pvf(one, atm, 0.5, [1.0])!.t / boiling[i] - 1.0))
		worst = math.max(worst, math.abs(flash.tvf(one, boiling[i], 0.5, [1.0])!.p / atm - 1.0))
	}
	println('bubble, dew and vapour-fraction points: ${worst:.1e}')
	// Measured, 2.2e-13, and Newton stops at 1e-12 in ln f.
	assert worst <= 2e-12
}

fn test_ph_and_ps_by_raoults_law() {
	m := raoult([0, 1])
	fl := fluid([0, 1])!
	z := [0.5, 0.5]
	mut phases := map[int]int{}
	mut worst := 0.0
	for t in [300.0, 340.0, 362.0, 365.0, 368.0, 371.0, 420.0] {
		h, s := whole(m, t, atm, z)
		by_h := flash.ph(fl, atm, h, z)!
		by_s := flash.ps(fl, atm, s, z)!
		phases[by_h.phases.len]++
		assert by_h.phases.len == by_s.phases.len, '${t} K'
		worst = math.max(worst, math.abs(by_h.t / t - 1.0))
		worst = math.max(worst, math.abs(by_s.t / t - 1.0))
	}
	// A pure component inside its jump at its boiling point.
	one := fluid([0])!
	hv := formation[0] + heat_capacities[0] * (boiling[0] - thermo.t_ref)
	r := flash.ph(one, atm, hv - 0.75 * thermo.r_gas * slopes[0], [1.0])!
	assert r.phases.len == 2 && math.abs(r.beta - 0.25) <= 1e-10, '${r.beta}'
	assert math.abs(r.t / boiling[0] - 1.0) <= 1e-11
	println('PH and PS: T to ${worst:.1e}, ${phases}')
	assert (phases[1] or { 0 }) >= 2 && (phases[2] or { 0 }) >= 2
	// Measured, 1.8e-12. etol is 1e-10 of RT in H, and Cp is 10R here, more
	// in a split, so T is good to 1e-11 of itself.
	assert worst <= 1e-11
}

// Fluid.props takes the speed of sound and the Joule–Thomson coefficient
// from what the seam gives: an ideal gas's √(γRT/M) and 0.
fn test_props_on_a_model_without_f() {
	fl := fluid([0])!
	t := 400.0
	gas := fl.props(fl.model.phase(t, atm, [1.0], .vapor)!)
	cp, r := heat_capacities[0], thermo.r_gas
	assert math.abs(gas.w / math.sqrt(cp / (cp - r) * r * t / benzene.mw) - 1.0) <= 1e-14
	assert math.abs(gas.mu_jt) <= 1e-18
	assert math.abs(gas.cv / (cp - r) - 1.0) <= 1e-14
	liquid := fl.props(fl.model.phase(300.0, atm, [1.0], .liquid)!)
	assert math.abs(liquid.cv / cp - 1.0) <= 1e-14 && math.is_inf(liquid.w, 1)
	assert math.abs(liquid.mu_jt / (-liquid_volumes[0] / cp) - 1.0) <= 1e-14
}

// The UV flash works in F(T, V, n), which this model hasn't: refused, with
// why (D3). flash.tv asks for a Model, so its caller meets helmholtz's error.
fn test_tv_and_uv_refuse_a_model_without_f() {
	fl := fluid([0, 1])!
	if _ := flash.uv(fl, 5e4, 1e-3, [0.5, 0.5]) {
		assert false, 'a UV flash with no F(T, V, n)'
	} else {
		assert err is flash.FlashError, err.msg()
		assert (err as flash.FlashError).kind == .bad_input
		assert err.msg().contains('F(T, V, n)'), err.msg()
	}
	if _ := fl.helmholtz() {
		assert false, 'a Helmholtz energy from Raoult'
	} else {
		assert err is thermo.NoHelmholtzError, err.msg()
	}
}

// What the model refuses reaches the caller as the flash's typed error,
// with the model's reason.
fn test_the_models_refusals_reach_the_caller() {
	hot := thermo.PhaseModel(raoult([0, 1], hot: 360.0))
	if _ := flash.pt(hot, 400.0, atm, [0.5, 0.5]) {
		assert false, 'a flash where the model has no phase'
	} else {
		assert err is flash.FlashError && (err as flash.FlashError).kind == .no_volume_root, err.msg()
		assert err.msg().contains("Raoult's vapour pressures end"), err.msg()
	}
	// No derivatives is a typed refusal too, with the model's reason, which
	// numeric.newton passes on (M6 step 6).
	flat := thermo.PhaseModel(raoult([0, 1], flat: true))
	if _ := flash.bubble_p(flat, 360.0, [0.5, 0.5]) {
		assert false, 'a Newton with no derivatives'
	} else {
		assert err is flash.FlashError, err.msg()
		assert err.msg().contains('Raoult has no derivatives here'), err.msg()
	}
}

// A phase's residual can be NaN, where the model has none (PhaseModel), and
// a NaN energy balance mustn't pass for a met one: these bracket the edge
// where the NaN begins, which Brent closes in on.
fn test_a_residual_that_is_none_is_no_answer() {
	m := raoult([0, 1])
	z := [0.5, 0.5]
	h, s := whole(m, 280.0, atm, z)
	assert math.abs(flash.ph(fluid([0, 1], cold: 270.0)!, atm, h, z)!.t / 280.0 - 1.0) <= 1e-10
	broken := fluid([0, 1], cold: 300.0)!
	if r := flash.ph(broken, atm, h, z, estimate: 340.0) {
		assert false, 'PH answered ${r.t} K, on a NaN enthalpy'
	} else {
		assert err is flash.FlashError, err.msg()
	}
	if r := flash.ps(broken, atm, s, z, estimate: 340.0) {
		assert false, 'PS answered ${r.t} K, on a NaN entropy'
	} else {
		assert err is flash.FlashError, err.msg()
	}
	hl := formation[0] + heat_capacities[0] * (280.0 - thermo.t_ref) - thermo.r_gas * slopes[0]
	// Whether Brent's last point lies on the NaN's side of the edge or the
	// finite one turns on where the edge is, so several.
	for edge in 0 .. 40 {
		cold := 285.0 + 1.37 * f64(edge)
		pure := fluid([0], cold: cold)!
		if r := flash.ph(pure, atm, hl, [1.0]) {
			assert false, 'a pure PH answered ${r.t} K, on a NaN enthalpy below ${cold} K'
		} else {
			assert err is flash.FlashError, err.msg()
		}
	}
}
