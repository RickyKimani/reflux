module water

import math
import thermo
import transport
import water.if97

// IAPWS-IF97 as a thermo.PhaseModel for water (M5b's D12, D13, D14), so that
// the flash, and a flowsheet, give the steam tables' states.
//
// - A phase at (T, P): the stable one is the state R7-97's region map gives,
//   region 4 deciding liquid or vapour below the critical point; the liquid
//   and the vapour are each side's equation across the line, as far as the
//   flash's paths for one compound ask (if97.liquid_pt, if97.vapour_pt).
// - Its ideal gas is region 2's ideal-gas part γ°, nine powers of
//   τ = 540/T, a thermo.Polynomial Cp (D14), and every region's residual is
//   against it, so that a stream's h and s are IF97's plus constants. IF97's
//   gas constant is carried as IAPWS-95's is (M5a's D4): the residual
//   entropy holds (R − R₉₇)·ln(P/P°), so that ln φ's P-derivative is
//   (V − RT/P)/RT with reflux's R.
// - Its saturation line is region 4's (D13, thermo.Saturates), not where its
//   regions' Gibbs energies meet, which is up to 3.5 mK away.
// - Outside IF97's range, 273.15 K to 1073.15 K to 100 MPa and to 2273.15 K
//   to 50 MPa, it has no state (thermo.NoState), so a flash's search steps
//   back into it (M5a's at_or_closer).
//
// It has no Helmholtz energy in (T, V) for regions 1, 2 and 5, so the TV
// and UV flashes refuse it.
pub struct If97Phases {
	comps []thermo.Component
}

// r97 is IF97's gas constant per mole, with IAPWS's molar mass.
pub const r97 = if97.r * molar_mass

// if97_phases makes the model: water, with IAPWS's critical constants.
pub fn if97_phases() If97Phases {
	return If97Phases{
		comps: iapws95().comps
	}
}

// components implements thermo.PhaseModel.
pub fn (m If97Phases) components() []thermo.Component {
	return m.comps
}

// ideal_cp_if97 is region 2's ideal-gas Cp per mole, −R₉₇τ²γ°_ττ =
// −R₉₇Σ nᵢJᵢ(Jᵢ − 1)(540/T)^Jᵢ, as a polynomial in T (D14): its terms with J
// of 0 and 1 have none.
pub fn ideal_cp_if97() thermo.CpForm {
	mut terms := []thermo.PowerTerm{}
	for k, j in if97.ideal2_j {
		if j == 0 || j == 1 {
			continue
		}
		c := -r97 * if97.ideal2_n[k] * f64(j * (j - 1)) * math.pow(540.0, f64(j))
		terms << thermo.PowerTerm{c, -j}
	}
	return thermo.Polynomial{
		terms: terms
	}
}

// if97_fluid is IF97 water as a thermo.Fluid, its H and S on reflux's
// reference through water's ΔH°f and S°, as IAPWS-95's (water.fluid).
pub fn if97_fluid(dhf ?f64, s0 ?f64) !thermo.Fluid {
	return thermo.new_phase_fluid(if97_phases(), [
		thermo.IdealGas{
			cp:  ideal_cp_if97()
			dhf: dhf
			s0:  s0
		},
	])
}

fn one(n []f64) !f64 {
	if n.len != 1 {
		return error('${n.len} mole numbers for water alone')
	}
	if !(n[0] > 0) || math.is_inf(n[0], 1) {
		return error('${n[0]} mol of water')
	}
	return n[0]
}

// edge_hair is how far past an edge of IF97's range in T a phase is the
// edge's own equation's: a flash's search can probe that far past an answer
// on the edge (Brent's tolerance step went 4.8e-7 K past 1073.15 K at
// 100 MPa), as IF97's backward functions take a value within its rounding
// of an edge as on it.
const edge_hair = 1e-9

// phase implements thermo.PhaseModel.
pub fn (m If97Phases) phase(t f64, p f64, n []f64, root thermo.Root) !thermo.State {
	n0 := one(n)!
	if p > 0 && p <= if97.p_max {
		if t > if97.t5 && p > if97.p5_max && t <= if97.t5 * (1 + edge_hair) {
			return state_of(if97.region2(p, t), t, p, n0)
		}
		if t > if97.t_max && t <= if97.t_max * (1 + edge_hair) {
			return state_of(if97.region5(p, t), t, p, n0)
		}
		if t < if97.t_min && t >= if97.t_min * (1 - edge_hair) {
			region := if p >= if97.psat97(if97.t_min) { 1 } else { 2 }
			return state_of(if region == 1 { if97.region1(p, t) } else { if97.region2(p, t) }, t,
				p, n0)
		}
	}
	st := match root {
		.stable { if97.state_pt(p, t)! }
		.liquid { if97.liquid_pt(p, t)! }
		.vapor { if97.vapour_pt(p, t)! }
	}
	return state_of(st, t, p, n0)
}

// state_of is an IF97 state as a phase of n0 mol: its volume, Z with
// reflux's R, and ln φ = Gʳ/RT against region 2's ideal gas (D14).
fn state_of(st if97.State, t f64, p f64, n0 f64) thermo.State {
	ideal := if97.region2_ideal(p, t)
	gr := molar_mass * ((st.h - t * st.s) - (ideal.h - t * ideal.s)) - t * (thermo.r_gas - r97) * math.log(p / thermo.p_ref)
	v := n0 * molar_mass * st.v
	return thermo.State{
		t:      t
		p:      p
		n:      [n0]
		v:      v
		z:      p * v / (n0 * thermo.r_gas * t)
		ln_phi: [gr / (thermo.r_gas * t)]
	}
}

// state_at is the IF97 state of a phase of water at T, P and density ρ
// (kg/m³), the region the phase came from: by ρ against ρc below Tc (region
// 1 or 3's liquid side, region 2 or 3's vapour side by B23), by the region
// map above it. Region 3's is at ρ, the rest at (P, T).
pub fn (m If97Phases) state_at(t f64, p f64, rho f64) if97.State {
	if t > if97.t5 {
		return if97.region5(p, t)
	}
	if t <= if97.t13 {
		return if rho > if97.rhoc { if97.region1(p, t) } else { if97.region2(p, t) }
	}
	if t < if97.tc {
		if rho > if97.rhoc || p > if97.p23(t) {
			return if97.region3(rho, t)
		}
		return if97.region2(p, t)
	}
	if t <= if97.t23_max && p > if97.p23(t) {
		return if97.region3(rho, t)
	}
	return if97.region2(p, t)
}

fn (m If97Phases) state_of_phase(st thermo.State) if97.State {
	return m.state_at(st.t, st.p, st.n[0] * molar_mass / st.v)
}

// ln_phi_derivs implements thermo.PhaseModel: ∂ln φ/∂T = −Hʳ/RT²,
// ∂ln φ/∂P = (V − RT/P)/RT, and n·∂ln φ/∂n = 0 for one compound.
pub fn (m If97Phases) ln_phi_derivs(st thermo.State) !thermo.LnPhiDerivs {
	one(st.n)!
	rs := m.state_of_phase(st)
	ideal := if97.region2_ideal(st.p, st.t)
	hr := molar_mass * (rs.h - ideal.h)
	rt := thermo.r_gas * st.t
	return thermo.LnPhiDerivs{
		dt: [-hr / (rt * st.t)]
		dp: [(st.v / st.n[0] - rt / st.p) / rt]
		dn: [0.0]
	}
}

// phase_label implements thermo.PhaseModel: region 1 liquid, regions 2 and
// 5 vapour, region 3 by ρc (D12).
pub fn (m If97Phases) phase_label(st thermo.State) thermo.PhaseLabel {
	rs := m.state_of_phase(st)
	return match rs.region {
		1 { .liquid }
		3 {
			if rs.rho() > if97.rhoc { .liquid } else { .vapor }
		}
		else { .vapor }
	}
}

// phase_residual implements thermo.PhaseModel: each property less region
// 2's ideal gas's, extensive (D14). Cvʳ is Cv less the ideal gas's Cp − R,
// as thermo.Fluid adds it back.
pub fn (m If97Phases) phase_residual(st thermo.State) thermo.PhaseResidual {
	n0 := st.n[0]
	rs := m.state_of_phase(st)
	ideal := if97.region2_ideal(st.p, st.t)
	return thermo.PhaseResidual{
		h:     n0 * molar_mass * (rs.h - ideal.h)
		s:     n0 * (molar_mass * (rs.s - ideal.s) + (thermo.r_gas - r97) * math.log(st.p / thermo.p_ref))
		cp:    n0 * molar_mass * (rs.cp - ideal.cp)
		cv:    n0 * (molar_mass * (rs.cv - ideal.cp) + thermo.r_gas)
		dv_dt: n0 * molar_mass * rs.dvdt
		dp_dv: -1.0 / (n0 * molar_mass * rs.v * rs.v * rs.drhodp)
	}
}

// saturation_p implements thermo.Saturates: region 4's Pˢᵃᵗ(T), below Tc.
pub fn (m If97Phases) saturation_p(t f64) !f64 {
	if !(t >= if97.t_min && t < if97.tc) {
		return error("IF97's saturation line runs from ${if97.t_min} K to its critical point, ${if97.tc} K: not ${t} K")
	}
	return if97.psat97(t)
}

// saturation_t implements thermo.Saturates: region 4's Tˢᵃᵗ(P), below Pc.
pub fn (m If97Phases) saturation_t(p f64) !f64 {
	if !(p >= if97.psat97(if97.t_min) && p < if97.pc) {
		return error("IF97's saturation line runs from ${if97.psat97(if97.t_min)} Pa to its critical point, ${if97.pc} Pa: not ${p} Pa")
	}
	return if97.tsat97(p)
}

// saturated implements thermo.Saturates: the liquid and the vapour at (T, P)
// on the line, regions 1 and 2 to 623.15 K, region 3's two sides above.
pub fn (m If97Phases) saturated(t f64, p f64, n []f64) ![]thermo.State {
	n0 := one(n)!
	liq, vap := if97.saturated_at(p, t)!
	return [state_of(liq, t, p, n0), state_of(vap, t, p, n0)]
}

// saturation_gap implements thermo.Saturates: how far apart the two
// phases' ln φ may be on the line, R7-97's Prague value for the saturation
// pressure's consistency, 0.05%, which a gap of Δ ln φ moves Pˢᵃᵗ by (the
// largest is 5.4e-5, M5b step 0).
pub fn (m If97Phases) saturation_gap() f64 {
	return 5e-4
}

// no_state implements thermo.NoState: outside IF97's range.
pub fn (m If97Phases) no_state(t f64, p f64) bool {
	return !(t >= if97.t_min && t <= if97.t_max) || !(p > 0 && p <= if97.p_max)
		|| (t > if97.t5 && p > if97.p5_max)
}

// TransportIf97 is a stream of IF97 water's transport (M5b's D10): R12-08's
// and R15-11's industrial forms on the phase's IF97 state, and R1-76's
// surface tension for a liquid.
pub struct TransportIf97 {}

// phase implements transport.Model.
pub fn (_ TransportIf97) phase(label thermo.PhaseLabel, t f64, p f64, x []f64, v f64) transport.PhaseTransport {
	if x.len != 1 || !(v > 0) || math.is_inf(v, 1) {
		why := if x.len != 1 {
			'${x.len} mole fractions for water alone'
		} else {
			'V = ${v} m³/mol'
		}
		return transport.PhaseTransport{
			viscosity:       transport.Result{
				why: why
			}
			conductivity:    transport.Result{
				why: why
			}
			surface_tension: transport.Result{
				why: why
			}
		}
	}
	rho := molar_mass / v
	st := if97_phases().state_at(t, p, rho)
	mu := if m := viscosity_industrial(t, rho) {
		transport.Result{
			value: m
		}
	} else {
		transport.Result{
			why: err.msg()
		}
	}
	lambda := if c := conductivity_industrial(st) {
		transport.Result{
			value: c.lambda
		}
	} else {
		transport.Result{
			why: err.msg()
		}
	}
	mut sigma := transport.Result{
		why: 'a gas has no surface tension'
	}
	if label != .vapor {
		sigma = if s := surface_tension(t) {
			transport.Result{
				value: s
			}
		} else {
			transport.Result{
				why: err.msg()
			}
		}
	}
	return transport.PhaseTransport{
		viscosity:       mu
		conductivity:    lambda
		surface_tension: sigma
	}
}

// methods implements transport.Model: M5a's industrial forms (D10).
pub fn (_ TransportIf97) methods() []transport.MethodUsed {
	return [
		transport.MethodUsed{'viscosity', 'any', "IAPWS R12-08's industrial form, at the phase's T and IF97 density"},
		transport.MethodUsed{'thermal_conductivity', 'any', "IAPWS R15-11's industrial form, with its critical enhancement on the phase's IF97 state"},
		transport.MethodUsed{'surface_tension', 'liquid', 'IAPWS R1-76(2014), at T'},
	]
}
