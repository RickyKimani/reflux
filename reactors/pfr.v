module reactors

import flowsheet
import kinetics
import math
import numeric
import thermo

// PressureDrop is how a PFR's pressure falls along it.
pub enum PressureDrop {
	no_drop // isobaric
	given   // a given drop, spread evenly over the volume
	ergun   // a packed bed's, by Ergun's equation
}

// Integrator is which of numeric's integrators a PFR runs on (D8).
pub enum Integrator {
	rodas4 // stiff: Rosenbrock, on the rates' own Jacobian
	dopri5 // not: Dormand–Prince, which refuses a stiff problem
}

// Catalyst is a packed bed's: its loading, kg of catalyst per m³ of
// reactor, which a rate per mass of catalyst is multiplied by, and for
// Ergun's drop its void fraction and its particles' diameter, m. The void
// fraction enters Ergun alone, as in DWSIM (step 0).
pub struct Catalyst {
pub:
	loading  f64
	void     f64
	diameter f64
}

// Pfr is a plug-flow reactor or packed bed (D10): the flows of the species
// along its volume V,
//
//     dFᵢ/dV = Sᵢ = Σⱼ sⱼ νᵢⱼ rⱼ(T, P, y),
//
// isothermal, adiabatic or with a duty spread over its volume, isobaric or
// with a given drop or Ergun's,
//
//     dP/dV = −[150 μ(1 − ε)²u/(ε³d²) + 1.75 ρ(1 − ε)u²/(ε³d)]/A,
//
// u the gas's superficial velocity and A = V/L the cross-section. Where T
// moves, H·F rises by the duty per volume q:
//
//     dT/dV = (q − Σᵢ h̄ᵢSᵢ + Cp·μ_JT·dP/dV)/Cp,
//
// h̄ᵢ each species' partial molar enthalpy and Cp the flow's, from the
// fluid's model. One phase, gas: an outlet that would split is refused.
pub struct Pfr {
pub:
	name          string
	k             kinetics.Kinetics
	reactions     []int // by index into the kinetics' reactions
	volume        f64   // m³
	length        f64   // m: 0 for a reactor whose length isn't given
	catalyst      Catalyst
	thermal       Thermal
	t             f64 // K, if isothermal: 0 for the inlet's
	duty          f64 // W, for .duty
	drop          PressureDrop
	pressure_drop f64 // Pa, for .given
	viscosity     ?Viscosity
	integrator    Integrator
	rtol          f64
}

// PfrParams are a PFR's specification.
@[params]
pub struct PfrParams {
pub:
	k             kinetics.Kinetics @[required]
	reactions     []int             @[required]
	volume        f64               @[required]
	length        f64
	catalyst      Catalyst
	thermal       Thermal
	t             f64
	duty          f64
	drop          PressureDrop
	pressure_drop f64
	viscosity     ?Viscosity
	integrator    Integrator
	rtol          f64 = 1e-8
}

// new_pfr is a PFR with its specification.
pub fn new_pfr(name string, p PfrParams) !Pfr {
	check_thermal(name, p.thermal, p.t, p.duty, p.pressure_drop)!
	if !(p.volume > 0) || math.is_inf(p.volume, 1) {
		return spec_error(name, 'a volume of ${p.volume} m³')
	}
	if !(p.length >= 0) || math.is_inf(p.length, 1) {
		return spec_error(name, 'a length of ${p.length} m')
	}
	if !(p.rtol > 1e-14 && p.rtol < 1) {
		return spec_error(name, 'a tolerance of ${p.rtol}')
	}
	if p.reactions.len == 0 {
		return spec_error(name, 'no reactions')
	}
	for i in p.reactions {
		if i < 0 || i >= p.k.reactions.len {
			return spec_error(name, 'no reaction ${i}')
		}
		if i !in p.k.rated {
			return spec_error(name, '${p.k.reactions[i].name} has no rate law')
		}
	}
	if p.drop == .ergun {
		c := p.catalyst
		if !(p.length > 0) {
			return spec_error(name, "Ergun's pressure drop needs the reactor's length")
		}
		if !(c.void > 0 && c.void < 1) {
			return spec_error(name, 'a void fraction of ${c.void}')
		}
		if !(c.diameter > 0) || math.is_inf(c.diameter, 1) {
			return spec_error(name, "a particle's diameter of ${c.diameter} m")
		}
		if p.viscosity == none {
			return spec_error(name, "Ergun's pressure drop needs the gas's viscosity")
		}
	}
	return Pfr{
		name:          name
		k:             p.k
		reactions:     p.reactions
		volume:        p.volume
		length:        p.length
		catalyst:      p.catalyst
		thermal:       p.thermal
		t:             p.t
		duty:          p.duty
		drop:          p.drop
		pressure_drop: p.pressure_drop
		viscosity:     p.viscosity
		integrator:    p.integrator
		rtol:          p.rtol
	}
}

// inlet_ports is in.
pub fn (r Pfr) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (r Pfr) outlet_ports() []string {
	return ['out']
}

// PfrOde is a PFR's equations, for the integrators: the state is the flows,
// then T if it moves, then P under Ergun's drop. A given drop's P is a line
// in V.
struct PfrOde {
	r     Pfr
	fl    thermo.Fluid
	n     int
	has_t bool
	has_p bool
	t_iso f64 // K, where T doesn't move
	p_in  f64 // Pa
	area  f64 // m², for Ergun
	q     f64 // W/m³, the duty per volume
	mw    []f64
	ideal []thermo.IdealGas
mut:
	src    Source
	s      []f64
	dflows []f64
	dt     []f64
	dp     []f64
	f1     []f64
	f2     []f64
	y1     []f64
}

// tp is T and P at V and the state.
fn (o PfrOde) tp(v f64, y []f64) (f64, f64) {
	t := if o.has_t { y[o.n] } else { o.t_iso }
	p := if o.has_p {
		y[o.n + int(o.has_t)]
	} else if o.r.drop == .given {
		o.p_in - o.r.pressure_drop * v / o.r.volume
	} else {
		o.p_in
	}
	return t, p
}

// ergun is dP/dV at the phase just read: −(viscous + inertial)/A.
fn (mut o PfrOde) ergun(t f64, p f64, flows []f64) !f64 {
	c := o.r.catalyst
	visc := o.r.viscosity or { return error('no viscosity') }
	mu := visc.viscosity(t, p, o.src.y)!
	qv := o.src.volume // m³/s
	mut mass := 0.0
	for i, f in flows {
		mass += math.max(f, 0) * o.mw[i]
	}
	rho := mass / qv
	u := qv / o.area
	e3 := c.void * c.void * c.void
	viscous := 150 * mu * (1 - c.void) * (1 - c.void) * u / (e3 * c.diameter * c.diameter)
	inertial := 1.75 * rho * (1 - c.void) * u * u / (e3 * c.diameter)
	return -(viscous + inertial) / o.area
}

// eval is the PFR's right-hand side at V.
fn (mut o PfrOde) eval(v f64, y []f64, mut f []f64) ! {
	t, p := o.tp(v, y)
	if !(t > 0) || !(p > 0) {
		return error('no state at ${t} K and ${p} Pa')
	}
	flows := y[..o.n]
	st := o.src.state(t, p, flows)!
	o.src.at(st, mut o.s)!
	for i in 0 .. o.n {
		f[i] = o.s[i]
	}
	mut dpdv := 0.0
	if o.has_p {
		dpdv = o.ergun(t, p, flows)!
		f[o.n + int(o.has_t)] = dpdv
	} else if o.r.drop == .given {
		dpdv = -o.r.pressure_drop / o.r.volume
	}
	if o.has_t {
		props := o.fl.props(o.src.st)
		lp := thermo.ln_phi_derivs(o.fl.model, o.src.st)
		mut absorbed := 0.0
		for i in 0 .. o.n {
			if o.s[i] == 0 {
				continue
			}
			g := o.ideal[i]
			dhf := g.dhf or { return error('species ${i} has no ΔH°f, which the energy balance needs') }
			hbar := dhf + g.cp.int_cp(thermo.t_ref, t) - thermo.r_gas * t * t * lp.dt[i]
			absorbed += hbar * o.s[i]
		}
		f[o.n] = (o.q - absorbed + props.cp * props.mu_jt * dpdv) / props.cp
	}
}

// jacobian is the rates' block exactly, from the kinetics (D5); the T and P
// rows, and a given drop's ∂f/∂V, by central differences of eval.
fn (mut o PfrOde) jacobian(v f64, y []f64, mut jac []f64, mut dfdt []f64) ! {
	m := y.len
	n := o.n
	t, p := o.tp(v, y)
	flows := y[..n]
	st := o.src.state(t, p, flows)!
	o.src.at(st, mut o.s)!
	mut total := 0.0
	for x in flows {
		total += math.max(x, 0)
	}
	o.src.derivs(st, total, mut o.dflows, mut o.dt, mut o.dp)
	for i in 0 .. n {
		for k in 0 .. n {
			jac[i * m + k] = o.dflows[i * n + k]
		}
		if o.has_t {
			jac[i * m + n] = o.dt[i]
		}
		if o.has_p {
			jac[i * m + n + int(o.has_t)] = o.dp[i]
		}
	}
	// The rows of T and P by central differences, each column's step a
	// cube root of ε of its scale.
	if o.has_t || o.has_p {
		for k in 0 .. m {
			h := 6e-6 * math.max(math.abs(y[k]), if k < n { 1e-6 * total } else { 1.0 })
			for i in 0 .. m {
				o.y1[i] = y[i]
			}
			o.y1[k] = y[k] + h
			o.eval(v, o.y1, mut o.f1)!
			o.y1[k] = y[k] - h
			o.eval(v, o.y1, mut o.f2)!
			for i in n .. m {
				jac[i * m + k] = (o.f1[i] - o.f2[i]) / (2 * h)
			}
		}
	}
	for i in 0 .. m {
		dfdt[i] = 0
	}
	if o.r.drop == .given {
		// Only P depends on V, through its line: ∂f/∂V = ∂f/∂P·dP/dV.
		dpdv := -o.r.pressure_drop / o.r.volume
		for i in 0 .. n {
			dfdt[i] = o.dp[i] * dpdv
		}
		if o.has_t {
			h := 6e-6 * math.max(math.abs(v), o.r.volume)
			o.eval(v + h, y, mut o.f1)!
			o.eval(v - h, y, mut o.f2)!
			dfdt[n] = (o.f1[n] - o.f2[n]) / (2 * h)
		}
	}
}

// calculate integrates the PFR from its inlet and finds its outlet: at the
// integrated T and P where it's isothermal, and otherwise by a PH flash on
// the energy balance, held to the integrated T. The profile is the inlet
// and every point the integrator stepped to.
pub fn (r Pfr) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	n := r.k.species.len
	inlet := one_inlet(r.name, fl, inlets, n)!
	src := new_source(r.name, r.k, r.reactions, r.catalyst.loading, fl.model)!
	has_t := r.thermal != .isothermal
	has_p := r.drop == .ergun
	m := n + int(has_t) + int(has_p)
	mut o := PfrOde{
		r:      r
		fl:     fl
		n:      n
		has_t:  has_t
		has_p:  has_p
		t_iso:  if r.t > 0 { r.t } else { inlet.t }
		p_in:   inlet.p
		area:   if r.length > 0 { r.volume / r.length } else { 0.0 }
		q:      if r.thermal == .duty { r.duty / r.volume } else { 0.0 }
		mw:     fl.model.components().map(it.mw)
		ideal:  fl.ideal
		src:    src
		s:      []f64{len: n}
		dflows: []f64{len: n * n}
		dt:     []f64{len: n}
		dp:     []f64{len: n}
		f1:     []f64{len: m}
		f2:     []f64{len: m}
		y1:     []f64{len: m}
	}
	if r.drop == .given && !(inlet.p - r.pressure_drop > 0) {
		return spec_error(r.name, 'a drop of ${r.pressure_drop} Pa from ${inlet.p} Pa leaves no pressure')
	}
	mut y0 := inlet.flows.clone()
	total := inlet.total()
	mut atols := []f64{len: n, init: 1e-6 * r.rtol * total}
	if has_t {
		y0 << inlet.t
		atols << r.rtol * inlet.t
	}
	if has_p {
		y0 << inlet.p
		atols << r.rtol * inlet.p
	}
	params := numeric.OdeParams{
		rtol:   r.rtol
		atols:  atols
		record: true
	}
	sol := if r.integrator == .rodas4 {
		numeric.rodas4(mut o, 0, y0, r.volume, params) or {
			return failed(r.name, 'its integration failed: ${err.msg()}')
		}
	} else {
		numeric.dopri5(mut o, 0, y0, r.volume, params) or {
			return failed(r.name, 'its integration failed: ${err.msg()}')
		}
	}
	end := sol.y
	t_end, p_end := o.tp(r.volume, end)
	flows := end[..n].map(math.max(it, 0))
	outlet := match r.thermal {
		.isothermal {
			at_tp(r.name, fl, flows, t_end, p_end)!
		}
		.adiabatic, .duty {
			h := (enthalpy_flow(inlet) + if r.thermal == .duty { r.duty } else { 0.0 }) / numeric.sum(flows)
			s := flowsheet.feed_ph(fl, flows, p_end, h) or { return named(r.name, err) }
			// The flash's T and the integrated one agree to the integration's
			// tolerance, or the energy equation and the fluid disagree.
			if math.abs(s.t - t_end) > 100 * r.rtol * t_end + 1e-9 * t_end {
				return failed(r.name, 'its outlet is at ${s.t} K by its energy balance and ${t_end} K by its integration')
			}
			s
		}
	}
	one_phase(r.name, outlet)!
	energy := match r.thermal {
		.isothermal {
			[flowsheet.EnergyStream{'duty', enthalpy_flow(outlet) - enthalpy_flow(inlet)}]
		}
		.adiabatic { []flowsheet.EnergyStream{} }
		.duty { [flowsheet.EnergyStream{'duty', r.duty}] }
	}
	mut points := []flowsheet.ProfilePoint{}
	for k, v in sol.ts {
		yk := sol.ys[k]
		tk, pk := o.tp(v, yk)
		points << flowsheet.ProfilePoint{
			volume: v
			length: if r.length > 0 { v / r.volume * r.length } else { 0.0 }
			t:      tk
			p:      pk
			flows:  yk[..n].clone()
		}
	}
	return flowsheet.Outcome{
		outlets:   [outlet]
		energy:    energy
		conserved: element_rows(r.k.species)
		profile:   flowsheet.Profile{
			points: points
		}
	}
}
