module reactors

import flowsheet
import kinetics
import math
import numeric
import thermo

// Cstr is a continuous stirred-tank reactor (D10): its outlet is its
// contents, so with V its volume,
//
//     Fᵢ,in − Fᵢ + V·Sᵢ(T, P, y) = 0,
//
// Sᵢ = Σⱼ sⱼνᵢⱼrⱼ as the PFR's, at a temperature, or where its energy
// balance puts it, F·H out = F·H in + Q. Solved by Newton (numeric.newton),
// in the outlet's flows and, where T moves, T, from the inlet's.
pub struct Cstr {
pub:
	name      string
	k         kinetics.Kinetics
	reactions []int // by index into the kinetics' reactions
	volume    f64   // m³
	loading   f64   // kg of catalyst per m³, for rates per mass
	thermal   Thermal
	t         f64 // K, if isothermal: 0 for the inlet's
	duty      f64 // W, for .duty
	drop      f64 // Pa
}

// CstrParams are a CSTR's specification.
@[params]
pub struct CstrParams {
pub:
	k         kinetics.Kinetics @[required]
	reactions []int             @[required]
	volume    f64               @[required]
	loading   f64
	thermal   Thermal
	t         f64
	duty      f64
	drop      f64
}

// new_cstr is a CSTR with its specification.
pub fn new_cstr(name string, p CstrParams) !Cstr {
	check_thermal(name, p.thermal, p.t, p.duty, p.drop)!
	if !(p.volume > 0) || math.is_inf(p.volume, 1) {
		return spec_error(name, 'a volume of ${p.volume} m³')
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
	return Cstr{
		name:      name
		k:         p.k
		reactions: p.reactions
		volume:    p.volume
		loading:   p.loading
		thermal:   p.thermal
		t:         p.t
		duty:      p.duty
		drop:      p.drop
	}
}

// inlet_ports is in.
pub fn (c Cstr) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (c Cstr) outlet_ports() []string {
	return ['out']
}

// CstrSystem is a CSTR's equations for Newton, scaled: the balances over
// the inlet's total flow, and the energy balance over F·RT. Each evaluation
// makes its own workspace, numeric.System's receiver being immutable; a
// CSTR's Newton takes a handful.
struct CstrSystem {
	c       Cstr
	fl      thermo.Fluid
	inlet   flowsheet.Stream
	p       f64
	t_fixed f64 // K, where T doesn't move
	target  f64 // W: the outlet's F·H, where T moves
	scale   f64 // mol/s
	src     Source
	// fraction is the share of the volume solved for: 1 but in the
	// continuation that solve uses when Newton from the inlet fails.
	fraction f64 = 1.0
}

fn (sys CstrSystem) eval(x []f64, mut f []f64, mut jac []f64) ! {
	n := sys.inlet.flows.len
	has_t := x.len > n
	t := if has_t { x[n] } else { sys.t_fixed }
	if !(t > 0) {
		return error('no state at ${t} K')
	}
	flows := x[..n]
	for v in flows {
		if v < 0 {
			return error('a flow below zero')
		}
	}
	mut src := sys.src
	st := src.state(t, sys.p, flows)!
	mut s := []f64{len: n}
	src.at(st, mut s)!
	mut total := 0.0
	for v in flows {
		total += v
	}
	mut dflows := []f64{len: n * n}
	mut dt := []f64{len: n}
	mut dp := []f64{len: n}
	src.derivs(st, total, mut dflows, mut dt, mut dp)
	m := x.len
	v := sys.c.volume * sys.fraction
	for i in 0 .. n {
		f[i] = (sys.inlet.flows[i] - flows[i] + v * s[i]) / sys.scale
		for k in 0 .. n {
			jac[i * m + k] = (v * dflows[i * n + k] - if i == k { 1.0 } else { 0.0 }) / sys.scale
		}
		if has_t {
			jac[i * m + n] = v * dt[i] / sys.scale
		}
	}
	if has_t {
		// F·H(T, P, F) − target, and its derivatives: each species' partial
		// molar enthalpy, and Cp.
		props := sys.fl.props(src.st)
		lp := thermo.ln_phi_derivs(sys.fl.model, src.st)
		e := sys.scale * thermo.r_gas * sys.inlet.t
		f[n] = (props.h - sys.target) / e
		for k in 0 .. n {
			g := sys.fl.ideal[k]
			dhf := g.dhf or { 0.0 }
			hbar := dhf + g.cp.int_cp(thermo.t_ref, t) - thermo.r_gas * t * t * lp.dt[k]
			jac[n * m + k] = hbar / e
		}
		jac[n * m + n] = props.cp / e
	}
}

// continuation_start is the share of the volume continuation starts from,
// where the outlet is nearly the inlet; first_step and least_step are its
// steps in ln V, which double on success and halve on failure.
const continuation_start = 1e-6
const first_step = math.ln10 / 4
const least_step = 1e-4

// solve finds the CSTR's outlet by Newton from x0, the inlet. Where that
// fails, it solves by continuation in the volume: from a millionth of it,
// where the inlet is nearly the answer, up to the whole, each solve starting
// from the last, the step in ln V doubled when Newton converges and halved
// when it fails. Fast rates near equilibrium make the balances steep enough
// that Newton from the inlet can't take a step: the owner's reforming at
// 1000 K (M9 step 10). The last solve is the whole volume's, to the same
// tolerance.
fn (c Cstr) solve(sys CstrSystem, x0 []f64) ![]f64 {
	direct := numeric.newton(sys, x0, ftol: 1e-12, xtol: 1e-14) or { numeric.NewtonSolution{} }
	if direct.x.len > 0 {
		return direct.x
	}
	mut x := x0.clone()
	mut solved := math.log(continuation_start) - first_step // ln of the share solved
	mut step := first_step
	mut first := true
	for {
		at := if first { math.log(continuation_start) } else { math.min(0.0, solved + step) }
		part := CstrSystem{
			...sys
			fraction: math.exp(at)
		}
		if sol := numeric.newton(part, x, ftol: 1e-12, xtol: 1e-14) {
			x = sol.x.clone()
			solved, first = at, false
			if at == 0.0 {
				return x
			}
			step = math.min(2 * step, math.ln10)
		} else {
			step /= 2
			if step < least_step {
				return error('Newton failed from the inlet, and continuation in the volume stalled at ${kinetics.sig(math.exp(solved))} of it: ${err.msg()}')
			}
		}
	}
	return x
}

// calculate solves the CSTR by Newton from its inlet, or by continuation in
// its volume where that fails, and finds its outlet at the T it reaches.
pub fn (c Cstr) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	n := c.k.species.len
	inlet := one_inlet(c.name, fl, inlets, n)!
	p := pressure_out(c.name, inlet, c.drop)!
	has_t := c.thermal != .isothermal
	if has_t {
		for i in c.reactions {
			for j in c.k.reactions[i].species {
				if fl.ideal[j].dhf == none {
					return spec_error(c.name, '${c.k.species[j].name} has no ΔH°f, which the energy balance needs')
				}
			}
		}
	}
	sys := CstrSystem{
		c:       c
		fl:      fl
		inlet:   inlet
		p:       p
		t_fixed: if c.t > 0 { c.t } else { inlet.t }
		target:  enthalpy_flow(inlet) + if c.thermal == .duty { c.duty } else { 0.0 }
		scale:   inlet.total()
		src:     new_source(c.name, c.k, c.reactions, c.loading, fl.model)!
	}
	mut x0 := inlet.flows.clone()
	if has_t {
		x0 << inlet.t
	}
	x := c.solve(sys, x0) or { return failed(c.name, err.msg()) }
	flows := x[..n].clone()
	t := if has_t { x[n] } else { sys.t_fixed }
	outlet := at_tp(c.name, fl, flows, t, p)!
	one_phase(c.name, outlet)!
	energy := match c.thermal {
		.isothermal {
			[flowsheet.EnergyStream{'duty', enthalpy_flow(outlet) - enthalpy_flow(inlet)}]
		}
		.adiabatic { []flowsheet.EnergyStream{} }
		.duty { [flowsheet.EnergyStream{'duty', c.duty}] }
	}
	return flowsheet.Outcome{
		outlets:   [outlet]
		energy:    energy
		conserved: element_rows(c.k.species)
	}
}
