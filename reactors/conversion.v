module reactors

import expr
import flowsheet
import kinetics
import math
import numeric
import thermo

// Conversion is one step of a conversion reactor: a reaction, the reactant
// it's measured by, and the share of what enters the step that reacts: a
// number, or an expression in the reactor's T (D10, as amended before step
// 10), in the reactor's context and with its definitions.
pub struct Conversion {
pub:
	reaction   kinetics.Reaction
	key        int // the species, by index
	conversion f64 // 0 to 1
	expression ?expr.Expr
}

// ConversionReactor runs its reactions in order (D10), each converting its
// share of its key reactant as it reaches that step, so a later reaction
// sees what the earlier ones left. The outlet is at a temperature, or where
// the energy balance puts it.
pub struct ConversionReactor {
pub:
	name    string
	species []kinetics.Species
	steps   []Conversion
	thermal Thermal
	t       f64              // K, if isothermal: 0 for the inlet's
	duty    f64              // W, for .duty
	drop    f64              // Pa
	ctx     expr.Context     // what a step's expression is bound in, T its first slot
	defs    expr.Definitions // and the definitions it may use
}

// ConversionParams are a conversion reactor's specification. A step with an
// expression needs the context and definitions it was read in.
@[params]
pub struct ConversionParams {
pub:
	species []kinetics.Species @[required]
	steps   []Conversion       @[required]
	thermal Thermal
	t       f64
	duty    f64
	drop    f64
	ctx     expr.Context
	defs    expr.Definitions
}

// new_conversion_reactor is a conversion reactor with its specification:
// each step's key a reactant of its reaction, and its conversion from 0 to 1.
pub fn new_conversion_reactor(name string, p ConversionParams) !ConversionReactor {
	check_thermal(name, p.thermal, p.t, p.duty, p.drop)!
	if p.steps.len == 0 {
		return spec_error(name, 'no reactions')
	}
	for s in p.steps {
		r := s.reaction
		k := r.species.index(s.key)
		if k < 0 || r.nu[k] >= 0 {
			what := if s.key >= 0 && s.key < p.species.len {
				p.species[s.key].name
			} else {
				'${s.key}'
			}
			return spec_error(name, '${what} is not a reactant of ${r.name}')
		}
		if e := s.expression {
			if p.ctx.variables.len == 0 || p.ctx.variables[0].name != 'T' {
				return spec_error(name, "${r.name}'s conversion is an expression with no context to read T from")
			}
			if !e.dim.is_dimensionless() {
				return spec_error(name, "${r.name}'s conversion must be a pure number")
			}
			if what := beyond_t(e, &p.defs) {
				return spec_error(name, "${r.name}'s conversion may depend on T alone, and this uses ${what}")
			}
		} else if !(s.conversion >= 0 && s.conversion <= 1) {
			return spec_error(name, "${r.name}'s conversion is ${s.conversion}, not from 0 to 1")
		}
	}
	return ConversionReactor{
		name:    name
		species: p.species
		steps:   p.steps
		thermal: p.thermal
		t:       p.t
		duty:    p.duty
		drop:    p.drop
		ctx:     p.ctx
		defs:    p.defs
	}
}

// beyond_t is what an expression uses besides T, through its definitions
// too, if anything: P, or a member such as p[CH4].
fn beyond_t(e expr.Expr, defs &expr.Definitions) ?string {
	for n in e.nodes {
		match n.kind {
			.var {
				if n.index != 0 {
					return n.text
				}
			}
			.member {
				return n.text
			}
			.def {
				if what := beyond_t(defs.exprs[n.index], defs) {
					return what
				}
			}
			else {}
		}
	}
	return none
}

// conversion_at is a step's conversion at T: its number, or its expression's
// value there, which must be from 0 to 1.
fn (c ConversionReactor) conversion_at(s Conversion, t f64) !f64 {
	e := s.expression or { return s.conversion }
	mut x := []f64{len: c.ctx.slots()}
	x[0] = t
	v := e.eval(&c.defs, x)
	if !(v >= 0 && v <= 1) {
		return spec_error(c.name, "${s.reaction.name}'s conversion is ${kinetics.sig(v)} at ${kinetics.sig(t)} K, not from 0 to 1")
	}
	return v
}

// inlet_ports is in.
pub fn (c ConversionReactor) inlet_ports() []string {
	return ['in']
}

// outlet_ports is out.
pub fn (c ConversionReactor) outlet_ports() []string {
	return ['out']
}

// flows_out runs the steps on the inlet's flows at the outlet's T: each
// step's extent is its conversion of its key as it enters the step, over the
// key's coefficient. A step that would take a reactant below zero beyond
// rounding, 1e-12 of the flow, is refused, with how far; within it, the
// reactant is zero.
fn (c ConversionReactor) flows_out(flows []f64, t f64) ![]f64 {
	total := numeric.sum(flows)
	mut n := flows.clone()
	for s in c.steps {
		r := s.reaction
		before := n[s.key]
		x := c.conversion_at(s, t)!
		extent := x * before / -r.nu[r.species.index(s.key)]
		for q, j in r.species {
			n[j] += r.nu[q] * extent
		}
		for j in r.species {
			if n[j] >= 0 {
				continue
			}
			if n[j] >= -1e-12 * total {
				n[j] = 0
				continue
			}
			return spec_error(c.name, '${r.name} needs more ${c.species[j].name} than there is: converting ${x} of its ${c.species[s.key].name} would leave ${n[j]} mol/s')
		}
	}
	return n
}

// calculate finds the outlet at T or by the energy balance, the steps run at
// the outlet's T, and warns of any reaction it takes past equilibrium.
pub fn (c ConversionReactor) calculate(fl thermo.Fluid, inlets []flowsheet.Stream) !flowsheet.Outcome {
	inlet := one_inlet(c.name, fl, inlets, c.species.len)!
	at := ConversionAt{c, fl, inlet.flows, pressure_out(c.name, inlet, c.drop)!}
	o := thermal_outcome(c.name, at, inlet, c.thermal, c.t, c.duty, element_rows(c.species))!
	return flowsheet.Outcome{
		...o
		warnings: c.past_equilibrium(o.outlets[0])
	}
}

// ConversionAt is a conversion reactor's outlet at a temperature: the steps
// run on the inlet's flows at that T, whose conversions may depend on it.
struct ConversionAt {
	c     ConversionReactor
	fl    thermo.Fluid
	flows []f64 // the inlet's
	p     f64
}

fn (a ConversionAt) outlet_at(t f64) !flowsheet.Stream {
	if !(math.is_finite(t) && t > 0) {
		return failed(a.c.name, 'no outlet at ${t} K')
	}
	return at_tp(a.c.name, a.fl, a.c.flows_out(a.flows, t)!, t, a.p)
}

// equilibrium_tol is how far past equilibrium, in ln(Q/K), rounding may put
// an outlet: Q from a flash's ln φ, and K from a quadrature, each good to
// about 1e-12.
const equilibrium_tol = 1e-8

// past_equilibrium warns of each step's reaction that the outlet is past
// equilibrium for (D10, as amended before step 10): its Q, from the outlet's
// fugacities in bar, above its K at 1 bar from the Gibbs energies (D7). A
// reaction whose species lack formation data is noted as unchecked.
fn (c ConversionReactor) past_equilibrium(outlet flowsheet.Stream) []string {
	mut out := []string{}
	st := outlet.result.phases[0].state
	total := numeric.sum(st.n)
	ln_p := math.log(st.p / 1e5)
	mut seen := []string{}
	for s in c.steps {
		r := s.reaction
		if r.name in seen {
			continue
		}
		seen << r.name
		std := r.standard(c.species, st.t) or {
			out << "${r.name}'s equilibrium wasn't checked: ${err.msg()}"
			continue
		}
		mut ln_q := 0.0
		for q, j in r.species {
			ln_q += r.nu[q] * (math.log(st.n[j] / total) + st.ln_phi[j] + ln_p)
		}
		gap := ln_q - std.ln_k
		if gap > equilibrium_tol {
			ratio := if gap < 700 { kinetics.sig(math.exp(gap)) } else { 'beyond 1e304' }
			out << '${r.name} is past equilibrium at the outlet: Q/K = ${ratio} at ${kinetics.sig(st.t)} K and ${kinetics.sig(st.p / 1e5)} bar, so no reactor could convert this much'
		}
	}
	return out
}
