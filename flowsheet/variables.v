module flowsheet

import expr
import kinetics
import math
import quantity

// Energies is a unit that names its energy streams before it calculates
// (M8's D3), so that an expression can name them: reactor.duty, pump.work.
// Each name is one its Outcome.energy gives. A unit without it has none to
// name, as a mixer has none.
pub interface Energies {
	energy_ports() []string
}

// Measure is what a flowsheet variable measures.
enum Measure {
	t      // a stream's T, K
	p      // its P, Pa
	f      // its total molar flow, mol/s
	h      // its molar enthalpy, J/mol
	vf     // its vapour fraction
	n      // a component's molar flow, mol/s
	y      // a component's mole fraction, of the whole stream
	energy // a unit's energy stream, W
}

// Slot is where one of an expression's inputs comes from.
struct Slot {
	stream    string
	unit      int = -1
	measure   Measure
	component int
	energy    string
}

// Variables is what an expression over a flowsheet may name (M8's D3), and
// where each name's value comes from in a solution:
// - each stream's T, P, F (its molar flow), h (its molar enthalpy) and vf
//   (its vapour fraction), and each component's flow n[…] and mole
//   fraction y[…], product.n[CH4];
// - each unit's energy streams by their names, reactor.duty, pump.work.
pub struct Variables {
pub:
	context expr.Context
mut:
	slots []Slot
}

// variables are a flowsheet's names for expressions. members name the
// components as an expression writes them back, and keys give every name a
// component is known by.
pub fn (f Flowsheet) variables(members []string, keys map[string]int) Variables {
	t := quantity.Dimension.temperature.dim()
	pressure := quantity.Dimension.pressure.dim()
	flow := quantity.Dimension.molar_flow.dim()
	enthalpy := quantity.Dimension.molar_enthalpy.dim()
	power := quantity.Dimension.power.dim()
	mut variables := []expr.Variable{}
	mut families := []expr.Family{}
	mut scalar := []Slot{}
	mut scopes := []expr.Scope{}
	for name in f.stream_names() {
		scopes << expr.Scope{name, 'stream'}
		for v in [
			expr.Variable{'${name}.T', t},
			expr.Variable{'${name}.P', pressure},
			expr.Variable{'${name}.F', flow},
			expr.Variable{'${name}.h', enthalpy},
			expr.Variable{'${name}.vf', quantity.Dim{}},
		] {
			variables << v
		}
		for m in [Measure.t, .p, .f, .h, .vf] {
			scalar << Slot{
				stream:  name
				measure: m
			}
		}
		families << expr.Family{'${name}.n', flow, members, keys}
		families << expr.Family{'${name}.y', quantity.Dim{}, members, keys}
	}
	for i, p in f.units {
		u := p.unit
		name := u.name
		k := scopes.map(it.name).index(name)
		if k >= 0 {
			scopes[k] = expr.Scope{name, 'stream and unit'}
		} else {
			scopes << expr.Scope{name, 'unit'}
		}
		if u is Energies {
			for e in u.energy_ports() {
				variables << expr.Variable{'${name}.${e}', power}
				scalar << Slot{
					unit:    i
					measure: .energy
					energy:  e
				}
			}
		}
	}
	mut slots := scalar.clone()
	for fam in families {
		measure := if fam.name.ends_with('.n') { Measure.n } else { Measure.y }
		for c in 0 .. members.len {
			slots << Slot{
				stream:    fam.name.all_before_last('.')
				measure:   measure
				component: c
			}
		}
	}
	return Variables{
		context: expr.Context{
			variables: variables
			families:  families
			scopes:    scopes
		}
		slots:   slots
	}
}

// stream_names is every stream once: the feeds, then each unit's outlets,
// unit by unit, as a report lists them.
fn (f Flowsheet) stream_names() []string {
	mut out := []string{}
	for feed in f.feeds {
		out << feed.name
	}
	for p in f.units {
		for name in p.outlets {
			if name !in out {
				out << name
			}
		}
	}
	return out
}

// member_names are the components as expressions name them (M8's D3, as
// kinetics.rate_context has them): each by its formula where no other
// component shares it, else by the file's name; and each answering to the
// file's name, its compound's name, CAS number, aliases and that formula.
pub fn member_names(names []string, species []kinetics.Species) ([]string, map[string]int) {
	mut hill := map[string]int{}
	for s in species {
		if f := s.formula {
			hill[f.str()]++
		}
	}
	mut members := []string{}
	mut keys := map[string]int{}
	for i, name in names {
		mut member := name
		mut all := [name]
		if i < species.len {
			s := species[i]
			all << s.name
			all << s.cas
			all << s.aliases
			if f := s.formula {
				if hill[f.str()] == 1 {
					member = f.str()
					all << member
				}
			}
		}
		members << member
		for k in all {
			if k != '' && k !in keys {
				keys[k] = i
			}
		}
	}
	return members, keys
}

// values are each slot's value in a solution, NaN where there's none, with
// why.
pub fn (v Variables) values(sol Solution) ([]f64, []string) {
	return v.fill(sol.streams, sol.units)
}

// fill is each slot's value from these streams and units' last results.
fn (v Variables) fill(streams map[string]Stream, runs []UnitResult) ([]f64, []string) {
	mut x := []f64{len: v.slots.len}
	mut why := []string{len: v.slots.len}
	for k, s in v.slots {
		x[k], why[k] = slot_value(s, streams, runs)
	}
	return x, why
}

fn slot_value(s Slot, streams map[string]Stream, runs []UnitResult) (f64, string) {
	if s.measure == .energy {
		r := runs[s.unit] or { return math.nan(), 'it was not calculated' }
		for e in r.energy {
			if e.name == s.energy {
				return e.rate, ''
			}
		}
		return math.nan(), "it wasn't calculated, or reported no ${s.energy}"
	}
	st := streams[s.stream] or { return math.nan(), '${s.stream} was not calculated' }
	return match s.measure {
		.t {
			st.t, ''
		}
		.p {
			st.p, ''
		}
		.f {
			st.total(), ''
		}
		.h {
			st.h, ''
		}
		.vf {
			st.beta, ''
		}
		.n {
			st.flows[s.component], ''
		}
		.y {
			total := st.total()
			if total > 0.0 {
				st.flows[s.component] / total, ''
			} else {
				math.nan(), '${s.stream} has no flow'
			}
		}
		.energy {
			math.nan(), ''
		}
	}
}

// expression reads an expression over the file's flowsheet (M8's D3): a
// spec's side, a study's output.
pub fn (l Loaded) expression(text string) !expr.Expr {
	return expr.parse(l.variables.context, &no_definitions, text)
}

// value is an expression's value in a solution, or why it has none (M8's
// D9): never ±∞ or NaN.
pub fn (l Loaded) value(e expr.Expr, sol Solution) !f64 {
	x, why := l.variables.values(sol)
	return e.value(&no_definitions, x, why)
}

const no_definitions = expr.Definitions{}
