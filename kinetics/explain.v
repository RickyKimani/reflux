module kinetics

import expr
import math
import quantity

// The factor table `reflux kinetics explain` prints (D12, kinetics.md): at
// a state, each rate as reflux read it, and every part of it with its
// value, so that which term dominates shows at once.

// Row is one part of a rate, or a definition it uses, with its value at the
// state, in SI.
pub struct Row {
pub:
	label string // what it is
	text  string // as written, in Unicode
	value f64
	unit  string // its SI unit, or '' for a pure number
}

// Explanation is one reaction's rate at a state: how it's written, its
// parts and the definitions it uses with their values, and the rate.
pub struct Explanation {
pub:
	reaction int
	name     string
	equation string
	law      string // the form, its basis and what it's per
	plain    string // the rate as plain text, which reads back
	unicode  string
	latex    string
	parts    []Row
	uses     []Row // every definition it uses, in the order defined
	rate     f64
	unit     string
}

// describe says a rate law's form, basis and what it's per.
fn describe(law RateLaw) string {
	form := match law {
		PowerLaw {
			'a power law'
		}
		Lhhw {
			over := if law.terms != none {
				"its own adsorption term, in place of site ${law.site}'s,"
			} else {
				"site ${law.site}'s adsorption term"
			}
			power := match law.power {
				1.0 { '' }
				2.0 { ' squared' }
				3.0 { ' cubed' }
				else { ' to the power ${exponent_text(law.power)}' }
			}
			'LHHW over ${over}${power}'
		}
		EleyRideal {
			'Eley–Rideal on site ${law.site}'
		}
		Expression {
			'an expression'
		}
	}
	basis := match basis_of(law) {
		.pressure { 'partial pressures' }
		.fugacity { 'fugacities' }
		.concentration { 'concentrations' }
		.fraction { 'mole fractions' }
	}
	per := if per_of(law) == .mass { 'per mass of catalyst' } else { 'per volume' }
	eq := equilibrium_of(law) or { Equilibrium{} }
	rev := match eq.kind {
		.irreversible {
			if law is Expression { '' } else { ', irreversible' }
		}
		.thermo { ', reversible, K from the Gibbs energies' }
		.correlation { ', reversible, K fitted' }
	}
	return '${form} on ${basis}, ${per}${rev}'
}

// The dimensions a rate's parts are written over, each with its unit's
// numerator and denominator: a pure number, a rate per mass and per
// volume, a concentration and a pressure; a part's other factors are
// pascals or moles per cubic metre to a power.
const unit_bases = [
	UnitBase{quantity.Dim{}, '', ''},
	UnitBase{per_mass, 'mol', 'kg·s'},
	UnitBase{per_volume, 'mol', 'm³·s'},
]!

struct UnitBase {
	dim  quantity.Dim
	up   string
	down string
}

// unit_of writes a part's unit as a reader would: a rate per mass or per
// volume, or a pure number, times pascals or moles per cubic metre to a
// power, mol·Pa^0.5/(kg·s) rather than SI's base units; and the base
// units where it's none of those.
fn unit_of(e expr.Expr) string {
	pressure := quantity.Dimension.pressure.dim()
	concentration := quantity.dim(0, -3, 0, 0, 1, 0)
	for x in [pressure, concentration] {
		name := if x.eq(pressure) { 'Pa' } else { 'mol/m³' }
		for twice in -8 .. 9 {
			power := quantity.Ratio{i64(twice), 2}
			reduced := quantity.ratio(power.num, power.den) or { continue }
			rest := e.dim.div(x.pow(reduced) or { continue }) or { continue }
			for b in unit_bases {
				if !rest.eq(b.dim) {
					continue
				}
				return compose(b, name, reduced)
			}
		}
	}
	s := e.dim.str()
	return if s == '1' { '' } else { s }
}

// compose writes a base's unit times a unit to a power.
fn compose(b UnitBase, name string, power quantity.Ratio) string {
	mut up := if b.up == '' { []string{} } else { [b.up] }
	mut down := if b.down == '' { []string{} } else { [b.down] }
	if power.num != 0 {
		a := if power.num < 0 { power.neg() } else { power }
		w := if a.num == 1 && a.den == 1 {
			name
		} else if name.contains('/') {
			'(${name})^${a.value()}'
		} else {
			'${name}^${exponent_text(a.value())}'
		}
		if power.num > 0 {
			up << w
		} else {
			down << w
		}
	}
	n := if up.len == 0 { '1' } else { up.join('·') }
	if down.len == 0 {
		return if n == '1' { '' } else { n }
	}
	d := down.join('·')
	return if down.len > 1 || d.contains('·') { '${n}/(${d})' } else { '${n}/${d}' }
}

// explain evaluates every rate and its parts at a state.
pub fn (k &Kinetics) explain(s State) ![]Explanation {
	mut pt := k.point()
	k.fill(s, mut pt)!
	mut out := []Explanation{}
	for j, rate in k.rates {
		i := k.rated[j]
		r := k.reactions[i]
		mut parts := []Row{}
		for part in k.parts[j] {
			e := expr.parse(k.ctx, &k.defs, part.text) or { continue }
			parts << Row{part.label, e.unicode(), e.eval(&k.defs, pt.x), unit_of(e)}
		}
		mut uses := []Row{}
		for d in k.defs.order {
			if !uses_definition(k.defs, rate, d) {
				continue
			}
			e := k.defs.exprs[d]
			uses << Row{k.defs.names[d], e.unicode(), e.eval(&k.defs, pt.x), unit_of(e)}
		}
		out << Explanation{
			reaction: i
			name:     reaction_name(k.reactions, i)
			equation: r.str()
			law:      describe(k.laws[j])
			plain:    rate.plain()
			unicode:  rate.unicode()
			latex:    rate.latex()
			parts:    parts
			uses:     uses
			rate:     rate.eval(&k.defs, pt.x)
			unit:     unit_of(rate)
		}
	}
	return out
}

// uses_definition is whether an expression uses definition d, directly or
// through others.
fn uses_definition(defs expr.Definitions, e expr.Expr, d int) bool {
	mut seen := []int{}
	mut todo := e.uses()
	for todo.len > 0 {
		x := todo.pop()
		if x == d {
			return true
		}
		if x in seen {
			continue
		}
		seen << x
		todo << defs.exprs[x].uses()
	}
	return false
}

// value writes a value with its unit.
fn value(v f64, unit string) string {
	return if unit == '' { sig(v) } else { '${sig(v)} ${unit}' }
}

// str writes an explanation as `reflux kinetics explain` prints it.
pub fn (e Explanation) str() string {
	mut lines := ['${e.name}: ${e.equation}', '  ${e.law}', '  r = ${e.unicode}']
	if e.uses.len > 0 {
		lines << '  where'
		width := arrays_max(e.uses.map(it.label.runes().len))
		for u in e.uses {
			lines << '    ${u.label}${' '.repeat(width - u.label.runes().len)} = ${u.text} = ${value(u.value, u.unit)}'
		}
	}
	if e.parts.len > 0 {
		lines << '  its parts'
		width := arrays_max(e.parts.map(it.label.runes().len))
		for p in e.parts {
			lines << '    ${p.label}${' '.repeat(width - p.label.runes().len)}  ${p.text} = ${value(p.value, p.unit)}'
		}
	}
	lines << '  r = ${value(e.rate, e.unit)}'
	return lines.join('\n')
}

fn arrays_max(xs []int) int {
	mut m := 0
	for x in xs {
		m = math.max(m, x)
	}
	return m
}
