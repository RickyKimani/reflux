module transport

import math
import thermo

// Result is a value, or none and why (D3, D16).
pub struct Result {
pub:
	value f64 = math.nan() // in SI; NaN when there's none
	why   string // why there's none; '' with a value
}

// ok says whether it has a value.
pub fn (r Result) ok() bool {
	return r.why == ''
}

// value_of is a mixing rule's value, or none where its arithmetic gave one
// that isn't a positive finite number: a pair's G large enough to overflow
// Grunberg and Nissan's exp, say.
fn value_of(v f64) Result {
	if !(v > 0) || math.is_inf(v, 1) {
		return none_because('the mixing rule gives ${v}, not a positive finite number')
	}
	return Result{
		value: v
	}
}

fn none_because(why string) Result {
	return Result{
		why: why
	}
}

// PhaseTransport is a phase's viscosity (Pa·s), thermal conductivity
// (W/(m·K)) and surface tension (N/m), each a value or why there's none,
// and the notes on how they were had: compounds left out of a liquid (D9),
// values taken beyond their correlations' ranges (D3).
pub struct PhaseTransport {
pub:
	viscosity       Result
	conductivity    Result
	surface_tension Result
	notes           []string
}

// Model is what gives a phase its transport: a mixture's correlations
// (Mixture) or a formulation's own, IAPWS's for water's (water's). A phase
// is its label, T (K), P (Pa), mole fractions x and molar volume v (m³/mol),
// as a flash gave it.
pub interface Model {
	phase(label thermo.PhaseLabel, t f64, p f64, x []f64, v f64) PhaseTransport
	// methods names how each property is had, for a report (D16).
	methods() []MethodUsed
}

// MethodUsed is how a model has one of a phase's properties, for a report
// (D16): the property, the phases it's for, and the method.
pub struct MethodUsed {
pub:
	property string // 'viscosity', 'thermal_conductivity' or 'surface_tension'
	phase    string // 'vapour', 'liquid' (liquid2 too) or 'any'
	method   string
}

// Pair is a liquid pair's Grunberg–Nissan Gᵢⱼ (D6), with where it came
// from, by the compounds' indices.
pub struct Pair {
pub:
	i      int
	j      int
	g      Parameter
	source string
}

// Mixture is a fluid's compounds for transport, by their correlations
// (D3 to D9), and any liquid pairs' Gᵢⱼ, every other pair's 0.
pub struct Mixture {
pub:
	pures []Pure
	pairs []Pair
}

// dissolved_limit is the share of a liquid's moles that its compounds
// above their Tc may have and leave it a value (D9): left out, the rest
// renormalised, and noted; more, and it has none.
pub const dissolved_limit = 0.10

// phase gives the phase's transport (D9):
// - a vapour, or a phase whose every compound is above its Tc, by the
//   gas's methods (D4, D5);
// - a liquid by the liquid's (D6, D7, D8), a compound above its Tc left out
//   and noted, unless such compounds are more than dissolved_limit of it.
pub fn (m Mixture) phase(label thermo.PhaseLabel, t f64, p f64, x []f64, v f64) PhaseTransport {
	n := m.pures.len
	if x.len != n {
		return all_none('${x.len} mole fractions for ${n} compounds', [])
	}
	if !(t > 0) || !(p > 0) || !(v > 0) || math.is_inf(t, 1) || math.is_inf(p, 1) || math.is_inf(v, 1) {
		return all_none('T = ${t} K, P = ${p} Pa and V = ${v} m³/mol must be positive', [])
	}
	mut total := 0.0
	for xi in x {
		if !(xi >= 0) || math.is_inf(xi, 1) {
			return all_none('the mole fractions ${x} are not a composition', [])
		}
		total += xi
	}
	if !(total > 0) {
		return all_none('the mole fractions are all 0', [])
	}
	mut present := []int{}
	mut super := 0
	for i in 0 .. n {
		if x[i] > 0 {
			present << i
			if !(t < m.pures[i].tc) {
				super++
			}
		}
	}
	y := x.map(it / total)
	if label == .vapor {
		return m.gas(t, p, y, v, [])
	}
	if super == present.len {
		return m.gas(t, p, y, v, ["every compound is above its Tc, so the gas's methods (D9)"])
	}
	return m.liquid(t, p, y)
}

// methods implements Model: D4 to D9's.
pub fn (m Mixture) methods() []MethodUsed {
	left_out := "a compound above its Tc left out, up to ${dissolved_limit} of the liquid's moles"
	gas_too := 'also a phase whose every compound is above its Tc'
	return [
		MethodUsed{'viscosity', 'vapour', "Wilke's rule on each compound's vapour viscosity, times Lucas's (1980) high-pressure factor taken against its value at zero pressure; ${gas_too}"},
		MethodUsed{'thermal_conductivity', 'vapour', "Wassiljewa's form with Mason and Saxena's Aᵢⱼ on each compound's vapour conductivity, plus Stiel and Thodos' (1964) residual on Yorizane's pseudo-critical constants; ${gas_too}"},
		MethodUsed{'viscosity', 'liquid', "Grunberg and Nissan's rule, with each pair's G, on each compound's liquid viscosity times Lucas's (1981) compressed-liquid factor; ${left_out}"},
		MethodUsed{'thermal_conductivity', 'liquid', "Li's rule on each compound's liquid conductivity times Missenard's (1970) compressed-liquid factor, its volume by its liquid density; ${left_out}"},
		MethodUsed{'surface_tension', 'liquid', "Winterfeld, Scriven and Davis's rule on each compound's surface tension, its volume by its liquid density; ${left_out}"},
	]
}

// reasons is why the phase has none of each property it has none of, each
// reason once with the properties it's for: "no viscosity or surface
// tension: …" (D16). A vapour's surface tension, which no vapour has, is
// left out.
pub fn (pt PhaseTransport) reasons(label thermo.PhaseLabel) []string {
	what := ['viscosity', 'thermal conductivity', 'surface tension']
	mut whys := []string{}
	mut whats := [][]string{}
	for k, r in [pt.viscosity, pt.conductivity, pt.surface_tension] {
		if r.ok() || (k == 2 && label == .vapor) {
			continue
		}
		at := whys.index(r.why)
		if at < 0 {
			whys << r.why
			whats << [what[k]]
		} else {
			whats[at] << what[k]
		}
	}
	mut out := []string{}
	for i, why in whys {
		w := whats[i]
		listed := if w.len == 1 { w[0] } else { '${w[..w.len - 1].join(', ')} or ${w.last()}' }
		out << 'no ${listed}: ${why}'
	}
	return out
}

fn all_none(why string, notes []string) PhaseTransport {
	return PhaseTransport{
		viscosity:       none_because(why)
		conductivity:    none_because(why)
		surface_tension: none_because(why)
		notes:           notes
	}
}

// sub is the compounds at indices idx, and y renormalised over them.
fn (m Mixture) sub(idx []int, y []f64) ([]Pure, []f64) {
	mut cs := []Pure{}
	mut ys := []f64{}
	mut total := 0.0
	for i in idx {
		cs << m.pures[i]
		ys << y[i]
		total += y[i]
	}
	return cs, ys.map(it / total)
}

// gas is a gas's transport (D4, D5), its surface tension none.
fn (m Mixture) gas(t f64, p f64, y []f64, v f64, notes_in []string) PhaseTransport {
	mut notes := notes_in.clone()
	idx := []int{len: y.len, init: index}.filter(y[it] > 0)
	cs, ys := m.sub(idx, y)
	mut mus := []Correlated{}
	mut mw := []f64{}
	mut visc_why := ''
	for c in cs {
		r := c.value(.vapour_viscosity, t) or {
			visc_why = err.msg()
			break
		}
		if !r.in_range {
			notes << c.outside(.vapour_viscosity, t)
		}
		mus << c.sources[int(Quantity.vapour_viscosity)].c
		mw << c.mw
	}
	mut viscosity := none_because(visc_why)
	mut conductivity := none_because(visc_why)
	if visc_why == '' {
		w := Wilke{mus, mw}
		mut k := w.room()
		w.at(t, mut k) or { return all_none(err.msg(), notes) }
		mu0 := w.mixed(k, ys)
		viscosity = if f := lucas_factor(cs, ys, t, p) {
			value_of(mu0 * f)
		} else {
			none_because(err.msg())
		}
		mut lambda := []f64{}
		mut why := ''
		for c in cs {
			r := c.value(.vapour_conductivity, t) or {
				why = err.msg()
				break
			}
			if !r.in_range {
				notes << c.outside(.vapour_conductivity, t)
			}
			lambda << r.value
		}
		conductivity = if why != '' {
			none_because(why)
		} else {
			lambda0 := w.conductivity(k, ys, lambda)
			if residual := residual(cs, ys, v) {
				value_of(lambda0 + residual)
			} else {
				none_because(err.msg())
			}
		}
	} else {
		conductivity = none_because("Mason and Saxena's Aᵢⱼ need each compound's viscosity, and ${visc_why}")
	}
	return PhaseTransport{
		viscosity:       viscosity
		conductivity:    conductivity
		surface_tension: none_because('a gas has no surface tension')
		notes:           unique(notes)
	}
}

fn residual(cs []Pure, ys []f64, v f64) !f64 {
	return stiel_thodos(yorizane(cs, ys)!, v)!
}

// liquid is a liquid's transport (D6 to D9).
fn (m Mixture) liquid(t f64, p f64, y []f64) PhaseTransport {
	mut notes := []string{}
	mut kept := []int{}
	mut left := []string{}
	mut out := 0.0
	for i, c in m.pures {
		if y[i] == 0 {
			continue
		}
		if t < c.tc {
			kept << i
		} else {
			left << c.name
			out += y[i]
		}
	}
	if left.len > 0 {
		above := if left.len == 1 { 'above its Tc' } else { 'above their Tc' }
		if out > dissolved_limit {
			are := if left.len == 1 { 'is' } else { 'are' }
			why := "${left.join(', ')}, ${above}, ${are} ${shown(out)} of the liquid's moles, more than ${dissolved_limit} (D9)"
			return all_none(why, notes)
		}
		notes << "${left.join(', ')}, ${above}, left out: ${shown(out)} of the liquid's moles (D9)"
	}
	cs, xs := m.sub(kept, y)
	g := m.g_matrix(kept, t) or {
		return PhaseTransport{
			viscosity:       none_because(err.msg())
			conductivity:    m.liquid_conductivity(cs, xs, t, p, mut notes)
			surface_tension: m.liquid_surface_tension(cs, xs, t, mut notes)
			notes:           unique(notes)
		}
	}
	viscosity := m.liquid_viscosity(cs, xs, g, t, p, mut notes)
	conductivity := m.liquid_conductivity(cs, xs, t, p, mut notes)
	surface := m.liquid_surface_tension(cs, xs, t, mut notes)
	return PhaseTransport{
		viscosity:       viscosity
		conductivity:    conductivity
		surface_tension: surface
		notes:           unique(notes)
	}
}

// unique is notes without repeats, in their first order: two values can
// read the same correlation (a liquid's density, for Li's rule and the
// surface tension's).
fn unique(notes []string) []string {
	mut out := []string{}
	for s in notes {
		if s !in out {
			out << s
		}
	}
	return out
}

// g_matrix is the Gᵢⱼ among the compounds at idx, n×n, at T.
fn (m Mixture) g_matrix(idx []int, t f64) ![]f64 {
	n := idx.len
	mut g := []f64{len: n * n}
	for pair in m.pairs {
		a := idx.index(pair.i)
		b := idx.index(pair.j)
		if a < 0 || b < 0 {
			continue
		}
		v := pair.g.value(t) or {
			return error('the Grunberg–Nissan G of ${m.pures[pair.i].name} and ${m.pures[pair.j].name} at ${shown(t)} K: ${err.msg()}')
		}
		if !math.is_finite(v) {
			return error('the Grunberg–Nissan G of ${m.pures[pair.i].name} and ${m.pures[pair.j].name} at ${shown(t)} K is ${v}')
		}
		g[a * n + b] = v
		g[b * n + a] = v
	}
	return g
}

// readings is each compound's q at T, or why one has none, noting each
// taken beyond its range.
fn readings(cs []Pure, q Quantity, t f64, mut notes []string) ![]f64 {
	mut out := []f64{}
	for c in cs {
		r := c.value(q, t)!
		if !r.in_range {
			notes << c.outside(q, t)
		}
		out << r.value
	}
	return out
}

// volumes is each compound's molar volume, m³/mol, from its liquid density.
fn volumes(cs []Pure, t f64, mut notes []string) ![]f64 {
	return readings(cs, .liquid_density, t, mut notes)!.map(1.0 / it)
}

// liquid_viscosity is the Grunberg–Nissan rule on each compound's viscosity
// times Lucas's factor (D6). One compound alone is its own value.
fn (m Mixture) liquid_viscosity(cs []Pure, xs []f64, g []f64, t f64, p f64, mut notes []string) Result {
	mu := readings(cs, .liquid_viscosity, t, mut notes) or { return none_because(err.msg()) }
	psat := readings(cs, .vapour_pressure, t, mut notes) or {
		return none_because("Lucas's compressed-liquid factor needs each vapour pressure: ${err.msg()}")
	}
	mut each := []f64{}
	for i, c in cs {
		each << mu[i] * lucas_liquid(t, p, c.tc, c.pc, c.omega, psat[i])
	}
	if cs.len == 1 {
		return value_of(each[0])
	}
	return value_of(grunberg_nissan(xs, each, g))
}

// liquid_conductivity is Li's rule on each compound's conductivity times
// Missenard's factor (D7). One compound alone is its own value.
fn (m Mixture) liquid_conductivity(cs []Pure, xs []f64, t f64, p f64, mut notes []string) Result {
	lambda := readings(cs, .liquid_conductivity, t, mut notes) or { return none_because(err.msg()) }
	mut each := []f64{}
	for i, c in cs {
		each << lambda[i] * missenard(t, p, c.tc, c.pc)
	}
	if cs.len == 1 {
		return value_of(each[0])
	}
	v := volumes(cs, t, mut notes) or {
		return none_because("Li's rule needs each liquid's molar volume: ${err.msg()}")
	}
	return value_of(li(xs, v, each))
}

// liquid_surface_tension is Winterfeld, Scriven and Davis' rule (D8). One
// compound alone is its own value.
fn (m Mixture) liquid_surface_tension(cs []Pure, xs []f64, t f64, mut notes []string) Result {
	sigma := readings(cs, .surface_tension, t, mut notes) or { return none_because(err.msg()) }
	if cs.len == 1 {
		return value_of(sigma[0])
	}
	v := volumes(cs, t, mut notes) or {
		return none_because("Winterfeld, Scriven and Davis' rule needs each liquid's molar volume: ${err.msg()}")
	}
	return value_of(winterfeld(xs, v, sigma))
}
