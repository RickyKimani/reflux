module kinetics

import expr
import math

// The checks `reflux kinetics check` runs (D12, kinetics.md), on what has
// already been read and built. What build refuses (an equation that doesn't
// balance, a rate in the wrong units, a name that isn't there) never gets
// this far; these are what's legal but worth a look.

// Severity is how much a finding asks of its reader: a note to know, or a
// warning to look into.
pub enum Severity {
	note
	warning
}

// Finding is one check's result: how much it matters, what it's about (a
// reaction by its index, a site, or neither, for the set), and what it
// says.
pub struct Finding {
pub:
	severity Severity
	reaction int = -1
	site     string
	message  string
}

// str writes a finding as `reflux kinetics check` prints it.
pub fn (f Finding) str() string {
	return '${f.severity}: ${f.message}'
}

// k_limit is how far a fitted K may be from the Gibbs energies' before the
// check warns rather than notes: a factor of 2.
const k_limit = math.ln2

// The temperatures a fitted K is held to the Gibbs energies' at, K.
const k_temperatures = [300.0, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500]

// gap_limit is how far a reaction's Σν ΔG°f may be from what its ΔH°f and
// S° give at 298.15 K before it's warned of, J/mol: 1 kJ/mol moves K by a
// factor of 1.5 there.
const gap_limit = 1000.0

// check runs every check on the kinetics, and, given a state, says how far
// each typed constant moves the rates there.
pub fn (k &Kinetics) check(s ?State) []Finding {
	mut out := []Finding{}
	out << k.check_set()
	for i, r in k.reactions {
		out << k.check_formation(i, r)
	}
	for j, law in k.laws {
		out << k.check_law(j, law, s)
	}
	for site in k.sites {
		for t in site.terms {
			out << check_adsorption(t, -1, site.name)
		}
	}
	out << k.check_look_alikes(s)
	return out
}

fn name_of(k &Kinetics, i int) string {
	return reaction_name(k.reactions, i)
}

// check_set notes each reaction that combines others.
fn (k &Kinetics) check_set() []Finding {
	d := dependence(k.reactions, k.species.len) or { return [] }
	mut out := []Finding{}
	for c in d.combinations {
		out << Finding{
			severity: .note
			reaction: c.reaction
			message:  '${c.describe(k.reactions)}: fine for kinetics, each having its own rate, but an equilibrium over these reactions needs them independent (${d.rank} of ${k.reactions.len} are)'
		}
	}
	return out
}

// check_formation warns of a reaction whose data's ΔG°f disagree with
// their ΔH°f and S°.
fn (k &Kinetics) check_formation(i int, r Reaction) []Finding {
	gap := r.formation_gap(k.species) or { return [] }
	if math.abs(gap) <= gap_limit {
		return []
	}
	return [
		Finding{
			severity: .warning
			reaction: i
			message:  "${name_of(k, i)}: its species' ΔG°f sum to ${sig(gap / 1e3)} kJ/mol from what their ΔH°f and S° give at 298.15 K; K comes from ΔH°f and S° (D7), so one of the three is wrong"
		},
	]
}

// molar_energy is an energy per mole written with its unit, in J/mol.
fn molar_energy(text string) ?f64 {
	e := expr.parse(expr.Context{}, &expr.Definitions{}, text) or { return none }
	return e.eval(&expr.Definitions{}, []f64{})
}

// check_law checks a rate law's constants and, for a fitted K, how far it
// is from the Gibbs energies'.
fn (k &Kinetics) check_law(j int, law RateLaw, s ?State) []Finding {
	i := k.rated[j]
	mut out := []Finding{}
	kf := match law {
		PowerLaw, Lhhw, EleyRideal { ?TempFn(law.k) }
		Expression { ?TempFn(none) }
	}
	if f := kf {
		out << check_constant(f, "${name_of(k, i)}'s rate constant", i)
		if f is Arrhenius {
			if ea := molar_energy(f.ea) {
				if ea < 0 {
					out << Finding{
						severity: .warning
						reaction: i
						message:  "${name_of(k, i)}'s rate constant has a negative activation energy, ${sig(ea / 1e3)} kJ/mol: an apparent one can be, from constants lumped together, but an elementary step's can't"
					}
				}
			}
		}
	}
	if law is Lhhw {
		if terms := law.terms {
			out << k.check_override(i, law, terms)
			for t in terms {
				out << check_adsorption(t, i, '')
			}
		}
	}
	if eq := equilibrium_of(law) {
		if eq.kind == .correlation {
			out << k.check_fitted(i, eq, basis_of(law), s)
		}
	}
	return out
}

// check_constant warns of a constant's A at or below zero.
fn check_constant(f TempFn, what string, i int) []Finding {
	a := match f {
		Arrhenius { f.a }
		VantHoff { f.a }
		Fixed { f.value }
		Written { '' }
	}
	if a == '' {
		return []
	}
	e := expr.parse(expr.Context{}, &expr.Definitions{}, a) or { return [] }
	if e.eval(&expr.Definitions{}, []f64{}) <= 0 {
		return [
			Finding{
				severity: .warning
				reaction: i
				message:  '${what} is ${a}: a rate or an adsorption constant is above zero'
			},
		]
	}
	return []
}

// check_adsorption warns of an adsorption term's constant that rises with
// T: adsorption is exothermic, so its ΔH is usually below zero.
fn check_adsorption(t Monomial, i int, site string) []Finding {
	what := if t.name != '' { t.name } else { 'a term' }
	mut out := check_constant(t.k, what, i)
	if t.k is VantHoff {
		if dh := molar_energy(t.k.dh) {
			if dh > 0 {
				out << Finding{
					severity: .warning
					reaction: i
					site:     site
					message:  '${what} rises with T: its ΔH, ${sig(dh / 1e3)} kJ/mol, is above zero, where adsorption is exothermic (unless the term lumps more than an adsorption)'
				}
			}
		}
	}
	return out
}

// check_override notes a reaction's own adsorption term: the site's belongs
// to every reaction on it (kinetics.md), so one of a reaction's own is the
// exception, and this says how it differs.
fn (k &Kinetics) check_override(i int, law Lhhw, terms []Monomial) []Finding {
	site := k.sites.filter(it.name == law.site)
	if site.len == 0 {
		return []
	}
	key := fn (t Monomial) string {
		return t.powers.map('${it.species}^${it.exponent}').join(' ')
	}
	own := terms.map(key(it))
	theirs := site[0].terms.map(key(it))
	mut missing := []string{}
	for x, t in site[0].terms {
		if theirs[x] !in own {
			missing << if t.name != '' { t.name } else { theirs[x] }
		}
	}
	mut extra := []string{}
	for x, t in terms {
		if own[x] !in theirs {
			extra << if t.name != '' { t.name } else { own[x] }
		}
	}
	mut differs := []string{}
	if missing.len > 0 {
		differs << 'without ${missing.join(', ')}'
	}
	if extra.len > 0 {
		differs << 'with ${extra.join(', ')} too'
	}
	how := if differs.len > 0 { ', ${differs.join(' and ')}' } else { ', the same terms' }
	return [
		Finding{
			severity: .note
			reaction: i
			site:     law.site
			message:  "${name_of(k, i)} has its own adsorption term in place of site ${law.site}'s${how}: a site's term belongs to every reaction on it"
		},
	]
}

// check_fitted says how far a fitted K is from the Gibbs energies', in the
// rate's activities: at the state's T, where there's a state, which it's
// judged by, and from 300 to 1500 K, the range it can be compared over. A
// correlation is fitted over a narrower range, so a distance far from where
// the reaction runs says little.
fn (k &Kinetics) check_fitted(i int, eq Equilibrium, basis Basis, s ?State) []Finding {
	r := k.reactions[i]
	name := if eq.name != '' { eq.name } else { 'K${i + 1}' }
	d := k.defs.find(name)
	if d < 0 {
		return []
	}
	mut x := []f64{len: k.ctx.slots()}
	distance := fn [k, r, d, basis, mut x] (t f64) !f64 {
		st := r.standard(k.species, t)!
		x[0] = t
		x[1] = 1e5
		fitted := k.defs.exprs[d].eval(&k.defs, x)
		dn := r.dn().value()
		standard := match basis {
			.pressure, .fugacity { dn * math.log(1e5) }
			.concentration { dn * math.log(1e5 / (gas_constant * t)) }
			.fraction { 0.0 }
		}
		return math.log(fitted) - (st.ln_k + standard)
	}
	mut lo, mut hi := math.inf(1), math.inf(-1)
	for t in k_temperatures {
		v := distance(t) or {
			return [
				Finding{
					severity: .note
					reaction: i
					message:  "${name_of(k, i)}'s fitted K can't be held to the Gibbs energies': ${err.msg()}"
				},
			]
		}
		lo, hi = math.min(lo, v), math.max(hi, v)
	}
	span := 'from 300 to 1500 K its ln K is ${signed(sig(lo))} to ${signed(sig(hi))} from theirs'
	if st := s {
		v := distance(st.t) or { return [] }
		sev := if math.abs(v) > k_limit { Severity.warning } else { Severity.note }
		return [
			Finding{
				severity: sev
				reaction: i
				message:  "${name_of(k, i)}'s fitted K, ${name}, is a factor of ${sig(math.exp(v))} of the Gibbs energies' at ${sig(st.t)} K (ln K ${signed(sig(v))}); ${span}"
			},
		]
	}
	worst := if math.abs(lo) > math.abs(hi) { lo } else { hi }
	sev := if math.abs(worst) > k_limit { Severity.warning } else { Severity.note }
	return [
		Finding{
			severity: sev
			reaction: i
			message:  "${name_of(k, i)}'s fitted K, ${name}: ${span}"
		},
	]
}

// check_look_alikes flags each number typed where a built-in constant was
// likely meant, 8.314 J/(mol·K) for R, and, given a state, how far it moves
// each rate there.
fn (k &Kinetics) check_look_alikes(s ?State) []Finding {
	mut out := []Finding{}
	mut pt := k.point()
	mut have := false
	if st := s {
		k.fill(st, mut pt) or {}
		have = pt.x[0] > 0
	}
	for d, e in k.defs.exprs {
		for la in e.look_alikes() {
			out << k.look_alike(la, k.defs.names[d], d, -1, have, pt)
		}
	}
	for j, e in k.rates {
		for la in e.look_alikes() {
			out << k.look_alike(la, "${name_of(k, k.rated[j])}'s rate", -1, j, have, pt)
		}
	}
	return out
}

fn (k &Kinetics) look_alike(la expr.LookAlike, where string, d int, j int, have bool, pt Point) Finding {
	c := expr.constants.filter(it.name == la.constant)[0]
	direction := if la.relative < 0 { 'below' } else { 'above' }
	mut msg := '${where} has ${sig(la.value)} ${c.unit}, ${sig(math.abs(la.relative))} ${direction} ${c.name}, ${c.about}: write ${c.name} unless the paper fitted with this number'
	o := k.owners[if d >= 0 { k.defs.names[d] } else { '' }] or { Owner{} }
	reaction := if j >= 0 { k.rated[j] } else { o.reaction }
	if have {
		moved := k.moved(la, d, j, pt)
		if moved.len > 0 {
			msg += '; at ${sig(pt.x[0])} K it moves ${moved.join(', ')} from what ${c.name} gives'
		}
	}
	return Finding{
		severity: .warning
		reaction: reaction
		site:     o.site
		message:  msg
	}
}

// moved says how far each rate moves at a point when the look-alike is the
// constant itself, by evaluating the trees with the number swapped.
fn (k &Kinetics) moved(la expr.LookAlike, d int, j int, pt Point) []string {
	c := expr.constants.filter(it.name == la.constant)[0]
	mut defs := k.defs
	mut rates := k.rates.clone()
	if d >= 0 {
		mut exprs := k.defs.exprs.clone()
		exprs[d] = swapped(exprs[d], la.node, c.value)
		defs = expr.Definitions{
			names: k.defs.names
			exprs: exprs
			order: k.defs.order
		}
	} else {
		rates[j] = swapped(rates[j], la.node, c.value)
	}
	mut out := []string{}
	for q, e in k.rates {
		before := e.eval(&k.defs, pt.x)
		after := rates[q].eval(&defs, pt.x)
		if before != after && before != 0 {
			out << "${name_of(k, k.rated[q])}'s rate by ${signed(sig((before / after - 1) * 100))}%"
		}
	}
	return out
}

fn swapped(e expr.Expr, node int, value f64) expr.Expr {
	mut nodes := e.nodes.clone()
	nodes[node] = expr.Node{
		...e.nodes[node]
		value: value
	}
	return expr.Expr{
		...e
		nodes: nodes
	}
}

fn signed(s string) string {
	return if s.starts_with('-') { s } else { '+${s}' }
}

// sig writes x to five significant figures: its digits rounded by
// math.round to a whole number and written out, so that the text is the
// same whichever C compiler built reflux, in plain notation from 1e-4 to
// 1e6 and with an exponent beyond, 1.3155e10.
pub fn sig(x f64) string {
	if x == 0 || !math.is_finite(x) {
		return '${x}'
	}
	a := math.abs(x)
	mut e := int(math.floor(math.log10(a)))
	mut m := i64(math.round(a / math.pow(10, e - 4)))
	if m < 10000 {
		e--
		m = i64(math.round(a / math.pow(10, e - 4)))
	}
	if m >= 100000 {
		e++
		m = i64(math.round(a / math.pow(10, e - 4)))
	}
	sign := if x < 0 { '-' } else { '' }
	digits := '${m}'.trim_right('0')
	if e >= 6 || e < -4 {
		frac := digits[1..]
		return '${sign}${digits[0..1]}${if frac == '' { '' } else { '.' + frac }}e${e}'
	}
	if e >= digits.len - 1 {
		return '${sign}${digits}${'0'.repeat(e - digits.len + 1)}'
	}
	if e >= 0 {
		return '${sign}${digits[..e + 1]}.${digits[e + 1..]}'
	}
	return '${sign}0.${'0'.repeat(-e - 1)}${digits}'
}
