module kinetics

import expr
import quantity

// Kinetics is a set of reactions with their rate laws, read, bound and
// checked, and compiled onto one tape: each rate law's expression, the
// definitions they use (each rate constant, each site's adsorption term
// and its constants, each fitted K, and the user's own), and the context
// they're bound in.
pub struct Kinetics {
pub:
	species   []Species
	reactions []Reaction
	sites     []Site
	laws      []RateLaw // one per rate, in the order of rated
	rated     []int     // the reaction each rate is
	ctx       expr.Context
	defs      expr.Definitions
	rates     []expr.Expr
	tape      expr.Tape
	k_used    []int // the reactions whose Keq[...] some rate uses
}

// Rated is a reaction's rate law, for build.
pub struct Rated {
pub:
	reaction int     @[required]
	law      RateLaw @[required]
}

// The rate's dimension on each basis: mol/(kg·s), and mol/(m³·s).
const per_mass = quantity.dim(-1, 0, -1, 0, 1, 0)
const per_volume = quantity.dim(0, -3, -1, 0, 1, 0)

// rate_context is a rate's context: T and P; each species' partial
// pressure p, fugacity f, concentration c and mole fraction y; and each
// reaction's K at 1 bar from the Gibbs energies, Keq[...], a pure number
// (not K, which would make every number in kelvins ambiguous). A species
// is a member by the name it's written with in the equations,
// and answers to that, its name, its aliases, its CAS number and its
// formula in Hill's order where no other species shares it.
pub fn rate_context(species []Species, reactions []Reaction) expr.Context {
	mut members := []string{len: species.len}
	for r in reactions {
		for j, k in r.species {
			if members[k] == '' {
				members[k] = r.written[j]
			}
		}
	}
	mut hill_count := map[string]int{}
	for s in species {
		if f := s.formula {
			hill_count[f.str()]++
		}
	}
	mut keys := map[string]int{}
	for i, s in species {
		if members[i] == '' {
			members[i] = s.name
			if f := s.formula {
				if hill_count[f.str()] == 1 {
					members[i] = f.str()
				}
			}
		}
		mut names := [members[i], s.name, s.cas]
		names << s.aliases
		if f := s.formula {
			if hill_count[f.str()] == 1 {
				names << f.str()
			}
		}
		for n in names {
			if n != '' && n !in keys {
				keys[n] = i
			}
		}
	}
	for r in reactions {
		for j, k in r.species {
			if r.written[j] !in keys {
				keys[r.written[j]] = k
			}
		}
	}
	mut rkeys := map[string]int{}
	mut rnames := []string{}
	for i, r in reactions {
		name := if r.name != '' { r.name } else { 'r${i + 1}' }
		rnames << name
		rkeys[name] = i
		rkeys['r${i + 1}'] = i
	}
	pressure := quantity.Dimension.pressure.dim()
	return expr.Context{
		variables: [
			expr.Variable{'T', quantity.Dimension.temperature.dim()},
			expr.Variable{'P', pressure},
		]
		families:  [
			expr.Family{'p', pressure, members, keys},
			expr.Family{'f', pressure, members, keys},
			expr.Family{'c', quantity.dim(0, -3, 0, 0, 1, 0), members, keys},
			expr.Family{'y', quantity.Dim{}, members, keys},
			expr.Family{'Keq', quantity.Dim{}, rnames, rkeys},
		]
	}
}

// KineticsError is a failure to build kinetics: which reaction or site, and
// why, with an expression's own error inside where one failed.
pub struct KineticsError {
	Error
pub:
	what   string
	reason string
}

pub fn (e KineticsError) msg() string {
	return '${e.what}: ${e.reason}'
}

fn refuse(what string, reason string) IError {
	return KineticsError{
		what:   what
		reason: reason
	}
}

// Builder holds what build writes: the definitions' sources, in order.
struct Builder {
	species   []Species
	reactions []Reaction
	ctx       expr.Context
mut:
	sources []expr.Source
	names   map[string]string // a definition's name, to what it's for
	texts   map[string]string // and its text
}

fn (mut b Builder) define(name string, text string, what string) ! {
	if name in b.names {
		return refuse(what, '${name} already names ${b.names[name]}: give this one a name of its own')
	}
	b.names[name] = what
	b.texts[name] = text
	b.sources << expr.Source{name, text}
}

// share defines a site term's constant, or uses the one already defined by
// that name if its text is the same: an override may share its site's
// constants.
fn (mut b Builder) share(name string, text string, what string) ! {
	if b.texts[name] or { '' } == text && name in b.names {
		return
	}
	b.define(name, text, what)!
}

// member writes a species, as written anywhere find takes it, as the member
// of a family: p[CH4], or p['n-butane'].
fn (b &Builder) member(family string, species string, what string) !string {
	k := find(b.species, species) or { return refuse(what, err.msg()) }
	m := b.ctx.families[0].members[k]
	quoted := if m.bytes().all(it.is_letter() || it.is_digit() || it == `_`) {
		m
	} else {
		"'${m}'"
	}
	return '${family}[${quoted}]'
}

// build reads a set of reactions' rate laws and the sites they use, with
// the user's own definitions, into Kinetics: every name bound, every
// rate's units checked against its basis, and all compiled onto one tape.
pub fn build(species []Species, reactions []Reaction, sites []Site, laws []Rated, own []expr.Source) !Kinetics {
	ctx := rate_context(species, reactions)
	mut b := Builder{
		species:   species
		reactions: reactions
		ctx:       ctx
	}
	for s in own {
		b.define(s.name, s.text, 'a definition of your own')!
	}
	mut site_symbols := map[string]string{}
	for site in sites {
		symbol := if site.symbol != '' { site.symbol } else { 'theta_${site.name}' }
		if site.name in site_symbols {
			return refuse('site ${site.name}', 'two sites are called ${site.name}')
		}
		site_symbols[site.name] = symbol
		theta := b.terms(site.terms, 'site ${site.name}')!
		b.define(symbol, theta, "site ${site.name}'s adsorption term")!
	}
	mut texts := []string{}
	mut rated := []int{}
	mut kept := []RateLaw{}
	for l in laws {
		if l.reaction < 0 || l.reaction >= reactions.len {
			return refuse('a rate law', 'no reaction ${l.reaction + 1}')
		}
		if l.reaction in rated {
			return refuse(reaction_name(reactions, l.reaction), 'two rate laws')
		}
		texts << b.rate(l.reaction, l.law, site_symbols, sites)!
		rated << l.reaction
		kept << l.law
	}
	defs := expr.define(ctx, b.sources) or {
		e := err as expr.ExprError
		return refuse('definition ${e.definition}', '${b.names[e.definition]}: ${e.msg()}')
	}
	// The fitted K of each reaction must be in its rate's activities,
	// checked before the rates, which would otherwise refuse it less clearly.
	for i, l in kept {
		eq := equilibrium_of(l) or { continue }
		if eq.kind != .correlation {
			continue
		}
		name := b.k_eq_name(rated[i], eq)
		k := defs.find(name)
		r := reactions[rated[i]]
		want := activity_dim(basis_of(l)).pow(r.dn()) or {
			return refuse(reaction_name(reactions, rated[i]), err.msg())
		}
		defs.exprs[k].must_be(want, "its K, in its rate's activities,") or {
			return refuse("${reaction_name(reactions, rated[i])}'s K", err.msg())
		}
	}
	mut rates := []expr.Expr{}
	for i, text in texts {
		what := reaction_name(reactions, rated[i])
		e := expr.parse(ctx, &defs, text) or { return refuse("${what}'s rate", err.msg()) }
		per := per_of(kept[i])
		d := if per == .mass { per_mass } else { per_volume }
		e.must_be(d, if per == .mass { 'a rate per mass of catalyst' } else { 'a rate per volume' }) or {
			return refuse("${what}'s rate", err.msg())
		}
		rates << e
	}
	tape := expr.compile(ctx, &defs, rates)
	mut k_used := []int{}
	k_slot := ctx.slot(4, 0)
	for i in 0 .. reactions.len {
		if uses_input(tape, k_slot + i) {
			k_used << i
		}
	}
	return Kinetics{
		species:   species
		reactions: reactions
		sites:     sites
		laws:      kept
		rated:     rated
		ctx:       ctx
		defs:      defs
		rates:     rates
		tape:      tape
		k_used:    k_used
	}
}

fn uses_input(t expr.Tape, slot int) bool {
	for n in t.instrs {
		if n.op == .input && n.a == slot {
			return true
		}
	}
	return false
}

fn reaction_name(reactions []Reaction, i int) string {
	return if reactions[i].name != '' { reactions[i].name } else { 'reaction ${i + 1}' }
}

// dn is Δn, the sum of the coefficients, exactly.
pub fn (r Reaction) dn() Ratio {
	mut s := zero
	for c in r.exact {
		s = s.add(c) or { return s }
	}
	return s
}

fn activity_dim(b Basis) quantity.Dim {
	return match b {
		.pressure, .fugacity { quantity.Dimension.pressure.dim() }
		.concentration { quantity.dim(0, -3, 0, 0, 1, 0) }
		.fraction { quantity.Dim{} }
	}
}

fn per_of(l RateLaw) Per {
	return match l {
		PowerLaw, Lhhw, EleyRideal, Expression { l.per }
	}
}

fn basis_of(l RateLaw) Basis {
	return match l {
		PowerLaw, Lhhw, EleyRideal, Expression { l.basis }
	}
}

fn equilibrium_of(l RateLaw) ?Equilibrium {
	return match l {
		PowerLaw, Lhhw, EleyRideal { l.equilibrium }
		Expression { none }
	}
}

// terms writes a site's adsorption term, 1 + Σ K·Π a^e, defining each K.
fn (mut b Builder) terms(terms []Monomial, what string) !string {
	mut parts := ['1']
	for t in terms {
		if t.powers.len == 0 {
			return refuse(what, 'a term with no species')
		}
		mut factors := []Factor{}
		for p in t.powers {
			factors << Factor{b.member('p', p.species, what)!, p.exponent}
		}
		name := if t.name != '' {
			t.name
		} else {
			'K_' + b.ctx.families[0].members[find(b.species, t.powers[0].species)!]
		}
		if !identifier(name) {
			return refuse(what, "${name} can't name a constant: give the term a name, letters, digits and _")
		}
		b.share(name, t.k.text(), 'a term of ${what}')!
		up, down := product(factors)
		parts << ratio(if up == '' { name } else { '${name}*${up}' }, down)
	}
	return parts.join(' + ')
}

fn identifier(s string) bool {
	return s.len > 0 && s[0].is_letter() && s.bytes().all(it.is_letter() || it.is_digit() || it == `_`)
}

// k_eq_name is the name of a reaction's fitted K.
fn (b &Builder) k_eq_name(i int, eq Equilibrium) string {
	return if eq.name != '' { eq.name } else { 'K${i + 1}' }
}

// k_eq writes a reaction's K in its rate's activities: K at 1 bar times
// the standard state's activity to the Δn, or its fitted correlation.
fn (mut b Builder) k_eq(i int, eq Equilibrium, basis Basis) !string {
	r := b.reactions[i]
	what := reaction_name(b.reactions, i)
	if eq.kind == .correlation {
		name := b.k_eq_name(i, eq)
		b.define(name, eq.correlation.text(), "${what}'s K")!
		return name
	}
	m := b.ctx.families[4].members[i]
	key := if identifier(m) { m } else { "'${m}'" }
	dn := r.dn()
	if dn.num == 0 {
		return 'Keq[${key}]'
	}
	power := if dn.den == 1 && dn.num > 0 { decimal(dn) } else { '(${decimal(dn)})' }
	standard := match basis {
		.pressure, .fugacity { 'p_std' }
		.concentration { '(p_std/(R*T))' }
		.fraction { '(p_std/P)' }
	}
	return '(Keq[${key}]*${standard}^${power})'
}

// orders writes a rate's orders as factors on its basis.
fn (b &Builder) orders(family string, orders []Power, what string) ![]Factor {
	mut out := []Factor{}
	for o in orders {
		out << Factor{b.member(family, o.species, what)!, o.exponent}
	}
	return out
}

// quotient is the reaction quotient's factors, Π aᵢ^νᵢ.
fn (b &Builder) quotient(family string, i int) ![]Factor {
	r := b.reactions[i]
	mut out := []Factor{}
	for j, k in r.species {
		out << Factor{b.member(family, b.species[k].name, reaction_name(b.reactions, i))!, r.nu[j]}
	}
	return out
}

// driving writes k·(driving force): a power law, irreversible, as a
// difference, or as an approach to equilibrium.
fn (mut b Builder) driving(i int, k string, family string, orders []Power, eq Equilibrium, form Form, basis Basis, what string) !string {
	ord := b.orders(family, orders, what)!
	if eq.kind == .irreversible {
		up, down := product(ord)
		return if up == '' { ratio(k, down) } else { ratio('${k}*${up}', down) }
	}
	keq := b.k_eq(i, eq, basis)!
	q := b.quotient(family, i)!
	if form == .approach {
		up, down := product(ord)
		qu, qd := product(q)
		pre := if up == '' { ratio(k, down) } else { ratio('${k}*${up}', down) }
		return '${pre}*(1 - ${ratio(qu, qd)}/${keq})'
	}
	// The difference: the prefactor takes the negative orders, the forward
	// term the positive ones, and the backward term those times Q.
	mut fwd := []Factor{}
	mut pre := []Factor{}
	mut exps := map[string]f64{}
	mut order_of := []string{}
	for f in ord {
		if f.exponent < 0 {
			pre << f
		} else {
			fwd << f
			if f.member !in exps {
				order_of << f.member
			}
			exps[f.member] = exps[f.member] + f.exponent
		}
	}
	for f in q {
		if f.member !in exps {
			order_of << f.member
		}
		exps[f.member] = exps[f.member] + f.exponent
	}
	mut back := []Factor{}
	for m in order_of {
		back << Factor{m, exps[m]}
	}
	_, pdown := product(pre)
	fu, fd := product(fwd)
	bu, bd := product(back)
	prefix := ratio(k, pdown)
	return '${prefix}*(${ratio(fu, fd)} - ${ratio(bu, bd)}/${keq})'
}

// rate writes reaction i's rate law as expression text, defining its k.
fn (mut b Builder) rate(i int, law RateLaw, symbols map[string]string, sites []Site) !string {
	what := reaction_name(b.reactions, i)
	match law {
		Expression {
			return law.text
		}
		PowerLaw {
			k := if law.k_name != '' { law.k_name } else { 'k${i + 1}' }
			b.define(k, law.k.text(), "${what}'s rate constant")!
			return b.driving(i, k, law.basis.family(), law.orders, law.equilibrium, law.form,
				law.basis, what)!
		}
		Lhhw {
			k := if law.k_name != '' { law.k_name } else { 'k${i + 1}' }
			b.define(k, law.k.text(), "${what}'s rate constant")!
			mut theta := symbols[law.site] or {
				return refuse(what, 'no site ${law.site}')
			}
			if terms := law.terms {
				theta = if law.symbol != '' { law.symbol } else { '${theta}_${i + 1}' }
				text := b.terms(terms, "${what}'s own adsorption term")!
				b.define(theta, text, "${what}'s own adsorption term, overriding site ${law.site}'s")!
			}
			d := b.driving(i, k, law.basis.family(), law.orders, law.equilibrium, law.form,
				law.basis, what)!
			n := if law.power == 1 { theta } else { '${theta}^${exponent_text(law.power)}' }
			return '${d}/${n}'
		}
		EleyRideal {
			k := if law.k_name != '' { law.k_name } else { 'k${i + 1}' }
			b.define(k, law.k.text(), "${what}'s rate constant")!
			theta := symbols[law.site] or { return refuse(what, 'no site ${law.site}') }
			site := sites.filter(it.name == law.site)[0]
			a := find(b.species, law.adsorbed) or { return refuse(what, err.msg()) }
			mut ka := ''
			for t in site.terms {
				if t.powers.len == 1 && t.powers[0].exponent == 1
					&& find(b.species, t.powers[0].species) or { -1 } == a {
					ka = if t.name != '' {
						t.name
					} else {
						'K_' + b.ctx.families[0].members[a]
					}
				}
			}
			if ka == '' {
				return refuse(what, '${law.adsorbed} has no term of its own on site ${law.site}, K·a¹, for an Eley–Rideal rate')
			}
			orders := [Power{law.adsorbed, 1}, Power{law.gas, 1}]
			d := b.driving(i, '${k}*${ka}', law.basis.family(), orders, law.equilibrium,
				.approach, law.basis, what)!
			return '${d}/${theta}'
		}
	}
}
