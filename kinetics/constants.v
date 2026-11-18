module kinetics

import expr

// NamedConstant is one of a set of rates' named constants, as build names
// them: a rate constant (k1, or the law's own k_name), a site term's
// constant (K_a, or K_ and its species), or a fitted K (K1, or its own
// name); and how it depends on T. A fit names its parameters by them, k1.A
// or K_a.dH (M8's D14).
pub struct NamedConstant {
pub:
	name  string
	value TempFn
	what  string
}

// constants are the rates' named constants, each once, in the order build
// defines them: the sites' terms, then each rate's constant, its own terms
// and its K.
pub fn (k &Kinetics) constants() []NamedConstant {
	mut out := []NamedConstant{}
	mut seen := map[string]bool{}
	for site in k.sites {
		for t in site.terms {
			name := k.term_name(t)
			if name !in seen {
				seen[name] = true
				out << NamedConstant{name, t.k, 'a term of site ${site.name}'}
			}
		}
	}
	for i, law in k.laws {
		r := k.rated[i]
		what := reaction_name(k.reactions, r)
		match law {
			PowerLaw, EleyRideal {
				name := if law.k_name != '' { law.k_name } else { 'k${r + 1}' }
				out << NamedConstant{name, law.k, "${what}'s rate constant"}
			}
			Lhhw {
				name := if law.k_name != '' { law.k_name } else { 'k${r + 1}' }
				out << NamedConstant{name, law.k, "${what}'s rate constant"}
				if terms := law.terms {
					for t in terms {
						tn := k.term_name(t)
						if tn !in seen {
							seen[tn] = true
							out << NamedConstant{tn, t.k, "a term of ${what}'s own adsorption term"}
						}
					}
				}
			}
			Expression {}
		}
		if eq := equilibrium_of(law) {
			if eq.kind == .correlation {
				name := if eq.name != '' { eq.name } else { 'K${r + 1}' }
				out << NamedConstant{name, eq.correlation, "${what}'s K"}
			}
		}
	}
	return out
}

// term_name is a site term's constant's name, as build gives it.
fn (k &Kinetics) term_name(t Monomial) string {
	if t.name != '' {
		return t.name
	}
	i := find(k.species, t.powers[0].species) or { return 'K_${t.powers[0].species}' }
	return 'K_' + k.ctx.families[0].members[i]
}

// with_constants is the same rates built again with the constants named in
// replace changed, and more variables in the context after T and P, which a
// replacement may use: a fit writes its parameters in as variables, so that
// the tape differentiates by them (M8's D14).
pub fn (k &Kinetics) with_constants(replace map[string]TempFn, extra []expr.Variable) !Kinetics {
	mut sites := []Site{}
	for site in k.sites {
		sites << Site{
			...site
			terms: k.replaced_terms(site.terms, replace)
		}
	}
	mut laws := []Rated{}
	for i, law in k.laws {
		r := k.rated[i]
		key := 'k${r + 1}'
		laws << Rated{
			reaction: r
			law:      match law {
				PowerLaw {
					RateLaw(PowerLaw{
						...law
						k:           replace[if law.k_name != '' { law.k_name } else { key }] or {
							law.k
						}
						equilibrium: replaced_eq(law.equilibrium, r, replace)
					})
				}
				EleyRideal {
					RateLaw(EleyRideal{
						...law
						k:           replace[if law.k_name != '' { law.k_name } else { key }] or {
							law.k
						}
						equilibrium: replaced_eq(law.equilibrium, r, replace)
					})
				}
				Lhhw {
					mut own := ?[]Monomial(none)
					if terms := law.terms {
						own = k.replaced_terms(terms, replace)
					}
					RateLaw(Lhhw{
						...law
						k:           replace[if law.k_name != '' { law.k_name } else { key }] or {
							law.k
						}
						equilibrium: replaced_eq(law.equilibrium, r, replace)
						terms:       own
					})
				}
				Expression {
					RateLaw(law)
				}
			}
		}
	}
	return build_with(k.species, k.reactions, sites, laws, k.own, extra)
}

fn (k &Kinetics) replaced_terms(terms []Monomial, replace map[string]TempFn) []Monomial {
	return terms.map(Monomial{
		...it
		k: replace[k.term_name(it)] or { it.k }
	})
}

fn replaced_eq(eq Equilibrium, r int, replace map[string]TempFn) Equilibrium {
	if eq.kind != .correlation {
		return eq
	}
	name := if eq.name != '' { eq.name } else { 'K${r + 1}' }
	return Equilibrium{
		...eq
		correlation: replace[name] or { eq.correlation }
	}
}
