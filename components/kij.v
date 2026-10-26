module components

// KijSource says where a matrix entry came from.
pub enum KijSource {
	file     // the bank's first listing of the pair
	override // given by the caller
	missing  // no data: 0, as ChemSep also defaults it
}

// PairKij records how one pair i < j of a component list got its kᵢⱼ.
pub struct PairKij {
pub:
	i      int
	j      int
	k      f64
	source KijSource
	used   ?Kij  // the listing used, when source is .file
	others []Kij // the bank's listings of this pair that weren't used
}

// KijMatrix is a kᵢⱼ matrix for a list of components, with how each pair
// got its value.
pub struct KijMatrix {
pub:
	bank  string // the file the values came from
	cas   []string
	k     []f64 // nc×nc, row-major, symmetric, zero diagonal: thermo.new_cubic's kij
	pairs []PairKij
}

@[params]
pub struct KijOptions {
pub:
	// overrides set pairs by hand, in place of whatever the bank lists.
	overrides []Kij
}

// matrix builds the kᵢⱼ matrix for components given by CAS number. A pair
// the bank lists takes its first listing (M2 decision: DECHEMA's repeated
// fits for one pair overlap in temperature, so nothing better picks one),
// and the others are kept in the result. A pair with no data gets 0.
pub fn (b KijBank) matrix(cas []string, opts KijOptions) !KijMatrix {
	nc := cas.len
	for i in 0 .. nc {
		for j in 0 .. i {
			if cas[i] == cas[j] {
				return error('${cas[i]} is in the component list twice')
			}
		}
	}
	mut overrides := map[string]Kij{}
	for o in opts.overrides {
		i, j := cas.index(o.cas1), cas.index(o.cas2)
		if i < 0 || j < 0 {
			return error("an override for ${o.cas1}/${o.cas2}, which isn't a pair in the list")
		}
		if i == j {
			return error('an override pairs ${o.cas1} with itself')
		}
		key := pair_key(i, j)
		if key in overrides {
			return error('two overrides for ${o.cas1}/${o.cas2}')
		}
		overrides[key] = o
	}
	mut k := []f64{len: nc * nc}
	mut pairs := []PairKij{}
	for i in 0 .. nc {
		for j in i + 1 .. nc {
			listed := b.for_pair(cas[i], cas[j])
			mut p := PairKij{
				i:      i
				j:      j
				source: .missing
				others: listed
			}
			if o := overrides[pair_key(i, j)] {
				p = PairKij{
					...p
					k:      o.k
					source: .override
				}
			} else if listed.len > 0 {
				p = PairKij{
					...p
					k:      listed[0].k
					source: .file
					used:   listed[0]
					others: listed[1..]
				}
			}
			k[i * nc + j] = p.k
			k[j * nc + i] = p.k
			pairs << p
		}
	}
	return KijMatrix{
		bank:  b.file
		cas:   cas.clone()
		k:     k
		pairs: pairs
	}
}

fn pair_key(i int, j int) string {
	return if i < j { '${i},${j}' } else { '${j},${i}' }
}

// missing returns the pairs that had no data, and so are 0.
pub fn (m KijMatrix) missing() []PairKij {
	return m.pairs.filter(it.source == .missing)
}

// report describes every pair's value: where it came from, and any listings
// of it that weren't used, so the choice can be checked and overridden.
pub fn (m KijMatrix) report() string {
	mut s := ['kij from ${m.bank} for ${m.cas.len} components:']
	for p in m.pairs {
		pair := '${m.cas[p.i]}/${m.cas[p.j]}'
		mut line := match p.source {
			.missing {
				'  ${pair}: 0, no data'
			}
			.override {
				'  ${pair}: ${p.k}, given by hand'
			}
			.file {
				u := p.used or { Kij{} }
				'  ${pair}: ${p.k}, line ${u.line} (${u.comment})'
			}
		}
		if p.others.len > 0 {
			alts := p.others.map('${it.k} on line ${it.line}')
			line += '; also listed: ${alts.join(', ')}'
		}
		s << line
	}
	return s.join('\n') + '\n'
}
