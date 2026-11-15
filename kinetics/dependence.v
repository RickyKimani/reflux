module kinetics

// Share is one reaction's part in a combination: factor times it.
pub struct Share {
pub:
	reaction int
	factor   Ratio
}

// Combination is a reaction that adds nothing to those before it: its
// equation is the sum of theirs, each times its factor.
pub struct Combination {
pub:
	reaction int
	of       []Share
}

// Dependence is a set of reactions' independence: its rank, how many are
// independent, and each reaction that is a combination of earlier ones.
pub struct Dependence {
pub:
	rank         int
	combinations []Combination
}

// dependence finds which of a set of reactions are combinations of those
// before them, by exact elimination on their coefficients, one reaction at
// a time: each is reduced by the independent ones so far, and if nothing
// is left, the factors that reduced it say what it's made of. A set with
// combinations is fine for kinetics, each reaction having its own rate,
// but not for an equilibrium over reactions, which needs them independent.
pub fn dependence(reactions []Reaction, species int) !Dependence {
	for r in reactions {
		if r.species.any(it >= species) {
			return error('${r.name} has a species beyond the ${species} given')
		}
	}
	n := reactions.len
	mut rows := [][]Ratio{} // the independent reactions, reduced
	mut pivots := []int{}
	mut makeup := [][]Ratio{} // each reduced row as a combination of reactions
	mut combinations := []Combination{}
	for k, r in reactions {
		mut row := []Ratio{len: species, init: zero}
		for j, s in r.species {
			row[s] = r.exact[j]
		}
		mut made := []Ratio{len: n, init: zero}
		made[k] = one
		for b, p in pivots {
			if is_zero(row[p]) {
				continue
			}
			f := div(row[p], rows[b][p])!
			for j in 0 .. species {
				row[j] = sub(row[j], f.mul(rows[b][j])!)!
			}
			for j in 0 .. n {
				made[j] = sub(made[j], f.mul(makeup[b][j])!)!
			}
		}
		p := first_nonzero(row)
		if p >= 0 {
			rows << row
			pivots << p
			makeup << made
			continue
		}
		// 0 = Σ made[j]·νⱼ with made[k] = 1, so νₖ = Σ −made[j]·νⱼ.
		mut of := []Share{}
		for j in 0 .. k {
			if !is_zero(made[j]) {
				of << Share{j, made[j].neg()}
			}
		}
		combinations << Combination{k, of}
	}
	return Dependence{rows.len, combinations}
}

fn first_nonzero(row []Ratio) int {
	for j, x in row {
		if !is_zero(x) {
			return j
		}
	}
	return -1
}

// describe writes a combination with the reactions' names, "overall =
// reforming + shift" or "r3 = 2 × r1 − r2".
pub fn (c Combination) describe(reactions []Reaction) string {
	mut s := '${reactions[c.reaction].name} ='
	for i, sh in c.of {
		negative := sh.factor.num < 0
		a := if negative { sh.factor.neg() } else { sh.factor }
		sign := if i == 0 {
			if negative { ' −' } else { '' }
		} else if negative {
			' −'
		} else {
			' +'
		}
		times := if a == one { '' } else { '${decimal(a)} × ' }
		s += '${sign} ${times}${reactions[sh.reaction].name}'
	}
	return s
}
