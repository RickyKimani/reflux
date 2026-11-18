module expr

// Relation is how an equation's two sides are to compare.
pub enum Relation {
	eq // =
	le // <=, ≤
	ge // >=, ≥
}

// Equation is two expressions and how they compare: a spec's target,
// "drum_vapour.y[CH4] = 5 %" (M8's D5), or an optimisation's constraint,
// "product.T <= 1100 K" (D13). Both sides have one dimension. Each side's
// nodes have their line and column in the whole text.
pub struct Equation {
pub:
	text     string
	lhs      Expr
	rhs      Expr
	relation Relation
	line     int // where the relation is, from 1, in characters
	column   int
}

// equation reads an equation: two expressions with =, <= or >= between
// them (≤ and ≥ too), bound in ctx with these definitions, whose sides must
// have the same units.
pub fn equation(ctx Context, defs &Definitions, text string) !Equation {
	runes := text.runes()
	p := new_parser(ctx, defs.names, text)
	at, width, relation := find_relation(&p, runes)!
	line, column := p.where(at)
	mut left := runes.clone()
	mut right := runes.clone()
	for k in 0 .. runes.len {
		if k >= at && runes[k] != `\n` {
			left[k] = ` `
		}
		if k < at + width && runes[k] != `\n` {
			right[k] = ` `
		}
	}
	lhs := side(ctx, defs, text, left.string(), at, 'left')!
	rhs := side(ctx, defs, text, right.string(), at, 'right')!
	if !lhs.dim.eq(rhs.dim) {
		sign := runes[at..at + width].string()
		return ExprError{
			text:   text
			line:   line
			column: column
			reason: 'the two sides of ${sign} must have the same units, but the left is ${describe(lhs.dim)} and the right ${describe(rhs.dim)}${reading_hint(lhs.dim, rhs.dim)}${reading_hint(rhs.dim, lhs.dim)}'
		}
	}
	return Equation{
		text:     text
		lhs:      lhs
		rhs:      rhs
		relation: relation
		line:     line
		column:   column
	}
}

// side parses one side of an equation, the rest of the text blanked so that
// its places are the whole text's.
fn side(ctx Context, defs &Definitions, text string, masked string, at int, which string) !Expr {
	if masked.trim_space() == '' {
		p := new_parser(ctx, defs.names, text)
		line, column := p.where(at)
		return ExprError{
			text:   text
			line:   line
			column: column
			reason: 'the ${which} side of the equation is empty'
		}
	}
	e := parse(ctx, defs, masked) or {
		if err is ExprError {
			x := err as ExprError
			return ExprError{
				...x
				text: text
			}
		}
		return err
	}
	return Expr{
		...e
		text: text
	}
}

// find_relation is where the equation's relation is, its width in
// characters, and which it is. It looks outside brackets and quotes, and
// refuses none, two, or a strict inequality.
fn find_relation(p &Parser, runes []rune) !(int, int, Relation) {
	mut depth := 0
	mut quote := rune(0)
	mut found := -1
	mut width := 0
	mut relation := Relation.eq
	mut k := 0
	for k < runes.len {
		c := runes[k]
		if quote != 0 {
			if c == quote {
				quote = 0
			}
			k++
			continue
		}
		next := if k + 1 < runes.len { runes[k + 1] } else { rune(0) }
		mut w := 0
		mut r := Relation.eq
		match c {
			`'`, `"` {
				quote = c
			}
			`(`, `[` {
				depth++
			}
			`)`, `]` {
				depth--
			}
			`=` {
				if next == `=` {
					return p.fail(k, 'an == here: write = once')
				}
				if next == `<` || next == `>` {
					return p.fail(k, 'write ${next}= with the ${next} first')
				}
				w = 1
			}
			`<`, `>` {
				if next != `=` {
					return p.fail(k, 'a strict ${c}: write ${c}=, since a bound is met where its two sides are equal')
				}
				w = 2
				r = if c == `<` { Relation.le } else { Relation.ge }
			}
			`≤`, `≥` {
				w = 1
				r = if c == `≤` { Relation.le } else { Relation.ge }
			}
			else {}
		}
		if w > 0 {
			if depth != 0 {
				return p.fail(k, 'a relation inside brackets: an equation has it between its two sides')
			}
			if found >= 0 {
				return p.fail(k, 'a second relation: an equation has one, between its two sides')
			}
			found, width, relation = k, w, r
			k += w
			continue
		}
		k++
	}
	if found < 0 {
		return p.fail(0, 'an equation needs =, <= or >= between its two sides')
	}
	return found, width, relation
}

// plain is the equation in its own syntax, which reads back to the same
// sides.
pub fn (q &Equation) plain() string {
	sign := match q.relation {
		.eq { '=' }
		.le { '<=' }
		.ge { '>=' }
	}
	return '${q.lhs.plain()} ${sign} ${q.rhs.plain()}'
}

// unicode is the equation as a terminal shows it.
pub fn (q &Equation) unicode() string {
	sign := match q.relation {
		.eq { '=' }
		.le { '≤' }
		.ge { '≥' }
	}
	return '${q.lhs.unicode()} ${sign} ${q.rhs.unicode()}'
}

// latex is the equation in LaTeX.
pub fn (q &Equation) latex() string {
	sign := match q.relation {
		.eq { '=' }
		.le { '\\le' }
		.ge { '\\ge' }
	}
	return '${q.lhs.latex()} ${sign} ${q.rhs.latex()}'
}
