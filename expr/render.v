module expr

// Style is how an expression is written back: plain, in its own syntax,
// which reads back to the same tree; Unicode, for a terminal; or LaTeX, for
// a notebook.
pub enum Style {
	plain
	unicode
	latex
}

// plain is the expression in its own syntax, with no more parentheses than
// it needs. It reads back to the same tree.
pub fn (e &Expr) plain() string {
	return e.form(e.root, .plain)
}

// unicode is the expression as a terminal shows it: k·p_CH₄/(1 + K·p_CH₄)².
pub fn (e &Expr) unicode() string {
	return e.form(e.root, .unicode)
}

// latex is the expression in LaTeX: \frac{k \cdot p_{\mathrm{CH_{4}}}}{…}.
pub fn (e &Expr) latex() string {
	return e.form(e.root, .latex)
}

// precedence: + and − bind loosest, then · and /, then a minus sign, then
// ^; everything else is an atom.
fn precedence(n Node) int {
	return match n.kind {
		.add, .sub { 1 }
		.mul, .div { 2 }
		.neg { 3 }
		.pow { 4 }
		else { 5 }
	}
}

// prec is a node's precedence as a style writes it. Unicode and LaTeX write
// a number in scientific notation as a product, 4.22·10¹⁵, or as a power,
// 10⁻³, so it binds as one: x/(2.5·10⁻¹) needs its parentheses there,
// where plain text's 2.5e-1 is one token.
fn (e &Expr) prec(i int, s Style) int {
	n := e.nodes[i]
	if s != .plain && n.kind == .num && n.text.to_lower().contains('e') {
		mantissa := n.text.to_lower().all_before('e')
		return if mantissa == '1' || mantissa == '1.0' { 4 } else { 2 }
	}
	return precedence(n)
}

fn (e &Expr) form(i int, s Style) string {
	n := e.nodes[i]
	return match n.kind {
		.num { number_form(n, s) }
		.var, .def { name_form(n.text, s) }
		.member { e.member_form(n, s) }
		.constant { constant_form(n.index, s) }
		.neg { e.neg_form(n, s) }
		.add, .sub { e.sum_form(n, s) }
		.mul { e.product_form(n, s) }
		.div { e.quotient_form(n, s) }
		.pow { e.power_form(n, s) }
		.call { e.call_form(n, s) }
		.reading { e.reading_form(n, s) }
		.scaled { e.scaled_form(n, s) }
	}
}

fn paren(text string, s Style) string {
	return if s == .latex { '\\left(${text}\\right)' } else { '(${text})' }
}

// has_unit is whether a node is written with a unit after it: a number's,
// or a group's.
fn has_unit(n Node) bool {
	return (n.kind == .num && n.unit != '') || n.kind == .scaled
}

// operand writes a child, in parentheses if it binds looser than need. In
// Unicode a number with a unit gets them too, in a product or a quotient,
// since its unit's · and / would run into the operator's.
fn (e &Expr) operand(i int, need int, s Style) string {
	text := e.form(i, s)
	n := e.nodes[i]
	if e.prec(i, s) < need || (s == .unicode && has_unit(n)) {
		return paren(text, s)
	}
	return text
}

// reading_form writes a reading in parentheses always, so that it reads as
// one value wherever it stands: (T in degC), (T in °C).
fn (e &Expr) reading_form(n Node, s Style) string {
	inner := e.operand(n.a, 5, s)
	return match s {
		.plain {
			unit := if n.braced { '[${n.unit}]' } else { n.unit }
			'(${inner} in ${unit})'
		}
		.unicode {
			'(${inner} in ${unit_unicode(degrees(n.unit))})'
		}
		.latex {
			'\\left(${inner}\\ \\mathrm{in}\\ \\mathrm{${unit_latex(degrees(n.unit))}}\\right)'
		}
	}
}

// degrees writes degC and degF as °C and °F.
fn degrees(unit string) string {
	return match unit {
		'degC' { '°C' }
		'degF' { '°F' }
		else { unit }
	}
}

// scaled_form writes a group with its unit after it: (…) [%].
fn (e &Expr) scaled_form(n Node, s Style) string {
	inner := e.form(n.a, s)
	return match s {
		.plain { '(${inner}) [${n.unit}]' }
		.unicode { '(${inner}) ${unit_unicode(n.unit)}' }
		.latex { '\\left(${inner}\\right)\\,\\mathrm{${unit_latex(n.unit)}}' }
	}
}

fn (e &Expr) neg_form(n Node, s Style) string {
	child := e.nodes[n.a]
	inner := e.form(n.a, s)
	minus := if s == .unicode { '−' } else { '-' }
	// In Unicode a number with a unit keeps its parentheses here too, since
	// the minus sign may itself be divided: −(240.1 kJ/mol)/(R·T).
	if e.prec(n.a, s) < 3 || child.kind == .neg || (s == .unicode && has_unit(child)) {
		return minus + paren(inner, s)
	}
	return minus + inner
}

fn (e &Expr) sum_form(n Node, s Style) string {
	left := e.form(n.a, s)
	right_node := e.nodes[n.b]
	mut right := e.form(n.b, s)
	// A sum on the right keeps its parentheses: a − (b + c) needs them, and
	// a + (b + c), though equal to a + b + c, is another tree, which the
	// rendering must give back. A minus sign on the right gets them to be
	// read easily, a + (−b).
	if e.prec(n.b, s) <= 1 || right_node.kind == .neg {
		right = paren(right, s)
	}
	op := match true {
		n.kind == .add { ' + ' }
		s == .unicode { ' − ' }
		else { ' - ' }
	}
	return left + op + right
}

fn (e &Expr) product_form(n Node, s Style) string {
	left := e.operand(n.a, 2, s)
	// A product or quotient on the right keeps its parentheses: a·(b/c) is
	// another tree than a·b/c. In LaTeX a quotient is a \frac, grouped
	// already.
	right_node := e.nodes[n.b]
	need := if s == .latex && right_node.kind == .div { 2 } else { 3 }
	mut right := e.operand(n.b, need, s)
	if right_node.kind == .neg {
		right = paren(e.form(n.b, s), s)
	}
	op := match s {
		.plain { ' * ' }
		.unicode { '·' }
		.latex { ' \\cdot ' }
	}
	return left + op + right
}

fn (e &Expr) quotient_form(n Node, s Style) string {
	if s == .latex {
		return '\\frac{${e.form(n.a, s)}}{${e.form(n.b, s)}}'
	}
	left := e.operand(n.a, 2, s)
	// a/(b·c) and a/(b/c) need their parentheses: the right binds tighter.
	mut right := e.operand(n.b, 3, s)
	return left + if s == .plain { ' / ' } else { '/' } + right
}

fn (e &Expr) power_form(n Node, s Style) string {
	base_node := e.nodes[n.a]
	mut base := e.form(n.a, s)
	// a^b^c is a^(b^c), so a power as the base needs parentheses; and a
	// number's unit would take the ^ for its own, so a number with a unit
	// does too, as does a call in LaTeX, where \exp(x)^{2} reads badly.
	if e.prec(n.a, s) <= 4 || has_unit(base_node) || (s == .latex && base_node.kind == .call) {
		base = paren(base, s)
	}
	exp_node := e.nodes[n.b]
	if s == .latex {
		return '${base}^{${e.form(n.b, s)}}'
	}
	if s == .unicode {
		if sup := superscript_of(e, exp_node) {
			return base + sup
		}
	}
	exponent := e.form(n.b, s)
	if e.prec(n.b, s) < 4 || (exp_node.kind == .num && exp_node.unit != '') {
		return '${base}^${paren(exponent, s)}'
	}
	return '${base}^${exponent}'
}

// superscript_of is a whole power as superscripts, ² or ⁻¹, if it is one.
fn superscript_of(e &Expr, n Node) ?string {
	if n.kind == .num && n.unit == '' && is_digits(n.text) {
		return superscript_digits(n.text)
	}
	if n.kind == .neg {
		inner := e.nodes[n.a]
		if inner.kind == .num && inner.unit == '' && is_digits(inner.text) {
			return '⁻' + superscript_digits(inner.text)
		}
	}
	return none
}

fn is_digits(s string) bool {
	if s == '' {
		return false
	}
	for c in s {
		if !c.is_digit() {
			return false
		}
	}
	return true
}

const superscript_table = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']!
const subscript_table = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']!

fn superscript_digits(s string) string {
	mut out := ''
	for c in s {
		out += if c == `-` {
			'⁻'
		} else if c.is_digit() {
			superscript_table[c - `0`]
		} else {
			c.ascii_str()
		}
	}
	return out
}

fn subscript_digits(s string) string {
	mut out := ''
	for c in s.runes() {
		out += if c >= `0` && c <= `9` { subscript_table[int(c - `0`)] } else { c.str() }
	}
	return out
}

fn (e &Expr) call_form(n Node, s Style) string {
	arg := e.form(n.a, s)
	name := functions[n.index]
	if s == .plain {
		return '${name}(${arg})'
	}
	if s == .latex {
		return match name {
			'exp' { '\\exp\\left(${arg}\\right)' }
			'ln' { '\\ln\\left(${arg}\\right)' }
			'log10' { '\\log_{10}\\left(${arg}\\right)' }
			'sqrt' { '\\sqrt{${arg}}' }
			else { '\\left|${arg}\\right|' }
		}
	}
	return match name {
		'sqrt' {
			if e.prec(n.a, s) == 5 && e.nodes[n.a].unit == '' { '√${arg}' } else { '√(${arg})' }
		}
		'abs' {
			'|${arg}|'
		}
		'log10' {
			'log₁₀(${arg})'
		}
		else {
			'${name}(${arg})'
		}
	}
}

// greek is the Greek letters a name may spell out, rendered as letters.
const greek = {
	'alpha':   ['α', '\\alpha']
	'beta':    ['β', '\\beta']
	'gamma':   ['γ', '\\gamma']
	'delta':   ['δ', '\\delta']
	'epsilon': ['ε', '\\varepsilon']
	'eta':     ['η', '\\eta']
	'theta':   ['θ', '\\theta']
	'kappa':   ['κ', '\\kappa']
	'lambda':  ['λ', '\\lambda']
	'mu':      ['μ', '\\mu']
	'nu':      ['ν', '\\nu']
	'rho':     ['ρ', '\\rho']
	'sigma':   ['σ', '\\sigma']
	'tau':     ['τ', '\\tau']
	'phi':     ['φ', '\\phi']
	'chi':     ['χ', '\\chi']
	'psi':     ['ψ', '\\psi']
	'omega':   ['ω', '\\omega']
	'Gamma':   ['Γ', '\\Gamma']
	'Delta':   ['Δ', '\\Delta']
	'Theta':   ['Θ', '\\Theta']
	'Phi':     ['Φ', '\\Phi']
	'Omega':   ['Ω', '\\Omega']
}

// name_form writes a name: a Greek letter where it spells one, its trailing
// digits as a subscript (k1 is k₁), and what follows a _ as a subscript too.
fn name_form(name string, s Style) string {
	if s == .plain {
		return name
	}
	base_full := name.all_before('_')
	suffix := if name.contains('_') { name.all_after('_') } else { '' }
	mut letters := base_full
	mut digits := ''
	for letters.len > 1 && letters[letters.len - 1].is_digit() {
		digits = letters[letters.len - 1..] + digits
		letters = letters[..letters.len - 1]
	}
	shown := if g := greek[letters] { g[if s == .unicode { 0 } else { 1 }] } else { letters }
	if s == .unicode {
		mut out := shown + subscript_digits(digits)
		if suffix != '' {
			out += if is_digits(suffix) {
				subscript_digits(suffix)
			} else {
				'_' + subscript_digits(suffix)
			}
		}
		return out
	}
	mut sub := []string{}
	if digits != '' {
		sub << digits
	}
	if suffix != '' {
		sub << if suffix.len == 1 || is_digits(suffix) {
			suffix
		} else {
			'\\mathrm{${formula_latex(suffix)}}'
		}
	}
	return if sub.len > 0 { '${shown}_{${sub.join(',')}}' } else { shown }
}

// formula_latex writes a formula's digits as subscripts: CH4 is CH_{4}.
fn formula_latex(s string) string {
	mut out := ''
	mut i := 0
	for i < s.len {
		if s[i].is_digit() {
			mut j := i
			for j < s.len && s[j].is_digit() {
				j++
			}
			out += '_{${s[i..j]}}'
			i = j
		} else {
			out += s[i].ascii_str()
			i++
		}
	}
	return out
}

fn (e &Expr) member_form(n Node, s Style) string {
	family := e.family_name(n)
	return match s {
		.plain {
			if n.text.len > 0 && n.text.bytes().all(it.is_letter() || it.is_digit() || it == `_`) {
				'${family}[${n.text}]'
			} else {
				"${family}['${n.text}']"
			}
		}
		.unicode {
			'${name_form(family, s)}_${subscript_digits(n.text)}'
		}
		.latex {
			'${name_form(family, s)}_{\\mathrm{${formula_latex(n.text)}}}'
		}
	}
}

// family_name is the name of the family a member belongs to.
fn (e &Expr) family_name(n Node) string {
	return e.families[n.a] or { '?' }
}

fn constant_form(i int, s Style) string {
	c := constants[i]
	if s == .plain {
		return c.name
	}
	return match c.name {
		'pi' {
			if s == .unicode { 'π' } else { '\\pi' }
		}
		'p_std' {
			if s == .unicode { 'p°' } else { 'p^{\\circ}' }
		}
		'g_n' {
			if s == .unicode { 'gₙ' } else { 'g_{n}' }
		}
		else { name_form(c.name, s) }
	}
}

fn number_form(n Node, s Style) string {
	num := match s {
		.plain { n.text }
		.unicode { scientific(n.text, '·10', true) }
		.latex { scientific(n.text, ' \\times 10', false) }
	}
	if n.unit == '' {
		return num
	}
	return match s {
		.plain {
			if n.braced { '${num} [${n.unit}]' } else { '${num} ${n.unit}' }
		}
		.unicode { '${num} ${unit_unicode(n.unit)}' }
		.latex { '${num}\\,\\mathrm{${unit_latex(n.unit)}}' }
	}
}

// scientific writes 4.22e15 as 4.22·10¹⁵ or 4.22 \times 10^{15}, and 1e-3
// as 10⁻³ or 10^{-3}.
fn scientific(text string, times string, unicode bool) string {
	lower := text.to_lower()
	if !lower.contains('e') {
		return text
	}
	mantissa := lower.all_before('e')
	mut power := lower.all_after('e').trim_left('+')
	negative := power.starts_with('-')
	power = power.trim_left('-').trim_left('0')
	if power == '' {
		power = '0'
	}
	if negative {
		power = '-' + power
	}
	ten := if unicode { '10' + superscript_digits(power) } else { '10^{${power}}' }
	if mantissa == '1' || mantissa == '1.0' {
		return ten
	}
	return mantissa + times + ten.all_after('10')
}

// unit_unicode writes a unit's * as · and its whole powers as superscripts.
fn unit_unicode(u string) string {
	mut out := ''
	runes := u.runes()
	mut i := 0
	for i < runes.len {
		c := runes[i]
		if c == `*` {
			out += '·'
			i++
		} else if c == `^` && i + 1 < runes.len && (runes[i + 1] == `-` || is_digit(runes[i + 1])) {
			mut j := i + 1
			if runes[j] == `-` {
				j++
			}
			for j < runes.len && is_digit(runes[j]) {
				j++
			}
			if j < runes.len && runes[j] == `.` {
				out += runes[i..j].string()
			} else {
				out += superscript_digits(runes[i + 1..j].string())
			}
			i = j
		} else if is_digit(c) && i > 0 && is_letter(runes[i - 1]) {
			mut j := i
			for j < runes.len && is_digit(runes[j]) {
				j++
			}
			out += superscript_digits(runes[i..j].string())
			i = j
		} else {
			out += c.str()
			i++
		}
	}
	return out
}

// unit_latex writes a unit inside \mathrm{}: * as \cdot, powers raised.
fn unit_latex(u string) string {
	mut out := ''
	runes := u.runes()
	mut i := 0
	for i < runes.len {
		c := runes[i]
		if c == `*` || c == `·` || c == `⋅` {
			out += '\\cdot '
			i++
		} else if c == `^` {
			mut j := i + 1
			if j < runes.len && runes[j] == `(` {
				for j < runes.len && runes[j] != `)` {
					j++
				}
				j++
			} else {
				if j < runes.len && (runes[j] == `-` || runes[j] == `+`) {
					j++
				}
				for j < runes.len && (is_digit(runes[j]) || runes[j] == `.`) {
					j++
				}
			}
			out += '^{${runes[i + 1..j].string().trim('()')}}'
			i = j
		} else if superscript_digit_value(c) >= 0 || c == `⁻` {
			mut j := i
			mut power := ''
			for j < runes.len && (superscript_digit_value(runes[j]) >= 0 || runes[j] == `⁻`) {
				power += if runes[j] == `⁻` {
					'-'
				} else {
					superscript_digit_value(runes[j]).str()
				}
				j++
			}
			out += '^{${power}}'
			i = j
		} else if is_digit(c) && i > 0 && is_letter(runes[i - 1]) {
			mut j := i
			for j < runes.len && is_digit(runes[j]) {
				j++
			}
			out += '^{${runes[i..j].string()}}'
			i = j
		} else if c == `µ` || c == `μ` {
			out += '\\mu '
			i++
		} else if c == `%` {
			out += '\\%'
			i++
		} else if c == `°` {
			out += '{}^{\\circ}'
			i++
		} else {
			out += c.str()
			i++
		}
	}
	return out
}

fn superscript_digit_value(c rune) int {
	for i, t in superscript_table {
		if c.str() == t {
			return i
		}
	}
	return -1
}
