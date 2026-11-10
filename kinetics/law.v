module kinetics

import math

// The structured forms of a rate law (D4, D6, docs/kinetics.md). Each
// keeps its structure, for a site shared by every reaction on it and for
// its own checks, and writes the expression a user would type for it, which
// expr parses: so a structured rate and the same rate typed give the same
// tree, and render, evaluate and differentiate alike.
//
// Numbers are kept as written, with their units, "4.225e15 kmol*bar^0.5/(kg*h)",
// so that a rate renders as its paper prints it.

// Arrhenius is a rate constant, A·exp(−Eₐ/RT), or centred on T_ref,
// k_ref·exp(−Eₐ/R·(1/T − 1/T_ref)), which fits better conditioned.
pub struct Arrhenius {
pub:
	a     string @[required] // A, or k at T_ref: a number with its unit
	ea    string @[required] // Eₐ, an energy per mole
	t_ref string // '' for A·exp(−Eₐ/RT)
}

// VantHoff is an adsorption or equilibrium constant, A·exp(−ΔH/RT), or
// centred, K_ref·exp(−ΔH/R·(1/T − 1/T_ref)). An adsorption's ΔH is below
// zero, so its K falls as T rises.
pub struct VantHoff {
pub:
	a     string @[required]
	dh    string @[required]
	t_ref string
}

// Fixed is a constant that doesn't depend on T, with its unit.
pub struct Fixed {
pub:
	value string @[required]
}

// Written is a constant as any expression, in T and the built-in constants
// and definitions: exp(30.114 − 26830 K/T)·1 bar^2, say.
pub struct Written {
pub:
	text string @[required]
}

// TempFn is how a constant depends on temperature.
pub type TempFn = Arrhenius | Fixed | VantHoff | Written

// exponential writes c·exp(−e/(R·T)), or centred on t_ref, with a minus
// sign written out of the energy where it has one: exp(70.65 kJ/mol/(R·T))
// for an adsorption's ΔH of −70.65 kJ/mol, as papers print it.
fn exponential(c string, e string, t_ref string) string {
	energy := e.trim_space()
	negative := energy.starts_with('-') || energy.starts_with('−')
	magnitude := if negative { energy.runes()[1..].string().trim_space() } else { energy }
	sign := if negative { '' } else { '-' }
	if t_ref.trim_space() == '' {
		return '${c} * exp(${sign}${magnitude} / (R*T))'
	}
	return '${c} * exp(${sign}${magnitude} / R * (1/T - 1/(${t_ref})))'
}

// text is the constant as expression text.
pub fn (f TempFn) text() string {
	return match f {
		Arrhenius { exponential(f.a, f.ea, f.t_ref) }
		VantHoff { exponential(f.a, f.dh, f.t_ref) }
		Fixed { f.value }
		Written { f.text }
	}
}

// Power is a species raised to an exponent: an order, or a site term's
// part. species is the species as written, by any name find takes.
pub struct Power {
pub:
	species  string @[required]
	exponent f64    @[required]
}

// Monomial is one term of a site's adsorption term, K(T)·Π aᵢ^eᵢ: K_A·p_A,
// √K_H2·p_H2^0.5, or K·p_H2O/p_H2. name names its K, as a definition.
pub struct Monomial {
pub:
	name   string
	k      TempFn  @[required]
	powers []Power @[required]
}

// Site is a catalyst's site: its adsorption term, θ = 1 + Σ terms, which
// every reaction on it shares, each to its own power. name is the site's,
// and symbol names θ as a definition (theta_ and the name if not given).
pub struct Site {
pub:
	name   string @[required]
	symbol string
	terms  []Monomial
}

// Basis is what a rate's activities are: partial pressures p, fugacities
// f, concentrations c or mole fractions y.
pub enum Basis {
	pressure
	fugacity
	concentration
	fraction
}

// family is the basis's family in a rate's context.
fn (b Basis) family() string {
	return match b {
		.pressure { 'p' }
		.fugacity { 'f' }
		.concentration { 'c' }
		.fraction { 'y' }
	}
}

// Per is what a rate is per: a mass of catalyst or a volume of reactor.
pub enum Per {
	mass
	volume
}

// Equilibrium is how a reversible rate knows its K: from the Gibbs
// energies (D7), K at 1 bar from each species' ideal-gas H and S, or from
// a correlation the paper fitted, in the rate's activities (bar², say).
// An irreversible rate has none.
pub enum Reversibility {
	irreversible
	thermo
	correlation
}

pub struct Equilibrium {
pub:
	kind        Reversibility
	correlation TempFn = Fixed{'1'}
	name        string // the correlation's definition, K and the reaction's number if not given
}

// Form is how a reversible driving force is written: as a difference,
// k·(Π aᵢ^αᵢ − Π aᵢ^(αᵢ+νᵢ)/K), as Xu and Froment print theirs, or as an
// approach to equilibrium, k·Π aᵢ^αᵢ·(1 − Q/K), as Vanden Bussche and
// Froment print theirs. They're the same function; each renders as its
// paper does.
pub enum Form {
	difference
	approach
}

// PowerLaw is r = k(T)·Π aᵢ^αᵢ, times (1 − Q/K) or as a difference where
// it's reversible.
pub struct PowerLaw {
pub:
	k           TempFn @[required]
	k_name      string // k and the reaction's number if not given
	orders      []Power @[required]
	equilibrium Equilibrium
	form        Form = .approach
	basis       Basis
	per         Per
}

// Lhhw is a Langmuir–Hinshelwood–Hougen–Watson rate: a power law over a
// site's adsorption term to a power, r = k·(driving force)/θⁿ. terms, if
// given, override the site's for this reaction alone, and are marked as
// an override.
pub struct Lhhw {
pub:
	k           TempFn @[required]
	k_name      string
	orders      []Power @[required]
	equilibrium Equilibrium
	form        Form = .approach
	basis       Basis
	per         Per
	site        string @[required]
	power       f64 = 1
	terms       ?[]Monomial
	symbol      string // the override's θ, as a definition
}

// EleyRideal is a gas-phase species reacting with an adsorbed one, r =
// k·K_A·a_A·a_B/θ, times (1 − Q/K) where it's reversible: adsorbed is A,
// whose term of its own on the site gives K_A, and gas is B.
pub struct EleyRideal {
pub:
	k           TempFn @[required]
	k_name      string
	adsorbed    string @[required]
	gas         string @[required]
	equilibrium Equilibrium
	basis       Basis
	per         Per
	site        string @[required]
}

// Expression is a rate as any expression in the rate's context, with the
// same checks as the structured forms.
pub struct Expression {
pub:
	text  string @[required]
	basis Basis
	per   Per
}

// RateLaw is a reaction's rate law. A reactor takes reactions and doesn't
// care which form each has.
pub type RateLaw = EleyRideal | Expression | Lhhw | PowerLaw

// exponent_text writes an exponent as reflux reads it back: a whole number
// without a point, any other as the shortest double that reads back.
fn exponent_text(x f64) string {
	if x == math.trunc(x) && math.abs(x) < 1e15 {
		return '${i64(x)}'
	}
	return '${x}'
}

// Factor is a family member to a power, written as text by product.
struct Factor {
	member   string // the member as text, p[CH4]
	exponent f64
}

// product writes Π factors: those with exponents above zero multiplied,
// those below divided, each to its power; 1 if there are none above.
fn product(factors []Factor) (string, string) {
	mut up := []string{}
	mut down := []string{}
	for f in factors {
		if f.exponent == 0 {
			continue
		}
		a := math.abs(f.exponent)
		w := if a == 1 { f.member } else { '${f.member}^${exponent_text(a)}' }
		if f.exponent > 0 {
			up << w
		} else {
			down << w
		}
	}
	return up.join('*'), down.join('*')
}

// ratio writes a numerator over a denominator, each a product's text.
fn ratio(up string, down string) string {
	n := if up == '' { '1' } else { up }
	if down == '' {
		return n
	}
	return if down.contains('*') { '${n}/(${down})' } else { '${n}/${down}' }
}
