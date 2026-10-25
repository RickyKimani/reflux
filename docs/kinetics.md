# Reaction kinetics design

## What went wrong in DWSIM

Setting up an LHHW model for a heterogeneous catalytic reactor meant:

- typing each rate expression into a text box, with no way to check that
  the program understood it the way the paper meant it;
- entering the adsorption (denominator) terms one by one for every
  reaction, and leaving out by hand the species that weren't involved.

reflux treats kinetics as **structured data that the program understands**,
not as strings it evaluates. The structure is what makes checking possible.

## The LHHW form, as data

Every Langmuir–Hinshelwood–Hougen–Watson rate has the same shape:

    r = k(T) · (driving force) / (adsorption term)ⁿ

| Part | General form | Stored as |
|---|---|---|
| kinetic factor | k(T) = k_ref · exp(−Eₐ/R · (1/T − 1/T_ref)) | Arrhenius, centred on T_ref (better conditioned when fitting) |
| driving force | Π aᵢ^αᵢ · (1 − Q/K) | forward orders, plus reversibility (see below) |
| adsorption term | (1 + Σ terms)ⁿ, where each term is Kⱼ(T)·Π aₖ^eₖ | a list of monomials; the power n per reaction |

Here aᵢ is the activity basis: partial pressure, fugacity or concentration,
chosen per reaction and stated explicitly.

The monomial form covers every adsorption term in the literature without
special cases:

| Physical picture | Term | Monomial |
|---|---|---|
| molecular adsorption | K_A p_A | K_A · p_A¹ |
| dissociative adsorption | √(K_H₂ p_H₂) | K′ · p_H₂^0.5 |
| quasi-equilibrated surface step | K p_H₂O / p_H₂ | K · p_H₂O¹ · p_H₂⁻¹ |

### The adsorption term belongs to the site, not to the reaction

In a site-balance derivation, the denominator comes from counting what
occupies the catalyst's active sites. It contains a term for **every
species or surface intermediate that holds a significant share of those
sites**, whether or not it takes part in the reaction at hand. Reactions on
the same site therefore share the same term, though each reaction may raise
it to its own power.

Methanol synthesis on Cu/ZnO (Vanden Bussche & Froment, 1996) is the
classic example. Both reactions share one site term:

    θ = 1 + Kₐ·p_H₂O/p_H₂ + √(K_b·p_H₂) + K_c·p_H₂O

    r_MeOH = k₁ · p_CO₂ · p_H₂ · (1 − Q₁/K_eq,1) / θ³
    r_RWGS = k₂ · p_CO₂ · (1 − Q₂/K_eq,2) / θ

H₂O and H₂ appear in the RWGS denominator. CO₂, CO and CH₃OH appear in
neither, because the derivation keeps only the most abundant surface
species. Whether a species appears depends on *site coverage*, not on
*participation in the reaction*.

So in reflux the term is defined **once per site**, and every reaction on
that site uses it automatically. Mechanisms with more than one kind of site,
or a published model whose denominators really do differ by reaction, get
one site per denominator, or a per-reaction override. The override is the
exception, and it is marked as one.

### Reversibility comes from thermodynamics

The driving force is written as Π aᵢ^αᵢ · (1 − Q/K), where Q is the
reaction quotient. By default K(T) comes from the Gibbs energies of
formation in the component data, so it's consistent with the rest of the
flowsheet by construction. You can still give K(T) as a correlation when a
paper fitted its own, and reflux then reports how far that correlation is
from the thermodynamic value.

## Input format

TOML, like flowsheets. The numbers below are placeholders, not the
published parameters.

```toml
[[site]]
name = "Cu"
terms = [
  { k = { k_ref = 1.0, dh = -40e3, t_ref = 500.0 }, a = { H2O = 1.0, H2 = -1.0 } },
  { k = { k_ref = 1.0, dh = -20e3, t_ref = 500.0 }, a = { H2 = 0.5 } },
  { k = { k_ref = 1.0, dh = -60e3, t_ref = 500.0 }, a = { H2O = 1.0 } },
]

[[reaction]]
name = "methanol"
equation = "CO2 + 3 H2 = CH3OH + H2O"    # parsed; element balance checked
basis = "partial_pressure"
units = { rate = "mol/(kg_cat*s)", pressure = "bar" }
k = { k_ref = 1.0, ea = 36e3, t_ref = 500.0 }
orders = { CO2 = 1.0, H2 = 1.0 }
reversible = "thermo"                      # K from ΔG°f
site = "Cu"
site_power = 3

[[reaction]]
name = "rwgs"
equation = "CO2 + H2 = CO + H2O"
basis = "partial_pressure"
units = { rate = "mol/(kg_cat*s)", pressure = "bar" }
k = { k_ref = 1.0, ea = 90e3, t_ref = 500.0 }
orders = { CO2 = 1.0 }
reversible = "thermo"
site = "Cu"
site_power = 1
```

Power-law and Eley–Rideal rates use the same building blocks. A free-form
expression type (a small parsed expression language, never `eval`) is the
escape hatch for anything else, and it gets the same checks.

## Checking that the model is what you meant

`reflux kinetics check reactions.toml` runs the checks below before any
reactor uses the model. `reflux kinetics explain` prints the rendering and
the factor table.

1. **Render it back.** Every rate is printed as a formula with the real
   species names (Unicode in the terminal, LaTeX for notebooks), exactly as
   reflux understood it, so you can hold it against the paper.
2. **Element balance.** Every equation must balance its elements, using the
   formulas in the component data.
3. **Units.** Every adsorption term must be dimensionless, given the stated
   pressure unit. The rate constant's units must turn the driving force into
   the stated rate unit. bar versus Pa errors are caught here.
4. **Thermodynamic consistency.** At an equilibrium composition from Gibbs
   minimisation (roadmap E6), every reversible rate is zero, and its sign
   always points toward equilibrium.
5. **Physical sanity (warnings).** Negative constants; an adsorption
   enthalpy that's positive (adsorption is normally exothermic, so K should
   fall with T); a negative activation energy; a species in an adsorption
   term that isn't in the component list.
6. **Factor table.** At a chosen (T, P, composition), print k(T), the
   driving force, each adsorption term, θ and r for every reaction. This
   shows at once which term dominates.
7. **Derivatives.** ∂r/∂aⱼ and ∂r/∂T are analytic. Stiff PFR integration
   and equation-oriented solving need them, and tests compare them with
   finite differences.

## Fitting parameters

Fitting LHHW constants to lab data is as common as using them, so it is
planned alongside them (roadmap H4). There are two parts:

- k_ref and Eₐ are fitted around a reference temperature, which removes the
  strong correlation between A and Eₐ that makes plain Arrhenius fits
  ill-conditioned.
- Each fit reports confidence intervals and parameter correlations. The
  Python API (M10) lets you drive it from NumPy and pandas.

## How it maps onto V

```v
pub type TempFn = Arrhenius | Constant | VantHoff

pub struct Monomial {
	k     TempFn
	exps  []SpeciesExp // (component index, exponent)
}

pub struct Site {
	name  string
	terms []Monomial  // θ = 1 + Σ terms
}

pub type RateLaw = EleyRideal | Expression | Lhhw | PowerLaw
```

A rate law is a sum type, because the module owns the list of forms. A
reactor takes `[]Reaction` and doesn't care which form each one uses.

## Dependencies

- Rate evaluation, rendering, unit checks and the factor table need only
  component identities. They can be built any time after **M2**.
- Thermodynamic reversibility and the equilibrium check need ideal-gas Gibbs
  energies (M1, M2) and Gibbs minimisation (E6).
- Reactors (CSTR, PFR, packed bed) need the flowsheet (M4) and, for PFRs,
  an ODE integrator (B5).
