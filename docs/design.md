# Thermo design

## Everything is F(T, V, n)

Every thermodynamic model in reflux supplies one function, the reduced
residual Helmholtz energy

    F(T, V, n) = Aʳ(T, V, n) / RT        (V total volume, n mole numbers)

and its derivatives. That is the whole contract. It is the formulation of
Michelsen & Mollerup, *Thermodynamic Models: Fundamentals and Computational
Aspects*, and it's what modern open libraries such as teqp, Clapeyron.jl,
FeOs and thermopack use.

The textbook generic cubic (σ, ε, Ω, Ψ, α(Tr, ω)) describes only one family
of models. These are the ones it can't express, and where each fits here:

| Model | What it changes | Where it goes |
|---|---|---|
| PRSV, Twu, Mathias–Copeman | α has per-component parameters | a variant of `Alpha`, chosen per component (PRSV is done) |
| PRBS (Babalola–Susu) | the attractive term itself | its own `Model`, or a new cubic term |
| CPA | SRK plus a Wertheim association term; V is no longer a cubic root | `F = F_SRK + F_assoc`; the generic volume solver already handles it |
| CPC, PC-SAFT | sums of chain, dispersion and association terms | their own `Model`s |
| Péneloux volume shift | V → V + Σ nᵢcᵢ | a wrapper `Model` around any other one |

The `Model` interface is in [thermo/model.v](../thermo/model.v):

```v
pub interface Model {
	components() []Component
	min_volume(n []f64) f64                            // packing limit, e.g. Σnᵢbᵢ
	volume_derivs(t f64, v f64, n []f64) VolumeDerivs  // F, F_V, F_VV (cheap)
	residual(t f64, v f64, n []f64) Residual           // F, F_V, F_VV, F_T, F_n
}
```

Everything else is derived, once, in model-independent code:

| Property | From F |
|---|---|
| P | RT(n/V − F_V) |
| ∂P/∂V | −RT(F_VV + n/V²) |
| ln φᵢ | F_nᵢ − ln Z |
| Gʳ/RT | F + n(Z − 1) − n ln Z |
| Hʳ | −RT²F_T + PV − nRT |
| Sʳ (at T, P) | −RT F_T − R F + nR ln Z |

## Finding the volume

`solve_tp` doesn't solve a polynomial, because CPA and SAFT don't have one.
It runs two safeguarded Newton iterations:

- **Liquid:** in V, starting just above the packing limit and moving right.
  P(V) is convex there, so the iteration can't overshoot.
- **Vapor:** in ρ = n/V, starting from the ideal gas (capped at half the
  packing limit). P(ρ) is concave there.

Steps are limited so that neither iteration can jump over the van der Waals
loop onto the other branch. Each returns `none` (a V Option) when its branch
doesn't reach the requested pressure. For cubics the tests check both
against the polynomial roots from `Cubic.z_roots`, over a grid of T and P.

## Derivatives

Each model writes its derivatives analytically. The tests make that safe:

- every first derivative is compared against central finite differences;
- Euler's theorem: F is homogeneous of degree 1 in (V, n), so
  F = V·F_V + Σ nᵢ·F_nᵢ;
- Σ nᵢ ln φᵢ = Gʳ/RT.

**Why not automatic differentiation?** It's the modern choice (teqp and
Clapeyron.jl use it), but in V 0.5 operator overloading only works between
operands of the same type, so `2.0 * x` doesn't compile when `x` is a dual
number. Model code written over a generic number type would be cluttered with
conversions. If V lifts that restriction this is worth revisiting.
Complex-step differentiation (V has `math.complex`) is a cheap way to
double-check analytic derivatives in the meantime.

Next to add: F_TT, F_TV, F_nᵢV, F_nᵢT and F_nᵢnⱼ. These give Cp, Newton-based
flash, and the Jacobians a flowsheet solver needs.

## How CPA will fit

```v
pub struct Cpa {
	srk   Cubic     // the physical part
	sites []Sites   // association scheme per component (2B, 4C, …)
	// association energy ε and volume β per component, combining rules
}
```

`residual` returns the SRK residual plus the association term, where the site
fractions X_A come from the site balance at each (T, V, n). Michelsen &
Hendriks (2001) showed that the association term is stationary with respect
to X, so its *first* derivatives need no ∂X/∂(T, V, n). `min_volume` is the
SRK covolume. `model.v`, `saturation.v` and every flash routine stay as
they are.

## Performance

These notes respond to the slow solves and slow sensitivity runs in DWSIM.

- Analytic derivatives are fast, and they give Newton methods their
  Jacobians for free.
- Right now each residual call allocates a few small arrays. That's fine
  until a profile says otherwise. The fix then is to cache the
  temperature-dependent parts (αᵢ, √aᵢ) and pass in scratch buffers.
- Sensitivity analysis: warm-start each case from its neighbour's converged
  solution, run independent cases in parallel with `spawn`, and later get
  exact local sensitivities from the converged Jacobian (implicit function
  theorem) instead of re-solving the flowsheet.
- Benchmark only `-prod` builds made with gcc or clang. tcc doesn't optimise.
