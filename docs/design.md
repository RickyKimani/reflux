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
	residual2(t f64, v f64, n []f64) Residual2         // all of the above + every second derivative
}
```

There are three tiers so that each caller pays only for what it uses. The
volume solver calls the cheapest one on every iteration. The composition
Hessian in `residual2` costs O(nc²) and is needed only by Cp, Newton flash
and Jacobians.

### Properties from F

Everything else is derived once, in model-independent code
([thermo/properties.v](../thermo/properties.v)), following Michelsen &
Mollerup chapter 2. n is Σnᵢ. Residual properties are on the (T, P) basis,
Mʳ = M − Mⁱᵍ at the same T, P and n, and they're extensive.

| Property | From F |
|---|---|
| P | RT(n/V − F_V) |
| (∂P/∂V)_T,n | −RT(F_VV + n/V²) |
| (∂P/∂T)_V,n | P/T − RT·F_TV |
| (∂P/∂nᵢ)_T,V | RT(1/V − F_Vnᵢ) |
| V̄ᵢ | −(∂P/∂nᵢ)/(∂P/∂V) |
| ln φᵢ | F_nᵢ − ln Z |
| (∂ln φᵢ/∂T)_P,n | F_Tnᵢ + 1/T − V̄ᵢ(∂P/∂T)/RT |
| (∂ln φᵢ/∂P)_T,n | V̄ᵢ/RT − 1/P |
| n(∂ln φᵢ/∂nⱼ)_T,P | n·F_nᵢnⱼ + 1 + n(∂P/∂nᵢ)(∂P/∂nⱼ)/(RT·∂P/∂V) |
| Gʳ/RT | F + n(Z − 1) − n ln Z  (= Σ nᵢ ln φᵢ) |
| Hʳ | −RT²F_T + PV − nRT |
| Sʳ | −RT F_T − R F + nR ln Z |
| Uʳ | −RT²F_T |
| Aʳ | Uʳ − TSʳ = RT(F − n ln Z) |
| Cvʳ | −RT(T·F_TT + 2F_T) |
| Cpʳ | Cvʳ − T(∂P/∂T)²/(∂P/∂V) − nR |

The n·∂ln φᵢ/∂nⱼ matrix is symmetric, and Gibbs–Duhem makes each of its
columns sum to zero when weighted by nᵢ. Both are tested.

### Ideal gas, reference state and total properties

`Model` is purely residual. Total properties come from a `Fluid`
([thermo/fluid.v](../thermo/fluid.v)), which pairs a `Model` with each
component's `IdealGas` data ([thermo/ideal.v](../thermo/ideal.v)).

**Reference state:** elements in their standard states at T° = 298.15 K
and P° = 1 bar. For component i:

    Hᵢⁱᵍ(T) = ΔH°f,i + ∫_T°^T Cp,i dT
    Sᵢⁱᵍ(T, P) = S°ᵢ + ∫_T°^T Cp,i/T dT − R ln(P/P°)

With this convention a reactor's energy balance needs no separate heat of
reaction. ΔH°f and S° are Options; if either is missing, `Props.absolute`
is false. Differences between states of the same composition are still
correct then.

**Cp forms** (`CpForm` sum type), each with closed-form integrals:

| Form | Cp | ∫Cp dT | ∫Cp/T dT |
|---|---|---|---|
| `Polynomial` | Σ cₖTᵖᵏ | Σ cₖTᵖᵏ⁺¹/(pₖ+1); ln T when pₖ = −1 | Σ cₖTᵖᵏ/pₖ; ln T when pₖ = 0 |
| `AlyLee` (DIPPR 107) | A + B(x/sinh x)² + D(y/cosh y)², x = C/T, y = E/T | AT + BC·coth x − DE·tanh y | A ln T + B(x coth x − ln sinh x) − D(y tanh y − ln cosh y) |

sinh and cosh overflow above about 710, which C/T reaches at low T. The
Aly–Lee terms are therefore written in terms of e^(−x), and ln cosh x uses
log1p(2 sinh²(x/2)) for small x, where cosh x is too close to 1.

**Totals** for a state (T, P, V, n), with xᵢ = nᵢ/n and m = Σ nᵢMᵢ:

| | |
|---|---|
| H | Σ nᵢHᵢⁱᵍ + Hʳ |
| S | Σ nᵢ(Sᵢⁱᵍ − R ln xᵢ) + Sʳ |
| G, U, A | H − TS, H − PV, U − TS |
| Cp | Σ nᵢCp,i + Cpʳ |
| Cv | Σ nᵢCp,i − nR + Cvʳ |
| speed of sound | w² = −(V²/m)·(Cp/Cv)·(∂P/∂V)_T |
| Joule–Thomson | μ_JT = [T(∂V/∂T)_P − V]/Cp, with (∂V/∂T)_P = −(∂P/∂T)_V/(∂P/∂V)_T |

μ_JT does *not* go to zero at zero pressure. It tends to (T·dB/dT − B)/Cp,
where B is the second virial coefficient, which is why the low-pressure test
checks only Cp − Cv and w.

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

- every first and second derivative is compared against central finite
  differences;
- Euler's theorem: F is homogeneous of degree 1 in (V, n), so
  F = V·F_V + Σ nᵢ·F_nᵢ. Its derivatives in V and nᵢ are homogeneous of
  degree 0, and F_T of degree 1, which gives three more identities:
  V·F_VV + Σ nᵢF_Vnᵢ = 0, V·F_Vnᵢ + Σⱼ nⱼF_nᵢnⱼ = 0, and
  V·F_TV + Σ nᵢF_Tnᵢ = F_T;
- Σ nᵢ ln φᵢ = Gʳ/RT.

**Why not automatic differentiation?** It's the modern choice (teqp and
Clapeyron.jl use it), but in V 0.5 operator overloading only works between
operands of the same type, so `2.0 * x` doesn't compile when `x` is a dual
number. Model code written over a generic number type would be cluttered with
conversions. If V lifts that restriction this is worth revisiting.
Complex-step differentiation (V has `math.complex`) is a cheap way to
double-check analytic derivatives in the meantime.

## Second derivatives of the cubic

These follow Michelsen & Mollerup's structure. The code in
[thermo/cubic.v](../thermo/cubic.v) uses the same names.

The cubic's F is built from two functions of (V, B) and one of (T, n):

    F = −n·g(V, B) − D(T, n)/(RT) · f(V, B)
    g = ln(1 − B/V)
    f = ln(d₁/d₂) / (B(δ₁ − δ₂)),  d₁ = V + δ₁B,  d₂ = V + δ₂B
    D = A = ΣᵢΣⱼ nᵢnⱼ(1 − kᵢⱼ)·sᵢsⱼ,  sᵢ = √aᵢ(T),  B = Σ nᵢbᵢ

When δ₁ = δ₂ = δ, f takes its limit 1/(V + δB). Van der Waals is δ = 0.

**g and f.** Every needed derivative has a closed form:

| | g | f |
|---|---|---|
| ∂/∂V | 1/(V−B) − 1/V | −1/(d₁d₂) |
| ∂²/∂V² | 1/V² − 1/(V−B)² | (d₁ + d₂)/(d₁²d₂²) |
| ∂/∂B | −1/(V−B) | −(f + V·f_V)/B |
| ∂²/∂B∂V | 1/(V−B)² | (δ₁d₂ + δ₂d₁)/(d₁²d₂²) |
| ∂²/∂B² | −1/(V−B)² | −(2f_B + V·f_BV)/B |

The two f_B entries come from homogeneity: f is of degree −1 in (V, B), and
f_B of degree −2. The general formulas are also correct in the limit
δ₁ = δ₂.

**D.** With qᵢ = Σⱼ(1 − kᵢⱼ)sⱼnⱼ and rᵢ = Σⱼ(1 − kᵢⱼ)s′ⱼnⱼ, where a prime
is d/dT:

| | |
|---|---|
| D | Σᵢ nᵢsᵢqᵢ |
| Dᵢ = ∂D/∂nᵢ | 2sᵢqᵢ |
| Dᵢⱼ | 2(1 − kᵢⱼ)sᵢsⱼ |
| D_T | 2Σᵢ nᵢs′ᵢqᵢ |
| D_iT | 2(s′ᵢqᵢ + sᵢrᵢ) |
| D_TT | 2Σᵢ nᵢ(s″ᵢqᵢ + s′ᵢrᵢ) |

The symmetric forms need kᵢⱼ = kⱼᵢ, which `new_cubic` enforces. Since
sᵢ² = a_c,i·αᵢ(Tr), the s derivatives come from the α function:

    s′ = a_c·α_T / (2s),   s″ = (a_c·α_TT/2 − s′²) / s,
    α_T = α′/Tc,   α_TT = α″/Tc²

where α′ and α″ are taken with respect to Tr.

**F.** Writing RT for R·T:

| | |
|---|---|
| F_V | −n·g_V − D/RT·f_V |
| F_VV | −n·g_VV − D/RT·f_VV |
| F_T | −f/RT·(D_T − D/T) |
| F_TV | −f_V/RT·(D_T − D/T) |
| F_TT | −f/RT·(D_TT − 2D_T/T + 2D/T²) |
| F_nᵢ | −g − n·g_B·bᵢ − f·Dᵢ/RT − D/RT·f_B·bᵢ |
| F_Vnᵢ | −g_V − n·g_BV·bᵢ − f_V·Dᵢ/RT − D/RT·f_BV·bᵢ |
| F_Tnᵢ | −f/RT·(D_iT − Dᵢ/T) − f_B·bᵢ/RT·(D_T − D/T) |
| F_nᵢnⱼ | −g_B(bᵢ + bⱼ) − n·g_BB·bᵢbⱼ − [Dᵢⱼ·f + f_B(Dᵢbⱼ + Dⱼbᵢ) + D·f_BB·bᵢbⱼ]/RT |

F_TT comes from d²(D/T)/dT² = D_TT/T − 2D_T/T² + 2D/T³. F_nᵢnⱼ is
visibly symmetric in i and j.

**α.** Each α function returns α, α′ and α″ in Tr. For the Soave form
α = k², with k = 1 + m(1 − √Tr):

    α′ = 2k·k′,   α″ = 2k′² + 2k·k″,   k′ = −m/(2√Tr),   k″ = m/(4·Tr^1.5)

PRSV is the same, with m replaced by κ(Tr) = κ₀ + κ₁h, h = (1 + √Tr)(0.7 − Tr):

    k′ = κ′u + κu′,   k″ = κ″u + 2κ′u′ + κu″,   u = 1 − √Tr
    h′ = (0.7 − Tr)/(2√Tr) − (1 + √Tr),   h″ = −1/√Tr − (0.7 − Tr)/(4·Tr^1.5)

**Twu (1991):** α = Tr^(N(M−1))·exp(L(1 − Tr^(NM))). This is easiest
through ln α:

    (ln α)′ = N(M−1)/Tr − L·NM·Tr^(NM−1)
    (ln α)″ = −N(M−1)/Tr² − L·NM(NM−1)·Tr^(NM−2)
    α′ = α·(ln α)′,   α″ = α·[(ln α)″ + (ln α)′²]

**Mathias–Copeman (1983):** α = k², with k = 1 + c₁u + c₂u² + c₃u³ and
u = 1 − √Tr, so k′ = k_u·u′ and k″ = k_uu·u′² + k_u·u″. Above Tc the
default (`truncated`) keeps only c₁, as in the original paper and in
`thermo`. At Tr = 1, u = 0, so α and α′ are continuous, but α″ drops by
exactly c₂, and Cp steps with it. `truncated: false` gives the smooth form.

**A caveat.** At very high Tr the Soave form passes through k = 0, which
makes α = 0, and s′ divides by s there. It's unphysical territory anyway
(roughly Tr > 5 for typical m), and the Twu α is the physical fix.

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
