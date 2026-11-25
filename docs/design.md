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
	min_volume(n []f64) f64                                 // packing limit, e.g. Σnᵢbᵢ
	volume_function(t f64, n []f64) VolumeFunction          // F at fixed T and n, as a function of V
	residual(t f64, v f64, n []f64) Residual                // F, F_V, F_VV, F_T, F_n
	residual2(t f64, v f64, n []f64) Residual2              // all of the above + every second derivative
}

pub interface VolumeFunction {
	at(v f64) VolumeDerivs  // F, F_V, F_VV: what the volume solver needs
	residual(v f64) Residual // the state at the root it found
}
```

There are three tiers so that each caller pays only for what it uses. The
composition Hessian in `residual2` costs O(nc²) and is needed only by Cp,
Newton flash and Jacobians.

**The volume solver's tier fixes T and n.** A volume root is found at fixed
T and n by iterating on V alone, ten or twenty steps. So `volume_function`
computes the parts that depend only on T and n once, and each step costs
only the terms in V: a cubic's A(T, n) = ΣΣ nᵢnⱼ(1 − kᵢⱼ)√(aᵢaⱼ) and
B(n) = Σ nᵢbᵢ, CPA's association strengths, SAFT's segment diameters. Its
`residual` then gives the state at the root without computing them again.
Until M3 step 10 this tier was `volume_derivs(t, v, n)`, which recomputed
every αᵢ(T) and the double sum on every step; the step-10 benchmarks found
the volume roots taking 88% of a PT flash. (That amends M1's D2.)
`thermo.volume_derivs(m, t, v, n)` remains, for one volume.

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
([thermo/fluid.v](../thermo/fluid.v)), which pairs a phase model (below)
with each component's `IdealGas` data ([thermo/ideal.v](../thermo/ideal.v)).

**Reference state:** elements in their standard states at T° = 298.15 K
and P° = 1 bar. For component i:

    Hᵢⁱᵍ(T) = ΔH°f,i + ∫_T°^T Cp,i dT
    Sᵢⁱᵍ(T, P) = S°ᵢ + ∫_T°^T Cp,i/T dT − R ln(P/P°)

**R** here and everywhere is `quantity.molar_gas_constant`, N_A·k_B =
8.31446261815324 J/(mol·K), the 2019 SI's exact value and the double
nearest it. reflux has one R: `thermo.r_gas`, `activity.r_gas`,
`kinetics.gas_constant` and an expression's `R` are that constant, and a
test holds them so (`fit/gas_constant_test.v`).

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

## The phase seam

Everything above is how an equation of state gives a phase. A γ–φ liquid
gives one with no F at all, so the flash asks for less than a `Model`
(M6's D3, [thermo/phases.v](../thermo/phases.v)): a **`PhaseModel`**, a
phase at (T, P, n).

| A phase model gives | A Helmholtz model's (`VolumeRoots`) |
|---|---|
| `components()` | the model's |
| `phase(T, P, n, root)`: V, Z and ln φᵢ of the liquid-like, the vapour-like or the phase of lower G | `solve_tp`'s volume root |
| `ln_phi_derivs(state)`: ∂ln φᵢ/∂T, ∂ln φᵢ/∂P, n·∂ln φᵢ/∂nⱼ, or why there are none | from F's second derivatives |
| `phase_label(state)` | Π, with the four-covolume rule |
| `phase_residual(state)`: Hʳ, Sʳ, Cpʳ, Cvʳ, (∂V/∂T)_P, (∂P/∂V)_T, each NaN where there's none | one `residual2` |

A phase's mass density is Σ nᵢMᵢ/V, of its state and the components' molar
masses. `thermo.by_roots(model)` puts a `Model` behind the seam, and each of
its methods gives the function's bits. `new_fluid` takes a `Model` and does
so itself; `new_phase_fluid` takes any phase model.

**What works in F(T, V, n) itself** stays on `Model`: the TV and UV
flashes, the reactors, `equilibrium`, Psat, and a fit of a cubic's kᵢⱼ.
Each takes it from `Fluid.helmholtz()`, which gives the model behind the
roots, or a `NoHelmholtzError` saying why a fluid has none.

Since a residual can be NaN, a PH or PS flash takes an energy balance for
met only when |r| ≤ etol holds, which a NaN fails, and never when |r| > etol
doesn't, which a NaN passes. Without that, Brent closes in on the edge of
where the NaN begins, as on a change of sign, and the flash returned it.

**A workspace for one flash** (M6 step 5). A model is immutable and shared,
yet a γ–φ liquid's Pˢᵃᵗ, Vᴸ, φˢᵃᵗ and τ depend on T alone, and a flash asks
for many phases at one T. A model that keeps them is also a
`thermo.Scoped`, whose `scoped()` is a copy with a workspace of its own;
`thermo.for_one_thread(m)` is that copy, or m, and every flash calls it
first. The seam itself is unchanged, as M8's `Warm` is optional to a unit.
A copy gives the bits of a model that makes its own each call.

**`numeric.newton` passes on why a system can't be evaluated** at its
start: "the system cannot be evaluated at x0" is followed by the system's
own error, or by the residual or Jacobian entry that wasn't finite. A
phase model's refusal there reaches the caller with its reason (M6 step 6).

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

**Near a critical point** the iterations stop on the residual instead. At
(Tc, Pc), P(V) − p has a triple root: dP/dV and d²P/dV² vanish there too.
Newton then converges only linearly, and rounding limits V to about ε^(1/3),
roughly 1e-5 relative. Close to the critical point dP/dV is merely tiny,
with the same effect. So once |P − p| is down to a few rounding errors of
P's own terms, 16ε·RT(n/V + |F_V|), and the Newton step still isn't below
1e-12, the volume is as converged as it can get and is returned. Elsewhere
the step test is met first, and nothing changes. M3 found this: before it,
the solver failed on 276 of 528 requests within 1e-6 of a pure component's
Tc and 1e-8 of its Pc. At the critical point itself, the liquid iteration
now stops just below Vc and the vapour one just above.

**A cubic's workspace** (M8 step 7). `solve_tp` allocates 944 bytes a call
on a five-component cubic: 704 for the volume function, its √aᵢ arrays and
the struct boxed as an interface, and 240 for the state's n and ln φ. A
reactor reads its phase about twenty times an integration step, so a
study of the owner's PFR allocated 11.8 MB a point, and threads waited on
the collector's lock. `Cubic.work()` keeps √aᵢ and ∂√aᵢ/∂T at the last T,
Aᵢ, ln φ, and one volume function, boxed once, which each call moves to
its composition in place. It does `solve_tp`'s arithmetic in `solve_tp`'s
order, `attraction`, `fix`, the volume solvers and `state_at`'s ∂F/∂nᵢ − ln Z
term by term, so it gives `solve_tp`'s bits, and allocates nothing. A
workspace is mutable: one per thread, one per reactor's calculation.

## Phase identification

`thermo.pip` ([thermo/phase.v](../thermo/phase.v)) is the phase
identification parameter of Venkatarathnam and Oellrich (2011):

    Π = V·[(∂²P/∂T∂V)/(∂P/∂T) − (∂²P/∂V²)/(∂P/∂V)]

`phase_label` calls a single phase liquid when Π > 1 and vapour otherwise,
as Caleb Bell's `thermo` does, but for a phase more than four covolumes
large, less dense than its critical density, which is vapour whatever Π
says (amended in M9 step 9). It needs no saturation curve or critical
point, so it works for any mixture, supercritical ones included.

**Why the covolumes.** At low density

    Π ≈ 1 + (B − T·dB/dT)/V,

B the second virial coefficient, and B − T·dB/dT is above zero for a gas
above its Joule–Thomson inversion temperature at zero pressure: hydrogen's
is about 200 K, nitrogen's 620 K, methane's 950 K. So Π alone calls such a
gas liquid however dilute, and M9's reformer gas at 1000 K came out of the
flash a liquid, which a drum would have sent out of its liquid port. A
cubic's critical volume is 3 (van der Waals) to 3.95 (Peng–Robinson)
covolumes, and Kay's rule on the components' is the mixture's covolume times
the same ratio, so four covolumes is less dense than the pseudo-critical
density on any cubic: no liquid is that dilute. Of M3's 1,327 fixture roots,
17 supercritical ones, 1.3 to 2.5 times their Tc with Π from 1.004 to 1.17,
are vapour where `thermo` calls them liquid; every other label is
`thermo`'s.

∂P/∂V and ∂P/∂T are analytic. The second derivatives are third derivatives
of F, which the `Model` contract doesn't have, so they're central
differences of the analytic first ones (D6). ∂²P/∂V² is differenced in V,
with a step of 1e-5·(V − V_min) so that a liquid squeezed against its
covolume still gets a small enough step. ∂²P/∂T∂V is differenced in T, with
a step of 1e-5·T, as the T-derivative of ∂P/∂V. Π then agrees with
`thermo`'s analytic value to about 5e-9.

At a pure component's exact critical point ∂P/∂V and ∂²P/∂V² are both zero,
so Π is 0/0 and no label means anything. In a two-phase split the flash
labels by mass density instead, the less dense phase being the vapour,
because Π can call both phases liquid-like. (`thermo` then calls both
liquids.)

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
sᵢ = √a_c,i·√αᵢ(Tr), the s derivatives come from the square root of α,
with no division:

    s′ = √a_c·(√α)′/Tc,   s″ = √a_c·(√α)″/Tc²

where the primes on √α are taken with respect to Tr (`Alpha.root`; see
"Through the Soave zero" below).

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

**Through the Soave zero.** At Tr = (1 + 1/m)² the Soave form reaches
k = 0, so α = 0. With PR that's about 1,400 K for nitrogen and 2,400 K for
methane, so a furnace or a reformer gets there. Building s as √(a_c·α) and
dividing by it for s′ and s″ would give 0/0 at that point. So each α
function also returns its square root, `Alpha.root`. For the forms written
as α = k² (Soave, PRSV, Mathias–Copeman) that's k itself, signed, and

    s = √a_c·k,   s′ = √a_c·k′/Tc,   s″ = √a_c·k″/Tc²,

with no division at all. The other forms take √α directly: Tr^(−1/4) for
1/√Tr, and exp(½ ln α) for Twu, with (√α)′ = √α·½(ln α)′.

Past the zero k is negative and α = k² rises again, as the published form
does, and a pure component is unchanged. A mixture's cross term
(1 − kᵢⱼ)sᵢsⱼ turns negative while only one of the pair has passed its
zero. Taking |k| instead, as `thermo` does, would put a kink in a(T) and a
step in Hʳ at the zero; the signed root keeps everything smooth. Neither is
physical: Soave's α isn't meant for such Tr. The Twu α, which stays positive
and keeps falling, is the one to use for components that far above their
Tc.

## Volume shift

`Translated` ([thermo/translated.v](../thermo/translated.v)) wraps any
`Model`. With C = Σ nᵢcᵢ and W = V + C, it is the inner model evaluated at
W:

    F′(T, V, n) = F(T, W, n) + n·ln(V/W),   so   P′(T, V, n) = P(T, W, n)

In the derivatives below, F and its subscripts are the inner model's at
(T, W, n). Because ∂W/∂nᵢ = cᵢ, each n-derivative picks up a V-derivative:

| | |
|---|---|
| F′_V | F_V + n/V − n/W |
| F′_VV | F_VV − n/V² + n/W² |
| F′_T, F′_TT, F′_TV | F_T, F_TT, F_TV |
| F′_nᵢ | F_nᵢ + cᵢF_V + ln(V/W) − n·cᵢ/W |
| F′_Vnᵢ | F_Vnᵢ + cᵢF_VV + 1/V − 1/W + n·cᵢ/W² |
| F′_Tnᵢ | F_Tnᵢ + cᵢF_TV |
| F′_nᵢnⱼ | F_nᵢnⱼ + cⱼF_Vnᵢ + cᵢF_Vnⱼ + cᵢcⱼF_VV − (cᵢ + cⱼ)/W + n·cᵢcⱼ/W² |

What follows, and is tested:

- V′ = V − C at the same (T, P), and ln φ′ᵢ = ln φᵢ − cᵢP/RT;
- so Psat and phase splits don't change;
- G′ʳ = Gʳ − PC, H′ʳ = Hʳ − PC, and Sʳ and Cpʳ are unchanged (for constant
  c);
- the packing limit moves to B − C, which is why `translate` rejects any
  cᵢ ≥ bᵢ.

## Rachford–Rice

`numeric.rachford_rice` ([numeric/rachford_rice.v](../numeric/rachford_rice.v))
solves, for fixed K-values,

    f(β) = Σᵢ zᵢ(Kᵢ − 1) / (1 + β(Kᵢ − 1)) = 0,   xᵢ = zᵢ / (1 + β(Kᵢ − 1)),   yᵢ = Kᵢxᵢ.

The pole-free window runs from c₁ = 1/(1 − K_max) < 0 to cₙ = 1/(1 − K_min) > 1,
taking K_max and K_min over the components present. f falls from +∞ to −∞
across it, so the root is unique. It may lie outside [0, 1]: that's the
negative flash, and the compositions are still positive.

The numerical danger is the denominator eᵢ = 1 + β(Kᵢ − 1). It is tiny for
a heavy component when the liquid is a trace, and computing it from β
subtracts two numbers near 1. So β is written as a position in the window,
β = c₁ + t(cₙ − c₁) with s = 1 − t. Substituting gives

    eᵢ = s·uᵢ + t·vᵢ,   uᵢ = (K_max − Kᵢ)/(K_max − 1),   vᵢ = (Kᵢ − K_min)/(1 − K_min),

and both uᵢ and vᵢ are non-negative for every component present. (Check
the ends: at t = 0, β = c₁ and the K_max component's e is 0; at t = 1,
β = cₙ and the K_min component's is.) A sum of non-negative terms has no
cancellation, so every eᵢ has full relative precision anywhere in the
window, given t and s to full relative precision. The smaller of the two is
stored and the larger is 1 minus it. So a root 1e-30 from an end still has
all its digits, and so do β and L near their own small values: β from
c₁ + (cₙ − c₁)t, L from (1 − cₙ) + (cₙ − c₁)s.

The equation solved is Leibovici and Neoschil's, f multiplied by the
distances to both poles:

    Φ(t) = t·s·Σᵢ zᵢ(Kᵢ − 1)/(s·uᵢ + t·vᵢ).

It's finite at both ends, Φ(0) = Σ over the K_max components of zᵢ(Kᵢ − 1)/vᵢ
and Φ(1) = Σ over the K_min ones of zᵢ(Kᵢ − 1)/uᵢ, and close to linear near
each. For two components it is exactly linear, so the first guess, the
straight line between Φ(0) and Φ(1), is the root. Newton runs on Φ in t,
and bisection takes over whenever a step would leave the bracket.

What β can be trusted to is set by the K-values, not the method. Given eᵢ
to a few ulps, β = (eᵢ − 1)/(Kᵢ − 1) carries an absolute error of about
ε·(|β| + min(|c₁|, cₙ)). That's about ε for an ordinary mixture, but about
2e-8 when every K is within 1e-8 of 1. Mass balance doesn't suffer, because
L·xᵢ + β·yᵢ = xᵢ·eᵢ = zᵢ holds whatever β's error.

## Stability

`flash.stability` ([flash/stability.v](../flash/stability.v)) is Michelsen's
(1982) tangent-plane test. For a trial phase of mole numbers W, with
w = W/ΣW and dᵢ = ln zᵢ + ln φᵢ(z),

    tm*(W) = 1 + Σᵢ Wᵢ(ln Wᵢ + ln φᵢ(w) − dᵢ − 1).

Its stationary points satisfy ln Wᵢ + ln φᵢ(w) = dᵢ. There tm* = 1 − ΣW,
and the normalised distance tm = Σ wᵢ(ln wᵢ + ln φᵢ(w) − dᵢ) equals
−ln ΣW. The feed is unstable when some stationary point has tm < 0, and
then that w is the first estimate of the new phase. Every phase takes its
volume root of lower Gibbs energy.

**Trials.** One starts vapour-like (Wᵢ = zᵢKᵢ) and one liquid-like
(Wᵢ = zᵢ/Kᵢ), with Wilson's K. Then comes one nearly pure phase of each
component present: Wᵢ = 1, with the others at 1e-10·zⱼ. A caller can add
more: the phases of a split, to re-check it.

The pure trials were added in M3 step 4. Wilson's two trials can miss a
liquid–liquid split entirely. With water, n-decane and methane (PR, a
water–hydrocarbon kᵢⱼ of 0.5) at 350 K and 5 bar, they call the methane-rich
phase stable, yet a nearly pure water phase lies at tm = −2.2 against it.
On 6,000 random PR mixtures the pure trials find 29 instabilities that
Wilson's missed, and once they're there, dropping either Wilson trial
changes no verdict. They cost about 2.4 times the work on the reference
set. `pure_trials: false` turns them off where no second liquid can form.

**Successive substitution,** ln Wᵢ ← dᵢ − ln φᵢ(w), lowers tm* at every
step. Every fifth step Michelsen's dominant-eigenvalue method extrapolates:
with Δ the last step and λ = Δ·Δ/(Δ_prev·Δ), substitution behaves like a
geometric series of ratio λ, and adding λ/(1 − λ) times the last step sums
the rest. Michelsen and Mollerup keep the extrapolated point only if it
lowers tm*. On 6,000 random states that test changed no verdict, so reflux
doesn't make it. Without the extrapolation the reference set takes a third
more steps and nearly three times as many second-order finishes.

**The second-order finish.** After 20 steps the trial switches to Newton
in αᵢ = 2√Wᵢ, minimising tm* itself. The gradient and Hessian are

    gᵢ = √Wᵢ·rᵢ,   Hᵢⱼ = δᵢⱼ(1 + rᵢ/2) + √(WᵢWⱼ)·(n∂ln φᵢ/∂nⱼ)/ΣW,

with rᵢ = ln Wᵢ + ln φᵢ − dᵢ. The δᵢⱼ term keeps H well conditioned even
with traces. Two choices mattered, both found near the methane/n-butane
critical locus:
- **Minimise tm*, don't solve g = 0.** On the way to the trivial solution
  near a critical point, tm* has a long valley whose curvature is nearly
  zero, and ‖g‖ has dips there that aren't stationary points. A root-finder
  (numeric.newton) stopped in them and handed the trial back to
  substitution, which then took up to 465 steps. One state failed outright.
- **Know f's rounding.** tm* = 1 + Σ… is a small difference of terms of
  size 1, so its rounding is about ε, not ε·|tm*|. numeric.minimize takes
  that as `fscale`. Where f stops showing progress it finishes on the
  gradient (the polish), and it counts that as convergence only if the
  unshifted Newton model agrees.

On 33,600 states 0.1 bar apart near that critical locus, no trial takes
more than 42 steps and none fails.

**The trivial solution.** A trial that ends within 1e-5 of the feed in
every ln(Wᵢ/zᵢ) is the feed itself: tm = 0, and it says nothing. It's
flagged trivial and never offered as a new phase.

**The verdict.** The feed is unstable when a non-trivial trial reaches
tm < −1e-10. A trial that runs out of steps still proves instability that
way, since any W with tm < 0 is a proof. If no trial does and one didn't
converge, the result is a FlashError, not "stable".

## PT flash

`flash.pt` ([flash/pt.v](../flash/pt.v)) tests the feed's stability first.
A stable feed is one phase, labelled by its PIP, with β = 1 for a vapour and
0 for a liquid. An unstable feed is split into two phases.

**Starting K.** From an unstable trial phase W, stationarity,
Wᵢφᵢ(w) = zᵢφᵢ(z), gives Kᵢ = Wᵢ/zᵢ. The trials are tried in order of their
tm, lowest first. Whether W is the lighter phase doesn't matter: K and 1/K
describe the same split with the roles swapped. Rachford–Rice then gives
β′ = 1 − β, x′ = y and y′ = x, every later step mirrors the other, and the
phases are labelled by density at the end.

**Successive substitution,** ln Kᵢ ← ln φᵢᴸ(x) − ln φᵢⱽ(y), with x, y and β
from Rachford–Rice at each step. Every fifth step it uses the same
dominant-eigenvalue extrapolation as stability. Rachford–Rice may give β
outside [0, 1] on the way (the negative flash); only the converged split
must have 0 < β < 1.

Unlike stability's, the split's extrapolation is kept only if it lowers the
residual, max |ln fᵢⱽ − ln fᵢᴸ|, below that of the split it was made from.
Otherwise the plain step it replaced is taken, and the rejected one counts
as a step. Michelsen and Mollerup's test is that G falls, but G has no
meaning at a β outside [0, 1], which the negative flash passes through, and
the residual does. M10a step 8 found the case. A natural gas (PR: 80%
methane, ethane to n-pentane) at 265 K and 81 bar is well inside its
two-phase region, with K from 1.7 to 0.07. There the extrapolation at the
sixth step threw β from 0.14 to 0.87 and the residual from 0.005 to 1.5.
The split then slid off towards the feed and collapsed onto it, and the
flash failed. On 16,160 PT flashes of that gas, from 150 to 350 K and 1 to
160 bar, the test undoes 392 extrapolations. It removes the one failure and
leaves the total steps the same: 67,754 against 67,751, with that point's
21 included.

**Newton on G.** After 20 steps it switches to minimising the Gibbs energy
of the split over the vapour mole numbers v (one mole of feed, l = z − v):

    G/RT = Σᵢ vᵢ ln fᵢ(y) + Σᵢ lᵢ ln fᵢ(x),   gᵢ = ln fᵢ(y) − ln fᵢ(x),
    Hᵢⱼ = (δᵢⱼ/yᵢ − 1 + Φⱽᵢⱼ)/V + (δᵢⱼ/xᵢ − 1 + Φᴸᵢⱼ)/L,

where Φ is each phase's n·∂ln φᵢ/∂nⱼ. This uses numeric.minimize, with
fscale set to the size of the terms G is summed from. Near the critical
locus, 44 of the reference set's 166 splits end this way.

A component with almost no moles in a phase, 1e-305 say, puts 1/(V·yᵢ) ≈
1e305 on H's diagonal. Its step leaves the domain unless the shift μ is as
large as that, and then no step moves anything else. So μ climbs to
overflow, and minimize fails, as it must. Until step 9 it looped forever
there, because its limit on μ, 1e20 times H's largest diagonal, had
overflowed to +inf first. A failed Newton leaves substitution to carry on.

**Accepted** only with every |ln fᵢⱽ − ln fᵢᴸ| ≤ tol (1e-12) and 0 < β < 1,
and a Newton result is held to the same test. A split that slides onto the
feed (every |ln K| < 1e-4) is refused: there the fugacities are equal too,
and it would pass for two identical phases.

**Extreme K-values.** A K is held within e^±700, so that it and the mole
fractions it makes stay finite doubles. A component whose K would lie
beyond that has a mole fraction below 1e-300 in one phase, which is zero
to double precision: tripalmitin's in the vapour below about 150 K, or
mercury's in PR with ChemSep's constants. Below that floor (`trace_floor`)
its fugacities aren't compared, since they can't be represented well
enough to agree. The robustness sweep found 49 such points that the
unbounded K turned into 0 or +inf for Rachford–Rice.

**More than one start.** The starts are tried in turn until one gives a
verified split. A split with an unstable phase can be a local minimum of G
that another start avoids. If none does better, the answer is
`more_phases`, with that split's proof, rather than a later start's
failure to converge.

**The feed's test stops early** (M4 step 9). It stops at the first proof of
instability once Wilson's two trials have run (`first_proof`): one proof is
all a split needs, and the nearly pure trials, nc more, were most of a split
flash's stability work. Both of Wilson's run, because the better of the two is
the better start: stopping at the first proof alone took 7% more split steps
on the reference set. If no split from those trials stands, the test runs all
of them and the starts are tried as above. All 166 of the reference set's splits
start from Wilson's trials alone; 11 of the robustness sweep's 10,000 points
need the rest.

**Verified** (`verify`, on by default): each phase is re-tested for
stability, with the other phase as a trial. An unstable phase means a third
phase or a better split. That's a FlashError (`more_phases`) rather than an
answer, since M3 handles two phases. The water, n-decane and methane
mixture above is the example. There, a two-phase flash, reflux's without
the re-test or `thermo`'s FlashVL, returns a split whose fugacities are
equal but which isn't the equilibrium.

**The order of the phases.** By increasing *mass* density, ΣnᵢMᵢ/V (D6,
amended in step 6). Molar volume isn't a substitute. In methane with
n-decane at 300 K and 200 bar, the methane-rich phase has the smaller molar
volume and yet weighs 168 kg/m³, against the decane-rich phase's 585. Step
4 compared molar volumes and called the decane-rich phase the vapour.

**Their names** (M6's D9, `labelled` in [flash/multi.v](../flash/multi.v)).
Each phase is named by its model's `phase_label` (an equation of state's
by Π and the four-covolume rule, above; a γ–φ fluid's by the model that
made it); of the phases called vapour only the least dense stays so, the
rest being liquids; and of the liquids the less dense is the liquid and
the denser liquid2. A third liquid has no name and is refused. β is the
vapour's fraction, or 0 with none. Until M6 the less dense phase of a split
was the vapour whatever Π said, so a liquid–liquid split by an equation of
state was called vapour and liquid: of M3's sweep's 3,344 splits by PR,
2,289 are liquids by Π, now liquid and liquid2.

## Three phases

M6's D8 ([flash/multi.v](../flash/multi.v)), Michelsen's stage-wise flash.
A two-phase split whose verification finds a phase unstable takes that
phase's trial composition as a third phase: the split's two and the trial,
the trial with 1% of the feed and the others 99% of theirs. With no moles
the trial would hold a component that lives almost only in it in no phase
at all.

**The split.** Successive substitution: at each step each phase's φᵢⱼ at its
composition, then the fractions and compositions from Michelsen's (1994)
multiphase Rachford–Rice (`numeric.multiphase_rr`,
[numeric/multiphase.v](../numeric/multiphase.v)) on Kᵢⱼ = 1/φᵢⱼ:

    Q(β) = Σⱼ βⱼ − Σᵢ zᵢ ln Eᵢ,   Eᵢ = Σⱼ βⱼKᵢⱼ,   β ≥ 0,
    xᵢⱼ = zᵢKᵢⱼ/Eᵢ,   ∂Q/∂βⱼ = 1 − Σᵢ xᵢⱼ,   ∂²Q/∂βⱼ∂βₖ = Σᵢ zᵢKᵢⱼKᵢₖ/Eᵢ².

Q is convex, so its minimum over β ≥ 0 is unique; there each phase present
has Σᵢ xᵢⱼ = 1 (so Σⱼ βⱼ = 1 and Σⱼ βⱼxᵢⱼ = zᵢ), and a phase held at 0 has
Σᵢ xᵢⱼ ≤ 1, stable absent. Each component's K-values are scaled by its
largest, so none overflows. Newton runs on the free phases, a phase at 0
whose step would go below it held there and the step taken again; the step
is cut short to keep β ≥ 0, and halved until Q falls or rises by no more
than its rounding (near the minimum a Newton step moves Q by about the
gradient squared). The phase that cuts the step short goes to 0 exactly:
β + (−β/d)·d rounds to about β·ε, not 0, and a phase leaving by such
residues shrinks by ε a step until −β/d underflows and nothing moves, as two
components in three phases, with a singular Hessian, did after 20 steps
(0.6.0). With two phases it is a Rachford–Rice with the negative flash cut
off at 0 and 1: the same root inside, to 2e-15.

A phase held at 0 three steps running leaves the split, and of two phases
that come within 1e-4 in every ln x, one does. After 30 steps the split
goes on by Newton on G.

**Newton on G over the phases.** Each component present has a dependent
phase rᵢ, the one holding the most of it, which takes the rest of the feed;
its moles in every other phase are the variables. Then

    gᵢⱼ = ln fᵢⱼ − ln fᵢᵣᵢ,
    Hᵢⱼ,ₖₘ = (δⱼₘ − δⱼᵣₖ)Aʲᵢₖ − (δᵣᵢₘ − δᵣᵢᵣₖ)Aʳⁱᵢₖ,
    Aʲᵢₖ = (δᵢₖ/xᵢⱼ − 1 + Φʲᵢₖ)/Nⱼ,

the two-phase form, H = δⱼₘAʲ + Aˡ, when one phase l is every component's.
With one dependent phase for all, a trace such as water's 1e-29 of decane is
z − Σ of numbers near 0.3, rounding, and Newton is lost; as a variable of
its own it keeps its digits. And the variables are scaled by √n: a trace's
Hessian entry is 1/n, 1e29 beside entries of order 1, which Cholesky can't
factor, and scaled every diagonal entry is of order 1 (Michelsen and
Mollerup). The scaled gradient is s·g, so its test is held at the smallest
s, and the stop is left to a step too small to move anything or a fall G's
rounding can't show; the split holds the answer to tol in every ln f.

**The answer** has every |ln fᵢⱼ − ln fᵢᵣ| ≤ tol and every fraction above 0,
and with `verify` each phase is tested for stability, the others as
trials. An unstable phase of a split of two grows it to three; of a split
of three, it is a fourth, `.more_phases`. A phase that leaves on the way
leaves a split of two, verified again.

**PH and PS** (D10) search first on unverified PT flashes, which stay
two-phase, and if that misses the energy balance, again on verified ones:
in a region of three phases an unverified flash's H(T) is a false split's.
Newton's finish in (v, ln T) is two-phase, so a split of three is finished
by Brent alone, to etol. **TV and UV** stay two-phase: an answer of three
is refused, `.more_phases`, though the search may pass through three.

## Vapour-fraction flashes

`flash.tvf(model, T, β, z)` finds P, and `flash.pvf(model, P, β, z)`
finds T ([flash/vf.v](../flash/vf.v)). `bubble_p`, `dew_p`, `bubble_t`
and `dew_t` are the β = 0 and β = 1 cases. The answer is a FlashResult
with both phases even when one has no moles.

**The equations** (Michelsen & Mollerup, chapter 12). The unknowns are ln Kᵢ
for each component present, and τ = ln P (T given) or ln T (P given):

    Fᵢ = ln Kᵢ + ln φᵢ(y) − ln φᵢ(x) = 0,
    F_R = Σᵢ (yᵢ − xᵢ) = 0,
    xᵢ = zᵢ/dᵢ,  yᵢ = Kᵢxᵢ,  dᵢ = 1 + β(Kᵢ − 1),

with each phase solved at its normalised composition. ln φ depends on mole
numbers only through the mole fractions, so with Φ = n·∂ln φᵢ/∂nⱼ at one
mole and S the sum of each phase's unnormalised x or y:

    ∂Fᵢ/∂ln Kq = δᵢq + Φⱽᵢq·yq(1 − β)/(dq·Sy) + Φᴸᵢq·βKqxq/(dq·Sx),
    ∂Fᵢ/∂τ = P(∂ln φᵢⱽ/∂P − ∂ln φᵢᴸ/∂P)   or   T(∂ln φᵢⱽ/∂T − ∂ln φᵢᴸ/∂T),
    ∂F_R/∂ln Kq = zqKq/dq²,   ∂F_R/∂τ = 0.

Newton stops with every |F| ≤ tol, so the mass balance holds to tol. Its
last step is quadratic, so in practice F_R ends near round-off.

**The roots.** Phase x takes its liquid root and phase y its vapour root
throughout. Giving each its lower-Gibbs root instead let an incipient phase
start on the feed's kind of root, so the phases looked alike and the
iteration slid to the trivial solution K = 1, which solves these equations
at any τ. At the end each phase must be on its lower-Gibbs root, or the
boundary is between metastable phases and is refused.

**Starts.**
- Newton (with a line search) starts from Wilson's K-values at the τ where they
  satisfy Rachford–Rice. With Kᵢ = aᵢ/P that's a bisection in ln P.
  For P given, ln Kᵢ is linear in 1/T, so it's a bisection in 1/T.
- If Newton fails from there, it starts again from the K-values of the
  nearest two-phase PT flash it can find: τ is stepped out both ways, 10%
  in P or 2% in T, up to 20 steps each way.
- If no PT flash in that range splits, the failure is `outside_envelope`.
  That's the evidence for "no such point here", such as a dew point above
  the cricondenbar.
- A caller's `estimate` starts from the PT split there, if there is one.
- On the reference set:
  - Newton takes 163 steps for `thermo`'s 40 bubble and dew points, and
    the dogleg 182. (While Rachford–Rice's equation was weighted by 100,
    the line search's merit function was distorted and the dogleg did
    better; the weight bought nothing and is gone.)
  - successive substitution first, as in the PT flash, slid to K = 1 from
    Wilson's poor starts near a critical point.

**Two answers.** At a temperature between a mixture's critical point and its
cricondentherm there are two dew pressures, and near the critical point both
ends of an isobar can be bubble points. The flash returns the one Newton
reaches from its start, and `estimate` chooses another start. At 90 bar,
methane/n-butane (40/60) splits from 327 to 379 K, and both ends are bubble
points: there is no dew point at that pressure at all.

**Checks.**
- The fugacity residual must be ≤ tol (1e-12), and the iteration must not
  have collapsed onto the feed.
- Phase y, of fraction β, must come out the less dense by mass. Otherwise
  the answer's vapour fraction is 1 − β, and it is refused
  (`phases_swapped`) rather than relabelled.
- With `verify`, each phase is re-tested for stability, as in the PT flash.
  In both flashes the re-test calls a phase unstable at
  tm < −max(1e-8, 10·tol): an answer converged only to a loose tol has tm
  of about that size against its own phases.

**A pure component** has K = 1 at any answer, so it gets its own Newton in τ
alone: ln φ on the liquid root equals ln φ on the vapour root. It first
moves τ until both roots exist. With only a vapour root, P is too low or T
too high, and the reverse with only a liquid root. Any β gives the same T
and P.

**Azeotropes** (M6 step 6, [flash/azeotrope.v](../flash/azeotrope.v)).
`flash.azeotrope_t(model, P, x0)` finds the azeotrope near x0 at P, and
`azeotrope_p(model, T, x0)` at T. A vapour-fraction flash can't: its two
phases are the same composition there, which it takes for the trivial
solution. The unknowns are x over the components present but the last,
which is 1 − Σ of the rest, and τ = ln T (or ln P); the equations are equal
fugacities with y = x,

    gᵢ = ln φᵢᴸ(x) − ln φᵢⱽ(x) = 0,   one for each component present,

with Φ = n·∂ln φᵢ/∂nⱼ and l the last component,

    ∂gᵢ/∂xₖ = (Φᴸᵢₖ − Φᴸᵢₗ) − (Φⱽᵢₖ − Φⱽᵢₗ),
    ∂gᵢ/∂τ = T(∂ln φᵢᴸ/∂T − ∂ln φᵢⱽ/∂T)   or   P(∂ln φᵢᴸ/∂P − ∂ln φᵢⱽ/∂P).

Newton starts from x0's bubble point. Where x0 is an azeotrope already, or
within about 1e-4 of one, that bubble point collapses onto K = 1, the
trivial solution it refuses, so the start is the bubble point 2% from x0
towards each component in turn, Newton still starting at x0. The answer is
shown: each |gᵢ| ≤ tol, every xᵢ in (0, 1), the liquid and the vapour two
phases, not one root twice (where an equation of state has one root, every
gᵢ is 0 at any x), and the liquid stable with the vapour as a trial. Which
phase has the lower G says nothing here: at an azeotrope the two are
equal. A liquid that would split is a
heterogeneous azeotrope, refused as `.more_phases` (three phases are step
7's). The result is the vapour and the liquid, each of composition x, at
β = 0. A caller finds where to look by where y − x changes sign along the
bubble points, as the tests do.

## PH and PS flashes

`flash.ph(fluid, P, H, z)` and `flash.ps(fluid, P, S, z)`
([flash/ph.v](../flash/ph.v)) find T and the split from a molar H or S, on
the Fluid's reference state. They take a Fluid, not a phase model, because
H and S need the ideal gas.

**The Q-function** (Michelsen 1999, D5). At fixed P and H the equilibrium
is the stationary point of

    Q(T, v) = (G(T, P, v) − H_spec)/RT      (PS: (G + T·S_spec)/R),

a minimum in the vapour mole numbers v and a maximum in T. Minimising over
v at fixed T is the PT flash, and what is left, Q*(T), has

    dQ*/dT = −(H(T) − H_spec)/RT²           (PS: S_spec − S(T)),

where H(T) is the equilibrium enthalpy. H(T) at fixed P never falls as T
rises, latent heat included, so the sign of H(T) − H_spec brackets the
answer, however steep or broken H(T) is.

**The steps.**
- **Start:** the temperature at which the feed as an ideal gas has
  H_spec (or S_spec), by bisection, but no lower than half the feed's mean
  Tc, or the caller's `estimate`. A liquid's H is far below the ideal
  gas's, so that temperature can be near absolute zero, where PR's liquids
  split in ways no answer needs and the PT flash refuses.
- **Bracket:** T is stepped out by 50% at a time, with a PT flash at
  each step, until the residual changes sign. (Steps of 1.5× and 2× take
  the fewest flashes on the reference set, and 1.25× takes 3% more.)
- **Narrow:** Brent's method on the residual, to 1e-3 in T.
- **Finish:** if the PT flash there splits, Newton on Q's stationary
  conditions in (v, ln T) together, equal fugacities and the energy
  balance, takes it to round-off.
- **Fallback:** otherwise, or if Newton fails, Brent goes on.
- **The answer** is a verified PT flash at the final T, with
  |H − H_spec| ≤ etol·RT or |S − S_spec| ≤ etol·R (etol 1e-10).
- **Only the answer is verified** (M4 step 9). The PT flashes of the bracket
  and Brent's steps only steer T, and verifying a split's phases costs two
  more stability tests each. A split there that isn't the equilibrium can
  move the search but not the answer, which is verified wherever it ends.
  One phase has nothing to verify, so a single-phase answer is taken as it is.

**Q's Hessian**, Newton's Jacobian. With gᵢ = ln fᵢⱽ − ln fᵢᴸ,
r_E = (H − H_spec)/RT (or (S − S_spec)/R) and
cᵢ = T·∂(ln φᵢⱽ − ln φᵢᴸ)/∂T at fixed composition:

    ∂gᵢ/∂vⱼ = the PT flash's Hessian (above),   ∂gᵢ/∂ln T = cᵢ,
    ∂r_E/∂vⱼ = −cⱼ             (PS: −gⱼ − cⱼ),
    ∂r_E/∂ln T = Cp/R − r_E    (PS: Cp/R),

with Cp the two phases' heat capacity at fixed composition. They follow
from H = Σᵢ(vᵢh̄ᵢⱽ + lᵢh̄ᵢᴸ) and h̄ᵢⱽ − h̄ᵢᴸ = −RT²·∂(ln φᵢⱽ − ln φᵢᴸ)/∂T.
For S there is also the mixing entropy:
s̄ᵢⱽ − s̄ᵢᴸ = −R·gᵢ − RT·∂(ln φᵢⱽ − ln φᵢᴸ)/∂T.

**Brent's method** ([numeric/roots.v](../numeric/roots.v)) is inverse
quadratic interpolation or the secant where they make progress, and
bisection where they don't (Brent 1973). On a smooth function it converges
superlinearly, and it is never much slower than bisection, so a jump in
H(T) only narrows the bracket onto it. It's slow on a multiple root, as
Brent's method is, but a flash's residual has none. Its tolerance is
relative (4ε|x|, plus xtol/2), so a root at 0 needs an xtol.

**A pure component's** H jumps at Tsat by the heat of vaporisation, so no
T gives an H inside the jump. There the Q-function's maximum is at the
kink, T = Tsat, with the vapour fraction the energy balance gives:

    β = (H_spec − Hᴸ(Tsat))/(Hⱽ(Tsat) − Hᴸ(Tsat)).

The flash finds Tsat with pvf, then answers this way for an H_spec
between the saturated liquid's and vapour's, within etol of either end.
Outside that range the answer is a single phase on its own side of Tsat,
found by Brent on that root alone.

**A model that refuses a T** (M6 step 6). A γ–φ liquid is refused outside
its correlations' ranges (D4), and a flash says so as `no_volume_root`. The
bracket then steps back from that T towards the last, by halves, up to 40
times, and goes on from the farthest T the model takes. Only if the model
takes none is the specification outside its range. An equation of state
never refuses a phase there, so its brackets don't move.

**Newton's own state** (M6 step 6). Within a hundredth of a kelvin of an
azeotrope, H rises so steeply with T that the PT flash at the nearest double
to Newton's T can miss the balance: methanol–water at 1.39e-10 against
1e-10. Newton's own (v, ln T) meets it. So where the answer's balance is
off and the PT flash split, the answer is Newton's state from Brent's last
split, if it is one: each phase on its phase of lower G, the fugacities
equal to tol, the balance met to etol, both phases' amounts positive, and
with `verify` each phase stable with the other as a trial. The vapour is
the phase of lower mass density, as in a PT split.

## TV and UV flashes

`flash.tv(model, T, V, z)` finds P and the split from a molar volume, and
`flash.uv(fluid, U, V, z)` finds T, P and the split from a molar internal
energy and volume ([flash/tv.v](../flash/tv.v)). TV needs only a Model,
since the ideal gas's part of A is the same for any split at fixed T, and
nothing less, since it works in F(T, V, n): a caller with a fluid passes
`fluid.helmholtz()!`. UV refuses a fluid without one, saying why.

**Q-functions** (D5). At fixed T and V the equilibrium minimises A, and
Q(P) = G − P·V_spec has dQ*/dP = V(P) − V_spec, with V(P) the equilibrium
volume, which falls as P rises. At fixed U and V it maximises S, and
Q(T) = (A − U_spec)/T has dQ*/dT = −(U(T) − U_spec)/T², with U(T) the TV
flash's U, which rises with T. So each has a monotone residual to bracket,
as the PH flash does.

**TV, step by step.**
- **The feed as one phase:** its pressure P(T, V) comes from the model
  directly, and a PT flash there confirms the phase is stable and has this
  volume. If that PT flash splits instead, its split starts the Helmholtz
  minimisation below.
- **Otherwise,** ln P is stepped by ln 2 with PT flashes until
  V(P) − V_spec changes sign, and Brent narrows that to 0.1%.
- **The Helmholtz minimisation,** in the model's own variables. Over the
  vapour's mole numbers nⱽ and volume Vⱽ, with nᴸ = z − nⱽ and
  Vᴸ = V − Vⱽ, it minimises

      A/RT = Σ over phases [F(T, V, n) + Σᵢ nᵢ(ln(nᵢ/V) − 1)] + a constant,
      ∂/∂nᵢⱽ = (F_nᵢ + ln(nᵢ/V))ⱽ − (…)ᴸ = (μᵢⱽ − μᵢᴸ)/RT,
      ∂/∂Vⱽ = (F_V − n/V)ⱽ − (…)ᴸ = (Pᴸ − Pⱽ)/RT,

  with the Hessian the two phases' sums of F_nn + δᵢⱼ/nᵢ, F_Vn − 1/V and
  F_VV + n/V². It needs no volume root: only F(T, V, n) and its
  derivatives.
- **The start:** the PT split, the liquid keeping its volume, since it
  hardly compresses, and the vapour taking the rest of V.
- **The answer** is a verified PT flash at the pressure found, held to
  |V − V_spec| ≤ vtol·V_spec (vtol 1e-10).

**UV, step by step.**
- **Bracket:** ln T is stepped by ln 1.5 with TV flashes until
  U − U_spec changes sign, from the ideal gas's T for U_spec. That T can
  be far too low for a liquid's U, low enough that no TV flash at this V
  has an answer; then the search starts from 300 K.
- **Narrow:** Brent, to 0.1%.
- **Finish:** Newton.
  - For a split it works in (nⱽ, Vⱽ, ln T): the Helmholtz conditions plus
    r_U = (U − U_spec)/RT, with
    U = Σᵢ zᵢ(ΔH°f,ᵢ + ∫Cpᵢ dT − RT) − RT²·Σ F_T. With
    cᵢ = T(F_Tnᵢⱽ − F_Tnᵢᴸ) and c_V = T(F_TVⱽ − F_TVᴸ), the extra
    derivatives are
    ∂gᵢ/∂ln T = cᵢ, ∂g_V/∂ln T = c_V, ∂r_U/∂nᵢⱽ = −cᵢ, ∂r_U/∂Vⱽ = −c_V
    and ∂r_U/∂ln T = Cv/R − r_U, with Cv = Σᵢ zᵢ(Cpᵢ − R) − RT·Σ(T·F_TT +
    2F_T).
  - For one phase it works on U(T) at this V, from F alone.
- **The answer** is a verified TV flash at the final T, held to
  |U − U_spec| ≤ etol·RT.

**Why not start with the one-phase Newton.** Newton on the one phase's U(T)
from the start would save the bracket, but it can run far. From a split's
U it found a one-phase answer at 5035 K, where the test data's quadratic
Cp has turned negative (n-butane's at 3290 K) and U(T) falls again. The
bracket finds the nearest change of sign instead.

**A pure component** inside its saturation dome has P = Psat(T), from tvf,
and β = (V − Vᴸ)/(Vⱽ − Vᴸ). Outside it, the one phase's pressure is P(T, V),
on the root of that volume.

**Conditioning.** A liquid's P moves V little, so V moves P a lot: by the
bulk modulus over P, a factor of a thousand here and there. `thermo`'s V
differs from reflux's by about 1e-11, which is 1e-8 in a liquid's P, so
the round trips start from reflux's own V and U.

## The robustness sweep

`sweep` runs PT flashes of random PR mixtures (D8) and checks each answer
from scratch, without trusting the flash's own tests.

**The points** come from splitmix64 (Steele, Lea & Flood 2014), seeded, so
they're the same on every platform:
- a u64 to f64 in (0, 1) is (⌊u/2¹¹⌋ + ½)/2⁵³, never 0;
- 2 to 10 compounds, drawn without replacement from the 640 with Tc, Pc,
  ω, M and an ideal-gas Cp;
- T = T_min·(T_max/T_min)^u, and P likewise: log-uniform over 100–1000 K
  and 0.1–1000 bar;
- zᵢ = eᵢ/Σe with eᵢ = −ln uᵢ, exponential draws, which makes z uniform
  over the simplex.

**Stability is judged from scratch.** `flash.stability` only searches, from
its own trials and five random compositions drawn the same way. At every
composition w its trials reach, the sweep computes

    tpd(w) = Σᵢ wᵢ(ln wᵢ + ln φᵢ(w) − ln xᵢ − ln φᵢ(x))

from the model's phases alone (`thermo.solve_tp`'s, for a cubic), each at
its root of lower G, and the
lowest tpd is the verdict. A negative tpd at any w shows that x lies above
some phase's tangent plane, so x is unstable. A phase of a converged split
is no use as a trial for the other: it sits on that phase's tangent plane,
a stationary point with tm = 0.

**One phase** must have no tpd below −1e-10, the flash's own tolerance.

**A split** must satisfy criterion 3:
- 0 < β < 1;
- the phases differ, some |ln(yᵢ/xᵢ)| ≥ 1e-4, since two copies of a stable
  feed pass every other check;
- |β·yᵢ + (1 − β)·xᵢ − zᵢ| ≤ 1e-12;
- |ln yᵢ + ln φᵢ(y) − ln xᵢ − ln φᵢ(x)| ≤ 1e-10, except where xᵢ or yᵢ is
  below `trace_floor`;
- |Σx − 1| and |Σy − 1| ≤ 1e-14;
- neither phase with a tpd below −1e-8.

**A third phase** (`more_phases`) is proved, not taken on trust. Unverified,
the flash returns the first split it found, which is the one the verified
flash found first and refused. One of its phases must have a tpd below
−1e-8.

## Units of measure

`quantity` is where numbers with units come in (M4's D10). Since M9 (its
D3) it has a dimension algebra, so a unit can be compound and every unit
knows its dimension.

**A dimension** is the exponent of each SI base unit, kg, m, s, K, mol and
A: a pressure is kg·m⁻¹·s⁻². The exponents are exact rationals p/q in
lowest terms, with |p| and q at most 10⁶, so that bar^0.5 squared is bar
again exactly, and no sum or product of two in range overflows 64 bits; one
out of range is refused. The ampere is there for the coulomb, and so for the
Faraday constant. Each of M4's named dimensions is one of these, and a
difference has its absolute quantity's: ΔT is K, ΔP is kg·m⁻¹·s⁻². A
dimension is written in base units, its positive powers over its negative
ones, as kg/(m·s²), mol/(kg·s), kg^(1/2)/(m^(1/2)·s) or 1, and that text
reads back as the same dimension with scale 1.

**A unit** is a scale on a dimension. Each symbol is a row of one table,

    SI = (x + pre)·scale + post,        x = (SI − post)/scale − pre

A scaled unit has pre = post = 0. Celsius and Fahrenheit shift before they
scale, and gauge pressures add a standard atmosphere after:

    T/K  = t/°C + 273.15
    T/K  = (t/°F + 459.67)·5/9
    T/K  = (T/°R)·5/9
    P/Pa = (P/barg)·10⁵ + 101325
    P/Pa = (P/psig)·6894.757… + 101325

A temperature or pressure *difference* uses the same scales with no
offset: 9 degF of ΔT is 5 K, and there are no gauge differences.

The scales that aren't SI prefixes, minutes or hours are built from the
exact definitions in NIST SP 811 (2008), appendix B. None is a rounded
number typed in:

| Unit | Definition |
|---|---|
| pound | 0.45359237 kg |
| inch | 0.0254 m, so the foot is 0.3048 m |
| standard gravity g₀ | 9.80665 m/s² |
| psi | lbf/in², with lbf = pound·g₀: 6894.757293… Pa |
| atm | 101325 Pa; the torr is 1/760 of it, 133.322368… Pa |
| conventional mmHg | 13.5951 g/cm³ · g₀ · 1 mm = 133.322387415 Pa, 1.4e-7 above the torr |
| IT Btu | 1055.05585262 J, so Btu/lb is exactly 2326 J/kg |
| IT kilocalorie | 4186.8 J |
| US gallon | 231 in³ |
| horsepower | 550 ft·lbf/s = 745.699871… W |
| pound-mole | 453.59237 mol: a pound's worth of moles |
| IT calorie | 4.1868 J |
| metric ton (t) | 1000 kg |
| day (d) | 86400 s |
| poise (P) | 0.1 Pa·s, for the centipoise |

SP 811 prints some of the factors built from them (psi, mmHg, torr, lb/h,
gal/min, Btu, Btu/h, hp) to 7 figures, and reflux's values round to them.

The symbols are the SI units m, kg and g, s, K, mol, A, C, L, N, J, W, Pa
and P; in, ft, t, lb, min, h, d, lbmol, gal, cal, kcal, Btu, hp, bar, atm,
psi, psia, mmHg and torr; the temperatures degC, °C, degF, °F and R; the
gauge pressures barg and psig; and %. Of these, m, g, s, K, mol, A, C, L,
N, J, W, Pa, bar and P take the prefixes G, M, k, h, d, c, m, µ (or μ, or
u) and n. The kilogram is a symbol of its own, so that it is exactly 1. Pico and smaller are left out, so
that "ps" can't read as a picosecond where psi was meant.

**A compound unit** is read by

    unit     := factor (("*" | "·") factor | "/" factor)*
    factor   := ("1" | symbol | "(" unit ")") [exponent]
    exponent := digits right after a symbol (m3), superscripts (m³, s⁻¹),
                or "^" and an integer, a decimal or "(p/q)" (bar^-2, bar^0.5)

- **Its scale** is the product of its numerator's scales over the product of
  its denominator's, in one division at the end, as M4's table computed
  kmol/h = 1e3/3600. A whole power is a repeated product, so cm³ is
  0.01·0.01·0.01, the same bits on every platform; a fractional one is
  `math.pow`, the C library's. (From the fourth power, `math.pow` is
  sometimes an ulp nearer the exact value, but Windows' and Linux's
  libraries differ.) So every unit M4 read, read now through the algebra,
  gives M4's bits.
- **Division goes left to right,** so J/mol/K is J/(mol·K). A * after a / is
  refused as ambiguous: kmol/kg*h might mean kmol/(kg·h) or kmol·h/kg.
- **Offsets only standing alone.** °C, °F, barg and psig keep their offsets
  only as a whole unit of an absolute temperature or pressure. In a compound
  unit °C is a kelvin of difference and °F 5/9 K, as in kJ/(kg·°C), and a
  gauge pressure is refused.
- A scale must be finite and above zero, and parentheses nest at most 16
  deep.

Reading is strict:
- the text is a number, whitespace and one unit;
- the number is what `strconv.atof64` reads, and must be finite;
- only a fraction may be a bare number;
- an absolute T or P must be above zero.

A unit that isn't one gets the closest of the dimension's common units
(M4's): those within two edits (Levenshtein, ignoring case), nearest first,
at most three. A unit of another dimension is named by what it is a unit
of: a named dimension, a common one (energy, viscosity, molar entropy, …),
or else its base units, as in "kmol/(kg*h) is a unit of mol/(kg·s), not
molar flow". A unit that doesn't parse says why and at which character.

## Expressions

`expr` is reflux's one expression language (M9's D1, D2 and A9): a rate
law, a rate constant, an adsorption term, later a pressure drop or a design
specification. An expression is parsed into a tree and never evaluated as
text. Every name in it is bound, every node has a dimension, and the tree
writes itself back.

**The grammar,** loosest first:

    expr   := term (("+" | "-") term)*
    term   := factor (("*" | "·" | "×" | "/") factor)*
    factor := ("-" | "+") factor | power
    power  := atom ["in" unit] ["^" factor]
    atom   := number [unit] | name | name "[" key "]" | function "(" expr ")"
              | "(" expr ")" ["[" unit "]"]
    name   := word ["." word] | quoted "." word

So ^ binds tighter than a minus sign and to the right, and · and / go left
to right, as in Python and SymPy: −x^2 is −(x²), x^-2 is x⁻², a^b^c is
a^(b^c), and a/b·c is (a/b)·c. − (U+2212) is a minus too. Whitespace,
line breaks included, goes anywhere between tokens.

**A number's unit** follows it after a space and runs to the next space, or
to a ) or , that isn't its own: `exp(-240.1 kJ/mol / (R*T))`. Or it's in
brackets, `240.1 [kJ/mol]`. The number becomes its SI value as it's read,
by `quantity`'s scale, and its node has the unit's dimension. A word after a
number that is also a name here is refused rather than guessed: "2 R" could
be two rankines or twice the gas constant, so it asks for 2*R or 2 [R]. A
unit with an offset (°C, °F, barg, psig) after a number is refused, since in
arithmetic its offset would be wrong; in a compound, °C is a kelvin of
difference and reads. Two values side by side, with no operator, are
refused: "2T", "2 T", "T (P)".

**A reading** (D2, as amended before M9 step 10) is a value's number on a
unit's scale, a pure number, so that a correlation goes in as its source
prints it: `0.85 - 0.0018*(T in degC)` is the tutorial's X = 0.85 − 0.0018 T
with T in °C, and `log10(P in mmHg)` is Antoine's. With `quantity`'s
SI = (x + pre)·scale + post, the reading of a value v is

    x = (v − post)/scale − pre,

T in °C being T − 273.15 K over 1 K, and P in barg (P − 101325 Pa)/10⁵ Pa.
It's the one place an offset applies, and these rules keep it there:
- `in` takes the name, member or bracketed group just before it, as a unit
  binds to its number: `0.0018*T in degC` is `0.0018*(T in degC)`, and
  `T in degC ^ 2` squares the reading. A number before `in` is refused,
  as "2 in degC" (after a number, `in` alone is still the inch).
- The unit must have the value's dimension.
- A unit with an offset reads only a context variable, T or P, which is
  absolute: `(T - 300 K) in degC` is refused, since a difference reads in
  K, and so are a family member and a definition, which may be differences.
- A reading renders as `(T in degC)` in plain text, `(T in °C)` in Unicode
  and `\left(T\ \mathrm{in}\ \mathrm{{}^{\circ}C}\right)` in LaTeX,
  always in parentheses, so that it reads as one value wherever it stands.
- `in` is reserved, so no definition can take it. A reading that repeats
  is named once by a definition, `t = "T in degC"`; there is no declaration
  for a whole expression, which would rebind T and let `T/300` through.
- A sum of a pure number and a temperature or pressure, or a pure number
  asked for and a temperature given, is refused with the reading
  suggested: likeliest, a correlation typed as printed without one.

The tape computes a reading as the tree does, the three operations in
order, skipping only x − 0 and x/1, which change no bit; its derivative is
1/scale.

**A group may take a unit** in brackets after it, `(85 - 0.18*(T in degC))
[%]`: the group's value times the unit's scale, with the group's dimension
times the unit's, as a number's unit does for a number. Only in brackets,
so that a word after a group is still two values side by side; a unit with
an offset is refused there as after a number.

**Names** are letters, digits and _, starting with a letter, and bind to:
- the context's variables, such as T and P, and its families' members, by
  any key the context gives a member, so p[CH4], p[Methane] and
  p['74-82-8'] are the same input. A key with other characters is quoted;
- a definition, by name;
- a built-in constant, at the 2019 SI's exact value:

| Name | Value | Unit |
|---|---|---|
| R | N_A·k_B = 8.31446261815324, exactly | J/(mol·K) |
| N_A | 6.02214076e23 | 1/mol |
| k_B | 1.380649e-23 | J/K |
| h | 6.62607015e-34 | J·s |
| F | N_A·e = 96485.3321233100184…, exactly | C/mol |
| c | 299792458 | m/s |
| g_n | 9.80665, standard gravity | m/s² |
| p_std | 1e5, the standard pressure, 1 bar | Pa |
| pi | π | |

- the functions exp, ln, log10, sqrt and abs, each of one argument. log is
  refused as ambiguous, with ln and log10 offered.

A name that binds to nothing is refused, with the closest names within two
edits.

**A flowsheet's names** (M8's D3) have one dot: a stream's or a unit's
name, then one of its variables, as the file names them. A scope that
isn't a word is quoted, `'let down'.T`, and no space goes either side of
the dot. Each is an input slot, filled from a solution:

| Name | Dimension | Its value |
|---|---|---|
| `s.T`, `s.P` | temperature, pressure | the stream's |
| `s.F` | molar flow | Σ nᵢ |
| `s.n[X]` | molar flow | X's flow nₓ |
| `s.y[X]` | none | nₓ / Σ nᵢ, of the whole stream; none where Σ nᵢ = 0 |
| `s.h` | molar enthalpy | the stream's, on the fluid's reference |
| `s.vf` | none | its vapour fraction, β |
| `u.duty`, `u.work` | power | the unit's energy stream of that name |

A component answers to the file's name, its compound's name, CAS number and
aliases, and its formula where no other component has it, which is what it
renders as; otherwise the file's name. A unit names its energy streams
before it calculates (`Energies`), so an adiabatic reactor has no duty to
name, and a mixer nothing. A name that isn't there is refused by what's
wrong: a scope that isn't there ("no stream or unit prodct here", with
product.T offered), a variable its scope hasn't ("stream product has no TT:
it has T, P, F, h, vf, n[…], y[…]"), a scope alone ("product is a stream:
name one of its variables"), or a variable without its scope ("n is a
stream's: say whose, as feed.n[CH4]"). Plain text writes a quoted scope
back quoted; LaTeX writes a scope upright, `\mathrm{let\ down}.T`.

**A value, or why there's none** (M8's D9). `Expr.value` computes what `eval`
does, in its order and to its bits, but refuses, at the node and in words,
where `eval` would give ±∞ or NaN:

| Node | Refused when | Saying |
|---|---|---|
| an input | it's NaN | it has no value here, and why if the caller knows: "gas has no flow" |
| a/b | b = 0 | it divides by b, which is zero here |
| ln a, log10 a | the result isn't finite (a ≤ 0) | a, and that a logarithm needs a number above zero |
| sqrt a | a < 0 | a, which is negative here |
| a^b | 0^(b < 0), or a < 0 and b not whole | which |
| anything else | the result isn't finite | that it's too large for a number |

Tutorial 2's product conversion, (in − out)/in with none in, is the first
division's refusal, not −∞.

**Equations** are two expressions with `=`, `<=` or `>=` between them (`≤`
and `≥` too): a spec's target (D5) or an optimisation's constraint (D13).
The relation is found outside brackets and quotes; none, two, a strict `<`
or `>`, or `==` is refused. Each side is parsed with the other blanked, so
its lines and columns are the whole text's, and the two must have one
dimension.

**Definitions** name expressions for reuse, in any order. Each may use the
others but not itself through any chain, found by a depth-first walk that
reports the chain: "u uses itself: u → v → w → u". A definition can't take
the name of a function, a constant, a variable or a family. Definitions are
checked in an order where each follows what it uses.

**Dimensions,** node by node, in the tree's order, where every operand
comes before its use:

| Node | Dimension |
|---|---|
| number | its unit's, or none |
| a + b, a − b | a's, which must be b's |
| a·b, a/b | a's times or over b's |
| −a, abs(a) | a's |
| sqrt(a) | a's to the ½ |
| exp(a), ln(a), log10(a) | none, and a must have none |
| a^b | b must have none. If a has none, none. If a has one, b must be a fixed value, folded from numbers, constants and fixed definitions, and a fraction p/q with q ≤ 10⁶ within 1e-12 (by continued fractions: 1/3, 2.5, 1.37), and the result is a's to the p/q |

A place that needs a dimension asks with `must_be`: "a rate per mass of
catalyst must be a mol/(kg·s), but this is a mol/(m³·s)". DWSIM's
`exp(-240100/8.314/T)` is refused here: its argument is in 1/K.

**A number that looks like a constant** is flagged, not refused: one with a
constant's dimension within 1% of it but not equal, such as 8.314
J/(mol·K) for R, with how far it is. A paper fitted with 8.314 may still
write it (D2). The standard pressure is left out, since a pressure near
1 bar, 1 atm say, usually means what it says.

**Writing back.** Three styles from one walk, with the same precedence:
- **plain,** in the language's own syntax, which reads back to the same
  tree, node for node. So it keeps the tree's grouping, not just its value:
  a + (b + c) and a·(b/c) keep their parentheses, being other trees than
  a + b + c and a·b/c;
- **Unicode,** for a terminal: · for a product, superscripts for a whole
  power (x², x⁻¹), a trailing number in a name as a subscript (k₁), a
  formula's digits as subscripts (p_CH₄), Greek names as letters (θ), π, p°
  and gₙ. 4.22e15 is 4.22·10¹⁵, so it binds as a product and keeps its
  parentheses as a divisor or a base. A number with a unit keeps its
  parentheses as an operand of ·, / or a minus sign, so that its unit never
  runs into the operator after it: (4.22·10¹⁵ kmol/(kg·h))·exp(−(240.1
  kJ/mol)/(R·T));
- **LaTeX,** for a notebook: \frac, \cdot, ^{}, \left( \right), \exp, \ln,
  \log_{10}, \sqrt, \left| \right|, 4.22 \times 10^{15}, units in
  \mathrm{}, species as p_{\mathrm{CH_{4}}}.

**Powers** are a^b by `raise`: a whole power up to 16 as a repeated
product, a negative one as 1 over that, and any other by `math.pow`. The
product's bits are the same on every platform, where `math.pow` is the C
library's, and the tree and the tape both raise this way.

### The tape

Expressions are compiled for speed into a tape (D5, as step 0 amended it):
a flat array of operations on slots, one slot per operation, every operand
before its use. Several expressions compile onto one tape, so what they
share is computed once.

**Compiling,** node by node:
- **Folding:** an operation whose operands are all constants is computed
  as eval would, and becomes a constant.
- **Merging:** an operation identical to one already there (its kind, its
  operands, its constants) is that one. + and · put their operands in
  order first, which IEEE arithmetic allows exactly. A definition is
  compiled once, whatever uses it, so a chain of definitions each using the
  next twice compiles in linear time.
- **Fusing,** each keeping the tree's arithmetic to the bit:

| Tree | Tape | Why the bits are the same |
|---|---|---|
| x + c, c + x | addc(x; c) | + commutes exactly |
| x − c | addc(x; −c) | x − c is x + (−c), and negation is exact |
| c − x | csub(x; c) | the same operation |
| c·x, x·c | mulc(x; c) | · commutes exactly |
| x/c, c/x | divc(x; c), cdiv(x; c) | the same operations |
| exp(d/x) | aexp(x; 1, d) | 1·y is y |
| c·exp(d/x) | aexp(x; c, d), if the exponential isn't scaled already | the same operations |
| exp(c + d/x), exp(c − d/x) | expl(x; c, ±d) | c − d/x is c + (−d)/x, exactly |
| a^c, c fixed | powc(a; c) | `raise` both ways |

  So Arrhenius' A·exp(−E/(R·T)) is one operation on the slot of R·T, which
  every such term shares, and an equilibrium constant's exp(30.42 −
  27106 K/T) is one operation on T's.
- **Pruning:** what no output uses is dropped. Then the constants go
  first, the inputs used next, and the operations after, in order.

**Evaluating** writes each slot's value into a buffer the caller owns,
whose constants were put in place once when it was made. Nothing is
allocated. The tape's value is eval's, to the bit.

**The gradient** comes from one reverse sweep after an evaluation, carrying
an adjoint per output in each slot, so every output's derivative with
respect to every input comes at once. Each operation adds its adjoint times
its local derivatives to its operands':

| Operation, v its value | ∂v/∂a | ∂v/∂b |
|---|---|---|
| a + b, a − b | 1 | 1, −1 |
| a·b | b | a |
| a/b | 1/b | −v/b |
| a^b | b·v/a (b·a^(b−1) at a = 0) | v·ln a (0 where v = 0) |
| −a | −1 | |
| exp(a), ln(a), log10(a) | v, 1/a, 1/(a·ln 10) | |
| √a, \|a\| | 1/(2v), sign(a) | |
| a^c | c·v/a (c·a^(c−1) at a = 0) | |
| a + c, c − a, c·a, a/c, c/a | 1, −1, c, 1/c, −v/a | |
| c·exp(d/a), exp(c + d/a) | −v·d/a² | |

An input no output uses has a derivative of 0.

## Reactions

M9 step 4: a compound's formula, equations read and balanced, a set's
independence, and each reaction's standard properties and K(T).

### Formulas

[chem/](../chem/chem.v) knows the elements: their symbols, and CIAAW's
standard atomic weights of 2021, abridged, with the conventional value for
an element whose weight is an interval. A formula is element counts, by
atomic number, read from:

- **a formula,** element symbols each with a whole count after it if not 1,
  and groups in ( ) or [ ] with a count after them: CH3OH, C2H5OH,
  CH3(CH2)2CH3. Case matters, so Co is cobalt and CO carbon monoxide. A
  fractional count is refused: no compound has a fraction of an atom.
- **a SMILES string** ([chem/smiles.v](../chem/smiles.v)): atoms in the
  organic subset (B C N O P S F Cl Br I, aromatic b c n o p s) and in
  brackets (isotope, element, chirality, hydrogens, charge, class), bonds,
  branches, ring bonds (`1`, `%12`) and dots. An organic-subset atom takes
  OpenSMILES' implicit hydrogens: enough to bring the sum of its bond orders
  v to its lowest normal valence that holds them (B 3; C 4; N 3, 5; O 2;
  P 3, 5; S 2, 4, 6; halogens 1), and none if none does. An aromatic atom
  gives one bond to its ring, so it takes max(0, V₀ − v − 1), V₀ its lowest
  valence: benzene's c one, pyridine's n none, furan's o and thiophene's s
  none. A bracket atom has the hydrogens it says. The net charge is the
  bracket atoms' charges summed.

A formula is written in Hill's order: carbon, then hydrogen, then the rest
by symbol; without carbon, every element by symbol, hydrogen among them.

**A compound's formula** ([components/formula.v](../components/formula.v))
comes from its file's three sources: its structure (ChemSep's, with the -
that marks a ring's ends and the spaces between groups dropped), its SMILES
string, and its molar mass. It's the formula two of them agree on: the
structure's and the SMILES's if they're equal, otherwise whichever of them
weighs the file's molar mass, to

    |Σ nₑ Aₑ − M| ≤ 2·10⁻⁴ M,

which is the most ChemSep's older atomic weights move a mass (sulfur's
32.065 against 32.06). A hydrogen more or less is beyond that below
5,000 g/mol. A compound with no formula two agree on, a mixture's
pseudo-compound, can't react.

### Equations

[kinetics/equation.v](../kinetics/equation.v) reads an equation, `CH4 +
H2O = CO + 3 H2`: terms joined by + on each side of one =, each a species
with a coefficient before it if not 1. An arrow is refused: whether a
reaction runs both ways is its rate law's to say.

- **A coefficient** is a whole number, a decimal or a fraction, read
  exactly, and may touch its species (3H2) unless the whole term is a
  species' name (1-butene).
- **A species** is found ([kinetics/species.v](../kinetics/species.v)) by
  its name or an alias, ignoring case, by its CAS number, or by its
  formula, read as above, so CH3OH, CH4O and HOCH3 are methanol. A name
  that two species answer to, or a formula two share (ethanol's and
  dimethyl ether's), is refused, and so is a name not found, with the
  closest names and formulas. Quotes hold a name with a + or an = in it.
- Each species appears once: twice on one side, or on both sides, is
  refused.

The coefficients ν are signed, below zero on the left, and exact.

**Balance.** For each element e, and for the charge z where any species
has one,

    Σ_left |νᵢ| nₑᵢ = Σ_right νᵢ nₑᵢ,

exactly. An equation that doesn't balance is refused with each element's
counts on each side, and a suggestion from the element matrix's null
space. With A the matrix of counts (an element per row, the charge last, a
species per column), the balanced coefficients are the vectors ν with
A·ν = 0, found by exact elimination to reduced row echelon form:

- **one direction:** scaled to the smallest whole numbers, signed so that
  the first species keeps its side, and suggested. A species with 0 is left
  out, and one whose sign changed is moved, each said;
- **none:** no coefficients would balance it. An element on one side only
  is named;
- **more than one:** the species balance in that many independent ways, so
  there's no one equation to suggest.

### Independence

A set of reactions is independent if its stoichiometric matrix (a reaction
per row) has full rank. [kinetics/dependence.v](../kinetics/dependence.v)
eliminates one reaction at a time, exactly, against the independent ones
before it, carrying each row's makeup, the combination of reactions it is.
A reaction reduced to nothing is a combination of earlier ones:

    0 = Σⱼ mⱼ νⱼ,  mₖ = 1   ⇒   νₖ = Σ_{j<k} (−mⱼ) νⱼ

so the owner's overall reforming reaction is steam reforming plus the
shift. That's fine for kinetics, each reaction having its own rate, and
refused for an equilibrium over reactions (step 7).

### Standard properties and K(T)

At a standard state of the ideal gas at 1 bar (D7),
[kinetics/standard.v](../kinetics/standard.v) takes each species' H and S
as thermo's `Fluid` does,

    H°ᵢ(T) = ΔH°f,i + ∫_T°^T Cp,i dT,    S°ᵢ(T) = S°ᵢ + ∫_T°^T Cp,i/T dT,

with T° = 298.15 K, and

    ΔH° = Σ νᵢ H°ᵢ,   ΔS° = Σ νᵢ S°ᵢ,   ΔG° = ΔH° − TΔS°,
    ln K = −ΔG°/RT,   d ln K/dT = ΔH°/RT²,

with R = N_A·k_B = 8.31446261815324 J/(mol·K), expr's built-in constant,
since K appears in rates that use it. K is on activities p/1 bar.

**What a rate's workspace keeps** (M8's D11,
[kinetics/evaluate.v](../kinetics/evaluate.v)). The two integrals are
quadratures for most compounds: 521 of the database's 640 give Cp as
ChemSep's equation 16, A + exp(B/T + C + D·T + E·T²), which has no
antiderivative in closed form. So `Kinetics.fill` computes each K only
when T has moved since the workspace's last fill, which an isothermal
reactor's never does, and at a new T it integrates each species once,
however many reactions it's in. Each reaction then sums its species' H° and
S° exactly as above, in the same order, so a K kept or shared is the K
`Reaction.standard` gives, to the bit. A fill that fails keeps nothing.

**The check against ΔG°f.** In a balanced equation the elements' entropies
cancel, so a data set whose ΔG°f agrees with its ΔH°f and S° gives

    Σ νᵢ ΔG°f,i = ΔH°(T°) − T°ΔS°(T°),

and `formation_gap` is the difference. For a compound alone
([components/formation.v](../components/formation.v)) the gap is

    ΔG°f − [ΔH°f − T°(S° − Σₑ nₑ S°ₑ)],

with the elements' entropies in their reference states from CODATA's key
values (H₂, N₂, O₂, F₂ and Cl₂ gases, Br₂ liquid, I₂ solid, each per atom;
graphite, rhombic sulfur, silicon, liquid mercury, and the noble gases).
The data-quality report lists every compound over 1 kJ/mol.

**ChemSep's standard states are mixed.** Its S° for H₂, H₂O, CO, CO₂, N₂ and
NH₃ are 0.10 to 0.11 J/(mol·K) below thermo's 1-bar values, as 1-atm values
are (R ln 1.01325 = 0.109), while methane's and methanol's are at 1 bar.
reflux takes the files as they are, so reforming's K is about 4% below
what consistent 1-bar data would give. The differences from thermo's data
are recorded in the tests.

### Rate laws and sites

M9 step 5 ([kinetics/law.v](../kinetics/law.v),
[kinetics/kinetics.v](../kinetics/kinetics.v)), as
[kinetics.md](kinetics.md) designs it (D4, D6).

**A rate's context** has T and P; each species' partial pressure p,
fugacity f, concentration c and mole fraction y; and each reaction's K at
1 bar from its Gibbs energies (D7), `Keq[reforming]`, a pure number. Not
`K`, which would make every number in kelvins ambiguous by D2's rule. A
species is a member by the name it's written with in the equations, so
methanol written CH3OH renders as p_CH₃OH, and answers to that, its name,
its aliases, its CAS number and its formula in Hill's order where no other
species shares it.

**Structured forms write text.** Each form keeps its structure, for its
checks and for a site every reaction on it shares, and writes the
expression a user would type for it, which `expr` parses. So a structured
rate and the same rate typed are the same tree (D4), and the constants keep
the numbers and units as written, rendering as the paper prints them. Each
constant becomes a definition:

| Form | Text |
|---|---|
| Arrhenius | A·exp(−Eₐ/(R·T)), or centred: k_ref·exp(−Eₐ/R·(1/T − 1/T_ref)) |
| van 't Hoff | A·exp(−ΔH/(R·T)), or centred; a ΔH below zero written as exp(\|ΔH\|/(R·T)) |
| fixed, written | the number, or the expression, as given |
| a site's term | K·Π aᵢ^eᵢ, a negative power divided: K_H2O·p_H2O/p_H2 |
| a site | θ = 1 + Σ terms, defined once (`theta_` and its name, or its symbol) |

**The driving force,** with orders α and coefficients ν:

| Form | k·(driving force) |
|---|---|
| irreversible | k·Π aᵢ^αᵢ |
| approach (Vanden Bussche and Froment) | k·Π aᵢ^αᵢ·(1 − Q/K), Q = Π aᵢ^νᵢ |
| difference (Xu and Froment) | k·Π_{αᵢ<0} aᵢ^αᵢ·(Π_{αᵢ>0} aᵢ^αᵢ − Π aᵢ^(αᵢ⁺+νᵢ)/K) |

The two reversible forms are the same function: the difference's backward
term is its forward one times Q. **K in the rate's activities** is
K°·p°^Δn on partial pressures or fugacities, K°·(p°/RT)^Δn on
concentrations and K°·(p°/P)^Δn on mole fractions, with K° = Keq[...] and
p° = `p_std`; or a correlation the paper fitted, refused unless its
dimension is the activity's to the Δn (bar² for steam reforming).

**The forms:**
- **power law:** k·(driving force);
- **LHHW:** k·(driving force)/θⁿ. A reaction may override its site's term
  with terms of its own, defined as its own θ (`DEN_2`), sharing the site's
  constants where their texts are the same and refused where a name would
  mean two;
- **Eley–Rideal:** k·K_A·a_A·a_B·(1 − Q/K)/θ, K_A the adsorbed species' own
  term on the site, K_A·a_A¹;
- **expression:** any expression in the context.

Every rate must be a mol/(kg·s), per mass of catalyst, or a mol/(m³·s),
per volume, as its law says.

**Evaluating.** A point holds a state's inputs, made once and filled at
each state ([kinetics/evaluate.v](../kinetics/evaluate.v)):

    pᵢ = yᵢP,   fᵢ = φᵢyᵢP,   cᵢ = yᵢP/(ZRT),   yᵢ,   Keq_r = exp(ln K_r(T)),

K only for the reactions some rate uses. The rates are the tape's outputs,
into the caller's buffer. Their derivatives with respect to T, P and each
yⱼ, with φ and Z held fixed (exact for the ideal gas; a reactor that knows
their derivatives adds their terms), follow from the tape's gradient g by
the chain rule:

    dr/dT  = g_T − Σᵢ g_cᵢ cᵢ/T + Σ_r g_Keq_r Keq_r ΔH°_r/(RT²)
    dr/dP  = g_P + Σᵢ (g_pᵢ + g_fᵢ φᵢ) yᵢ + g_cᵢ yᵢ/(ZRT)
    dr/dyⱼ = (g_pⱼ + g_fⱼ φⱼ) P + g_cⱼ P/(ZRT) + g_yⱼ

### Checks and the factor table

M9 step 6 ([kinetics/check.v](../kinetics/check.v),
[kinetics/explain.v](../kinetics/explain.v)), `reflux kinetics check` and
`explain` (D12). Reading a file refuses what's wrong: an equation that
doesn't balance, a rate or a constant in the wrong units, a name that isn't
there. The checks are on what's legal but worth a look, each a note or a
warning:

| Finding | When |
|---|---|
| a reaction that combines others (note) | `dependence`: fine for kinetics, not for an equilibrium over reactions |
| ΔG°f against ΔH°f and S° (warning) | a reaction's Σν ΔG°f more than 1 kJ/mol from ΔH°(T°) − T°ΔS°(T°) |
| a fitted K against the Gibbs energies' | ln K_fit − (ln K° + Δn ln a°), a° the basis's standard state (p°, p°/RT or 1), at the state's T if there is one, else from 300 to 1500 K; a warning past a factor of 2 |
| a negative activation energy (warning) | Eₐ < 0 in an Arrhenius constant: an apparent one can be, from lumped constants |
| an adsorption enthalpy above zero (warning) | ΔH > 0 in a site term's van 't Hoff constant: K rising with T |
| a constant at or below zero (warning) | A, or a fixed value, ≤ 0 |
| a reaction's own adsorption term (note) | an override, with the terms it lacks or adds |
| a look-alike constant (warning) | a number within 1% of a built-in constant of its dimension (expr's lint), and, at a state, how far it moves each rate: the trees evaluated with the number swapped for the constant, r_typed/r_constant − 1 |

**The factor table** evaluates each rate's parts at a state, as the rate
law recorded them when it was built: its constant; K in its activities, Q
and Q/K for a reversible rate; its site's term and each of that term's
terms. Then every definition the rate uses, directly or through others,
with its value, and the rate. Values are in SI, a part's unit written over
the rate's bases (mol/(kg·s), mol/(m³·s), or none) times Pa or mol/m³ to a
power, mol·Pa^0.5/(kg·s) rather than its base units. Numbers are written to
five significant figures from their digits rounded to a whole number by
`math.round`, so the text is the same whichever compiler built reflux.

## Chemical equilibrium

M9 step 7 ([equilibrium/rand.v](../equilibrium/rand.v),
[equilibrium/equilibrium.v](../equilibrium/equilibrium.v)), D9 and D10.

**What's minimised,** at T and P, in one phase (D9):

    G/RT = Σⱼ nⱼ μⱼ/RT,   μⱼ/RT = gⱼ + ln(P/p°) + ln(nⱼ/n) + ln φⱼ,
    gⱼ = μ°ⱼ/RT = (H°ⱼ − TS°ⱼ)/RT,

with H° and S° along Cp as kinetics' standard properties take them, R
exact and p° = 1 bar, subject to C·n = b, b = C·n⁰. A Gibbs minimisation's
C is the element matrix, a row per element of the feed's species, with a
row of its own for each species held inert or without a formula. Before
anything else, what the feed can't make goes: a row whose coefficients
share a sign, with nothing of it in the feed (an element the feed hasn't),
fixes its species at zero, and they and it are dropped. What's left is
reduced to independent rows.

**The RAND method in CEA's form** (Gordon and McBride, NASA RP-1311, 1994):
Newton on the constraints' potentials π and on ln n, the total amount a
variable of its own, every species in log space, so that a vanishing one is
a large negative ln nⱼ and never a negative nⱼ. Each step solves

    Σᵢ (Σⱼ cₖⱼcᵢⱼnⱼ) πᵢ + (Σⱼ cₖⱼnⱼ) Δln n = bₖ − Σⱼ cₖⱼnⱼ + Σⱼ cₖⱼnⱼ μⱼ/RT
    Σᵢ (Σⱼ cᵢⱼnⱼ) πᵢ + (Σⱼ nⱼ − n) Δln n = n − Σⱼ nⱼ + Σⱼ nⱼ μⱼ/RT

then Δln nⱼ = −μⱼ/RT + Σᵢ cᵢⱼπᵢ + Δln n, damped by CEA's control factor

    λ = min(1, 2/max(5|Δln n|, |Δln nⱼ| over majors), λ₂),
    λ₂ = min over traces rising of |(−ln(nⱼ/n) − 9.2103404)/(Δln nⱼ − Δln n)|,

a trace being below 1e-8 of the total, so no trace rises past 1e-4 of it in
one step. Every species starts at an equal share of the feed's total. A
full step that moves no species' ln n by more than 1e-10 has converged,
Newton's error after it being of the order of its square; 500 steps
haven't. A step that small is never damped.

**A cubic** (D9): φ at the answer's amounts, from the model's stable root,
held fixed inside RAND and updated outside it until ln φ moves by no more
than 1e-12, at most 100 times. The answer's phase is then tested by
Michelsen's stability test, and one that would split is refused: a second
phase is M6's.

**Over chosen reactions** (D10, as amended in step 7), the same
minimisation with other constraints and potentials. C is what the
reactions conserve, a basis of the null space of their coefficients ν (the
vectors c with c·νᵣ = 0 for every reaction), by reduction to row echelon
form. The potentials are stand-ins,

    g = −νᵀ(ννᵀ)⁻¹ ln K,

the least-norm g with Σⱼ νᵣⱼ gⱼ = −ln Kᵣ for every reaction. At a
stationary point μⱼ/RT = Σᵢ cᵢⱼπᵢ, so each reaction's Σⱼ νᵣⱼ μⱼ/RT = 0:
ln Q = ln K, Q on fugacities at 1 bar. A species no reaction touches is a
conserved quantity of its own, held at its feed. The null space's basis
mixes the elements, so the element rows, which every balanced equation
conserves, are looked at too for what the feed can't make. The reactions
must be independent (`kinetics.dependence`).

**Shown right before it's returned:** every conserved row holds the feed's
to 1e-12 of its largest term; a Gibbs minimisation's every species present
has μⱼ/RT = Σᵢ cᵢⱼπᵢ to 1e-9, G at a minimum over what the elements allow;
an equilibrium over reactions has each ln Q its ln K to 1e-9; and a
cubic's phase is stable. A species held alone by a row, an inert one, is
set to the feed's amount exactly. The answers match mpmath's at 50 digits
by another method (Newton on ln n, `tools/gen_m9_equilibrium.py`) to
6.0e-13 in ln y, combustion's traces at 1e-29 included.

### The consistency check

Kinetics.md's check 4 ([equilibrium/consistency.v](../equilibrium/consistency.v)):
from a state's T, P and mole fractions, taken as amounts, the equilibrium
over the kinetics' independent reactions, K from the Gibbs energies. Then
each reversible rate, a structured law with a K or an expression, is
followed along its own reaction from there,

    n(ξ) = n_eq + νᵣ ξ,

to where it's zero: probing outwards both ways, 2⁻⁵⁰ of the way to where a
species runs out and doubling, then halving what's left, the way it runs
first, until its sign changes, then bisecting to rounding. Both ways,
because a rate with its sign reversed is zero at equilibrium but runs away
from it. Its reaction's ΔG/RT = ln(Q/K) there, Q on fugacities, is where it
stops:

| Rate | Stops at ΔG/RT |
|---|---|
| consistent with the Gibbs energies | 0, within 1e-10 (2.7e-15 at worst on Xu and Froment's) |
| a fitted K, in the ideal gas | ln K_fit − (ln K° + Δn ln a°), check_fitted's distance |
| on partial pressures in a gas that isn't ideal | Σν ln φ, where Q on pressures is K |
| on concentrations in such a gas | Σν ln φ + Δn ln Z |
| on mole fractions without (p°/P)^Δn | Δn ln(P/p°) |

A consistent rate is zero where its own Q is K, whatever the equilibrium's
accuracy, so the tolerance is rounding's. A stop past 1e-10 is a note, past
a factor of 2 in Q/K a warning, as a fitted K's distance is. At the state,
a rate that runs against its reaction's ΔG (r·ΔG > 0, |ΔG/RT| past 1e-10),
and not because it stops beyond equilibrium as a fitted K can, is warned
of. A rate never zero along its reaction is a warning for a structured law
and a note for an expression, irreversible as written. A reaction whose
species the state can't make isn't judged, with a note saying so.
`reflux kinetics check` runs it in the ideal gas when given a state.

## The integrators

M9 step 8 ([numeric/ode.v](../numeric/ode.v),
[numeric/dopri5.v](../numeric/dopri5.v),
[numeric/rodas4.v](../numeric/rodas4.v)), D8: `dopri5` for what isn't stiff,
`rodas4` for what is. Both follow Hairer's codes, DOPRI5 and RODAS, in
their steps, error norm and step control, and take the problem as an
interface: `Ode`, y' = f(t, y), and `StiffOde`, with ∂f/∂y and ∂f/∂t.

**The error norm** is Hairer's: a step is accepted when

    ‖e‖ = √(Σᵢ (eᵢ/sᵢ)²/n) ≤ 1,   sᵢ = atolᵢ + rtol·max(|y₀ᵢ|, |y₁ᵢ|),

e the error estimate, y₀ and y₁ the step's ends. rtol must be above 10
units of roundoff and every atol above zero, as Hairer's codes ask.

**The first step** is Hairer's HINIT: with d₀ = ‖y₀‖ and d₁ = ‖f(t₀, y₀)‖
in the norm above, h₀ = 0.01·d₀/d₁ (1e-6 where either is below 1e-5), one
Euler step of h₀ estimates d₂ = ‖f(t₀ + h₀, y₀ + h₀f₀) − f₀‖/h₀, and

    h = min(100·h₀, (0.01/max(d₁, d₂))^(1/p), hmax),

p the method's order, 5 or 4.

**Dormand and Prince's 5(4) pair** (Hairer's CDOPRI coefficients, SciPy's
RK45's too): seven stages, the last f at the new point and the next step's
first (FSAL). The fifth-order solution goes on; the error estimate is its
difference from the fourth-order one, e = h·Σ eᵢkᵢ. The step control is
Hairer's PI controller, with β = 0.04:

    h_new = h / max(1/10, min(5, ‖e‖^(0.2 − 0.75β) / ‖e_prev‖^β / 0.9)),

‖e_prev‖ the last accepted step's, at least 1e-4. A rejected step shrinks
by min(5, ‖e‖^(0.2−0.75β)/0.9) and the step after a rejection can't grow.
**Stiffness:** every 1000 accepted steps, and on each step while it
suspects stiffness, |hλ| ≈ |h|·‖k₇ − k₆‖₂/‖y₁ − y₆‖₂, the last two stages
at the same t. Above 3.25 for 15 steps in a row is stiff, and refused for
rodas4; 6 steps below clear the suspicion.

**RODAS4** (Hairer and Wanner II, §IV.7; the coefficients ROCOE's, METH =
1): six stages of

    (I/(hγ) − J)·kᵢ = f(t + cᵢh, y + Σⱼ<ᵢ aᵢⱼkⱼ) + Σⱼ<ᵢ (cᵢⱼ/h)kⱼ + h·dᵢ·∂f/∂t,

γ = 1/4, J and ∂f/∂t at the step's start, one LU of I/(hγ) − J for all six.
y₅ = y + Σ a₅ⱼkⱼ, y₆ = y₅ + k₅, y_new = y₆ + k₆: stiffly accurate, so
L-stable, with y₆ the embedded solution of order 3 and k₆ the error
estimate. The step control is Gustafsson's predictive controller, as
RODAS's:

    h_new = h / max(1/6, min(5, ‖e‖^(1/4)/0.9)),
    and after the first step, h / max(that, (h_prev/h)·(‖e‖²/‖e_prev‖)^(1/4)/0.9),

‖e_prev‖ at least 0.01. A singular I/(hγ) − J halves the step, five times
in all at most, as RODAS counts them.

**Failures:** an evaluation that fails, or isn't finite, rejects the step
and halves it, so a problem can refuse a point outside its domain (a
negative flow) and the integrator steps back. A step that falls to 10
units of roundoff in t, or more than max_steps, is an `OdeError`, with
the last accepted point.

**What the tests show** ([numeric/ode_test.v](../numeric/ode_test.v),
[numeric/stiff_test.v](../numeric/stiff_test.v)):
- Each solution's order, by one step at halving sizes: RODAS4's 4 and its
  embedded 3, Dormand–Prince's 5 and 4, within 0.35 (0.14 and 0.29 at
  worst).
- RODAS4's |R(z)| ≤ 1 over the left half-plane, and R(−1e12) = 8.8e-12.
- Closed forms: a linear system with eigenvalues −1, −1e3 and −1e6, and
  Prothero and Robinson's problem, within 0.16 tolerance units; Kepler's
  orbit and Arenstorf's, whose errors grow along the orbit (730 to 24,000
  units), each recorded.
- The Test Set for IVP Solvers' HIRES, Robertson, van der Pol, Oregonator
  and pollution problems within 0.13 tolerance units of the reference, the
  Oregonator's within 2.9, at rtol 1e-4 to 1e-10, the steps recorded.
- **E5 is beyond RODAS4 in double precision.** y₂ − y₃ − y₄ is conserved
  exactly, and any Rosenbrock method keeps a linear invariant in exact
  arithmetic. But a Rosenbrock stage solves for its whole increment from
  f's values, so its rounding is relative to f's terms, where an implicit
  method's Newton solves for small corrections. While E5's y₂ and y₃ are
  1e-10 their rates' terms are 1e-7 and cancel, and the 1e-22 lost then
  stays when they have fallen to 9e-23 at t = 1e13. RODAS4's y₂ and y₃ end
  51% from the reference under tcc and a factor of 20 under gcc, whose last
  bits differ, the drift 7.9e-23 and 1.9e-21; SciPy's Radau 3.4e-4; the
  reference, in 128-bit arithmetic, keeps the invariant to 1.3e-25. The
  test holds y₁'s and y₄'s decay and the drift below 1e-20. A PFR's flows
  don't fall a dozen orders below their rates' terms.

## The reactors

M9 step 9 ([reactors/](../reactors/pfr.v)), D10: five flowsheet units with
makers for the file (D1), each with one inlet and one outlet, its outlet a
flash at its state, in one phase.

**Balances on the elements.** A unit that reacts reports what it conserves
(`Outcome.conserved`): a row per element of the file's components, with a
row of its own for a component without a formula, which can't react. The
flowsheet's solver checks such a unit on those rows instead of each
component, to 1e-10 of its flow, and the whole flowsheet too where it has a
reactor, a tear's gap in one component allowed to move an element's
balance by the most atoms a component has. The energy needs nothing new: H
is referred to the elements (§ Ideal gas), so each stream's H holds its
heats of formation.

**Temperature.** Every reactor is isothermal (at its `T`, or its inlet's),
adiabatic, or at a given duty. An isothermal reactor's duty is F·H out
less F·H in. For an equilibrium, conversion or Gibbs reactor out of
isothermal, the outlet's T solves F·H(T) = F·H_in + Q by Brent, bracketed
outwards from the inlet's T by doubling steps from 50 K, between 0.5 and
10,000 K, the outlet recomputed at each T; the CSTR and the PFR carry T
themselves.

**The reactors:**

| Reactor | Its outlet |
|---|---|
| Gibbs | `equilibrium.gibbs` over every component (D9) |
| equilibrium | `equilibrium.over_reactions`, K from the Gibbs energies, its reactions independent |
| conversion | each step in order: extent ξ = X·n_key/(−ν_key), on the key as it enters the step; a reactant taken below zero beyond 1e-12 of the flow is refused, and within it is zero. X is a number or an expression in T, evaluated at the outlet's T (an adiabatic or a duty's found with the energy balance), which must give 0 to 1 there |
| CSTR | Fᵢ,in − Fᵢ + V·Sᵢ = 0, and where T moves F·H(T, P, F) = F·H_in + Q, by `numeric.newton` from the inlet; where that fails, by continuation in V, from 10⁻⁶ of it to the whole, each solve from the last, the step in ln V doubled on success and halved on failure (fast rates near equilibrium, as in reforming at 1000 K, make the balances too steep for Newton from the inlet) |
| PFR | dFᵢ/dV = Sᵢ, with T and P where they move, by `rodas4` (or `dopri5`) |

**A conversion reactor past equilibrium** (D10, as amended before step
10) is warned of, not refused, since its conversion is the user's to give:
for each step's reaction, at the outlet,

    ln(Q/K) = Σᵢ νᵢ (ln yᵢ + ln φᵢ + ln(P/1 bar)) − ln K(T) > 10⁻⁸,

K at 1 bar from the Gibbs energies (D7) and φ from the outlet's flash, the
tolerance being what the two's rounding allows. The warning gives Q/K, T
and P; a reaction whose species lack formation data is noted as unchecked.
The owner's Tutorial 2 at 320 °C holds 27.4% at 30 bar, where Q/K = 4.56.

**A reactor's source of each species,** from the kinetics' rates
([reactors/source.v](../reactors/source.v)):

    Sᵢ = Σⱼ sⱼ νᵢⱼ rⱼ(T, P, y),   mol/(m³·s),

sⱼ the catalyst loading (kg/m³) for a rate per mass of catalyst and 1 for
one per volume. φ, Z and the phase's volume come from the fluid's model at
each state. A flow that rounding takes below zero is read as zero in the
rates, as DWSIM reads it. Its derivatives by the flows follow from the
rates' by the chain rule through yₗ = Fₗ/ΣF:

    ∂Sᵢ/∂Fₖ = Σⱼ sⱼνᵢⱼ [∂rⱼ/∂yₖ − Σₗ yₗ ∂rⱼ/∂yₗ]/ΣF,

with φ and Z held fixed, as kinetics' `jacobian` holds them.

**The PFR** ([reactors/pfr.v](../reactors/pfr.v)): the state is the flows,
then T where it moves, then P under Ergun's drop; a given drop's P is a
line in V. Where T moves, F·H rises by the duty per volume q, and with
h̄ᵢ = H°ᵢ(T) − RT²(∂ln φᵢ/∂T) each species' partial molar enthalpy, and
(∂H/∂P)_T = −Cp·μ_JT,

    dT/dV = (q − Σᵢ h̄ᵢSᵢ + Cp·μ_JT·dP/dV)/Cp.

Ergun's drop, with u = Q_v/A the superficial velocity, A = V/L, ρ the gas's
density and μ its viscosity:

    dP/dV = −[150 μ(1 − ε)²u/(ε³d²) + 1.75 ρ(1 − ε)u²/(ε³d)]/A.

`rodas4` takes the flows' block of the Jacobian exactly, from the rates;
the rows of T and P, and a given drop's ∂f/∂V, are central differences of
the right-hand side, each column's step 6e-6 of its value (about ε^(1/3)).
The absolute tolerances are 1e-6·rtol of the inlet's flow for each flow,
rtol·T_in and rtol·P_in. At the outlet, an isothermal PFR is flashed at its
T and P; otherwise the outlet is a PH flash on the energy balance, and the
integrated T must agree with its T to 100·rtol of T (plus 1e-9 of it), or
the energy equation and the fluid disagree and the PFR fails. The profile
is the inlet and every accepted step.

**The gas viscosity** ([reactors/viscosity.v](../reactors/viscosity.v)),
for Ergun (D10, as amended in step 0): each component's from its file's
vapour-viscosity correlation (ChemSep's equation 102) at T, mixed by
Wilke's rule,

    μ = Σᵢ yᵢμᵢ / Σⱼ yⱼΦᵢⱼ,   Φᵢⱼ = [1 + (μᵢ/μⱼ)^½(Mⱼ/Mᵢ)^¼]² / [8(1 + Mᵢ/Mⱼ)]^½,

pressure ignored. It matches `chemicals`' Wilke and EQ102 to 2.2e-16.

Everything in it but the mole fractions depends on T alone: each μᵢ, and
each pair's 1 + (μᵢ/μⱼ)^½(Mⱼ/Mᵢ)^¼ and [8(1 + Mᵢ/Mⱼ)]^½. A PFR keeps those
(`WilkeAt`) until T moves (M8's D11), so that along an isothermal reactor
each state costs the sums alone, with no correlation, power or square
root. The sums are the rule's own arithmetic in its own order, so the
viscosity is the same number, to the bit, kept or not. A viscosity that
isn't Wilke's, as a test's, is asked whole at every state.

**What the tests show:** an isomerisation's exp(−kPV/F) and 1/(1 + kPV/F)
to 3.2e-11 and 1e-12, and a dehydrogenation's −X − 2 ln(1 − X) = kPV/F₀ to
1e-11; the Gibbs and equilibrium reactors on mpmath's equilibria to 1e-11
in ln y; Ergun's drop within 3.3e-4 of its ideal-gas closed form; a long
PFR on fugacities within 2.8e-14 of the Gibbs reactor's equilibrium; and
every reactor in a chain closing its elements to 1e-10.

## Streams and units

A material stream holds its component molar flows nᵢ (mol/s) and the flash
that gave its state (M4's D2), so it's always fully specified. Its molar
properties are its phases', weighted by their shares of its moles:

    h = Σₖ βₖ hₖ,    s = Σₖ βₖ sₖ,    v = Σₖ βₖ vₖ

A feed is given by (T, P), (P, h), (T, β) or (P, β), each by its own flash.
An empty stream, which a split fraction of 0 makes, has no flows but keeps
its composition in its flash, so it has a state and passes through a valve
or heater.

A unit takes its inlets and gives its outlets and any energy streams. A
rate Q above zero goes into the unit's material, so with F the total flow,
each unit balances as

    Σ_in nᵢ = Σ_out nᵢ,    Σ_in F h + ΣQ = Σ_out F h

- **Mixer.** The flows and enthalpy flows add, then a PH flash at the outlet
  pressure:

      nᵢ = Σⱼ nᵢⱼ,    h = Σⱼ Fⱼ hⱼ / F,    P = minⱼ Pⱼ

  P may be given, but not above the lowest inlet's. An empty inlet adds
  nothing and doesn't set P. The flash starts from the flow-weighted mean T.
- **Splitter.** Each outlet k has the inlet's state and a fraction of its
  flows, nᵢₖ = fₖ nᵢ with Σₖ fₖ = 1. Nothing is flashed.
- **Valve.** Isenthalpic: a PH flash at h_out = h_in and the outlet pressure,
  from the inlet's T. It can't raise the pressure. Throttling is
  irreversible, so s can only rise.
- **Heater.** With P_out = P_in − ΔP, the specification picks the flash:

  | Specification | Outlet |
  |---|---|
  | outlet T | PT at (T, P_out) |
  | ΔT | PT at (T_in + ΔT, P_out) |
  | outlet β | the vapour-fraction flash at (P_out, β) |
  | duty Q | PH at (P_out, h_in + Q/F), from the inlet's T |

  The duty it reports is always Q = F (h_out − h_in). With a duty given,
  that's the one the flash reached, within its energy tolerance.
- **Pump.** With v the inlet's molar volume,

      h_out = h_in + v (P_out − P_in) / η,    W = F (h_out − h_in)

  then a PH flash. At η = 1 this holds v at the inlet's, where dh = v dP at
  constant s has v follow the liquid, so it leaves
  T Δs ≈ ½ (v_in − v_out) ΔP. A pump takes no vapour. An empty inlet has
  none to refuse.
- **Compressor and expander.** A PS flash finds the isentropic outlet,
  h_s at (P_out, s_in). Then

      compressor:  h_out = h_in + (h_s − h_in) / η
      expander:    h_out = h_in + η (h_s − h_in)

  and a PH flash finds the outlet, which can be two-phase, from T_s. The work
  W = F (h_out − h_in) is below 0 for an expander.
- **Flash drum.** The inlets blend as in a mixer, then a flash at the drum's
  pressure: PT at a given T, or PH at h = (Σⱼ Fⱼ hⱼ + Q) / F. A duty of 0 is
  adiabatic. Each phase k leaves by its label's outlet (M6's D9 and D11):
  the vapour by the vapour port, and every liquid, liquid2 too, by the
  liquid port, a stream of both at their shares. Its flows are

      nᵢₖ = F βₖ xᵢₖ

  saturated at the flash's T and P. That needs no second flash: its result is
  the first's, with the phase's composition as z and the other phase at no
  moles, the form a bubble or dew point's answer has. A single phase leaves by
  its own outlet, and the other outlet is empty. A liquid–liquid split, D9's
  liquid and liquid2, leaves by the liquid port, and the drum's outcome warns
  that its liquid has two phases.
- **Decanter** (M6's D11). A flash drum with an outlet for each phase: the
  same blend and flash, and each phase k by its label, the vapour by vapour,
  liquid by light and liquid2 by heavy, each nᵢₖ = F βₖ xᵢₖ as above. A
  label the flash didn't give has an empty outlet, so a single liquid, D9's
  liquid, leaves by light.
- **Component separator.** A fraction fᵢ of each component goes to the first
  outlet and the rest to the second, each flashed at the given T and P. It's
  a black box, so Q is whatever closes the balance, and its pressure can rise.
  An empty outlet takes the inlet's composition, so it has a state.

## The flowsheet graph

A flowsheet's units are nodes, and the streams between them are edges (M4's
D5). Two units can be joined by parallel edges, and a stream can leave and
enter one unit. `Flowsheet.plan` turns that into blocks to solve in order.

**Partition.** Tarjan's algorithm (1972) finds the strongly connected
components, each after the ones it needs. A component with a cycle in it is
a *net*: its units have to be solved together. Everything else is solved
once, in order. A net with nothing coming into it, from a feed or another
unit, has no material in it at all, and is refused.

**Cycles.** Johnson's algorithm (1975) lists a net's elementary cycles, as
lists of edges, so that parallel streams are separate cycles.

**Tears.** Tearing a set of streams T breaks a cycle c if T meets it. With
b_c = |T ∩ c| the times c is torn, the choice is the best T by, in turn:

1. fewest tears, |T|, with every b_c ≥ 1: a minimum hitting set;
2. fewest second breaks, Σ_c (b_c − 1). A set with none is what Upadhye and
   Grens (1975) call nonredundant, and it's what converges best by direct
   substitution. It isn't always there at the minimum size: the sulfuric
   acid plant's 5 tears can't do better than 36;
3. earliest streams, by the unit they leave, so a plan is reproducible.

Among streams in exactly the same cycles, which all do equally well, it takes
one coming back to where material enters the block: the recycle, in other
words. A tear there can start empty, so the first pass is once-through with
real material, and no guess reaches a pump as a gas.

This order is the reverse of FOQUS's and Pyomo's, which keep each cycle's
breaks down first and count tears second. Every tear here costs a block of
iteration variables, so reflux takes the fewest.

The search is branch and bound, which is exact when it finishes:
- Streams in exactly the same cycles, as those in series are, form a class,
  and only one class member can be chosen. These are Upadhye and Grens's
  families, whose members behave alike under direct substitution.
- It branches on the shortest cycle not yet torn, one branch per class. A
  class tried in one branch is excluded from its siblings, so that no set is
  reached twice.
- Cycles that are untouched and share no class need a tear each, which
  bounds what's left.

Its budget is a count of cycles looked at, not a time, so plans are
reproducible. If it runs out, or a net has too many cycles to list (50,000),
the plan says so in its notes. In the second case a greedy cover takes over,
tearing the stream in the most cycles until none is left.

**Order.** Within a block, with its tears cut, the units go in topological
order, the earliest in the flowsheet first.

## Converging the tears

Each tear has the variables x = (n₁ … n_c, P, h), with h molar (M4's D6). A
pass runs a block's units once from x, and gives g(x), the tears' new values.
The solution is a fixed point, g(x) = x. A tear is rebuilt from x by a PH
flash, starting from its last T, unless x is the pass's own answer, which
needs no flash.

**Convergence** (D7). Every variable must move by at most tol = 1e-9 of its
scale in a pass:

    |g_i − x_i| ≤ tol · { F for a flow, P for P, RT for h }

with F the tear's total flow, or for an empty tear, what enters the block. h
is against RT because h itself can be near zero. Each pass's largest moves,
tear by tear, are the history.

**The first values.** A tear starts from its estimate. Without one, what
enters the block is mixed as a mixer would, and a tear whose unit has another
inlet starts empty, in that mixture's state, so the first pass is once-through.
Any other starts as the mixture itself.

**Direct substitution:** x ← g(x). On a loop of gain λ the error falls by λ a
pass, so it takes about ln(tol)/ln(λ) passes.

**Bounded Wegstein.** Variable by variable, from the last two passes,

    s = (g(x) − g(x₀)) / (x − x₀),    q = s / (s − 1),    x ← q·x + (1 − q)·g(x)

with q held to [−5, 0], so that it only ever accelerates. For a linear loop of
gain s that's exact. Held at q = −5, it leaves −5 + 6s of the error. It's
applied every other pass, from the third: the first pass starts from a guess,
so a slope through it says nothing about the loop, and taking one nearly
doubled the passes on a nonlinear loop.

**Broyden.** One pass, for the same reason, and then Broyden's "good" method
on all a block's tear variables together, in unknowns scaled by the flow into
the block, P and RT, with residuals (g − x) over the scales above. With
J₀ = −I its first step is a direct one. With s the step and y the change in
the residual, the inverse Jacobian updates as

    H ← H + (s − H y) (sᵀ H) / (sᵀ H y)

**A loop that's kinder than it looks.** Where a recycle is a drum's liquid
brought back to the drum's own T and P, as in the reference flowsheet, adding
it changes nothing about the equilibrium. The liquid just grows, the loop is
affine with a gain of the split fraction, and Wegstein's one secant lands on the
answer in 4 passes. Put the cooler on a duty, so that T moves with the recycle,
and it takes 14.

**Balances** (D8). A unit must close to 1e-10 of its flow and of F·RT on the
inlets it was calculated from. The whole can't close better than the tears
allow, and it must close to exactly that: the gap between the last pass's start
and end, in each tear's flows and enthalpy flow, plus rounding. A solution that
doesn't is an error, never a result.

## Varying a number

A study, a spec or an optimiser varies a number the file gives a feed or a
unit, by its key (M8's D4): `unit.cooler.outlet_T`, `feed.f.T`,
`feed.f.molar_flows.methane`, `unit.split.fractions[0]`. Nothing is set on
a unit: the value replaces the file's in its tree, and the feed or unit is
made again by its own maker, which checks it as it checks the file.

**Readings.** A `Section` notes every number its maker reads, with its key,
its SI value and the dimension it was read as, M4's named one where there
is one (a temperature, or a temperature difference). A key can be varied
only if its maker read it as a number: so `27 degC` is 300.15 K, as the
maker had it, and `type`, `in`, or a conversion that's an expression in T
are refused, with the keys that can be.

**Writing.** A new value goes into the tree as a pure number, or as
`value.str()` with the dimension's SI unit (K, Pa, mol/s, W, J/mol, or base
units such as kg/m^3), each of scale 1. V's shortest decimal reads back to
the same double (200,000 random doubles, measured), so the maker reads the
value given to the bit, and `with` checks that it did. Only the tables and
lists on the key's path are copied, and the rest of the file is shared.

**What it costs** (gcc `-prod`, `tools/remake_bench.v`): 8 to 80 µs for a
unit, against 0.3 µs to 8 ms for its own calculation, and 0.3 to 0.6 ms for
a feed, which is a flash.

## Specs

A spec varies an input within its bounds until an equation holds (M8's D5):

    r = (lhs − rhs) / s,    s = max(|lhs|, |rhs|, s₀),    |r| ≤ tol

where s₀ is the block's own scale for the equation's dimension: the flow
into the block for a flow, F·R·T_max for a power, R·T_max for a molar
enthalpy, 1 for a pure number, and 0 otherwise, a temperature or pressure
being far from zero. So a target of 0 W is met to tol of F·RT.

**Planning.** A spec is a link from each unit that makes what its equation
reads to the unit it varies, or the units the feed it varies enters, as if
a stream joined them. Tarjan's components then put both in one block, and
the link is torn already: the spec's value is the block's unknown, and the
tears are the fewest streams that break the cycles left. A spec none of
whose links closes a loop is refused, since nothing it reads depends on
what it varies. A link that doesn't close one still orders the blocks, its
reading made before its block.

**One spec, no tears: Brent** (1973) within the bounds, each evaluation the
spec's unit remade and the block calculated once. The bracket is the bounds
if r changes sign between them; else nine values between, evenly spaced; and
a target no pair brackets is refused with what the values reached ("gas.F
goes from 41.99 to 76.65 mol/s, never 200 mol/s"). A value Brent asks for
again is the sampled one, not calculated again.

**Several specs, or tears: Broyden** (1965) on the tears' scaled variables
and the specs' values together, each evaluation one pass, its residuals the
tears' changes and each spec's r. The specs' first Jacobian J is a pass's,
each value stepped by 1e-4 of its bounds' width from the tears' start, and
a spec's value is

    v = v₀ + m·u,    m = −J⁻¹,

so that the first step in u, with Broyden's H₀ = −I, is Newton's on the
specs, as it's direct substitution on the tears. A value is held within
its bounds. Held at a bound three evaluations running, the target is likely
out of reach: Broyden is stopped, and the spec is tried again nested, the
specs' values outside (Brent for one, Broyden for more) around the tears
converged from where the last evaluation left them. Specs whose J is
singular can't be met together from there, and are refused.

**The answer** is a solution like any other, every balance checked (M4's
D8), and each spec checked again from the solution's streams and units:
|r| ≤ tol, and the value within its bounds.

## Warm starts

A point of a study starts from its neighbour's solution (M8's D6), and its
answer is checked as a cold one is. What it carries was measured (M8 step
5), and is what pays:

- **Units** that can, through `Warm`: the CSTR's Newton from its last
  outlet, falling back to the inlet and continuation; an energy balance's
  T bracketed out from its last T by steps of 5 K that double, where cold
  it starts from the inlet's by 50 K; the PFR's first step.
- **Specs,** from their last values. One spec alone is met or bracketed
  near there first, by secant steps out from it, each at most four times
  the last, the first as large as the last two points' change, before the
  cold search over its bounds.
- **Tears,** with direct substitution and in a block with specs only.
  Wegstein's and Broyden's own starts beat a neighbour's tears.

With two points behind it, a tear's flows, P and h and a spec's value
start from the secant through them,

    x = b + r·(b − a),    r = (this step)/(the last step),

where both steps are along one input in one direction; a flow is held at
zero, and a tear is a PH flash from b's T. Methods' own state, Broyden's
inverse Jacobian or Wegstein's slopes, isn't carried: it took more passes.

## Studies

A study (M8's D7, D9) solves a grid of one or two inputs' values, each
point remade from the file (D4), and evaluates its outputs there (D3, D9).

**Order.** Along the first input's values; with two, along each row of the
second's, alternating direction, so that each point's neighbour is the last
solved. Points go in chunks of 8, or as the file says: each chunk's first
point is cold and the rest warm, so the answers don't depend on how chunks
are shared out (D8's workers, step 7).

**Each point** is solved, solved with an output that has no value (its
reason kept, D9), or failed (its reason kept). The study goes on past a
failure, and the next point starts cold. Evenly spaced values are
from + (to − from)·i/(n − 1), the last set to `to` exactly.

**Workers** (D8): threads, made with `spawn`, share the loaded file, which
nothing changes, and take whole chunks from a queue as they come free. A
chunk starts cold, so its points are the same bits whichever thread solves
them, and in whichever order the chunks are taken. `reflux study` runs as
many as the processors, `NUMBER_OF_PROCESSORS` on Windows and
`sysconf(_SC_NPROCESSORS_ONLN)` elsewhere, or `--workers`.

## Sensitivities

A study's point gives the derivative of each output f with respect to each
input v it's asked for (M8's D10), by central differences of whole solves:
nothing about a unit's equations is needed, so every unit has them, specs
met again at each step included.

    h = 10⁻⁴·|v|,    D(s) = (f(v + s) − f(v − s)) / (2s),
    f′ ≈ D(h/2) + (D(h/2) − D(h))/3                      (Richardson)

The step balances the difference's truncation, h²f‴/6, against the solves'
noise over h, ε·f/h, at ε^(1/3) for ε = 10⁻¹² (the solves' tolerance, from
the point's own solution). Each difference is over the steps the doubles
actually are, (v + s) − (v − s), not 2s. The error's estimate has two
parts:

    e = |D(h/2) − D(h)| + |Δ⁴|/h,
    Δ⁴ = f(v − h) − 4f(v − h/2) + 6f(v) − 4f(v + h/2) + f(v + h)

The first is three times D(h/2)'s truncation error, which the extrapolation
removes to O(h⁴). The second measures the solves' noise: Δ⁴'s truncation is
(h/2)⁴f⁗, so what's left is noise, about √70 times one solve's, and the
extrapolation carries about 1.9 times one solve's over h. A spec met at
10⁻¹² of its scale, or a unit's own tolerance (a PFR's rtol, which isn't
tightened), shows there. An input at 0 has no scale for a step, and is
refused with that reason; an output with no value at the point has no
derivative there, which says so. Exact sensitivities, by the converged Jacobian,
are M13's.

**What the estimate can't see:** a solve started from the point's own
solution stops where its tears move by the tolerance, which leaves an error
in proportion to how far it started, and so to the step. That's a bias
linear in the step: it cancels from both parts of the estimate, and only
the tolerance bounds it, to about ε/10⁻⁴ of the derivative, 10⁻⁸ at
ε = 10⁻¹². At a study's 10⁻⁹ it would be 10⁻⁵, which is why the tolerance
is tightened (M8 step 8 found it, by a mistake planted there).

## The optimiser

`numeric.sqp` (M8's D12) minimises f(x) subject to equalities c(x) = 0,
inequalities c(x) ≥ 0 and bounds l ≤ x ≤ u, by Han and Powell's sequential
quadratic programming as Nocedal and Wright give it (*Numerical
Optimization*, ch. 18). The problem is an interface, `Constrained`: its values,
then its gradients at the same x, so that a flowsheet differentiates from
the solve it just made.

**The QP.** At each iterate the step d solves

    minimise ½·dᵀBd + ∇fᵀd
    subject to c + ∇c·d = 0 (equalities), ≥ 0 (inequalities),
               l − x ≤ d ≤ u − x

and its multipliers λ̂ (the constraints'), μ̂ (the bounds') are the next
estimates. `numeric.qp` solves it by Goldfarb and Idnani's dual active-set
method (Math. Programming 27, 1983): from the unconstrained minimum −B⁻¹∇f,
it adds each equality not met and then the most violated inequality, by its
distance s/‖row‖. Adding row p moves d along z and the active multipliers
along −r, where

    B·z + N·r = nₚ,   Nᵀ·z = 0,

N the active rows. The step is the smaller of t₂ = −sₚ/(zᵀnₚ), which meets
p, and t₁, the least uⱼ/rⱼ over active inequalities with rⱼ > 0, which takes
one to zero and drops it. z = 0 with no t₁ proves p can't be met: the QP has
no point. The problems are small, so z and r come from the whole system,
factorised at each step.
z is taken as zero where zᵀnₚ ≤ 10⁻¹²·nₚᵀB⁻¹nₚ, its size with nothing
active.

**A QP with no point.** A linearisation can have no point though the
problem has; then the QP is relaxed, as Powell's and Schittkowski's codes
do: a variable δ ∈ [0, 1] scales each equality and violated inequality,

    ∇cᵢ·d + (1 − δ)·cᵢ = 0 or ≥ 0,   at a cost ½·ρδ²,  ρ = 10⁴·max(1, maxᵢ Bᵢᵢ),

and δ = 1, d = 0 always meets every row. A δ that comes out 1 means no step
lowers the constraints' violation, and the problem is refused: it may have
no feasible point near x. So is one whose relaxed QP fails too, which only
a badly conditioned one does.

**The answer.** x is returned once, with λ̂ and μ̂ from the QP at x:

    max |cᵢ| (equalities), max(−cᵢ) (inequalities)   ≤ tol
    ‖∇f − Σλᵢ∇cᵢ − μ_l + μ_u‖∞                         ≤ tol·max(1, ‖∇f‖∞)
    λᵢ (inequalities), μ_l, μ_u                       ≥ −tol·max(1, ‖∇f‖∞)
    |λᵢcᵢ|, μ_l·(x − l), μ_u·(u − x)                    ≤ tol·max(1, ‖∇f‖∞), the bounds' × max(1, |x|)

tol is absolute on c, so the caller scales its constraints; 10⁻⁸ by default.

**The step.** α·d is accepted by backtracking on the ℓ1 merit function

    φ(x) = f(x) + Σ νᵢ·vᵢ(x),   vᵢ = |cᵢ| (equalities), max(0, −cᵢ) (inequalities),
    φ(x + αd) ≤ φ(x) + 10⁻⁴·α·Dφ,   Dφ = ∇fᵀd − (1 − δ)·Σ νᵢvᵢ,

with Powell's weights νᵢ = max(|λ̂ᵢ|, (νᵢ + |λ̂ᵢ|)/2), at least |λ̂ᵢ|, which
makes Dφ ≤ −dᵀBd < 0. Each shorter α is the quadratic's minimum through
φ(0), Dφ and φ(α), kept within [α/10, α/2]; a step where the problem can't
be evaluated is halved. Every iterate is clamped into the bounds, which the
QP's d already keeps it within, so the problem is never evaluated outside
them.

**B.** The Lagrangian's Hessian, from the identity, by BFGS on
s = x₊ − x and y = ∇L(x₊, λ̂) − ∇L(x, λ̂), the bounds' part linear and so
cancelling, with Powell's damping:

    θ = 1 if sᵀy ≥ 0.2·sᵀBs, else 0.8·sᵀBs/(sᵀBs − sᵀy),
    r = θy + (1 − θ)Bs,   B ← B − BssᵀB/(sᵀBs) + rrᵀ/(sᵀr),

which keeps B positive definite. The first update scales the identity by
yᵀy/sᵀy (Shanno and Phua) where sᵀy > 0.

## Optimising a flowsheet

An `[optimize.NAME]` (M8's D13) is `numeric.sqp` on the flowsheet's
inputs, every evaluation a converged solve at 10⁻¹², specs and all,
started from the last evaluation's solution (D6); its first evaluation, at
the start, is cold. The problem the SQP sees is scaled:

    u = (x − lo)/(hi − lo) ∈ [0, 1],
    f = ±F/s_F,   s_F = maxₖ |∂F/∂uₖ| at the start (|F| there, or 1, if that's 0),
    c = (lhs − rhs)/sᵢ for = and ≥,  (rhs − lhs)/sᵢ for ≤,
    sᵢ = max(|lhs|, |rhs|) at the start, or 1 in SI if both are 0,

with − to maximise. The objective's gradient so starts at about 1, and
`optimize_tol`, 10⁻⁶, is the SQP's tol: well above the gradients' noise,
about 10⁻⁸ of their scale. The equalities go first, as the SQP takes them.

**Gradients** are the sensitivities' differences (above), of the objective
and of each constraint's two sides, with a step of 10⁻⁴ of the largest
value the input's bounds allow, max(|x|, |lo|, |hi|), so that an input at
0 has one, and at most a quarter of its range. Where a bound is within the
step, they're one-sided, at x + σ·h·k/4, k = 0 to 4, σ = ±1 into the range:

    D(s) = (f(x + σs) − f(x))/(σs),
    R₁ = 2·D(h/2) − D(h) = f′ − h²f‴/12,   R₂ = 2·D(h/4) − D(h/2) = f′ − h²f‴/48,
    f′ ≈ R₂,   e = 4/3·|R₁ − R₂| + |Δ⁴|/h,

the first part R₁'s truncation, four times R₂'s, since 4/3·|R₁ − R₂| alone
would be R₂'s estimate with no margin. R₂ carries about 10ε/h of the
solves' noise, twice R₁'s.

`--workers` takes the inputs at once, each from the point's solution, so
the bits are the same.

**The multipliers,** unscaled, are how fast the optimum's objective
improves (falls when minimised, rises when maximised) as each is relaxed,
per unit of its own:

    an inequality:   λᵢ·s_F/sᵢ   (≥ 0)
    an equality:    −λᵢ·s_F/sᵢ   (as its right side grows, either sign)
    a bound:         μ·s_F/(hi − lo)

so that, by the envelope theorem, a bound's is ∂F/∂x at the optimum, and
a constraint's ∂F/∂x over its sides' ∂/∂x there, for one input.

## Least squares

`numeric.least_squares` (M8's D14) minimises ½‖r(θ)‖² by Levenberg and
Marquardt's method as Moré gives it (Lecture Notes in Mathematics 630,
1978), over the interface `Residuals`: r and its Jacobian J at θ.

**The step** minimises the linear model with a damping term,

    ‖J·δ + r‖² + μ·‖D·δ‖²,   by Householder QR on [J; √μ·D]·δ = [−r; 0],

so JᵀJ is never formed and its squared conditioning never paid. Dₖ is the
largest norm of column k of J seen so far, which makes the steps
independent of the parameters' scales; μ starts at 10⁻³.

**Taking it.** With ρ the residuals' fall over the model's,
‖r‖² − ‖r + J·δ‖², a step is taken where ρ > 0, and then
μ ← μ·max(1/3, 1 − (2ρ − 1)³); otherwise μ ← μ·ν, ν doubling each refusal
and back to 2 at each step taken (Nielsen). Where the model predicts a fall
below 10⁻¹⁰·‖r‖², Σr²'s rounding, tens of ε of it, can't resolve it and ρ
is noise: there a step is taken if the gradient test falls. Without that, a
fit near its minimum stalls, every step refused, before the test can be met
(M8 step 11's coverage test found it, in 1 data set of 400; the gcc check
found a window of 8ε·‖r‖² on the actual fall too narrow).

**The answer.** θ is returned once the gradient test holds,

    maxₖ |Jₖᵀr| / (‖Jₖ‖·‖r‖) ≤ gtol = 10⁻¹⁰,

the largest cosine between r and a column of J, zero at a minimum whatever
the scales; or r is zero. (JᵀJ)⁻¹ comes from the QR of J with its columns
scaled to unit length, J·S⁻¹ = Q·R:

    (JᵀJ)⁻¹ = S⁻¹·R⁻¹·R⁻ᵀ·S⁻¹,

and a diagonal entry of R below 10⁻¹⁰ of the largest is a parameter the data
don't determine, refused by its index.

**Student's t,** for the intervals, is the root of its CDF by bisection to
the doubles' resolution,

    P(T ≤ t) = 1 − ½·I_x(ν/2, ½),   x = ν/(ν + t²),   1 − x = t²/(ν + t²),

with the regularized incomplete beta I by its continued fraction (Lentz's
method, Numerical Recipes §6.4) and its front factor's
ln(Γ(a + b)/(Γ(a)Γ(b))) by Stirling's series with log1p where a is large:
the difference of two numbers near 6·10⁶ at ν = 10⁶ otherwise. It's good to
1.2·10⁻¹⁴ of t at ν ≤ 1000, and 3.2·10⁻¹² at 10⁶.

## Fitting

A `[fit.NAME]` (M8's D14, D15) is least squares on weighted residuals,
measurement less model over its error σ (σ = w·|measured| for a relative
weight w), and the file's form of each parameter at the end:

    s² = Σrᵢ²/(n − p),   C = s²·M·(JᵀJ)⁻¹·Mᵀ,   seₖ = √Cₖₖ,
    θₖ ± t₀.₉₇₅(n − p)·seₖ,   ρₖₗ = Cₖₗ/(seₖ·seₗ),

M = ∂(the file's form)/∂(θ as fitted), the delta method, which is exact for
the linearised covariance (Bates and Watts, *Nonlinear Regression Analysis
and Its Applications*, 1988).

**Rate constants.** The kinetics are built again with the fitted numbers as
variables in the rates' context, after T and P, so the tape's reverse sweep
gives ∂r/∂θ exactly. Where an Arrhenius constant's A and Eₐ are both
fitted, and it isn't centred already, it's fitted centred on the data's
mean T, T̄:

    k = k_ref·exp(−Eₐ/R·(1/T − 1/T̄)),   A = k_ref·exp(Eₐ/(R·T̄)),
    ∂A/∂k_ref = exp(Eₐ/(R·T̄)),   ∂A/∂Eₐ = A/(R·T̄),

which takes most of A's and Eₐ's correlation out of the solve: in the file's
form, the ethylene fit's k1.A and k1.Ea correlate at 0.996. A van 't Hoff
constant's A and ΔH are centred alike. Each row's state is the ideal gas, pᵢ = yᵢP, fᵢ = pᵢ,
cᵢ = pᵢ/RT, as rate data from a differential reactor are.

**kᵢⱼ** come from bubble points: at each row's T and liquid x, the model's
bubble pressure and vapour composition (`flash.bubble_p`) less the
measured P and y, each over its error. J is by central differences,
kᵢⱼ ± 10⁻⁶, good to about 10⁻⁷: kᵢⱼ is of order 0.01 and a bubble point's
rounding about 10⁻¹².

**An activity model's pairs** (M6's D16) come from bubble points the same
way, on the γ–φ fluid with the fitted pairs in place of its own. A pair's
energy Aᵢⱼ, J/mol, is NRTL's τᵢⱼ = Aᵢⱼ/RT or UNIQUAC's τᵢⱼ = exp(−Aᵢⱼ/RT), R
the SI's exact value, and NRTL's α a number; each evaluation builds the
model again from those expressions. Its start is the model's own pair at
T₀ = 298.15 K,

    A = τ(T₀)·R·T₀ (NRTL),   A = −R·T₀·ln τ(T₀) (UNIQUAC),   α = α(T₀),

which a bank's pair, A₁₂ cal/mol over DECHEMA's R, gives at any T. J is by
central differences of an energy ± 10⁻⁴·max(|A|, 2500 J/mol) and α ± 10⁻⁴:
the bubble points' rounding, 10⁻¹², over the step is then about 10⁻⁸ of J,
and the difference's truncation about as much, which moves the answer by
as little of its standard error and stops no convergence where the data
determine the parameters. With steps of 10⁻⁶ the rounding kept NRTL's α
and both energies from the gradient test. Where they don't (α with both
energies on one isobaric binary), whether the test is met or Levenberg–
Marquardt stops at the rounding, saying so, is the platform's. The pairs fitted are written back as
`[[model.pairs]]`, each fitted energy as `"A J/mol / (R*T)"` (UNIQUAC's in
an `exp`).

**The point test** (M6's criterion 8) checks a binary's VLE data for
Gibbs–Duhem consistency, as Fredenslund, Gmehling and Rasmussen give Van
Ness, Byer and Gibbs' test (*Vapor–Liquid Equilibria Using UNIFAC*, 1977).
Gᴱ is a Legendre series, fitted to the bubble pressures alone by Barker's
method, the vapour the ideal gas:

    g = Gᴱ/RT = x₁x₂ Σₖ aₖ Lₖ(z),   z = x₁ − x₂,   k = 0 … 3,
    ln γ₁ = g + x₂ g′,   ln γ₂ = g − x₁ g′,   g′ = dg/dx₁,
    g′ = Σₖ aₖ [(x₂ − x₁) Lₖ(z) + 2x₁x₂ Lₖ′(z)],
    P_model = x₁γ₁P₁ˢᵃᵗ + x₂γ₂P₂ˢᵃᵗ,   rᵢ = (P_model − P)/P,

Lₖ by Bonnet's recursion, kLₖ = (2k − 1)zLₖ₋₁ − (k − 1)Lₖ₋₂, and
Lₖ′ = Lₖ₋₂′ + (2k − 1)Lₖ₋₁, each point at its own T for isobaric data, Gᴱ's
T-dependence over the set neglected. Then y₁ = x₁γ₁P₁ˢᵃᵗ/P_model, and the
data are consistent where the mean |y₁ − y₁,measured| ≤ 0.01. ln γ is
linear in the aₖ, so J is exact. Barker's fit can have two minima where the
data cover part of the range, so it starts from a = 0 and from the Gᴱ the
measured vapour implies, gᵢ = x₁ ln(y₁P/x₁P₁ˢᵃᵗ) + x₂ ln(y₂P/x₂P₂ˢᵃᵗ)
fitted linearly, and keeps the lower Σr².

## Activity models

M6's liquid models ([activity](../activity/model.v), D2 in
[m6-plan.md](m6-plan.md)) give each component's activity coefficient γᵢ at
a temperature and mole fractions, and nothing of phases: the γ–φ liquid is
built on them. A model holds nothing about an evaluation, so threads share
one, each with its own `Work`, which keeps what depends on T alone (each
τ, and UNIFAC's ψ and pure-group terms) until T moves, as M8's D11 keeps a
rate's K. A T the model refuses leaves nothing kept, since a refusal can
come with some of that T's parameters written.

**An evaluation answers or refuses,** with an `ActivityError` that names
what's wrong:
- **The composition:** each xᵢ ≥ 0 and finite, and |Σxᵢ − 1| ≤ 1e-12, which
  is far above normalising's rounding. All zero has no answer, and
  UNIQUAC's, UNIFAC's and Wilson's formulas give an unnormalised x a wrong,
  finite one (UNIQUAC's ln(rᵢ/R̄) moves by −ln Σxᵢ).
- **A parameter that's divided by or whose log is taken,** NRTL's G =
  exp(−ατ), τ in UNIQUAC, Λ, ψ: at least the smallest normal double,
  2.2250738585072014e-308, with finite T-derivatives. A subnormal one has
  lost digits, so its ln γ would be finite and wrong, and at 0 it's NaN.
- **The answer:** whatever was asked for, ln γ, its T-derivatives or dn, is
  finite, in every component.

**The models,** each with τᵢᵢ = 0 (NRTL) or 1 (UNIQUAC, Wilson's Λᵢᵢ):

- **Wilson** (1964): ln γᵢ = 1 − ln Sᵢ − Σₖ xₖΛₖᵢ/Sₖ, with Sᵢ = Σⱼ xⱼΛᵢⱼ.
- **NRTL** (Renon and Prausnitz, 1968): ln γᵢ = Aᵢ/Bᵢ + Σⱼ (xⱼGᵢⱼ/Bⱼ)(τᵢⱼ −
  Aⱼ/Bⱼ), with Aⱼ = Σₖ xₖτₖⱼGₖⱼ, Bⱼ = Σₖ xₖGₖⱼ and Gᵢⱼ = exp(−αᵢⱼτᵢⱼ),
  αᵢⱼ = αⱼᵢ.
- **UNIQUAC** (Abrams and Prausnitz, 1975), in DECHEMA's form, q in both
  terms and a coordination number of 10:

      ln γᵢ = ln(rᵢ/R̄) + 5qᵢ ln((qᵢ/Q̄)/(rᵢ/R̄)) + lᵢ − (rᵢ/R̄)Σⱼ xⱼlⱼ
              + qᵢ(1 − ln Sᵢ − Σⱼ θⱼτᵢⱼ/Sⱼ),

  R̄ = Σ xⱼrⱼ, Q̄ = Σ xⱼqⱼ, θⱼ = xⱼqⱼ/Q̄, Sᵢ = Σⱼ θⱼτⱼᵢ and lᵢ = 5(rᵢ − qᵢ) −
  (rᵢ − 1). ln(rᵢ/R̄) is ln(Φᵢ/xᵢ), written so that an absent component's
  xᵢ = 0 divides nothing.
- **UNIFAC**, each component i a set of subgroups k counted νₖ⁽ⁱ⁾, with rᵢ =
  Σₖ νₖ⁽ⁱ⁾Rₖ and qᵢ = Σₖ νₖ⁽ⁱ⁾Qₖ:

      ln γᵢ = ln γᵢᶜ + Σₖ νₖ⁽ⁱ⁾(ln Γₖ − ln Γₖ⁽ⁱ⁾),
      ln Γₖ = Qₖ(1 − ln Σₘ Θₘψₘₖ − Σₘ Θₘψₖₘ/Σₙ Θₙψₙₘ),
      Θₘ = QₘXₘ/Σₙ QₙXₙ,   Xₘ = Σᵢ νₘ⁽ⁱ⁾xᵢ / Σᵢ Σₙ νₙ⁽ⁱ⁾xᵢ,
      ψₘₙ = exp(−(aₘₙ + bₘₙT + cₘₙT²)/T),
      ln γᵢᶜ = 1 − Vᵢ′ + ln Vᵢ′ − 5qᵢ(1 − Vᵢ/Fᵢ + ln(Vᵢ/Fᵢ)),
      Vᵢ = rᵢ/Σⱼ xⱼrⱼ,   Fᵢ = qᵢ/Σⱼ xⱼqⱼ,

  ln Γₖ⁽ⁱ⁾ being ln Γₖ in pure i. In the original and LLE forms Vᵢ′ = Vᵢ and
  b = c = 0; in Dortmund's, Vᵢ′ = rᵢ^(3/4)/Σⱼ xⱼrⱼ^(3/4). ψ is taken between
  main groups, 1 within one. Each ln Γₖ⁽ⁱ⁾ is computed by the mixture's own
  arithmetic on xᵢ = 1, so that a pure component's γ is exactly 1.

**The derivatives, by forward differentiation.** Each formula is written
once, generic over f64 and `numeric.Jet`: a value with its first two
derivatives along one direction, whose arithmetic carries them,

    (ab)′ = a′b + ab′,   (ab)″ = a″b + 2a′b′ + ab″,
    (a/b)′ = (a′ − qb′)/b,   (a/b)″ = (a″ − 2q′b′ − qb″)/b,   q = a/b,
    f(a)′ = f′(a)a′,   f(a)″ = f″(a)a′² + f′(a)a″.

A logarithm's f″ = −k/a² (k = 1, or 1/ln 10 for log₁₀) overflows for a
below about 1e-154, where the derivative itself needn't. There, and only
there, its second derivative is taken in its other form,

    (k ln a)″ = k·(a″/a − (a′/a)²).

The f64 instance gives ln γ; the jet's value is the same bits.
- **In T,** the parameters as jets in T, each τ's own derivatives from
  `expr`'s jet pass over its tree (below), and x constant: ∂ln γᵢ/∂T and
  ∂²ln γᵢ/∂T² at constant composition.
- **In composition,** one pass for each xⱼ seeded, the parameters constant:
  ∂ᵢⱼ, the formula's derivative with the x's taken as independent. Since
  xₖ = nₖ/n, ∂xₖ/∂nⱼ = (δₖⱼ − xₖ)/n, so whatever the formula does off the
  simplex,

      n·∂ln γᵢ/∂nⱼ = ∂ᵢⱼ − Σₖ xₖ∂ᵢₖ.

  Gibbs–Duhem makes Σᵢ xᵢ·n∂ln γᵢ/∂nⱼ = 0, and n∂ln γᵢ/∂nⱼ is symmetric, n
  times the Hessian of nGᴱ/RT; the tests check both from scratch.

**The excess properties,** per mole, with R as `thermo.r_gas`:

    Gᴱ = RT Σ xᵢ ln γᵢ,   Hᴱ = −RT² Σ xᵢ ∂ln γᵢ/∂T,   Sᴱ = (Hᴱ − Gᴱ)/T,
    Cpᴱ = ∂Hᴱ/∂T = −RT(2Σ xᵢ ∂ln γᵢ/∂T + TΣ xᵢ ∂²ln γᵢ/∂T²).

`thermo.r_gas` is `quantity.molar_gas_constant`, reflux's one R: N_A·k_B =
8.31446261815324 J/(mol·K), exactly, as `expr`'s `R`, the kinetics' and the
Python `thermo`'s are. Until 0.6.0 `thermo.r_gas` and `activity.r_gas` were
8.314462618, the exact value rounded, 1.84e-11 below it, so a Gibbs reactor
or a PFR mixed two R's. ln γ doesn't involve R.

**Parameters** are expressions in T, dimensionless, parsed by `expr` (D7):
any paper's form goes in as printed. Their T-derivatives come from the jet
pass: `Expr.jet` evaluates the tree with each node on jets, its value
`eval`'s bits (a whole power by repeated multiplication, as `raise`), so a
τ's τ′ and τ″ are exact. ChemSep's banks give A in cal/mol and state no R.
Each is written with DECHEMA's, 1.98721 cal/(mol·K), the calorie cancelling.
That's confirmed for NRTL's bank (M6 step 0), and inferred for UNIQUAC's and
Wilson's, as is Wilson's form, until M6's step 12:

| Model | ChemSep's bank, written as | |
|---|---|---|
| NRTL | τᵢⱼ = `Aᵢⱼ K / (1.98721*T)`, α as printed | 352 lines, 329 pairs |
| UNIQUAC | τᵢⱼ = `exp(Aᵢⱼ K / (-1.98721*T))` | 291 lines, 273 pairs |
| Wilson | Λᵢⱼ = `Vⱼ / Vᵢ * exp(Aᵢⱼ K / (-1.98721*T))`, V each compound's Wilson volume | 366 lines, 362 pairs |

A line names its pair either way round; A₁₂ is its first compound's with the
second, turned to the model's order. A pair listed twice takes its first
listing. A pair the bank doesn't list is refused by name, unless the caller
lets it have no energy of interaction. With NRTL's τ = 0 or Wilson's Λ = 1
that pair is ideal, γ = 1 in their binary. UNIQUAC's τ = 1 removes only the
residual term: the combinatorial one, in r and q, stays, so the pair is
reported as having no residual interaction, not as ideal.

**UNIFAC's tables** are ChemSep's base files with their extensions
([components/activity.v](../components/activity.v)): main groups and their
parameters' rows, then subgroups. ChemSep numbers C=C as subgroup 9, as the
1987 *Properties of Gases and Liquids* did, so its 10 to 70 are the
published 9 to 69. Where a file contradicts itself the reader repairs it,
by the file's own evidence, and reports it:
- **A subgroup** whose main group's name is another main group's in the
  list, and neither name begins the other (so not a spelling, as Morph for
  Morpholine), takes that main group: furfural's LLE subgroup, which
  pointed to DMFA, and AC-CN's and AC-Br's, swapped.
- **A parameter row** whose name isn't its number's: an a row takes its
  names' pair when its numbers' pair is listed by another row too, or a
  number isn't in the list; a b or c row takes whichever pair the a rows
  have. Dortmund's c rows copy a name from the row above, and the numbers
  stand there.
- **A number** the list hasn't, named nothing it names, is a new main
  group, named by its rows.

And where the published table (DDBST's public pages) shows a typo,
`unifac_corrections` sets the published value, applied only where the file
still holds the shipped one: the original table's 9 parameters and 2
subgroups, the LLE table's 2 subgroups, and Dortmund's quaternary C, whose Q
is 0, with six ionic-liquid subgroups. [data-quality.md](data-quality.md)
lists every repair, and each table against `thermo`'s copy.

## The γ–φ liquid

A γ–φ fluid ([thermo/gamma_phi.v](../thermo/gamma_phi.v), M6's D4 to D6) is
a phase model whose liquid comes from an activity model and each
component's vapour pressure and liquid density, and whose vapour is the
ideal gas or an equation of state's. It has no F(T, V, n) behind its
liquid, so the flash takes it as it takes a cubic, and the TV and UV
flashes and the reactors refuse it. `Database.gamma_phi(names, options)`
([components/gamma_phi.v](../components/gamma_phi.v)) builds one from
ChemSep's banks, tables and component files.

**A solvent's fugacity** is xᵢγᵢPᵢˢᵃᵗφᵢˢᵃᵗ exp(Vᵢᴸ(P − Pᵢˢᵃᵗ)/RT), so

    ln φᵢᴸ = ln γᵢ(x′) + ln Pᵢˢᵃᵗ + ln φᵢˢᵃᵗ + Vᵢᴸ(P − Pᵢˢᵃᵗ)/RT − ln P,

- **γ** from the activity model over the solvents alone, at x′, their own
  mole fractions, x′ᵢ = xᵢ/s with s = Σ xⱼ over the solvents;
- **Pˢᵃᵗ and Vᴸ** from the component's own correlations, Vᴸ = 1/ρᴸ, the
  liquid incompressible over the Poynting integral;
- **φˢᵃᵗ** the vapour model's φ of the pure component at the correlation's
  Pˢᵃᵗ, on its vapour root, and 1 for the ideal gas;
- **Poynting's term** on by default, and off if asked, as Aspen's NRTL is.

**A Henry's-law solute's** is xᵢHᵢ,mix, so ln φᵢᴸ = ln Hᵢ,mix − ln P, with

    ln Hᵢ,mix = Σⱼ xⱼ ln Hᵢⱼ / Σⱼ xⱼ,

the sums over the solvents that have constants for it, an absent one
weighing nothing, which is `thermo`'s 'solvents_with_parameters' mode. Each
Hᵢⱼ, in Pa, is in one of two forms:
- **IAPWS's, for a gas in water** that its guideline covers (G7-04;
  Fernández-Prini, Alvarez and Harvey, 2003; 14 gases), with A, B and C
  from its Table 2:

      ln(kH/p*) = A/Tr + B·τ^0.355/Tr + C·Tr^−0.41·e^τ,   Tr = T/Tc, τ = 1 − Tr,
      ln(p*/pc) = (Tc/T)·Σ aᵢτ^bᵢ,

  on water's critical point (647.096 K, 22.064 MPa), p* being Wagner and
  Pruss's vapour pressure (IAPWS, 1992), whose six aᵢ and bᵢ the guideline
  gives. It takes the place of ChemSep's row (D6 as amended after step 6):
  ChemSep's hydrogen is 5.95 times it at 298 K and its neon 4.5e32 times
  (data-quality.md). Equation 121 can't hold it, so it's this form, on jets.
- **ChemSep's equation 121,** ln H = A + B/T + C ln T + D·Tᴱ, for every
  other gas and solvent.

A solute is outside the activity model, so its γ* is 1 (the
unsymmetric convention with no correction), it adds no volume, and it has
no Poynting term. ChemSep's banks pair water, the one solvent with
constants, with none of the gases, and its UNIFAC tables have no group for
hydrogen, CO, CO₂, nitrogen, oxygen or the noble gases, so there is no γ∞
to correct by for those (the light hydrocarbons, made of UNIFAC's own
groups, could take theirs from it, a later piece), and no gas's volume at
infinite dilution. (`thermo`
gives a solute the symmetric γ from the whole model, and Poynting on any
volume it's given; the two agree on an ideal solvent, which is how they're
compared.) ChemSep's bank has 41 gases, all in water, and IAPWS 14, all
but sulfur hexafluoride among them.

**Which compounds are solutes:** those named, and any with a constant in
another compound present whose Tc is below the lowest T at which a
compound that can only be a solvent has a vapour pressure: wherever the
liquid can be had, such a gas is above its critical point. A fluid where no
T is inside every solvent's correlations is refused when it's built, naming
the solvent above its critical point there.

**The liquid's volume** is Σ nᵢVᵢᴸ over the solvents, with no excess volume
(D4), and (∂V/∂T)_P = Σ nᵢ dVᵢᴸ/dT. It is incompressible, so (∂P/∂V)_T is
−∞, and it has no Cv and no speed of sound: NaN, as the seam allows.

**Which phase is which.** `phase(T, P, n, root)` gives the liquid, the
vapour model's vapour (the ideal gas, or a cubic's vapour root), or the one
of lower Σ nᵢ ln φᵢ, which is G less the ideal gas's at the same (T, P, n).
A state's label and derivatives come from the model that made it: a
liquid's volume is Σ nᵢVᵢᴸ to the bit, and a vapour whose volume ever
equalled it would be refused.

**Its derivatives**, all exact but φˢᵃᵗ's second:
- **In T,** through jets (`numeric.Jet`, as the activity models use):
  `tdep.Correlation.jet` gives each vapour pressure and density with its
  first two T-derivatives, for the forms a γ–φ liquid is made of (10, 16,
  101, 105 and 106); the activity model gives ∂ln γ/∂T and ∂²ln γ/∂T², and
  each ln Hᵢⱼ is a jet.
- **φˢᵃᵗ's:** d ln φˢᵃᵗ/dT = ∂ln φ/∂T + ∂ln φ/∂P·dPˢᵃᵗ/dT exactly, at
  (T, Pˢᵃᵗ(T)), and its second by Richardson's extrapolation of central
  differences of the first, steps h = 10⁻⁴T and h/2,
  (4D(h/2) − D(h))/3: an equation of state has no third derivatives of F.
- **In P:** ∂ln φᵢ/∂P = −1/P + Vᵢᴸ/RT, the second term for a solvent with
  Poynting's term only.
- **In n:** a solvent's n·∂ln φᵢ/∂nⱼ is the activity model's over x′,
  divided by s, and 0 in a solute's amount; a solute's is
  (ln Hᵢⱼ − ln Hᵢ,mix)/Σₖ xₖ in a solvent j with a constant, the sum over
  those, and 0 otherwise.

**Calorics from the fugacities** (D5), with R as `thermo.r_gas`:

    Hʳ = −RT² Σ nᵢ ∂ln φᵢ/∂T,   Gʳ = RT Σ nᵢ ln φᵢ,   Sʳ = (Hʳ − Gʳ)/T,
    Cpʳ = ∂Hʳ/∂T = −RT(2Σ nᵢ ∂ln φᵢ/∂T + T Σ nᵢ ∂²ln φᵢ/∂T²).

So a pure liquid's heat of vaporisation is Clapeyron's, T·ΔV·dPˢᵃᵗ/dT on
its own Pˢᵃᵗ and both phases' volumes, and the mixture's Hᴱ is the activity
model's: no ΔHvap correlation enters.

**Why φˢᵃᵗ at the correlation's Pˢᵃᵗ.** Only there is a pure liquid's
fugacity its vapour's at its own vapour pressure: at Pˢᵃᵗ, xγ = 1 and the
Poynting term vanishes, so f = Pˢᵃᵗφ(Pˢᵃᵗ), the vapour's. `thermo` takes
φˢᵃᵗ at the equation of state's own Pˢᵃᵗ instead, which moves a pure
liquid's ln f off its vapour's by the difference of the two ln φ's: up to
3.3e-4 for methanol and water and 9.0e-4 for ethanol under PR, from 300 to
355 K. The tests move `thermo`'s answers onto reflux's φˢᵃᵗ before
comparing.

**What it refuses**, with its reason (D4): a T outside a solvent's vapour
pressure or density range, or a solute's Henry's-law range, for the liquid
and for the phase of lower G alike, never the vapour in the liquid's place;
a liquid of solutes alone; a solute with a constant in none of the solvents
present; and a T where the vapour model has no vapour root at Pˢᵃᵗ, so no
φˢᵃᵗ. A flash reports a refused phase as `no_volume_root`, the kind that
says the model has no such phase there.

**Kept until T moves** (D4): each solvent's ln Pˢᵃᵗ, Vᴸ and ln φˢᵃᵗ, each
ln Hᵢⱼ as jets, and the activity model's own `Work`, in a `GammaWork` that a
model's `scoped()` copy holds (the seam's workspace, above). A model without
one makes its own for each call, with the same bits. A refusal leaves
nothing kept.

## Reports

The report (M4's D11) is the solution in SI, with what a reader works out from
it written beside it ([report-format.md](report-format.md)). For a stream with
molar flows nᵢ, the flash's composition zᵢ (which an empty stream keeps) and
molar masses Mᵢ:

| | |
|---|---|
| mass flows | mᵢ = nᵢMᵢ, and ṁ = Σ mᵢ |
| mass fractions | wᵢ = zᵢMᵢ / Σⱼ zⱼMⱼ |
| volume flow | F·v, at the stream's T and P |
| each phase | its own state, which is one mole of it: x, h, s and v |

Report format 4 (M6's D13) adds how the fluid was built: an equation of
state's kᵢⱼ for every pair, or a γ–φ fluid's activity model, vapour and
solutes, each binary pair's parameters as the expressions the model
evaluates, and each pair's source, a bank's line with its DECHEMA reference,
the file, or none listed.

A unit's balance is the one `solve` checks (D8), on the inlets it was
calculated from: its worst component against its flow F, and its energy
against F·RT at its hottest inlet,

    maxᵢ |Σ_in nᵢ − Σ_out nᵢ| / F,    |Σ_in F h + Σ q − Σ_out F h| / (F·R·T_max)

and 0 for a unit with nothing coming in. The whole flowsheet's is the same
over its feeds and products, with T the feeds' mean weighted by flow.

**Numbers for reading** are rounded to n significant figures by
`math.round`, half away from zero, and printed as integers, so tcc and gcc
print the same. With e = ⌊log₁₀|x|⌋,

    m = round(|x| · 10^(n−1−e)),    and if m = 10ⁿ, e + 1 and again

For n − 1 − e below 0 it divides by 10^(e−n+1), which is exact, rather than
multiply by 10^(n−1−e), which isn't: 846064.9999999999 × 0.1 rounds to the
tie 84606.5, and up, where the quotient by 10 is below it. A log₁₀ that rounds
up just short of a power of ten needs nothing: m rounds up to 10ⁿ⁻¹ with it.

## The C API

M10a's library ([capi/](../capi), declared by
[capi/reflux.h](../capi/reflux.h)) is the one place reflux meets another
language, so its conventions are the same for every function (M10a's D2):

- **A result code,** 0 or negative for a failure of the whole call, with the
  reason written into a buffer the caller passes, NUL-terminated and cut at
  the start of a UTF-8 character. A per-point status comes in step 4.
- **Fixed widths:** `int32_t`, `int64_t`, `double` and pointers. Every array is
  allocated by the caller, so no memory of V's reaches it but handles.
- **Handles** are `int64_t` ids into a registry, a map behind a lock, from a
  counter that never repeats. An id at or below the counter that the map
  lacks was freed; any other isn't a handle. The registry is also what keeps
  an object alive: the collector can't see the caller's references.
- **Threads:** every function registers the calling thread with the collector
  for the length of the call, and unregisters it after, so a thread the
  collector never knew (Python's) can allocate and collect. The library
  starts itself on its first call, which mustn't race, so one call comes
  first, from one thread: the Python package makes it at import.
- **The first call's thread** is where the collector starts, and the
  collector keeps it registered as its main thread. On Linux a registered
  thread that has ended makes the next collection abort, so the library
  unregisters it at the end of that call, unless it's the process's first
  thread, which lives as long as the process. Then it registers per call,
  like any other. Windows drops an exited thread itself.
- **Exports** are exactly the header's functions: a `.def` file on Windows and
  a linker version script on Linux, both written from the header by
  [tools/python.vsh](../tools/python.vsh). On Windows the header is compiled
  with the module, so a declaration that disagrees with the code fails the
  build. That's also why the header has no `const`: V can't write it, so
  the definitions it generates take `char *`, and a `const char *` in the
  header would conflict.
- **Text** goes in as NUL-terminated UTF-8, checked, and comes out into the
  caller's buffer with the room it takes in `*needed`. Text that doesn't fit
  is `REFLUX_TOO_SMALL`, with nothing written, and the Python package asks
  again with the room it was told. What comes out is in formats already
  documented: a compound is its component file
  ([data-format.md](data-format.md)), the index is `index.toml`'s format, and
  a fluid's setup is TOML written by `components.Setup.to_toml`. Floats are
  V's shortest round-trip form, so they read back exactly.
- **Codes** say whose the failure is: `REFLUX_BAD_ARGUMENT` for anything the
  caller could get wrong, checked before the V API is called; `REFLUX_BAD_HANDLE`
  for an id that isn't a live handle of its kind (a fluid's where a database's
  belongs, say); `REFLUX_NOT_FOUND` for a compound lookup, whose reason names
  the closest; and `REFLUX_DATA_ERROR` for what the data can't do, which is
  what's left of `components`' failures once the arguments are checked.
- **Many points in a call.** A function over points takes their number n,
  and each per-point input with a count of 1 or n, so a scalar or one
  composition needn't be repeated. The Python package passes an array that's
  already of the points' shape, C-contiguous `float64`, as it is, and
  broadcasts the rest by NumPy's rules. Outputs are property-major,
  `out[k·n + i]`, so each property's values are contiguous: a NumPy array of
  the points' shape with no copy.
- **A status per point** (`REFLUX_POINT_*`, whose flash kinds follow
  `flash.Failure`'s order). A point without an answer gets NaN, and
  `REFLUX_NO_PHASE` for its label, and leaves the others as they would be
  alone. `err` holds the first failed point's reason, so a point's reason is
  the same call made for it alone (M10a's D7).
- **A point's inputs are checked first:** finite, with T and P above zero and
  z not negative nor all zero. thermo's own check is `t <= 0.0`, which NaN
  and ∞ pass.
- **Properties per mole** are the V API's extensive ones divided by Σz, so a
  V program doing the same gets the same bits (criterion 3). A compound's
  vapour pressure is the fluid's cubic for it alone, with its own α.
- **One flash function, one layout** (M10a's D6). `reflux_flash` takes a
  spec saying what a point's two values are, (T, P) for the PT flash. It
  writes the whole's β, T, P, H, S and V per mole of feed, and each phase's
  fraction, V, H, S and composition in the place of its label: the
  liquid's (0) then the vapour's (1). A single phase takes its label's
  place, the other NaN, so the vapour's composition is always the same
  column of the arrays. An answer with two liquids, which M6's D9 names
  liquid and liquid2, or with three phases, has no place in the layout, so
  it is `REFLUX_POINT_MORE_PHASES`, its reason naming the phases, until ABI 2
  (M6's D14 as settled in step 7). The whole's
  H, S and V are its phases', summed in the flash's order, as
  `flowsheet.Stream` sums them.
- **Held to the V API bit for bit.** `tools/capi_reference.v` computes a
  reference set of states with `components`, `thermo` and `flash` alone, as
  the header documents each point, and the Python tests compare the
  library's answers with it to the bit. Both are built by one compiler with
  one set of flags, because bits differ between compilers and maths
  libraries (M10a's D12). The flags include `-ffp-contract=off`. Where the
  target has FMA, gcc fuses a·b + c into one rounding wherever its inlining
  lets it, and the library (`-shared`, without LTO) and the reference (with
  it) would differ. Ubuntu's CI found it (M10a's exit review).
- **The specs** are PT, PH, PS, TV, UV, and the vapour-fraction flashes TVF
  and PVF, whose β of 0 or 1 makes a bubble or dew point. Each checks its own
  two values before the feed. A vapour-fraction flash gives both phases,
  even one with no moles, whose composition is the incipient phase's.

### Where the data comes from

A **data root** is a directory with `components/`, a file per compound, and
`vendor/chemsep/`, ChemSep's kᵢⱼ banks (`pr.ipd`, `srk.ipd`) and
Mathias–Copeman constants (`mc-PR.lib`). `components.open(root, extra)` opens
one, with directories of the caller's own compounds after its own, and the
`Database` it returns keeps its root's `vendor/chemsep/`, so a fluid's kᵢⱼ and
α come from the same copy of the data as its compounds (M10a's D9). The
default root, `components.data_root()`, is `$REFLUX_DATA` if that's set, and
otherwise the source tree's `data/`, whose path V fills in when it compiles
(`@VMODROOT`). So a library or a binary built from the source tree finds the
data wherever it's moved to, as long as the tree stays, and a copy of the
data can go anywhere with `REFLUX_DATA` pointing at it. Before M10a,
`REFLUX_DATA` named the compounds' directory alone, and the kᵢⱼ always came
from the source tree.

## How equation-oriented solving will fit

M13 solves a whole flowsheet's equations at once, by sparse Newton (roadmap
J1 and J2), where M4 calculates unit by unit. Each unit then needs its model as
residuals, r(x_in, x_out) = 0, in the variables the tears already use: each
stream's x = (n₁ … n_c, P, h) (D6). A converged sequential-modular solution
satisfies them, so it's the equation-oriented solve's first point.

- **Balances,** for every unit: Σ_in nᵢ − Σ_out nᵢ, and
  Σ_in F h + Σ q − Σ_out F h. `UnitResult.balance` evaluates both already. A
  splitter's outlets are fractions of its inlet, and a separator's of each
  component.
- **Specifications:** the outlet's P or ΔP; then its T, β or a duty for a heater
  or a drum; h_out = h_in for a valve; s(T_s, P_out, n) = s_in at a machine's
  isentropic point, with h_out from η; h_out − h_in = v·ΔP/η for a pump.
- **Equilibrium** for a two-phase outlet: ln Kᵢ + ln φᵢⱽ − ln φᵢᴸ = 0 and
  Rachford–Rice, whose derivatives in T, P and n `thermo` already gives. A
  phase that appears or vanishes changes the equations themselves, which
  Newton can't see. So an equation-oriented solve keeps each outlet's number of
  phases, and tests stability at the end, as `flash.pt` verifies its split.

`Unit` would gain a `residuals` method beside `calculate`, with its Jacobian
blocks, and the sequential-modular solve stays the way in.

The same Jacobian gives M8's sensitivities exactly, where its D10 takes them
by differences of solves. At a converged point r(x, p) = 0, so
dx/dp = −(∂r/∂x)⁻¹ ∂r/∂p, one factorisation for every input. An output
y(x, p) then has dy/dp = ∂y/∂p + ∂y/∂x · dx/dp. An optimiser could also take
the stream variables as its own unknowns, the infeasible path, and meet the
residuals only at its optimum (M8's D13).

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
SRK covolume. `volume_function(t, n)` fixes what depends on T alone (the
association strengths Δ_AB(T) and SRK's a(T, n)) and solves the site
balance at each V it's asked for.
`model.v`, `saturation.v` and every flash routine stay as they are.

## Performance

These notes respond to the slow solves and slow sensitivity runs in DWSIM.

- Analytic derivatives are fast, and they give Newton methods their
  Jacobians for free.
- **What M3 step 10 measured** (`docs/benchmarks.md`, with a profile of the
  same flashes):
  - The volume roots took 88% of a PT flash, because every Newton step on
    V recomputed A(T, n). Fixing T and n first (`volume_function`, above)
    and sharing A with the residual made every flash 2.3–5.7 times faster,
    timed back to back against the old code.
  - Each `solve_tp` then still computed A twice, once for the root and
    once for the state there. `VolumeFunction.residual` shares it.
  - `Correlated`'s quadrature of the ideal-gas Cp is 0.6% of a PH flash,
    so it stays as it is.
  - A PH flash is about 12 PT flashes, each with its stability test and
    each split verified. That, not H, is where its time goes, and it's
    what is left between reflux and DWSIM on PH.
  - Each residual call still allocates a few small arrays (√aᵢ, ∂A/∂nᵢ,
    F_n). With the above fixed they're no longer near the top of the
    profile.
- **What M4 step 9 measured** (`docs/benchmarks.md`, and a profile of the
  reference flowsheet):
  - PH flashes were 95% of a flowsheet's solve, and stability tests 84%: 16
    PH flashes a solve, each about 8 PT flashes, each about 2.2 stability
    tests of about 7.5 trials.
  - So only a PH or PS flash's answer is verified, and a split's feed test
    stops at Wilson's two trials when they prove it unstable (above). Timed
    back to back against the old code, PH's median is 2.2, 2.3 and 3.7 times
    faster for 2, 5 and 10 components, and the reference flowsheet's solve 2.5
    times. In the benchmark, beside DWSIM: PH in 1.6, 4.0 and 5.9 ms against
    DWSIM's 3.5, 4.8 and 5.2, and the flowsheet in 26 ms against 559.
- **No closures in the solvers.** V 0.5.2 never frees a closure. When the
  solvers took function values, every PH flash leaked 18 KB and every
  bubble point 3 KB, a few closures each, which a sensitivity study would
  have turned into gigabytes. `numeric`'s solvers take interfaces
  (`Function`, `FallibleFunction`, `System`, `Objective`), and each caller
  passes a small struct with one method. The heap now stays flat over
  thousands of flashes.
- Sensitivity analysis: warm-start each case from its neighbour's converged
  solution, run independent cases in parallel with `spawn`, and later get
  exact local sensitivities from the converged Jacobian (implicit function
  theorem) instead of re-solving the flowsheet.
- Benchmark only `-prod` builds made with gcc or clang. tcc doesn't optimise.
