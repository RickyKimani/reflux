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

## Phase identification

`thermo.pip` ([thermo/phase.v](../thermo/phase.v)) is the phase
identification parameter of Venkatarathnam and Oellrich (2011):

    Π = V·[(∂²P/∂T∂V)/(∂P/∂T) − (∂²P/∂V²)/(∂P/∂V)]

`phase_label` calls a single phase liquid when Π > 1 and vapour otherwise,
as Caleb Bell's `thermo` does. It needs no saturation curve or critical
point, so it works for any mixture, supercritical ones included.

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

**Newton on G.** After 20 steps it switches to minimising the Gibbs energy
of the split over the vapour mole numbers v (one mole of feed, l = z − v):

    G/RT = Σᵢ vᵢ ln fᵢ(y) + Σᵢ lᵢ ln fᵢ(x),   gᵢ = ln fᵢ(y) − ln fᵢ(x),
    Hᵢⱼ = (δᵢⱼ/yᵢ − 1 + Φⱽᵢⱼ)/V + (δᵢⱼ/xᵢ − 1 + Φᴸᵢⱼ)/L,

where Φ is each phase's n·∂ln φᵢ/∂nⱼ. This uses numeric.minimize, with
fscale set to the size of the terms G is summed from. Near the critical
locus, 44 of the reference set's 166 splits end this way.

**Accepted** only with every |ln fᵢⱽ − ln fᵢᴸ| ≤ tol (1e-12) and 0 < β < 1,
and a Newton result is held to the same test. A split that slides onto the
feed (every |ln K| < 1e-4) is refused: there the fugacities are equal too,
and it would pass for two identical phases.

**Verified** (`verify`, on by default): each phase is re-tested for
stability, with the other phase as a trial. An unstable phase means a third
phase or a better split. That's a FlashError (`more_phases`) rather than an
answer, since M3 handles two phases. The water, n-decane and methane
mixture above is the example. There, a two-phase flash, reflux's without
the re-test or `thermo`'s FlashVL, returns a split whose fugacities are
equal but which isn't the equilibrium.

**Which phase is the vapour.** The phase of lower *mass* density, ΣnᵢMᵢ/V
(D6, amended in step 6). Molar volume isn't a substitute. In methane with
n-decane at 300 K and 200 bar, the methane-rich phase has the smaller molar
volume and yet weighs 168 kg/m³, against the decane-rich phase's 585. Step
4 compared molar volumes and called the decane-rich phase the vapour.

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

## PH and PS flashes

`flash.ph(fluid, P, H, z)` and `flash.ps(fluid, P, S, z)`
([flash/ph.v](../flash/ph.v)) find T and the split from a molar H or S, on
the Fluid's reference state. They take a Fluid, not a Model, because H and
S need the ideal gas.

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
  H_spec (or S_spec), by bisection, or the caller's `estimate`.
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
