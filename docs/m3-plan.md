# M3 plan — flash

M3 makes reflux answer the question every unit operation asks: given a
mixture and two specifications, how many phases are there, and what are
they? It's planned like M1 and M2: research first, then exit criteria,
decisions and steps. The roadmap entry is in [roadmap.md](roadmap.md).

## Goal

Equilibrium flash calculations for any `thermo.Model`:

- **stability analysis:** does a phase split?
- **PT flash:** the split at given temperature and pressure;
- **PH, PS, TV and UV flashes**, and **vapour-fraction flashes**: the split
  at other pairs of specifications;
- **bubble and dew points**;
- a **Newton solver** with line search and trust region (roadmap track B2),
  and a **benchmark harness** (A5).

The flash is vapour–liquid only: two phases at most. Three-phase VLLE and
free water are M6. Right after M3 comes M10a, thermo and flash from Python,
so the API is designed for that from the start.

## What the research found

**What M1 and M2 give the flash.**
- `thermo.ln_phi_derivs` returns ln φᵢ and its T, P and composition
  derivatives at a solved state. Those are exactly what a second-order
  (Newton) flash and stability analysis need.
- `Root.stable` picks the volume root with the lower Gibbs energy.
- `Residual2` carries every second derivative of F.
- `numeric` has a dense LU solver.
- M2's `db.fluid` builds PR or SRK fluids with real kᵢⱼ for 642 compounds,
  640 of them with an ideal-gas Cp, which gives the robustness sweep real
  mixtures to draw from.

**The algorithms.** They're standard and well documented:
- Michelsen & Mollerup, *Thermodynamic Models: Fundamentals and
  Computational Aspects*, chapters 9–11, already the reference in
  [design.md](design.md);
- Michelsen 1982, *The isothermal flash problem*, Part I (stability) and
  Part II (phase split), in Fluid Phase Equilibria 9;
- Whitson & Michelsen 1989, the negative flash;
- Michelsen 1999, *State function based flash specifications*, Fluid Phase
  Equilibria 158–160, 617–626.

What they say:
- **Stability:** minimise the tangent-plane distance (TPD) over a trial
  phase W:

      tm(W) = 1 + Σ Wᵢ (ln Wᵢ + ln φᵢ(W) − dᵢ − 1),  dᵢ = ln zᵢ + ln φᵢ(z)

  - Start from Wilson's K, with one vapour-like trial (Wᵢ = zᵢKᵢ) and one
    liquid-like (zᵢ/Kᵢ).
  - Iterate by successive substitution on ln W, accelerated by Michelsen's
    dominant-eigenvalue method (DEM).
  - Finish with a second-order step in αᵢ = 2√Wᵢ, whose Hessian is
    well-conditioned.
  - tm < 0 means unstable, and the trial phase is a good first guess for
    the split. A trial that converges back to z is the trivial solution
    and tells you nothing, so it must be detected.
- **PT flash:**
  - Solve Rachford–Rice for the vapour fraction β inside the window where
    it has no poles (Leibovici & Neoschil). Nichita & Leibovici 2013 give a
    convex transformation that converges monotonically, a robust
    alternative.
  - Iterate ln Kᵢ = ln φᵢᴸ − ln φᵢⱽ by successive substitution with DEM.
  - Switch to Newton on the Gibbs energy, in vapour mole numbers with a
    line search, when successive substitution slows.
  - The negative flash, which lets β lie outside [0, 1] while iterating,
    is what bubble and dew points and phase boundaries need.
- **PH, PS and UV flashes:** Michelsen 1999 writes each specification as
  finding a stationary point of a state function, his Q-function. This
  handles **pure components and narrow-boiling mixtures**. The textbook
  nested loop (guess T, PT-flash, correct T) fails on those: at fixed P,
  a pure component's H jumps at Tsat, so no T gives a specified H inside
  the jump.
- **TV and UV flashes** are natural for reflux. Its models are Helmholtz
  functions F(T, V, n), so these need no volume-root solve at all.

**References.**
- **`thermo` 0.6.1** (Caleb Bell), checked by inspecting the installed
  library:
  - `FlashVL` does TP, PH, PS, TV and vapour-fraction flashes with PR
    phases;
  - it labels phases by the phase identification parameter (PIP);
  - `chemicals.rachford_rice` has Leibovici–Neoschil and other
    Rachford–Rice solvers;
  - its ideal-gas Cp takes a polynomial (`poly_fit`), so reflux and
    `thermo` can share one Cp and agree to 1e-9 on H and S, as in M1.
- **DWSIM 9.0.5** is installed here, on .NET Framework 4.6.2.
  - DWSIM's modern Python route, the FluentAPI, needs DWSIM 10 or later,
    which you couldn't use.
  - A probe showed that `pythonnet` 3.0.5, on Python 3.13, loads DWSIM 9's
    `DWSIM.Automation.Automation3` from the install. It offers
    `CreateFlowsheet`, `AvailableCompounds`, `AvailablePropertyPackages`
    and `CalculateFlowsheet`.
  - A scripted capture of reference flashes is therefore possible, but
    unproven beyond loading the assembly, which is why step 0 is a spike.
  - DWSIM's compounds come from its own, older ChemSep copy, whose
    constants differ from M2's for some compounds. So the comparison
    records DWSIM's constants and kᵢⱼ and uses those: the same model with
    the same numbers.

**Phase labels.** A single phase has to be called liquid or vapour, and
unit operations and users both depend on the label.
- `thermo` uses the PIP of Venkatarathnam & Oellrich (2011):
  V·[(∂²P/∂T∂V)/(∂P/∂T) − (∂²P/∂V²)/(∂P/∂V)], liquid above 1.
- It needs ∂²P/∂V² and ∂²P/∂T∂V, which are third derivatives of F that the
  `Model` contract doesn't have. Finite differences of the analytic first
  derivatives are plenty for a label.

**Carried into M3:**
- from M1, the Soave α → 0 guard (supercritical components are where it
  bites), the benchmark harness, and the allocation work (M1 D8);
- from M2, the cost of `thermo.Correlated`, which integrates ChemSep's
  ideal-gas Cp numerically on every H and S. A PH flash evaluates H many
  times.

## Exit criteria

All checked by tests in CI, except the full sweep and DWSIM's capture,
which are tools whose committed output tests check.

1. **Stability is right.** On a reference set (binaries, ternaries, a
   10-component gas condensate, points near a critical point), stability
   analysis agrees with `thermo`'s on stable or unstable. The trivial
   solution is detected, never reported as a result.
2. **PT flash matches `thermo`.** On the same set, β and both phases'
   compositions agree to 1e-8 with `thermo`'s `FlashVL`, run with the same
   PR, constants and kᵢⱼ.
3. **Every converged split is an equilibrium**, in every test and sweep:
   - fugacities equal, |ln fᵢⱽ − ln fᵢᴸ| ≤ 1e-10;
   - each component's mass balance to 1e-12;
   - Σx = Σy = 1;
   - neither phase unstable (TPD ≥ −1e-8 on each).
4. **Every specification round-trips.** From each reference PT flash, take
   its H, S, V, U and β, then run the PH, PS, TV, UV and vapour-fraction
   flashes from them. Each must recover T, P and the split, to 1e-8
   relative. Bubble and dew points are β = 0 and β = 1 splits that satisfy
   criterion 3. The set includes pure components and a narrow-boiling
   mixture.
5. **Robust.**
   - A sweep of 10,000 random (T, P, z) points, PR with 2 to 10 components
     drawn from the database with its kᵢⱼ, has **zero unexplained
     failures**. D8 defines unexplained, and every failure is classified.
   - The sweep's report is committed. A fixed, seeded 1,000-point subset
     runs in CI.
6. **Agrees with DWSIM 9.** On about 20 cases (PT, PH, PS, bubble, dew),
   using DWSIM's own constants and kᵢⱼ, reflux agrees within DWSIM's
   convergence tolerance. The target is β and compositions to 1e-4, T to
   0.01 K and P to 0.01%. The deviations are recorded and can't grow.
7. **The Newton solver works.** It solves a standard set of test problems
   (from Moré, Garbow & Hillstrom) and reports non-convergence as an error,
   never a wrong answer.
8. **Speed is measured.** The benchmark harness times PT, PH and
   bubble-point flashes for 2-, 5- and 10-component mixtures, with gcc
   `-prod`. It compares reflux with `thermo` and with DWSIM on the same
   flashes, the last because slow flowsheets were a reason to start
   reflux. The results are kept in `docs/benchmarks.md`. There's no speed
   target yet: M3 makes the first honest measurement.
9. **Soave α is safe at any Tr.** α and its derivatives stay finite through
   the point where Soave's form reaches zero, and M1's fixtures are
   unchanged.
10. CI is green, and the docs are updated: README, roadmap, design notes
    and V notes.

## Decisions

| # | Decision | Choice | Why |
|---|---|---|---|
| D1 | Where it lives | A new module `flash`, importing `thermo` and `numeric`. | `thermo` stays about models, as `components` is about data. |
| D2 | Which models | Any `thermo.Model`: the flash uses only ln φ, its derivatives and volume roots. | CPA and SAFT (M12) must drop in without touching the flash. |
| D3 | Stability | Michelsen's TPD: Wilson vapour- and liquid-like trials; successive substitution with DEM acceleration every fifth step; a second-order finish in α = 2√W. Re-checking a split uses each phase as a trial too. | The standard, well-tested formulation, and the one `thermo` follows. |
| D4 | PT flash | Initial K from stability's unstable trial. Rachford–Rice in the Leibovici–Neoschil window, with Newton safeguarded by bisection. Successive substitution with DEM, then Newton on G in vapour mole numbers with a line search when it slows. Negative flash internally. | Fast where easy, robust where hard (near-critical), and the negative flash is what phase boundaries need. |
| D5 | PH, PS, TV, UV | Michelsen's state-function (Q-function) formulation. TV and UV work directly in F(T, V, n). | It handles pure components and narrow-boiling mixtures, where nested T-loops fail. |
| D6 | Phase labels | PIP, with finite-difference second derivatives of P. In a split, the less dense phase is the vapour. | `thermo`'s choice, so comparisons line up. No change to the `Model` contract. |
| D7 | The result | `FlashResult`: T, P, β, the phases (each a label, a `thermo.State` and its fraction), iterations, and which path converged. It allows N phases for M6. Failure is a typed `FlashError` (no convergence, bad specification, …), never a panic. | M10a wraps this for Python next: no panics across the C boundary, and failures carry a reason. |
| D8 | The robustness sweep | Seeded. Components drawn from the database among those with Tc, Pc, ω, M and an ideal-gas Cp. T and P log-uniform over wide ranges; z random. **Unexplained:** no convergence, a broken invariant (criterion 3), or a contradiction such as a single-phase answer that stability calls unstable. **Explained, and counted separately:** a specification outside the model's range, or a point within 1e-6 of a critical point, where the two phases can't be told apart. | "Zero failures" needs a definition, or it means nothing. |
| D9 | References | `thermo` 0.6.1 fixtures with the same PR, constants and kᵢⱼ, and a shared polynomial Cp, to 1e-8. DWSIM 9.0.5 through `pythonnet` and `DWSIM.Automation`, run by hand like the CoolProp fixtures, with DWSIM's own constants and kᵢⱼ recorded and used. **Fallback** if automation fails: capture by hand in the DWSIM GUI, following a written procedure. | `thermo` checks the implementation; DWSIM checks what M4 will be compared with. |
| D10 | Benchmarks | A `benchmark` module (roadmap A5), built with gcc `-prod`. Results, with machine and compiler, in `docs/benchmarks.md`. DWSIM's timings come from the same automation script. | Speed claims have to be real, and kept so that regressions show. |
| D11 | The Soave α guard | Compute Soave's α through its square root, s = 1 + m(1 − √Tr), so the zero at Tr = (1 + 1/m)² is harmless instead of a division by zero. Past it, keep the published form, which rises again, document that, and recommend Twu for supercritical components such as H₂. | The fix is exact and changes nothing below the zero, so M1's fixtures hold. |
| D12 | Who writes what | As before, I write it all unless you claim a step. | |

## Steps

```mermaid
flowchart LR
  S0[0 DWSIM spike] --> S5
  S1[1 Newton, Rachford–Rice] --> S3[3 Stability]
  S2[2 PIP, Soave guard] --> S3
  S3 --> S4[4 PT flash]
  S4 --> S5[5 Reference fixtures]
  S4 --> S6[6 Vapour fraction, bubble, dew]
  S4 --> S7[7 PH, PS]
  S7 --> S8[8 TV, UV]
  S6 --> S9[9 Robustness sweep]
  S8 --> S9
  S5 --> S9
  S9 --> S10[10 Benchmarks]
  S10 --> S11[11 Docs and exit review]
```

### 0. DWSIM automation spike
- From Python, through `pythonnet` and `Automation3`:
  - build a flowsheet with a methane/n-butane stream on PR;
  - run one PT flash and one PH flash;
  - read back β, the compositions and H;
  - read the constants and kᵢⱼ DWSIM used.
- Go or no-go for D9, written up in `docs/spike-dwsim.md`, as M1 did for
  the Python spike.

### 1. Newton solver and Rachford–Rice
- `numeric`: Newton with line search and trust region (B2), and a
  two-phase Rachford–Rice solver that's safe for extreme K (a heavy
  component's K can be 1e-20).
- Tests:
  - a Moré–Garbow–Hillstrom subset;
  - Rachford–Rice against `chemicals`' solvers and an mpmath reference, on
    random and extreme K;
  - β outside [0, 1] for the negative flash.

### 2. Phase identification and the Soave guard
- `thermo`: PIP (D6), and Soave α through its square root (D11).
- Tests: PIP labels match `thermo`'s; α stays finite and smooth through its
  zero; M1's fixtures are unchanged.

### 3. Stability analysis
- `flash.stability(model, t, p, z)`: TPD with the trials, successive
  substitution with DEM, the second-order finish, and trivial-solution
  detection (D3).
- Tests: agreement with `thermo`'s stability test on the reference set;
  TPD ≥ 0 at converged stable points; the trivial solution caught.

### 4. PT flash
- `flash.pt(model, t, p, z)`: stability, then the split (D4), returning a
  `FlashResult` (D7).
- Tests: criterion 3's invariants at every split; the hard cases (near
  critical, retrograde, a component far below its Tc beside one far above).

### 5. Reference fixtures
- `tools/gen_m3_fixtures.py` writes `testdata/m3_thermo.json`: `thermo`'s
  PT, PH, PS, TV and vapour-fraction flashes and bubble and dew points on
  the reference set, with a shared polynomial Cp.
- `tools/gen_m3_dwsim.py`, grown from the spike, writes
  `testdata/m3_dwsim.json`: about 20 DWSIM cases, with DWSIM's constants,
  kᵢⱼ and version.
- Tests: criteria 1, 2 and 6.

### 6. Vapour-fraction flashes, bubble and dew points
- β given with T, or with P; bubble (β = 0) and dew (β = 1) points.
  Newton on ln K together with T or P, started from Wilson.
- Retrograde cases have two dew points at one T, so the API says which one
  it returns.
- Tests: consistency with the PT flash; retrograde cases; clear errors near
  the cricondenbar and cricondentherm.

### 7. PH and PS flashes
- Michelsen's Q-function (D5), including pure components (H jumps at Tsat)
  and a narrow-boiling mixture.
- Tests: criterion 4's round trips; agreement with `thermo`.

### 8. TV and UV flashes
- The Helmholtz-natural specifications (D5).
- Tests: round trips from PT; agreement with `thermo` (TV).

### 9. Robustness sweep
- `tools/flash_sweep.vsh` runs the 10,000-point sweep (D8) and writes
  `docs/flash-robustness.md`, with every failure classified.
- A test runs the seeded 1,000-point subset in CI.

### 10. Benchmarks and a performance pass
- The `benchmark` module and `docs/benchmarks.md` (D10), with reflux,
  `thermo` and DWSIM on the same flashes.
- Then the carried items, measured before anything is changed:
  allocation in hot paths, and the cost of `Correlated` in PH and PS.

### 11. Docs and exit review
- README, roadmap, design notes on the flash, the V notes, and a
  criterion-by-criterion review, as in M1 and M2.

## Risks

| Risk | Effect | Fallback |
|---|---|---|
| DWSIM's automation can't be driven headless from Python | no DWSIM references | capture by hand in the GUI, following a written procedure (D9). The spike in step 0 finds out first. |
| Near-critical points converge slowly or to the trivial solution | sweep failures | the second-order finish, trivial-solution detection, and D8's explained class for true criticality |
| Pure and narrow-boiling PH flashes | wrong T or no convergence | the Q-function formulation (D5), with pure components in the tests from the start |
| `thermo` itself fails at hard points | no reference there | compare only where `thermo` converges, and record where it doesn't |
| The 10,000-point sweep is too slow under tcc | CI time | the full sweep is a tool with a committed report; CI runs 1,000 points |
| `Correlated`'s quadrature makes PH and PS slow | slow flashes | measure in step 10 first; then cache or tabulate Cp integrals if it matters |
| The Soave change perturbs M1's results | regressions | the square-root form is exact below the zero; M1's fixtures must still pass |

## Not in M3

- Three-phase VLLE and free-water handling: M6, with activity models.
- Phase envelope tracing (Michelsen 1980). Bubble and dew points come
  point by point here, and tracing a whole envelope comes later.
- Solids (SLE), electrolytes, and reactive or chemical-equilibrium flashes.
- γ–φ flashes with activity models: M6.
- The Python API: M10a, straight after M3.
