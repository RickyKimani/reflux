# M8 plan — specifications, sensitivity studies and optimisation

M8 makes a flowsheet answer questions, not just solve once.
- "Vary the cooler's duty until the drum's vapour is 5% methane."
- "Sweep the feed from 700 to 1000 °C and report the methane conversion."
- "Find the feed temperature and steam ratio that maximise hydrogen, with the
  outlet below 1100 K."
- "Fit these LHHW constants to my rate data, with their confidence
  intervals."

Its benchmark is the owner's. DWSIM 9 took 30 minutes in the app for a
50-point study of the reforming PFR (M9's Tutorial 1). M9 solved that study as
50 cold solves in 2.1 s. M8 adds warm starts, parallel cases, design
specifications, optimisation and fitting, and the speed the owner's other
models need, the published reforming model's K(T) above all. It's planned like
M1 to M4, M9 and M10a: research first, then exit criteria, decisions and steps.
The roadmap entry is in [roadmap.md](roadmap.md).

Status: **approved** on 2026-11-16, after a draft that split M8 in two: it's
one milestone, as the owner decided, with its optimiser written in V (D1,
D12). Step 0 is done and amended D6 and D18, and step 1 is done and amended
D11, both approved on 2026-11-17. Whether equation 16 is integrated by a
fixed rule (step 1) is still the owner's to decide. Step 2 is next.

## Goal

- **Design specifications** (F4): "vary X until Y = target", with X any number
  the file gives and Y an expression over the flowsheet's streams and units.
  Each spec is converged with its block's tears, not in a loop around the
  whole solve, and the answer is checked against the target like every
  balance.
- **Flowsheet variables in expressions:** a stream's T, P, flows and
  composition, and a unit's duty or work, by name, units checked. Specs, study
  outputs, objectives and constraints are expressions, as M9 promised (A9).
- **Sensitivity studies** (H1, H2): one or two inputs over a range, any
  expressions reported at each point.
  - Each point has its own status. A point that fails says why, and the rest
    stand.
  - A quantity with no value at a point, such as a product's conversion
    when the feed has none of it, is refused with its reason, never plotted
    as −∞ (Tutorial 2's lesson).
- **Warm starts** (F5): each point starts from its neighbour's answer, its
  tears and its units' own solvers alike, and is checked as a cold one is.
- **Parallel cases** with `spawn`. A study's answers don't depend on how many
  workers ran it.
- **Local sensitivities:** the derivative of each output with respect to each
  input at a point, with its error estimated.
- **Optimisation** (H3): an objective, minimised or maximised over chosen
  inputs within bounds, subject to constraints. The answer is returned only
  once it's shown to be an optimum, by its first-order conditions checked
  from scratch.
- **Fitting constants to data** (H4):
  - LHHW and power-law constants to rate data, on M9's kinetics with exact
    derivatives;
  - kᵢⱼ to VLE data.

  Each fit comes with its confidence intervals and its parameters'
  correlations.
- **What reactors recompute though their T hasn't moved** (M9 step 11): each
  reaction's K(T), 93% of the published reforming model's solve, and each
  compound's viscosity.
- **Against DWSIM 9:**
  - the owner's study against its sensitivity analysis;
  - a spec against its Adjust;
  - an optimisation against its optimiser;
  - a kᵢⱼ fit against its regression;

  every difference explained.

## What the research found

**What M8 builds on:**
- **`flowsheet.solve`** converges each block's tears by direct substitution,
  Wegstein or Broyden on x = (n₁ … n_c, P, h) (M4's D6, D7). Tears take first
  values from `SolveParams.estimates`. Nothing is returned unless every
  balance closes (D8).
- **A `Unit` holds nothing about a solve,** so calculating it twice gives the
  same answer. Units are therefore safe to share between threads. But a unit
  can't remember its last answer, which a warm start needs: that has to be
  passed in.
- **`Loaded.root` keeps the file's TOML tree,** so a study can replace one
  value and remake that feed or unit through its own maker. Every unit,
  including later ones, would then be variable without each writing setters.
- **`expr`'s names** are a context's variables (`T`, `P`) and family members
  (`p[CH4]`). Nothing names a stream's or a unit's variable yet. The tape
  differentiates with respect to its inputs, and a rate's constants are
  folded into it, so a fit must make them inputs.
- **`numeric.minimize`** is Newton with a Levenberg–Marquardt shift, for an
  unconstrained objective with its Hessian. `numeric` has no least squares
  and no constrained optimiser.

**Where a reactor's time goes** (M9 step 11, and checked for this plan):
- On the published model, 93% of a solve is `Kinetics.fill`. It calls
  `Reaction.standard`, which integrates each species' ideal-gas Cp from
  298.15 K at every evaluation of the rates.
- **Four of the reforming case's five species give Cp as ChemSep's eq. 16,**
  Y = A + exp(B/T + C + D·T + E·T²), which has no closed-form integral. So
  each integral is a Gauss–Kronrod quadrature at 1e-13 (`tdep`). Only water's
  eq. 107 integrates in closed form.
- **Across the database,** 521 of 640 compounds' ideal-gas Cp are eq. 16,
  116 eq. 107, 2 eq. 100 and 1 eq. 1.
- **So keeping K at the last T** fixes an isothermal reactor, whose T never
  moves. An adiabatic one moves T every step, and would need eq. 16
  integrated faster. Every PH flash integrates the same Cp, but M3 measured
  that at 0.6% of a flash.
- On DWSIM's typed rates, Wilke's viscosity is about a quarter of a solve.
  Each compound's correlation is evaluated again at every step, though T
  doesn't move.

**How DWSIM does it** (its source on GitHub, `DWSIM.FlowsheetSolver` and
`Adjust.vb`; its API reference and help):
- **An Adjust** varies one variable until another equals a target, by
  bisection or secant. Its defaults are a tolerance of 1e-4 and at most 10
  iterations.
- **Simultaneous adjusts** run Newton on every Adjust marked for it:
  - the Jacobian by central differences of ±1% of each variable;
  - each evaluation a whole flowsheet solve, so 2n + 1 solves an iteration;
  - at most 25 iterations, each adjust's own tolerance.

  Either way a spec is a loop around the flowsheet, not converged with its
  recycles.
- **Its sensitivity analysis** varies up to two variables, each over equally
  spaced points, and records dependent variables or expressions. Each point
  recalculates the flowsheet from the state the last left. That's a warm
  start of a sort, but an unchecked one: a point's answer depends on the
  order of the points.
- **Its optimiser:**
  - its methods are Brent, L-BFGS, truncated Newton, Nelder–Mead (each with
    a bounded form) and IPOPT;
  - the objective is a flowsheet variable or an expression;
  - its help describes bounds on the varied inputs and no other
    constraints;
  - it stops when the objective's change falls below its tolerance, so its
    answer isn't checked to be an optimum.
- **Its parallel mode:** the solver documents a "background parallel"
  mode. In the code read here, its parallel queue (`Parallel.ForEach`) is
  never reached, and every mode takes the serial path. Step 0 checks that in
  DWSIM 9.0.5 itself.
- **Its −∞:** the owner's Tutorial 2 swept a product's conversion,
  (in − out)/in with nothing in. DWSIM plotted −∞ at every point, with no
  message (M9 step 10).
- **Its regression:** a study's results can go to a data-regression plugin.
  Step 0 finds out whether that fits kᵢⱼ headless.

**Threads or processes,** measured for this plan:
- **The setup:** the owner's case (`reformer_dwsim.toml`), 48 feed
  temperatures, reflux built with gcc `-prod`, on the owner's laptop (Ryzen 5
  PRO 3500U, 4 cores).
- **One thread:** 1.67 to 1.78 s.
- **Threads with `spawn`,** sharing one `Loaded`: 1.11 s with 2 (1.5×) and
  0.81 s with 4 (2.1×).
- **Four processes** of 12 points each: 0.48 to 0.64 s of solving each
  (2.6×), and 1.09 s from start to finish, each process reading the file and
  the data.
- **What limits threads:** the collector's lock, which every allocation
  takes (M10a measured 1.3 to 2.3× with 8 threads), and the cores
  themselves.
- **What it suggests:** threads are worth having, and an allocation less per
  step would raise their ceiling. Processes gain a little more, at the cost
  of reading everything once each and merging reports.

**The methods** (independent of reflux):
- **Specs and tears.** In sequential-modular simulation a spec is converged
  either nested, by a 1-D method around the recycle loops, or together with
  the tears by a quasi-Newton method (Westerberg, Hutchison, Motard and
  Winter, *Process Flowsheeting*, 1979; Biegler, Grossmann and Westerberg,
  *Systematic Methods of Chemical Process Design*, 1997, ch. 4).
  - Together needs far fewer passes, since each Broyden evaluation is one
    pass, not a converged solve.
  - Nested is sturdier far from the answer. Aspen Plus offers both.
- **Warm starts along a sweep** are natural-parameter continuation. A secant
  predictor through the last two answers usually halves the corrector's
  work (Allgower and Georg, *Numerical Continuation Methods*, 1990).
- **A derivative by differences** has an error that two steps estimate:
  Richardson's extrapolation. A derivative of a solve can't be closer than
  the solve's tolerance over the step allows, so its tolerance must tighten.
- **Optimising a flowsheet:**
  - A flowsheet optimisation has few inputs, a handful to a few dozen, and
    each evaluation is a converged solve.
  - For problems of that size, sequential quadratic programming (SQP) is
    the standard: a quadratic model from a BFGS update, kept positive
    definite by Powell's damping, a QP subproblem on the linearised
    constraints, and a merit function to accept steps (Nocedal and Wright,
    *Numerical Optimization*, 2006, ch. 18).
  - Converging each evaluation fully is the "feasible path". Making the
    tears optimisation variables too, the "infeasible path", takes far fewer
    passes (Biegler and Hughes, AIChE J. 28, 1982), at the cost of the
    optimiser seeing every tear.
  - IPOPT is built for thousands of variables with sparse Jacobians: an
    equation-oriented flowsheet's (M13), not this one's.
  - Hock and Schittkowski's collection (*Test Examples for Nonlinear
    Programming Codes*, 1981) gives small constrained problems with their
    published optima.
- **Fitting:**
  - Nonlinear least squares by Gauss–Newton with a Levenberg–Marquardt
    shift.
  - Confidence intervals by linearising at the optimum, s²(JᵀJ)⁻¹, whose
    coverage is approximate for a nonlinear model (Bates and Watts,
    *Nonlinear Regression Analysis and Its Applications*, 1988).
  - Arrhenius constants centred on a reference T in the data's range, which
    removes most of A's and Eₐ's correlation ([kinetics.md](kinetics.md)).

**References for the tests:**
- **Closed forms:**
  - an isothermal first-order PFR, X = 1 − exp(−kPV/F): its volume for a
    target conversion, and dX/dT;
  - a heater's duty for an outlet T, against a PH flash;
  - series reactions A → B → C in a PFR, whose B peaks at a space time of
    ln(k₂/k₁)/(k₂ − k₁) (Levenspiel), for the optimiser.
- **mpmath,** for those derivatives and optima at high precision.
- **Hock and Schittkowski's problems,** for the SQP.
- **SciPy** (`least_squares`, with `thermo`'s PR for kᵢⱼ), fitting the same
  data in a uv script, and seeded synthetic data from known constants.
- **reflux's own cold solves,** for every warm one.
- **DWSIM 9:** its Adjust, its optimiser, its sensitivity analysis and its
  regression, as far as its API allows them headless (step 0).

## Exit criteria

All checked by tests in CI, as before.

1. **Specs are met.** A spec converges with its block's tears, and the
   solution meets its target to the solver's tolerance, checked from the
   solution itself. It holds on a block with a recycle and on one without:
   - an isothermal PFR's volume for a conversion matches its closed form;
   - a heater's duty for an outlet T matches the PH flash.

   A target outside what the varied value's bounds can reach is refused,
   with the range they reach. Every balance is still checked (M4's D8).
2. **Flowsheet variables bind.** Every name in a spec, an output, an
   objective or a constraint binds to a stream's or a unit's variable, with
   its dimension, or is refused with its line, its column and a suggestion.
   Each renders back.
3. **Studies report.** One or two inputs over values with units, and every
   output at every point, each point with its status.
   - A failed point gives its reason, and every other point stands.
   - An output with no value at a point is NaN with its reason, and
     Tutorial 2's methanol conversion is one.
4. **Warm starts are checked.** Every warm-started point agrees with its cold
   solve within the solver's tolerance, point by point. The passes and
   integrator steps each takes are recorded, and can't grow.
5. **Workers don't change answers.** A study run by 1, 2, 4 or 8 workers gives
   the same bits. The throughput is recorded, with no target.
6. **Sensitivities are honest.** Each derivative matches its closed form
   within its own error estimate. Each estimate's ratio to the true error is
   recorded, so that a loose estimate shows.
7. **The optimiser is right.** Hock and Schittkowski's problems reach their
   published optima within tolerance, from their published starts. Each
   answer's first-order conditions are checked from scratch: feasibility,
   the Lagrangian's gradient, the multipliers' signs and complementarity. A
   problem with no feasible point, or one that doesn't converge, is refused
   with its reason. Iterations are recorded, and can't grow.
8. **Flowsheet optimisation is right.** The series reactions' PFR reaches
   its closed-form optimum. One optimum sits on a bound and one on a
   constraint, and each is verified as in criterion 7. The same problem in
   DWSIM's optimiser is compared, and the difference explained.
9. **Fits are right.**
   - Rate data and VLE data fitted by reflux match SciPy's fits of the same
     data, parameters and intervals, to the optimum's tolerance.
   - Over many seeded synthetic data sets from known constants, the 95%
     intervals cover the true values about as often as they claim, within
     binomial bounds.
   - A kᵢⱼ fit is compared with DWSIM's regression, if step 0 can drive it.
10. **Nothing a reactor keeps changes an answer.** K(T) and the viscosity
    kept at the last T leave every golden report exactly as it was. The
    published model's solve is recorded before and after.
11. **The owner's study against DWSIM 9.**
    - The 50-point study is timed with gcc `-prod`, warm, cold and on
      workers, against DWSIM's own sensitivity analysis. If step 0 finds that
      can't be driven, it's timed against M9's loop of solves instead.
    - As the roadmap's exit asks, it must be faster, which M9 already shows,
      and the ratio is recorded.
    - A spec is compared with DWSIM's Adjust at a tight tolerance, every
      difference explained.
12. **Files and the CLI.**
    - `[spec.NAME]`, `[study.NAME]`, `[optimize.NAME]` and `[fit.NAME]` in
      the file, every key checked, every refusal at its line.
    - `reflux solve` converges specs. `reflux study`, `reflux optimize` and
      `reflux fit` each print a table and write JSON and CSV.
13. CI is green, and the docs are updated: README, roadmap, design notes,
    kinetics.md, the file and report formats, and the V notes. v0.4.0 is
    released.

## Decisions

| # | Decision | Choice | Why |
|---|---|---|---|
| D1 | Scope | One milestone, the roadmap's M8 (discussed on 2026-11-16): specs, warm starts, studies, workers, sensitivities, optimisation and fitting. The roadmap is amended in three places:<br>- B6: no C optimiser in M8, its own SQP (D12) instead, with IPOPT moving to M13;<br>- H1: exact sensitivities move to M13 (D10);<br>- H4: data reconciliation moves to M13. | The optimiser in V removes the C dependencies, NLopt and IPOPT, which were the case for splitting M8. A flowsheet optimisation is small and dense. Reconciliation needs measured data against the whole flowsheet's equations, M13's. |
| D2 | Where it lives | <br>- **`flowsheet`:** specs, since they change how a block converges.<br>- **A new module, `study`:** studies, workers, sensitivities and flowsheet optimisation, over `flowsheet`.<br>- **`numeric`:** the SQP and least squares, each over an interface.<br>- **A new module, `fit`:** fitting, over `kinetics`, `flash`, `thermo` and `numeric`.<br>- **`main`:** `reflux study`, `optimize` and `fit`.<br>- **`kinetics` and `reactors`:** caching K and the viscosity.<br><br>vlib has no module named `study` or `fit`. | A study and an optimisation are many solves; the solve itself stays one. Fitting is many evaluations of a model, not of a flowsheet. |
| D3 | Naming flowsheet variables | A name with one dot:<br>- for a stream, `product.T`, `product.P`, `product.F` (total molar flow), `product.n[CH4]` (a component's flow), `product.y[CH4]`, `product.h` and `product.vf` (vapour fraction);<br>- for a unit, its energy streams by their names, `reactor.duty` and `pump.work`.<br><br>The grammar's atom gains `name "." name`. A stream or unit that isn't there, or a variable it lacks, is refused with the closest names. Each variable has its dimension, so `product.T - 300` is refused as M9's sums are. | It reads like the file, where a stream is `[feed.syngas]` and a unit `[unit.reactor]`, and keeps M9's members for components. Families keyed by stream, `T[product]`, would need two keys for a composition, `y[product, CH4]`. |
| D4 | What can be varied | Any number the file gives a feed or a unit, by its key: `"feed.syngas.T"`, `"unit.cooler.duty"`, `"unit.reactor.volume"`. Each point replaces that value in the file's tree, with its unit, and remakes the feed or unit through its maker, which checks it as it checks the file. The value's dimension must be the key's. | Every unit, today's and later ones, is variable at once, with its own checks, and no unit writes setters. Remaking a unit costs microseconds against a solve's milliseconds; step 3 measures it. |
| D5 | Specs | `[spec.NAME]` has three keys:<br>- `vary`, a key, as D4;<br>- `between`, its bounds;<br>- `until`, an equation of two expressions whose sides share a dimension, `"drum_vapour.y[CH4] = 5 %"`.<br><br>In planning, a spec is an edge from the unit that makes what it measures back to what it varies (a unit, or the units its feed enters), so they land in one block. That block's unknowns are its tears and its specs' values, and Broyden converges them together (F4). The residuals are the tears' changes and each spec's two sides' difference, scaled.<br><br>A block with one spec and no tears is solved by Brent within the bounds instead. A spec that fails in a block with tears is tried once more nested, Brent around the block's tears, before it's refused.<br><br>The answer must meet each target to the solver's tolerance and stay within its bounds. A unit's value set equal to an expression of others, DWSIM's Spec block, is a spec that varies that value until it equals the expression. | Together is a pass an evaluation, where DWSIM's simultaneous adjusts take a whole solve. Brent can't fail inside a bracket, and a single spec is the common case. Bounds let a target out of reach be refused with what the bounds reach, rather than wandering. |
| D6 | Warm starts | A study's or an optimiser's point starts its tears from the last point's, with a secant predictor through the last two where both exist.<br><br>A unit that can use its last answer says so with an optional interface, `Warm`, whose `calculate_from(fl, inlets, last Outcome)` takes it. `Unit` is unchanged, and a unit without it starts cold:<br>- the CSTR starts Newton from its last outlet, skipping continuation;<br>- a reactor at a duty or adiabatic brackets its T around the last;<br>- a PFR starts from its last first step.<br><br>A warm answer is checked exactly as a cold one is.<br><br>**Amended in step 0, for step 5 to settle.** A tear started from its neighbour's answer, as `SolveParams.estimates` allows today, took more passes than a cold start on the reference flowsheet with Wegstein (11.0 a point against 4.0) and with Broyden (11.8 against 8.1), and fewer only with direct substitution (24 against 29). The predictor brought Broyden's to 7.8 and Wegstein's to 8.4. So a tear's value alone isn't the warm start. Step 5 first measures what is, Broyden's Jacobian carried along with the tears above all, and a block is started warm only where that takes no more passes than cold. | `Unit`'s promise, nothing about a solve, keeps units shareable between threads, so the last answer comes in from outside. A predictor is cheap and halves the work along a smooth sweep.<br><br>Step 0 measured the tears' part by hand, and it doesn't hold as written. |
| D7 | Determinism | A study's points go in chunks of a fixed size, 8 unless the file says otherwise. Each chunk starts cold at its first point and warm-starts the rest in order. Workers take whole chunks. A grid of two inputs runs along its rows, alternating direction, so each point's neighbour is the last one solved. | The same chunks give the same bits whatever the number of workers (criterion 5). DWSIM's answers depend on the order of its points. |
| D8 | Workers | Threads, with `spawn`, sharing one loaded file. The number is the CPU count unless `--workers` says otherwise. Step 7 profiles a solve's allocations and removes those in the hot path that it can. | Measured: 2.1× with 4 threads against 2.6× with 4 processes, which would each read the file and the data and need their reports merged. One process keeps one report. |
| D9 | A point's outcome | Each point has a status (solved, failed, or solved with an output that has no value) and reasons in words, as M10a's points have. An expression that divides by zero, or takes the log of zero, at a point is NaN with its reason ("methanol's conversion: the feed has no methanol"), never ±∞. A failed point doesn't stop the study, and the next point starts cold. | Tutorial 2's −∞ is what a sweep must never do. M10a's C API already works this way, so M10b can hand it to Python unchanged. |
| D10 | Sensitivities | The derivative of each output with respect to each input, at any point, by central differences of warm-started solves.<br>- The solve's tolerance is tightened to 1e-12 for them.<br>- Each step is chosen from the input's scale.<br>- The error is estimated from steps h and h/2 by Richardson's extrapolation, and reported beside the derivative.<br><br>The same derivatives are the optimiser's gradients. Exact sensitivities, by the implicit function theorem on the converged Jacobian (H1), need each unit's Jacobian, which comes with M13's residuals (J1). | Differences work for every unit, and give an error estimate that a user can judge. An optimiser needs gradients good to about 1e-6, which a solve at 1e-12 gives. The PFR's rows in T and P are themselves differences (M9 step 9), so an exact path would be partial today. |
| D11 | Keeping K(T) and the viscosity | A reactor's rates keep each reaction's K, and each compound's viscosity, at the T they were last computed at. A new T computes them as now, so every answer keeps its bits. Whether eq. 16's Cp should integrate faster, for reactors whose T moves, is measured in step 1 and decided there.<br><br>**Amended in step 1.**<br>- Wilke's rule keeps each pair's two factors with each compound's viscosity: its powers and square roots, not the correlations, were most of its cost.<br>- At a T that moves, each species' two integrals are taken once and shared by the reactions it's in, each summing them as before: 10 integrals for the reforming model where its three reactions took 24, and the same bits. | Exact by construction. An isothermal reactor, the owner's case, never moves T. |
| D12 | The optimiser | `numeric.sqp`, written in V, for small dense problems: minimise f(x) subject to bounds, equalities c(x) = 0 and inequalities c(x) ≥ 0.<br>- **The method:** Han and Powell's SQP. The Hessian of the Lagrangian comes from a BFGS update with Powell's damping, and each step's QP is solved by an active-set method on the bounds and the linearised constraints. Steps are accepted by an ℓ1 merit function with a backtracking line search (Nocedal and Wright, ch. 18).<br>- **What it takes:** the problem as an interface, with gradients, exact where the caller has them and D10's differences for a flowsheet.<br>- **What it returns:** an answer only once the first-order conditions hold to their tolerances, each multiplier with it. A QP with no feasible point, a line search that can't descend, or the iteration budget spent is a typed `OptimizeError` with its reason. | Every flowsheet optimisation M8 meets has a handful of inputs. SQP converges superlinearly there, takes constraints directly, and needs no dependency. Checking the conditions is how the optimiser shows its answer, as every other solver here does. DWSIM's stops when its objective stops changing. |
| D13 | Optimising a flowsheet | `[optimize.NAME]` has four keys:<br>- `vary`, keys with bounds, as D4;<br>- `minimize` or `maximize`, an expression;<br>- `subject_to`, equations and inequalities over flowsheet variables, `"product.T <= 1100 K"`;<br>- `start`, optionally.<br><br>The path is feasible: each evaluation is a converged solve, specs included, warm-started from the last (D6), and the constraints are the SQP's. The equation reader of D5 takes `<=` and `>=` too. The report gives the optimum, the multipliers, which bounds and constraints are active, and each evaluation. | The feasible path keeps every evaluation a checked solve. The infeasible path's fewer passes come from giving the optimiser the tears, which M13's equation-oriented solve does anyway. |
| D14 | Fitting kinetics | `[fit.NAME]` in a reactions file has four keys:<br>- `data`, a CSV whose header gives each column with its unit (`T [K]`, `P [bar]`, `y[CH4]`, `r[reforming] [mol/(kg*s)]`);<br>- `parameters`, by name: constants and definitions such as `k1.A`, `k1.Ea` or `K_CO.dH`;<br>- `start`, if not the file's values;<br>- `weights`, each rate's absolute or relative error.<br><br>The residuals are the measured rates less the model's. The fit is least squares by Gauss–Newton with a Levenberg–Marquardt shift (`numeric.least_squares`), on ∂r/∂θ from the tape, exactly: the fitted constants become the tape's inputs. An Arrhenius or van 't Hoff constant is fitted centred on the data's mean T, and reported in the file's own form.<br><br>The report gives each value, its standard error and 95% interval by linearisation, the correlation matrix, the residuals, and the fitted section as file text. Rate data only, as from a differential reactor; outlet data from an integral reactor, a PFR solved per evaluation, comes later. | M9's exact derivatives make the fit exact to rounding in its Jacobian. Centring removes most of A's and Eₐ's correlation (kinetics.md), and the correlation matrix shows what's left. Rate data is how kinetics is measured for fitting, and needs no reactor model. |
| D15 | Fitting kᵢⱼ | `[fit.NAME]` with VLE data (T, P, x, y), the kᵢⱼ of named pairs as parameters, and residuals in bubble pressure and vapour composition, each weighted. The derivatives are by differences, and the solver and report are D14's. | The data DWSIM's regression and thermo's users fit. kᵢⱼ is one number a pair, so differences are cheap and enough. |
| D16 | The file | Format 3 adds `[spec.NAME]`, `[study.NAME]`, `[optimize.NAME]` and `[fit.NAME]`.<br><br>A study has:<br>- `vary`, one or two, each a key with `from`, `to` and `points`, or `values`;<br>- `report`, named expressions;<br>- optionally `chunk` and `sensitivities` (which outputs, with respect to which inputs).<br><br>Formats 1 and 2 read as before. | A new section is a new format, by the CHANGELOG's rules. |
| D17 | The CLI and reports | `reflux solve FILE` converges the file's specs, and its report gains each spec's value and residual.<br><br>`reflux study FILE [NAME]`, `reflux optimize FILE [NAME]` and `reflux fit FILE [NAME]` each print a table, with `--json OUT` and `--csv OUT`, and the first two `--workers N`. The JSON gives every point or iteration, with its passes and reasons, and the CSV one row a point, NaN where an output has no value.<br><br>The exit status:<br>- 0 when the command ran, even with failed points, which the table says;<br>- 1 when it couldn't run, or found no optimum or fit;<br>- 2 for a wrong command.<br><br>Each report is a format of its own, documented beside the solve's. | A study that ran with some failed points is still a result. Each point's reason is in the report. |
| D18 | Comparing with DWSIM | <br>- The owner's 50-point study against DWSIM 9's own sensitivity analysis, if step 0 can drive it headless, else against M9 step 11's loop of solves.<br>- A spec on the reference flowsheet (M4) against its Adjust, its tolerance set to 1e-10.<br>- An optimisation within bounds against its optimiser.<br>- A kᵢⱼ fit against its regression, if step 0 can drive it.<br><br>Each is compared by state, as M4 compared.<br><br>**Amended in step 0** ([spike-dwsim-study.md](spike-dwsim-study.md)): a script can run each of its forms, unshown.<br>- The study: against its own form.<br>- The spec: the cooler's outlet T until the first drum's vapour is 60 mol/s, against a single Adjust's secant and Brent's at 1e-10, with the recycle's tolerances scaled by 1e-8. Not the simultaneous solver at that tolerance, which leaves the flowsheet 1% off its answer.<br>- The optimisation: one variable within bounds, against bounded Brent's, the only one of its 11 methods that found an optimum. The pressure between two compression stages, which step 0 ran, and criterion 8's series reactions if its PFR takes them.<br>- The kᵢⱼ fit: its regression ran and fitted nothing, so the comparison is recorded as not possible, unless the owner's run of the same nine points in the app gives the kᵢⱼ. | The comparison is of the same question asked of both, and its differences are each explained. |
| D19 | Versions | M8 ends with v0.4.0: format 3 and the new reports are additions, a minor release. | The CHANGELOG's rules. |
| D20 | Who writes what | As before, I write it all unless you claim a step. | |

## Steps

```mermaid
flowchart LR
  S0[0 Spikes, test_memory] --> S1[1 Keeping K and viscosity]
  S0 --> S2[2 Flowsheet variables]
  S0 --> S3[3 Remaking a feed or unit]
  S2 --> S4[4 Specs]
  S3 --> S4
  S3 --> S5[5 Warm starts]
  S2 --> S6[6 Studies, the CLI]
  S5 --> S6
  S6 --> S7[7 Workers]
  S6 --> S8[8 Sensitivities]
  S0 --> S9[9 SQP]
  S8 --> S10[10 Optimising a flowsheet]
  S9 --> S10
  S9 --> S11[11 Fitting kinetics]
  S11 --> S12[12 Fitting kij]
  S4 --> S13[13 DWSIM references]
  S7 --> S13
  S10 --> S13
  S12 --> S13
  S13 --> S14[14 Speed]
  S14 --> S15[15 Release v0.4.0]
  S15 --> S16[16 Docs, exit review]
```

Steps 1 to 8 meet the roadmap's exit, the study faster than DWSIM's, before
the optimiser begins.

### 0. Spikes, and `test_memory.py` — done: go, with D6 and D18 amended
- **DWSIM:**
  - can its sensitivity analysis be built and run headless, and how long
    does the owner's study take in it;
  - its Adjust and simultaneous adjusts on a spec of the reference
    flowsheet: their iterations, their whole-flowsheet solves, and their
    answers at their default and tight tolerances;
  - its optimiser on a bounded problem: which methods run headless, and
    what they stop at;
  - its regression plugin: whether it fits kᵢⱼ headless;
  - what its study does with a failed point.

  Written up as `docs/spike-dwsim-study.md`.
- **Warm starts by hand:** the reference flowsheet's recycle swept by its
  estimates, to measure the passes a warm start saves, with and without a
  predictor. This sets criterion 4's expectations.
- **Workers:** each thread's answers against one thread's, bit for bit, on
  the owner's case, to show that sharing a `Loaded` is safe. Then a profile
  of a solve's allocations.
- **`test_memory.py`,** flaky on Ubuntu since M10a: compare the bytes live
  after a full collection, not the collector's heap size. `reflux_collect`
  reports what it measures, so the header, `docs/python.md` and the next
  CHANGELOG entry say so.

**What was built:**
- **[spike/dwsim_study/run.py](../spike/dwsim_study/run.py)** and its
  write-up, [spike-dwsim-study.md](spike-dwsim-study.md). None of DWSIM's
  study, single Adjust, optimiser or regression has an API: each is a form
  of the classic app, and its work is its Run button's handler. The app's
  first automation class makes its main window without showing it, so the
  script fills each form's own case and calls the handler: what ran is
  DWSIM's code.
- **[spike/study/spike.v](../spike/study/spike.v):** warm starts by hand,
  threads' bits, what a solve allocates, and a checksum of every number of
  four studies, to compare two builds of reflux by.
  [spike/study/probe/reactors_alloc.v](../spike/study/probe/reactors_alloc.v)
  is a test of `reactors` kept there, to copy in and run: where a PFR's
  solve allocates.
- **`reflux_collect`** ([capi/capi.v](../capi/capi.v)) reports
  `gc_memory_use()`, the heap less its free blocks, after its collection,
  and its parameter is `live_bytes`. `reflux.collect()` returns it, and
  `docs/python.md` has it, with an example its test runs.

**What DWSIM does** (the write-up has each number):
- **The study** ran by its own form: 50 points of the owner's case in 727 s,
  and every row the answer of a loop of solves, to the bit (639 s). So
  criterion 11 has its like-for-like comparison.
- **A failed point records the last point's numbers,** or the solve
  before the study's if it's the first, with "Run completed" and no
  message.
- **A single Adjust isn't run by the flowsheet solver,** only by its control
  panel, which also runs unshown: a secant or Brent's method around whole
  flowsheet solves, 6 to 10 of them for the spec here (the cooler's outlet
  T until the first drum's vapour is 60 mol/s: 219.7546705409 K). With the
  recycle at DWSIM's default tolerances, Brent's ends 0.051 K off and says
  it succeeded.
- **The simultaneous adjust solver** takes 19 flowsheet solves at its
  default tolerance, and at 1e-10 leaves the flowsheet at its derivative's
  last point, the varied value 1% low and the target missed by 1.3 mol/s
  in 60, with no error.
- **The optimiser:** bounded Brent's is the one method of its 11 that finds
  an optimum, inside the bounds or on one.
  - An objective that's an expression of other variables fails in the
    nine that aren't Brent's, on an index of −1 in DWSIM's objective
    function.
  - Unbounded Brent sets the variable to 0, and the L-BFGS pair stop after
    one evaluation of one variable.
  - The DotNumerics methods and IPOPT end at the same point asked to
    maximise as asked to minimise.
  - It has "constraint" variables, a barrier added to the objective, which
    the research above didn't find.
- **The regression** ran and fitted nothing: on nine exact P-x-y points
  from PR with kᵢⱼ = 0.13, none of its four methods reported 0.13. Whether
  the app does better isn't known; the owner can try the same nine points
  there.
- **"Background parallel" is serial:** three PFRs take 2.6 to 3.1 times
  one, in every mode.

**Warm starts by hand** (the reference flowsheet, its cooler swept down
from 230 K, each point cold, then with its tear started from the last
point's answer through `SolveParams.estimates`, then from 2b − a of the
last two; passes a point, and the most any point took):

| method | points | cold | warm | predicted | warm against cold |
|---|---|---|---|---|---|
| Wegstein, tol 1e-9 | 21, 1 K apart | 4.0 | 11.0 (11) | 8.4 (9) | 6.2e-11 |
| direct | 21, 1 K apart | 29.0 | 24.0 (24) | 19.2 (20) | 1.7e-9 |
| Broyden | 21, 1 K apart | 8.1 | 11.8 (13) | 7.8 (9) | 3.0e-9 |
| Wegstein | 21, 0.1 K apart | 4.0 | 9.7 (10) | 5.4 (7) | 6.7e-11 |
| Wegstein | 9, 5 K apart | 4.0 | 12.1 (15) | 10.0 (11) | 4.7e-11 |
| Broyden | 9, 5 K apart | 7.7 | 11.1 (13) | 10.4 (12) | 1.2e-9 |
| Wegstein, tol 1e-12 | 21, 1 K apart | 5.0 | 15.5 (16) | 12.5 (13) | 8.2e-14 |

- **Every warm answer is the cold one** within the solver's tolerance, as
  criterion 4 asks.
- **But a tear's value alone doesn't save passes:** it costs them with
  Wegstein and Broyden, even 0.1 K from the last answer, and saves a sixth
  with direct substitution. The predictor helps each, and only brings
  Broyden level with cold.
- My reading, not yet checked: from an empty recycle every tear variable
  follows one geometric sequence, half the liquid coming back each pass,
  for which Wegstein's secant is exact, so cold converges in one step of
  it. A neighbour's answer is close but off in several directions at once,
  which a secant per variable handles badly, and Broyden starts each solve
  with no Jacobian, however near it starts.
- **So D6 is amended:** step 5 first measures what a warm start must carry,
  Broyden's Jacobian with the tears above all, and a block is started warm
  only where that takes no more passes than cold. Criterion 4's "can't
  grow" then has a floor: never more passes than cold.
- The owner's case has no tear, so its warm start is its units' own
  (step 5).

**Workers** (the owner's study, 48 feed temperatures, gcc `-prod`, four
cores):
- **Every bit is one thread's:** 218,200 numbers (each product's flows, T,
  P and H, each duty and every profile point) from 2, 4 and 8 threads
  sharing one `Loaded`, three runs each, and one thread run three times.
- **Throughput,** over three runs after step 1: 2 threads 1.25 to 1.33
  times one, 4 threads 1.55 to 1.62, 8 threads 1.33 to 1.41. Less than the
  2.1 this plan measured for 4, and the solves are shorter now: the
  collector's lock is a larger share.
- **What a solve allocates:** 9.7 MB for either reforming case, 20.6 kB a
  step over 471 steps. All but a few bytes of it is the phase's state,
  `Source.state`: 1,025 bytes a call (the volume root's `State` and a copy
  of the flows), about 20 calls a step. The rates, their derivatives and,
  since step 1, the viscosity allocate nothing. Step 7 starts there.

**`test_memory.py`:**
- What the library holds, by `reflux.collect()`, over ten rounds of 10,000
  calls: 2,836 KiB at every round in one run, 2,844 to 2,848 in another,
  948 to 1,016 in a third. The second half's most was never above the
  first's. The test allows 256 KiB.
- The Python suite is 83 tests, in four and a half minutes here.
- **For 0.4.0's CHANGELOG entry:** `reflux_collect` and `reflux.collect()`
  report the bytes live after the collection, not the heap's size. The
  ABI's version stays 1: the function's types are the same.

**Probed in V 0.5.2:** a value held as one interface can be tested for
another with `is`, under tcc and gcc, which D6's `Warm` needs. Inside the
branch it's the second interface only ([v-for-go.md](v-for-go.md), §5).
And a struct update on a smart-cast value copied nothing again, in the
spike's own code: an adiabatic copy of a PFR was all zeros until it was
cast with `as` first.

**Planted mistakes,** in a scratch copy:
- `reflux_collect` reporting the heap's size, or its free bytes: each
  caught by `test_collect_gives_the_live_bytes`.
- `reflux_collect` without its collection survived at first. The test now
  holds the collector's bytes since its last collection under 64 KiB after
  a call, and catches it.
- A freed fluid kept in the registry, built into the library:
  `test_memory.py`'s count climbs 7 MB a round, 36 MB over its 256 KiB.

### 1. Keeping K(T) and the viscosity — done, with D11 amended
- D11, in `kinetics` and `reactors`. The golden reports are unchanged to the
  bit.
- The published model timed before and after.
- Eq. 16's integral measured: whether a faster form holds tdep's 1e-13 and
  pays for an adiabatic reactor. Decided here, with the owner.

**What was built:**
- **[kinetics/evaluate.v](../kinetics/evaluate.v):** a `Point` keeps its K
  and d ln K/dT with the T they're at, and `fill` computes them only when T
  has moved. At a new T it integrates each species' Cp once and shares the
  two integrals among the reactions it's in. A fill that fails keeps
  nothing.
- **[kinetics/standard.v](../kinetics/standard.v):** `Reaction.standard` in
  three parts, the refusal of a T, a species' H° and S°, and ΔG° and K from
  the sums, so that `fill` sums the same numbers in the same order.
- **[reactors/viscosity.v](../reactors/viscosity.v):** `WilkeAt`, what
  Wilke's rule needs at a T whatever the composition: each compound's
  viscosity and each pair's two factors. `Wilke.at` fills one unless it
  holds that T already, and `Wilke.mixed` sums it at a composition.
  `viscosity` is the two together.
- **[reactors/pfr.v](../reactors/pfr.v)** keeps a `WilkeAt` for a reactor
  whose viscosity is Wilke's. Any other viscosity is asked whole.
- **[spike/eq16/spike.v](../spike/eq16/spike.v):** what K(T) costs, and
  equation 16 by a fixed rule.

**D11 is amended, in two ways:**
- **Wilke's pair factors are kept too.** With only each compound's
  viscosity kept, the owner's case barely moved: the rule's 25 powers and
  50 square roots a call, which depend on T alone, were most of its cost.
- **Each species is integrated once at a new T.** The reforming model's
  three reactions took 24 integrals, two for each species of each, where
  its five species need 10. That's for a T that moves, which keeping can't
  help.

**Nothing moved** (criterion 10):
- **Every bit, against the code before:** the spike built from the last
  commit and from this one gives the same checksum for the owner's study
  (218,200 numbers), the published model's, and each of them adiabatic,
  where T moves at every state: four checksums, over 48 feed temperatures
  each.
- `v run tools/golden.vsh` leaves every golden report as it was.
- In the tests: a kept K against a new point's and against
  `Reaction.standard`'s, at nine states; the viscosity kept and mixed
  against the whole, 148 mixtures; and a PFR keeping Wilke's values against
  one asking the viscosity whole, isothermal and adiabatic, 472 and 432
  profile points.

**Times** (gcc `-prod`, three runs of each build, taking turns; the
laptop's clock moves between runs):

| | before | after | |
|---|---|---|---|
| the published model, a solve | 407 to 467 ms | 18 to 23 ms | 20× |
| the published model, 48 points | 24.4 to 31.4 s | 1.20 to 1.50 s | 17 to 22× |
| the owner's case, 48 points | 1.58 to 2.05 s | 1.21 to 1.47 s | 1.2 to 1.4× |
| the published model, adiabatic, 48 points | 31.9 to 39.9 s | 17.7 to 21.0 s | 1.5 to 2.0× |
| the owner's case, adiabatic, 48 points | 13.4 to 16.9 s | 14.2 to 16.7 s | none |

- M9 measured the published model at 393 ms, 93% of it K(T). It now costs
  what the owner's typed rates do.
- A solve allocates 9.7 MB, where the published model's took 78.3: the
  quadrature's panels.
- An adiabatic solve of the owner's case takes about 300 ms against an
  isothermal one's 25, with no K in it. Not profiled yet: its Jacobian's T
  rows are central differences of the whole right-hand side, the likely
  place. That's step 14's.

**Equation 16's integral** (`spike/eq16`, gcc `-prod`):
- **Now:** a fill at a new T takes 16.8 µs for the reforming model, and
  0.05 µs at an unmoved one. One reaction's standard properties take
  11.4 µs, so the three took about 34 µs before the species were shared.
- **Where it goes:** each equation-16 species' two integrals take 3.4 to
  4.1 µs, 15 to 75 evaluations of the integrand each by the adaptive
  Gauss–Kronrod, and 540 to 740 bytes of panels. Water's closed forms take
  0.27 µs.
- **A fixed Gauss–Legendre rule from 298.15 K,** both integrals in one
  pass: over the database's 521 equation-16 compounds, each to nine
  temperatures from 200 to 1500 K, against 32 nodes on each of 16 pieces:

  | nodes | compounds within 1e-13 at every T | the worst |
  |---|---|---|
  | 16 | 519 | 2.9e-7 |
  | 20 | 520 | 1.4e-10 |
  | 24 | 521 | 2.1e-14 |
  | 32 | 521 | 4.7e-15 |

  tdep's adaptive integrals are within 3.6e-15 of the same reference. The
  worst compound at every size is para-hydrogen.
- **What it would buy:** a 16-node pass takes 0.25 to 0.41 µs a species, so
  24 nodes about 0.5: a fill at a new T in about 2 µs, eight times less,
  with no allocation.
- **What it would cost:** every answer through equation 16 would move in
  its last digits, a few parts in 1e15, so the bits of every flash fixture
  and golden report would have to be looked at again, or the rule kept to
  `kinetics` alone and two integrals of one Cp kept in the code.
- **Not decided yet:** it's the owner's call. My recommendation is to leave
  it: the owner's cases are isothermal, where K is now free, and an
  adiabatic PFR spends most of its time elsewhere (above). If step 14 finds
  K(T) in an adiabatic reactor's profile, the 24-node rule is measured and
  ready.

**Tests,** five new or extended:
- `test_k_is_kept_until_t_moves` and
  `test_a_fill_integrates_each_species_once`
  ([kinetics/law_test.v](../kinetics/law_test.v)): the bits above; a K
  overwritten by hand left at an unmoved T and replaced at another; a fill
  at −1 K refused twice and nothing kept; and a fill at a new T allocating
  2,752 bytes where its reactions one by one allocate 6,048, held under
  0.6 of it, and nothing at an unmoved T.
- `test_kept_and_mixed_is_the_whole` and `test_refused`
  ([reactors/viscosity_test.v](../reactors/viscosity_test.v)): the same for
  `WilkeAt`.
- `test_the_owners_pfr_against_dwsim`
  ([reactors/dwsim_test.v](../reactors/dwsim_test.v)): the two PFRs, kept
  and whole.

**Planted mistakes:** 22, in a scratch copy, one at a time.
- **Caught, 18:** a K or a viscosity kept when T moves, or never kept (the
  overwritten value); a failed fill still marked kept; the species'
  integrals not cleared at a new T, taken at another T, or S° summed from
  H°; ln K with R typed as 8.314; the pair factors' indices or molar
  masses swapped; the PFR keeping the viscosity at its inlet's T, with no
  room for it, or a rounding off.
- **Two survived at first:** a fill integrating every species for every
  reaction again, by either of two mistakes, changes no number.
  `test_a_fill_integrates_each_species_once` is for them, and catches
  both.
- **Two survive, explained:** S° off by one part in 1e15, and Wilke's sum
  regrouped as y·(a·a). Each moves a last digit, below every reference's
  tolerance (1e-13 for the viscosity), and the code they'd be compared with
  bit for bit is the code before this step, which the spike's checksums
  were compared with once.

### 2. Flowsheet variables in expressions
- D3: the grammar's dotted names, a context over a solution's streams and
  units, dimensions, refusals and renderings.
- Evaluating an output at a solution, with D9's NaN and reason.

### 3. Remaking a feed or unit with a value replaced
- D4: a key into the file's tree, a value with its unit put in its place,
  and the feed or unit remade by its maker, its checks included.
- What remaking costs, measured on every unit type.

### 4. Specs
- D5: `[spec.NAME]` read and checked; the planning edge; Broyden over tears
  and spec values; Brent for one spec alone; the nested retry; bounds and the
  reachable range.
- Closed forms (criterion 1) and the reference flowsheet with a spec.

### 5. Warm starts
- D6: the tears' predictor, `Warm`, and the CSTR, the reactors with an
  energy balance, and the PFR using it.
- Every warm point against its cold solve.

### 6. Studies and the CLI
- D7, D9, D16 and D17: `[study.NAME]`, chunks, each point's outcome, the
  table, JSON and CSV, and `reflux study`.
- Tutorial 2's methanol conversion as a study output: NaN, with its reason.

### 7. Workers
- D8: chunks on threads, bit for bit for any number of workers, and the
  throughput recorded. Allocations in the hot path removed where step 0's
  profile finds them.

### 8. Sensitivities
- D10: central differences at a tight tolerance, Richardson's error
  estimate, and closed forms against mpmath.

### 9. The SQP
- D12: `numeric.sqp`, the damped BFGS update, the active-set QP, the merit
  function, and the first-order conditions checked at the answer.
- Hock and Schittkowski's problems from a uv script that writes their
  published optima and starts, with each problem written again in the test.

### 10. Optimising a flowsheet
- D13: `[optimize.NAME]`, inequalities in the equation reader, the feasible
  path with warm starts, and `reflux optimize`.
- The series reactions' closed-form optimum against mpmath, an optimum on a
  bound and one on a constraint.

### 11. Fitting kinetics
- D14: `numeric.least_squares`, fitted constants as the tape's inputs, the
  centring, the data's CSV with units, `[fit.NAME]` and `reflux fit`.
- Against SciPy's fits of the same data, and the intervals' coverage over
  seeded synthetic data.

### 12. Fitting kᵢⱼ
- D15: VLE data, bubble points as residuals, and the same report.
- Against SciPy with `thermo`'s PR on the same data.

### 13. DWSIM references
- D18: the study, the spec, the optimisation and, if it runs, the kᵢⱼ fit
  against DWSIM 9, every difference explained.

### 14. Speed
- The owner's study warm, cold and on workers; the published model with K
  kept; an optimisation's solves; against DWSIM, in `docs/benchmarks.md`.

### 15. Release v0.4.0
- The CHANGELOG entry, the versions, and the release gate, as for v0.3.0.

### 16. Docs and exit review
- README, roadmap (D1's amendments of B6, H1 and H4), design notes,
  kinetics.md's fitting, the file and report formats, the V notes, and a
  criterion-by-criterion review.

## Risks

| Risk | Effect | Fallback |
|---|---|---|
| A unit hides state that threads would share | wrong answers in parallel, or a crash | step 0 checks each thread's bits against one thread's on the owner's case; `Unit`'s contract already forbids it, and every unit type is checked in step 7 |
| The collector's lock caps threads near 2× | parallel studies little faster | fewer allocations in the hot path (step 7); processes, which step 0 measured at 2.6×, if threads stay far behind |
| A warm start lands on another steady state (an adiabatic CSTR has up to three) | a smooth-looking sweep with a jump in it, or the wrong branch | every chunk starts cold, so a jump shows at a chunk's edge; a point whose warm and cold answers disagree beyond tolerance is reported, not hidden |
| Broyden over tears and specs fails far from the answer | a spec that DWSIM's nested loop would meet | D5's nested retry: Brent around the block's tears, within the bounds |
| Gradients by differences are too rough for the SQP near an optimum | slow or stalled convergence | D10's tight solves give about 1e-6; the first-order conditions are checked to what the gradients allow, and a stall is reported with how far it got |
| An SQP written from scratch has bugs a library wouldn't | a wrong optimum | Hock and Schittkowski's problems, the conditions checked at every answer, and planted mistakes, as every solver here has had |
| A fit's parameters are too correlated to estimate | huge intervals, or no convergence | centred Arrhenius forms; the correlation matrix and each interval reported, so the user sees it; a singular JᵀJ refused with the parameters it can't separate |
| DWSIM's sensitivity analysis, optimiser or regression can't be driven headless | no like-for-like comparison | M9 step 11's loop of solves for the study; for the rest, the comparison is recorded as not possible, with why |

## Not in M8

- Exact sensitivities by the implicit function theorem, and the infeasible
  path: with M13's residuals (D10, D13).
- IPOPT, for large sparse problems, and data reconciliation: M13 (D1).
- Fitting to an integral reactor's outlets: later, on M8's optimiser and a
  PFR per evaluation (D14).
- Monte Carlo uncertainty propagation (H5): after M8.
- Studies, optimisation and fits from Python: M10b, over D9's per-point
  outcomes.
- Faster PH for 10 components, `flash.uv` on cold liquids, and V's code
  running slower inside the library: still carried, not M8's.
