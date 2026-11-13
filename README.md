# reflux

An open-source steady-state chemical process simulator written in
[V](https://vlang.io). The long-term goal is to be faster and nicer to use
than DWSIM. The short-term goal is to learn V by building it.

```bash
v run . solve testdata/flowsheets/reference.toml --json results.json   # a flowsheet with a recycle: its stream table, and its results as JSON
```

```bash
v run . kinetics check testdata/kinetics/xu_froment.toml   # a file's reactions: every check (docs/flowsheet-format.md)
```

```bash
v run . demo          # n-butane at 350 K and 9.4573 bar, Mathias–Copeman water, a mixture, H, S, Cp, w, μ_JT, a flash
```

```bash
v run . demo 300 2.5  # T in K, P in bar
```

```bash
v -W test .           # all tests, warnings treated as errors
```

```bash
v run tools/python.vsh build   # the C API's library, for reflux from Python: docs/python.md
```

- [docs/roadmap.md](docs/roadmap.md) is what gets built, in what order, and
  how each piece is tested. [docs/m1-plan.md](docs/m1-plan.md) is M1's
  plan, step by step, with its exit review at the end. The same goes for
  [docs/m2-plan.md](docs/m2-plan.md), component data,
  [docs/m3-plan.md](docs/m3-plan.md), the flash, and
  [docs/m4-plan.md](docs/m4-plan.md), the first flowsheet, and
  [docs/m10a-plan.md](docs/m10a-plan.md), thermo and the flash from Python.
  [docs/m9-plan.md](docs/m9-plan.md) is M9's, in progress: equations,
  kinetics and reactors.
- [docs/python.md](docs/python.md) is reflux from Python: NumPy arrays in
  and out, many states in a call. Every example on it runs in CI.
- [docs/flash-robustness.md](docs/flash-robustness.md) is the flash on
  10,000 random mixtures from the database, every answer checked, and
  [docs/benchmarks.md](docs/benchmarks.md) is reflux timed against
  `thermo` and DWSIM 9 on the same flashes, and against DWSIM on the
  reference flowsheet.
- [CHANGELOG.md](CHANGELOG.md) is every release, and what a version
  promises.
- [docs/design.md](docs/design.md) explains how the thermodynamics is built,
  and why, with every formula the code uses.
- [docs/kinetics.md](docs/kinetics.md) is the design for reaction kinetics:
  LHHW and friends as checkable data, not text.
- [docs/spike-python-ffi.md](docs/spike-python-ffi.md) is what we learned
  calling reflux from Python and NumPy, and the rules M10a follows.
- [docs/flowsheet-format.md](docs/flowsheet-format.md) is the flowsheet file
  format, and [docs/report-format.md](docs/report-format.md) the results
  `reflux solve` writes: JSON and CSV.
- [docs/spike-dwsim.md](docs/spike-dwsim.md) is what we learned driving
  DWSIM 9 from Python, and how M3 compares reflux with it.
  [docs/spike-dwsim-flowsheet.md](docs/spike-dwsim-flowsheet.md) does the
  same for its flowsheets, for M4, and
  [docs/spike-dwsim-pfr.md](docs/spike-dwsim-pfr.md) for its plug-flow
  reactor, for M9.
- [docs/v-for-go.md](docs/v-for-go.md) is V for Go programmers. It points
  back into this code for every feature.

## Why another simulator

These goals come from what went wrong with DWSIM:

| DWSIM problem | What reflux does instead |
|---|---|
| The UI gets in the way | CLI and plain-text flowsheet files come first. They can be diffed, reviewed and scripted. A GUI comes later, as a separate layer over a stable core. |
| Docs are years behind the code (8.x docs for 9.x/10.x) | Docs are generated from the code (`v doc`). Every example in the docs is compiled and run by the test suite, so it can't drift from the release. |
| 10.x didn't work, so we had to go back to 9.x | A regression suite of reference cases runs on every commit. Flowsheet files carry a format version and are migrated, never silently reinterpreted. |
| Flowsheets, especially sensitivity analyses, solve slowly | Analytic derivatives throughout, Newton/Broyden recycle convergence, warm-started and parallel sensitivity cases, and later exact sensitivities from the Jacobian. |
| Every equation is a text box: LHHW rates, equilibrium constants, specifications. Nothing checks it, a rate names its species by their place in the reaction, and every adsorption term is typed in again for each reaction | Equations are data reflux understands: parsed, never evaluated as text, every name bound to the model, units checked, rendered back as a formula and differentiated exactly. Kinetics is the first user ([docs/kinetics.md](docs/kinetics.md)): adsorption terms are defined once per catalyst site and shared by every reaction on it, and `reflux kinetics check` renders each rate back and checks element balance, units and consistency with equilibrium. |
| Scripting is awkward, and you can't bring NumPy, pandas or Jupyter to your own analysis | A Python package over reflux's C API, available from the first flash onwards. NumPy arrays go in and come out without copying, calls are vectorised over many states, and flowsheets are Python objects. The CLI also writes JSON and CSV that pandas reads directly. |

## Layout

```
v.mod               module manifest (Go's go.mod), with reflux's version
CHANGELOG.md        every release, what a version promises, how to make one
main.v              the command line: reflux solve, kinetics, demo, version, help
kinetics_command.v  reflux kinetics check and explain
demo.v              reflux demo, on database components
main_test.v         the command line, run in the test's own process
numeric/
  functions.v       what the solvers take: interfaces, not closures (V never
                    frees a closure)
  cubic.v           real roots of a cubic
  roots.v           bracketed root finding: bisection, Brent's method
  lu.v              dense LU with partial pivoting
  cholesky.v        Cholesky factorisation, the positive-definiteness test
  vec.v             small vector helpers
  integrate.v       adaptive Gauss–Kronrod quadrature
  newton.v          Newton for nonlinear systems: line search or dogleg
  broyden.v         Broyden's quasi-Newton method, with no Jacobian
  minimize.v        Newton minimisation with a Levenberg–Marquardt shift
  rachford_rice.v   two-phase Rachford–Rice, negative flash, extreme K
  ode.v             initial-value problems (M9): the interfaces, Hairer's
                    error norm and first step
  dopri5.v          Dormand and Prince's 5(4) pair, with stiffness detection
  rodas4.v          RODAS4, Hairer and Wanner's L-stable Rosenbrock method,
                    on the problem's own Jacobian
tdep/
  tdep.v            temperature-dependent correlations (ChemSep and DIPPR
                    forms): value, derivative and both integrals
components/
  record.v          a compound as a component file holds it (docs/data-format.md)
  file.v, index.v   the TOML component files and index: deterministic writer,
                    strict reader
  database.v        load and look up compounds, with suggestions for misses
  build.v           records → thermo.Component, thermo.IdealGas, and
                    correlations in SI per mol
  params.v, kij.v   ChemSep's kᵢⱼ (PR, SRK) and Mathias–Copeman (PR) files,
                    and kᵢⱼ matrices with a report of where each value came from
  formula.v         a compound's formula, where two of its structure, SMILES
                    and molar mass agree (M9)
  formation.v       ΔG°f against ΔH°f and S°, with the elements' entropies
  quality/          the sweep behind docs/data-quality.md
  chemsep/          the ChemSep XML importer: strict, with every repair reported
thermo/
  model.v           the Model contract F(T, V, n), with F at fixed T and n for
                    the volume solver; P, Z, ln φ
  properties.v      P derivatives, ln φ derivatives, residual H, S, G, U, A, Cv, Cp
  cubic.v           cubic EOS (vdW, RK, SRK, PR), exact Ωa and Ωb, vdW mixing with kᵢⱼ
  alpha.v           per-component α (Soave, PRSV, Twu, Mathias–Copeman)
  translated.v      Péneloux volume shift, as a wrapper around any Model
  ideal.v           ideal-gas Cp forms (polynomial, Aly–Lee) and formation data
  fluid.v           Model + ideal gas → total H, S, G, Cp, Cv, speed of sound, μ_JT
  saturation.v      pure-component Psat, for any Model
  phase.v           phase identification parameter and liquid/vapour labels
flash/
  stability.v       Michelsen's tangent-plane stability test, for any Model
  pt.v              the PT flash: stability, then the two-phase split
  vf.v              vapour-fraction flashes: bubble and dew points, T or P given
  ph.v              PH and PS flashes, pure components included
  tv.v              TV and UV flashes, in the model's own variables (T, V, n)
  wilson.v          Wilson's K-values
  error.v           FlashError: every failure, typed, never a panic
chem/               the elements (M9): symbols and standard atomic weights,
                    formulas (CH3(CH2)2CH3) and SMILES strings read to
                    element counts, Hill's order
quantity/           numbers with units ("30 bar", "100 kmol/h") read into SI,
                    by NIST SP 811's definitions, and the dimension algebra
                    (M9): compound units such as kmol/(kg*h) and bar^0.5,
                    each with its dimension in exact rational exponents
expr/               the expression language (M9): every equation a user gives
                    parsed into a tree, its names bound, its units checked,
                    written back as plain text, Unicode and LaTeX
  parse.v           the grammar, numbers with units, names and members,
                    readings on a unit's scale, (T in degC), and a group's
                    unit, (…) [%]
  check.v           every node's dimension, and fixed powers as fractions
  render.v          plain text that reads back, Unicode and LaTeX
  expr.v            nodes, contexts, built-in constants, definitions
  eval.v, lint.v    the tree's value; numbers that look like a constant
  tape.v            expressions compiled to a tape: folded, merged, fused;
                    evaluation into the caller's buffer, and the gradient by
                    one reverse sweep
kinetics/           reactions (M9)
  species.v         a compound as reactions know it, found by name, CAS
                    number or formula
  equation.v        equations read and balanced by elements and charge, with
                    the balanced equation suggested
  dependence.v      which reactions of a set combine others
  standard.v        ΔH°, ΔS°, ΔG° and K(T) from the ideal gas
  law.v             rate laws and sites: Arrhenius, van 't Hoff, monomials,
                    power law, LHHW, Eley–Rideal, expressions, as data
  kinetics.v        rate laws to definitions and expressions, bound, checked
                    and compiled onto one tape
  evaluate.v        rates at a state, and their derivatives by the chain rule
  check.v           reflux kinetics check: notes and warnings on what's legal
  explain.v         the factor table: each rate's parts at a state
  exact.v           exact arithmetic: null spaces, whole coefficients
equilibrium/        chemical equilibrium (M9), one phase, ideal gas or a cubic
  rand.v            the RAND method in CEA's form: Newton on the
                    constraints' potentials, every species in log space
  equilibrium.v     Gibbs minimisation over the elements, and the
                    equilibrium over chosen reactions, each shown right
  consistency.v     each reversible rate held to the equilibrium from a
                    state: where it stops along its reaction
reactors/           the reactors (M9): flowsheet units with makers for the file
  reactor.v         what every reactor shares: its elements, its errors, and
                    the temperature an energy balance sets
  gibbs.v           the Gibbs reactor
  equilibrium_reactor.v  equilibrium over chosen reactions
  conversion.v      the conversion reactor: its reactions in order, each
                    conversion a number or an expression in T, and a warning
                    where the outlet is past equilibrium
  source.v          a reactor's rates as a source of each species, with
                    their derivatives
  cstr.v            the CSTR, by Newton, with continuation in the volume
                    where Newton from the inlet fails
  pfr.v             the PFR and packed bed: isothermal, adiabatic or with a
                    duty, Ergun's drop, its profile
  viscosity.v       the gas viscosity Ergun needs: Wilke's rule
  makers.v          the reactors' tables in a flowsheet file (format 2)
flowsheet/
  stream.v          material streams, each with the flash that gave its state
  unit.v            the Unit interface, and energy streams
  mixer.v           mixer: flows and enthalpy add, then a PH flash
  splitter.v        splitter: fractions of the flow, the inlet's state
  pressure.v        outlet P or ΔP, for whatever changes pressure
  valve.v           valve: isenthalpic, by a PH flash
  heater.v          heater and cooler: outlet T, ΔT, β or duty
  pump.v            pump: v·ΔP/η, liquid only
  compressor.v      compressor and expander: isentropic by PS, then η
  drum.v            flash drum: mixed inlets, saturated vapour and liquid out
  separator.v       component separator: a fraction of each component
  sheet.v           Flowsheet: units joined by named streams, and its plan
  graph.v           strongly connected components, elementary cycles, order
  tear.v            tear selection: fewest streams, fewest second breaks
  solve.v           Solution: tears converged, every balance checked
  file.v, section.v the TOML file format: load, with every field checked
  makers.v          a unit type's table in a file, to the unit
  reactions.v       a file's [definitions], [site.NAME] and [reaction.NAME],
                    to kinetics (format 2)
  write.v           a file back, in a fixed order
  report.v          a solve's report, or its failure's, as JSON and CSV
                    (docs/report-format.md)
  table.v           the stream table and the trace, for reading
  error.v           FlowsheetError: every failure, typed, never a panic
sweep/              the flash robustness sweep: seeded database mixtures,
                    every answer checked, every failure classified
benchmarks/         the benchmark harness: reflux, thermo and DWSIM timed on
                    the same flashes, and reflux and DWSIM on the reference
                    flowsheet (docs/benchmarks.md)
release/            what must hold before a version is tagged, for
                    tools/release.vsh
capi/               the C API (M10a): reflux.h, and the functions it declares,
  capi.v            each registering its thread with the collector
  registry.v        handles: ids into a registry behind a lock
  database.v        component data and fluids, their text as TOML
  props.v           a phase's properties and vapour pressure, many points
                    at once, each with its status
  flash.v           stability, and the flashes in one layout: each phase in
                    the place of its label
  gc_threads.h      registering a thread the collector doesn't know
python/             the Python package over the C API (M10a): reflux/
                    (Database, Fluid, and points.py's results on NumPy arrays,
                    over _native's ctypes), with the library built into
                    reflux/_lib/; tests/; and examples/phase_diagrams.py, M10a's
                    exit script (docs/python.md is its guide)
docs/               roadmap, milestone plans, design notes, data format,
                    data-quality report, spike write-ups, V notes, the
                    Python guide (python.md, its figure in images/)
tools/              check.vsh (the CI checks), gen_m*_fixtures.py (references),
                    gen_m3_dwsim.py, gen_m4_dwsim.py, gen_m3_bench.py and
                    gen_m4_bench.py (DWSIM 9's answers and timings, by hand on
                    Windows), golden.vsh (the golden flowsheets' results),
                    release.vsh (the release gate), python.vsh (the C API's
                    library, and the Python tests), capi_reference.v (the V
                    API's answers, for those tests to match bit for bit),
                    gen_m10a_bench.py (the cost of a call from Python),
                    gen_m9_fixtures.py and expr_render.v (expressions against
                    SymPy, and the constants against NIST's), tape_bench.v
                    (the tape against rates written in V), gen_m9_reactions.py
                    (formulas against RDKit and chemicals, and K(T) against
                    mpmath and thermo), gen_m9_kinetics.py (two papers'
                    kinetics by mpmath), gen_m9_equilibrium.py (equilibria
                    by mpmath, ideal gas and Peng–Robinson), gen_m9_odes.py
                    (stiff test problems from the Test Set for IVP Solvers,
                    checked by SciPy), gen_m9_viscosity.py (gas viscosities
                    by chemicals), gen_m9_dwsim.py (DWSIM 9's reactors, and
                    the exact answers of its equations),
                    import_chemsep.vsh, data_quality.vsh, flash_sweep.vsh and
                    benchmarks.vsh (generated data and reports), bench_case.v
                    (reflux's side of the benchmarks)
testdata/           reference fixtures from thermo, CoolProp, chemicals, mpmath,
                    RDKit and DWSIM 9, with provenance, and flowsheets/ (TOML files,
                    and the golden flowsheets' results, those with reactors in
                    golden/reactors/), and kinetics/ (reaction files: two
                    papers' and a DWSIM case's)
data/vendor/        the ChemSep 8.50 database, verbatim, with its licence
data/components/    the same, converted: one TOML file per compound
                    (generated by tools/import_chemsep.vsh)
spike/data_format/  the TOML component-file trial (M2 step 1)
spike/dwsim/        driving DWSIM 9 from Python, headless (M3 step 0)
spike/dwsim_flowsheet/  DWSIM 9's flowsheets and recycle, headless (M4 step 0)
spike/dwsim_pfr/    DWSIM 9's plug-flow reactor on the owner's case (M9 step 0)
spike/tape/         a tape of rates against the same written in V (M9 step 0)
.v-version          the V commit CI builds
```

Units: **SI everywhere inside** (K, Pa, m³, mol, J). Conversion happens only
at the edges, i.e. input files and printed output.

## Status

**M0 (foundations), M1 (thermo core), M2 (component data), M3 (flash), M4
(the first flowsheet) and M10a (thermo and the flash from Python) are done.**
The exit reviews are at the ends of [docs/m1-plan.md](docs/m1-plan.md),
[docs/m2-plan.md](docs/m2-plan.md), [docs/m3-plan.md](docs/m3-plan.md),
[docs/m4-plan.md](docs/m4-plan.md) and [docs/m10a-plan.md](docs/m10a-plan.md).
What exists now:

- [x] Model contract built on the residual Helmholtz energy F(T, V, n), with
      first and second derivatives
- [x] Generic volume solver that works for non-cubic models; checked against
      cubic polynomial roots
- [x] Cubic EOS for mixtures (vdW, RK, SRK, PR) with exact critical
      constants and van der Waals mixing with kᵢⱼ
- [x] Per-component α functions: Soave, PRSV, Twu, Mathias–Copeman
- [x] Péneloux volume shift for any model
- [x] Fugacity coefficients and their T, P and composition derivatives;
      pure-component Psat
- [x] Residual and total H, S, G, U, A, Cp, Cv, speed of sound and
      Joule–Thomson coefficient, on a formation-enthalpy reference state
- [x] Ideal-gas Cp: polynomial and Aly–Lee (DIPPR 107)
- [x] Dense LU solver for the Newton flash to come
- [x] Agrees with Caleb Bell's `thermo` to 1e-9. PR's deviations from
      CoolProp's reference equation for n-butane are recorded and can't
      grow.
- [x] Proven callable from Python and NumPy on Windows and Linux (a spike;
      the real API is M10a)

And from M2:

- [x] The official ChemSep 8.50 database, vendored verbatim with its
      licence, checksums verified in CI
- [x] A strict importer to reflux's own TOML format: 642 compounds, with
      every defect repaired or reported
      ([data/components/CHANGES.md](data/components/CHANGES.md)), never
      fixed silently
- [x] 24 correlation forms (ChemSep's and DIPPR's), each with its value,
      derivative and both integrals, checked against mpmath and `chemicals`.
      Two rare forms aren't evaluated yet: ChemSep's eq. 200 for 5 vapour
      pressures, which fall back to Antoine, and DIPPR 124 for one liquid Cp.
- [x] Lookup by name, CAS number, InChIKey or alias, with suggestions for a
      miss. A directory of your own adds to or overrides the built-in data.
- [x] `db.fluid(.pr, ['methane', 'n-butane'])`: a working fluid, with kᵢⱼ
      from ChemSep's DECHEMA tables (and a report of where each came from)
      and optional Mathias–Copeman α. It reproduces M1's n-butane results
      with ChemSep's constants.
- [x] ChemSep checked against CoolProp for the 85 fluids both have. The
      typical deviations are 0.03% in Tc, 0.3% in ideal-gas Cp and 0.8% in
      vapour pressure, and none may grow.
- [x] Every correlation swept for sanity. What fails is listed, not hidden,
      in [docs/data-quality.md](docs/data-quality.md).

And from M3, for any `thermo.Model`:

- [x] Michelsen's stability test, with Wilson's trials and a nearly pure
      phase of each component. It agrees with `thermo` and with brute-force
      scans on 513 states, near critical points included.
- [x] The PT flash. β and both compositions match `thermo`'s to 1e-8, and
      every split is verified: a phase that's itself unstable is reported
      as a third phase (`more_phases`), not returned.
- [x] Bubble and dew points with T or P given, and PH, PS, TV and UV
      flashes, pure components and narrow-boiling mixtures included.
      Each specification round-trips to its PT state to 1e-8.
- [x] Agrees with DWSIM 9 on 25 cases, with every difference traced to
      something in DWSIM itself: its rounded constants, its saturation of a
      pure component on a vapour-pressure correlation, its ideal mixing
      entropy.
- [x] 10,000 random mixtures of 2 to 10 database compounds, with zero
      unexplained failures; every answer checked from scratch
      ([docs/flash-robustness.md](docs/flash-robustness.md)).
- [x] Timed against `thermo` and DWSIM 9
      ([docs/benchmarks.md](docs/benchmarks.md)): faster than DWSIM on PT
      flashes (1.4–5.6×) and bubble points (3.9–10×), as fast on a binary's
      PH, and 2–4× slower on PH for 5 and 10 components. M4 step 9 made PH
      2–4 times faster: twice as fast as DWSIM's for the binary, a little
      faster for 5 components and a little slower for 10.
- [x] Every failure a typed error, never a panic and never an unproven
      answer.

And from M4, the first flowsheet:

- [x] `reflux solve plant.toml --json results.json --csv streams.csv`: a
      flowsheet from a TOML file with a unit on every number, every field
      checked, and a wrong one refused with its line and a suggestion
      ([docs/flowsheet-format.md](docs/flowsheet-format.md)). A stream table
      on the screen, and every stream, unit, balance and pass as JSON and CSV
      ([docs/report-format.md](docs/report-format.md)).
- [x] Units of measure read into SI by SP 811's exact definitions: 75 units,
      with Celsius, Fahrenheit and gauge pressure only ever explicit.
- [x] The mixer, splitter, valve, heater and cooler, pump, compressor,
      expander, flash drum and component separator, each closing its
      balances to 1e-10 and checked against physics that doesn't go through
      the flash: μ_JT, Cp, Clapeyron, ∫v dP.
- [x] A sequential-modular solver: Tarjan's blocks, Johnson's cycles, and the
      fewest tears found exactly (Gundersen and Hertzberg's nine graphs), then
      direct substitution, Wegstein or Broyden. No answer is returned until
      every unit and the whole flowsheet balance.
- [x] The reference flowsheet, a recycle with two drums, matches DWSIM 9 on
      DWSIM's constants, and every difference over the targets is reproduced
      from something in DWSIM: its pump's liquid density, its mixing entropy,
      its rounded constants.
- [x] It solves in 26 ms, against DWSIM's 0.56 s
      ([docs/benchmarks.md](docs/benchmarks.md)), and five golden flowsheets
      are solved again in every check.
- [x] A first release, v0.1.0, with a release gate
      ([CHANGELOG.md](CHANGELOG.md)).

And from M10a, thermo and the flash from Python
([docs/python.md](docs/python.md)):

- [x] A C API, [capi/reflux.h](capi/reflux.h): component data and fluids,
      properties, vapour pressures, stability and every flash, many points
      in a call, each with its own status. The library exports exactly the
      header's functions, on Windows and Linux.
- [x] A Python package over it, with `ctypes`: NumPy arrays in and out,
      broadcast by NumPy's rules, an input that's already contiguous
      `float64` never copied, and NaN with a status and its reason where a
      point has no answer.
- [x] The same numbers as the V API, bit for bit, checked in CI against a V
      program built with the library's compiler and flags.
- [x] Safe from many Python threads at once, with the serial answers bit for
      bit, memory flat over 100,000 calls, and seeded fuzzing of every
      function with no crash.
- [x] An exit script that draws a T–x–y diagram and a phase-fraction map
      with matplotlib. It found a PT-flash bug, now fixed.
- [x] A second release, v0.2.0, tagged on `main` after M10a's merge
      ([CHANGELOG.md](CHANGELOG.md)).

By the roadmap, next is M10b: flowsheets from Python. M5 to M9 build on M4:
steam, activity models, columns, specifications and reactors.

## Roadmap

The full roadmap is in [docs/roadmap.md](docs/roadmap.md). It covers
milestones with exit criteria, headless testing, C numeric backends and
how they get replaced, and CAPE-OPEN / FMI interop. The critical path:

**M0** foundations ✓ → **M1** thermo core ✓ → **M2** component data
(ChemSep) ✓ → **M3** flash ✓ → **M4** first flowsheet, compared with DWSIM 9
✓ → **M10a** thermo and flash from Python ✓. Next, **M10b** makes the
flowsheets callable from Python.

## Development

```bash
v run tools/check.vsh            # what CI runs, with tcc (fast)
```

```bash
v run tools/check.vsh -cc gcc    # also builds and runs an optimised -prod binary
```

CI ([.github/workflows/ci.yml](.github/workflows/ci.yml)) runs the same
script on Windows and Linux, with tcc and with gcc. It uses the V commit
pinned in [.v-version](.v-version). After `v up`, put the new commit there
so CI and your machine stay on the same V.

Releases are tags on `main`, `vX.Y.Z`, as source: build one with
`v -prod -cc gcc -o reflux .`. [CHANGELOG.md](CHANGELOG.md) has each release,
what the version numbers promise, and how a release is made:
`v run tools/release.vsh` refuses until everything holds, then prints the tag
commands.

## Property data

The plan is ChemSep's pure-component database as the main source. It's what
DWSIM uses, so results can be compared number for number. The `chemicals` /
`thermo` Python libraries and CoolProp serve as references to check against
in tests, and NIST data is for spot checks. Commercial databanks (DIPPR,
DDB) can't be redistributed.

## Exercises

These are left for you to do while learning V. Each is small. (Earlier
exercises, the Twu α and residual H and S, are now part of M1.)

1. **Boston–Mathias α.** Boston & Mathias (1980) replace the Soave form
   above Tc with an exponential that stays positive. Add it as a new `Alpha`
   variant. The compiler will list every `match` that needs a branch, and
   [thermo/alpha_test.v](thermo/alpha_test.v) shows what to check: α(1) = 1,
   α′ and α″ against finite differences, and continuity with Soave at
   Tr = 1.
2. **Second virial coefficient.** Write
   `fn second_virial(m Model, t f64) f64` for a pure component, using
   B = lim_(V→∞) V·F/n² (evaluate at a very large V). Then check something
   step 4 found: at zero pressure the Joule–Thomson coefficient tends to
   (T·dB/dT − B)/Cp⁰, not to zero.
3. **Benchmark.** Time one million `solve_tp` calls under tcc, then under
   `v -prod -cc gcc`. Compare the same loop written in Go.
