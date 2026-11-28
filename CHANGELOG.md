# Changelog

Every release of reflux, newest first. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the versions
[semantic versioning](https://semver.org/spec/v2.0.0.html), as below.

## Versions

From 0.1.0 on:

- **The flowsheet file and the report** each carry a `format` number of their
  own ([flowsheet-format.md](docs/flowsheet-format.md),
  [report-format.md](docs/report-format.md)). A change to either is a new
  format, migrated or read beside the old, and never a silent change of meaning.
- **The command line's** commands, options and exit statuses change only in a
  minor release before 1.0, and only in a major one after.
- **The V modules' API** may change in a minor release before 1.0, as semantic
  versioning allows for 0.x, and never in a patch.
- **A golden flowsheet's results** (`testdata/flowsheets/golden/`) change only
  with an entry here that names them and says why.

From 0.2.0 on:

- **The C API** (`capi/reflux.h`) carries an ABI version of its own,
  `reflux_abi_version()`. A change that would break a caller is a new ABI
  version, and the Python package refuses a library of a version it wasn't
  written for.
- **The Python package's API** changes as the V modules' does: in a minor
  release before 1.0, and never in a patch.

## Making a release

1. Put the version in `v.mod`, and give its entry below today's date.
2. Run `v run tools/release.vsh`. It refuses unless the tree is clean, the
   entry is v.mod's version and dated, the version isn't tagged and comes after
   the last release, every golden result changed since that release is named in
   the entry, and the full check passes. Then it prints the tag command.
3. Tag on `main`, and push the tag.

## [0.7.0]

Transport of mixtures, and IF97 behind the flash: the roadmap's M5b. Every
phase of every stream now reports its viscosity and thermal conductivity,
and a liquid its surface tension, from each compound's ChemSep correlations
by mixing rules, with pressure corrections for dense gases and compressed
liquids, or why it has none; water by IAPWS-IF97 goes through every flash
but TV and UV, in a flowsheet file and its Rankine cycle; and a phase's
transport comes from C and Python. The methods were chosen against measured
mixtures from NIST's ThermoML archive, whose 27 sets are recorded with each
alternative's deviations beside reflux's; DWSIM 9's transport is
reproduced from its own formulas state by state, every difference put down
to its methods. Like 0.6.0, it's a source release
(`v -prod -cc gcc -o reflux .`).

### Added

- **Transport of mixtures** (`transport`): a gas's viscosity by Wilke's rule
  times Lucas's (1980) high-pressure factor, taken against its value at
  zero pressure; its conductivity by Mason and Saxena's rule plus Stiel and
  Thodos' dense-gas residual; a liquid's viscosity by Grunberg and Nissan's
  rule with each pair's G (0 unless a file gives it) and Lucas's (1981)
  compressed-liquid factor; its conductivity by Li's rule with Missenard's
  factor; its surface tension by Winterfeld, Scriven and Davis. A compound
  above its Tc is left out of a liquid up to 10% of its moles, and a phase
  whose every compound is above its Tc takes the gas's methods. A value is
  a correlation's, never an estimate in place of a missing one.
- **Water's own transport** in a stream of water by IAPWS-95 (R12-08,
  R15-11, R1-76) or IF97 (their industrial forms).
- **Report format 5:** each phase's `viscosity`, `thermal_conductivity` and
  `surface_tension`, a number or null, and `transport_notes` (why each null,
  the compounds a liquid left out, the values beyond their correlations'
  fitted ranges); `model.transport`, each method and every pair's G with
  its source. The CSV and the printed table are unchanged.
- **IF97 behind the flash** (`water.If97Phases`): IF97's region map as a
  phase model, its saturation region 4's own (through the new optional
  `thermo.Saturates`, which the flash's paths for one compound take), its
  range IF97's (`thermo.NoState`), its H and S IF97's on reflux's reference.
- **File format 5:** `[model] equation = "IF97"` (water alone), and
  `[[transport.pairs]]`, a liquid pair's Grunberg–Nissan G as an expression
  in T. Formats 1 to 4 read as before.
- **A golden flowsheet,** `rankine_if97.json`: M5a's Rankine cycle on IF97
  (`testdata/flowsheets/rankine_if97.toml`), within 4.9e-7 J/kg of IF97's
  own functions worked through the same unit models, and DWSIM 9's Steam
  Tables cycle by state.
- **`reflux_transport`** in the C API and **`Fluid.transport`** in Python: a
  phase's three values at many points, each with a status of its own
  (`REFLUX_POINT_NO_VALUE`, new); bit for bit the V API's. The ABI version
  stays 1.
- **Speed** (docs/benchmarks.md): a phase's transport in 2.8 to 8.7 µs with
  gcc `-prod`, 3.4 to 22 times DWSIM 9's own functions and 11 to 21 times
  thermo's from Python; an IF97 PT flash in 1.3 µs and a PH flash in 26 µs,
  eight and seven times faster than IAPWS-95's.

### Changed

- **The golden flowsheets' reports** were written again, all nine
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json`, `saturated_loop.json`, `decanter.json`,
  `rankine.json`, and `reactors/reformer.json` and
  `reactors/methanol.json`), in report format 5 and for this version. With
  each phase's transport and `model.transport` taken out, each is its
  format-4 report value for value: no other number moved.
- **Wilke's rule moved** from `reactors` to `transport` (`transport.Wilke`,
  `transport.WilkeAt`; `components.wilke` builds it); the PFR's every
  answer keeps its bits.
- **A fits file reads format 5,** leaving `[transport]` to the flowsheet.
- **Water's transport range messages** give T and ρ to 6 significant
  figures, as transport's notes do, so that a report's notes read the same
  from every compiler.

### Fixed

- **IF97's backward functions** (`water.if97`, 0.6.0's): a T within
  rounding of a closed bracket's root no longer halves the bracket without
  end ("T not found", met at 25 MPa just below 273.16 K's h), and a value
  within rounding above region 2's top at 1073.15 K and above 50 MPa is
  region 2's, not region 5's search, which has no state there (it gave
  1072.9 K at 100 MPa).

## [0.6.0] — 2026-11-27

Water and steam: the roadmap's M5a. reflux now has water by IAPWS-95, the
scientific formulation, through every flash, in a flowsheet and from Python,
and by IAPWS-IF97 with its region-3 supplements as steam-table functions;
water's viscosity, thermal conductivity and surface tension by IAPWS's
releases; and a steam table on the command line, `reflux water`. Every
verification table IAPWS publishes for what it implements is reproduced to
its printed digits, and DWSIM 9's two water packages are held to it state
by state, every difference explained (its Steam Tables package fails
silently in IF97's region 3 and 5, and above 25 MPa in its PH flash). It
begins with a fix to reflux's gas constant. Like 0.5.0, it's a source
release (`v -prod -cc gcc -o reflux .`).

### Added

- **Water by IAPWS-95** (`water`): IAPWS R6-95(2018) as a `thermo.Model`,
  with its own gas constant so that its tables match, through every flash,
  in a flowsheet file (format 4's `[model] equation = "IAPWS-95"`, water
  alone) and from the C API and Python (`REFLUX_IAPWS95`,
  `Database.fluid("IAPWS-95", ["water"])`). Its model has no state below
  251.165 K.
- **IAPWS-IF97** (`water.if97`), with its region-3 supplements SR3-03 and
  SR5-05, as steam-table functions from (p, T), (p, h), (p, s), (T, x) and
  (p, x), each backward answer finished on its basic equation.
- **Water's transport** (`water`): R12-08's viscosity and R15-11's thermal
  conductivity with their critical enhancements on IAPWS-95, their
  industrial forms on IF97, and R1-76's surface tension.
- **`reflux water`**: a steam table on the command line, a state from two of
  T, P, h, s, x and ρ with units, by IAPWS-95 or IF97 (`--if97`), per kg and
  per mole on IAPWS's reference, as a table and as JSON (the water report,
  format 1).
- **A golden flowsheet**, `rankine.json`: a Rankine cycle on IAPWS-95
  (`testdata/flowsheets/rankine.toml`), held against a hand calculation from
  the `iapws` package's states to 1.8e-10 kJ/kg, and against DWSIM 9's same
  cycle on its CoolProp package to 3.7e-7 kJ/kg, its PH flash's tolerance.
- **Speed** (docs/benchmarks.md): an IAPWS-95 PT flash takes 17.5 µs and a
  PH flash 251 µs, IF97's (p, T) 0.7 µs and (p, h) 5.5 µs, with gcc
  `-prod`; 470 to 1400 times faster than DWSIM 9's water packages, and
  1.4 to 3 times CoolProp's AbstractState from Python.
- **IF97's states carry (∂ρ/∂p)_T** (`water.if97.State.drhodp`), which
  R15-11's industrial thermal conductivity needs.

### Changed

- **A pure component's PH and PS flash above 1.01 Pc** no longer looks for
  a saturation temperature, which can't exist there and took the
  vapour-fraction flash 8 ms to fail for IAPWS-95's water at 30 MPa. Every
  answer is the same, bit for bit (1,680 pure PH and PS flashes of 40
  compounds by PR and water by IAPWS-95, from 0.5 Pc to 3 Pc); the flash
  is 30 times faster there.
- **`Database.fluid` in Python** takes `kij=None` by default: ChemSep's
  kᵢⱼ for a cubic, as before, and none for IAPWS-95, which refuses one.

- **The golden flowsheets' reports** were written again, all eight
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json`, `saturated_loop.json`, `decanter.json`, and
  `reactors/reformer.json` and `reactors/methanol.json`), for the gas
  constant above and for this version. Each number moved by at most 3.2e-11
  of itself, beyond the balances' own rounding; no label, warning, pass or
  step count changed. The reformer's PFR takes its 472 steps at volumes up
  to 7.9e-8 of themselves away, and its outlet moves by 2e-11.

### Fixed

- **reflux has one gas constant,** N_A·k_B = 8.31446261815324 J/(mol·K),
  the 2019 SI's exact value, in `quantity.molar_gas_constant`.
  `thermo.r_gas` and `activity.r_gas` were 8.314462618, that value rounded,
  1.84e-11 below the R that expressions, the kinetics and the fits used,
  though their comment called it exact; a Gibbs reactor and a PFR used
  both. They are now the one constant, as `kinetics.gas_constant` and an
  expression's `R` are. An answer that depends on R moves by about that
  fraction, and no more, except where a flash refuses at the edge of
  converging: at five of the flash robustness sweep's 10,000 states, cold
  mixtures of five to ten heavy compounds, the refusal's reason moves
  between a fourth phase and a three-phase split that doesn't converge.
  The sweep has 4 states not solved, from 5, and none unexplained.
- **A phase leaving a split of three reaches 0** (`numeric.multiphase_rr`).
  The step that took a leaving phase to its bound left a residue of about
  β·ε, so the phase shrank by ε a step; where Newton's steps were long, as
  for two components in three phases, its step underflowed to 0 after
  about 20 and the flash refused, saying Q fell no further. Water and
  1-butanol a step from their heteroazeotrope's T did so once the gas
  constant above moved them there. The phase that bounds the step now goes
  to 0 exactly; no golden result and no outcome of the flash sweep moves.

## [0.5.0] — 2026-11-25

Activity models and liquid–liquid equilibrium: the roadmap's M6. reflux now
flashes the liquids a cubic describes badly, polar and hydrogen-bonding
mixtures, by Wilson, NRTL, UNIQUAC and three UNIFACs on ChemSep's
parameters, with Henry's law for gases above their critical point; and it
finds up to three phases, two liquids and a vapour, for any model, Peng–
Robinson's included. A flowsheet takes an activity model in format 4, with
a decanter to separate the liquids. NRTL and UNIQUAC parameters can be
fitted to data, checked first for Gibbs–Duhem consistency. Like 0.4.0, it's
a source release (`v -prod -cc gcc -o reflux .`).

### Added

- **Activity models** (`activity`): Wilson, NRTL, UNIQUAC, and UNIFAC in
  its original form for VLE, its LLE table (Magnussen et al.) and Dortmund's
  modification. Each gives ln γ with its T-derivatives and its composition
  derivatives exactly, from one formula over doubles and forward jets, and
  holds Gibbs–Duhem by construction; a pure component's γ is exactly 1. Each
  answers finitely or refuses, saying why.
- **ChemSep's banks and tables,** vendored byte for byte: `nrtl.ipd`,
  `uniquac.ipd`, `wilson.ipd`, `henrys_law.ipd` and the three UNIFAC
  tables with their extensions. A bank's A is in cal/mol over DECHEMA's R,
  1.98721, which it was fitted with. The UNIFAC files contradict themselves
  in eight places, repaired by their own evidence, and DDBST's published
  tables settle 26 typos; every repair and every difference from `thermo`'s
  copy is listed in [data-quality.md](docs/data-quality.md).
- **The γ–φ fluid** (`components.Database.gamma_phi`, `thermo.GammaPhi`): a
  liquid of fugacity xᵢγᵢPᵢˢᵃᵗφᵢˢᵃᵗ with the Poynting correction, each
  compound's vapour pressure and density from its component file; the
  vapour the ideal gas or a cubic. Its enthalpy, entropy and heat capacity
  come from the same fugacities, so a pure liquid's latent heat is
  Clapeyron's and a drum's duty agrees with its phase equilibrium. It
  matches `thermo`'s `GibbsExcessLiquid` to 4.7e-13 in ln φ.
- **Henry's law:** a gas above its critical point in a liquid is a solute,
  with its Henry's constant mixed over the solvents. In water, a gas
  IAPWS's guideline G7-04 covers (14 gases) takes IAPWS's form and
  constants; others take ChemSep's bank. A compound above its Tc with
  neither is refused, by name.
- **Flashes on any phase model:** `flash` takes a `thermo.PhaseModel`, a
  phase at (T, P, n), which an equation of state gives by its volume roots
  and a γ–φ fluid by its liquid and vapour. PT, PH, PS, bubble and dew
  points and the vapour-fraction flashes work on both, and
  `flash.azeotrope_t` and `azeotrope_p` find a model's azeotrope.
- **Three phases:** a PT, PH or PS flash finds up to three phases,
  stage-wise, each phase shown stable against every trial
  (`flash.max_phases`). The fractions come from Michelsen's multiphase
  Rachford–Rice (`numeric.multiphase_rr`), then Newton on G over every
  phase. A fourth phase, or three liquids, is refused with its reason. It
  matches `thermo`'s `FlashVLN` to 1.9e-13 in the fractions.
- **Phase labels:** each phase is vapour, liquid or liquid2, in order of
  mass density: an equation of state's by its phase identification
  parameter, a γ–φ fluid's by the model that made it; at most one vapour,
  and of two liquids the denser liquid2.
- **Flowsheet file format 4** ([flowsheet-format.md](docs/flowsheet-format.md)):
  `[model]` may give `liquid`, an activity model, with `vapour`, `poynting`,
  `henry`, `missing` and `[[model.pairs]]`, each pair's parameters typed as
  expressions with units (`tau_ab = "3317.08 J/mol / (R*T)"`). A pair the
  bank lacks is refused by name, or made ideal by `missing = "ideal"`, with
  a warning. Formats 1 to 3 read as before.
- **The decanter,** a flash drum with outlets for the vapour and the light
  and heavy liquids. A flash drum sends two liquids out of its liquid port,
  and warns that it does.
- **Report format 4:** `model`, how the fluid was built, each pair's
  parameters with their source (a bank's line and its DECHEMA reference,
  the file, or none listed) and any warnings; and a phase `"liquid2"`.
- **Published data and the point test:** `fit.point_test`, Fredenslund et
  al.'s test of a binary's VLE data for Gibbs–Duhem consistency. The models
  are held against data sets from NIST's ThermoML archive
  (`testdata/m6_published.json`): methanol–water and ethanol–water VLE at
  1 atm, water–1-butanol's LLE and heteroazeotrope, and a ternary LLE, each
  model's deviations recorded.
- **Fitting NRTL and UNIQUAC:** a fit's `A[a, b]`, a pair's energy in J/mol,
  and `alpha[a, b]`, on bubble points, written back as `[[model.pairs]]`.
  The fits match SciPy's with `thermo`'s models to 1.5e-6 of their standard
  errors.
- **Compared with DWSIM 9** (`testdata/m6_dwsim.json`), by state, every
  difference explained: its NRTL and UNIFACs are reflux's to 1.7e-9; its
  Poynting term loses its sign below the vapour pressure; its UNIQUAC for
  methanol–water mixes two records; its Wilson takes other volumes and R;
  its UNIFAC-LL is its VLE UNIFAC with wrong group fractions; its default
  three-phase flash is right, two of its others aren't.
- **Speed** ([benchmarks.md](docs/benchmarks.md)): an NRTL flash of
  methanol–water in 75 µs, against 1.3 ms in `thermo` and 10.5 ms in DWSIM;
  half the time of reflux's own Peng–Robinson flash.

### Changed

- **A liquid–liquid split by an equation of state is liquid and liquid2.**
  Until now the less dense phase of any split was the vapour. By the phase
  identification parameter, 2,289 of the flash robustness sweep's 3,344
  Peng–Robinson splits are two liquids, and a few near critical points are
  too. `flash.FlashResult.beta` is the vapour's fraction, 0 with none.
- **The C API refuses two liquids,** as it refuses three phases:
  `REFLUX_POINT_MORE_PHASES`, with the phases named in its reason. Its
  layout of one liquid and one vapour can't hold them; γ–φ fluids and three
  phases come to it with M10b, as ABI 2. So the Python package's natural-gas
  map has 125 points of 6,480 refused near the critical point, where it had
  answers. The ABI version is still 1.
- **The report is format 4.** It has `model`, and a phase may be
  `"liquid2"`.
- **The V API,** as 0.x allows:
  - `flash` takes a `thermo.PhaseModel` where it took a `thermo.Model`;
    `thermo.by_roots(model)` makes any `Model` one, and `tv` and `uv` still
    take a `Model`;
  - `thermo.Fluid.model` is a `PhaseModel`, and `Fluid.helmholtz()` gives
    its `Model` where it has one;
  - `thermo.PhaseLabel` has `liquid2`, and a `FlashResult` may have three
    phases;
  - `flash.Failure.more_phases` is a fourth phase, three liquids, or more
    phases than a two-phase flash can give;
  - `flowsheet.Loaded` has `model`, and `flowsheet.FitFile` the γ–φ
    fluid's setup;
  - there's a new module, `activity`.
- **The reactors and the equilibrium refuse a γ–φ fluid,** saying they work
  in one phase by an equation of state.
- **The golden flowsheets' reports** were written again, all seven
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json`, `saturated_loop.json`, and `reactors/reformer.json`
  and `reactors/methanol.json`), in report format 4 and for this version.
  Each gains `model`. No number in them moved, nor any phase's label.
  `decanter.json` is new.

### Fixed

- **PH and PS took a NaN energy balance as met.** A balance that isn't a
  number is now unmet.
- **`numeric.newton` dropped a system's reason** when it refused at its
  start; it now gives it.

### Known limitations

- **No γ–φ fluids or three phases from Python or the C API** until M10b,
  ABI 2; two liquids are refused there too.
- **A γ–φ fluid's vapour-fraction flashes, TV and UV** keep two phases:
  a bubble point whose liquid would split is refused, and TV and UV refuse a
  γ–φ fluid and three phases.
- **Reactors work in one phase,** on an equation of state.
- **No electrolytes, Gᴱ mixing rules or associating vapours** (Hayden and
  O'Connell's); a Henry's-law solute's γ* is 1.
- **The Wilson and UNIQUAC banks' R** is inferred to be DECHEMA's, as
  NRTL's is shown to be.
- **A fit of parameters the data don't determine,** NRTL's α with both
  energies on one isobaric binary, says so in its interval, or stops at its
  residuals' rounding, which differs by platform.
- **Windows and Linux,** with V 0.5.2 (eab77e7).

## [0.4.0] — 2026-11-18

Specifications, sensitivity studies, optimisation and fitting: the roadmap's
M8. A flowsheet now answers questions as well as solving once. It can vary a
unit until a target is met, sweep its inputs with each point's derivatives,
find an optimum under constraints, and fit rate constants and kᵢⱼ to data
with confidence intervals. Each answer is shown to be right, as a solve's
is. The owner's 50-point reforming study takes 0.70 s, where DWSIM 9's own
form takes 15 minutes. Like 0.3.0, it's a source release
(`v -prod -cc gcc -o reflux .`).

### Added

- **Flowsheet variables in expressions:** a stream's `product.T`, `.P`,
  `.F`, `.n[CH4]`, `.y[CH4]`, `.h` and `.vf`, and a unit's energy streams,
  `reactor.duty`. A name whose stream or unit isn't one word is quoted,
  `'let down'.T`. A name that isn't there is refused with the closest names,
  and every one has its dimension.
- **Equations,** `=`, `<=` and `>=` between two expressions of one dimension.
  An expression with no value at a state says why, never ±∞: "the feed has
  no methanol".
- **Specs,** `[spec.NAME]` ([flowsheet-format.md](docs/flowsheet-format.md)).
  A spec varies any number the file gives a feed or a unit, by its key,
  within bounds, until an equation holds.
  - The key is `unit.cooler.outlet_T` or `feed.syngas.T`, and the feed or
    unit is remade by its own maker, with its own checks.
  - A spec is met with its block's tears, by Broyden on both, or by Brent
    alone where it has none. A nested retry takes over where that fails.
  - Each spec is checked again from the solution. A target the bounds can't
    reach is refused, with the range they reach.
  - `reflux solve` meets a file's specs, and its report gives each one.
- **Studies,** `[study.NAME]`, and `reflux study FILE [NAME]`.
  - One or two inputs, over values with units, and named outputs.
  - Each point is solved, or failed with its reason, or solved with an output
    that has no value, with why. The table, `--json` and `--csv` give every
    point.
  - Each point starts from its neighbour's answer. A unit that can use its
    last answer does: the CSTR's Newton starts from its last outlet, 10
    times faster, and an energy balance brackets its T around the last.
  - `--workers N` runs the points on threads, the processors by default. The
    answers are the same bits on any number of workers.
- **Sensitivities:** a study's `sensitivities`, each point's derivatives by
  central differences of solves at 1e-12, extrapolated by Richardson's rule.
  Each derivative comes with an estimate of its error, which the tests hold
  above the true error against closed forms.
- **Optimisation,** `[optimize.NAME]`, and `reflux optimize FILE [NAME]`.
  - It minimises or maximises an expression, over inputs with bounds, subject
    to equations and inequalities over flowsheet variables.
  - The optimiser is `numeric.sqp`, Han and Powell's SQP, written in V, on
    `numeric.qp`, Goldfarb and Idnani's dual active-set QP. Every evaluation
    is a converged solve.
  - An optimum is returned only once its first-order conditions hold. It
    comes with each bound's and constraint's multiplier, as how fast the
    optimum improves as each is relaxed, and every iterate.
  - 21 of Hock and Schittkowski's problems reach their published optima.
    Series reactions in a PFR reach their closed-form optima inside the
    bounds, on a bound and on a constraint.
- **Fitting,** `[fit.NAME]` with a CSV of data, and `reflux fit FILE [NAME]`.
  - Rate constants are fitted to rate data on the tape's exact derivatives,
    an Arrhenius or van 't Hoff pair centred on the data's mean T. kᵢⱼ are
    fitted to bubble points.
  - The solver is `numeric.least_squares`, Levenberg–Marquardt on QR.
  - Each constant comes with its standard error, its 95% interval by
    Student's t (`numeric.student_t_quantile`) and the correlations, all in
    the file's own form. The fitted section is written back as file text.
  - It matches SciPy's fits of the same data, and over 400 seeded data sets
    the 95% intervals held the truth 381 and 382 times, against 380.
- **Flowsheet file format 3:** specs, studies, optimisations and fits. A
  file of reactions alone may be format 2 or 3. Formats 1 and 2 read as
  before.
- **Report format 3:** the file's `specs`, each block's `specs` and
  `spec_residuals`, and the failure `"spec_not_met"`. The study, the
  optimisation and the fit each have a report format of their own, format 1
  ([report-format.md](docs/report-format.md)).
- **Compared with DWSIM 9** (`testdata/m8_dwsim.json`), every difference
  explained by a number a test checks:
  - the reference flowsheet's spec is 2.8 mK from its Adjust's, which is
    reflux's vapour 1.6e-3 mol/s short at DWSIM's answer, over the slope;
  - a two-stage compression's least work is 1.5 Pa from its optimiser's;
  - the series reactions in DWSIM's PFR give an optimum 1.5% off, from the
    noise in its outflow.
- **Speed** ([benchmarks.md](docs/benchmarks.md)): the owner's study
  in 1.24 s cold, 1.14 s warm and 0.70 s on 2 workers, against 904 s by DWSIM's
  sensitivity analysis. The spec is 130 times faster than DWSIM's Adjust, and
  the optimisation 4.8 times faster than its optimiser.

### Changed

- **The published reforming model solves in 13.4 ms,** from 393. A rate
  keeps each reaction's K until T moves, and at a new T integrates each
  species' Cp once for all its reactions. The PFR keeps what Wilke's rule
  needs at a T, and its states' volume roots allocate nothing
  (`thermo.Cubic.work()`): a PFR solve allocated 9.72 MB, and now 0.57 MB.
  Every answer keeps its bits.
- **The report is format 3.** It has `specs` beside format 2's keys, and
  each block has `specs` and `spec_residuals`.
- **The V API,** as 0.x allows. Everything is added and nothing removed:
  - `flowsheet.Loaded` has `input`, `with`, `solve_from`, `expression` and
    `value`, and the file's `specs`, `studies` and `optimizations`, and
    `flowsheet.load_fits` reads a file's fits;
  - a `Section` notes each number its maker reads;
  - a unit may implement `Warm` and `Energies`;
  - `flowsheet.Solution` has `specs`;
  - there are new modules: `study` and `fit`.
- **`reflux_collect` gives the bytes live after its collection,** in the C
  API and as `reflux.collect()` in Python, where it gave the collector's heap
  size. The heap's size grows in large steps when the collector chooses to,
  and is never given back, so it said little about what the library holds.
  The call is the same, so the ABI version is still 1.
- **The golden flowsheets' reports** were written again, all seven
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json`, `saturated_loop.json`, and `reactors/reformer.json`
  and `reactors/methanol.json`), in report format 3 and for this version.
  Each gains an empty `specs` and its blocks empty `specs` and
  `spec_residuals`. No number in them moved.

### Fixed

- **`python/tests/test_memory.py` failed now and then on Ubuntu.** It held
  the collector's heap size, which doubled once in a step and failed it. It
  now holds the bytes live after a full collection within 256 KiB: they
  stay within 68.

### Known limitations

- **Sensitivities and gradients are by differences.** Exact sensitivities,
  from the converged Jacobian, need each unit's Jacobian, which comes with
  M13's equation-oriented solve.
- **The optimiser is for small, dense problems,** a handful of inputs, on
  the feasible path: each evaluation is a whole solve. IPOPT, for large
  sparse problems, and data reconciliation are M13's.
- **Fits take rate data,** as from a differential reactor, and VLE as
  bubble points. Outlet data from an integral reactor comes later. A kᵢⱼ fit
  isn't compared with DWSIM's regression, which fitted nothing when a script
  ran it.
- **Workers are bounded by the collector's lock,** which every allocation
  takes: the owner's study runs about 1.6 times faster on 2 or 4 threads than
  on one.
- **Reactors work in one phase,** as in 0.3.0. There are no heat exchangers,
  and no flowsheets, studies or fits from Python (M10b).
- **Two phases in a flash,** and cubic equations of state only, as in 0.3.0.
- **Windows and Linux,** with V 0.5.2 (eab77e7).

## [0.3.0] — 2026-11-15

Equations, kinetics and reactors: the roadmap's M9. Every equation a user
gives, a rate law, a rate constant, an adsorption term or a conversion, is
read as an expression reflux understands, its names bound and its units
checked, never taken as text. On it stand reaction kinetics, chemical
equilibrium and five reactors, compared with DWSIM 9. Like 0.2.0, it's a
source release (`v -prod -cc gcc -o reflux .`).

### Added

- **The expression language,** `expr` ([design.md](docs/design.md),
  "Expressions").
  - Numbers carry their units, `240.1 kJ/mol`, and every node's dimension is
    checked: a pressure added to a temperature is refused, with why.
  - Names bind to the context (T, P, each species' p, f, c and y by name,
    formula or CAS number), to a file's definitions, and to built-in
    constants at the 2019 SI's exact values: `R`, not a typed 8.314.
  - A value can be read on a unit's scale, `0.85 - 0.0018*(T in degC)`, and
    a group can take a unit, `(85 - 0.18*(T in degC)) [%]`, so a correlation
    goes in as its source prints it.
  - Each expression writes itself back as text, Unicode and LaTeX, and is
    compiled onto a tape that evaluates allocating nothing and gives every
    derivative exactly, by one reverse sweep.
- **Compound units,** in `quantity`: SI prefixes and exact rational
  exponents, `kmol*bar^0.5/(kg*h)`.
- **Chemistry,** `chem`: elements, formulas and SMILES. Each component's
  formula is the one its structure, SMILES and molar mass agree on.
- **Reactions and kinetics,** `kinetics` ([kinetics.md](docs/kinetics.md)).
  - Equations balanced by element and charge, with the balanced equation
    suggested; a set's dependence; ΔH°, ΔS°, ΔG° and K(T) from the ideal gas.
  - Rate laws as data: power law, LHHW, Eley–Rideal and the user's own
    expression, with sites and their adsorption terms, K from the Gibbs
    energies or the paper's correlation.
  - `reflux kinetics check` flags what's wrong with a file's rates (a typed
    8.314, a denominator cut short, a K that disagrees with the Gibbs
    energies), and `reflux kinetics explain` prints each rate's factors at a
    state.
- **Chemical equilibrium,** `equilibrium`: Gibbs minimisation by the RAND
  method, on the ideal gas or a cubic, and the equilibrium over chosen
  reactions; and a check of where each reversible rate stops against it.
- **Integrators,** in `numeric`: `dopri5`, and `rodas4` for stiff problems,
  each returning an answer only once every step has met its tolerance.
- **Reactors,** `reactors`: Gibbs, equilibrium, conversion, CSTR and PFR or
  packed bed, each isothermal, adiabatic or at a duty, each checked on its
  elements.
  - The PFR carries the real gas's energy balance and Ergun's pressure drop
    on the gas's viscosity (Wilke's rule), and reports its profile.
  - A conversion may be an expression in T, and a conversion reactor warns
    when its outlet is past equilibrium, with Q/K.
  - The CSTR solves by Newton, with continuation in its volume where that
    fails.
- **Flowsheet file format 2:** reactions, sites, definitions and the reactor
  units ([flowsheet-format.md](docs/flowsheet-format.md)). Format 1 files
  read as before.
- **Report format 2:** a unit's `profile` and `warnings`, balances on the
  elements where something reacts, and the failure `"reaction_failed"`
  ([report-format.md](docs/report-format.md)).
- **Compared with DWSIM 9** (`testdata/m9_dwsim.json`): the owner's reforming
  PFR is within 2.3e-10 of the exact answer of DWSIM's own equations, and
  within 1.0e-4 of DWSIM's answer at its finest sections; each other reactor
  alone, every difference taken apart.
- **Speed** ([benchmarks.md](docs/benchmarks.md)): the owner's reforming PFR
  in 28 ms against DWSIM's 10.7 s, and its 50-point study in 2.1 s against
  11.7 minutes.

### Changed

- **A single phase more than four covolumes large is vapour,** whatever its
  phase identification parameter says. Π alone called any gas above its
  Joule–Thomson inversion temperature liquid, however dilute. 17 of M3's
  1,327 reference roots, supercritical, now differ from `thermo`'s labels.
- **The report is format 2.** Every unit has `profile` (null but for a PFR)
  and `warnings` (empty for most units), beside format 1's keys.
- **The V API,** as 0.x allows:
  - `flowsheet.load` takes the makers of other modules' units (`makers:`),
    and a unit's maker reads its table through `Section`, which now has
    `context`, `measure` and `expression`;
  - `flowsheet.Outcome` and `UnitResult` have `conserved`, `profile` and
    `warnings`.
- **The golden flowsheets' reports** were written again, all five
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json` and `saturated_loop.json`), in report format 2 and
  for this version. No number in them moved. Two flowsheets with reactors
  join them, in `testdata/flowsheets/golden/reactors/`: `reformer.json` and
  `methanol.json`. `tools/golden.vsh` takes the names of those to write.

The C API and the Python package are unchanged but for their version: ABI
version 1, as in 0.2.0.

### Fixed

- **A hot gas was labelled liquid,** and a drum would have sent it out of its
  liquid port: a reformer's product at 1000 K and 0.95 bar, for one. See the
  change above.

### Known limitations

- **Reactors work in one phase.** An outlet that would split is refused, a
  second phase being M6's.
- **K from the Gibbs energies is slow in a rate,** recomputed by quadrature
  at every evaluation: the published reforming model takes 393 ms a solve,
  14 times the same case on typed K. M8 keeps it, and the viscosity Ergun's
  drop needs, at the last T.
- **ChemSep's data,** which the Gibbs energies come from: its S° mixes 1 atm
  and 1 bar, which puts reforming's K about 4% low, and 142 compounds' ΔG°f
  disagree with their ΔH°f and S° by over 1 kJ/mol
  ([data-quality.md](docs/data-quality.md)).
- **No heat exchangers,** so no PFR cooled by a stream (a later milestone);
  no design specifications, sensitivity studies or warm starts (M8); no
  flowsheets from Python (M10b).
- **Two phases in a flash,** and cubic equations of state only, as in 0.2.0.
- **Windows and Linux,** with V 0.5.2 (eab77e7).

## [0.2.0] — 2026-11-08

Thermodynamics and flashes from Python: the roadmap's M10a. A C API that any
language can load, and a Python package over it that takes and gives NumPy
arrays, a whole array of states in each call. Like 0.1.0, it's a source
release: build the library with gcc (`v run tools/python.vsh build`), then
install the package from the source tree (`pip install -e python`).

### Added

- **The C API,** `capi/reflux.h`, ABI version 1. It covers component
  databases and fluids, a phase's properties, vapour pressures, the stability
  test, and every flash: PT, PH, PS, TV, UV and at a vapour fraction, with
  bubble and dew points.
  - Each call takes many points, and each point gets its own status.
  - It's plain functions of `int32_t`, `int64_t`, `double` and pointers, and
    nothing V allocates crosses it. Handles are ids that are never reused.
  - The library exports exactly the header's functions, on Windows and Linux.
    It's safe to call from many threads at once, and seeded fuzzing of every
    function finds no input that crashes it.
- **The Python package,** `python/`, `ctypes` over the C API.
  - `reflux.Database` and `Fluid`, with `props`, `psat`, `stability` and the
    flashes.
  - Inputs broadcast by NumPy's rules, and one already contiguous `float64`
    isn't copied.
  - A point without an answer has NaN, a `Status`, and `reason(index)` in
    words.
  - Python 3.11 or later, NumPy 2. Its guide is [python.md](docs/python.md),
    and every example there runs in CI.
- **The same numbers as the V API, bit for bit.** CI checks a reference set
  that covers every function against a V program built with the library's own
  compiler and flags.
- **The exit script,** `python/examples/phase_diagrams.py`: a T–x–y diagram and
  a phase-fraction map from Python, drawn with matplotlib.
- **Speed from Python** ([benchmarks.md](docs/benchmarks.md), "From Python").
  A call costs about 0.1 to 0.2 ms more than V's own loop. Per point, the C API
  adds a few µs to properties and nothing measurable to a PT flash.

### Changed

- **`REFLUX_DATA` names a data root,** with `components/` and
  `vendor/chemsep/` in it, not the `components/` directory itself. ChemSep's
  kᵢⱼ and Mathias–Copeman constants come from that root now too, where 0.1.0
  always read them from the source tree. A flowsheet file's own `data`
  directories are read after its database's root, and take that root's kᵢⱼ.
- **The V API of `components`:**
  - `components.load(dirs)` is now `components.open(root, extra)`;
  - `load_kij` and `load_mathias_copeman_pr` take a directory;
  - `Database` has `root` and `vendor`, and `Setup` has `kind` and
    `to_toml()`.
- **The golden flowsheets' reports** were written again, all five
  (`reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json` and `saturated_loop.json`), for the version each
  carries. No number in them moved.

### Fixed

- **A PT flash that failed inside the two-phase region.** A natural gas at 265
  K and 81 bar failed. There, the split's extrapolation overshot, and the split
  slid off towards the feed and collapsed onto it. An extrapolation is now kept
  only if it lowers the fugacity residual.
  - On M3's robustness sweep of 10,000 flashes, one refusal (still three
    phases) is proved by another tm.
  - M3's reference splits take 2,248 steps, against 2,246.
- **An empty compound name** matched the 35 compounds without an InChIKey. It
  now finds none.

### Known limitations

- **Two phases.** A state with three is refused (`more_phases`). Three-phase
  flashes come in M6.
- **Cubic equations of state only.** Activity models, steam tables, CPA and
  SAFT come later.
- **No columns, reactors or heat exchangers** (M7, M9), and no design
  specifications or warm starts (M8).
- **Flowsheets from Python** come in M10b.
- **No wheels.** The library needs gcc to build, and the package installs from
  the source tree.
- **Threads are safe but don't scale far:** with 8, 1.3 to 2.3 times one
  thread's rate, since every point allocates and the collector's lock is shared.
- **V's code runs slower inside the library:** 1.3 to 1.8 times slower in a
  Python process than in a V program, for a reason not yet known (M8).
- **A build finds its data in its source tree's `data/`,** unless `REFLUX_DATA`
  or the caller names another data root.
- **Windows and Linux,** with V 0.5.2 (eab77e7). reflux builds with tcc or gcc,
  and the library with gcc.

## [0.1.0] — 2026-11-03

The first release: the roadmap's M0 to M4. A steady-state flowsheet with a
recycle, solved from a file, on cubic equations of state and ChemSep's data,
and compared with DWSIM 9. It's a source release: build it with V 0.5.2
(`v -prod -cc gcc -o reflux .`).

### Added

- **Thermodynamics (M1).** Any model as its residual Helmholtz energy F(T, V,
  n) and its derivatives, from which the volume roots, fugacities, H, S, G, U,
  A, Cp, Cv, the speed of sound, the Joule–Thomson coefficient and a pure
  component's Psat follow. Cubic equations of state (van der Waals, Redlich–
  Kwong, SRK, Peng–Robinson) with kᵢⱼ, per-component α (Soave, PRSV, Twu,
  Mathias–Copeman), and the Péneloux volume shift. They agree with Caleb Bell's
  `thermo` to 1e-9.
- **Component data (M2).** ChemSep 8.50, vendored with its licence and converted
  to a TOML file per compound, 642 of them, with every repair reported. Lookup
  by name, CAS number, InChIKey or alias, with suggestions for a miss. kᵢⱼ from
  ChemSep's banks, and files of your own that add compounds or replace them
  ([data-format.md](docs/data-format.md)).
- **The flash (M3).** Stability, PT, bubble and dew points, and PH, PS, TV and
  UV flashes, for any model, each answer verified before it's returned. 10,000
  random mixtures of 2 to 10 compounds with no unexplained failure
  ([flash-robustness.md](docs/flash-robustness.md)), and agreement with `thermo`
  and DWSIM 9, every difference from DWSIM traced to DWSIM.
- **Flowsheets (M4).** Material and energy streams, and the mixer, splitter,
  valve, heater and cooler, pump, compressor, expander, flash drum and
  component separator. A sequential-modular solver: the flowsheet split into
  blocks (Tarjan), each block's cycles listed (Johnson) and the fewest tears
  chosen exactly, then direct substitution, Wegstein or Broyden. No answer is
  returned until every unit and the whole flowsheet balance.
- **The flowsheet file, format 1** ([flowsheet-format.md](docs/flowsheet-format.md)):
  TOML with a unit on every number, every field checked, and a wrong one refused
  with its line and a suggestion for a misspelt name.
- **The command line:** `reflux solve plant.toml [--json out.json] [--csv
  streams.csv] [--trace]`, which prints a stream table, and `reflux demo`,
  `reflux version` and `reflux help`.
- **The report, format 1** ([report-format.md](docs/report-format.md)): every
  stream, unit, balance and pass as JSON, and the streams as CSV.
- **DWSIM 9.0.5 as a reference:** each unit alone and the reference flowsheet,
  on DWSIM's constants. Every deviation has a bound that can't grow, and those
  over the targets are reproduced from DWSIM itself: its rounded constants, its
  pump's liquid density and its ideal mixing entropy.
- **Golden flowsheets:** five whose results are kept, and solved again in every
  check to 1e-9: `reference.json`, `every_unit.json`, `no_recycle.json`,
  `interlocking.json` and `saturated_loop.json`.
- **Speed** ([benchmarks.md](docs/benchmarks.md)): the reference flowsheet solves
  in 26 ms, against DWSIM's 0.56 s, and a PH flash of 5 components in 4.0 ms.

### Known limitations

- **Two phases.** A state with three is refused (`more_phases`); three-phase
  flashes come in M6.
- **Cubic equations of state only.** Activity models, steam tables, CPA and
  SAFT come later.
- **No columns, reactors or heat exchangers** (M7, M9), and no design
  specifications or warm starts (M8).
- **A build finds its component data where it was built**
  (`data/components/` in the source tree), or wherever `REFLUX_DATA` points.
- **Windows and Linux,** with V 0.5.2 (eab77e7), tcc and gcc.
