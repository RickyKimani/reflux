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
