# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What reflux is, and for whom

An open-source steady-state chemical process simulator in V. The aim is to
be faster and nicer to use than DWSIM. The README's "Why another simulator"
table lists the DWSIM problems that drive the design: the UI, stale docs,
DWSIM 10 not working, slow flowsheets and sensitivity runs, no NumPy in
scripting, and LHHW kinetics typed as text.

It's also how the owner learns V. He is a chemical-engineering student,
strong in thermodynamics, who knows Go. Explain V by contrast with Go, and
put lasting lessons in `docs/v-for-go.md`, checked by compiling a probe
before they're written. He works on Windows with nushell. His Go package
`C:\Users\Ricky\dev\zfactor` is a reference only: don't copy from it.

## Where things stand

Keep this section current: update it whenever a step is committed.

- **M0 to M4 are done** and merged into `main` (PRs #1 to #4).
  `docs/roadmap.md` has every milestone (M0–M14) and its exit criteria.
  Each milestone's plan, `docs/mN-plan.md`, ends with its exit review, and
  each finished step there has its write-up: the record of what it did.
- **M3, the flash,** carried forward (in its review): VLLE is M6; reports
  to make upstream to V, `thermo` and DWSIM. Its PH speed was M4 step 9's,
  which made it 2–4× faster.
- **M4, the first flowsheet,** carried forward (in its review):
  - `flash.uv` fails on cold liquids at or below 100 K, in its TV flash. No
    unit uses UV yet; fix it before one does.
  - PH for 10 components is still a little slower than DWSIM's (5.9 ms
    against 5.2).
  - Warm starts across solves (F5), and a solver callback for `--trace` as
    it happens: M8, if its PFR needs them.
  - Reports to make upstream: to V, the TOML parser's panic on a `\u`
    escape and a smart-cast interface's field inside `${}`; to `thermo`, its
    flash at a vapour fraction; to DWSIM, its two-inlet drum and its silent
    unit failures.
- **v0.1.0** is released: M0 to M4, a source release, tagged on
  2026-11-03 at `45d4ce6` on `main` (an annotated tag, pushed). The next
  release starts with a new CHANGELOG entry and its version in `v.mod`, and
  a changed golden report goes in that entry by name.
- **M10a (thermo and the flash from Python) is planned** on branch
  `m10a-python`. Its plan, `docs/m10a-plan.md`, was approved on 2026-11-04:
  exit criteria 1–11, decisions D1–D14 and steps 1–10, ending in v0.2.0.
  - Its probes found that V 0.5.2's collector has thread support on Windows
    too, so the spike's thread registration is needed on both OSes (the
    spike's write-up says otherwise, and step 1 corrects it). 8 Python threads
    at once get the serial answers bit for bit, at about 2× the serial rate.
  - D9 changes `components`: the kᵢⱼ banks and `mc-PR.lib` load from the
    source tree whatever `REFLUX_DATA` says, and must follow the data root.
  - Next is step 1: the library skeleton (`capi`, `capi/reflux.h`, the
    registry, thread registration, restricted exports, `tools/python.vsh`).
  - M10b (flowsheets from Python) follows; M5 to M9 build on M4.
- **The owner's own benchmark,** for M8 and M9: a one-PFR methane steam
  reforming flowsheet, 15 to 18 s a solve in DWSIM, in `docs/roadmap.md`.

## How we work

- **Plan a milestone before coding it.** Research, then write
  `docs/mN-plan.md` with exit criteria, decisions (D1, D2, …) and steps. It
  needs the owner's approval before any code.
- **One step at a time.** The owner says "Commit step N and start step N+1".
  - Commit only when asked, and push only when asked.
  - Open the PR into `main` only after the milestone's docs and exit review,
    and only when asked.
  - Tag a release only when asked, on `main`, once `v run tools/release.vsh`
    passes ("Making a release" in `CHANGELOG.md`). A change to a golden
    report goes in the CHANGELOG's next entry, by name.
  - Commit subjects look like `M3 step 4: PT flash`, with a body below.
- **Finishing a step** means all of:
  - marking it "— done" in the plan, with what was built, the tests,
    measured numbers, bugs found and the planted-mistake results;
  - updating `docs/design.md` (every formula the code uses), the README
    layout block, `testdata/README.md` (one row per fixture) and
    `docs/v-for-go.md` (V lessons);
  - running the tcc check, and the gcc check in the background.
- **Changing an approved decision** happens in the plan's table, marked
  "Amended in step N", and is flagged to the owner.

## Commands

```bash
v run tools/check.vsh            # what CI runs, with tcc: about 4 min
v run tools/check.vsh -cc gcc    # with gcc, plus a -prod build: about 8 min
v -W test flash/pt_test.v        # one test file
v -W -stats test flash/pt_test.v # per-test results, and println output from passing tests
v run . solve testdata/flowsheets/reference.toml --json out.json   # the CLI (docs/report-format.md)
v run . demo                     # the demo; `v run . demo 300 2.5` sets T (K) and P (bar)
v run tools/release.vsh          # the release gate: why this commit can't be tagged, or the tag commands
```

- **`check.vsh`** runs `v fmt -verify .`, `v vet .`, `v -cc <cc> -W test .`,
  then builds reflux and runs `reflux demo`, solves the reference flowsheet
  with both reports, and checks that a missing file exits 1. With gcc it
  does the same with a `-prod` build.
- **Test builds are the slow part:** each `_test.v` file is compiled into a
  program of its own. `check.vsh` runs up to 4 at once (V's runner would
  allow one per 8 GiB of RAM: one here), and builds with `-nocache` under
  gcc, where V 0.5.2's module cache made each build about three times
  slower. Run the gcc check in the background, and don't edit or regenerate
  anything it compiles or reads until it finishes.
- **CI** (`.github/workflows/ci.yml`) runs `check.vsh` on Ubuntu and Windows
  with both tcc and gcc. The gcc jobs also run the Python FFI spike. V is
  pinned to the commit in `.v-version` (0.5.2, eab77e7).
- **Generated data:**
  - `v run tools/import_chemsep.vsh` regenerates `data/components/` from the
    vendored ChemSep XML. A test fails if it's stale.
  - `v run tools/data_quality.vsh` writes `docs/data-quality.md`.
  - `tools/gen_m3_bench.py` (by hand, with DWSIM) times reflux, thermo and
    DWSIM case by case and writes `testdata/m3_bench.json`; `v run
    tools/benchmarks.vsh` then writes `docs/benchmarks.md`.
    `tools/gen_m4_bench.py` does the same for the reference flowsheet,
    into `testdata/m4_bench.json`. Rerun them after a change that could
    make flashes faster or slower. `v -prod -cc gcc run
    tools/benchmarks.vsh --time` times reflux alone, against the fixtures,
    without DWSIM.
  - `v -prod -cc gcc run tools/flash_sweep.vsh` writes
    `docs/flash-robustness.md`, 10,000 PT flashes (about three minutes,
    the build included; about twenty under tcc). Regenerate it after
    changing the flash or the sweep's checks.
  - `v run tools/golden.vsh` writes the golden flowsheets' reports
    (`testdata/flowsheets/golden/`). Run it only when a change is meant to
    move the answers, and explain why in the commit.
- **Reference fixtures** come from uv scripts whose headers pin their
  dependencies: `thermo` 0.6.1, `chemicals` 1.5.2, `mpmath` 1.3.0,
  `CoolProp` 8.0.0, and for `tools/gen_m4_graphs.py`, `networkx` 3.4.2 and
  `scipy` 1.14.1.
  - On this machine uv's managed-Python link fails in the sandbox, so pass
    the interpreter explicitly, and set `PYTHONIOENCODING=utf-8` (the
    console is cp1252):

    ```bash
    PYTHONIOENCODING=utf-8 uv run --python "C:/Users/Ricky/AppData/Roaming/uv/python/cpython-3.13.17-windows-x86_64-none/python.exe" tools/gen_m3_fixtures.py
    ```

  - `gen_m3_fixtures.py` takes the names of the fixtures to write (`thermo`,
    `stability`, `pip`, `rachford_rice`); with none, it writes them all.
  - A generator must reproduce its committed output byte for byte on the
    same day. Check with `cmp` after regenerating.
- **DWSIM 9.0.5** is installed at `%LOCALAPPDATA%\DWSIM`. pythonnet 3.0.5
  stops at Python 3.13, so pass the same interpreter as above.
  - `tools/gen_m3_dwsim.py` writes `testdata/m3_dwsim.json`. It imports
    `gen_m3_fixtures.py` for its `thermo` helpers.
  - `tools/gen_m4_dwsim.py` writes `testdata/m4_dwsim.json`, DWSIM's
    constants as component files, and `reference_dwsim.toml`. It imports
    `gen_m3_dwsim.py`, and needs V on the PATH: it solves that file with
    reflux for the open loop. About five minutes.
  - thermo 0.6.1's flash at a vapour fraction strictly between 0 and 1
    fails on a mixture (`dew` unbound): bisect on its PT flash instead.
  - `spike/dwsim/run.py` is M3 step 0's spike, and
    `spike/dwsim_flowsheet/run.py` M4 step 0's: DWSIM's flowsheets and
    recycle (`docs/spike-dwsim-flowsheet.md`).
  - DWSIM's PR uses rounded constants, saturates pure components on
    ChemSep's vapour-pressure correlation, and has a wrong ideal mixing
    entropy. "What step 5 found" in `docs/spike-dwsim.md` has the details.
    Compare with it by state, and explain each difference.
  - In DWSIM's flowsheets, a failed unit is silent: check every unit's
    `Calculated` and `ErrorMessage`. Its recycle needs its tolerances scaled
    by 1e-8 to close the balances, its pump uses a liquid density that
    isn't PR's, and its flash drum is wrong with two inlets.
  - Its separator's `EnergyImb` is what goes in less what comes out, and its
    wrong mixing entropy moves a two-phase isentropic outlet: 0.42 K in
    step 8's expander. "What step 8 found" in
    `docs/spike-dwsim-flowsheet.md` has the rest.

## Architecture

The modules depend in one direction: `numeric` ← `tdep` ← `thermo` ←
`components`, and `flash` uses `thermo` and `numeric` but not
`components`. `quantity`, units of measure (M4), depends on nothing, and
so does `release`, the release gate's logic.
`sweep`, the flash robustness sweep, is the one module that uses both
`components` and `flash`. `flowsheet` (M4) uses `flash`, `thermo` and
`numeric`, and its tests use `components`.

- **`thermo` is built on one contract.** `Model` (thermo/model.v) supplies
  the reduced residual Helmholtz energy F(T, V, n) = Aʳ/RT and its
  derivatives, with V the total volume and n the mole numbers.
  - Everything else is derived generically from F: the volume roots, P,
    ln φ and its derivatives, H, S, Cp, PIP and Psat.
  - Cubics (`Cubic`, with per-component `Alpha` functions) and the Péneloux
    wrapper `Translated` implement it. CPA and SAFT will too (M12), with no
    change to the flash.
  - The volume solver's tier is `volume_function(t, n)`: F at fixed T and n
    as a function of V, so that the parts in T and n (a cubic's A and B) are
    computed once per root, not on every Newton step.
  - `Fluid` adds the ideal-gas part, for total properties.
- **`flash` works for any `thermo.Model`** (D2).
  - `stability` is Michelsen's tangent-plane test. Its trials are Wilson's
    vapour- and liquid-like ones plus a nearly pure phase of each component.
    Substitution with DEM is followed by a second-order finish in α = 2√W.
  - `pt` tests stability, then splits: substitution with Rachford–Rice and
    DEM, then Newton on G through `numeric.minimize`, then each phase
    re-verified. Its feed's test stops at the first proof of instability
    once Wilson's two trials have run, and runs the rest only if no split
    from those stands.
  - `tvf` and `pvf` (`vf.v`) find P or T for a given β; bubble and dew
    points are β = 0 and 1. Newton runs on ln K and ln P (or ln T)
    from Wilson's start, with a PT-split start as fallback. Two answers
    are possible near critical points, and `estimate` picks between them.
  - `ph` and `ps` (`ph.v`) take a `thermo.Fluid` and find T from a molar H
    or S, by Michelsen's Q-function. T is bracketed by PT flashes on the
    sign of H(T) − H_spec, narrowed by Brent, and finished by Newton in
    (v, ln T). Only the answer's PT flash is verified. A pure component
    inside its jump at Tsat gets β from the energy balance.
  - `tv` and `uv` (`tv.v`) find P (and T) from a molar V (and U). The finish
    minimises A over the vapour's mole numbers and volume, in F(T, V, n)
    alone, with no volume roots (D5); UV adds ln T and the energy balance to
    that Newton.
  - In a split, the vapour is the phase of lower *mass* density (D6, as
    amended in step 6), not the larger molar volume.
  - Results are `FlashResult`; failures are `FlashError` with a
    `Failure` kind.
- **`numeric`** holds the solvers the flash rests on: LU, Cholesky, `newton`
  for systems (line search or dogleg), `broyden` for systems with no Jacobian,
  and `minimize` (a Levenberg–Marquardt-shifted Newton that knows f's rounding
  through `fscale`). Also `rachford_rice` (accurate for K from 1e-30 to 1e25),
  QAG integration and bisection. They take what they solve as an interface
  (`Function`, `FallibleFunction`, `System`, `Equations`, `Objective`), never
  a closure: V 0.5.2 never frees a closure, and one per flash leaked kilobytes.
- **`components`** runs the data pipeline: ChemSep XML, a strict importer,
  one TOML file per compound in `data/components/`, then the database API.
  The API yields `thermo.Component`, `IdealGas` and kᵢⱼ matrices from the
  vendored `pr.ipd`/`srk.ipd`, taking the first entry listed for a pair.
- **`flowsheet`** is M4's: streams, units, the graph, the solver, the file
  and the reports.
  - A `Stream` is its molar flows plus the flash that gave its state, so
    it's never half specified. Its H, S and V are its phases'.
  - A `Unit` is `calculate(fluid, inlets) !Outcome` and holds nothing about
    a solve. An energy rate above zero goes into the unit's material.
  - A `Flowsheet` is units with a stream named on each port, and `plan`
    checks it and returns blocks to solve in order. Tarjan gives the blocks
    and Johnson the cycles, and `tear.v` picks the fewest tears exactly, by
    branch and bound, then the fewest second breaks.
  - A file is read by `load` through a `Section` per table, and a `Maker` per
    unit type. Later unit families add makers, so `flowsheet` needn't import
    them (D1).
  - `solve` iterates each block's tears, whose variables are flows, P and
    molar H, by direct substitution, Wegstein or Broyden. A `Solution` is only
    returned with every balance checked: each unit to 1e-10 on the inlets it
    used, and the whole to the gap the tears left.
  - A `Report` is a solution, or a failure, in SI: `Loaded.report`,
    `Loaded.failure` or `load_failure`. Its JSON and CSV come from that one
    struct, and its stream table rounds by `math.round`, never printf.
    `main.v` is the command line over it, tested in process by
    `main_test.v`.
- **Units are SI inside** (K, Pa, m³, mol, J), with conversion only at the
  edges.
- **No panics reach a caller.** Every failure is a typed error
  (`FlashError`, `NewtonError`, `InputError`, `LookupError`,
  `QuantityError`, `FlowsheetError`), because the Python API (M10a) can't
  catch panics.
- **A solver returns a result only when it has shown the result to be
  right**, meaning its tolerance is met and the split verified. Anything
  else is an error, never a plausible-looking answer.

## Testing practice

- **References and tolerances:**
  - References come from independent sources (`thermo`, `chemicals`,
    CoolProp, mpmath at high precision, brute-force scans). They are stored
    as JSON in `testdata/` with provenance. Tests never call Python.
  - Measure first, then assert a tolerance justified by conditioning or
    rounding, and say why in the test's comment.
  - Model-accuracy deviations and iteration budgets are recorded as
    numbers that can't grow.
- **Check the claim, not the code path.** For example, prove an "unstable"
  verdict with a TPD computed from scratch, and check every split's
  fugacities, mass balance and phase stability independently.
- **Planted-mistake checks for every step:**
  - Mutate one snippet at a time in a scratch copy of the modules (never in
    the repo while a check is running), run the relevant tests, restore,
    and confirm the files match with `cmp`.
  - Each survivor gets a new test, gets its code deleted, or is explained
    in the plan.
- **Keep sensitive states exact.** Where a test state is sensitive to an
  ulp, store the exact double: build it from integers, for instance
  `(915 + k) * 1e4`.

## V notes that bite here

`docs/v-for-go.md` has the full set. The ones that most often break a build
or a test:

- `-W` turns warnings into errors, unused variables among them.
- **Shadowing is uneven:** a local named like a const is always an error
  ("duplicate of a const name"). One named like a function is only a
  notice, but the checks are kept free of notices, so rename it.
- **Results:** `return f()` forwards a result, but inside a `match` or `if`
  arm the call needs `!`. An Option can't be propagated with `?` from a
  function that returns `!T`.
- **Arrays:** a `mut` variable can't alias an immutable array, so use
  `.clone()`. Two `mut` variables share, as Go slices do. There's no `a + b`;
  `a << b` appends every element of b.
- **`v vet` doesn't type-check**, so compile what you change, `.vsh`
  scripts included.
- **On Windows, `import runtime` defeats V's default compiler** (tcc can't
  parse `psapi.h`). Only the fallback compiler builds such a program.
- **V reports its own failed builds** to bugs.vlang.io when its fallback
  compiler succeeds: build options and C errors, with local paths. Set
  `V_C_ERROR_BUG_REPORT_DISABLED=1` when probing V's edges.
- **Structs:** a struct field of function type needs `@[required]` or a
  default. Writing `= false` on a `bool` field is an error, because fields
  are zeroed already.
- **Strings:** `'${x}'` and `"${x}"` both interpolate, so V source built as
  text needs raw strings, `r'...'`.
- **Writing files from the shell:** bash heredocs mangle backslashes, so
  write such files with the Write tool or a Python script.
- **Formatting:** printf-style formats round ties differently under tcc
  and gcc, so round with `math.round` before formatting anything compared
  byte for byte.
- **Shared libraries** need gcc: tcc's Windows DLLs export nothing.
- **Closures are never freed**, so nothing that runs per flash may make one:
  pass a struct with a method, as `numeric`'s interfaces ask.
- **A panic can't be caught.** There's no `recover`, so what a library might
  panic on must be kept from it: `toml.parse_text` on a `\u` escape, for one.
- **Module names:** not like a vlib module's. `benchmark` broke external
  tests, because V's test runner imports vlib's.
- **Leftovers:** V can leave `.*.v3cc.*/` build folders next to a binary.
  They're git-ignored; delete them.
- **Smart App Control** can refuse V's programs, and V reports it as `failed
  CreateProcess ...: The operation completed successfully`. A refused run
  proves nothing, so a mutation loop must retry it, not count it as caught.
  It's the owner's security setting: don't work around it. See
  `docs/v-for-go.md`, §1.
- **`json2`** leaves a `none` field out unless it has `@[json_null]`, and
  writes a map's keys sorted, so anything ordered goes in an array. In
  0.5.2 `Any.arr()` is deprecated, an error under `-W`: use `as_array()`.
- **A `!` in a test function that doesn't return a Result** can turn a typed
  error into a `V panic` with no FAIL line, so a panic in a test's output
  isn't always a panic in the code.
