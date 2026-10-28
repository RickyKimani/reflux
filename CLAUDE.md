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

- **M0, M1 and M2 are done** and merged into `main` (PRs #1 and #2).
  `docs/roadmap.md` has every milestone (M0–M14) and its exit criteria.
- **M3 (flash) is done** on branch `m3-flash`, which is pushed and has no
  PR yet. Its exit review, at the end of `docs/m3-plan.md`, finds all ten
  criteria met. The PR into `main` waits for CI on the last commit, and for
  the owner to ask.
  - The plan, `docs/m3-plan.md`, holds exit criteria 1–10, decisions
    D1–D12 and steps 0–11, each step "— done" with its write-up: the
    record of what each step did.
  - What M3 built: Newton and Rachford–Rice, PIP and the Soave guard,
    stability, the PT flash (`flash/pt.v`), the reference fixtures
    (`testdata/m3_thermo.json`, `testdata/m3_dwsim.json`), the
    vapour-fraction flashes with bubble and dew points (`flash/vf.v`), the
    PH and PS flashes (`flash/ph.v`), the TV and UV flashes (`flash/tv.v`),
    the robustness sweep (`sweep/`, `docs/flash-robustness.md`), and the
    benchmarks with a performance pass (`benchmarks/`,
    `docs/benchmarks.md`).
  - Amended along the way: D3 in step 4 (pure-component stability trials),
    D6 in step 6 (the vapour is the phase of lower mass density),
    criterion 6 after step 5 (a DWSIM deviation may exceed a target only
    when the fixture reproduces it from DWSIM itself, and reflux still
    matches the exact model to 1e-8), D8 in step 9 (a third phase, proved
    by a negative tangent-plane distance, is an explained failure: VLLE is
    M6), and in step 10 D10 (the module is `benchmarks`, since vlib has a
    `benchmark`) and M1's D2 (the volume solver's tier on `Model` is
    `volume_function(t, n)`).
  - Carried forward (the review's list): PH for 5 and 10 components takes
    2–4× as long as DWSIM's, about 12 verified PT flashes per answer; VLLE
    is M6; upstream reports to V, `thermo` and DWSIM.
- **Next is M10a** (thermo and the flash from Python, with NumPy arrays),
  then M4, the first flowsheet. Each starts with its plan,
  `docs/mN-plan.md`, approved before any code.

## How we work

- **Plan a milestone before coding it.** Research, then write
  `docs/mN-plan.md` with exit criteria, decisions (D1, D2, …) and steps. It
  needs the owner's approval before any code.
- **One step at a time.** The owner says "Commit step N and start step N+1".
  - Commit only when asked, and push only when asked.
  - Open the PR into `main` only after the milestone's docs and exit review,
    and only when asked.
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
v run .                          # the demo CLI; `v run . 300 2.5` sets T (K) and P (bar)
```

- **`check.vsh`** runs `v fmt -verify .`, `v vet .`, `v -cc <cc> -W test .`,
  then builds and runs `main.v`. With gcc it also builds and runs `-prod`.
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
    tools/benchmarks.vsh` then writes `docs/benchmarks.md`. Rerun both
    after a change that could make flashes faster or slower. `v -prod -cc
    gcc run tools/benchmarks.vsh --time` times reflux alone, against the
    fixture, without DWSIM.
  - `v -prod -cc gcc run tools/flash_sweep.vsh` writes
    `docs/flash-robustness.md`, 10,000 PT flashes (about three minutes,
    the build included; about twenty under tcc). Regenerate it after
    changing the flash or the sweep's checks.
- **Reference fixtures** come from uv scripts whose headers pin their
  dependencies: `thermo` 0.6.1, `chemicals` 1.5.2, `mpmath` 1.3.0,
  `CoolProp` 8.0.0.
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
  - `spike/dwsim/run.py` is step 0's spike.
  - DWSIM's PR uses rounded constants, saturates pure components on
    ChemSep's vapour-pressure correlation, and has a wrong ideal mixing
    entropy. "What step 5 found" in `docs/spike-dwsim.md` has the details.
    Compare with it by state, and explain each difference.

## Architecture

The modules depend in one direction: `numeric` ← `tdep` ← `thermo` ←
`components`, and `flash` uses `thermo` and `numeric` but not
`components`. `sweep`, the flash robustness sweep, is the one module that
uses both `components` and `flash`.

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
    re-verified.
  - `tvf` and `pvf` (`vf.v`) find P or T for a given β; bubble and dew
    points are β = 0 and 1. Newton runs on ln K and ln P (or ln T)
    from Wilson's start, with a PT-split start as fallback. Two answers
    are possible near critical points, and `estimate` picks between them.
  - `ph` and `ps` (`ph.v`) take a `thermo.Fluid` and find T from a molar H
    or S, by Michelsen's Q-function. T is bracketed by PT flashes on the
    sign of H(T) − H_spec, narrowed by Brent, and finished by Newton in
    (v, ln T). A pure component inside its jump at Tsat gets β from the
    energy balance.
  - `tv` and `uv` (`tv.v`) find P (and T) from a molar V (and U). The finish
    minimises A over the vapour's mole numbers and volume, in F(T, V, n)
    alone, with no volume roots (D5); UV adds ln T and the energy balance to
    that Newton.
  - In a split, the vapour is the phase of lower *mass* density (D6, as
    amended in step 6), not the larger molar volume.
  - Results are `FlashResult`; failures are `FlashError` with a
    `Failure` kind.
- **`numeric`** holds the solvers the flash rests on: LU, Cholesky, `newton`
  for systems (line search or dogleg) and `minimize` (a
  Levenberg–Marquardt-shifted Newton that knows f's rounding through
  `fscale`). Also `rachford_rice` (accurate for K from 1e-30 to 1e25),
  QAG integration and bisection. They take what they solve as an interface
  (`Function`, `FallibleFunction`, `System`, `Objective`), never a closure:
  V 0.5.2 never frees a closure, and one per flash leaked kilobytes.
- **`components`** runs the data pipeline: ChemSep XML, a strict importer,
  one TOML file per compound in `data/components/`, then the database API.
  The API yields `thermo.Component`, `IdealGas` and kᵢⱼ matrices from the
  vendored `pr.ipd`/`srk.ipd`, taking the first entry listed for a pair.
- **Units are SI inside** (K, Pa, m³, mol, J), with conversion only at the
  edges.
- **No panics reach a caller.** Every failure is a typed error
  (`FlashError`, `NewtonError`, `InputError`, `LookupError`), because the
  Python API (M10a) can't catch panics.
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
- **Module names:** not like a vlib module's. `benchmark` broke external
  tests, because V's test runner imports vlib's.
- **Leftovers:** V can leave `.*.v3cc.*/` build folders next to a binary.
  They're git-ignored; delete them.
