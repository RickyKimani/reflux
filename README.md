# reflux

An open-source steady-state chemical process simulator written in
[V](https://vlang.io). The long-term goal is to be faster and nicer to use
than DWSIM. The short-term goal is to learn V by building it.

```bash
v run .               # n-butane at 350 K and 9.4573 bar, plus PRSV and mixture examples
```

```bash
v run . 300 2.5       # T in K, P in bar
```

```bash
v -W test .           # all tests, warnings treated as errors
```

- [docs/roadmap.md](docs/roadmap.md) is what gets built, in what order, and
  how each piece is tested.
- [docs/design.md](docs/design.md) explains how the thermodynamics is built,
  and why.
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
| Scripting is awkward, and you can't bring NumPy, pandas or Jupyter to your own analysis | A Python package over reflux's C API, available from the first flash onwards. NumPy arrays go in and come out without copying, calls are vectorised over many states, and flowsheets are Python objects. The CLI also writes JSON and CSV that pandas reads directly. |

## Layout

```
v.mod               module manifest (Go's go.mod)
main.v              demo CLI
numeric/            plain numerics: cubic roots, small vector helpers
thermo/
  model.v           the Model contract F(T, V, n); volume solver; P, Z, ln φ
  cubic.v           cubic EOS (vdW, RK, SRK, PR) with vdW mixing and kᵢⱼ
  alpha.v           per-component α functions (Soave, PRSV, …)
  saturation.v      pure-component Psat, for any Model
docs/               roadmap, design notes and V notes
tools/              check.vsh (the CI checks) and, later, fixture generators
testdata/           test fixtures, each value with its source
.v-version          the V commit CI builds
```

Units: **SI everywhere inside** (K, Pa, m³, mol, J). Conversion happens only
at the edges, i.e. input files and printed output.

## Status

- [x] Model contract built on the residual Helmholtz energy F(T, V, n)
- [x] Generic volume solver that works for non-cubic models; checked against
      cubic polynomial roots
- [x] Cubic EOS for mixtures (vdW, RK, SRK, PR), van der Waals mixing with kᵢⱼ
- [x] Per-component α functions: Soave, PRSV
- [x] Fugacity coefficients, pure-component Psat
- [x] Derivatives checked against finite differences, Euler homogeneity and
      Σnᵢ ln φᵢ = Gʳ/RT
- [x] Reproduces Smith, Van Ness & Abbott Example 3.9 for all four EOS

## Roadmap

The full roadmap is in [docs/roadmap.md](docs/roadmap.md). It covers
milestones with exit criteria, headless testing, C numeric backends and
how they get replaced, and CAPE-OPEN / FMI interop. The critical path:

**M0** foundations (git, CI, gcc) → **M1** thermo core (second derivatives,
H, S, Cp) → **M2** component data (ChemSep) → **M3** flash →
**M4** first flowsheet, compared with DWSIM 9. As soon as M3 lands,
**M10a** makes thermo and flash callable from Python with NumPy arrays.

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

## Property data

The plan is ChemSep's pure-component database as the main source. It's what
DWSIM uses, so results can be compared number for number. The `chemicals` /
`thermo` Python libraries and CoolProp serve as references to check against
in tests, and NIST data is for spot checks. Commercial databanks (DIPPR,
DDB) can't be redistributed.

## Exercises

These are left for you to do while learning V. Each is small.

1. **Add the Twu (1991) α.** Add `Twu{l f64, m f64, n f64}` to the `Alpha`
   sum type, with α = Tr^(N(M−1))·exp(L(1 − Tr^(NM))). The compiler will
   list every `match` that needs a branch, and the finite-difference test
   pattern in `thermo/cubic_test.v` will check your dα/dTr.
2. **Residual enthalpy and entropy.** Write `fn residual_hs(m Model, s State)`
   using Hʳ = −RT²F_T + PV − nRT and Sʳ = −RT F_T − RF + nR ln Z. To test it,
   use the fact that for a pure component ∂(ln φ)/∂T at constant P is
   −Hʳ/(RT²).
3. **Volume translation as a wrapper.** Write a struct `Translated` with
   fields `inner Model` and `c []f64`, and make it implement `Model`.
   With C = Σnᵢcᵢ: F′(T, V, n) = F(T, V + C, n) + n·ln(V/(V + C)). This is
   an interface that holds an interface, and none of the solvers will
   notice the difference.
4. **Benchmark.** Time one million `solve_tp` calls under tcc, then under
   `v -prod -cc gcc`. Compare the same loop written in Go.
