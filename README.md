# reflux

An open-source steady-state chemical process simulator written in
[V](https://vlang.io). The long-term goal is to be faster and nicer to use
than DWSIM. The short-term goal is to learn V by building it.

```bash
v run .               # n-butane at 350 K and 9.4573 bar, Mathias–Copeman water, a mixture, H, S, Cp, w, μ_JT
```

```bash
v run . 300 2.5       # T in K, P in bar
```

```bash
v -W test .           # all tests, warnings treated as errors
```

- [docs/roadmap.md](docs/roadmap.md) is what gets built, in what order, and
  how each piece is tested. [docs/m1-plan.md](docs/m1-plan.md) is M1's
  plan, step by step, with its exit review at the end, and
  [docs/m2-plan.md](docs/m2-plan.md) is the plan for M2, component data.
- [docs/design.md](docs/design.md) explains how the thermodynamics is built,
  and why, with every formula the code uses.
- [docs/kinetics.md](docs/kinetics.md) is the design for reaction kinetics:
  LHHW and friends as checkable data, not text.
- [docs/spike-python-ffi.md](docs/spike-python-ffi.md) is what we learned
  calling reflux from Python and NumPy, and the rules M10a follows.
- [docs/spike-dwsim.md](docs/spike-dwsim.md) is what we learned driving
  DWSIM 9 from Python, and how M3 compares reflux with it.
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
| LHHW kinetics is a text box: no way to check the equation, and every adsorption term typed in again for each reaction | Kinetics is structured data ([docs/kinetics.md](docs/kinetics.md)). Adsorption terms are defined once per catalyst site and shared by every reaction on it. `reflux kinetics check` renders each rate back as a formula and checks element balance, units and consistency with equilibrium. |
| Scripting is awkward, and you can't bring NumPy, pandas or Jupyter to your own analysis | A Python package over reflux's C API, available from the first flash onwards. NumPy arrays go in and come out without copying, calls are vectorised over many states, and flowsheets are Python objects. The CLI also writes JSON and CSV that pandas reads directly. |

## Layout

```
v.mod               module manifest (Go's go.mod)
main.v              demo CLI, on database components
numeric/
  cubic.v           real roots of a cubic
  roots.v           bracketed root finding: bisection, Brent's method
  lu.v              dense LU with partial pivoting
  cholesky.v        Cholesky factorisation, the positive-definiteness test
  vec.v             small vector helpers
  integrate.v       adaptive Gauss–Kronrod quadrature
  newton.v          Newton for nonlinear systems: line search or dogleg
  minimize.v        Newton minimisation with a Levenberg–Marquardt shift
  rachford_rice.v   two-phase Rachford–Rice, negative flash, extreme K
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
  quality/          the sweep behind docs/data-quality.md
  chemsep/          the ChemSep XML importer: strict, with every repair reported
thermo/
  model.v           the Model contract F(T, V, n); volume solver; P, Z, ln φ
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
sweep/              the flash robustness sweep: seeded database mixtures,
                    every answer checked, every failure classified
docs/               roadmap, milestone plans, design notes, data format,
                    data-quality report, spike write-ups, V notes
tools/              check.vsh (the CI checks), gen_m*_fixtures.py (references),
                    gen_m3_dwsim.py (DWSIM 9's, by hand on Windows),
                    import_chemsep.vsh, data_quality.vsh and flash_sweep.vsh
                    (generated data and reports)
testdata/           reference fixtures from thermo, CoolProp, chemicals and
                    mpmath, with provenance
data/vendor/        the ChemSep 8.50 database, verbatim, with its licence
data/components/    the same, converted: one TOML file per compound
                    (generated by tools/import_chemsep.vsh)
spike/python_ffi/   calling reflux from Python and NumPy (M1 step 9)
spike/data_format/  the TOML component-file trial (M2 step 1)
spike/dwsim/        driving DWSIM 9 from Python, headless (M3 step 0)
.v-version          the V commit CI builds
```

Units: **SI everywhere inside** (K, Pa, m³, mol, J). Conversion happens only
at the edges, i.e. input files and printed output.

## Status

**M0 (foundations), M1 (thermo core) and M2 (component data) are done.**
The exit reviews are at the ends of [docs/m1-plan.md](docs/m1-plan.md) and
[docs/m2-plan.md](docs/m2-plan.md). What exists now:

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

Next is M3: the flash.

## Roadmap

The full roadmap is in [docs/roadmap.md](docs/roadmap.md). It covers
milestones with exit criteria, headless testing, C numeric backends and
how they get replaced, and CAPE-OPEN / FMI interop. The critical path:

**M0** foundations ✓ → **M1** thermo core ✓ → **M2** component data
(ChemSep) ✓ → **M3** flash → **M4** first flowsheet, compared with DWSIM 9. As
soon as M3 lands, **M10a** makes thermo and flash callable from Python with
NumPy arrays.

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
