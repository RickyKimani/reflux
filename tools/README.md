# tools

Scripts for working on reflux. Nothing here is part of the simulator.

| Script | What it does |
|---|---|
| `check.vsh` | Runs the same checks as CI: format, vet, tests, build, then `reflux demo`, `reflux solve` on the reference flowsheet with both reports, and a solve that must fail. |
| `import_chemsep.vsh` | Regenerates `../data/components/` from the vendored ChemSep XML. Run it after changing the importer or the vendored files; a test fails until you do. |
| `data_quality.vsh` | Regenerates `../docs/data-quality.md`: the sweep over every compound's data, PR against ChemSep's vapour pressures, and the outliers against CoolProp. Run it after changing the data or anything the sweep evaluates; a test fails until you do. |
| `flash_sweep.vsh` | Runs M3's robustness sweep, 10,000 seeded PT flashes of database mixtures, and writes `../docs/flash-robustness.md`. Build it with `-prod -cc gcc`: then it takes about three minutes, the build included, and about twenty under tcc. `sweep/sweep_test.v` runs the first 1,000 points in CI. |
| `benchmarks.vsh` | Writes `../docs/benchmarks.md` from `../testdata/m3_bench.json` and `../testdata/m4_bench.json`: `v run tools/benchmarks.vsh`. With `--time`, built with `-prod -cc gcc`, it also times reflux now, the flashes and the reference flowsheet, and prints that beside the fixtures' timings, to spot a regression, without touching the report. |
| `bench_case.v` | Times reflux on one benchmark case and prints JSON, for `gen_m3_bench.py`, which builds it with gcc `-prod`; with `--flowsheet`, a flowsheet file's solves, for `gen_m4_bench.py`; with `--capi`, V's own loop over M10a's points, for `gen_m10a_bench.py`. A `.v` program rather than a script, because a script runs whenever V compiles it. |
| `remake_bench.v` | Times what remaking each unit type costs (M8 step 3): `Loaded.with` on one of its numbers, against the unit's own calculation on its solved inlets. `v -prod -cc gcc -o remake_bench.exe tools/remake_bench.v`, then run it. A `.v` program, for the same reason as `bench_case.v`. |
| `study_bench.v` | Times the owner's 50-point study on 1, 2, 4 and 8 workers, in chunks of 8 and of 2, with what it allocates a point, beside a loop that allocates nothing on as many threads (M8 step 7). `v -prod -cc gcc -o study_bench.exe tools/study_bench.v`, then run it. |
| `gen_m8_sensitivities.py` | Writes `../testdata/m8_sensitivities.json` (M8 step 8): an isothermal first-order PFR's conversion and its derivatives by volume and by T, by mpmath at 50 digits, each closed form checked against mpmath's numerical derivative. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m8_sensitivities.py`. |
| `gen_m8_hock_schittkowski.py` | Writes `../testdata/m8_hock_schittkowski.json` (M8 step 9): Hock and Schittkowski's problems in SymPy, each optimum solved again by mpmath at 40 digits on the published active set and checked against the published optimum. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m8_hock_schittkowski.py`. |
| `gen_m8_series.py` | Writes `../testdata/m8_series.json` (M8 step 10): series reactions' optima in a PFR, inside the bounds, on a bound and on a constraint, with their multipliers, by mpmath at 50 digits. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m8_series.py`. |
| `gen_m8_fits.py` | Writes `../testdata/m8_fits.json` and the two CSV files in `../testdata/m8_fit/` (M8 steps 11 and 12): rate data and VLE data from known constants, seeded, SciPy's fits of them, and Student's t quantiles by mpmath. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m8_fits.py`. |
| `gen_m8_dwsim.py` | Writes `../testdata/m8_dwsim.json` (M8 step 13), by hand with DWSIM 9: its answers to the reference flowsheet's spec, the two-stage compression's least work and the series reactions' most, by its own forms driven unshown and by SciPy on its solves. Takes the parts `spec`, `compression` and `series`; a part not named keeps what the file has. About twenty minutes. |
| `gen_m8_bench.py` | Writes `../testdata/m8_bench.json` (M8 step 14), by hand with DWSIM 9 and V: the owner's study, an optimisation and a spec, by DWSIM's forms and by reflux (`bench_m8.v`), reflux before and after DWSIM's runs. About twenty minutes. |
| `bench_m8.v` | reflux's side of M8's benchmark: a study cold, warm or on workers, a solve, an optimisation or a spec, timed, as JSON. Built with gcc -prod by `gen_m8_bench.py`. |
| `gen_m6_activity.py` | Writes `../testdata/m6_activity.json` and `../testdata/m6_unifac_tables.json` (M6 steps 1 to 3): Wilson, NRTL, UNIQUAC and the three UNIFACs by thermo and, for the binary models, mpmath at 40 digits, on ChemSep's banks and tables as shipped; and thermo's copies of the published UNIFAC tables. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m6_activity.py`. Reproduces its output byte for byte. |
| `gen_m6_henry.py` | Writes `../testdata/m6_henry_iapws.json` (M6, after step 6): IAPWS's Henry's constants for gases in water by the `iapws` package, for the test of reflux's IAPWS form. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m6_henry.py`. Reproduces its output byte for byte. |
| `gen_m6_three_phase.py` | Writes `../testdata/m6_three_phase.json` (M6 step 7): thermo's `FlashVLN` on NRTL's and PR's three-phase and two-liquid states, polished on thermo's phases. Imports `gen_m6_gamma_phi.py`. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m6_three_phase.py`. Reproduces its output byte for byte. |
| `gen_m6_gamma_phi.py` | Writes `../testdata/m6_gamma_phi.json` (M6 steps 5 and 6): thermo's γ–φ liquid on ChemSep's coefficients, on three bases; Henry's law in water and in a mixed solvent; PT, bubble and dew, PH and PS flashes and azeotropes, each solved on thermo's own liquid where its `FlashVL` can't be trusted. Imports `gen_m6_activity.py`. `PYTHONIOENCODING=utf-8 uv run --python <3.13> tools/gen_m6_gamma_phi.py`. Reproduces its output byte for byte. |
| `capi_reference.v` | The V API's answers for a reference set of states, every C API function's, written to JSON with each float as the hex of its bits (M10a's D12). `python.vsh test` builds it with the library's compiler and flags, runs it, and the Python tests match the library against it bit for bit. Never committed: bits differ between compilers and maths libraries. |
| `gen_m10a_bench.py` | Writes `../testdata/m10a_bench.json`: properties and a PT flash at 1, 100 and 10,000 points, timed in V's own loop (through `bench_case.v`) and from Python through the C API, alternately in eight rounds. With V and gcc on the PATH: `uv run --python <3.11 or later> tools/gen_m10a_bench.py`, then `v run tools/benchmarks.vsh`. A few minutes. Its timings can't be reproduced byte for byte. |
| `gen_m4_bench.py` | Writes `../testdata/m4_bench.json`: the reference flowsheet on DWSIM's constants, solved by DWSIM 9 at its default recycle tolerances and at them scaled by 1e-8, and by reflux (through `bench_case.v`) by Wegstein and by direct substitution, all within a minute. By hand on Windows with DWSIM 9 and V: `uv run --python 3.13 tools/gen_m4_bench.py`, then `v run tools/benchmarks.vsh`. Its timings can't be reproduced byte for byte. |
| `gen_m1_fixtures.py` | Writes the M1 reference fixtures in `../testdata` from Caleb Bell's `thermo` and CoolProp. |
| `gen_m2_fixtures.py` | Writes `../testdata/m2_tdep.json`: every correlation form on real ChemSep records, evaluated with mpmath and checked against `chemicals`. Also writes `../testdata/m2_coolprop.json`, ChemSep against CoolProp for the 85 fluids both have, and `../testdata/m2_pr_psat.json`, PR's vapour pressure from `thermo` against ChemSep's for 20 compounds. |
| `gen_m3_fixtures.py` | Writes the M3 fixtures: Rachford–Rice, phase identification, stability and the flash references (`m3_thermo.json`), from `thermo`, `chemicals` and mpmath. Name some of them (`thermo`, `stability`, `pip`, `rachford_rice`) to write only those. |
| `gen_m3_bench.py` | Writes `../testdata/m3_bench.json`: reflux's (through `bench_case.v`), thermo's and DWSIM 9's timings and answers on the benchmark flashes, with DWSIM's constants for all three. The three time each case one after the other, within seconds: on a laptop the same flash's median moves by up to 3× from one session to another. Windows only, with DWSIM 9 installed, by hand: `uv run --python 3.13 tools/gen_m3_bench.py`, then `v run tools/benchmarks.vsh`. A few minutes. Its timings differ from run to run, so unlike the fixture generators it can't reproduce its output byte for byte; its answers it can. |
| `gen_m3_dwsim.py` | Writes `../testdata/m3_dwsim.json` from DWSIM 9 through pythonnet. Windows only, with DWSIM 9 installed; run it by hand: `uv run --python 3.13 tools/gen_m3_dwsim.py`. |
| `gen_m4_dwsim.py` | Writes `../testdata/m4_dwsim.json`, DWSIM 9's unit operations and reference flowsheet, with DWSIM's constants as component files in `../testdata/m4_dwsim/components/` and the reference flowsheet on them, `../testdata/flowsheets/reference_dwsim.toml`. By hand on Windows, with DWSIM 9 installed and V on the PATH, since it solves that file with reflux for the open loop: `uv run --python 3.13 tools/gen_m4_dwsim.py`. Name `units` or `reference` to run only that part. About five minutes. |
| `golden.vsh` | Writes the golden flowsheets' reports to `../testdata/flowsheets/golden/`: `v run tools/golden.vsh`. Run it only when a change is meant to move the answers, or for a new version, which each report carries. Say why in the commit, and in the CHANGELOG's next entry. |
| `python.vsh` | Builds the C API (`../capi/`) into the Python package, `../python/reflux/_lib/reflux.dll` (or `.so`), with gcc `-prod`, exporting exactly the functions `../capi/reflux.h` declares: `v run tools/python.vsh build`, or `build -debug` for a quicker build. `test` builds it and the V reference (`capi_reference.v`), both with `-ffp-contract=off` so that they agree to the bit where gcc's target has FMA, then runs the package's tests with uv and pinned pytest, NumPy and matplotlib (for the exit script, `../python/examples/phase_diagrams.py`); `--python <path>` picks the interpreter. gcc only: tcc can't build a usable shared library. |
| `release.vsh` | The release gate: `v run tools/release.vsh`. It refuses, with every reason, unless this commit can be tagged as `../v.mod`'s version: on `main` with a clean tree, the CHANGELOG's newest entry that version and dated, the version untagged and after the latest release, and every golden report changed since then named in the entry. Then it runs `check.vsh` and prints the tag commands; it tags nothing. The logic is the `release` module's, and tested. `../CHANGELOG.md` has the whole procedure. |

```bash
v run tools/check.vsh
```

```bash
v run tools/check.vsh -cc gcc
```

The first form uses tcc, V's bundled compiler. It's fast, so use it while
working. The second also makes an optimised `-prod` build; gcc compiles much
more slowly than tcc. The script runs up to four test programs at once
(`VJOBS` overrides that), and builds with `-nocache` under gcc, where V
0.5.2's module cache makes each build three times slower.

## Fixture generators

These use Python, because that's where the reference libraries are. They
run by hand, never in tests: the V tests only read the JSON they write.

```bash
uv run tools/gen_m1_fixtures.py
```

```bash
uv run tools/gen_m2_fixtures.py
```

The script declares its dependencies inline, with pinned versions (PEP 723),
so `uv run` builds the environment itself and re-running gives the same
numbers. Regenerate when a fixture needs new points or a new library
version. The provenance block at the top of each file says which versions
produced it.
