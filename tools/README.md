# tools

Scripts for working on reflux. Nothing here is part of the simulator.

| Script | What it does |
|---|---|
| `check.vsh` | Runs the same checks as CI: format, vet, tests, build, demo. |
| `import_chemsep.vsh` | Regenerates `../data/components/` from the vendored ChemSep XML. Run it after changing the importer or the vendored files; a test fails until you do. |
| `data_quality.vsh` | Regenerates `../docs/data-quality.md`: the sweep over every compound's data, PR against ChemSep's vapour pressures, and the outliers against CoolProp. Run it after changing the data or anything the sweep evaluates; a test fails until you do. |
| `flash_sweep.vsh` | Runs M3's robustness sweep, 10,000 seeded PT flashes of database mixtures, and writes `../docs/flash-robustness.md`. Build it with `-prod -cc gcc`: then it takes about three minutes, the build included, and about twenty under tcc. `sweep/sweep_test.v` runs the first 1,000 points in CI. |
| `benchmarks.vsh` | Writes `../docs/benchmarks.md` from `../testdata/m3_bench.json`: `v run tools/benchmarks.vsh`. With `--time`, built with `-prod -cc gcc`, it also times reflux now and prints that beside the fixture's timings, to spot a regression, without touching the report. |
| `bench_case.v` | Times reflux on one benchmark case and prints JSON, for `gen_m3_bench.py`, which builds it with gcc `-prod`. A `.v` program rather than a script, because a script runs whenever V compiles it. |
| `gen_m1_fixtures.py` | Writes the M1 reference fixtures in `../testdata` from Caleb Bell's `thermo` and CoolProp. |
| `gen_m2_fixtures.py` | Writes `../testdata/m2_tdep.json`: every correlation form on real ChemSep records, evaluated with mpmath and checked against `chemicals`. Also writes `../testdata/m2_coolprop.json`, ChemSep against CoolProp for the 85 fluids both have, and `../testdata/m2_pr_psat.json`, PR's vapour pressure from `thermo` against ChemSep's for 20 compounds. |
| `gen_m3_fixtures.py` | Writes the M3 fixtures: Rachford–Rice, phase identification, stability and the flash references (`m3_thermo.json`), from `thermo`, `chemicals` and mpmath. Name some of them (`thermo`, `stability`, `pip`, `rachford_rice`) to write only those. |
| `gen_m3_bench.py` | Writes `../testdata/m3_bench.json`: reflux's (through `bench_case.v`), thermo's and DWSIM 9's timings and answers on the benchmark flashes, with DWSIM's constants for all three. The three time each case one after the other, within seconds: on a laptop the same flash's median moves by up to 3× from one session to another. Windows only, with DWSIM 9 installed, by hand: `uv run --python 3.13 tools/gen_m3_bench.py`, then `v run tools/benchmarks.vsh`. A few minutes. Its timings differ from run to run, so unlike the fixture generators it can't reproduce its output byte for byte; its answers it can. |
| `gen_m3_dwsim.py` | Writes `../testdata/m3_dwsim.json` from DWSIM 9 through pythonnet. Windows only, with DWSIM 9 installed; run it by hand: `uv run --python 3.13 tools/gen_m3_dwsim.py`. |

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
