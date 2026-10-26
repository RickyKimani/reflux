# tools

Scripts for working on reflux. Nothing here is part of the simulator.

| Script | What it does |
|---|---|
| `check.vsh` | Runs the same checks as CI: format, vet, tests, build, demo. |
| `import_chemsep.vsh` | Regenerates `../data/components/` from the vendored ChemSep XML. Run it after changing the importer or the vendored files; a test fails until you do. |
| `gen_m1_fixtures.py` | Writes the M1 reference fixtures in `../testdata` from Caleb Bell's `thermo` and CoolProp. |
| `gen_m2_fixtures.py` | Writes `../testdata/m2_tdep.json`: every correlation form on real ChemSep records, evaluated with mpmath and checked against `chemicals`. |

```bash
v run tools/check.vsh
```

```bash
v run tools/check.vsh -cc gcc
```

The first form uses tcc, V's bundled compiler. It's fast, so use it while
working. The second also makes an optimised `-prod` build; gcc compiles much
more slowly than tcc.

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
