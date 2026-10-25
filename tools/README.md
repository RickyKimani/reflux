# tools

Scripts for working on reflux. Nothing here is part of the simulator.

| Script | What it does |
|---|---|
| `check.vsh` | Runs the same checks as CI: format, vet, tests, build, demo. |

```bash
v run tools/check.vsh
```

```bash
v run tools/check.vsh -cc gcc
```

The first form uses tcc, V's bundled compiler. It's fast, so use it while
working. The second also makes an optimised `-prod` build; gcc compiles much
more slowly than tcc.

Fixture generators for [`../testdata`](../testdata) will also live here.
They may use Python (for `thermo` or CoolProp), because they run by hand,
not in tests.
