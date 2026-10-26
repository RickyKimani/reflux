# testdata

Fixtures that tests read, borrowed from Go's `testdata/` convention. V has
no convention of its own here.

Rules:

- Tests only **read** from here. Nothing in `v test` writes to it.
- Every file records where each value came from: the source (book and
  table, paper, database and version, or tool and version), and the date
  it was generated.
- Generated fixtures are produced by a script in [`../tools`](../tools) and
  committed. Tests never need Python, CoolProp or network access to run.
- JSON is the default format, read with V's `json2` module.

| File | Contents | Checked by |
|---|---|---|
| `m1_thermo.json` | PR and SRK, pure n-butane and a methane / n-butane / CO₂ ternary, computed by Caleb Bell's `thermo` (same equations and constants as reflux) | `thermo/fixtures_test.v`: agreement to 1e-9 |
| `m1_coolprop_butane.json` | n-butane from CoolProp's reference equation, with PR's deviations from it at generation time | `thermo/fixtures_test.v`: PR's deviations match and haven't grown |

Both are written by [`../tools/gen_m1_fixtures.py`](../tools/gen_m1_fixtures.py).
