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
| `m2_tdep.json` | Every correlation form `tdep` supports, on real ChemSep 8.50 records and a few synthetic sets: value, dY/dT, ∫Y dT and ∫Y/T dT at 40 digits (mpmath), and `chemicals` 1.5.2 values for the DIPPR forms | `tdep/fixtures_test.v`: values to 1e-13, derivatives and integrals to 1e-11, `chemicals` to 1e-12 |
| `m2_coolprop.json` | The 85 CoolProp 8.0.0 fluids that ChemSep also has: CoolProp's Tc, Pc, ω, M, ideal-gas Cp and vapour pressure, with ChemSep's deviation from each at generation. Also n-butane's Psat at four temperatures and its enthalpy at ten states. | `components/coolprop_test.v`: each deviation, recomputed through reflux, matches and hasn't grown. `components/fluid_test.v`: the end-to-end Psat and enthalpy changes. |
| `m2_pr_psat.json` | Peng–Robinson's vapour pressure (Caleb Bell's `thermo` 0.6.1, from ChemSep's Tc, Pc, ω) against ChemSep's own correlation, for 20 representative compounds at Tr 0.5–0.95 | `components/quality/quality_test.v`: reflux's PR matches `thermo` to 1e-9, and each deviation matches and hasn't grown |

| `m3_rachford_rice.json` | 87 seeded Rachford–Rice cases: 2 to 20 components, K from 1e-30 to 1e25, trace phases down to 1e-15, K within 1e-8 of 1, negative flashes and absent components. The reference is the exact root for each case's floats (bisection at 100 digits with mpmath 1.3.0), with what `chemicals` 1.5.2's own solvers return | `numeric/rachford_rice_test.v`: compositions to 1e-15, β and L to 4 ulps of what the K-values allow, and no larger an error than any of `chemicals`' float solvers |

The M1 files are written by
[`../tools/gen_m1_fixtures.py`](../tools/gen_m1_fixtures.py), the M2 files
by [`../tools/gen_m2_fixtures.py`](../tools/gen_m2_fixtures.py), and the M3
files by [`../tools/gen_m3_fixtures.py`](../tools/gen_m3_fixtures.py).
