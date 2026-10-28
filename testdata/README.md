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
| `m3_pip.json` | Phase identification from Caleb Bell's `thermo` 0.6.1: 1,264 PR and SRK states (1,327 roots), pure components and mixtures across the saturation curve and the critical point, plus liquids squeezed to within 0.02% of the covolume, each root with its volume and analytic PIP. Also pure CO₂ through the Soave zero, T₀ = Tc(1 + 1/m)², with the properties M1 checks | `thermo/phase_test.v`: volumes to 1e-9, PIP to 1e-7, every label the same; CO₂ through T₀ to 1e-9 |
| `m3_stability.json` | 513 PR (and some SRK) states for stability: binaries on grids of T, P and z, 62 near the methane/n-butane critical locus (33 of them 0.1 bar apart), a ternary and a 10-component gas condensate. For each binary, a brute-force scan of the tangent-plane distance over the whole composition range; for every state, Caleb Bell's `thermo` 0.6.1 FlashVL (agreeing with the scans everywhere), with each split's phases as a witness of instability | `flash/stability_test.v`: every verdict, the scans' global minima to 1e-11, every instability proved and every incipient phase stationary, the trivial solution caught, split phases stable, and a recorded iteration budget |
| `m3_rachford_rice.json` | 87 seeded Rachford–Rice cases: 2 to 20 components, K from 1e-30 to 1e25, trace phases down to 1e-15, K within 1e-8 of 1, negative flashes and absent components. The reference is the exact root for each case's floats (bisection at 100 digits with mpmath 1.3.0), with what `chemicals` 1.5.2's own solvers return | `numeric/rachford_rice_test.v`: compositions to 1e-15, β and L to 4 ulps of what the K-values allow, and no larger an error than any of `chemicals`' float solvers |
| `m3_thermo.json` | The flash references, from Caleb Bell's `thermo` 0.6.1 with a shared quadratic ideal-gas Cp: 83 PT states (binaries, a ternary, the gas condensate, near-critical points, pure n-butane and CO₂, and narrow-boiling n-butane/isobutane), each with its split, fugacity residual and molar V, H, S and U, and where `thermo`'s own PH, PS and TV flashes land from them; 36 pure-component saturation states (`FlashPureVLS`); 42 bubble and dew points, polished by Newton on `thermo`'s own ln φ and each checked to be a phase boundary | `flash/thermo_test.v`: splits to 1e-8, V, H, S and U to 1e-9, saturation and boundary fugacities to 1e-10, Psat to 1e-9, and reflux's PT flash sees each boundary. `flash/vf_test.v`: the bubble and dew points and saturation states to 1e-11, and every PT split's vapour-fraction round trip to 1e-8. `flash/ph_test.v`: every PT state's PH and PS round trip, and every saturation state's, to 1e-8. `flash/tv_test.v`: the same for TV and UV |
| `m3_bench.json` | `tools/gen_m3_bench.py`, by hand on Windows with DWSIM 9 | M3's benchmarks (D10, criterion 8): 12 cases (PT with two phases and one, PH by state, the bubble point at P, for 2, 5 and 10 components), 40 states each, with reflux's, thermo 0.6.1's and DWSIM 9's median and minimum times, timed case by case, and their answers at the first state. Peng–Robinson with DWSIM's constants and kᵢⱼ for all three programs. `benchmarks/benchmarks_test.v` checks that reflux's answers match; `tools/benchmarks.vsh` writes `docs/benchmarks.md` from it. The timings vary between runs, so this is the one fixture that can't be reproduced byte for byte. |
| `m3_dwsim.json` | 25 cases from DWSIM 9.0.5's Peng–Robinson package (PT, PH and PS by state, bubble and dew points) with DWSIM's constants, Cp and vapour-pressure coefficients, kᵢⱼ and flash settings. Each has DWSIM's answer and its own fugacity residual, the exact PR with the same constants (`thermo` 0.6.1, mpmath 1.3.0 for the Cp integrals), and the PR with DWSIM's rounded constants. Run by hand on Windows | `flash/dwsim_test.v`: reflux matches the exact PR to 1e-8 (split) and 1e-9 (H, S), is within 1e-4 of DWSIM and no further than the exact answer is, and DWSIM's differences are explained case by case. reflux's own bubble and dew points, and its PH and PS flashes by state, land on the exact PR's |

The M1 files are written by
[`../tools/gen_m1_fixtures.py`](../tools/gen_m1_fixtures.py), the M2 files
by [`../tools/gen_m2_fixtures.py`](../tools/gen_m2_fixtures.py), and the M3
files by [`../tools/gen_m3_fixtures.py`](../tools/gen_m3_fixtures.py),
except `m3_dwsim.json` and `m3_bench.json`, which
[`../tools/gen_m3_dwsim.py`](../tools/gen_m3_dwsim.py) and
[`../tools/gen_m3_bench.py`](../tools/gen_m3_bench.py) write.
