# M2 plan — component data

Status: **in progress.** Steps 0 and 1 are done: the official files are in
`data/vendor/chemsep/` (see their [NOTICE](../data/vendor/chemsep/NOTICE.md)),
and the format is TOML ([data-format.md](data-format.md)). The roadmap entry
is
[M2 in roadmap.md](roadmap.md#m2--component-data); this file is the detailed
plan.

## Goal

After M2, reflux can build any of the 642 compounds in ChemSep 8.50 by
name, CAS number or InChIKey. Each compound comes with:

- its constants, in SI per-mol units, each traceable to its source;
- the temperature-dependent correlations the flash and energy balances
  need (vapour pressure, ideal-gas Cp, liquid density, heat of
  vaporisation, liquid Cp), with the integrals H and S need;
- PR and SRK kᵢⱼ for the pairs ChemSep has.

One call should go from names to a ready `thermo.Fluid`:

```v
db := components.load_default()!
fl := db.fluid(.pr, ['methane', 'n-butane'])!   // Cubic + ideal gas + kij
```

M3's flash then starts from real data instead of hand-typed constants.

## What the research found

All of this comes from looking at the files themselves (October 2026).
The plan was first written from the copies in DWSIM's repository (v8.3,
2022). It was then checked against the official ChemSep LITE 8.50 install
(step 0), and the official files are what reflux uses. Where the two
differ, the official values are given here.

| Question | Answer |
|---|---|
| What is the data? | ChemSep pure-component databases `chemsep1.xml` (v8.49, **433 compounds**) and `chemsep2.xml` (v8.50, **209 compounds**). They don't overlap: 642 compounds in all. |
| Licence | Stated in each file's root element: *"ChemSep v8.49 pure component data 1 - Copyright (c) Harry Kooijman and Ross Taylor (2024)"*, with a link to the **Artistic License 2.0**, and the licence text ships in the same folder. |
| Defects | The XML is **not well-formed**: 43 SMILES attribute values end in a stray `<`. Strict XML parsers reject the files as shipped. |
| New in the official files | `InChIKey` and `PubChemID` per compound. Real revisions to Tc, Pc or ω for a few dozen compounds since v8.3; most other differences are reformatting. |
| What Artistic 2.0 allows | Verbatim redistribution with notices (§2); a converted copy is a "Modified Version", distributable if the changes are documented and it stays under Artistic 2.0 (§4(c)(i)); bundling with MIT code is fine (§7). |
| Identifiers | Name, CAS number (present, well-formed and unique in every record), InChIKey, PubChemID, SMILES, a structural formula, a family code. |
| Units | Per **kmol**: J/kmol, m³/kmol, kg/kmol. K and Pa elsewhere. |
| What M3 needs, and coverage | Tc, Pc, Vc, Zc, Tb, M, ω, ΔH°f, ΔG°f, S°: 431/431. Vapour pressure, ideal-gas Cp and heat of vaporisation: 431/431. Liquid density 424, liquid Cp 422. |
| Equation forms | Ideal-gas Cp is ChemSep **eq. 16**, A + exp(B/T + C + DT + ET²), for 429 compounds. It has no closed-form integral. Vapour pressure is DIPPR 101 (421), eq. 10 Antoine (8) and **eq. 200 (2, form unknown)**. Liquid density is DIPPR 105, heat of vaporisation DIPPR 106. Transport properties use 101, 102, 16 and others. |
| Extras | Mathias–Copeman c₁–c₃ in the XML (the EOS isn't stated), plus a separate `mc-PR.lib` explicitly fitted for Peng–Robinson; UNIQUAC r and q; UNIFAC, modified-UNIFAC and ASOG group assignments; Lennard-Jones parameters; a name-synonyms file. |
| kᵢⱼ | `pr.ipd` and `srk.ipd`: 292 pairs each, keyed by CAS pair, "DECHEMA Peng-Robinson [SRK] Parameters", "Copyright 2009 Harry Kooijman and Ross Taylor". **No licence named.** Also NRTL, UNIQUAC and Wilson files, for M6. |
| An oracle for tests | Caleb Bell's `chemicals` 1.5.2 implements DIPPR 100–107, 114–116 and 127, including integrals. It reproduces n-butane's Psat(350 K) = 9.467 bar and ρL(300 K) = 9.836 kmol/m³ from the ChemSep coefficients. |

The ChemSep-specific equations, from ChemSep's PCDmanager documentation:

| Eq. | Form | Eq. | Form |
|---|---|---|---|
| 1–5 | A + BT + CT² + DT³ + ET⁴, truncated to 1–5 terms | 10 | exp(A − B/(T + C)) (Antoine; K, Pa) |
| 11–15 | exp of the 1–5 polynomial | 16 | A + exp(B/T + C + DT + ET²) |
| 17 | A + exp(B + CT + DT² + ET³) | 45 | the integral of the eq. 5 polynomial |

## Exit criteria

All checked by `v -W test .`, so CI enforces them.

1. **Import is complete.** All 644 records import as 642 compounds,
   despite the 43 malformed SMILES values, which are repaired and listed.
   The import report lists zero unparsed records. The two compounds that
   `chemsep2.xml` lists twice are resolved by the documented rule (D8).
   Unsupported equation forms (such as eq. 200) are kept verbatim and
   reported, never dropped.
2. **The data files are reproducible.** Re-running the importer on the
   vendored XML regenerates `data/components/` byte for byte, and the
   vendored files still match the SHA-256 sums in their NOTICE. CI checks
   both, so neither the converted data nor its source can drift.
3. **Units are right.** Every scalar is SI per mol. Tests check a sample
   against hand conversions, for example n-butane's M of 58.123 kg/kmol
   becoming 0.058123 kg/mol, and ΔH°f in J/kmol becoming J/mol.
4. **Correlations are exact.** Every equation form used by the five
   M3-critical properties is implemented with its value, dY/dT,
   ∫Y dT and ∫Y/T dT:
   - DIPPR forms match `chemicals.dippr` to 1e-12 relative;
   - ChemSep-only forms match an independent evaluation to 1e-12;
   - numerical integrals (eq. 16) match high-resolution quadrature to 1e-10.
5. **The data is sane, and outliers are visible.** A sweep over all
   compounds checks that each correlation is finite, and positive where it
   must be, across its [Tmin, Tmax], and that Psat(Tb) is within 2% of
   1 atm. Failures go into a committed data-quality report rather than
   failing the build, but the report can't silently grow.
6. **Cross-checks are recorded.** Deviations are recorded at generation and
   can't grow, as in M1:
   - for fluids that also have a reference equation in CoolProp (about a
     dozen: light hydrocarbons, CO₂, N₂, water and others), ChemSep's Tc,
     Pc, ω and M, and its eq.-16 ideal-gas Cp, against CoolProp;
   - for 20 representative compounds, PR's Psat (from ChemSep's constants)
     against ChemSep's own vapour-pressure correlation, across Tr 0.5–0.95.
7. **Lookup works.** By CAS, by name (case-insensitive) or by alias. An
   unknown name is an error that lists the closest matches.
8. **kᵢⱼ works.** PR and SRK kᵢⱼ load from the ChemSep files. The matrix for
   any component list is symmetric, and pairs with no data get 0 and are
   listed in a report.
9. **End to end.** `db.fluid(.pr, […])` builds a working `thermo.Fluid`. An
   integration test reproduces M1's n-butane Psat and enthalpy changes using
   database constants.
10. CI is green, and the docs are updated.

## Decisions

| # | Decision | Choice | Why |
|---|---|---|---|
| D1 | Licence and layout | Verbatim ChemSep XML in `data/vendor/chemsep/`, with the Artistic 2.0 text and a NOTICE (source, version, SHA-256). The converted files in `data/components/` are a Modified Version, also under Artistic 2.0 (§4(c)(i)), with a CHANGES file. Code stays MIT. | It's exactly what the licence allows, and it keeps code and data licences visibly separate. |
| D2 | kᵢⱼ files | **Shipped, as DWSIM does (your call).** They go in `data/vendor/chemsep/` too. The NOTICE says plainly that the files name no licence and are assumed to fall under the database's Artistic 2.0, with a note to confirm with the authors. | Practical; the assumption is written down where anyone redistributing will see it. |
| D3 | Where the XML comes from | **The official ChemSep LITE download (your call).** Done with ChemSep LITE 8.50; the files come from its `pcd` and `ipd` folders. | The Standard Version, straight from the copyright holders. |
| D4 | reflux's own format | TOML, one file per compound (`data/components/<CAS>.toml`): readable, diffable, and the same format a user writes for their own compounds. Every value carries a source tag. | A per-compound diff is reviewable; one large file isn't. Step 1 checks V's `toml` module can do it, with JSON as the fallback. |
| D5 | Units of correlations | Coefficients are stored **verbatim**, along with the correlation's native output unit, and the evaluator converts the *result* to SI per mol. Scalars are converted at import. | Forms like A + exp(…) can't be rescaled coefficient by coefficient. Verbatim coefficients also compare directly with ChemSep and `chemicals`. |
| D6 | Where it lives | A new module `components`, which imports `thermo` and builds `thermo.Component`, `thermo.IdealGas` and `thermo.Fluid`. | `thermo` stays about models. Data handling is its own concern. |
| D7 | Outside [Tmin, Tmax] | Evaluate anyway (extrapolate) and say so: each evaluation returns the value plus an in-range flag. | A flash's trial temperatures often step outside a correlation's range briefly. Refusing would break the iteration, while silence would hide bad data. |
| D8 | Overlapping compounds | The first occurrence of a CAS number, in reading order (`chemsep1.xml`, then `chemsep2.xml`, top to bottom), wins. Each duplicate is listed in the import report with the elements that differ. | Simple and deterministic. In 8.50 the only duplicates are two compounds listed twice within `chemsep2.xml`. In both, the first copy is the more complete and more recently edited one (step 3). |
| D9 | Ideal-gas Cp from the database | A new `CpForm` variant, `thermo.Correlated`, wrapping a correlation, so eq. 16 drops into M1's `Fluid` unchanged. Integrals use adaptive Gauss–Kronrod quadrature (`numeric.integrate`) when there's no closed form. A failed integral is NaN rather than an error. | It reuses M1's whole total-property machinery, and `Fluid.props` keeps its signature. |
| D10 | Mathias–Copeman parameters | The XML's values are imported and stored but not used by default, since their EOS isn't stated. `mc-PR.lib`, fitted for PR, is the better source; vendoring it is a small addition to step 6. | Using parameters fitted for one EOS in another would be guessing. |
| D11 | Who writes what | As in M1, I write it all unless you claim a step. Good candidates are marked. | |

## Steps

```mermaid
flowchart LR
  S0[0 Install ChemSep LITE] --> S3
  S1[1 Data format spike] --> S4
  S2[2 Correlations] --> S3[3 XML importer]
  S3 --> S4[4 TOML writer and reader]
  S4 --> S5[5 Database API]
  S5 --> S6[6 kij banks]
  S2 --> S7[7 Reference fixtures]
  S5 --> S8[8 Sanity sweep]
  S5 --> S9[9 End to end]
  S6 --> S9
  S7 --> S10[10 Docs and exit review]
  S8 --> S10
  S9 --> S10
```

### 0. Get the official files — done
- ChemSep LITE **8.50** is installed at
  `%LOCALAPPDATA%\ChemSepL8v50`. `chemsep1.xml`, `chemsep2.xml`, `pr.ipd`,
  `srk.ipd` and ChemSep's own `Artistic_license_2_0.txt` were copied into
  `data/vendor/chemsep/`, and they're byte-identical to the installed files.
- [NOTICE.md](../data/vendor/chemsep/NOTICE.md) records the versions,
  SHA-256 sums, the licence statements quoted from the files, the kᵢⱼ
  licence assumption (D2), the 43 malformed SMILES values, and how the
  files differ from DWSIM's.
- `.gitattributes` marks `data/vendor/` as `-text`, so Git doesn't rewrite
  the CRLF line endings and the checksums stay true.

### 1. Data format spike — done: TOML
- n-butane was written by hand in the schema
  ([spike/data_format/106-97-8.toml](../spike/data_format/106-97-8.toml)) and
  read with V's `toml` module
  ([spike/data_format/spike.v](../spike/data_format/spike.v)). Everything the
  schema needs works: nested tables, arrays, inline-table overrides with
  their own source, quoted numeric group keys, comments, file order kept,
  and missing keys and malformed files reported.
- On all 53,977 distinct numbers in the ChemSep XML, the TOML parser and
  `strconv` agree exactly. V's `f64.str()` round-trips and always prints a
  TOML float, so the writer can be deterministic.
- Parsing takes about 1 ms per file under tcc. So the importer also writes an
  `index.toml` for name lookup, and components are parsed only when asked
  for.
- Two V gotchas found:
  - `toml.Any.f64()` on a string silently returns 0.0, so the reader checks
    types;
  - string interpolation of an inline `.map(…)` result prints the array's
    internals, a V bug (assign to a variable first).
- The schema is written down in [data-format.md](data-format.md).

### 2. Correlations — done
- A new module, `tdep`. It has `Correlation` (form, A–E, Tmin, Tmax, and Tc
  for the forms that need it) with `eval`, `deriv`, `integral` and
  `integral_over_t`. It's pure mathematics in the correlation's own units:
  the native unit and the conversion to SI belong to the `components`
  layer (D5).
- The forms: ChemSep 1–5, 10–17 and 45; DIPPR 100–102, 104–107 and 114–116.
  - Derivatives are closed-form for all of them.
  - ∫Y dT is closed-form for 1–5, 45, 100, 104, 107, 114 and 116, and
    ∫Y/T dT for 1–5, 45, 100, 104 and 107.
  - The rest use a new `numeric.integrate` (below).
- `new()` refuses forms it can't evaluate, so the data layer keeps them
  verbatim and reports them. ChemSep 8.50 uses five such forms:
  - 120 and 121, relative permittivity (324 and 5 records). These matter
    only for electrolytes.
  - 200, vapour pressure (5 records). Step 5 falls back to the Antoine
    record.
  - 118, liquid viscosity (1 record).
  - 124, liquid Cp (1 record).
- Findings that shaped the code:
  - 101 eq. 106 records (heats of vaporization, surface tensions, water's
    liquid density) end a little past their compound's own Tc, by up to
    48 K. n-butane's ends at 425.16 K against Tc = 425.12 K. Form 106 is for
    properties that vanish at Tc, so it returns 0 at and above Tc instead
    of NaN. 114 and 116 aren't defined above Tc, and `new()` refuses such
    ranges; none occur in the data.
  - Some solid-phase records start at 0 K, so Tmin = 0 is allowed.
- `numeric.integrate` is globally adaptive Gauss–Kronrod 7–15 (QUADPACK's
  QAG): it always splits the panel with the largest error, uses QUADPACK's
  error estimate and roundoff floor, and has a 2000-panel budget. My first
  version halved locally with a halved tolerance per side. When asked for
  1e-14, below what roundoff allows, it ran for two minutes; QAG settles in
  one panel. A divergent integral or noise is an error, not a hang.
- Fixtures: [tools/gen_m2_fixtures.py](../tools/gen_m2_fixtures.py) writes
  [testdata/m2_tdep.json](../testdata/m2_tdep.json).
  - It has 68 cases covering all 24 forms: up to three real ChemSep records
    per form, n-butane's main records, and synthetic sets for the forms
    ChemSep never uses. Two synthetic sets cover hard cases:
    - an Aly–Lee set run from 2 K to 10⁵ K, where sinh overflows a float;
    - an eq.-104 set where every term counts (see the mutation check).
  - Values, dY/dT and both integrals are mpmath at 40 digits, from formulas
    written separately from the V code.
  - For the DIPPR forms the script checks itself against `chemicals` 1.5.2,
    and stops rather than write a bad fixture: values to 1e-12, derivatives
    to 1e-9, closed-form integrals to 1e-10.
- Tests: values match mpmath to 1e-13, and derivatives and integrals to
  1e-11, which meets criterion 4. The DIPPR values also match `chemicals`
  to 1e-12. Unit tests cover:
  - the refusals;
  - the zero-fill of missing coefficients, `unused_nonzero` and the range ends;
  - closed forms against quadrature;
  - the behaviour at and above Tc;
  - Aly–Lee at 1 mK and 10⁹ K, with C = 0 and with negative C and E;
  - the hyperbolic helpers' series.
- Mutation check: of 26 planted mistakes, 25 were caught first time. The
  one that got through was a wrong D/T⁷ coefficient in the eq.-104
  integral. In ChemSep's second virial coefficients the D and E terms are
  about 1e-14 of the total over each record's range, so no test could see
  them. The synthetic eq.-104 set fixes that, and it and four more
  mutations of the D and E terms are now caught: 30 of 30.

### 3. ChemSep XML importer — done
- New module `components`, for now just the record model
  ([components/record.v](../components/record.v)). Its struct field names
  are the component-file keys. Values are set by key through V's
  compile-time reflection, so an unknown key is an error, and the step 4
  writer will take its key order from the same structs.
- New module `components.chemsep`, the importer.
  - [fields.v](../components/chemsep/fields.v) is the one table that maps
    each of the 68 ChemSep elements to a key, a kind, the unit it must carry,
    and whether it's per kmol.
  - [chemsep.v](../components/chemsep/chemsep.v) parses with V's
    `encoding.xml`.
  - The provenance is `ChemSep v8.49 (data/vendor/chemsep/chemsep1.xml)`
    and the like, with the version taken from the file itself.
- Strict, so a new ChemSep release can't be misread silently. These are all
  errors:
  - an unknown element or attribute;
  - a unit other than the expected one;
  - a malformed, NaN or infinite number;
  - a gap in a coefficient set;
  - an element twice in a compound;
  - text where elements belong;
  - a CAS number that fails its check digit. It becomes a file name, so
    this also keeps anything path-like out.
- Per-kmol values are converted by moving the decimal point in the text
  before parsing. Dividing by 1000 would be one unit in the last place off
  for 220 of the 2,637 distinct per-kmol values.
- The whole vendored database imports in about 4 s under tcc. What it found:
  - **644 records, 642 compounds.** Trifluoromethane and 3-chloropropene
    are each in `chemsep2.xml` twice. The first copy is the more complete
    and more recent one both times, so D8 now reads "first occurrence wins".
  - **48 repairs** (47 in step 3; step 4's round trip found the 48th).
    - 41 stray SMILES `<`s; the other 2 of the 43 are in the dropped copies.
    - Six ASOG group lists that name a group twice. Summing the two is right
      every time: diglyme's 2 + 1 ether groups for its three ether oxygens.
    - One zero group count, diphenyl disulfide's UNIFAC-LLE list. It's
      dropped: a zero count contributes nothing.
  - **Findings, kept verbatim and reported:**
    - 42 impossible correlation ranges, such as chloroform's solid Cp over
      [283, 0] K and several liquid viscosities ending at 1 K;
    - 334 correlations in forms reflux can't evaluate, 327 of them
      permittivities;
    - 8 records with suspect unused coefficients;
    - one correlation with no coefficients, para-hydrogen's RPP Cp;
    - 5 partial Mathias–Copeman sets;
    - 1,2-dichloroethylene carbonate, with no Tc, Pc or ω.
- Tests ([chemsep_test.v](../components/chemsep/chemsep_test.v)): the real
  import's counts, duplicates, repairs and findings. They were all
  cross-checked with an independent Python/ElementTree survey of the XML.
  Also:
  - n-butane against hand conversions, value by value;
  - the decimal shift;
  - CAS validation;
  - a synthetic file for each refusal, the repairs, and duplicates across
    files;
  - every table key against the record structs.
- Mutation check: 21 planted mistakes, all caught. They covered the decimal
  shift, each strictness check, the repairs, the duplicate rule, each kind
  of finding, the CAS check digit, a per-kmol flag in the table, and the
  reflection setters.

### 4. TOML writer and reader — done
- `components.write_toml` and `components.parse_toml`/`read_toml`
  ([components/file.v](../components/file.v)). Both walk the record structs
  by compile-time reflection, so key order and key names come from one
  place.
  - **The writer is deterministic:** struct order, empty tables left out,
    floats in V's shortest round-trip form, strings escaped. SMILES have
    `\` in them.
  - **It refuses** NaN, infinities, control characters, and a `sources`
    override that names no value.
  - **The reader is strict.** These are errors naming the file and key:
    - an unknown key;
    - a wrong type;
    - a missing required key;
    - a file not named after its CAS number;
    - a group number that isn't canonical, or a count that isn't positive.
  - V's toml parser reads `nan` as a string and `inf` as a huge integer.
    So a float field takes an integer only if it's exact (|i| ≤ 2⁵³), and
    both are refused.
- `index.toml` ([components/index.v](../components/index.v)) gives each
  compound's name, InChIKey and aliases, sorted by CAS number. Lookup by
  name reads only this file.
- `components.chemsep.generate` renders the whole of `data/components/`:
  - 642 compound files;
  - `index.toml`;
  - `README.md`: generated, don't edit, the licence;
  - `CHANGES.md`: how this Modified Version differs from ChemSep's Standard
    Version (D1). It's also criterion 1's import report: every repair, the
    dropped duplicates, and every finding.
- [tools/import_chemsep.vsh](../tools/import_chemsep.vsh) writes exactly
  that set: it rewrites what changed and removes what's no longer
  generated. 645 files, 2.9 MB.
- Tests:
  - **Criterion 2, in [generate_test.v](../components/chemsep/generate_test.v):**
    - the committed directory equals a fresh generation byte for byte,
      with no extra or missing files;
    - all five vendored files match the size and SHA-256 in their NOTICE.
    CI runs both as part of `v test .`.
  - Every committed file reads back as the record it was written from (V's
    `==` compares the structs deeply).
  - The index lists every file.
  - [file_test.v](../components/file_test.v) checks the writer's
    determinism and escaping, and gives each reader refusal its own case.
- The round trip found a data defect the importer had passed: a zero group
  count. It's now a reported repair (step 3).
- Mutation check: 25 planted mistakes, 23 caught on the first run.
  - **Missed, and fixed:** the index reader's unknown-key check. The test
    appended its stray key at the end of the file, where TOML puts it
    inside the last compound table, so a different check refused it. The
    test now puts the key at the top level, and each index refusal is
    matched to its message. That mutation and two new ones are caught.
  - **Missed, and left:** the reader's NaN/infinity check on floats. V's
    toml never produces a non-finite float (see above), so no file can
    reach it. It stays, documented, in case the parser changes.

### 5. Database API — done
- `components.load(dirs)` and `load_default()`
  ([components/database.v](../components/database.v)).
  - The built-in directory is `data/components/`, found through `@VMODROOT`
    at compile time, or `$REFLUX_DATA`.
  - A directory with an `index.toml` is read through it, and a stale index
    (one listing a missing file) is an error. A directory without one, such
    as your own, has its `*.toml` files read.
  - A later directory's compound replaces an earlier one with the same CAS
    number, whole: that's how user components override built-in ones.
  - Files are parsed only when a compound is asked for.
- `db.find(query)` matches a CAS number, an InChIKey, or a name or alias
  (ignoring case and spaces). Criterion 7:
  - A miss is a typed `LookupError` with the five closest names: those
    containing the query first, shortest first, then by edit distance.
    "butane" suggests N-butane and Isobutane, "watr" Water.
  - A query matching several compounds is a `LookupError` listing them all.
    That happens in the shipped data: ChemSep gives carbon monoxide
    formaldehyde's InChIKey. A name used again in a user directory
    triggers it too.
- What a record builds ([components/build.v](../components/build.v)):
  - `rec.component()` builds a `thermo.Component`. It fails, naming what's
    missing, for the one compound without Tc, Pc and ω.
  - `rec.ideal_gas()` builds a `thermo.IdealGas`. Its Cp is a new
    `thermo.Correlated` form (D9) wrapping the `ideal_gas_cp` correlation,
    scaled from J/(kmol·K).
  - `rec.property(key)` gives any correlation converted to SI per mol. Its
    `at(T)` returns the value and whether T is in the fitted range (D7).
- The step 3 findings, decided:
  - a correlation with no coefficients or an impossible range is refused,
    and the message points to CHANGES.md;
  - a vapour pressure that can't be used falls back to the Antoine record,
    with a note saying why. That covers 26 compounds: 5 in eq. 200, the rest
    without a usable eq.-101 record. Only the stub compound has no vapour
    pressure at all.
  - the reduced-temperature forms take the compound's Tc unless the
    correlation has its own.
- In thermo:
  - **`Correlated`** integrates numerically where it must. A failed
    integral, which takes a correlation that overflows, returns NaN, so
    `Fluid.props` stays infallible and M1's API is unchanged.
  - **`AlyLee`** now evaluates through tdep's eq. 107, so the duplicated
    overflow-safe helpers are gone from thermo. Their naive-form test moved
    to tdep. M1's fixtures still pass unchanged.
- Tests ([database_test.v](../components/database_test.v)): every lookup
  path and refusal, overrides and clashes in temporary user directories,
  n-butane's conversions against the formulas written out, and a sweep
  building every compound. Its counts were checked with Python's tomllib:
  - an ideal gas for 640 compounds (none for the stub or for R-1234yf);
  - a component for 641;
  - a vapour pressure for 641.
- Mutation check: 22 planted mistakes, all caught. They covered each lookup
  path, the suggestion ranking, directory order, the stale-index check,
  the unit scales, the fallback, each refusal, Tc and ΔH°f wiring, the
  in-range flag, and thermo's new Cp form.

### 6. kᵢⱼ banks
- An IPD parser; a symmetric matrix for a component list; a report of
  missing pairs.
- Tests: a few pairs against the file by hand; symmetry; reporting of
  missing pairs.

### 7. Reference fixtures
- The DIPPR values from `chemicals` were done in step 2. What's left is
  CoolProp's constants and Cp⁰ for the overlapping fluids, added to
  `tools/gen_m2_fixtures.py`.
- The V tests read them (criterion 6).

### 8. Sanity sweep and data-quality report
- Over every compound: finiteness, sign, Psat(Tb) ≈ 1 atm, and the PR-Psat
  comparison for the 20 representative compounds.
- `docs/data-quality.md` lists what fails, with numbers, generated by a
  tool and checked by a test. It also lists what step 2 found:
  - the records in forms reflux doesn't evaluate;
  - the 101 eq.-106 records whose range runs past Tc;
  - the 8 records whose unused coefficients aren't zero (`unused_nonzero`):
    four eq.-2 solid densities with a C, and four eq.-4 amine liquid Cps
    with an E. Those suggest a wrong equation number.
- And what step 5 found:
  - an InChIKey shared by two compounds: carbon monoxide has formaldehyde's;
  - compounds with no usable ideal-gas Cp: the stub record, and R-1234yf
    (754-12-1);
  - the 26 compounds whose vapour pressure falls back to Antoine.
    ChemSep's eq. 200 looks like a Wagner equation, with Pc as its fifth
    coefficient and Tmax at Tc. Supporting it would turn 5 of those back
    into their own fits.

### 9. End to end
- `db.fluid(.pr, names)`; the integration test from criterion 9; the demo
  updated to use database components.

### 10. Docs and exit review
- README, roadmap, the data-format doc, the V notes (`encoding.xml`,
  `toml`, file walking), and a criterion-by-criterion review, as in M1.

## Risks

| Risk | Effect | Fallback |
|---|---|---|
| V's `toml` module can't handle the schema | step 1 changes format | JSON via `json2`, which M1 already uses |
| `encoding.xml` is slow, or rejects the malformed files | import fails | **Happened as expected, and handled.** V's parser rejects the stray `<`, so the importer escapes that one pattern first and lists each repair. Speed is fine: 4 s for both files under tcc, and import runs only when the vendored XML changes. |
| ChemSep LITE's files differ from DWSIM's copies | different numbers | **Happened**, and harmlessly: the official files win (D3), and the NOTICE records the differences. |
| The kᵢⱼ licence assumption (D2) turns out wrong | files must be removed | they're isolated in `data/vendor/chemsep/`, and the importer can read them from a local install instead |
| Bad entries in the data | wrong results downstream | the sanity sweep and data-quality report make them visible; no silent fixes |
| Repository size (about 5 MB vendored plus the converted files) | slower clones | acceptable at this scale; revisit with Git LFS if it grows |

## Not in M2

- Transport properties (viscosity, thermal conductivity, surface tension):
  their correlations are imported verbatim but not wired up until M5.
- UNIFAC, UNIQUAC and NRTL parameters: imported and stored raw, used in M6.
- Estimation methods (Joback and others): M6.
- Petroleum pseudo-components: after M7.
- The Soave α → 0 guard carried over from M1: still open, and not part of
  M2.
