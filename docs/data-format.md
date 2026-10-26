# Component data format

reflux keeps component data as plain-text TOML, **one file per compound**:
`data/components/<CAS number>.toml`. The importer writes them from the
vendored ChemSep XML (M2 step 4), and the same format is what you write by
hand for your own compounds. A per-compound file makes every data change a
small, reviewable diff.

This is format **1**. Every file says so in its first key, and the reader
refuses a format it doesn't know rather than guessing.

## An example, abbreviated

```toml
format = 1
name = "N-butane"
# Every value below comes from this source unless it says otherwise.
source = "ChemSep v8.49 (data/vendor/chemsep/chemsep1.xml)"

[ids]
cas = "106-97-8"
inchikey = "IJDNQMDRQITEOD-UHFFFAOYSA-N"
pubchem = 7843
smiles = "CCCC"
structure = "CH3(CH2)2CH3"
chemsep_index = 5
family = 1
aliases = ["butane"]

[constants]
tc = 425.12
pc = 3.796e+06
omega = 0.199
mw = 0.0581222
dhf = -125790.0
# A value from somewhere else carries its own source:
# tc = { value = 425.125, source = "CoolProp 8.0.0, n-Butane" }

[parameters]
rackett = 0.2728
mathias_copeman_c1 = 0.84209

[correlations.ideal_gas_cp]
eq = 16
coeffs = [44749.9, -338.141, 11.8145, 0.00097744, -3.35913e-07]
t_range = [80.0, 1500.0]
units = "J/kmol/K"

[groups.unifac_vle]
"1" = 2
"2" = 2
```

The step 1 spike's hand-written file,
[spike/data_format/106-97-8.toml](../spike/data_format/106-97-8.toml), shows
a whole compound. Writing the importer (step 3) changed three things since
then:
- `name` is now the source's own name;
- `ids.chemsep_name` is gone;
- Mathias–Copeman is three keys rather than an array.

The tables below are current.

The field names of the structs in
[components/record.v](../components/record.v) are these keys, and the
writer follows the fields' order.

## Top level

| Key | Type | Meaning |
|---|---|---|
| `format` | integer | the format version, 1 |
| `name` | string | the compound's name, as the source spells it (ChemSep's are capitalised oddly, like `N-butane`; lookup ignores case) |
| `source` | string | where every value in the file comes from, unless a value says otherwise |

## `[ids]`

| Key | Type | Notes |
|---|---|---|
| `cas` | string | required; must match the file name. A valid CAS number, check digit included. ChemSep adds a letter for two isomers it keeps apart, `1333-74-0p` and `1333-74-0o` (para- and ortho-hydrogen). |
| `inchikey` | string | when known |
| `pubchem` | integer | PubChem CID, when known |
| `smiles` | string | as given by the source, with the importer's repairs (for example a stray trailing `<` removed) |
| `structure` | string | ChemSep's structural formula, such as `CH3(CH2)2CH3` |
| `chemsep_index`, `family` | integer | ChemSep's own record index and family number |
| `aliases` | array of strings | other names to look the compound up by |

## `[constants]`: SI per mol

Each value is a number (integer or float), or an inline table
`{ value = …, source = "…" }` when it comes from somewhere other than the
file's `source`.

| Key | Unit | | Key | Unit |
|---|---|---|---|---|
| `tc` | K | | `dhf` | J/mol (ideal gas, 298.15 K) |
| `pc` | Pa | | `dgf` | J/mol (ideal gas, 298.15 K) |
| `vc` | m³/mol | | `s0` | J/(mol·K) (ideal gas, 298.15 K, 1 bar) |
| `zc`, `omega` | — | | `dh_fusion`, `dh_combustion` | J/mol |
| `mw` | kg/mol | | `v_liquid_tb` | m³/mol |
| `tb`, `tm`, `t_triple` | K | | `p_triple` | Pa |
| `dipole` | C·m | | `radius_of_gyration` | m |
| `solubility_parameter` | J^0.5/m^1.5 | | `charge` | — (an integer) |

The importer converts ChemSep's per-kmol values by moving the decimal
point. It parses "58.1222" written as "58.1222e-3", so the result is
exactly the double nearest 0.0581222 kg/mol. Dividing by 1000 instead
would round twice: for 220 of ChemSep's 2,637 distinct per-kmol values it
lands one unit in the last place away. Vc = 0.0417 m³/kmol, for example,
would become 4.1700000000000004e-05.

## `[parameters]`

Parameters of particular models, in SI where they have units:

| Key | Unit | Meaning |
|---|---|---|
| `rackett` | — | Rackett parameter |
| `costald_volume`, `costald_omega` | m³/mol, — | COSTALD's V* and ω (ChemSep calls the latter the SRK acentric factor) |
| `wilson_volume` | m³/mol | Wilson liquid volume |
| `uniquac_r`, `uniquac_q`, `uniquac_q_prime` | — | UNIQUAC r, q and q′ |
| `vdw_volume`, `vdw_area` | m³/mol, m²/mol | van der Waals (Bondi) volume and area |
| `lj_sigma`, `lj_epsilon_k` | m, K | Lennard-Jones diameter and energy/k |
| `fuller_volume` | — | Fuller diffusion volume |
| `parachor` | kg^0.25·m³/(s^0.5·mol) | parachor |
| `specific_gravity` | — | |
| `mathias_copeman_c1`, `_c2`, `_c3` | — | Mathias–Copeman α constants. The source doesn't say which EOS they were fitted for, so reflux doesn't use them by default. They're separate keys because five ChemSep compounds have only some of them. |
| `api_srk_s1`, `api_srk_s2` | — | API SRK parameters |
| `chao_seader_omega`, `chao_seader_solubility_parameter`, `chao_seader_liquid_volume` | —, J^0.5/m^1.5, m³/mol | Chao–Seader's own constants |

## `[correlations.<property>]`

| Key | Type | Meaning |
|---|---|---|
| `eq` | integer | the equation number, ChemSep's numbering (1–17, 45; DIPPR 100–127) |
| `coeffs` | array of numbers | A, B, C, … **verbatim from the source** |
| `t_range` | array of 2 numbers | [Tmin, Tmax] in K, as the source gives it. 42 ChemSep records have impossible ranges, such as chloroform's solid Cp over [283, 0]. They're kept as they are and listed in the import report. |
| `units` | string | the unit of the correlation's output, as the source gives it |
| `tc` | number | optional; Tc for forms that use reduced temperature, when it differs from `constants.tc` |
| `source` | string | optional, overrides the file's `source` |

Coefficients are never converted. reflux evaluates the correlation in its
native units and converts the **result** to SI per mol. A form like
A + exp(B/T + C + …) can't be rescaled coefficient by coefficient, and
verbatim coefficients can be compared directly with ChemSep and with the
`chemicals` library.

Property names, in file order: `vapor_pressure`, `antoine_vapor_pressure`,
`ideal_gas_cp`, `rpp_cp` (ideal-gas Cp from Reid, Prausnitz and Poling),
`liquid_cp`, `solid_cp`, `liquid_density`, `solid_density`,
`heat_of_vaporization`, `second_virial`, `liquid_viscosity`,
`liquid_viscosity_rps`, `vapor_viscosity`, `liquid_thermal_conductivity`,
`vapor_thermal_conductivity`, `surface_tension`, `relative_permittivity`.

## `[groups.<method>]`

Group counts for group-contribution methods: `unifac_vle`, `unifac_lle`,
`modified_unifac` (Dortmund), `asog`, `umr`, and `gc_method` (PPR78). The
keys are the method's group numbers as quoted strings, and the values are
integer counts. They're used from M6.

Six ChemSep compounds list an ASOG group twice, such as diglyme's ether
group 19 as 2 and then 1. The two always add up to the right count
(diglyme has three ether oxygens), so the importer sums them and reports
each sum.

## `data/components/index.toml`

The importer also writes one index file listing every compound's CAS
number, name, aliases and InChIKey. Lookup by name reads only the index,
and only the requested components' files are parsed. That matters at about
1 ms per file under tcc (below).

## Rules the reader enforces

- `format` must be a version it knows.
- **Types are checked.** V's `toml.Any.f64()` returns 0.0 for a string
  without complaint, so every value's type is checked before use. A wrong
  type is an error naming the file, the key and the type found.
- Unknown keys are errors, so a typo like `omgea` can't silently drop a
  value.
- `ids.cas` must match the file name.

## Rules the writer follows

Output is deterministic, so regenerating the files gives no diff
(M2 exit criterion 2):

- keys always in the order of the fields in
  [components/record.v](../components/record.v);
- floats printed with V's `f64.str()`, which the spike showed round-trips
  exactly and always reads back as a TOML float (`2.0`, not `2`);
- LF line endings and a final newline.

## Why TOML: the step 1 spike

[spike/data_format/spike.v](../spike/data_format/spike.v) runs these checks
against V's `toml` module:

| Check | Result |
|---|---|
| Nested tables, arrays, inline tables, quoted numeric keys, comments | all read correctly; tables iterate in file order |
| A constant given either as a number or as `{ value, source }` | both read |
| Missing key | `value()` gives `Null`; `value_opt()` gives an error |
| Malformed file | error with line and column |
| Every distinct number string in the ChemSep XML (53,977), read by the TOML parser and by `strconv.atof64` | 0 mismatches |
| `f64.str()` output read back | identical bits, always a TOML float |
| Speed | 642 parses of one component file in 641 ms (tcc, debug build) |

JSON via `json2` was the fallback; nothing required it. TOML wins on what
matters for a hand-edited, reviewed data file: comments, and readable
tables.
