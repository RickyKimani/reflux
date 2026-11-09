# The flowsheet file format

A flowsheet is a [TOML](https://toml.io) file. This is format 1.
[testdata/flowsheets/reference.toml](../testdata/flowsheets/reference.toml) is
a complete example:

```toml
format = 1

[components]
names = ["methane", "ethane", "propane", "n-butane", "n-pentane"]

[model]
equation = "PR"

[solver]
method = "wegstein"

[feed.feed]
T = "300 K"
P = "30 bar"
flow = "100 mol/s"
mole_fractions = { methane = 0.7, ethane = 0.1, propane = 0.08, n-butane = 0.07, n-pentane = 0.05 }

[unit.mixer]
type = "mixer"
in = ["feed", "recycle"]
out = "mixed"

[unit.cooler]
type = "heater"
in = "mixed"
out = "cold"
outlet_T = "230 K"
```

`flowsheet.load(path)` reads a file and checks all of it. Anything wrong is
refused with its line, its key, a reason, and suggestions for a misspelt name:

```
flowsheet: line 29: unit.cooler.type: no unit type heatr (did you mean heater?)
```

A key the format doesn't have is an error, not ignored, so a misspelling can't
silently leave a default in place.

`reflux solve plant.toml` solves a file and prints its stream table, and
`--json` and `--csv` write its results: [report-format.md](report-format.md)
describes them.

## Numbers and units

Every physical number is a string with its unit, such as `"230 K"` or
`"100 kmol/h"`. A bare number there is refused, because that's how a pressure in
bar gets read as pascals.

| Kind | Units |
|---|---|
| temperature | `K`, `degC`, `°C`, `degF`, `°F`, `R` |
| temperature difference | the same, with no offset: `"9 degF"` is 5 K |
| pressure | `Pa`, `kPa`, `MPa`, `bar`, `mbar`, `atm`, `psi`, `psia`, `mmHg`, `torr`, and the gauge `barg` and `psig`, above 1 atm |
| pressure difference | as pressure, without the gauge units |
| molar flow | `mol/s`, `mol/min`, `mol/h`, `kmol/s`, `kmol/min`, `kmol/h`, `lbmol/h` |
| mass flow | `kg/s`, `kg/min`, `kg/h`, `g/s`, `t/h`, `lb/s`, `lb/h` |
| power and duty | `W`, `kW`, `MW`, `kJ/s`, `kJ/h`, `MJ/h`, `GJ/h`, `kcal/h`, `Btu/h`, `hp` |
| molar enthalpy | `J/mol`, `kJ/mol`, `kJ/kmol`, `MJ/kmol`, `Btu/lbmol` |

A fraction, such as an efficiency, can be a plain number or a percentage: `0.75`
or `"75 %"`.

These are each kind's common units. Any unit of the right dimension reads too,
prefixed or compound, such as `"10 lbmol/min"`, `"3 kJ/s"` or `"5 m^3/h"`, with
no space inside the unit. A unit of another dimension is refused with what it's
a unit of. [design.md](design.md#units-of-measure) lists the symbols and
prefixes, and how a compound unit is written.

## `format`

Required, and 1. A later format will come with a way to migrate.

## `[components]`

| Key | |
|---|---|
| `names` | Required. A list of compounds, by name, alias or CAS number. |
| `data` | A list of directories of your own component files, which are read after the data root's compounds and override them. A relative path is from the flowsheet file. The data root is `$REFLUX_DATA`, or the source tree's `data/` ([data-format.md](data-format.md#your-own-components)). |

[data-format.md](data-format.md) describes a component file.

## `[model]`

| Key | |
|---|---|
| `equation` | Required: `"PR"`, `"SRK"`, `"RK"` or `"vdW"`. |
| `kij` | `true` by default: binary interaction parameters from ChemSep's bank. |
| `alpha` | `"equation"` by default, or `"mathias_copeman"`, for PR only. |
| `kij_overrides` | A list of pairs set by hand: `[{ a = "methane", b = "ethane", k = 0.01 }]`. |

## `[solver]`

Optional, as is each key.

| Key | |
|---|---|
| `method` | `"wegstein"` by default, `"direct"` or `"broyden"`. |
| `tol` | `1e-9` by default: how far a tear variable may move in a pass, as a fraction of its scale. |
| `max_iterations` | 100 by default. |
| `tears` | A list of streams to tear, whatever the search would choose. Others are added if the loops need them. |

## `[feed.NAME]`

A stream entering the flowsheet. It takes two state keys and a composition.

**State,** one of:

| Keys | |
|---|---|
| `T`, `P` | |
| `P`, `H` | `H` is the molar enthalpy, on the fluid's reference state. |
| `T`, `vapour_fraction` | |
| `P`, `vapour_fraction` | |

**Composition,** one of:

| Keys | |
|---|---|
| `flow`, `mole_fractions` | `flow` is a molar or a mass flow. |
| `flow`, `mass_fractions` | The same. |
| `molar_flows` | Each component's flow. |
| `mass_flows` | Each component's flow. |

Each is a table by component, and a component left out is 0: `mole_fractions = {
methane = 0.75, n-pentane = 0.25 }`. Fractions must sum to 1 within 1e-6. Mass
bases are converted to moles with the components' molar masses.

## `[estimate.NAME]`

A first value for a stream, with the same keys as a feed. The stream is then
torn, so it has to be in a loop. A tear needs one when its first pass would
otherwise fail, as when a gas guess would reach a pump.

## `[unit.NAME]`

Every unit has a `type` and its streams, and exactly one of each kind of
specification. A stream is named at both ends, on the outlet that makes it and
the inlet that takes it, and a stream that no inlet takes is a product.

`in` and `out` are a stream or a list, in port order. Two-outlet drums name their
outlets instead.

### `mixer`

Any number of inlets, at least 2.

| Key | |
|---|---|
| `in`, `out` | A list of inlets, and one outlet. |
| `P` | The outlet's pressure, at most the lowest inlet's. By default, the lowest inlet's. |

### `splitter`

| Key | |
|---|---|
| `in`, `out` | One inlet, and a list of outlets. |
| `fractions` | One for each outlet, summing to 1. |

### `valve`

One of `outlet_P` or `delta_P`, the drop. Isenthalpic.

### `heater`

For cooling too: a duty or `delta_T` below 0. One of:

| Key | |
|---|---|
| `outlet_T` | |
| `delta_T` | |
| `duty` | Heat added. |
| `outlet_vapour_fraction` | |

and optionally `pressure_drop`.

### `pump`, `compressor`, `expander`

One of `outlet_P` or `delta_P`, and `efficiency`, which is required.

A pump takes liquid only. Its work is v·ΔP/η, with v the inlet's. A compressor or
expander finds its isentropic outlet, then applies the efficiency.

### `flash_drum`

| Key | |
|---|---|
| `in` | One inlet, or a list. |
| `vapour`, `liquid` | The two outlets. With one phase, the other is empty. |
| `T` or `duty` | One of these. A duty of `"0 W"` is adiabatic. |
| `P` | The drum's pressure, at most the lowest inlet's. By default, the lowest inlet's. |

### `component_separator`

| Key | |
|---|---|
| `in` | One inlet. |
| `out` | Two outlets. |
| `fractions` | A table by component: how much of each goes to the first outlet. A component left out all goes to the second. |
| `T`, `P` | Both outlets' state. |

## What it refuses that TOML allows

- **A `\u` escape.** Write the character itself. V 0.5.2's TOML reader can crash
  on one, and nothing can catch that, so reflux refuses them first.
- **Nesting past 32 levels.** No flowsheet comes near it, and it's what would
  overflow the stack.
- **An escaped quote in a quoted table name,** such as `[unit."a\"b"]`: V's reader
  refuses that.

## Writing

`Loaded.write()` gives the file back in a fixed order, with its numbers and units
as they were written, and loading that gives the same flowsheet. Comments are
lost, because the reader doesn't tie them to anything.
