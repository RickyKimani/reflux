# The report format

`reflux solve plant.toml --json results.json` writes a flowsheet's results as
JSON, and `--csv streams.csv` its streams as CSV. This is format 1 of both. It
changes only with a new `format` number, so that what reads a report can rely
on it.

- **Every number is SI:** K, Pa, mol/s, kg/s, m³/s, J/mol, J/(mol·K), m³/mol
  and W.
- **A list by component,** such as `molar_flows`, is in the order of
  `components`, which is the order the flowsheet file lists them in.
- **Each number reads back as the same double:** it's the shortest decimal
  that does.

From Python, with pandas:

```python
import json
import pandas as pd

r = json.load(open("results.json"))
names = [c["name"] for c in r["components"]]
streams = pd.DataFrame(r["streams"]).set_index("name")
x = pd.DataFrame(streams["mole_fractions"].tolist(), index=streams.index, columns=names)

table = pd.read_csv("streams.csv", index_col="stream")
```

## The JSON

### The report

| Key | |
|---|---|
| `format` | 1. |
| `reflux` | The version of reflux that wrote it. |
| `file` | The flowsheet file, as the command line named it. |
| `status` | `"solved"` or `"failed"`. |
| `error` | Why it failed, or `null`. |
| `components` | The components. Empty if the file didn't load. |
| `streams` | Every stream: the feeds, then each unit's outlets, unit by unit. Empty on a failure. |
| `units` | Every unit, in the file's order. Empty on a failure. |
| `solver` | What it was solved with, or `null` if the file didn't load. |
| `blocks` | The plan: the parts of the flowsheet in the order they were solved, and how each loop converged. |
| `notes` | Anything the plan has to say about its choices. |
| `balance` | How far the whole flowsheet is from closing its balances, or `null` on a failure. |
| `worst_unit` | The worst unit's, for each balance, or `null` on a failure. |

### `components[]`

| Key | |
|---|---|
| `name` | As the file names it. |
| `compound` | As the database does. |
| `molar_mass` | kg/mol. |

### `streams[]`

| Key | |
|---|---|
| `name` | |
| `from` | The unit it leaves, or `null` for a feed. |
| `to` | The unit it enters, or `null` for a product. |
| `T` | K. |
| `P` | Pa. |
| `vapour_fraction` | The share of the moles in the less dense phase, or with one phase, 1 if it's labelled vapour and 0 if liquid. |
| `molar_flow` | mol/s. |
| `mass_flow` | kg/s. |
| `volume_flow` | m³/s, at the stream's T and P. |
| `H` | Molar enthalpy, J/mol, from the elements at 298.15 K and 1 bar ([design.md](design.md)). |
| `S` | Molar entropy, J/(mol·K), on the same reference. |
| `V` | Molar volume, m³/mol. |
| `molar_flows` | mol/s of each component. |
| `mass_flows` | kg/s of each component. |
| `mole_fractions` | |
| `mass_fractions` | |
| `phases` | One or two, the less dense first. |

An empty stream, such as a drum's liquid when there's no liquid, keeps the
composition it was given. So it has fractions and a phase though its flows are
0.

### `streams[].phases[]`

| Key | |
|---|---|
| `phase` | `"vapour"` or `"liquid"`, as the flash labels it by the phase identification parameter. |
| `fraction` | Its share of the stream's moles. |
| `mole_fractions` | |
| `H` | J/mol of the phase. |
| `S` | J/(mol·K). |
| `V` | m³/mol. |

### `units[]`

| Key | |
|---|---|
| `name` | |
| `type` | As the file gives it: `"mixer"`, `"heater"`, and so on. |
| `inlets` | The streams on its inlets, in port order. |
| `outlets` | The same for its outlets. |
| `energy` | Its duty or work. |
| `balance` | How far it closes its balances, on the inlets it was calculated from. Where an inlet is a tear, those are a pass behind the stream in `streams`, by no more than `tol`. |

### `units[].energy[]`

| Key | |
|---|---|
| `name` | `"duty"` for heat, `"work"` for shaft work. |
| `rate` | W, into the unit's material: a cooler's duty is below 0, and a pump's work above. |

### `solver`

| Key | |
|---|---|
| `method` | `"direct"`, `"wegstein"` or `"broyden"`. |
| `tol` | How far a tear variable may move in a pass, converged, as a fraction of its scale. |
| `max_iterations` | The most passes a loop may take. |

### `blocks[]`

| Key | |
|---|---|
| `units` | Its units, in the order they were calculated. |
| `tears` | The streams torn, or none for a block of one unit. |
| `cycles` | Its elementary cycles. |
| `residuals` | A row for each pass, holding each tear's largest change in that pass as a fraction of its scale: a component's flow against the stream's total flow, P against itself, and H against RT. It converged when the last row is within `tol`. Empty for a block with no tears, which is calculated once, and for one in which a failure came before a pass finished. |

### `balance`, `worst_unit`, `units[].balance`

| Key | |
|---|---|
| `component` | The worst component's imbalance, as a fraction of the flow in. |
| `energy` | The energy's imbalance, as a fraction of F·RT, with F the flow in. For a unit, T is its hottest inlet's, and for the whole flowsheet, the feeds' mean, weighted by flow. |

### `error`

| Key | |
|---|---|
| `kind` | `"bad_file"`, `"bad_flowsheet"`, `"bad_specification"`, `"flash_failed"`, `"no_convergence"` or `"unbalanced"`, and `"other"` for anything else. |
| `message` | All of it, as one line. |
| `unit` | The unit at fault, or `""`. |
| `reason` | Why, in words. |
| `cause` | For `"flash_failed"`, the flash's own failure, such as `"more_phases"`, or else `null`. |
| `line` | The line of the file at fault, or 0. |
| `key` | The key at fault, such as `"unit.cooler.type"`, or `""`. |
| `suggestions` | What a misspelt name may have meant. |

## The CSV

A row for each stream, in the JSON's order, and a column for each of a
stream's values in the JSON but its phases, named by its key with its unit:

    stream, from, to, T [K], P [Pa], vapour_fraction, molar_flow [mol/s],
    mass_flow [kg/s], volume_flow [m3/s], H [J/mol], S [J/(mol K)], V [m3/mol]

then for each component `molar_flows.NAME [mol/s]`, then each
`mass_flows.NAME [kg/s]`, each `mole_fractions.NAME` and each
`mass_fractions.NAME`. A feed's `from` and a product's `to` are empty. A field
with a comma, a quote, a line break or a leading space is quoted, as RFC 4180
has it.

## What it prints

`reflux solve` prints how the flowsheet was solved, a stream table, and each
unit's duty and work, in K, bar, mol/s, kg/s, kJ/mol and kW, to 5 significant
figures. `--trace` adds every pass of every loop. These are for reading, and
their layout may change. The JSON is for programs.

## Failures

On a failure the JSON has `status` `"failed"`, the `error`, and the `blocks`
with the passes as far as they got. It has no streams or units, since nothing
that hasn't been shown right is reported (D8 in [m4-plan.md](m4-plan.md)). A
CSV from an earlier run is removed, so that it can't be taken for this one's.
The message goes to stderr, its file and line first, as a compiler's would:

```
reflux: plant.toml:26: unit.cooler.type: no unit type heatr (did you mean heater?)
```

## Exit status

| Status | |
|---|---|
| 0 | It solved, and every report was written. |
| 1 | The flowsheet didn't load or didn't solve, or a report couldn't be written. |
| 2 | The command is wrong: an unknown option, say, or no file. |
