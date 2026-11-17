# The report format

`reflux solve plant.toml --json results.json` writes a flowsheet's results as
JSON, and `--csv streams.csv` its streams as CSV. This is format 3 of the
JSON, and of the CSV, which is format 1's. It changes only with a new
`format` number, so that what reads a report can rely on it. Format 2 (M9)
adds a reactor's profile to its unit, a unit's warnings, an element balance
where something reacts, and the failure `"reaction_failed"`. Format 3 (M8)
adds the file's specs, each block's specs and their residuals, and the
failure `"spec_not_met"`.

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
| `format` | 3. |
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
| `specs` | Each spec of the file as it was met, in the file's order. Empty without specs, and on a failure. |

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
| `profile` | A PFR's state along it, or `null` for a unit without one. |
| `warnings` | What its results should be read with, though they stand, each a sentence: a conversion reactor's outlet past equilibrium, say, with its Q/K. Empty for most units. The text report lists them under the unit table. |

### `units[].energy[]`

| Key | |
|---|---|
| `name` | `"duty"` for heat, `"work"` for shaft work. |
| `rate` | W, into the unit's material: a cooler's duty is below 0, and a pump's work above. |

### `units[].profile[]`

The PFR's inlet, then each point its integrator stepped to.

| Key | |
|---|---|
| `volume` | m³ from the inlet. |
| `length` | m from the inlet, for a reactor whose length is given, or 0. |
| `T` | K. |
| `P` | Pa. |
| `molar_flows` | mol/s. |

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
| `specs` | The specs met with the block, by name: those whose links put the block together (D5 in [m8-plan.md](m8-plan.md)). |
| `spec_residuals` | A row for each evaluation of the block's specs, each one's residual as a fraction of its scale (`specs[].residual`). Empty for a block without specs. |

### `specs[]`

| Key | |
|---|---|
| `name` | As the file names it, `[spec.NAME]`. |
| `vary` | The key it varied, such as `"unit.cooler.outlet_T"`. |
| `value` | The value the key took, in SI. |
| `unit` | That value's SI unit, such as `"K"`, or `""` for a pure number. |
| `until` | Its equation, as reflux writes it back. |
| `lhs` | The equation's left side in the solution, in SI. |
| `rhs` | Its right side. |
| `residual` | lhs − rhs as a fraction of the larger of the two, or of the block's own scale for that dimension where both are near zero (a flow's F, a duty's F·RT, a pure number's 1). It's within `tol`, checked from the solution itself. |

### `balance`, `worst_unit`, `units[].balance`

| Key | |
|---|---|
| `component` | The worst component's imbalance, as a fraction of the flow in. Where something reacts, the worst element's, in atoms: for a reactor, and for the whole flowsheet when it has one. |
| `energy` | The energy's imbalance, as a fraction of F·RT, with F the flow in. For a unit, T is its hottest inlet's, and for the whole flowsheet, the feeds' mean, weighted by flow. |

### `error`

| Key | |
|---|---|
| `kind` | `"bad_file"`, `"bad_flowsheet"`, `"bad_specification"`, `"flash_failed"`, `"no_convergence"`, `"unbalanced"`, `"reaction_failed"` (a reactor's equilibrium, integration or Newton) or `"spec_not_met"` (a spec's target out of reach of its bounds, or not met), and `"other"` for anything else. |
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

`reflux solve` prints how the flowsheet was solved, each spec's value, a
stream table, and each unit's duty and work, in K, bar, mol/s, kg/s, kJ/mol
and kW, to 5 significant figures. `--trace` adds every pass of every loop,
and each spec's residual at every evaluation. These are for reading, and
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

## The study report

`reflux study plant.toml --json study.json --csv study.csv` writes a
study's results (M8's D17): a format of its own, format 1, which changes
only with a new number. Every number is SI. A point that fails, or an
output with no value, is reported with its reason, and the study's other
points stand.

#### The study

| Key | |
|---|---|
| `format` | 1. |
| `reflux` | The version of reflux that wrote it. |
| `file` | The flowsheet file, as the command line named it. |
| `study` | The study's name. |
| `status` | `"ran"`, even with failed points, or `"failed"` if the file didn't load. |
| `error` | Why it failed, as the solve's `error`, or `null`. |
| `inputs` | What it varies: `key` and `unit`, the values' SI unit, `""` for a pure number. |
| `outputs` | What it reports: `name`, `expression` as reflux writes it back, and `unit`. |
| `specs` | The file's specs, whose values each point gives. |
| `chunk` | How many points ran from one cold start. |
| `points` | Every point, in the grid's order: the first input's values outer. |

#### `points[]`

| Key | |
|---|---|
| `index` | Its place in the grid's order, from 0. |
| `inputs` | Each input's value. |
| `status` | `"solved"`; `"no_value"`, solved, with an output that has no value; or `"failed"`. |
| `outputs` | Each output, `{ "value": x, "reason": "" }`, or `{ "value": null, "reason": "it divides by …" }`. |
| `error` | Why a failed point failed, or `""`. |
| `specs` | Each spec's value there. |
| `passes` | The passes of its blocks with tears. |
| `spec_evaluations` | Its specs' evaluations. |
| `steps` | Its integrators' accepted steps. |
| `warm` | Whether it started from its neighbour's solution, or cold. |

The CSV has a row for each point: `point`, each input as `KEY [unit]`,
each output as `NAME [unit]`, NaN where it has no value, and `status`.
`reflux study` prints the points to 5 significant figures, each output's
reasons once, and each failed point's. `--workers N` runs it on N threads,
as many as the processors by default; the report is the same bits whatever
N is.

## Exit status

| Status | |
|---|---|
| 0 | It solved, and every report was written. |
| 1 | The flowsheet didn't load or didn't solve, or a report couldn't be written. |
| 2 | The command is wrong: an unknown option, say, or no file. |

`reflux study` exits 0 when the study ran, failed points and all, which its
table lists; 1 when its file didn't load or has no studies, or a report
couldn't be written; and 2 for a wrong command, a study name that isn't
there included.
