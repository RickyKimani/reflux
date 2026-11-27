# The report format

`reflux solve plant.toml --json results.json` writes a flowsheet's results as
JSON, and `--csv streams.csv` its streams as CSV. This is format 4 of the
JSON, and of the CSV, which is format 1's. It changes only with a new
`format` number, so that what reads a report can rely on it. Format 2 (M9)
adds a reactor's profile to its unit, a unit's warnings, an element balance
where something reacts, and the failure `"reaction_failed"`. Format 3 (M8)
adds the file's specs, each block's specs and their residuals, and the
failure `"spec_not_met"`. Format 4 (M6) adds `model`, how the fluid was
built, and a phase `"liquid2"`, so that a stream may have three phases.

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
| `format` | 4. |
| `reflux` | The version of reflux that wrote it. |
| `file` | The flowsheet file, as the command line named it. |
| `status` | `"solved"` or `"failed"`. |
| `error` | Why it failed, or `null`. |
| `components` | The components. Empty if the file didn't load. |
| `model` | How the fluid was built, or `null` if the file didn't load (format 4). |
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

### `model`

How the file's fluid was built (M6's D13): an equation of state, or a γ–φ
fluid, with each pair's parameters and where they came from.

| Key | |
|---|---|
| `equation` | The equation of state, `"PR"` and so on, or `""` for a γ–φ fluid. |
| `liquid` | A γ–φ fluid's activity model, `"NRTL"` and so on, or `""`. |
| `vapour` | A γ–φ fluid's vapour, `"ideal gas"` or an equation, or `""`. |
| `poynting` | Whether a γ–φ liquid has the Poynting correction; `false` for an equation of state. |
| `solutes` | A γ–φ fluid's Henry's-law solutes, by the file's names. |
| `pairs` | Each pair of components the model has parameters for, in the file's order: an equation of state's kᵢⱼ, every pair; a binary activity model's, every pair of solvents. UNIFAC has none, its parameters being its groups'. |
| `notes` | Where other values came from: an α not the equation's own, a vapour pressure from its Antoine record, a gas's Henry's constant from IAPWS. |
| `warnings` | What the answers should be read with: a pair with no parameters, made ideal by `missing = "ideal"`. `reflux solve` prints them. |

### `model.pairs[]`

| Key | |
|---|---|
| `a` | The pair's first component, by the file's name: its parameters with `b` come first. |
| `b` | The second. |
| `parameters` | Each as `{ name, text }`: `kij`, a number; `tau_ab`, `tau_ba` and NRTL's `alpha`, or Wilson's `lambda_ab` and `lambda_ba`, each an expression in T as the file would write it. A bank's is in its own form: NRTL's `"792.8020 K / (1.98721*T)"` is A₁₂ in cal/mol over DECHEMA's R. |
| `source` | Where they came from: a bank's line and its comment, such as `"nrtl.ipd line 180: Methanol/Water p61 1/1a"`, the DECHEMA volume and page; `"the file"`; `"none listed: 0"` for a kᵢⱼ; or `"none listed: an ideal pair"`, `"none listed: no residual interaction"` for UNIQUAC. |

### `streams[]`

| Key | |
|---|---|
| `name` | |
| `from` | The unit it leaves, or `null` for a feed. |
| `to` | The unit it enters, or `null` for a product. |
| `T` | K. |
| `P` | Pa. |
| `vapour_fraction` | The vapour's share of the moles, 0 with no vapour: a liquid, or two (M6's D9). |
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
| `phase` | `"vapour"`, `"liquid"` or `"liquid2"`, as the flash labels it (M6's D9): by the phase identification parameter for an equation of state, at most one vapour, and of two liquids the denser `"liquid2"`. A stream may have up to three phases. Report format 4 (M6 step 9) makes `"liquid2"` its own. |
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
| `warnings` | What its results should be read with, though they stand, each a sentence: a conversion reactor's outlet past equilibrium, say, with its Q/K, or a flash drum's liquid outlet of two liquids. Empty for most units. The text report lists them under the unit table. |

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
γ–φ fluid's model and its warnings, a stream table, and each unit's duty and
work, in K, bar, mol/s, kg/s, kJ/mol
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
| `sensitivities` | The derivatives each point gives, each input's every output, in order: empty for a study without. |
| `points` | Every point, in the grid's order: the first input's values outer. |

#### `sensitivities[]`

| Key | |
|---|---|
| `output` | The output's name. |
| `input` | The key it's with respect to. |
| `unit` | The derivative's SI unit, the output's over the input's, `""` for a pure number. |

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
| `derivatives` | Each of `sensitivities` at this point, in its order. |

#### `points[].derivatives[]`

| Key | |
|---|---|
| `value` | The derivative, by central differences of solves at 1e-12 over steps of 1e-4 and 5e-5 of the input, Richardson's extrapolation of the two; or `null`. |
| `error` | Its error's estimate: the two differences' gap, and the solves' own noise over the step from a fourth difference; or `null`. |
| `reason` | Why it has no value, or `""`. |

The CSV has a row for each point: `point`, each input as `KEY [unit]`,
each output as `NAME [unit]`, NaN where it has no value, each derivative as
`d(NAME)/d(KEY) [unit]` and its error as `error of d(NAME)/d(KEY) [unit]`,
and `status`.
`reflux study` prints the points to 5 significant figures, each output's
reasons once, and each failed point's. `--workers N` runs it on N threads,
as many as the processors by default; the report is the same bits whatever
N is.

## The optimisation report

`reflux optimize plant.toml best --json optimum.json --csv iterations.csv`
writes an optimisation's outcome (M8's D13, D17): a format of its own,
format 1, which changes only with a new number. Every number is SI. An
optimisation with no optimum still writes its report, with its reason and
its iterates.

#### The optimisation

| Key | |
|---|---|
| `format` | 1. |
| `reflux` | The version of reflux that wrote it. |
| `file` | The flowsheet file, as the command line named it. |
| `optimize` | The optimisation's name. |
| `status` | `"optimal"`, with an optimum whose first-order conditions hold; or `"failed"`, with `reason`, or `error` if the file didn't load. |
| `error` | Why the file didn't load, as the solve's `error`, or `null`. |
| `reason` | Why there's no optimum, or `""`. |
| `goal` | `"minimize"` or `"maximize"`. |
| `objective` | What's minimised or maximised: `expression` as reflux writes it back, `unit`, and `value` at the optimum, or `null`. |
| `inputs` | What it varies, each as below, in the file's order. |
| `constraints` | What it's subject to, each as below, in the file's order. |
| `specs` | The file's specs at the optimum: each `name` and `value`. |
| `iterations` | Every iterate, the last the optimum, or the last before it failed. |
| `evaluations` | The solves at points the optimiser chose. |
| `solves` | Every solve, the gradients' included. |

#### `inputs[]`

| Key | |
|---|---|
| `key` | The input's key, as `vary` gives it. |
| `unit` | Its values' SI unit, `""` for a pure number. |
| `lower` | Its least value. |
| `upper` | Its most. |
| `start` | Where the optimisation started it. |
| `value` | Its value at the optimum, or at the last iterate. |
| `bound` | `"lower"` or `"upper"` where it's on a bound, else `""`. |
| `multiplier` | The bound's multiplier: how fast the objective improves (falls when minimised, rises when maximised) as that bound is relaxed, per unit of the input; 0 off a bound. |

#### `constraints[]`

| Key | |
|---|---|
| `equation` | The constraint, as reflux writes it back. |
| `unit` | Its sides' SI unit, `""` for a pure number. |
| `lhs` | Its left side at the optimum, or `null`. |
| `rhs` | Its right side, or `null`. |
| `active` | Whether it holds with equality there: always, for an equality. |
| `multiplier` | How fast the objective improves as the constraint is relaxed, per unit of its sides, at least 0 for an inequality; for an equality, as its right side grows, of either sign. `null` with no optimum. |

#### `iterations[]`

| Key | |
|---|---|
| `inputs` | Each input's value. |
| `objective` | The objective's value. |
| `violation` | The worst constraint's violation, as a fraction of its sides' size at the start. |
| `stationarity` | The Lagrangian's gradient's largest entry, in the inputs scaled over their bounds and the objective over its largest derivative at the start: below 1e-6 at an optimum. |
| `step` | The fraction of the QP's step the line search took from it, 0 at the last. |

The CSV has a row for each iterate: `iteration` from 0, each input as
`KEY [unit]`, the objective as `minimize [unit]` or `maximize [unit]`, and
`violation`, `stationarity` and `step`.
`reflux optimize` prints the optimum, each input with its bounds and any
bound it's on, each constraint, and each iterate, to 5 significant figures.
`--workers N` differentiates by N inputs at once, as many as the processors
by default; the report is the same bits whatever N is.

## The fit report

`reflux fit reactions.toml rates --json fit.json --csv residuals.csv`
writes a fit's outcome (M8's D14, D15, D17): a format of its own, format 1,
which changes only with a new number. Every number is SI, each parameter
in the file's own form: an Arrhenius constant fitted centred on the data's
mean T is reported back as A and Eₐ.

#### The fit

| Key | |
|---|---|
| `format` | 1. |
| `reflux` | The version of reflux that wrote it. |
| `file` | The file, as the command line named it. |
| `fit` | The fit's name. |
| `status` | `"fitted"`; or `"failed"`, with `reason`, or `error` if the file didn't load. |
| `error` | Why the file didn't load, as the solve's `error`, or `null`. |
| `reason` | Why there's no fit, or `""`. |
| `kind` | `"rates"`, rate constants to rates, or `"vle"`, kᵢⱼ or an activity model's pairs to bubble points. |
| `data` | The data's file, as the fit names it. |
| `parameters` | Each parameter, as below, in the fit's order. |
| `correlation` | The parameters' correlations, a row each, in the file's form. |
| `centred` | Each constant fitted centred: `constant` and `t_ref`, the data's mean T. |
| `measurements` | The residuals: each measured rate, or each bubble point's P and measured y. |
| `dof` | The measurements less the parameters. |
| `s` | The weighted residuals' standard deviation, √(Σr²/dof), about 1 where the weights are the measurements' errors; or `null`. |
| `t` | Student's t at 0.975 on `dof`, each interval's half-width over its standard error; or `null`. |
| `residuals` | Each measurement against the model at the answer, as below. |
| `section` | The fitted constants as the file writes them, in its units, to 10 significant figures. |
| `iterations` | Levenberg–Marquardt's iterations. |
| `evaluations` | Its evaluations of the residuals. |

#### `parameters[]`

| Key | |
|---|---|
| `name` | As the fit names it: `k1.A`, `K_E.dH`, `K_H` or `kij[ethane, propane]`. |
| `unit` | Its SI unit, `""` for a pure number. |
| `file_unit` | The unit the file writes it in, `""` for a pure number. |
| `start` | Where the fit started it. |
| `value` | Its fitted value, or `null`. |
| `standard_error` | Its standard error, by linearisation, s²·(JᵀJ)⁻¹, or `null`. |
| `lower` | Its 95% interval's lower end, value − t·standard error, or `null`. |
| `upper` | And its upper end, or `null`. |

#### `residuals[]`

| Key | |
|---|---|
| `row` | The measurement's row in the data, from 0. |
| `column` | Its column: `r[hydrogenation]`, `P` or `y[ethane]`. |
| `measured` | Its measured value. |
| `model` | The model's at the answer. |
| `weighted` | Their difference, model less measured, over the measurement's error. |

The CSV has a row for each residual: `row`, `column`, `measured`, `model`
and `weighted`. `reflux fit` prints each parameter in the file's unit with
its standard error and interval, their correlations, and the fitted
constants as the file writes them.

## The water report

`reflux water T=450degC P=40bar --json water.json` writes one state of water
(M5a's D13): a format of its own, format 1, which changes only with a new
number. Every number is SI, per kg on IAPWS's reference (u = s = 0 for the
liquid at the triple point), and per mole by IAPWS's molar mass,
18.015268 g/mol. A property a state hasn't (a mixture's cp, a viscosity
outside its correlation's range) is left out.

| Key | |
|---|---|
| `format` | 1. |
| `reflux` | The version of reflux that wrote it. |
| `formulation` | `"IAPWS-95"` or `"IAPWS-IF97"`. |
| `reference` | The reference of h, s and u, as a sentence. |
| `given` | The two SPECs, as typed: `name` (`T`, `P`, `h`, `s`, `x` or `rho`), `text` (the value as typed) and `value` (in SI, per kg). |
| `T`, `P` | K, Pa. |
| `phase` | `"liquid"`, `"vapour"`, `"supercritical"`, `"saturated liquid"`, `"saturated vapour"`, or `"liquid and vapour"` for a wet state. |
| `x` | The vapour's mass fraction, at a saturated or wet state; else left out. |
| `region` | IF97's region (4 for a saturated or wet state); left out by IAPWS-95. |
| `per_kg` | The state's properties per kg, as below; a wet state's are the mixture's. |
| `per_mole` | The same per mole: `rho` (mol/m³), `v` (m³/mol), `h`, `s`, `u`, `cv`, `cp`; no `w`, `mu` or `lambda`. |
| `liquid`, `vapour` | A wet state's saturated phases, each as `per_kg`; else left out. |
| `sigma` | The surface tension at T, N/m, at a saturated or wet state; else left out. |
| `notes` | What the state should be read with: a metastable state below the triple point, a transport property's range. |

#### `per_kg`, `liquid`, `vapour`

| Key | |
|---|---|
| `rho` | kg/m³. |
| `v` | m³/kg. |
| `h`, `u` | J/kg. |
| `s`, `cv`, `cp` | J/(kg·K). |
| `w` | The speed of sound, m/s. |
| `mu` | The viscosity, Pa·s (R12-08; on IF97 its industrial form). |
| `lambda` | The thermal conductivity, W/(m·K) (R15-11; on IF97 its industrial form). |

`reflux water` prints the same, rounded to seven figures: one phase per kg
and per mole, or a wet state's liquid, vapour and mixture per kg and the
mixture per mole, then σ and any notes.

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

`reflux optimize` exits 0 with an optimum; 1 when its file didn't load or
has no optimisations, when it found no optimum, or when a report couldn't be
written; and 2 for a wrong command.

`reflux fit` exits 0 with a fit; 1 when its file didn't load or has no
fits, when the fit failed, or when a report couldn't be written; and 2 for a
wrong command, `--workers` included.

`reflux water` exits 0 with a state; 1 when water has no state there, or
the formulation none (outside its range, above the melting pressure of
ice, no saturation above the critical point), or a report couldn't be
written; and 2 for a wrong command, a pair of SPECs that doesn't fix one
state included.
