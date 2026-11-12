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

## Reactions (format 2)

A file with reactions is format 2: `[definitions]`, `[site.NAME]` and
`[reaction.NAME]`, beside `[components]`. `reflux kinetics check` and
`reflux kinetics explain` read them (M9 step 6); the reactors that use them
in a flowsheet come with M9 step 9. A file of reactions alone has `format`,
`[components]` and its reactions, and a whole flowsheet's other sections
are left to `reflux solve`.
[testdata/kinetics/](../testdata/kinetics/xu_froment.toml) has three whole
files: Xu and Froment's reforming kinetics, Vanden Bussche and Froment's
methanol kinetics, and a DWSIM case's rates as DWSIM has them typed.

```toml
format = 2

[components]
names = ["methane", "water", "hydrogen", "carbon monoxide", "carbon dioxide"]

[site.Ni]
symbol = "DEN"
terms = [
  { name = "K_CO", K = { A = "8.23e-5 bar^-1", dH = "-70.65 kJ/mol" }, a = { CO = 1 } },
  { name = "K_H2O", K = { A = "1.77e5", dH = "88.68 kJ/mol" }, a = { H2O = 1, H2 = -1 } },
]

[reaction.reforming]
equation = "CH4 + H2O = CO + 3 H2"
law = "lhhw"
k = { A = "4.225e15 kmol*bar^0.5/(kg*h)", Ea = "240.1 kJ/mol" }
orders = { CH4 = 1, H2O = 1, H2 = -2.5 }
K = "gibbs"
driving = "difference"
site = "Ni"
site_power = 2
```

Every number has its unit, and every expression is reflux's expression
language ([design.md](design.md#expressions)): parsed, its names bound and
its units checked, never evaluated as text.

### Constants

A constant is an expression, as a string, or a table:

| Table | Constant |
|---|---|
| `{ A, Ea }` | a rate constant, A·exp(−Eₐ/RT) |
| `{ k_ref, Ea, T_ref }` | the same, centred: k_ref·exp(−Eₐ/R·(1/T − 1/T_ref)) |
| `{ A, dH }` | an adsorption constant, A·exp(−ΔH/RT) |
| `{ K_ref, dH, T_ref }` | the same, centred |

### `[definitions]`

Named expressions, `name = "text"`, which rates and constants may use. They
may use one another, in any order, but not themselves.

### `[site.NAME]`

A catalyst's site: its adsorption term, θ = 1 + Σ K·Π aᵢ^eᵢ, which every
reaction on it shares.

| Key | |
|---|---|
| `terms` | Required. A list of terms, each `{ name, K, a }`: `K` a constant, `a` a table of species and their exponents, and `name` what the constant is called (K_ and the first species' name if not given). |
| `symbol` | What θ is called in a rate: `theta_` and the site's name if not given. |

### `[reaction.NAME]`

| Key | |
|---|---|
| `equation` | Required. `CH4 + H2O = CO + 3 H2`: species by name, formula or CAS number, a coefficient whole, decimal or a fraction, written with `=` (an arrow is refused: whether a reaction runs both ways is its rate law's to say). It must balance every element, and the charge; one that doesn't is refused with each imbalance, and the balanced equation where its species allow just one. |
| `law` | `"power_law"`, `"lhhw"`, `"eley_rideal"` or `"expression"`. A reaction without one has no rate, and none of the keys below. |
| `k` | The rate constant: a constant, as above. |
| `k_name` | What it's called: k and the reaction's number if not given. |
| `orders` | Each species' order, a table: `{ CH4 = 1, H2 = -2.5 }`. Not for `eley_rideal`. |
| `K` | Required with a law but `expression`: `"gibbs"`, K at 1 bar from the Gibbs energies of the component data (each species' ideal-gas H and S along Cp); `"none"`, an irreversible rate; or the paper's correlation, in the rate's activities (bar² for steam reforming on partial pressures in bar). |
| `K_name` | What a correlation is called: K and the reaction's number if not given. |
| `basis` | `"partial_pressure"` by default, `"fugacity"`, `"concentration"` or `"mole_fraction"`. |
| `per` | `"catalyst_mass"` by default, or `"volume"`: the rate must be in mol/(kg·s) or mol/(m³·s), in whatever units it's written. |
| `driving` | `"approach"` by default, k·Π aᵢ^αᵢ·(1 − Q/K), or `"difference"`, k·(Π aᵢ^αᵢ − Π aᵢ^(αᵢ+νᵢ)/K): the same function, each written as its paper prints it. |
| `site`, `site_power` | `lhhw` and `eley_rideal`: the site, and its term's power, 1 by default. |
| `terms`, `symbol` | `lhhw`: the reaction's own adsorption term in place of its site's, as the site's `terms`, and what it's called. A site's term belongs to every reaction on it, so `check` notes one. |
| `adsorbed`, `gas` | `eley_rideal`: the adsorbed species, which has a term of its own on the site, K·a¹, and the one that comes from the gas. |
| `rate` | `expression`: the rate, in T, P, each species' `p[...]`, `f[...]`, `c[...]` and `y[...]`, each reaction's `Keq[...]` (K at 1 bar from the Gibbs energies), the definitions and the built-in constants. |

### `reflux kinetics check` and `explain`

```
reflux kinetics check FILE [T P NAME=y ...]
reflux kinetics explain FILE T P NAME=y ...
```

`check` reads the file, which refuses what's wrong (an equation that
doesn't balance, a rate in the wrong units, a name that isn't there), and
prints what else it finds, each a note or a warning: reactions that combine
others, a reaction whose ΔG°f disagree with its ΔH°f and S°, a fitted K's
distance from the Gibbs energies', a negative activation energy, an
adsorption enthalpy above zero, a constant at or below zero, a reaction's
own adsorption term, and a number typed where a built-in constant was
meant, 8.314 J/(mol·K) for R. Given a state, a fitted K is judged at its T,
each typed constant says how far it moves the rates there, and each
reversible rate is held to the equilibrium the reactions reach from it, in
the ideal gas: where along its reaction it's zero, which is at equilibrium
for a rate consistent with the Gibbs energies, and whether it runs towards
equilibrium at the state. `explain`
prints each rate as reflux read it, and its factor table at the state:
every constant, each term of its site's, Q, K and Q/K, with their values,
and the rate. A state is a temperature and a pressure with their units,
then each species' mole fraction, summing to 1:

```
reflux kinetics explain plant.toml "900 K" "15 bar" CH4=0.2 H2O=0.45 H2=0.2 CO=0.1 CO2=0.05
```

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
