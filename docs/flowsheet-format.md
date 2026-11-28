# The flowsheet file format

A flowsheet is a [TOML](https://toml.io) file, format 1, 2, 3, 4 or 5:
format 2 (M9) adds reactions and reactors, format 3 (M8) adds specs,
studies, optimisations and fits, format 4 (M6) adds γ–φ fluids, whose liquid
is an activity model, and the decanter, and (M5a) water by IAPWS-95, and
format 5 (M5b) adds water by IF97 and `[transport]`, a liquid pair's
viscosity parameter. A file of reactions alone may be format
2 or 3. Each reads every earlier format's sections as that format does.
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

Required: 1, 2 for a file with reactions or reactors, 3 for a file with
specs, 4 for a file with an activity model, a decanter or IAPWS-95, or 5
for a file with IF97 or `[transport]`. A later format will come with a way
to migrate.

## `[components]`

| Key | |
|---|---|
| `names` | Required. A list of compounds, by name, alias or CAS number. |
| `data` | A list of directories of your own component files, which are read after the data root's compounds and override them. A relative path is from the flowsheet file. The data root is `$REFLUX_DATA`, or the source tree's `data/` ([data-format.md](data-format.md#your-own-components)). |

[data-format.md](data-format.md) describes a component file.

## `[model]`

An equation of state:

| Key | |
|---|---|
| `equation` | Required: `"PR"`, `"SRK"`, `"RK"` or `"vdW"`. |
| `kij` | `true` by default: binary interaction parameters from ChemSep's bank. |
| `alpha` | `"equation"` by default, or `"mathias_copeman"`, for PR only. |
| `kij_overrides` | A list of pairs set by hand: `[{ a = "methane", b = "ethane", k = 0.01 }]`. |

Or in format 4 water alone by IAPWS-95 (M5a's D12), `equation =
"IAPWS-95"`, with `[components] names` water by any of its names or its CAS
number, and no other key: it has neither α nor kᵢⱼ, so `alpha`, `kij` and
`kij_overrides` are refused, and so are water with another compound (which
needs a mixture model reflux hasn't yet) and any other compound. Its H and S
are on reflux's reference, as every fluid's are, so its units balance with
the rest; `reflux water` gives its states on IAPWS's
([design.md](design.md#reflux-water)).
[testdata/flowsheets/rankine.toml](../testdata/flowsheets/rankine.toml) is a
Rankine cycle on it.

Or in format 5 water alone by IAPWS-IF97 (M5b's D15), `equation = "IF97"`,
refused as IAPWS-95 is: IF97's regions and region 4's saturation line behind
the flash, its range 273.15 K to 1073.15 K up to 100 MPa and to 2273.15 K up
to 50 MPa, its H and S IF97's plus reflux's reference
([design.md](design.md#if97-as-a-phase-model)). TV and UV flashes refuse
it. [testdata/flowsheets/rankine_if97.toml](../testdata/flowsheets/rankine_if97.toml)
is the same Rankine cycle on it.

Or in format 4 a γ–φ fluid (M6's D12), whose liquid's fugacity is
xᵢγᵢPᵢˢᵃᵗφᵢˢᵃᵗ times the Poynting correction, and whose vapour is the ideal
gas or an equation of state ([design.md](design.md#the-γφ-liquid)):

```toml
[model]
liquid = "NRTL"
vapour = "PR"
henry = ["carbon dioxide"]

[[model.pairs]]
a = "water"
b = "methanol"
tau_ab = "3317.08 J/mol / (R*T)"
tau_ba = "-790.97 J/mol / (R*T)"
alpha = "0.3"
```

| Key | |
|---|---|
| `liquid` | Required: `"Wilson"`, `"NRTL"`, `"UNIQUAC"`, `"UNIFAC"`, `"UNIFAC-LLE"` or `"UNIFAC-Dortmund"`. The binary models take ChemSep's banks, each pair's first listing; UNIFAC each compound's groups from its component file. |
| `vapour` | `"ideal gas"` by default, or an equation, `"PR"` and so on. |
| `kij` | For a cubic vapour: `true` by default, ChemSep's bank for its equation. |
| `poynting` | `true` by default: the Poynting correction, each solvent's liquid volume from its density correlation. |
| `henry` | Henry's-law solutes, by name, beyond those taken as solutes anyway: a gas with Henry's constants in a solvent present whose Tc is below every solvent's vapour-pressure range, such as hydrogen in water. In water, a gas IAPWS covers takes IAPWS's constant (G7-04), and the rest ChemSep's `henrys_law.ipd`. |
| `missing` | `"refuse"` by default: a pair of solvents the bank doesn't list is refused, by name. `"ideal"` gives it no energy of interaction, an ideal pair for NRTL and Wilson and no residual interaction for UNIQUAC, and the report warns of it. |
| `pairs` | `[[model.pairs]]`: a binary model's pairs, each in place of the bank's (below). |

Each `[[model.pairs]]` is a pair of solvents, `a` and `b`, and the model's
parameters for it, each an expression in T that comes out a pure number,
with units, the SI's constants (`R`) and `T` ([Reactions](#reactions-format-2)
has the language):

| Model | Keys |
|---|---|
| NRTL | `tau_ab`, `tau_ba`, `alpha`, all required |
| UNIQUAC | `tau_ab`, `tau_ba` |
| Wilson | `lambda_ab`, `lambda_ba` |

`tau_ab` is a's with b: NRTL's τ_ab in G_ab = exp(−α τ_ab), UNIQUAC's τ_ab,
Wilson's Λ_ab. UNIFAC takes none, its parameters being its groups'. A bank's
pair is written in the same language in the report: NRTL's methanol–water
is `tau_ab = "792.8020 K / (1.98721*T)"` for water with methanol, A₁₂ in
cal/mol over DECHEMA's R, which the bank was fitted with.

## `[transport]`

Optional, in format 5 (M5b's D15): what a phase's transport in the report
takes beyond the compounds' correlations ([report-format.md](report-format.md),
`model.transport`). Nothing in a solve uses it.

```toml
[[transport.pairs]]
a = "water"
b = "ethylene glycol"
G = "2.1 - 300 K / T"
```

| Key | |
|---|---|
| `pairs` | `[[transport.pairs]]`: a liquid pair's Grunberg–Nissan G, ln μ = Σxᵢ ln μᵢ + ½ΣΣxᵢxⱼGᵢⱼ ([design.md](design.md#transport-of-mixtures)). Every pair not listed has G = 0, the logarithmic rule, which leaves an aqueous mixture's viscosity 40% to 70% low ([m5b-plan.md](m5b-plan.md), D6). |

Each `[[transport.pairs]]` has `a` and `b`, two of the components, and `G`,
an expression in T that comes out a pure number, in the language of
`[[model.pairs]]`. A pair is listed once, either way round. Water by
IAPWS-95 or IF97 takes IAPWS's own transport, with no pairs.

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
| `vapour`, `liquid` | The two outlets. With one phase, the other is empty. Two liquids both leave by `liquid`, as from a vessel without a boot, and the unit's report warns of it. |
| `T` or `duty` | One of these. A duty of `"0 W"` is adiabatic. |
| `P` | The drum's pressure, at most the lowest inlet's. By default, the lowest inlet's. |

### `decanter`

A three-phase separator (format 4, M6's D11): a flash drum with an outlet
for each phase.

| Key | |
|---|---|
| `in` | One inlet, or a list. |
| `vapour`, `light`, `heavy` | The three outlets: the vapour, the lighter liquid and the denser, by mass density. Any of them may be empty. A single liquid leaves by `light`. |
| `T` or `duty` | One of these. A duty of `"0 W"` is adiabatic. |
| `P` | Its pressure, at most the lowest inlet's. By default, the lowest inlet's. |

### `component_separator`

| Key | |
|---|---|
| `in` | One inlet. |
| `out` | Two outlets. |
| `fractions` | A table by component: how much of each goes to the first outlet. A component left out all goes to the second. |
| `T`, `P` | Both outlets' state. |

### Reactors (format 2)

Five unit types run a file's reactions (M9, D10), each with one inlet `in`
and one outlet `out`, its outlet a flash at its state, and its balance
checked on its elements, which is what a reaction conserves. A reactor
works in one phase, by an equation of state: on a γ–φ fluid it's refused,
saying so, and an outlet that would split is refused too (M6's D11).
[testdata/flowsheets/reformer.toml](../testdata/flowsheets/reformer.toml)
is a whole file: a steam reformer as a packed bed, on Xu and Froment's
kinetics.

```toml
[unit.reformer]
type = "pfr"
in = "feed"
out = "product"
volume = "1 m^3"
length = "1 m"
catalyst = { loading = "0.386 kg/m^3", void_fraction = 0.4, particle_diameter = "2 mm" }
pressure_drop = "ergun"
```

Every reactor takes these, each optional:

| Key | |
|---|---|
| `T`, `duty`, `adiabatic` | At most one. `T` keeps the outlet (a PFR, all of it) at that temperature; `duty` adds that heat; `adiabatic = true` adds none. With none of them the reactor is isothermal at its inlet's T. An isothermal reactor reports its duty. |
| `pressure_drop` | The outlet's P below the inlet's. A PFR's may be `"ergun"` instead (below). |
| `reactions` | The file's reactions it runs, by name. By default every one it can: every one with a rate law for a PFR or a CSTR, every one for an equilibrium reactor. Not for a Gibbs reactor or a conversion reactor. |

A reactor's outlet is one phase: one that would split is refused, a second
phase being M6's.

#### `gibbs_reactor`

Its outlet is the equilibrium over every component the feed's elements
allow (D9), whatever reactions the file has, at its T or where its energy
balance puts it.

| Key | |
|---|---|
| `inert` | Components that pass through, by name. A component without a formula does anyway. |

#### `equilibrium_reactor`

Its reactions reach equilibrium, each one's Q its K from the Gibbs energies
(D7), and nothing else reacts. The reactions must be independent.

#### `conversion_reactor`

| Key | |
|---|---|
| `conversions` | Its steps, in order: each a table with a `reaction`, a `key` reactant by name, and its `conversion` of the key as it enters that step. A later step sees what earlier ones left. A step that needs more of a reactant than there is is refused. |

A `conversion` is a number from 0 to 1, a fraction with a unit such as
`"27.4 %"`, or an expression in `T`, the reactor's temperature, with the
file's definitions: a correlation typed as its source prints it, the
temperature read on the scale it gives.

```toml
conversions = [{ reaction = "synthesis", key = "carbon monoxide", conversion = "0.85 - 0.0018*(T in degC)" }]
# or, as a percentage
conversions = [{ reaction = "synthesis", key = "carbon monoxide", conversion = "(85 - 0.18*(T in degC)) [%]" }]
```

The expression is taken at the outlet's T: the given or the inlet's for an
isothermal reactor, and the one the energy balance finds for an adiabatic
one or one at a duty. It must give 0 to 1 there, or the reactor fails,
naming the T. It may use T alone. Typed without a reading, `0.85 -
0.0018*T` is refused, since T is a temperature, with the reading suggested.

A conversion reactor warns, in its report, of each reaction it takes past
equilibrium, Q at its outlet above K at 1 bar from the Gibbs energies,
with Q/K: no reactor could convert that much at that T and P.

#### `cstr`

Its outlet is its contents: each component's flow in less its flow out,
plus V times what its reactions make, is zero, solved by Newton.

| Key | |
|---|---|
| `volume` | The reactor's volume, `"0.2 m^3"`. |
| `catalyst_loading` | kg of catalyst per m³ of reactor, for rates per mass of catalyst. |

#### `pfr`

A plug-flow reactor or packed bed: the flows integrated along its volume
(D8), isothermal, adiabatic or with its duty spread evenly over its volume.

| Key | |
|---|---|
| `volume` | The reactor's volume. |
| `length` | Its length, which Ergun's drop needs, and which the profile then gives. |
| `catalyst` | `loading`, kg of catalyst per m³ of reactor, for rates per mass of catalyst; for Ergun's drop, `void_fraction` and `particle_diameter`. The void fraction enters Ergun alone, as in DWSIM. |
| `pressure_drop` | A drop, spread evenly over the volume, or `"ergun"`: Ergun's equation, on the gas's viscosity at low pressure, each component's from its file's correlation mixed by Wilke's rule. |
| `integrator` | `"rodas4"` (the default), stiff, or `"dopri5"`, which refuses a problem it finds stiff. |
| `rtol` | The integrator's relative tolerance, 1e-8 by default. |

The report gives a PFR's profile: its state at its inlet and at every point
the integrator stepped to (report format 2).

## Reactions (format 2)

A file with reactions is format 2: `[definitions]`, `[site.NAME]` and
`[reaction.NAME]`, beside `[components]`. `reflux kinetics check` and
`reflux kinetics explain` read them (M9 step 6), and so does `reflux solve`,
for its reactors. A file of reactions alone has `format`,
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

## Specs (format 3)

A spec varies a number the file gives until an equation over the
flowsheet holds (M8's D5): DWSIM's Adjust, met with the flowsheet's loops
rather than around them.

```toml
format = 3

[spec.vapour]
vary = "unit.cooler.outlet_T"
between = ["200 K", "260 K"]
until = "gas.F = 60 mol/s"
```

| Key | |
|---|---|
| `vary` | Required. A key the file gives a feed or a unit, as below. |
| `between` | Required. The least and the most it may take, each in a unit of its own (or a plain number for a fraction). The file's value must be between them: the spec starts from it. |
| `until` | Required. An equation, two expressions with `=` between them, of one dimension, over the flowsheet's variables, below. |

`reflux solve` meets every spec, and the report gives each one's value, its
equation's two sides, and how far apart they are
([report-format.md](report-format.md)). A target its bounds can't reach is
refused with what they do reach:

```
reflux: plant.toml:65: spec.vapour: its target is out of reach: with unit.cooler.outlet_T from 200 K to 260 K, gas.F goes from 41.9919 mol/s to 76.6534 mol/s, never 200 mol/s (at 10 values)
```

**A key to vary** names a feed's or unit's table and a number in it, as the
file writes them: `feed.syngas.T`, `unit.reactor.volume`,
`feed.syngas.molar_flows.methane`, `unit.split.fractions[0]`,
`unit.tube.catalyst.loading`. It must be a number the file gives and its
maker reads as one: `unit.cooler.duty` is refused where the cooler gives
`outlet_T`, and so are `type`, `in`, and a conversion that's an expression
in T. Each new value is checked as the file's own would be: a split whose
fractions don't sum to 1 is refused, so two fractions are varied together.

**The flowsheet's variables,** in a spec's equation and later a study's
outputs, a name with one dot:

| Name | |
|---|---|
| `gas.T`, `gas.P` | a stream's temperature and pressure |
| `gas.F` | its molar flow |
| `gas.n[methane]` | a component's molar flow |
| `gas.y[CH4]` | a component's mole fraction, of the whole stream |
| `gas.h` | its molar enthalpy |
| `gas.vf` | its vapour fraction |
| `cooler.duty`, `pump.work` | a unit's duty or work: an adiabatic reactor has none |

A stream whose name isn't one word is quoted: `'let down'.T`. A component
goes by the file's name, its formula, its CAS number or its aliases, quoted
where it isn't a word, `gas.n['n-butane']`. Every number has its unit, as in
any expression ([design.md](design.md#expressions)), so `gas.T = 300` is
refused, and `gas.T = 300 K` reads.

## Studies (format 3)

A study solves the flowsheet at each value of one input, or each pair of
two, and reports expressions over the flowsheet at each point (M8's D16).
`reflux study plant.toml [NAME]` runs it ([report-format.md](report-format.md)).

```toml
[study.temperature]
vary = [{ key = "feed.syngas.T", from = "200 degC", to = "320 degC", points = 13 }]
report = { CO_conversion = "(syngas.n[CO] - product.n[CO]) / syngas.n[CO]", T_out = "(product.T in degC)" }
```

| Key | |
|---|---|
| `vary` | Required. A list of one input or two, each `{ key, from, to, points }`, evenly spaced from `from` to `to`, or `{ key, values }`, a list. A key is as a spec's (above), and a key a spec varies can't be. Each value is read as the key's maker reads it, in a unit of its own. |
| `report` | Required. Named expressions over the flowsheet's variables, as a spec's `until`'s sides. |
| `chunk` | How many points run from one cold start, 8 by default: each chunk's first point is solved cold, and the rest from their neighbours. |
| `sensitivities` | `{ outputs = [...], inputs = [...] }`: at each point, the derivative of each of these outputs, by name, with respect to each of these keys, by central differences of solves at 1e-12, with its error's estimate ([design.md](design.md#sensitivities)). Each input costs five more solves a point. A key a spec varies can't be one. |

Two inputs make a grid: the first's values outer, the second's inner. An
expression that has no value at a point, such as a product's conversion
when the feed has none, is reported with its reason, and so is a point
that fails; the study goes on. [testdata/flowsheets/methanol_study.toml](../testdata/flowsheets/methanol_study.toml)
is Tutorial 2's.

## Optimisations (format 3)

An optimisation finds the inputs, within their bounds, that minimise or
maximise an expression over the flowsheet, subject to equations and
inequalities over it (M8's D13). Every evaluation is a converged solve,
specs and all; the optimum is returned only once its first-order
conditions hold ([design.md](design.md#optimising-a-flowsheet)).
`reflux optimize plant.toml [NAME]` runs it
([report-format.md](report-format.md)).

```toml
[optimize.best]
vary = [
    { key = "unit.tube.volume", between = ["0.1 m^3", "2 m^3"] },
    { key = "feed.f.T", between = ["380 K", "420 K"], start = "400 K" },
]
maximize = "product.n[isopentane]"
subject_to = ["product.n[neopentane] <= 0.2 mol/s"]
```

| Key | |
|---|---|
| `vary` | Required. A list of inputs, each `{ key, between }`, its least and most values, and optionally `start`, else the file's value, which must be within them. A key is as a spec's (above), and a key a spec varies can't be. |
| `minimize` or `maximize` | Required, one of the two. An expression over the flowsheet's variables, as a spec's `until`'s sides. |
| `subject_to` | A list of equations and inequalities over the flowsheet's variables, each with `=`, `<=` or `>=` (`≤` and `≥` too) between two sides of one dimension. |

An optimisation that finds no optimum, because no point meets its
constraints or its iterations run out, says why, and its report gives each
iterate. [testdata/flowsheets/series.toml](../testdata/flowsheets/series.toml)
has series reactions in a PFR, whose optimum sits on a bound and a
constraint at once.

## Fits (format 3)

A fit finds a file's constants that best match measured data, with each
constant's standard error and 95% interval (M8's D14, D15). It's in a
format-3 file with the reactions whose constants it fits, the `[model]`
whose kᵢⱼ it fits, or both, or in format 4 or 5 with the activity model
whose pairs it fits (M6's D16); a flowsheet file's fits read too.
`reflux fit file.toml [NAME]` runs it ([report-format.md](report-format.md)).

```toml
[fit.rates]
data = "ethylene.csv"
parameters = ["k1.A", "k1.Ea", "K_E.A", "K_E.dH", "K_H"]
weights = { "r[hydrogenation]" = "5 %" }

[fit.kij]
data = "ethane_propane.csv"
parameters = ["kij[ethane, propane]"]
weights = { P = "0.5 %", "y[ethane]" = "0.003" }
```

| Key | |
|---|---|
| `data` | Required. A CSV file, beside this one unless its path is absolute (below). |
| `parameters` | Required. What the fit varies, all of one kind. A reaction's named constant (`k1`, the rate constant of the first reaction, or the law's own `k_name`; a site term's constant by its `name`; a fitted K, `K1` or its own name) by its part: `NAME.A` and `NAME.Ea` of an Arrhenius constant, `NAME.A` and `NAME.dH` of a van 't Hoff one, `NAME` of a number with its unit. Or a pair's kᵢⱼ, `kij[a, b]`, which needs a `[model]` with an equation. Or an activity model's pair: `A[a, b]`, a's energy with b in J/mol, NRTL's τ_ab = A/RT or UNIQUAC's τ_ab = exp(−A/RT), and `alpha[a, b]`, NRTL's α, which need a `[model]` with `liquid = "NRTL"` or `"UNIQUAC"`. |
| `start` | Where to start a parameter, by name, a number with its unit; the file's value if not given. A pair's energy starts at the model's τ at 298.15 K, as A = τRT (NRTL) or −RT ln τ (UNIQUAC): a bank's A₁₂ in cal/mol over DECHEMA's R is that at any T. |
| `weights` | Each measured column's error: `"5 %"` of each measurement, or a number with the column's unit. A rate fit's columns are unweighted (an error of 1 in SI) where none is given; a kᵢⱼ fit needs P's and each y's, since they're of different kinds. |

**Rate data** have a row for each state: `T [K]` and `P [bar]` (any unit of
each), `y[species]`, each species' mole fraction (0 if its column is left
out), and `r[reaction] [mol/(kg*s)]`, each measured rate in its rate's
dimension, empty where it wasn't measured. Each row's rates are on the ideal
gas, p = yP, as from a differential reactor. An Arrhenius or van 't Hoff
constant whose A and energy are both fitted is fitted centred on the data's
mean T, and reported back in the file's own form.

**VLE data** have a row for each bubble point: `T`, `P`, `x[component]`
for every component of the liquid but one, which makes up the rest, and
`y[component]` for each measured vapour fraction. A fit of an activity
model's pairs is written back as the pairs' `[[model.pairs]]` tables, each
fitted energy as `"A J/mol / (R*T)"` (UNIQUAC's in an `exp`), the pair's
other parameters as they were. Data are worth giving the point test for
Gibbs–Duhem consistency first (`fit.point_test`, M6 step 10).

A fit that can't be made says why: data it can't read, a parameter the data
don't determine (an Arrhenius constant's A and Eₐ from data at one T), or a
solver that couldn't meet its test. [testdata/m8_fit/](../testdata/m8_fit/)
has one of each kind.

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
