# Spike: driving DWSIM 9 from Python (M3 step 0)

**Verdict: go**, with the conditions listed under "What step 5 must do".
DWSIM 9.0.5 runs from a Python script with no GUI. It builds a flowsheet,
flashes PT and PH both through the flowsheet solver and through the
property package directly, and hands back the constants, Cp coefficients
and kᵢⱼ it used. The manual GUI fallback in D9 isn't needed.

The spike's main result is that the whole gap between DWSIM and an exact
Peng-Robinson is explained, by three things:
- **DWSIM's PR uses rounded constants**: Ωa = 0.45724, Ωb = 0.0778 and
  R = 8.314. With those, a short PR in the script reproduces its enthalpy
  departures to 5e-16.
- **Its PT flash stops at a fugacity residual of 8e-7**, however tight
  the tolerance is set.
- **Its ideal-gas enthalpy integrates Cp approximately**, to about 1e-5.

On the way, the spike also found that `thermo`'s default flash tolerance is
too loose for criterion 2.

The spike is one script, [spike/dwsim/run.py](../spike/dwsim/run.py), a uv
script like the fixture generators:

```bash
uv run --python 3.13 spike/dwsim/run.py
```

It needs Windows and DWSIM 9 installed for the current user
(`%LOCALAPPDATA%\DWSIM`). It prints every number below.

## Results (Windows 11, DWSIM 9.0.5.0, Python 3.13.17, pythonnet 3.0.5, `thermo` 0.6.1)

The case is 40 % methane and 60 % n-butane at 250 K and 30 bar, on
DWSIM's "Peng-Robinson (PR)" package. The timing ranges cover four runs.
They vary a lot from run to run (the PH flash by 2×), so they're a first
look; step 10 measures properly.

| Check | Result |
|---|---|
| Load DWSIM, create a flowsheet, no GUI | yes: 2–3 s loading assemblies, 4–6.5 s `Automation3()`, 0.5–0.8 s `CreateFlowsheet()` |
| PT flash, flowsheet solver | β = 0.23926872, y(CH₄) = 0.97615638, x(CH₄) = 0.21879714, H = −18274.659 J/mol |
| PT flash, property package directly | the same numbers, bit for bit, in 6 iterations |
| PH flash at that H | back to 250 K within 3e-9 K (1.7e-13 K with tight tolerances) |
| Constants and kᵢⱼ read back | Tc, Pc, ω, M, ChemSep eq 16 Cp coefficients, kᵢⱼ = 0.0244 |
| DWSIM's fugacity residual, \|ln fᵥ − ln fₗ\| | 9.6e-5 at the default tolerances (1e-4); 8.0e-7 at 1e-6, 1e-8 or 1e-12 |
| DWSIM (tolerances 1e-12) minus `thermo`'s PR, same constants and kᵢⱼ | β 8.9e-7, y 6.1e-6, x 1.1e-6 |
| DWSIM minus a PR with DWSIM's rounded constants | β 4.1e-7, y 4.8e-8 |
| Enthalpy departures, DWSIM minus PR with its rounded constants | 5e-16 relative or less (against exact constants: 5.6e-5 in the vapour, 1.0e-4 in the liquid) |
| Ideal-gas enthalpy, DWSIM minus an exact quadrature of the same Cp | up to 1e-5 relative: 0.31 J/mol for n-butane from 298.15 K to 600 K |
| Timings, warm, through pythonnet | PT flash 1.0–1.6 ms, PH flash 6.7–15 ms, one-stream flowsheet solve 2.3–4 ms |

## What we learned

1. **Loading.**
   - pythonnet reaches DWSIM's classic build (.NET Framework 4.6.2) as
     `pythonnet.load("netfx")`. pythonnet 3.0.5 supports Python up to 3.13,
     not 3.14.
   - DWSIM finds its databases relative to the working directory, so the
     script changes into the install folder before `clr.AddReference`.
   - No window opens: `Automation3` runs headless. (DWSIM's FluentAPI
     needs version 10, so it's not an option here.)
   - One property package, ThermoC, fails to load because its `ThermoCS`
     assembly isn't installed. DWSIM prints a stack trace and carries on.
2. **Compounds are keyed by display name, and some names appear twice.**
   "N-butane" is ChemSep's (ω = 0.199) and "n-Butane" is CoolProp's
   (ω = 0.2008). Capture scripts must pick the ChemSep entry by name and
   check `OriginalDB`.
3. **pythonnet hands back interfaces.**
   - `CreateAndAddPropertyPackage` returns an `IPropertyPackage` with 30
     members. `.__implementation__` gives the concrete
     `PengRobinsonPropertyPackage`, with `RET_VKij`, `FlashSettings`,
     `DW_CalcFugCoeff` and the rest.
   - `AddObject` returns an `ISimulationObject`, and `.GetAsObject()` gives
     the `MaterialStream`.
4. **The property package needs a current stream, and every flowsheet
   solve clears it.**
   - `CurrentMaterialStream` must be set before any direct call:
     `CalculateEquilibrium`, `RET_VKij` or `RET_Hid`.
   - It must be set again after each `CalculateFlowsheet2`.
   - A missing stream shows up as a bare `NullReferenceException`.
5. **The direct flash is the flowsheet's flash.** It gives the same numbers
   bit for bit, without the stream-property work, so it's the faster way to
   capture fixtures. It reports H in kJ/kg; multiplied by the molar mass in
   g/mol, that's J/mol. The enthalpy reference is the ideal gas at
   298.15 K, without heats of formation.
6. **Tolerances, and a floor.**
   - DWSIM's defaults are 1e-4 for the PT and PH loops. They leave a
     fugacity residual of 9.6e-5, and β = 3.8e-5 away from DWSIM's own
     tight answer. That error alone is more than a third of criterion 6's 1e-4.
   - Tightening to 1e-6 brings the residual to 8.0e-7, and nothing tighter
     changes anything. The nested-loops flash stops after 10 iterations
     with the same answer at 1e-6, 1e-8 and 1e-12, with or without
     `NL_FastMode`.
   - Its compositions are within 5e-8 of the converged solution of its own
     equations, but β is 4e-7 away from it.
7. **DWSIM's Peng-Robinson uses rounded constants.**
   - The constants are Ωa = 0.45724, Ωb = 0.0778 and R = 8.314. The exact
     values are 0.45723553, 0.07779607 and 8.314462618.
   - With DWSIM's values, the script's own PR reproduces DWSIM's enthalpy
     departures to 5e-16 in both phases.
   - With exact values, which is what `thermo` and reflux use, the
     departures differ by 1e-4 relative: 2 J/mol in this liquid.
   - R cancels out of a PT flash, since A and B are each a ratio of R
     terms, but it scales every departure function.
8. **DWSIM's ideal-gas enthalpy is an approximate integral** of the same
   ChemSep eq 16. It's 1e-5 relative off an exact quadrature (0.3 J/mol at
   600 K). reflux's `Correlated` integrates to round-off, so that part of
   the difference belongs to DWSIM.
9. **DWSIM's copy of ChemSep isn't reflux's.**
   - n-butane's Cp coefficients carry one more digit in DWSIM: A is
     44749.95 against 44749.9 in reflux's copy (ChemSep LITE 8.50, whose
     `chemsep1.xml` says v8.49). Methane's B and E differ in the same
     way.
   - Tc, Pc, ω and M agree for both compounds.
   - The kᵢⱼ, 0.0244, is the first entry pr.ipd lists for the pair, which
     is the choice reflux made in M2 (D12).
10. **`thermo`'s default PT flash is too loose for criterion 2.**
    - FlashVL's `PT_SS_TOL = 1e-13` is a bound on a sum of squares. It
      leaves a fugacity residual of 1.1e-7 and β 2.7e-8 off: more than
      criterion 2's 1e-8 by itself.
    - With `PT_SS_TOL = 1e-30`, the residual is 1.8e-15 and β agrees with
      the script's independent PR to 1.4e-16.
    - This belongs to the `thermo` fixtures in step 5, not to DWSIM, but
      the spike is where it showed up.

## What step 5 must do

- **Capture DWSIM with its PT and PH tolerances at 1e-8.** Nothing past
  that changes the answer.
  - Record DWSIM's own fugacity residual with each case, so that its
    convergence floor shows in the fixture, not in reflux's error.
  - Use the property package's `CalculateEquilibrium`. Set
    `CurrentMaterialStream` first.
- **Use what DWSIM reports.** That means its Tc, Pc, ω, M, Cp coefficients
  and kᵢⱼ, read back by the generator and stored in the fixture, never
  reflux's database values.
- **Expect DWSIM's rounded constants to show, and record them as the
  explanation.**
  - They shift β by about 5e-7 here, which leaves criterion 6's 1e-4 for
    β and compositions plenty of room.
  - They shift enthalpy departures by about 1e-4 relative, which matters
    for PH. Give reflux the H that DWSIM computed, and a subcooled liquid
    lands 2 J/mol ÷ about 130 J/(mol K) ≈ 0.015 K away, over criterion 6's
    0.01 K.
  - So compare PH and PS cases by state. Give both programs the same
    (T, P) and z, let each compute its own H or S, flash back, and compare
    the T, β and compositions each one lands on.
  - Compare H and S at fixed (T, P) separately, against a recorded relative
    tolerance that the rounded constants explain.
- **Tighten `thermo` too.** Set `PT_SS_TOL` to 1e-30, with matching
  tolerances for the other flash types, and store each fixture's fugacity
  residual next to its values.
- **Run by hand on Windows**, like the CoolProp fixtures. CI never runs
  DWSIM, and the fixture records DWSIM's version.

## For the benchmarks (D10)

- DWSIM's timings come from the same automation route: warm, through
  pythonnet, from the property package. The per-call pythonnet overhead is
  microseconds, small against 1–10 ms flashes.
- Start-up costs 7–10 s per process. The first flowsheet solve takes
  0.5–0.8 s and the first direct flash about 30 ms, while .NET compiles
  them. A benchmark has to warm up first and keep
  one process running.

## Still open

- **DWSIM's other flash approaches.** The spike used the default:
  `UniversalFlash`, nested loops. Inside-out (`UseIOFlash`) or Gibbs
  minimisation might converge past the 8e-7 floor. That's worth one try in
  step 5 if the floor gets in the way. It doesn't block anything: the floor
  is well inside criterion 6.
- **Harder points.** The spike flashed one easy two-phase point. Near-
  critical points, single-phase points and pure components are step 5's
  20 cases. Since DWSIM's convergence there is unknown, each case records
  its residual.
