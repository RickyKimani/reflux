# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What reflux is, and for whom

An open-source steady-state chemical process simulator in V. The aim is to
be faster and nicer to use than DWSIM. The README's "Why another simulator"
table lists the DWSIM problems that drive the design: the UI, stale docs,
DWSIM 10 not working, slow flowsheets and sensitivity runs, no NumPy in
scripting, and LHHW kinetics typed as text.

It's also how the owner learns V. He is a chemical-engineering student,
strong in thermodynamics, who knows Go. Explain V by contrast with Go, and
put lasting lessons in `docs/v-for-go.md`, checked by compiling a probe
before they're written. He works on Windows with nushell. His Go package
`C:\Users\Ricky\dev\zfactor` is a reference only: don't copy from it.

## Where things stand

Keep this section current: update it whenever a step is committed.

- **M5b (transport of mixtures, IF97 behind the flash) is done,** not
  committed past the plan, its exit review at the end of its plan, on
  branch `claude/dazzling-lehmann-3ad069` from `main` at v0.6.0
  (`a1d1b91`), in the worktree `.claude/worktrees/dazzling-lehmann-3ad069`.
  Its plan, `docs/m5b-plan.md`, was approved on 2026-11-27 as drafted:
  criteria 1–12, decisions D1–D22, steps 0–13 in four bundles, A (1–4,
  transport), B (5–6, IF97 behind the flash), C (7–8, file format 5,
  report format 5, Python), D (9–13, data, DWSIM, speed, release, review).
  The plan and step 0's spikes are committed at `5c49cfd`.
  - Step 0 (`spike/transport/`): DWSIM 9.0.5's mixture transport
    reproduced from its own code to 2e-16 (a mole-mean liquid viscosity,
    Jossi–Stiel–Thodos at every density, its second liquid given the
    first's μ and λ); the candidates against ThermoML's measured mixtures,
    which chose each method (Lucas's factor 2.3–2.8% on dense hydrocarbon
    gases against JST's 8–12%); ChemSep's water within 2% of IAPWS's; IF97's
    regions 1 and 2 up to 5.4e-5 RT apart on region 4's line, hence D13.
  - **Bundle A (steps 1–4) is done, not committed:** the `transport`
    module (each compound's values, Wilke's rule moved from `reactors` with
    its bits, Mason–Saxena, Lucas's factor, Stiel–Thodos, Grunberg–Nissan,
    Li with Missenard, Winterfeld–Scriven–Davis, a phase's transport by
    D9), `components/transport.v`, `water.Transport95`, data-quality.md's
    transport section, `tools/gen_m5b_transport.py` and
    `testdata/m5b_transport.json`. Every transport correlation `chemicals`'
    within 2.8e-16 of its terms; each method within 1.3e-15 of `chemicals`
    or mpmath; 272 seeded phases with values or reasons. Planted mistakes:
    41, 38 caught, 1 after a test was added, 2 explained. The full tcc
    check passed (142 test files).
  - **Bundle B (steps 5–6, IF97 behind the flash) is done, not
    committed:** `water.If97Phases` (a `thermo.PhaseModel`: the region map,
    each side's equation across the line, residuals against region 2's
    ideal gas as a `thermo.Polynomial` Cp, `TransportIf97`),
    `thermo.Saturates` and `thermo.NoState`, the flash's `own_line` for a
    pure model's own saturation line, its step back for a phase model,
    `Database.if97`. The phases are IF97's states (h 1.2e-10 J/mol); region
    4 to the bit; PH and PS within 3e-11·T of IF97's own but on region
    boundaries. D14 and criterion 7 amended, approved on 2026-11-28: no new Cp form;
    on B13 and at 1073.15 K one h has two T (up to 58.7 mK). Found and
    fixed in IF97 (M5a's): a closed bracket's endless halving, region 5
    reached above 50 MPa, a search's probe past the range's edge. M5a's
    trace, with IAPWS-95 water added, 86,408 lines the same bytes before
    and after; planted mistakes 31, all caught, 2 after a test was added.
    The full tcc check passed (143 test files).
  - **Bundle C (steps 7–8, the file, the reports, Python) is done, not
    committed:** file format 5 (`equation = "IF97"`, `[[transport.pairs]]`
    with Grunberg–Nissan G as expressions), report format 5 (each phase's
    μ, λ and σ or null with `transport_notes`, `model.transport`),
    `Loaded.transport`, `testdata/flowsheets/rankine_if97.toml` and its
    golden; `reflux_transport` (a status per value, `REFLUX_POINT_NO_VALUE`,
    ABI still 1), `Fluid.transport`. The nine goldens, transport taken out,
    are their format-4 selves value for value; the IF97 cycle is IF97's by
    hand (4.9e-7 J/kg at worst) and DWSIM's Steam Tables' by state. Found:
    json2 knows a type by its short name, and a smart-cast spread into an
    interface variable copies garbage (both in v-for-go). Planted mistakes
    39: 38 caught, 1 after a test was added, 1 explained.
    `tools/python.vsh test` passed (94). The full tcc check ran 145 test
    files: 144 passed, and `reactors/flowsheet_test.v` (the report still at
    format 4) passed alone once fixed; fmt, vet, build and smoke by hand.
  - **Bundle D (steps 9–13) is done, not committed:**
    `tools/gen_m5b_published.py`, `testdata/m5b_published.json` (27 ThermoML
    sets, 1,662 points, two new nonaqueous liquids; M5b's deviations in
    `transport/published_test.v`, 89 numbers that can't grow, each
    alternative's beside them); `tools/gen_m5b_dwsim.py`,
    `testdata/m5b_dwsim.json` (DWSIM's PR and NRTL packages at ten states,
    every number its formulas' to 4.7e-16, its data the files', every
    difference the methods', 35 held signed in `transport/dwsim_test.v`);
    `tools/bench_m5b.v`, `tools/gen_m5b_bench.py`, `testdata/m5b_bench.json`,
    `benchmarks/m5b.v` (a phase's transport 3.4–22× DWSIM's own functions);
    v.mod and the package at 0.7.0, its CHANGELOG entry undated, the ten
    goldens written again (only the version moved); docs and the exit
    review. Found: DWSIM's σ order-dependent with a compound above its Tc,
    that compound's μ taken with its Henry's constant in water, its
    CoolProp package's water CoolProp's (not ChemSep's, as step 0 read),
    and the mass-fraction basis beating D6's on four liquid sets (for the
    owner). Planted mistakes 10, all caught (one only once step 10's
    differences were held two-sided). The full tcc check passed (148
    test files).
  - Found, for the owner: 622 usable liquid viscosities (not 637); 49
    compounds with no dipole moment, so no gas viscosity with one;
    Cantera's 0.6990 against polykin's and reflux's 6.990 in Lucas's Z₂
    below Tr = 1, for the owner's Poling; criterion 3's worked examples
    need that book too.
- **M0 to M4, M10a, M9 and M8 are done** and merged into `main` (PRs #1 to
  #7).
  `docs/roadmap.md` has every milestone (M0–M14) and its exit criteria.
  Each milestone's plan, `docs/mN-plan.md`, ends with its exit review, and
  each finished step there has its write-up: the record of what it did.
- **M5a (water and steam) is done,** merged in PR #10 (`4bf65d5`) and
  released as v0.6.0, with its exit review at the end of its plan; built on
  branch `m5a-water`, from
  `main` at v0.4.0. `main` (M6, v0.5.0, and the gas constant's fix that
  begins 0.6.0) was merged into it before bundle C, at the owner's word on
  2026-11-26, to meet M6's seam and format 4 before steps 7 and 8. Its plan,
  `docs/m5-plan.md`, was approved on 2026-11-24 as drafted: exit criteria
  1–12, decisions D1–D18 and steps 0–12, ending in the release `main` has
  begun, 0.6.0 (D17). Its "With M6" section has the merge, what broke and
  how it was put right.
  - D1 splits the roadmap's M5: M5a is IAPWS-95 as a `thermo.Model`, IF97
    with its region-3 supplements as steam-table functions, water's
    transport by IAPWS, `equation = "IAPWS-95"` in the file, `reflux water`
    and IAPWS-95 in the C API; M5b is mixtures' transport and IF97 behind
    the flash on M6's seam. CoolProp isn't linked, only a fixture source.
  - **Step 0 is done** (`spike/water/`, Python, porting `thermo`'s solvers
    line for line):
    - IAPWS-95's isotherms inside the dome have spurious roots with Gibbs
      energies up to 1e19 RT too low (910 of 2,508 grid states), so the
      stable root is only ever chosen between the two physical branches'.
      The generic solvers mostly do that, but the vapour start p/RT landed
      in the loops (300 K, 42 MPa gave 322 kg/m³), the liquid walk returned
      an artificial vmin as a root, and near Tc each solver could return the
      other branch's root. D7: an optional interface (IAPWS-95's alone)
      splitting the branches at ρc below Tc (a vapour-start bound by SR1-86's
      ρ'' too, which step 2 found redundant and took out), a liquid walk that
      must reach p, and a bracketed start for
      `saturation_pressure` (which failed from 646 K, and now converges to
      0.1 mK below Tc). A convexity guard, tried first, refused real liquid
      roots at 251–260 K.
    - IAPWS-95's R (8.314371357587) is carried in F (D4), so its tables
      match exactly; water's ln φ is then about 1e-4, not 0, at low P.
    - Table 5's 0/0 at δ = 1 is removable; only (Tc, ρc) is singular
      (φʳ_ττ → −∞, D8). The `iapws` package is GPL-3.0: fixtures only (D10).
    - DWSIM 9.0.5's CoolProp package is IAPWS-95 to 1e-9 kJ/kg; its Steam
      Tables (IF97) returns h = s = −1 silently in region 3, has no region 5,
      and its region-3 PH flash lands 1080 K off.
  - The steps go in four bundles (the plan, below its step graph), allowed
    on 2026-11-25: A, steps 1–2; B, 3–5; C, 6–8; D, 9–12. A full tcc check
    ends A, C and D; B touches only `water`, which nothing imports yet, so
    it runs `water`'s tests, fmt and vet instead.
  - The version is set at the release, as the next on `main` then (D17):
    0.6.0, whose CHANGELOG entry the gas constant's fix began, and M5a adds
    to it.
  - **Bundle A (steps 1 and 2) is done,** committed at `181e456`, its full
    tcc check passed (108 test files), and amends D7 and D8, approved on
    2026-11-26:
    - Step 1: `thermo.PlanckEinstein`; the `water` module, IAPWS-95 as a
      `Model` (`water.iapws95()`, coefficients transcribed from the release),
      its ideal part, `water.fluid`, Table 3 (`props95`) and D6's constants
      (`reference_offsets`). Tables 6 and 7 to their printed digits, by Table
      3 and through `thermo`; φʳ within 18.3 ε·scale of mpmath at 291 points.
      D8 amended: `Fluid.props`' w is NaN at (Tc, ρc) (Cp/Cv is ∞/∞ in the
      cubics' formula, which stays).
    - Step 2: `thermo.Branches` (`split`, `no_state`), the
      liquid walk that must reach p, the bracketed Psat start, labels by the
      split, and `flash.at_or_closer`, a search's step back from a T where
      the model has no state (never for a cubic). Table 8 through `tvf`; Psat to 0.1 mK below Tc; every flash back
      from 925 states. D7 amended: `min_volume` at δ = 4; no state below
      251.165 K (the liquid's H falls with T there, and a vapour gave UV a
      phantom); Π labelled liquid water vapour near its density maximum.
    - Every cubic answer keeps its bits: an 83,563-line trace of the
      sweep's mixtures and pure compounds through every flash, gcc `-prod`,
      before and after (`spike/water/trace`).
    - Found: the `iapws` package's densities are wrong at 640–647 K below a
      few hundred Pa; a water flash costs about 0.06 s under tcc (step 10).
  - **Bundle B (steps 3 to 5, IF97) is done** and committed, and amends
    criterion 6, approved on 2026-11-27. `water`'s tests, fmt and vet ran,
    not the full check (its exception, above).
    - Step 3: `water.if97`, a submodule (D2 allowed it): the basic equations
      of regions 1, 2 (and the metastable vapour), 3 and 5, region 4, B23
      and B2bc, per kg on IAPWS's reference. Its 1,555 coefficients are
      read from the three releases' PDFs by `tools/extract_if97.py` (their
      plain text loses exponents and the supplements' digits); each equals
      the `iapws` package's (`tools/gen_m5_if97.py`) but B23's n₅, which the
      package has with an 8 dropped.
    - Step 4: the backward equations as published (R7-97's, SR3-03's,
      SR5-05's 26 subregions) and the steam-table functions from (p, T),
      (p, h), (p, s), (T, x), (p, x), each answer finished on its basic
      equation (T back to 3.6e-15, region 3's 4.4e-14); region 3's
      saturated phases to 1 µK below Tc, each bracketed on its side of ρc.
      Every verification table to its digits but SR5-05's Table 13's last
      row, printed 0.53 units off its own equation; Tables 43 and 44
      reproduced, maxima and RMS, but for a slip at 623.15 K in Table 44
      (0.0046%, printed 0.0041%) and four maxima of Table 43 (0.4% of them).
    - Step 5: IF97 against IAPWS-95 at 4,712 states, each difference within
      R7-97's Figs. 3 to 6 but on the pseudo-critical ridge above the
      critical triangle (30–40 MPa), where cp differs by up to 0.99% against
      the figure's 0.5%: criterion 6 amended to report it apart.
    - Found besides: the `iapws` package's region-3 mixtures are SR5-05's
      unfinished densities (1% off at 22 MPa), and a NaN passed the first
      tests' `> bound` checks (now `!(… <= bound)`).
  - **Bundle C (steps 6 to 8) is done and committed,** and amends
    criterion 9, approved on 2026-11-27. `tools/python.vsh test` passed
    (88). The full tcc
    check ran 136 test files: 135 passed, and `capi/capi_test.v`'s header
    test, which didn't know `REFLUX_IAPWS95`, passed alone once it did;
    fmt, vet, the build and the smoke runs then passed by hand. Planted
    mistakes: 36, all caught, 4 after tests were added for them.
    - Step 6: `water/transport.v`, R12-08's μ and R15-11's λ with their
      critical enhancements on IAPWS-95, their industrial forms on IF97
      (IF97's state gained (∂ρ/∂p)_T), R1-76's σ; every verification table
      to its digits; the `iapws` package's 791 and 775 states to 6e-14 and
      2.9e-11 (`tools/gen_m5_transport.py`). D11's ranges as read: from
      250 K, to 1273.15 K, above it a dilute gas (ρ ≤ 0.01ρc).
    - Step 7: format 4's `equation = "IAPWS-95"` (`Database.iapws95`), its
      refusals at their lines; `testdata/flowsheets/rankine.toml`, a new
      golden (`rankine.json`, named in 0.6.0's entry), within 1.8e-10 kJ/kg
      of its hand calculation (`tools/gen_m5_rankine.py`); `reflux water`
      (`water_command.v`), six pairs (criterion 9 amended: the other nine
      refused with why), ice by R14-08's melting curves, the water report.
    - Step 8: `REFLUX_IAPWS95` (4), `reflux_psat` on a pure fluid's own
      model, `Database.fluid("IAPWS-95", ["water"])` (kij None by default),
      the reference program's water, bit for bit.
    - Found: the `iapws` package's (P, s) in the dome lands off the
      saturation line (0.15 K, 0.94 kJ/kg at 10 kPa), and its (P, x) stops
      1.5e-8 K high, so the hand calculation takes wet states by the lever
      rule at its own Newton's Tsat.
    - CI failed on bundle C's commit on Windows' gcc alone: `main`'s
      `fn finite` clashed with MinGW's math.h. Renamed `json_number`
      (`3b8f290`); CI passed, all four jobs (run 45984142988).
  - **Bundle D (steps 9 to 12) is done,** committed at `499d032`, CI
    passed on it, all four jobs (run 46032587308): M5a is finished, its
    exit review at the end of its plan. The owner approved criteria 6 and 9
    as amended on 2026-11-27. CI passed on the PR's tip, `7a0f3bb`, and it
    merged into `main` as `4bf65d5`.
  - Carried forward (in its review): M5b from `main` (mixtures' transport,
    IF97 behind the flash); `reflux water`'s other pairs (two of h, s and
    ρ); reports to DWSIM, the `iapws` package and V.
    - Step 9: `tools/gen_m5_dwsim.py`, `testdata/m5_dwsim.json`: DWSIM's
      CoolProp package is IAPWS-95 to 7.5e-7 kJ/kg, its PH flash 0.24–4%
      off at four near-critical states; its Steam Tables package is IF97
      where it answers, h = −1 in region 3 at 650 K, region 5 wrong, its PH
      flash above 25 MPa returning 1729/2000/1525 K whatever h (Brent on its
      own broken h, from its source); the Rankine cycle per kg to 3.7e-7.
    - Step 10: `tools/bench_m5.v`, `tools/gen_m5_bench.py`,
      `benchmarks/m5.v`: reflux 470–1400× DWSIM, 1.4–3× CoolProp from
      Python. Found and fixed: a pure PH/PS flash asked pvf above Pc, 8 ms
      to fail for water at 30 MPa; skipped above 1.01 Pc, the same bits
      (1,680 pure flashes before and after), PH 2.8 ms to 251 µs.
    - Step 11: 0.6.0's CHANGELOG entry written, undated; the release gate
      refuses for the tree, the branch and the date alone.
    - Step 12: docs, the roadmap (M5a done), and the exit review, which
      added criterion 2's Maxwell relations and ln φ = Gʳ/RT on IAPWS-95.
    - Planted mistakes: 6, all caught, 1 after a test was added. The full
      tcc check passed (138 test files).
- **M6 (activity models and liquid–liquid equilibrium) is done,** merged
  in PR #8 (`28724c1`) and released as v0.5.0, with its exit review at the
  end of its plan. Its plan, `docs/m6-plan.md`, was approved on 2026-11-18: exit
  criteria 1–13, decisions D1–D20 and steps 0–15, ending in v0.5.0.
  - **Step 0 is done** (committed at ea977c2), and amends D7 and criterion
    8, both approved on 2026-11-18:
    - The seam (`spike/seam`): `flash`'s PT flash behind a phase-model
      interface gives the same bits at all 10,000 sweep points, and the same
      time to the laptop's noise (median ratio 1.01).
    - DWSIM (`docs/spike-dwsim-activity.md`): its NRTL is τ = A/(1.98721·T)
      on ChemSep's banks to 4e-10, so the banks keep DECHEMA's R (D7),
      confirmed for NRTL alone and inferred for UNIQUAC and Wilson. Its
      activity packages take enthalpies from Lee–Kesler by default. Its
      UNIQUAC for methanol–water is 0.28 off in ln γ (explained in step
      12), and its three-phase flash seemed to miss NRTL's water–butanol
      split (corrected in step 12: the feed was outside it).
    - ThermoML (`spike/thermoml`) has a public search API and covers 2003
      to 2019, with nothing for methanol–water alone. So criterion 8 takes
      papers' tables where it has nothing.
  - **Steps 1 to 3 are done,** built together at the owner's word
    (committed at a70e0a9), and step 1 amends D7, approved on 2026-11-20.
    - Step 1: twenty ChemSep files vendored (the NRTL, UNIQUAC, Wilson and
      Henry's-law banks; each UNIFAC table with its extension), read in
      `components/activity.v`. The UNIFAC files contradict themselves in
      eight places (furfural's LLE subgroup on DMFA's main group, AC-CN's and
      AC-Br's swapped, a row's 61 for 62, Dortmund's unnumbered 230 and its
      `78 cy-CH2`), repaired by their own evidence; and DDBST's published
      tables settle 26 typos (D7 as amended). Each table against thermo's
      copy is in data-quality.md, every difference settled.
    - Step 2: `numeric.Jet` (forward differentiation), `Expr.jet` (τ(T)'s
      second derivative), and the `activity` module: Wilson, NRTL, UNIQUAC,
      each formula written once over f64 and jets, so its T- and
      composition derivatives are exact. Against mpmath at 40 digits: ln γ
      1.9e-15, dn 4.9e-14 of its row; against thermo, its own rounding.
    - Step 3: UNIFAC in three variants, a pure component's γ exactly 1.
    - Found: `thermo.r_gas` is N_A·k_B rounded, 1.84e-11 short of exact
      (made exact after M6, below); thermo's UNIQUAC divides by xᵢ.
    - Corrected after them (the plan, below step 3): a refused T left its
      parameters in the `Work` under the last T, so the next answer there
      was wrong; the comparison's rounding rule passed 342.4 for 300 (no
      count moved). D7 is amended, approved on 2026-11-22:
      UNIQUAC's unlisted pair is "no residual interaction", not ideal, and
      the banks' R is inferred for UNIQUAC and Wilson until step 12.
    - Then `activity.eval` was made to answer finitely or refuse: Σx within
      1e-12 of 1, G, τ, Λ and ψ normal doubles, every output asked for
      finite. `numeric.Jet`'s log no longer loses ∂² below 1e-154. No
      answer moved; newly refused only at 10 K and below.
  - **Step 4, the seam, is done,** and amends D3, approved on 2026-11-23
    (committed at e5c8055): `flash` takes a `thermo.PhaseModel`, a phase at
    (T, P, n), which an equation of state gives by `thermo.by_roots` and
    the γ–φ fluid gives itself.
    - Every cubic answer keeps its bits (criterion 4): the plan has the
      trace, the timings and what was compared.
    - What needs F(T, V, n) asks `Fluid.helmholtz()`: TV and UV, the
      reactors (D11, refused with why), equilibrium, Psat, a kij fit.
    - Found: PH and PS took a NaN energy balance as met (fixed); V's `as`
      to a type that can't be in the interface compiles and returns a zero
      value (v-for-go); `numeric.newton` drops a model's reason at x0
      (step 6's to pass on).
  - **Steps 5 and 6 are done,** built together at the owner's word
    (committed at 2871df2), and amend D4, D6 and criterion 3, approved on
    2026-11-23.
    - Step 5: `thermo.GammaPhi`, a γ–φ fluid behind the seam (a solvent by
      γ, Pˢᵃᵗ, φˢᵃᵗ at the correlation's Pˢᵃᵗ and Poynting; a Henry's-law
      solute by ln H mixed over the solvents with constants, γ* = 1),
      calorics from its fugacities, a workspace per flash
      (`thermo.Scoped`, `for_one_thread`); `Database.gamma_phi`;
      `tdep.Correlation.jet`. Against `thermo`'s `GibbsExcessLiquid`, ln φ
      to 4.7e-13, Cp to 3.6e-12, 7.9e-9 with φˢᵃᵗ; Henry's law in water and
      in a mixed solvent.
    - Step 6: the flashes on it (each split checked from scratch, against
      `FlashVL` β to 6.7e-11), `flash.azeotrope_t` and `azeotrope_p`,
      `numeric.newton`'s reason at x0, PH's step back from a refused T, and
      PH and PS from Newton's own state where a PT flash at its T misses
      the balance.
    - Every cubic answer keeps its bits but 868 trace lines: 865 refusals
      now with their reason, 3 refusals (one state's PH and PS) now its
      verified answer.
    - Found: ChemSep's hydrogen row in `henrys_law.ipd` is 5.9 times
      Sander's compilation at 298 K, slope 920 K against 490; `thermo`
      takes φˢᵃᵗ at its EOS's own Pˢᵃᵗ, and its solute's γ is the symmetric
      one.
  - **After step 6, Henry's constants by IAPWS,** at the owner's word to
    choose what's best for hydrogen, and amends D6, approved on
    2026-11-24: a gas in water that IAPWS's guideline covers
    (G7-04, 14 gases) takes IAPWS's own form and constants, in place of the
    bank's row. Against it, ChemSep's hydrogen is 5.95 times off at 298 K
    (62 at 623 K), its neon 4.5e32, the other eleven within 12% at 298 K
    (data-quality.md). Committed at 2648ebb.
  - **Steps 7 and 8 are done,** built together at the owner's word
    (committed at afb5cab). D9 and D14 are settled, approved on 2026-11-24 (by Π, 2,289
    of M3's sweep's 3,344 PR splits are two liquids, now liquid and liquid2;
    the C API refuses two liquids until ABI 2). D8 is amended, approved on
    2026-11-24: Michelsen's multiphase Rachford–Rice (Q over β ≥ 0) in
    place of Okuno's, Newton on G scaled by √n with each component's
    dependent phase its largest holder, a new phase begun at 1% of the
    feed, three liquids refused.
    - Step 7: `flash/multi.v`, the stage-wise split to three phases,
      `numeric.multiphase_rr`, D9's labels (`labelled`); against `FlashVLN`
      on 17 NRTL and PR states to 1.9e-13, each checked from scratch.
    - Step 8: PH and PS search again on verified flashes where the first
      misses the balance; TV and UV refuse three phases; the drum sends both
      liquids out of its liquid port; the C API and the sweep follow.
    - The sweep: 303 three-phase answers, 856 fourth phases or three
      liquids, 5 not solved (proved more than two, no convergence), 0
      unexplained.
    - From Python, D14 refuses 125 of the exit script's 6,480 map points
      (near the natural gas's critical point), one of M3's 39 splits by
      `thermo`, and the PT state of one TVF answer at 80.6 bar; their
      tests hold each so.
  - **Steps 9 to 11 are done,** built together at the owner's word
    (committed at d8a1e0c). D12 is amended, approved on 2026-11-24:
    Wilson's pairs are `lambda_ab` and `lambda_ba`.
    - Step 9: format 4 (a γ–φ `[model]`: `liquid`, `vapour`, `poynting`,
      `henry`, `missing`, `[[model.pairs]]` as expressions), the
      `decanter` (vapour, light, heavy; one liquid by light), a drum's
      two-liquid warning, report format 4's `model` (each pair's
      parameters and source). Every golden report written again, no number
      moved; `decanter.json` new.
    - Step 10: the ThermoML archive has methanol–water inside other papers'
      files after all (step 0 looked at two-compound files only):
      `testdata/m6_published.json`, 11 VLE sets at 101.3 kPa, water–1-butanol
      LLE and its heteroazeotrope, and a ternary LLE with p-xylene;
      `fit.point_test` (Fredenslund et al.'s, two starts for Barker's
      fit's two minima), 9 of 11 consistent; 74 model deviations recorded.
    - Step 11: a fit's `A[a, b]` (J/mol) and `alpha[a, b]` on NRTL and
      UNIQUAC, written back as `[[model.pairs]]`; against SciPy with
      thermo's models to 1.5e-6 of the standard errors. The published
      methanol–water set determines NRTL's two energies but not α nor
      UNIQUAC's, so those fits are on drawn ethanol–water data.
  - **Steps 12 and 13 are done,** built together (committed at 74787ff). No
    decision amended.
    - Step 12: DWSIM 9's activity packages by state
      (`testdata/m6_dwsim.json`, `tools/gen_m6_dwsim.py`), every difference
      explained, its γ worked out again from its own code where it isn't
      reflux's: NRTL and both UNIFACs are reflux's to 1.7e-9; its Poynting
      term has |P − Psat|; its methanol–water UNIQUAC mixes two records;
      its Wilson takes R = 1.9872 and volumes from the density at 298.15 K;
      its UNIFAC-LL is its VLE UNIFAC with wrong group fractions; its
      default three-phase flash is right, InsideOut and GibbsMinimization
      aren't. Four reports to make to DWSIM.
    - Step 13: γ–φ flashes timed (`testdata/m6_bench.json`,
      `tools/gen_m6_bench.py`, `tools/bench_m6.v`, `benchmarks/m6.v`):
      NRTL's PT flash 75 µs in reflux, 1.3 ms in thermo, 10.5 ms in DWSIM;
      half reflux's PR flash's time.
    - CI failed on d8a1e0c's Ubuntu jobs (both): the test of NRTL's α fitted
      where the data don't determine it, which stops at its rounding there
      (a gradient test of 2.3e-10) and converges on Windows. The test now
      takes either outcome; CI passed on 74787ff, all four jobs (run
      45477152064).
  - **Steps 14 and 15 are done,** built together (committed at f3fd1ca): M6
    is finished, its exit review at the end of its plan, every criterion
    met. CI passed on f3fd1ca (run 45496466774) and on the merge,
    `28724c1` (run 45501449501).
    - Step 14: v.mod and the package at 0.5.0; the CHANGELOG entry, undated
      until the tag on `main` after the merge; all eight golden reports
      written again, only the version moved. The release gate refuses for
      the tree, the branch and the date alone.
    - Step 15: the README, the roadmap (M6 done; C5, E4, H4's fits),
      data-format.md (the activity banks), the V note (no tuples in an
      array), and the criterion-by-criterion review: 317 planted mistakes,
      294 caught, 6 deleted, 17 explained.
    - Carried forward (in its review): the merge with `m5a-water`;
      `thermo.r_gas` (done after M6, below); γ–φ and three phases from
      Python (M10b, ABI 2); reports to DWSIM, `thermo`, ChemSep and V.
  - The data is ChemSep LITE's own: `nrtl.ipd`, `uniquac.ipd`,
    `wilson.ipd`, `henrys_law.ipd` and the UNIFAC tables, from
    `%LOCALAPPDATA%\ChemSepL8v50\ipd`, under the same licence assumption as
    `pr.ipd`.
  - `thermo` 0.6.1 is the reference: its activity models,
    `GibbsExcessLiquid` and `FlashVLN`. Experimental data come from NIST's
    ThermoML archive.
  - The roadmap gained Lee–Kesler–Plöcker in M12 on 2026-11-18,
    committed at ea977c2 with the plan.
- **M8 (specifications, sensitivity studies and optimisation) is done,**
  merged in PR #7 (`44fb1c7`) and released as v0.4.0.
  Its plan, `docs/m8-plan.md`, was
  approved on 2026-11-16: exit criteria 1–13, decisions D1–D20 and steps
  0–16, ending in v0.4.0.
  - One milestone, not split (D1): the optimiser is `numeric.sqp`, written
    in V (D12), so no NLopt or IPOPT. The roadmap is amended (step 16):
    IPOPT and data reconciliation move to M13 (B6, H4), and exact
    sensitivities too (H1, D10).
  - Its benchmark is the owner's 50-point reforming study; steps 1 to 8
    meet the roadmap's exit before the optimiser begins.
  - Research for the plan: four of the reforming case's five species give
    Cp as ChemSep's eq. 16 (521 of 640 compounds), which has no closed-form
    integral, so K(T) is a quadrature at every evaluation; threads with
    `spawn` ran the 48-point study 2.1× faster on 4 workers, processes 2.6×.
  - **Step 0 is done,** and amended D6 and D18, approved on 2026-11-17:
    the spikes (`docs/spike-dwsim-study.md`,
    `spike/study/spike.v`) and the `test_memory.py` fix.
    - DWSIM's study, Adjust, optimiser and regression have no API, but a
      script runs each form unshown ("DWSIM 9.0.5" below).
    - Warm starts by hand: a tear started from its neighbour's answer
      costs passes with Wegstein (11.0 against 4.0 cold) and Broyden (11.8
      against 8.1) on the reference flowsheet. So step 5 first measures
      what a warm start must carry, Broyden's Jacobian above all (D6).
    - Threads sharing one `Loaded` give one thread's bits: 218,200 numbers.
      4 threads run the owner's study about 1.6 times faster.
    - A PFR's solve allocates 20.6 kB a step, nearly all of it
      `Source.state`: 1,025 bytes a call. Step 7 starts there.
    - `reflux_collect` reports the bytes live after its collection
      (`gc_memory_use()`), for 0.4.0's CHANGELOG entry. `test_memory.py`
      holds them within 256 KiB: they stayed within 68.
    - V 0.5.2 tests an interface value for another interface with `is`,
      which D6's `Warm` needs.
  - **Step 1 is done,** and amended D11, approved on 2026-11-17: a rate's
    workspace keeps K until T moves, and at a new T integrates each species
    once for all its reactions; a PFR keeps what Wilke's rule needs at a T
    (`WilkeAt`).
    - Every answer keeps its bits: four checksums of studies, isothermal
      and adiabatic, against the code before, and the golden reports.
    - The published reforming model solves in about 20 ms, from 400; the
      owner's case 1.2 to 1.4 times faster.
    - Equation 16 by a fixed 24-node Gauss–Legendre rule holds 1e-13 for all
      521 compounds to 1500 K and would cut a fill at a new T from 17 µs to
      about 2 (`spike/eq16`). Whether to use it is the owner's to decide;
      the recommendation is to leave it until an adiabatic reactor needs it.
  - **Steps 2 to 4 are done,** built together, and amend D3, D4 and D5,
    approved on 2026-11-17 (committed at a43d8ea).
    - Step 2: dotted names in `expr`, `product.T`, `product.n[CH4]`,
      `reactor.duty`, `'let down'.T` quoted; `Expr.value`, a value or why
      there's none (D9); equations with `=`, `<=`, `>=`;
      `flowsheet/variables.v` and the optional `Energies` interface.
    - Step 3: a `Section` notes every number its maker reads;
      `Loaded.input(key)` and `Loaded.with(inputs, values)` remake a feed or
      unit by its maker, to the bit. A remake costs 8 to 80 µs a unit.
    - Step 4: format 3's `[spec.NAME]`, specs' links in planning, Brent for
      one spec alone, Broyden on tears and specs together with a nested
      retry, each spec checked again from the solution; report format 3.
      The reference flowsheet's D18 spec: 219.7574877219 K in 27 passes;
      an isothermal PFR's volume to 1.1e-11 of its closed form.
    - Found: two specs varying each other's targets need the whole first
      Jacobian; `with_history` dropped an error's line; a V bug, an empty
      `[]` as an `or {}` block's value.
  - **Steps 5 and 6 are done,** built together, and amend D6, approved on
    2026-11-17 (committed at 95743f7).
    - Step 5: `Loaded.solve_from(Start)` starts a point from its
      neighbour's solution; the optional `Warm` interface lets a unit start
      from its last answer (CSTR Newton from its last outlet, 10× faster;
      an energy balance's T bracket, 2×; the PFR's first step, nothing).
      Specs start from their last values (8 → 4.5 evaluations), tears only
      with direct substitution or specs: with Wegstein and Broyden a
      neighbour's tears cost passes, and carrying their Jacobian or slopes
      cost more (measured by `spike/warm`, then taken out).
    - Step 6: format 3's `[study.NAME]`, the `study` module (chunks of 8,
      a serpentine over two inputs, each point solved, no_value or failed
      with its reason), its own report format 1, and `reflux study`.
      Tutorial 2's methanol conversion reports why it has no value.
  - **Step 7 is done** (committed at 28150ec): `study.run(..., workers: N)`,
    threads sharing the file and taking whole chunks from a queue, the same
    bits on 1 to 8 workers; `reflux study --workers N`, the processors by
    default (`cpus_*.c.v`, not vlib's `runtime`, which tcc can't build on
    Windows).
    - `thermo.Cubic.work()`, a `solve_tp` for many states that allocates
      nothing, `solve_tp`'s bits; the reactors' `Source` reads its phase
      through it. A PFR solve allocated 9.72 MB and now 0.57; the owner's
      50-point study takes 0.71 s on 4 workers (1.74 s on one before).
      Step 0's four study checksums are unchanged.
  - **Step 8 is done,** and amends D10, approved on 2026-11-17
    (committed at 5166fc9): a study's `sensitivities`, each point's derivatives by
    central differences of five solves at 1e-12 (v, ±h, ±h/2, h = 1e-4·v),
    Richardson's extrapolation, and an error's estimate of the two
    differences' gap plus a fourth difference's noise over h (the
    amendment: the gap alone underestimated a spec's noisy derivative).
    Held against a PFR's closed forms by mpmath
    (`tools/gen_m8_sensitivities.py`), a heater's F·Cp from `thermo`, and
    Tutorial 2's dX/dT = −0.0018 /K; each estimate 3.7 to 10⁷ times its
    true error, recorded.
    - Found: a warm-started solve's error is in proportion to the step, a
      bias the estimate can't see; only the 1e-12 tolerance bounds it
      (1e-8 of the derivative). A recycle's derivative against cold solves
      at 1e-13 holds it.
  - Steps 1 to 8 meet the roadmap's exit.
  - **Steps 9 and 10 are done,** built together, and amend D12 and D13,
    approved on 2026-11-17 (committed at 05dd726).
    - Step 9: `numeric.qp` (Goldfarb and Idnani's dual active-set method)
      and `numeric.sqp` (Han and Powell's SQP: damped BFGS, an ℓ1 merit
      function, a QP with no point relaxed by δ, the first-order conditions
      checked at the answer). 21 of Hock and Schittkowski's problems reach
      their optima from their published starts, 206 iterations in all,
      each answer checked from scratch against mpmath's
      (`tools/gen_m8_hock_schittkowski.py`).
    - Step 10: format 3's `[optimize.NAME]`, `study.optimize` (the problem
      scaled; every evaluation a solve at 1e-12, warm; gradients by D10's
      differences, one-sided at a bound), its report format 1 and
      `reflux optimize`. Series reactions in a PFR reach their closed-form
      optima inside the bounds, on a bound and on a constraint, each
      multiplier within 1e-5 (`tools/gen_m8_series.py`).
    - Found: a one-sided difference's estimate, R₁'s, had no margin (1.01
      times its error), so its value is R₂; a relaxed QP at a badly
      conditioned iterate can fail, and is then reported as infeasible.
  - **Steps 11 and 12 are done,** built together at the owner's word, and
    amend D14 and D15, approved on 2026-11-18 (committed at 579c12c).
    - `numeric.least_squares` (Levenberg–Marquardt on QR, Moré's scaling,
      Nielsen's damping, the gradient test at 1e-10) and
      `numeric.student_t_quantile`; format 3's `[fit.NAME]`
      (`flowsheet/fits.v`), the `fit` module (rate constants on the tape's
      exact derivatives, centred on the data's mean T; kᵢⱼ on bubble
      points), its report format 1 and `reflux fit`.
    - Criterion 9: seven LHHW and power-law constants and a kᵢⱼ match
      SciPy's fits of the same data (`tools/gen_m8_fits.py`), values to
      1e-4 of their standard errors; over 400 seeded data sets the 95%
      intervals hold the truth 381 and 382 times, against 380.
    - Found: Levenberg–Marquardt stalled near a minimum, its fall below
      Σr²'s rounding, so a step there is judged by the gradient test;
      SciPy's t.ppf is off by 1.2e-12 at ν = 1, so the quantiles' reference
      is mpmath's.
  - **Steps 13 and 14 are done,** built together at the owner's word, no
    decision amended (committed at 6c2cf76, with steps 15 and 16).
    - Step 13 (`tools/gen_m8_dwsim.py`, `study/dwsim_test.v`): the
      reference flowsheet's spec on DWSIM's constants 2.8 mK from its
      Adjust's, explained to 1e-6 K by the vapour's 1.6e-3 mol/s at DWSIM's
      answer over its slope; the compression's least work 1.5 Pa from its
      optimiser's, explained by reflux's slope at DWSIM's optimum; the
      series reactions in DWSIM's PFR, whose optimum is 1.5% off from noise
      in its outflow. The kij fit against DWSIM's regression stays not
      possible from a script (step 0).
    - Step 14 (`tools/gen_m8_bench.py`, `tools/bench_m8.v`,
      `benchmarks/m8.v`): the owner's study 904 s by DWSIM's own form,
      0.70 to 1.24 s by reflux; the spec 130 times faster, the optimisation
      4.8; the published model 13.4 ms a solve.
  - **Steps 15 and 16 are done,** built together (committed at 6c2cf76):
    M8 is finished, with its exit review at the end of its plan. CI passed
    on 6c2cf76, all four jobs (run 44337209227).
    - CI first: the branch's one run (`05dd726`, run 44154785483) failed
      both gcc jobs. Step 1's test counted the collector's bytes to show a
      fill's integrals, and gcc's build of the collector counts small
      allocations in batches, so CI's counted 0 (a probe: 939 of 1,000 small
      arrays moved the count by nothing under gcc, none under tcc). A
      `kinetics.Point` now counts its `integrals` itself.
    - Step 15: v.mod and the package at 0.4.0; the CHANGELOG entry, dated
      2026-11-18 on `main` after the merge; all seven golden reports
      written again in report format 3, with no number moved.
    - Step 16: the README, the roadmap (D1's amendments: B6's IPOPT, H1's
      exact sensitivities and H4's data reconciliation to M13, which lists
      them), design.md's exact sensitivities to come, the V note, and the
      criterion-by-criterion review: every criterion met but CI and the tag,
      and the kᵢⱼ fit against DWSIM's regression, not possible from a script.
  - Carried forward (in its review): equation 16 by a fixed rule, the
    owner's to decide; M13's exact sensitivities, infeasible path, IPOPT
    and data reconciliation; fits to an integral reactor's outlets; a kᵢⱼ
    fit against DWSIM's regression if the owner's app run gives one;
    workers bounded by the collector's lock; studies, optimisations and
    fits from Python (M10b); reports to V, DWSIM and SciPy.
- **M3, the flash,** carried forward (in its review): VLLE is M6; reports
  to make upstream to V, `thermo` and DWSIM. Its PH speed was M4 step 9's,
  which made it 2–4× faster.
- **M4, the first flowsheet,** carried forward (in its review):
  - `flash.uv` fails on cold liquids at or below 100 K, in its TV flash. No
    unit uses UV yet; fix it before one does.
  - PH for 10 components is still a little slower than DWSIM's (5.9 ms
    against 5.2).
  - Warm starts across solves (F5) came in M8 (step 5). A solver callback
    for `--trace` as it happens is still to do.
  - Reports to make upstream: to V, the TOML parser's panic on a `\u`
    escape and a smart-cast interface's field inside `${}`; to `thermo`, its
    flash at a vapour fraction; to DWSIM, its two-inlet drum and its silent
    unit failures.
- **Released:** v0.1.0 (M0 to M4), tagged on 2026-11-03 at `45d4ce6`,
  v0.2.0 (M10a), tagged on 2026-11-08, v0.3.0 (M9), tagged on 2026-11-15,
  v0.4.0 (M8), tagged on 2026-11-18 at `18eb317`, v0.5.0 (M6), tagged
  on 2026-11-25 at `862cd2c`, and v0.6.0 (M5a), tagged on 2026-11-27, all
  on `main`: source releases, annotated tags, pushed. The
  next release starts with a new CHANGELOG entry and its version in `v.mod`
  and the package, and a changed golden report goes in that entry by name.
- **`test_memory.py` was flaky on Ubuntu,** and M8's step 0 fixed it, now
  on `main`. It held the collector's `heap_size`, free blocks included,
  which Boehm grows in one large step when it chooses to and never gives
  back: on `main`'s merge of M10a it doubled once and failed. It now holds
  the bytes live after a full collection.
- **M9 (equations, kinetics and reactors) is done,** merged in PR #6 and
  released as v0.3.0, with its exit review at the end of its plan. Its plan,
  `docs/m9-plan.md`, was approved on 2026-11-08: exit criteria 1–11,
  decisions D1–D16 and steps 0–13, ending in v0.3.0.
  - The owner's principle drives it: no equation is ever taken as text, not
    only LHHW rates (roadmap A9). Users type their own, with species by
    name, numbers with units, named definitions, and built-in constants at
    the SI's exact values (`R`, not a typed 8.314).
  - Its benchmark is the owner's reforming PFR (see the memory), whose DWSIM
    file gives the three reactions three different denominators.
  - Heat exchangers (G2) moved out of M9 to a later milestone (D14).
  - **Step 0 is done:** the spikes on DWSIM's PFR
    (`docs/spike-dwsim-pfr.md`) and on the tape (`spike/tape/spike.v`).
    It amended D5 (the tape kept, with fused operations: 3 to 5 times
    hand-written V, under 2 µs for three rates and their gradients), D10
    (Ergun needs a gas viscosity, so ChemSep's correlations mixed by Wilke's
    rule come forward from M5, and the PFR takes a viscosity through an
    interface) and D13, and set criterion 8's targets. The owner approved
    the amendments on 2026-11-09.
  - **Step 1 is done:** `quantity`'s dimension algebra. `Dim` (exact
    rational exponents of kg, m, s, K, mol, A), `unit(text)` for compound
    units with SI prefixes, and `parse` reading any unit of a dimension.
    Every M4 unit reads to M4's bits, against a frozen copy of M4's table,
    and the new units match SP 811 (`testdata/m9_units.json`). D3 amended:
    current as a base, rational exponents, symbols rather than M4's whole
    names, and a * after a / refused as ambiguous.
  - **Step 2 is done:** the `expr` module, the expression language
    (`docs/design.md`, "Expressions"). Expressions parse to arena trees,
    names bind to the context, definitions and the SI-exact constants, every
    node's dimension is checked, and each writes back as plain text (which
    reads back to the same tree), Unicode and LaTeX. Held against SymPy by
    `testdata/m9_expressions.json` from `tools/gen_m9_fixtures.py`, which
    needs V on the PATH (it runs `tools/expr_render.v`) and decodes reflux's
    renderings itself: SymPy's LaTeX readers misread function calls. D2
    amended: units after a space or in brackets, "2 R" refused as
    ambiguous, offset units and `log` refused, `·` `×` `−` accepted.
  - **Step 3 is done:** `expr/tape.v`, expressions compiled onto one tape
    (folded, merged, fused, pruned), evaluated into the caller's `Buffer`
    with nothing allocated, and differentiated by one reverse sweep for
    every output. Every fusion keeps the tree's arithmetic, so the tape
    gives eval's bits; derivatives match SymPy's (now in
    `testdata/m9_expressions.json`) to 6.5e-16. Xu and Froment's three
    rates: 55 operations, 3.2 times hand-written V for the rates and 6.7
    with gradients, about 1 µs in all (`tools/tape_bench.v`).
  - **Step 4 is done:** a new `chem` module (elements, formulas, SMILES;
    D1 amended, approved on 2026-11-10) and `kinetics`' first part:
    species by name, CAS number or formula; equations balanced by elements
    and charge, with the balanced equation suggested; a set's dependence; ΔH°, ΔS°, ΔG° and K(T)
    from the ideal gas (D7). A compound's formula is the one two of its
    structure, SMILES and molar mass agree on (640 of 642; five files
    disagree with themselves). Held against RDKit, chemicals, thermo and
    mpmath (`tools/gen_m9_reactions.py`). Found: ChemSep's S° mixes 1 atm
    and 1 bar (reforming's K about 4% low); 142 compounds' ΔG°f disagree
    with their ΔH°f and S° by over 1 kJ/mol (data-quality.md); and the
    owner's DWSIM case types an overall-reaction K that isn't the product
    of the other two (86 times too small at 700 K).
  - **Step 5 is done,** and amended criterion 5 (approved on 2026-11-10):
    `kinetics/law.v` (the forms as data), `kinetics.v` (each form writes the
    expression a user would type, parsed by `expr`: D4's same tree) and
    `evaluate.v` (rates at a state, allocating nothing, and the jacobian by
    the chain rule). A rate's context has `Keq[...]`, each reaction's K at
    1 bar (not `K`, which clashes with kelvins). Xu and Froment's and Vanden
    Bussche and Froment's rates match mpmath (`tools/gen_m9_kinetics.py`)
    to 9.4e-15, every constant checked against the papers' tables (the
    owner's copies, kept out of the repository). Xu and Froment's Table 5
    and Table 6 disagree by up to 1.7% (k₃); their paper gives no
    equilibrium constants, so criterion 5 compares Gibbs K with the
    literature's correlations instead.
  - **Step 6 is done,** and amended D11 (sites and reactions as tables by
    name, `[site.NAME]` and `[reaction.NAME]`, approved on 2026-11-10):
    `flowsheet/reactions.v` reads a file's reactions (format 2, alone or in a
    flowsheet file), `kinetics/check.v` the checks (notes and warnings),
    `kinetics/explain.v` the factor table, and `reflux kinetics check` and
    `explain` (`kinetics_command.v`). Three files in `testdata/kinetics/`: the two
    papers' kinetics and the owner's DWSIM case's rates as typed, whose
    check finds its seven 8.314s, its two cut denominators and its overall
    K a factor of 0.20 off at 900 K.
  - **Step 7 is done,** and amended D1 and D10, approved on 2026-11-11:
    a new `equilibrium` module. Gibbs minimisation by the RAND method in
    CEA's form (Gordon and McBride), species in log space, ideal gas or a
    cubic with φ updated outside it; and the equilibrium over chosen
    reactions by the same minimisation over what they conserve, with
    stand-in potentials −νᵀ(ννᵀ)⁻¹ ln K (D10: not Newton on extents). Every
    answer is checked (balance; stationarity, or each ln Q its ln K; a
    stable phase, two refused as M6's) and matches mpmath's
    (`tools/gen_m9_equilibrium.py`, 21 cases) to 6.0e-13 in ln y.
  - Step 7's consistency check (`equilibrium/consistency.v`; D1:
    `equilibrium` imports `kinetics` and `flash`) follows each reversible
    rate along its reaction from the equilibrium a state reaches to where
    it's zero, and `reflux kinetics check` runs it given a state. Xu and
    Froment's stop within ΔG/RT 2.7e-15 of equilibrium; the DWSIM case's
    overall rate stops at Q/K = 0.20. Found: an equilibrium over reactions
    from a feed without an element was a singular matrix, fixed by looking
    at the element rows too.
  - **Step 8 is done,** and amended criterion 6, approved on 2026-11-12:
    `numeric`'s integrators, `dopri5` (Dormand and Prince's 5(4) pair,
    with Hairer's stiffness detection, which refuses a stiff problem) and
    `rodas4` (Hairer and Wanner's L-stable Rosenbrock method, on the
    problem's own Jacobian), both as Hairer's codes step and control them,
    taking the problem as an interface (`Ode`, `StiffOde`). Measured: each
    solution's order, RODAS4's |R(z)| ≤ 1, closed forms, and the Test Set for
    IVP Solvers' stiff problems (`tools/gen_m9_odes.py`) within 0.13
    tolerance units, the Oregonator's within 2.9, steps recorded.
  - Found in step 8: E5's y₂ and y₃ are beyond RODAS4 in double precision
    (51% off under tcc, 20 times under gcc; a Rosenbrock stage's rounding is
    relative to f's terms, which cancel there), so its test holds its decay
    and its invariant's drift instead; and
    an orbit's error grows to many tolerances, so criterion 6's "within the
    tolerance asked" is read as errors in tolerance units, recorded.
  - **Step 9 is done,** and amended D11, approved on 2026-11-12: a new
    `reactors` module (Gibbs, equilibrium, conversion, CSTR, PFR) with makers
    for format 2 files, each reactor's `reactions` naming what it runs. Every
    reactor is isothermal, adiabatic or at a duty; the PFR's energy equation
    carries the real gas, Ergun's drop takes Wilke's viscosity (`chemicals`
    to 2.2e-16), and its profile goes in report format 2. A unit that reacts
    reports its elements (`Outcome.conserved`), and the solver checks it, and
    the whole flowsheet, on them. A long PFR on fugacities reaches the Gibbs
    reactor's equilibrium to 2.8e-14.
  - Found and fixed in step 9, at the owner's word: M3's PT flash labelled a
    hot gas liquid (Π − 1 ≈ (B − T·dB/dT)/V > 0 above the Joule–Thomson
    inversion temperature). M3's D6 is amended, approved on 2026-11-12: a
    single phase more than four covolumes large is vapour. 17 supercritical fixture roots now
    differ from `thermo`'s labels.
  - **Step 10 is done,** and amended D13, approved on 2026-11-13. The
    owner's Tutorial 2 (methanol in a conversion reactor; its −∞ sweeps
    plotted a product's conversion) amended D2 and D10 before it, approved
    on 2026-11-12:
    - `expr` reads `(T in degC)`, a value's number on a unit's scale, and
      a unit after a group, `(…) [%]`;
    - a conversion may be an expression in T, taken at the outlet's T;
    - a conversion reactor warns of an outlet past equilibrium, in
      `Outcome.warnings` and report format 2's `units[].warnings`.
  - Step 10's comparisons (`tools/gen_m9_dwsim.py`, `testdata/m9_dwsim.json`,
    `reactors/dwsim_test.v`):
    - The owner's PFR is within 2.3e-10 of the exact answer of DWSIM's
      equations, and within criterion 8's targets of DWSIM at dV = 0.001.
    - Tutorial 2's duties are explained to 3.2e-13.
    - DWSIM's equilibrium reactor puts fugacities over 1 atm against a
      1-bar K.
    - Its Gibbs reactor is good to about 1e-3.
    - Its CSTR reacts nothing in the vapour, so reflux's CSTR is held to
      the exact answer of DWSIM's equations instead (D13).
  - Found in step 10: reflux's CSTR couldn't start Newton on the owner's
    case, so it now solves by continuation in its volume where Newton from
    the inlet fails. The reactor goldens are in
    `testdata/flowsheets/golden/reactors/`, and `tools/golden.vsh` takes
    names, so that M4's goldens wait for v0.3.0.
  - **Step 11 is done:** speed, against DWSIM in one session, taking turns
    (`tools/gen_m9_bench.py`, `tools/bench_reactor.v`,
    `testdata/m9_bench.json`, `docs/benchmarks.md`). The owner's case
    (`testdata/flowsheets/reformer_dwsim.toml`): reflux 28 ms a solve against
    DWSIM's 10.7 s, and 2.1 s for the 50-point study against 11.7 min. The
    published model (`reformer.toml`) takes 393 ms, 93% of it in K(T), which
    `Kinetics.fill` recomputes by quadrature at every evaluation: left for
    M8, at the owner's word, with the viscosity Ergun's drop recomputes (the
    roadmap's M8 entry). The tape is about 10% of a solve, so it stays (D5).
  - **Steps 12 and 13 are done:** v.mod and the package at 0.3.0, the
    CHANGELOG entry undated until the tag, the seven golden reports written
    again (report format 2, no number moved), and the docs and the exit
    review. No source changed, so no full checks ran; CI passed on step 11.
  - Carried forward (in its review): `test_memory.py` (M8 step 0); K(T)
    and the viscosity at an unmoved T (M8); one-phase reactors (M6); the
    yield reactor; reports to DWSIM, `thermo` (a hot gas's label), SymPy and V.
  - Next: M8, chosen on 2026-11-16 (above).
- **M10a (thermo and the flash from Python) is done,** merged in PR #5:
  steps 1–10 of `docs/m10a-plan.md`, whose exit review is at its end. CI
  passed on its tip on Windows and Linux (run 42080776945), after the run
  before found Ubuntu's gcc fusing multiply-adds.
  - Its probes found that V 0.5.2's collector has thread support on Windows
    too, so the spike's thread registration is needed on both OSes. 8 Python
    threads at once get the serial answers bit for bit, at about 2× the serial
    rate.
  - **Step 1 is done:** the library skeleton (`capi`, `capi/reflux.h`, the
    registry, thread registration, exports restricted to the header's,
    `tools/python.vsh`, the package's skeleton). It amended D3 (a version
    script on Linux) and D5 (V starts the library on its first call, which
    mustn't race, so the package makes it at import).
  - **Step 2 is done:** databases and fluids through the C API and the
    package, and the data root (D9). `REFLUX_DATA` now names the data root
    (`components/` and `vendor/chemsep/`), and `components.load` became
    `open(root, extra)`: both go in 0.2.0's CHANGELOG entry, listed at the
    end of step 2's write-up. The package installs (D10), editable or not.
  - **Step 3 is done:** a phase's properties and each compound's Psat at many
    points, each with a status (`REFLUX_POINT_*`, 0 to 8). A point whose
    answer isn't finite has none (`NOT_FINITE`): above about 1e10 K the
    ideal-gas Cp overflows.
  - Carried forward from step 3: thermo's `saturation_pressure` (successive
    substitution) fails at Tc(1−1e-6) for n-butane and at or below 10 K.
  - **Step 4 is done:** stability and the PT flash at many points, through
    one `reflux_flash(fluid, spec, …)` whose layout every flash will share.
    D6 is amended: each phase is in the place of its label.
  - **Step 5 is done:** PH, PS, TV, UV and the vapour-fraction flashes, with
    bubble and dew points, as specs of `reflux_flash`. Every failure kind is
    reached and tested through the public API.
  - **Step 6 is done:** seeded fuzzing of every exported function
    (`capi/fuzz_test.v`), 8 threads on every vectorised function, and a flat
    heap over 100,000 calls. Throughput with 8 threads is 1.3 to 2.3 times
    serial, the collector's lock being shared.
  - **Step 7 is done:** `tools/capi_reference.v` (D12) writes the V API's
    answers, and the Python tests match them bit for bit through the C API.
    `docs/benchmarks.md` has the cost from Python (criterion 10). The spike's
    code and CI step are gone, every scenario living in the package's tests.
  - Carried forward from step 7: V's code runs 1.3 to 1.8 times slower inside
    the library in a Python process than in a V program, and fewer
    collections don't remove it. For M8.
  - **Step 8 is done:** the exit script, `python/examples/phase_diagrams.py`
    (a T–x–y and a phase-fraction map, matched to the V API bit for bit), and
    `docs/python.md`, whose every example `python/tests/test_docs.py` runs.
    The script found a bug in M3's PT flash, which is fixed: the split's
    extrapolation is kept only if it lowers the residual. It goes in 0.2.0's
    CHANGELOG entry.
  - **Step 9 is done:** v.mod and the package are 0.2.0, with its CHANGELOG
    entry undated until the tag, which comes on `main` after step 10 and the
    merge, as v0.1.0's did. The golden reports carry 0.2.0, and their test
    now holds each to v.mod's version. No number in them moved.
  - **Step 10 is done:** the README, roadmap and V notes, and the exit
    review. The review added a test of the header's types against ctypes'
    (criterion 1), and Mathias–Copeman constants to the install test
    (criterion 9).
  - Carried forward (in the review): V's code slower inside the library and
    threads that don't scale (M8); thermo's `saturation_pressure` near Tc
    and below 10 K; the stability test's unjudged extrapolation, to watch;
    wheels (I2); the committed figure, drawn by hand; a report to V on a
    shared library's first call.
  - M10b (flowsheets from Python) follows; M5 to M8 build on M4.
- **The owner's own benchmark,** for M8 and M9: a one-PFR methane steam
  reforming flowsheet, 15 to 18 s a solve in DWSIM, in `docs/roadmap.md`.
- **One gas constant, after M6,** at the owner's word on 2026-11-25, on
  branch `r-exact`, for its PR into `main`: R is `quantity.molar_gas_constant`,
  N_A·k_B = 8.31446261815324, the SI's exact value, and `thermo.r_gas`,
  `activity.r_gas`, `kinetics.gas_constant`, an expression's `R` and the
  fits' are that one constant (`fit/gas_constant_test.v`). The first two
  were 8.314462618, 1.84e-11 short, and a Gibbs reactor or PFR mixed both.
  - It begins 0.6.0's CHANGELOG entry, with v.mod and the package at 0.6.0:
    the next milestone adds to it.
  - The eight golden reports moved by at most 3.2e-11 of each number (the
    reformer's PFR step positions by 7.9e-8, the same 472 steps); no label,
    warning, pass or step count changed.
  - The flash sweep moved at five cold heavy mixtures whose three-phase
    split is at the edge of converging: 4 not solved, from 5; 857 fourth
    phases or three liquids, from 856; none unexplained.
  - The full check then found an M6 bug: `numeric.multiphase_rr` left a
    leaving phase at β·ε each step, not 0, until its step underflowed (20
    steps, water–1-butanol by its heteroazeotrope, `published_test.v`).
    The phase bounding the step now goes to exactly 0, held by a test on
    the failing doubles; no golden or sweep outcome moved.
  - Fixtures whose generators typed reflux's old R were made again on the
    exact one: `m6_activity.json` (mpmath's Gᴱ, Hᴱ, Cpᴱ, by 1.84e-11),
    `m6_unifac_tables.json`, `m1_coolprop_butane.json` (Sʳ's R ln Z),
    `m1_thermo.json` and `m9_equilibrium.json` (no number moved). Two
    Python tests typed it too.
  - Tests that held old rounding: `gamma_phi_test.v` counted only nonzero
    worsts, so its now-exact ∂P keys went missing (four limits tighten to
    1e-13); a vdW Cvʳ difference held to an exact 0 (now its rounding,
    ε|Uʳ|/h); the heater's sensitivity margin, its true error rounding
    (floor 5, measured 10.6).
- **The owner's CAT 1, an SMR–SOFC flowsheet,** offered on 2026-11-25 as a
  reference flowsheet, goes into the next milestone's plan (the owner's
  word), as a step. It has a Gibbs reformer, conversion reactors for the
  fuel cell and afterburner, and a water recycle. A probe solved it to
  within 0.3 kW of DWSIM's duties.
  - Found: a conversion reactor warns that a complete combustion is past
    equilibrium (Q/K infinite), a false alarm to fix in that step.
  - The owner's DWSIM file and the paper stay out of the repository. Whether
    the splitter takes a flow spec, as DWSIM's does, is the owner's to
    decide.

## How we work

- **Plan a milestone before coding it.** Research, then write
  `docs/mN-plan.md` with exit criteria, decisions (D1, D2, …) and steps. It
  needs the owner's approval before any code.
- **Steps in bundles.** Related steps are built together by default, each
  with its own write-up in the plan, and stop where a decision is the
  owner's. The owner says when to commit and what to start next.
  - Commit only when asked, and push only when asked.
  - Open the PR into `main` only after the milestone's docs and exit review,
    and only when asked.
  - Tag a release only when asked, on `main`, once `v run tools/release.vsh`
    passes ("Making a release" in `CHANGELOG.md`). A change to a golden
    report goes in the CHANGELOG's next entry, by name.
  - Commit subjects look like `M3 step 4: PT flash`, with a body below.
- **Finishing a step** means all of:
  - marking it "— done" in the plan, with what was built, the tests,
    measured numbers, bugs found and the planted-mistake results;
  - updating `docs/design.md` (every formula the code uses), the README
    layout block, `testdata/README.md` (one row per fixture) and
    `docs/v-for-go.md` (V lessons);
  - at the end of a bundle, before its commit and only if V source
    changed, one full tcc check. No local gcc check: CI runs both
    compilers, the `-prod` build and the Python tests on Ubuntu and Windows
    on every push, so its run is read after the push. Locally, gcc runs
    only what tcc can't show: `tools/python.vsh test` when the C API or
    what it calls changed, and a test that tcc and gcc may differ on.
    Inside a bundle, and for a fix between bundles, run the touched tests
    and their dependents, and the planted mistakes.
- **Changing an approved decision** happens in the plan's table, marked
  "Amended in step N", and is flagged to the owner.

## Commands

```bash
v run tools/check.vsh            # what CI runs, with tcc: about 4 min
v run tools/check.vsh -cc gcc    # with gcc, plus a -prod build: about 8 min
v -W test flash/pt_test.v        # one test file
v -W -stats test flash/pt_test.v # per-test results, and println output from passing tests
v run . solve testdata/flowsheets/reference.toml --json out.json   # the CLI (docs/report-format.md)
v run . study testdata/flowsheets/methanol_study.toml              # a study: its points, and why an output has no value
v run . optimize testdata/flowsheets/series.toml best              # an optimisation: its optimum, multipliers and iterates
v run . fit testdata/m8_fit/ethylene.toml                          # a fit: each constant with its interval, written back
v run . water T=450degC P=40bar [--if97] [--json out.json]         # a steam table: water at a state, by IAPWS-95 or IF97
v run . demo                     # the demo; `v run . demo 300 2.5` sets T (K) and P (bar)
v run tools/release.vsh          # the release gate: why this commit can't be tagged, or the tag commands
v run tools/python.vsh test      # the C API's library (gcc), then the Python tests (uv, pytest)
REFLUX_FUZZ_SEED=3 REFLUX_FUZZ_CALLS=20000 v -cc gcc -prod test capi/fuzz_test.v   # a longer fuzz
```

- **`check.vsh`** runs `v fmt -verify .`, `v vet .`, `v -cc <cc> -W test .`,
  then builds reflux and runs `reflux demo`, solves the reference flowsheet
  with both reports, and checks that a missing file exits 1. With gcc it
  does the same with a `-prod` build.
- **Test builds are the slow part:** each `_test.v` file is compiled into a
  program of its own. `check.vsh` runs up to 4 at once (V's runner would
  allow one per 8 GiB of RAM: one here), and builds with `-nocache` under
  gcc, where V 0.5.2's module cache made each build about three times
  slower: about an hour on this laptop by M6. Run any full check in the
  background, and don't edit or regenerate anything it compiles or reads
  until it finishes.
- **CI** (`.github/workflows/ci.yml`) runs `check.vsh` on Ubuntu and Windows
  with both tcc and gcc. The gcc jobs also run `tools/python.vsh test`. V is
  pinned to the commit in `.v-version` (0.5.2, eab77e7).
- **The Python tests** need the library, which only gcc can build, so
  `check.vsh` doesn't run them: `v run tools/python.vsh test`. It also builds
  `tools/capi_reference.v` with the library's compiler and flags, for the
  bit-for-bit tests (`REFLUX_REFERENCE`); pytest alone skips those. Both
  are built with `-ffp-contract=off`: where gcc's target has FMA, as
  Ubuntu's CI gcc's does, it fuses a·b + c by its inlining, and the two builds
  inline differently (M10a's exit review). Here uv needs
  the interpreter named, as below: `--python <path>`. `test_install.py`
  installs the package into a fresh environment, so uv needs setuptools and
  NumPy from its cache or the network. matplotlib (pinned) is for the exit
  script, whose test runs it. An editable install leaves
  `python/reflux.egg-info/`, which git ignores.
- **Generated data:**
  - `v run tools/import_chemsep.vsh` regenerates `data/components/` from the
    vendored ChemSep XML. A test fails if it's stale.
  - `v run tools/data_quality.vsh` writes `docs/data-quality.md`.
  - `tools/gen_m3_bench.py` (by hand, with DWSIM) times reflux, thermo and
    DWSIM case by case and writes `testdata/m3_bench.json`; `v run
    tools/benchmarks.vsh` then writes `docs/benchmarks.md`.
    `tools/gen_m4_bench.py` does the same for the reference flowsheet,
    into `testdata/m4_bench.json`. Rerun them after a change that could
    make flashes faster or slower. `v -prod -cc gcc run
    tools/benchmarks.vsh --time` times reflux alone, against the fixtures,
    without DWSIM.
  - `v -prod -cc gcc run tools/flash_sweep.vsh` writes
    `docs/flash-robustness.md`, 10,000 PT flashes (about three minutes,
    the build included; about twenty under tcc). Regenerate it after
    changing the flash or the sweep's checks.
  - `v run tools/golden.vsh` writes the golden flowsheets' reports
    (`testdata/flowsheets/golden/`). Run it only when a change is meant to
    move the answers, and explain why in the commit.
- **Reference fixtures** come from uv scripts whose headers pin their
  dependencies: `thermo` 0.6.1, `chemicals` 1.5.2, `mpmath` 1.3.0,
  `CoolProp` 8.0.0, and for `tools/gen_m4_graphs.py`, `networkx` 3.4.2 and
  `scipy` 1.14.1.
  - On this machine uv's managed-Python link fails in the sandbox, so pass
    the interpreter explicitly, and set `PYTHONIOENCODING=utf-8` (the
    console is cp1252):

    ```bash
    PYTHONIOENCODING=utf-8 uv run --python "C:/Users/Ricky/AppData/Roaming/uv/python/cpython-3.13.17-windows-x86_64-none/python.exe" tools/gen_m3_fixtures.py
    ```

  - `tools/gen_m6_henry.py` (the `iapws` package) writes
    `testdata/m6_henry_iapws.json`, IAPWS's Henry's constants for gases in
    water; `tools/gen_m6_gamma_phi.py` (thermo, chemicals, mpmath, SciPy)
    `testdata/m6_gamma_phi.json`, the γ–φ liquid, Henry's law and the
    flashes against `thermo`.
  - `tools/gen_m6_published.py` (SciPy) writes `testdata/m6_published.json`,
    data sets fetched from the ThermoML archive by DOI (cached in
    `$REFLUX_THERMOML_CACHE`), with each VLE set's point test;
    `tools/gen_m6_fits.py` (thermo, SciPy, mpmath) `testdata/m6_fits.json`
    and `testdata/m6_fit/*.csv`, SciPy's fits of NRTL's and UNIQUAC's pairs.
  - `tools/gen_m6_dwsim.py` (DWSIM, about two minutes) writes
    `testdata/m6_dwsim.json`, DWSIM's activity packages by state;
    `tools/gen_m6_bench.py` (DWSIM, thermo and V, a few minutes)
    `testdata/m6_bench.json`, the γ–φ timings, which can't reproduce byte
    for byte. Then `v run tools/benchmarks.vsh` writes `docs/benchmarks.md`.
  - `tools/gen_m6_activity.py` (thermo, mpmath) writes
    `testdata/m6_activity.json`, the activity models on ChemSep's banks and
    tables, and `testdata/m6_unifac_tables.json`, thermo's UNIFAC tables.
  - `tools/gen_m8_sensitivities.py` writes `testdata/m8_sensitivities.json`
    (mpmath only): a PFR's derivatives for M8's criterion 6.
  - `tools/gen_m8_hock_schittkowski.py` (SymPy and mpmath) writes
    `testdata/m8_hock_schittkowski.json`, and `tools/gen_m8_series.py`
    (mpmath) `testdata/m8_series.json`: criteria 7 and 8.
  - `tools/gen_m8_fits.py` (SciPy, NumPy, mpmath, thermo) writes
    `testdata/m8_fits.json` and the CSV files in `testdata/m8_fit/`:
    criterion 9's data, SciPy's fits of them, and t quantiles.
  - `tools/gen_m8_dwsim.py [spec compression series]` (DWSIM, about twenty
    minutes) writes `testdata/m8_dwsim.json`, DWSIM's answers for step 13;
    `tools/gen_m8_bench.py` (DWSIM and V, about twenty minutes)
    `testdata/m8_bench.json`, step 14's timings, which can't reproduce byte
    for byte. Then `v run tools/benchmarks.vsh` writes `docs/benchmarks.md`.
  - `tools/gen_m5_iapws95.py [coefficients derivatives states saturation]`
    (iapws, CoolProp, mpmath) writes `testdata/m5_iapws95.json`, and
    `tools/gen_m5_if97.py [coefficients states]` (iapws, with NumPy 2.2.6:
    under 2.5.5 its region-3 (p, h) raises) `testdata/m5_if97.json`.
  - `tools/gen_m5_transport.py` (iapws, NumPy 2.2.6) writes
    `testdata/m5_transport.json`, the package's IAPWS-95 and IF97
    transport at (T, P) states and σ; `tools/gen_m5_rankine.py` (iapws)
    `testdata/m5_rankine.json`, the Rankine cycle by hand.
  - `tools/gen_m5_dwsim.py` (DWSIM, about a minute) writes
    `testdata/m5_dwsim.json`, DWSIM's two water packages by state and the
    Rankine cycle on each; `tools/gen_m5_bench.py` (DWSIM, CoolProp 8.0.0,
    iapws and V) `testdata/m5_bench.json`, step 10's timings, which can't
    reproduce byte for byte. Then `v run tools/benchmarks.vsh` writes
    `docs/benchmarks.md`.
  - `tools/extract_if97.py R7-97.pdf SR3-03.pdf SR5-05.pdf` (PyMuPDF
    1.24.14) writes `water/if97/coefficients*.v` from IAPWS's release PDFs,
    which stay out of the repository; then `v fmt -w` them. Run it only to
    check or redo the transcription: the files are its output, as checked.
  - `gen_m3_fixtures.py` takes the names of the fixtures to write (`thermo`,
    `stability`, `pip`, `rachford_rice`); with none, it writes them all.
  - A generator must reproduce its committed output byte for byte on the
    same day. Check with `cmp` after regenerating.
- **DWSIM 9.0.5** is installed at `%LOCALAPPDATA%\DWSIM`. pythonnet 3.0.5
  stops at Python 3.13, so pass the same interpreter as above.
  - `tools/gen_m3_dwsim.py` writes `testdata/m3_dwsim.json`. It imports
    `gen_m3_fixtures.py` for its `thermo` helpers.
  - `tools/gen_m4_dwsim.py` writes `testdata/m4_dwsim.json`, DWSIM's
    constants as component files, and `reference_dwsim.toml`. It imports
    `gen_m3_dwsim.py`, and needs V on the PATH: it solves that file with
    reflux for the open loop. About five minutes.
  - `tools/gen_m9_dwsim.py [part ...]` writes `testdata/m9_dwsim.json`: the
    parts conversion, equilibrium, gibbs, cstr, pfr and ideal_gas. A part
    not named keeps what the file has; a run that fails writes nothing. pfr
    takes about a quarter of an hour, at dV = 0.001. It imports step 0's
    spike, `spike/dwsim_pfr/run.py`. `REFLUX_TUTORIAL2` names the owner's
    Tutorial 2 file, to check the script's case against it.
  - thermo 0.6.1's flash at a vapour fraction strictly between 0 and 1
    fails on a mixture (`dew` unbound): bisect on its PT flash instead.
  - `spike/dwsim/run.py` is M3 step 0's spike, and
    `spike/dwsim_flowsheet/run.py` M4 step 0's: DWSIM's flowsheets and
    recycle (`docs/spike-dwsim-flowsheet.md`).
  - DWSIM's PR uses rounded constants, saturates pure components on
    ChemSep's vapour-pressure correlation, and has a wrong ideal mixing
    entropy. "What step 5 found" in `docs/spike-dwsim.md` has the details.
    Compare with it by state, and explain each difference.
  - In DWSIM's flowsheets, a failed unit is silent: check every unit's
    `Calculated` and `ErrorMessage`. Its recycle needs its tolerances scaled
    by 1e-8 to close the balances, its pump uses a liquid density that
    isn't PR's, and its flash drum is wrong with two inlets.
  - Its separator's `EnergyImb` is what goes in less what comes out, and its
    wrong mixing entropy moves a two-phase isentropic outlet: 0.42 K in
    step 8's expander. "What step 8 found" in
    `docs/spike-dwsim-flowsheet.md` has the rest.
  - `spike/dwsim_pfr/run.py` is M9 step 0's spike: DWSIM's PFR on the
    owner's reforming case (`docs/spike-dwsim-pfr.md`). `REFLUX_TUTORIAL`
    names the owner's file, to check the script's case against it; the file
    stays out of the repository.
  - DWSIM's PFR holds each section's volumetric flow, Z and pressure at the
    section's start, so its error is first order in dV, and its integrators
    run at 1e-3 relative, which nothing sets. A rate names a reaction's
    species R1, R2, P1, … by their order in its compound list.
  - Ergun's drop needs the gas viscosity. DWSIM's is the mole-fraction mean
    of ChemSep's eq. 102 with Jossi–Stiel–Thodos' correction, 30% below
    Wilke's in hydrogen-rich gas. A PFR's failure is silent too, and a
    reaction above its Tmax is skipped with no error.
  - DWSIM's conversion reactor:
    - It gives every compound a conversion, (in − out)/in, so a product's
      is −∞. A sensitivity study of it plots −∞ with no message.
    - Its conversion is a function of T alone, and isn't checked against
      equilibrium.
    - Its "PR78" package's departures aren't its PR package's with PR78's
      α: they differ from `thermo`'s by 3e-5.
    - Its adiabatic conversion reactor evaluates the expression at the
      inlet's T, not the outlet's.
  - `spike/dwsim_study/run.py` is M8 step 0's spike
    (`docs/spike-dwsim-study.md`): DWSIM's sensitivity analysis, Adjust,
    optimiser and regression. None has an API; each is a form of the
    classic app, run by its button's handler.
    - `DWSIM.Automation.Automation` (the first automation class) makes the
      app's main window unshown. A flowsheet is built with `Automation3`,
      saved, and loaded there as a `FormFlowsheet`; a bare one can't take
      objects. The forms read `My.Application.ActiveSimulation`.
    - A study's failed point records the last point's numbers, with no
      message. After a study, the restored flowsheet doesn't solve again
      here: use one loaded fresh.
    - The flowsheet solver doesn't run an Adjust unless it's marked
      simultaneous; its control panel does. The simultaneous solver at a
      tolerance it can't reach leaves the flowsheet at its derivative's
      last point, 1% off in the varied value.
    - Of the optimiser's 11 methods only `AL_BRENT_B`, bounded Brent's,
      finds an optimum. An objective expression over auxiliary variables
      throws in the nine that aren't Brent's (an index of −1 in its
      `FunctionValue`), unbounded Brent's sets the variable to 0, and
      the DotNumerics methods and IPOPT minimise when asked to maximise.
    - The regression ran and fitted no kᵢⱼ, by any of its four methods.
    - "Background parallel" solves serially.
    - The script's process can take a minute to exit after its last line.
    - Its PFR isn't smooth in its volume (M8 step 13): its outflow matches
      an exact answer to 5e-11 at most volumes and is off by up to 1e-4 at
      others, its integrator's unset 1e-3 tolerance, enough to move a flat
      optimum by 1.5%. A PFR solve of three species at dV 0.01 takes 16 s.
    - A compressor's work differs from reflux's by 1e-5 to 5e-5 on DWSIM's
      own constants, scattered with T: its ideal-gas H and S, integrated
      approximately (M3).
  - DWSIM's activity packages (M6 step 12, `docs/m6-plan.md`):
    - Each is at its defaults an ideal-gas vapour, Poynting on, ChemSep's
      liquid densities and Lee–Kesler enthalpies; its vapour pressures carry
      the full ChemSep's constants, 5e-5 above ChemSep LITE's.
    - Its Poynting term is exp(V·|P − Psat|/RT): below Psat it raises φ.
    - Its UNIQUAC takes each ordered pair's τ's from that pair's own record,
      and holds two inconsistent records for methanol–water.
    - Its Wilson is (Vⱼ/Vᵢ)·exp(−A/(1.9872·T)), V from the liquid's density
      at 298.15 K.
    - Its UNIFAC-LL is its VLE UNIFAC (VLE groups, VLE a for most pairs),
      each group fraction divided by the Q of the compound's Dortmund group
      numbers read in the VLE table.
    - Its three-phase flash: NestedLoops (the default) is right; InsideOut
      adds a vapour where there's none; GibbsMinimization misses splits.
      Its three-phase numbers move from run to run.
  - DWSIM's other reactors (M9 step 10):
    - The equilibrium reactor puts fugacities over 1 atm against K at
      1 bar. Its default tolerances, 0.01 and 1e-3, leave its answer up to
      3e-2 off: set both to 1e-10.
    - The Gibbs reactor needs `CreateElementMatrix()` after the feed is
      calculated. Its first solve from a feed lacking some products
      refuses with "invalid initial estimates", and a second runs. Its
      default IPOPT answer can leave the elements unbalanced by 7e-4, and
      its two solvers agree only to 1e-3.
    - The CSTR reacts nothing in the vapour, kinetic or catalytic, and
      reports success. Its `IsothermalTemperature` defaults to 0 K, and it
      needs a `Headspace`.
    - DWSIM's copy of ChemSep gives water's ideal-gas Cp as eq. 16 with its
      own constants, where reflux's file has eq. 107: 72 J/mol apart at
      1000 K.

## Architecture

The modules depend in one direction: `numeric` ← `tdep` ← `thermo` ←
`components`, and `flash` uses `thermo` and `numeric` but not
`components`. `quantity`, units of measure (M4) and M9's dimension algebra
(compound units, exact rational exponents of kg, m, s, K, mol, A), depends
on nothing, and so do `release`, the release gate's logic, and `chem`
(M9), the elements, formulas and SMILES strings. `quantity` also holds
reflux's one R, `quantity.molar_gas_constant`, which `thermo`, `activity`,
`kinetics`, `expr` and `fit` take as theirs. `expr` (M9), the
expression language, uses `quantity`, and has M8's dotted names for a
flowsheet's variables; `kinetics` (M9), reactions, uses
`chem`, `expr`, `quantity` and `thermo`, and its tests `components`.
`flowsheet` reads a file's reactions into `kinetics` (step 6).
`equilibrium` (M9), chemical equilibrium and the rates' consistency with
it, uses `kinetics`, `flash`, `thermo` and `numeric`, and its tests
`components` and `flowsheet`; `main` uses it for `reflux kinetics check`.
`reactors` (M9), the reactors as flowsheet units, uses `flowsheet`,
`equilibrium`, `kinetics`, `components` (for the viscosity) and the rest;
`main` gives `flowsheet.load` its makers.
`water` (M5a) is water and steam: IAPWS-95 as a `thermo.Model` and water's
transport, over `thermo`, `numeric` and its submodule `water.if97` (IF97 as
steam-table functions, which uses only `math`); its tests use `flash`, and
`flash`'s step-back test uses it, so `water` itself can't import `flash`.
`components` builds IAPWS-95's fluid through it (`Database.iapws95`), and
`main`'s `reflux water` uses it, `water.if97` and `flash`.
`activity` (M6), liquid activity models, uses `numeric` and `expr`; `components`
reads ChemSep's banks and UNIFAC tables into it and builds its models.
`transport` (M5b), a phase's viscosity, thermal conductivity and surface
tension, uses `tdep`, `thermo` and `quantity`, and knows nothing of
component files: `components` builds its compounds (`transport_pure`,
`mixture`, `Setup.transport`, `wilke`), `water` gives water's own through
its `Model` (`water.Transport95`, `water.TransportIf97`), `reactors`' PFR
takes its `Wilke` for Ergun's drop, `flowsheet` keeps a file's model for
its report (`Loaded.transport`, report format 5) and `capi` serves it as
`reflux_transport`. Its tests are external, so they import `components`.
`components` uses `chem` for a compound's formula. `study` (M8) runs a file's studies
and optimisations over `flowsheet`, the latter on `numeric.sqp`, and `main`
uses it for `reflux study` and `reflux optimize`. `fit` (M8) fits a file's
constants to data over `flowsheet` (which reads `[fit.NAME]`), `kinetics`,
`flash`, `thermo` and `numeric.least_squares`; `main` uses it for
`reflux fit`. `capi` (M10a) is the C API
over the rest, built into a shared library for the Python package in
`python/`.
`sweep`, the flash robustness sweep, is the one module that uses both
`components` and `flash`. `flowsheet` (M4) uses `flash`, `thermo` and
`numeric`, and its tests use `components`.

- **`thermo` is built on one contract.** `Model` (thermo/model.v) supplies
  the reduced residual Helmholtz energy F(T, V, n) = Aʳ/RT and its
  derivatives, with V the total volume and n the mole numbers.
  - Everything else is derived generically from F: the volume roots, P,
    ln φ and its derivatives, H, S, Cp, PIP and Psat.
  - Cubics (`Cubic`, with per-component `Alpha` functions) and the Péneloux
    wrapper `Translated` implement it. CPA and SAFT will too (M12), with no
    change to the flash.
  - The volume solver's tier is `volume_function(t, n)`: F at fixed T and n
    as a function of V, so that the parts in T and n (a cubic's A and B) are
    computed once per root, not on every Newton step.
  - `Cubic.work()` (M8) is `solve_tp` for many states at no allocation: one
    volume function, boxed once, moved to each composition in place, with
    `solve_tp`'s bits. One per thread; the reactors' `Source` uses it.
  - `Fluid` adds the ideal-gas part, for total properties. Its `model` is
    a `PhaseModel` (M6's D3): never downcast it with `as`, which compiles
    and returns zeros. Ask `helmholtz()` for the `Model`, then test that.
  - A model may say what it knows of its isotherms' two branches through
    the optional `Branches` (M5a's D7): a split volume between them, and
    an isotherm with no state. The volume solvers keep to
    them, and `phase_label` labels by the split. IAPWS-95 needs it (its
    dome has spurious roots 1e19 RT below the real phase's); a cubic
    doesn't have it, and its paths are unchanged.
- **`water`** (M5a) holds IAPWS-95: `iapws95()`, the `Model`, F = n·ψ with
  ψ = (R₉₅/R)φʳ + (R₉₅/R − 1) ln δ (D4: its own gas constant, so its tables
  match); its ideal part as a `PlanckEinstein` Cp; `props95`, the release's
  Table 3 on IAPWS's reference, and `reference_offsets`, the constants
  between it and reflux's (D6). Below 251.165 K it has no state.
- **`water.if97`** (M5a) is IAPWS-IF97 with SR3-03 and SR5-05, per kg on
  IAPWS's reference: each basic and backward equation as published, and
  `state_pt`, `state_ph`, `state_ps`, `state_tx`, `state_px`, which choose
  the region, refuse what's outside the range (`If97Error`), and finish
  each backward answer by Newton on its basic equation (D9). Region 3's
  densities keep to their side of ρc below pc, as IAPWS-95's do.
- **`water.If97Phases`** (M5b) is IF97 as a `thermo.PhaseModel`: the
  region map's state, each side's equation across the line, residuals
  against region 2's ideal gas (a `thermo.Polynomial` Cp), so a stream's H
  and S are IF97's plus constants. Its saturation line is region 4's,
  through the optional `thermo.Saturates`, which the flash's paths for one
  compound take (`Path.own_line`); its range through `thermo.NoState`, for
  the searches' step back. It has no F(T, V): TV and UV refuse it.
- **`flash` works for any `thermo.PhaseModel`** (M6's D3): a phase at
  (T, P, n) with ln φ, their derivatives, a label and residual properties.
  `thermo.by_roots(model)` is any `Model` as one, and `tv` alone takes a
  `Model` itself.
  - `stability` is Michelsen's tangent-plane test. Its trials are Wilson's
    vapour- and liquid-like ones plus a nearly pure phase of each component.
    Substitution with DEM is followed by a second-order finish in α = 2√W.
  - `pt` tests stability, then splits: substitution with Rachford–Rice and
    DEM (an extrapolation kept only if it lowers the residual), then Newton
    on G through `numeric.minimize`, then each phase re-verified. Its feed's test stops at the first proof of instability
    once Wilson's two trials have run, and runs the rest only if no split
    from those stands.
  - `tvf` and `pvf` (`vf.v`) find P or T for a given β; bubble and dew
    points are β = 0 and 1. Newton runs on ln K and ln P (or ln T)
    from Wilson's start, with a PT-split start as fallback. Two answers
    are possible near critical points, and `estimate` picks between them.
  - `ph` and `ps` (`ph.v`) take a `thermo.Fluid` and find T from a molar H
    or S, by Michelsen's Q-function. T is bracketed by PT flashes on the
    sign of H(T) − H_spec, narrowed by Brent, and finished by Newton in
    (v, ln T). Only the answer's PT flash is verified. A pure component
    inside its jump at Tsat gets β from the energy balance.
  - `tv` and `uv` (`tv.v`) find P (and T) from a molar V (and U). The finish
    minimises A over the vapour's mole numbers and volume, in F(T, V, n)
    alone, with no volume roots (D5); UV adds ln T and the energy balance to
    that Newton.
  - In a split, the vapour is the phase of lower *mass* density (D6, as
    amended in step 6), not the larger molar volume.
  - Results are `FlashResult`; failures are `FlashError` with a
    `Failure` kind.
- **`numeric`** holds the solvers the flash rests on: LU, Cholesky, `newton`
  for systems (line search or dogleg), `broyden` for systems with no Jacobian,
  and `minimize` (a Levenberg–Marquardt-shifted Newton that knows f's rounding
  through `fscale`). Also `rachford_rice` (accurate for K from 1e-30 to 1e25),
  QAG integration and bisection, and M9's integrators, `dopri5` and `rodas4`
  (`Ode`, `StiffOde`), which return an answer only once every step has met
  the tolerance. M8's `qp` (Goldfarb and Idnani) and `sqp` (Han and Powell,
  `Constrained`, `OptimizeError`) return an optimum only once its first-order
  conditions hold, and `least_squares` (`Residuals`, `LeastSquaresError`) a
  minimum only once its gradient test does. They take what they solve as an interface
  (`Function`, `FallibleFunction`, `System`, `Equations`, `Objective`), never
  a closure: V 0.5.2 never frees a closure, and one per flash leaked kilobytes.
- **`components`** runs the data pipeline: ChemSep XML, a strict importer,
  one TOML file per compound in `data/components/`, then the database API.
  The API yields `thermo.Component`, `IdealGas` and kᵢⱼ matrices from the
  vendored `pr.ipd`/`srk.ipd`, taking the first entry listed for a pair.
  - A `Database` opens on a data root, `components/` and `vendor/chemsep/`,
    with `components.open(root, extra)`. The default root is `$REFLUX_DATA`,
    or else the source tree's `data/` (`@VMODROOT`). A fluid's kᵢⱼ and
    Mathias–Copeman α come from its database's root (M10a's D9), never from
    the source tree behind its back.
- **`flowsheet`** is M4's: streams, units, the graph, the solver, the file
  and the reports, and M8's specs.
  - A `Stream` is its molar flows plus the flash that gave its state, so
    it's never half specified. Its H, S and V are its phases'.
  - A `Unit` is `calculate(fluid, inlets) !Outcome` and holds nothing about
    a solve. An energy rate above zero goes into the unit's material.
  - A `Flowsheet` is units with a stream named on each port, and `plan`
    checks it and returns blocks to solve in order. Tarjan gives the blocks
    and Johnson the cycles, and `tear.v` picks the fewest tears exactly, by
    branch and bound, then the fewest second breaks.
  - A file is read by `load` through a `Section` per table, and a `Maker` per
    unit type. Later unit families add makers, so `flowsheet` needn't import
    them (D1).
  - `solve` iterates each block's tears, whose variables are flows, P and
    molar H, by direct substitution, Wegstein or Broyden. A `Solution` is only
    returned with every balance checked: each unit to 1e-10 on the inlets it
    used, and the whole to the gap the tears left.
  - A `Section` notes every number its maker reads (`Readings`), so a key
    such as `unit.cooler.outlet_T` can be varied: `Loaded.with` remakes the
    feed or unit by its maker (M8's D4). `Variables` names each stream's and
    unit's variables for expressions, `gas.F`, `cooler.duty` (D3).
  - A spec (`[spec.NAME]`, format 3) is a link in planning from what it
    reads to what it varies, so both land in one block; it's met by Brent
    alone, or by Broyden with the tears, and checked again from the
    solution (D5). `Loaded.solve` solves on the plan made at load.
  - `Loaded.solve_from(Start)` starts from a neighbour's solution (M8's
    D6); a unit implementing `Warm` starts from its last answer.
  - `[study.NAME]` and `[optimize.NAME]` are read into `StudyPlan` and
    `OptimizePlan` (`studies.v`, `optimize.v`); the `study` module runs
    them.
  - A `Report` is a solution, or a failure, in SI: `Loaded.report`,
    `Loaded.failure` or `load_failure`. Its JSON and CSV come from that one
    struct, and its stream table rounds by `math.round`, never printf.
    `main.v` is the command line over it, tested in process by
    `main_test.v`.
- **`capi`** is M10a's, in progress: the functions `capi/reflux.h` declares,
  and nothing else is exported. Each returns 0 or a negative code with its
  reason in the caller's buffer, takes `i32`, `i64`, `f64` and pointers only,
  and registers its thread with the collector for the call. Handles are ids
  into a registry behind a lock, never reused. The library starts on its
  first call, which mustn't race, so the Python package makes it at import.
  - Text comes out into the caller's buffer, with the room it takes in
    `*needed`, or `REFLUX_TOO_SMALL`. Its formats are components' own: a
    compound is its component file, and a fluid's setup `Setup.to_toml`.
  - Everything a caller can get wrong is checked in `capi` first
    (`REFLUX_BAD_ARGUMENT`), so that what `components` refuses after is the
    data's doing (`REFLUX_DATA_ERROR`). A lookup is `REFLUX_NOT_FOUND`.
  - The header has no `const`: V writes `char *`, and gcc refuses a header
    that disagrees.
  - Many points per call: each input's count is 1 or n, outputs are
    property-major (`out[k*n + i]`), and each point has a status
    (`REFLUX_POINT_*`), NaN when it has no answer. A point with a value that
    isn't finite has no answer (`NOT_FINITE`). `err` holds the first failed
    point's reason, so one point's reason is the call for it alone.
  - Every flash is `reflux_flash(fluid, spec, …)` with one layout (D6): the
    whole's values, then each phase's in the place of its label, liquid (0)
    then vapour (1). A single phase takes its label's place, the other NaN.
- **`equilibrium`** is M9's (D9, D10): one phase, ideal gas or a cubic.
  - `gibbs` minimises G over the amounts the feed's elements allow, and
    `over_reactions` over those chosen reactions allow, both by RAND
    (`rand.v`) in log space, φ updated outside it for a cubic. A row with
    every coefficient of one sign and none of it in the feed drops its
    species first.
  - An answer is returned only once its balance, its stationarity (or each
    reaction's ln Q = ln K) and, for a cubic, its phase's stability hold.
  - `consistency` judges each reversible rate by where it stops along its
    reaction from the equilibrium a state reaches: ΔG/RT there, zero for a
    rate consistent with the Gibbs energies.
- **`reactors`** is M9's (D10): Gibbs, equilibrium over chosen reactions,
  conversion, CSTR and PFR, each a `flowsheet.Unit` with a maker (D1).
  - A maker reads its table through a `Section` whose `context` holds the
    file's format, species, component records and reactions (format 2).
  - A reactor's `Outcome` gives its elements as `conserved`, and the solver
    checks it on them; a PFR's gives its `profile`, and a conversion
    reactor's its `warnings` (an outlet past equilibrium).
  - Rates become each species' source through `Source`, with exact
    derivatives; the PFR and CSTR use them in `rodas4` and `newton`. The
    CSTR continues in its volume where Newton from the inlet fails.
  - What depends on T alone is kept until T moves (M8's D11), and is the
    same number kept or not: each reaction's K, in the rates' `Point`, and
    what Wilke's rule needs (`WilkeAt`), in the PFR.
  - A conversion may be an expression in T, bound in the kinetics' rate
    context with the file's definitions, evaluated at the outlet's T.
- **`activity`** is M6's (D2): a liquid's activity coefficients at T and
  mole fractions, knowing nothing of phases.
  - `Model` is immutable and shared; each thread's `Work` keeps what
    depends on T alone (τ, UNIFAC's ψ and pure-group terms) until T moves,
    and nothing once a T is refused.
  - `eval` returns what was asked for only if it's finite, and refuses a
    composition that doesn't sum to 1 and a parameter that isn't a normal
    double, so a caller normalises x first.
  - Each formula is written once, generic over f64 and `numeric.Jet`, and
    `eval` takes its T-derivatives from the parameters as jets and its
    composition derivatives one xⱼ at a time, n·∂ln γᵢ/∂nⱼ = ∂ᵢⱼ − Σₖ xₖ∂ᵢₖ.
  - A parameter is an expression in T (`Param`), its τ″ from `Expr.jet`.
    ChemSep's banks are written with DECHEMA's R, `A K / (1.98721*T)`,
    confirmed for NRTL's bank and inferred for the others (step 12).
  - `components` reads the banks and UNIFAC tables (`load_pair_bank`,
    `load_unifac`), repairing each self-contradiction and correcting the
    published typos, all reported; `Database.nrtl`, `uniquac`, `wilson` and
    `unifac` build the models.
- **Units are SI inside** (K, Pa, m³, mol, J), with conversion only at the
  edges.
- **No panics reach a caller.** Every failure is a typed error
  (`FlashError`, `NewtonError`, `InputError`, `LookupError`,
  `QuantityError`, `FlowsheetError`), because the Python API (M10a) can't
  catch panics.
- **A solver returns a result only when it has shown the result to be
  right**, meaning its tolerance is met and the split verified. Anything
  else is an error, never a plausible-looking answer.

## Testing practice

- **References and tolerances:**
  - References come from independent sources (`thermo`, `chemicals`,
    CoolProp, mpmath at high precision, brute-force scans). They are stored
    as JSON in `testdata/` with provenance. Tests never call Python.
  - Measure first, then assert a tolerance justified by conditioning or
    rounding, and say why in the test's comment.
  - Model-accuracy deviations and iteration budgets are recorded as
    numbers that can't grow.
- **Check the claim, not the code path.** For example, prove an "unstable"
  verdict with a TPD computed from scratch, and check every split's
  fugacities, mass balance and phase stability independently.
- **Planted-mistake checks for every step:**
  - Mutate one snippet at a time in a scratch copy of the modules (never in
    the repo while a check is running), run the relevant tests, restore,
    and confirm the files match with `cmp`.
  - Each survivor gets a new test, gets its code deleted, or is explained
    in the plan.
  - A mutation that doesn't compile isn't caught: no test saw it. Under `-W`
    an import left unused or a `mut` never changed is an error, so write the
    mutation to keep both used (`if !utf8.validate_str(s) && false`). The
    scripts report it as NO BUILD.
  - A timeout must kill the whole process tree: killing `v` leaves the test
    program it started holding the pipes, and the script waits forever.
- **A NaN passes `x > bound`:** write a bound as `assert x <= bound`, or
  `if !(x <= bound)`. `math.max` carries a NaN through.
- **Keep sensitive states exact.** Where a test state is sensitive to an
  ulp, store the exact double: build it from integers, for instance
  `(915 + k) * 1e4`.

## V notes that bite here

`docs/v-for-go.md` has the full set. The ones that most often break a build
or a test:

- `-W` turns warnings into errors, unused variables among them.
- **Shadowing is uneven:** a local named like a const is always an error
  ("duplicate of a const name"). One named like a function is only a
  notice, but the checks are kept free of notices, so rename it.
- **Results:** `return f()` forwards a result, but inside a `match` or `if`
  arm the call needs `!`. An Option can't be propagated with `?` from a
  function that returns `!T`.
- **Arrays:** a `mut` variable can't alias an immutable array, so use
  `.clone()`. Two `mut` variables share, as Go slices do. There's no `a + b`;
  `a << b` appends every element of b.
- **`v vet` doesn't type-check**, so compile what you change, `.vsh`
  scripts included.
- **On Windows, `import runtime` defeats V's default compiler** (tcc can't
  parse `psapi.h`). Only the fallback compiler builds such a program.
- **V reports its own failed builds** to bugs.vlang.io when its fallback
  compiler succeeds: build options and C errors, with local paths. Set
  `V_C_ERROR_BUG_REPORT_DISABLED=1` when probing V's edges.
- **Structs:** a struct field of function type needs `@[required]` or a
  default. Writing `= false` on a `bool` field is an error, because fields
  are zeroed already.
- **Strings:** `'${x}'` and `"${x}"` both interpolate, so V source built as
  text needs raw strings, `r'...'`.
- **Writing files from the shell:** bash heredocs mangle backslashes, so
  write such files with the Write tool or a Python script. On Windows,
  Python's `write_text` writes CRLF: use `write_bytes`, or
  `newline='\n'`, and keep a script's own `\r` out of non-raw strings.
- **Formatting:** printf-style formats round ties differently under tcc
  and gcc, so round with `math.round` before formatting anything compared
  byte for byte.
- **Shared libraries** need gcc: tcc's Windows DLLs export nothing.
- **Closures are never freed**, so nothing that runs per flash may make one:
  pass a struct with a method, as `numeric`'s interfaces ask.
- **`gc_heap_usage().total_bytes` is batched under gcc:** small allocations
  are counted only when a free list refills, so one call's change can be 0.
  A test may hold an upper bound or "nothing allocated" (tcc's build
  catches it), never a count of work. Count the work itself.
- **A panic can't be caught.** There's no `recover`, so what a library might
  panic on must be kept from it: `toml.parse_text` on a `\u` escape, for one.
- **Module `main`'s functions keep their C names,** unprefixed: one named
  like a C library function (`finite`, in MinGW's `math.h`) builds under
  tcc and fails under Windows' gcc. Build `main` with gcc when you add one.
- **Module names:** not like a vlib module's. `benchmark` broke external
  tests, because V's test runner imports vlib's.
- **Leftovers:** V can leave `.*.v3cc.*/` build folders next to a binary.
  They're git-ignored; delete them.
- **Smart App Control** can refuse V's programs, and V reports it as `failed
  CreateProcess ...: The operation completed successfully`. A refused run
  proves nothing, so a mutation loop must retry it, not count it as caught.
  It's the owner's security setting: don't work around it. See
  `docs/v-for-go.md`, §1.
- **`json2`** leaves a `none` field out unless it has `@[json_null]`, and
  writes a map's keys sorted, so anything ordered goes in an array. In
  0.5.2 `Any.arr()` is deprecated, an error under `-W`: use `as_array()`.
  It knows a type by its short name: a struct named like another module's
  type that it decodes breaks the decode (M5b's `transport.MethodUsed`).
- **A `!` in a test function that doesn't return a Result** can turn a typed
  error into a `V panic` with no FAIL line, so a panic in a test's output
  isn't always a panic in the code.
- **Inside `or {}`, `err` is the error,** even where a parameter is called
  `err`: name a C buffer parameter `buf`. And `T{...err, …}` on a smart-cast
  `err` copies nothing: cast with `as` first, or name the fields. A
  smart-cast interface value spread straight into an interface variable
  (`m = T{...x, …}`) copies garbage: into a local first, or name them.
- **A fixed array of structs inside another struct's zero value** skips the
  element type's field defaults. Named struct fields keep theirs.
- **`v fmt -w` rewrites a file that doesn't parse,** and can drop lines.
  Compile first; mind an apostrophe in a single-quoted string.
- **`@[has_globals]`** allows `__global` in the file it heads, not its whole
  module.
- **A V shared library starts on its first call,** from whichever thread,
  behind a flag that isn't thread-safe: one call must come first, from one
  thread. A Python thread must be registered with the collector on Windows
  as on Linux (`capi/gc_threads.h`). The collector keeps the first call's
  thread registered for good, and on Linux it aborts once that thread has
  ended, so the library unregisters it unless it's the process's first.
