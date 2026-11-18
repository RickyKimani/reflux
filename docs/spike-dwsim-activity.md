# DWSIM 9's activity-coefficient packages, by script (M6 step 0)

What M6 needs to know before comparing with DWSIM's NRTL, UNIQUAC, Wilson and
UNIFAC packages (D17), and what they say about M6's own decisions. Two
scripts, run on 2026-11-18 with DWSIM 9.0.5 and `thermo` 0.6.1:
- [spike/dwsim_activity/run.py](../spike/dwsim_activity/run.py) drives
  DWSIM's packages on methanol–water and water–1-butanol;
- [spike/dwsim_activity/lle.py](../spike/dwsim_activity/lle.py) computes
  each model's own liquid–liquid split for water–1-butanol, by `thermo`'s
  NRTL and UNIQUAC with ChemSep's parameters, as the lower convex hull of the
  Gibbs energy of mixing on a grid of 1/20,000.

## Its settings

Every activity package has the same defaults, listed from its properties:

| Setting | Default |
|---|---|
| `VaporPhaseFugacityCalculationMode` | `Ideal`: the vapour is the ideal gas |
| `LiquidFugacity_UsePoyntingCorrectionFactor` | `True` |
| `LiquidDensityCalculationMode_Subcritical` | `Rackett_and_ExpData` |
| `EnthalpyEntropyCpCvCalculationMode` | `LeeKesler` |
| `FlashCalculationApproach` | `NestedLoops` |
| `ActivityCoefficientModels_IgnoreMissingInteractionParameters` | `False` |

## Its parameters are ChemSep's

DWSIM's NRTL and UNIQUAC packages hold ChemSep's banks, to the digit:
methanol–water's NRTL is A₁₂ = −189.0469, A₂₁ = 792.802 cal/mol and
α = 0.2999, and its UNIQUAC −337.1298 and 549.2958, each with its B and C
terms zero. Water–1-butanol's are ChemSep's too. Its UNIQUAC r and q are the
component files' (methanol 1.43 and 1.43, water 0.92 and 1.40).

## NRTL's form: τ = A/(1.98721·T)

DWSIM's liquid φ is γ·Pˢᵃᵗ·exp(Vᴸ(P − Pˢᵃᵗ)/RT)/P. With its Poynting term
taken out, on its own Rackett volumes (4.196e-4 for methanol and 5.610e-4
for water at 330 K and 1 atm), its γ at x = (0.3, 0.7) is NRTL's with
τ = A/RT for this R in calories:

| R, cal/(mol·K) | Gap in ln γ, methanol | Water |
|---|---|---|
| 1.98721 | −3.7e-10 | −5.0e-10 |
| 8.314462618/4.184 = 1.987204 | −5.3e-7 | −2.9e-7 |
| 8.314462618/4.1868 = 1.985875 | −1.2e-4 | −6.7e-5 |

So ChemSep's A's are in thermochemical calories, divided by DECHEMA's own
R, 1.98721 cal/(mol·K), with which DECHEMA regressed them. That R is
8.31449 J/(mol·K), not the SI's exact 8.314462618. The difference is 2.7e-6
of τ. **D7 is amended** (below): a bank's τ keeps the R its parameters were
fitted with, as `reflux kinetics check` already advises ("write R unless the
paper fitted with this number").

## Its enthalpies are Lee–Kesler's, not its activity model's

For a pure liquid at 330 K and 1 atm, DWSIM's H_L − H_ig is the same in all
six activity packages, and agrees neither with ChemSep's heat of
vaporization nor with Clapeyron on its own vapour pressure:

| | DWSIM's H_L − H_ig | ChemSep's ΔHvap (eq. 106) | Clapeyron, RT²·d ln Pˢᵃᵗ/dT |
|---|---|---|---|
| Methanol | −38,357.78 J/mol | 35,933.19 | 37,533.64 |
| Water | −45,414.76 | 42,652.86 | 42,801.09 |

Clapeyron's column takes the vapour as the ideal gas, as DWSIM does, and
neglects the liquid's volume. The difference is
`EnthalpyEntropyCpCvCalculationMode = LeeKesler`: the liquid's departure
comes from Lee and Kesler's corresponding states, whatever model gives its
fugacities. So a drum on DWSIM's NRTL takes its duty from one model and its
phase equilibrium from another. Methanol's latent heat there is 0.8 kJ/mol,
2%, above what its own vapour pressure implies. That's the inconsistency
D5 rules out. The other modes (step 12 lists them) decide how DWSIM is
compared.

## UNIQUAC for methanol–water: unexplained

At 330 K and x = (0.3, 0.7), DWSIM's UNIQUAC γ, Poynting taken out, is
(1.01196, 0.79960). DECHEMA's form, τᵢⱼ = exp(−Aᵢⱼ/RT), gives (1.34110,
1.09153) by `thermo` and by a hand-written UNIQUAC, to the last digit. That
agrees with NRTL's (1.37304, 1.09743), a fit to the same data. None of these
variants reproduces DWSIM's:

| Variant | γ |
|---|---|
| q in both terms (DECHEMA's) | (1.34110, 1.09153) |
| q′ in the residual term (0.96 and 1.00) | (1.32342, 1.08217) |
| A₁₂ and A₂₁ swapped | (1.36373, 1.06603) |
| A in kelvins, not cal/mol | (1.07410, 1.03957) |

Its PT flash of z = (0.5, 0.5) at 350 K and 1 atm finds one liquid, where
its NRTL finds 46% vapour. This pair is ChemSep's one methanol entry with no
DECHEMA page, "(Kojima+Kato)". It's for step 12.

## Its UNIQUAC splits water and butanol as the model does

| T | DWSIM's two liquids, water | The model's split (`lle.py`) |
|---|---|---|
| 298.15 K | 0.3808, 0.9901 | 0.3808, 0.9901 |
| 360 K | 0.4566, 0.9812 | 0.4567, 0.9812 |

The organic phase is 0.38 water by the model at 298 K, against about 0.51
measured: ChemSep's VLE-fitted UNIQUAC predicts the split, not well.

## Its NRTL misses the split its model has

ChemSep's NRTL for water–1-butanol (A₁₂ = 2633.6951, A₂₁ = 504.0381 cal/mol,
α = 0.4447) splits, 0.6008 and 0.9944 water at 298.15 K, and 0.5870 and
0.9837 at 360 K (`lle.py`).
- At 298.15 K, DWSIM's flash of z = (0.5, 0.5) returns one liquid at the feed,
  in 0 iterations.
- At 360 K it returns two "liquids" at 0.4998 and 0.5004 water, the trivial
  split.
- At 366 K both packages put the feed all vapour, which step 12 checks
  against each model's heterogeneous azeotrope.

So DWSIM's default three-phase flash can miss a liquid split, or return a
trivial one, with no message. Its other flash approaches are for step 12;
the algorithm classes weren't where the script looked for them.

## What it means for M6

- **D5 stands,** and DWSIM's default is the case against the alternative.
- **D7 is amended:** a bank's parameters keep DECHEMA's R.
- **D17's comparisons** are by state, as before: γ, a split's compositions
  and fugacities, and enthalpies under whichever of DWSIM's modes is closest
  to D5's, every difference explained.
