module water

// IAPWS-95's coefficients, transcribed from IAPWS R6-95(2018), "Revised
// Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of
// Ordinary Water Substance for General and Scientific Use" (Wagner and Pruß,
// J. Phys. Chem. Ref. Data 31, 387, 2002), Tables 1 and 2. IAPWS allows the
// release's publication in whole or in part with attribution. A test holds
// every number to the `iapws` package's (testdata/m5_iapws95.json), and
// Tables 6 to 8 hold the equations.

// Table 1: the ideal-gas part, Eq. (5). n₁° to n₃°, then n₄° to n₈° with
// their γᵢ°.
const ideal_n1 = -8.3204464837497
const ideal_n2 = 6.6832105275932
const ideal_n3 = 3.00632
const ideal_n = [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
const ideal_gamma = [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]

// Table 2, i = 1 to 7: the polynomial terms nᵢ δ^dᵢ τ^tᵢ.
const poly_d = [1, 1, 1, 2, 2, 3, 4]
const poly_t = [-0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0]
const poly_n = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1, 0.31802509345418,
	-0.26145533859358, -0.78199751687981e-2, 0.88089493102134e-2]

// Table 2, i = 8 to 51: the exponential terms nᵢ δ^dᵢ τ^tᵢ e^(−δ^cᵢ).
const expo_c = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]
const expo_d = [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6,
	6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
const expo_t = [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10, 10,
	6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50]
const expo_n = [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4, -0.19232721156002,
	-0.25709043003438, 0.16074868486251, -0.40092828925807e-1, 0.39343422603254e-6,
	-0.75941377088144e-5, 0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
	0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9, -0.10793600908932,
	0.17611491008752e-1, 0.22132295167546, -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
	-0.31358700712549e-1, -0.74315929710341, 0.47807329915480, 0.20527940895948e-1, -0.13636435110343,
	0.14180634400617e-1, 0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
	-0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2, 0.15870308324157e-3,
	-0.16388568342530e-4, 0.43613615723811e-1, 0.34994005463765e-1, -0.76788197844621e-1,
	0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
	0.31777497330738, -0.11841182425981]

// Table 2, i = 52 to 54: the Gaussian terms
// nᵢ δ^dᵢ τ^tᵢ exp(−αᵢ(δ − εᵢ)² − βᵢ(τ − γᵢ)²).
const gauss_d = [3, 3, 3]
const gauss_t = [0, 1, 4]
const gauss_n = [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]
const gauss_alpha = [20.0, 20.0, 20.0]
const gauss_beta = [150.0, 150.0, 250.0]
const gauss_gamma = [1.21, 1.21, 1.25]
const gauss_epsilon = [1.0, 1.0, 1.0]

// Table 2, i = 55 and 56: the non-analytic terms nᵢ Δ^bᵢ δ ψ, with
// Δ = θ² + Bᵢ[(δ − 1)²]^aᵢ, θ = (1 − τ) + Aᵢ[(δ − 1)²]^(1/2βᵢ) and
// ψ = exp(−Cᵢ(δ − 1)² − Dᵢ(τ − 1)²). Both terms share a, B, A and β, which
// the evaluation uses: their powers of (δ − 1)² are taken once.
const nonan_n = [-0.14874640856724, 0.31806110878444]
const nonan_b = [0.85, 0.95]
const nonan_c = [28.0, 32.0]
const nonan_dd = [700.0, 800.0]
const nonan_a = 3.5
const nonan_bb = 0.2
const nonan_aa = 0.32
const nonan_beta = 0.3
