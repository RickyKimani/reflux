// Internal test (module water): IAPWS-95 against its release's Tables 6 and
// 7, through the model's own Table 3 and through thermo's generic
// properties of the Model, and D4 to D8's claims.
module water

import math
import thermo

// half_unit is half a unit of the last figure of a value as the release
// prints it, "0.204797733e1": what reproducing it to its printed digits
// allows.
fn half_unit(printed string) f64 {
	mant := printed.all_before('e')
	exp := if printed.contains('e') { printed.all_after('e').int() } else { 0 }
	decimals := mant.all_after('.').len
	return 0.5 * math.pow(10.0, f64(exp - decimals))
}

fn printed_ok(got f64, printed string) bool {
	return math.abs(got - printed.f64()) <= half_unit(printed)
}

// Printed is a value computed and the release's printing of it. (A V array
// literal holds one type, unlike a Python list, so a row is a struct.)
struct Printed {
	got  f64
	want string
}

// Table 6: φ° and φʳ with their derivatives at 500 K and 838.025 kg/m³.
fn test_table_6() {
	tau := critical_t / 500.0
	delta := 838.025 / critical_rho_mass
	o := phi_ideal(tau, delta)
	r := phi_residual(tau, delta)
	rows := [
		Printed{o.f, '0.204797733e1'},
		Printed{o.d, '0.384236747'},
		Printed{o.dd, '-0.147637878'},
		Printed{o.t, '0.904611106e1'},
		Printed{o.tt, '-0.193249185e1'},
		Printed{r.f, '-0.342693206e1'},
		Printed{r.d, '-0.364366650'},
		Printed{r.dd, '0.856063701'},
		Printed{r.t, '-0.581403435e1'},
		Printed{r.tt, '-0.223440737e1'},
		Printed{r.dt, '-0.112176915e1'},
	]
	for i, row in rows {
		assert printed_ok(row.got, row.want), 'Table 6 row ${i}: ${row.got} against ${row.want}'
	}
	assert o.dt == 0.0
}

struct Row7 {
	t   f64
	rho string // kg/m³
	p   string // MPa
	cv  string // kJ/(kg·K)
	w   string // m/s
	s   string // kJ/(kg·K)
}

// Table 7: p, cv, w and s at eleven (T, ρ).
const table7 = [
	Row7{300, '0.9965560e3', '0.992418352e-1', '0.413018112e1', '0.150151914e4', '0.393062643'},
	Row7{300, '0.1005308e4', '0.200022515e2', '0.406798347e1', '0.153492501e4', '0.387405401'},
	Row7{300, '0.1188202e4', '0.700004704e3', '0.346135580e1', '0.244357992e4', '0.132609616'},
	Row7{500, '0.4350000', '0.999679423e-1', '0.150817541e1', '0.548314253e3', '0.794488271e1'},
	Row7{500, '0.4532000e1', '0.999938125', '0.166991025e1', '0.535739001e3', '0.682502725e1'},
	Row7{500, '0.8380250e3', '0.100003858e2', '0.322106219e1', '0.127128441e4', '0.256690919e1'},
	Row7{500, '0.1084564e4', '0.700000405e3', '0.307437693e1', '0.241200877e4', '0.203237509e1'},
	Row7{647, '0.3580000e3', '0.220384756e2', '0.618315728e1', '0.252145078e3', '0.432092307e1'},
	Row7{900, '0.2410000', '0.100062559', '0.175890657e1', '0.724027147e3', '0.916653194e1'},
	Row7{900, '0.5261500e2', '0.200000690e2', '0.193510526e1', '0.698445674e3', '0.659070225e1'},
	Row7{900, '0.8707690e3', '0.700000006e3', '0.266422350e1', '0.201933608e4', '0.417223802e1'},
]

// A per-mole value in J… as the table's per-kg one in kJ….
fn per_kg_k(x f64) f64 {
	return x / molar_mass * 1e-3
}

fn test_table_7_by_table_3() {
	for row in table7 {
		rho := row.rho.f64() / molar_mass
		pr := props95(row.t, rho)
		assert printed_ok(pr.p * 1e-6, row.p), 'p at ${row.t} K, ${row.rho}: ${pr.p * 1e-6} against ${row.p}'
		assert printed_ok(per_kg_k(pr.cv), row.cv), 'cv at ${row.t} K, ${row.rho}: ${per_kg_k(pr.cv)}'
		assert printed_ok(pr.w, row.w), 'w at ${row.t} K, ${row.rho}: ${pr.w} against ${row.w}'
		assert printed_ok(per_kg_k(pr.s), row.s), 's at ${row.t} K, ${row.rho}: ${per_kg_k(pr.s)}'
	}
}

fn test_table_7_through_thermo() {
	// The same table through the Model and thermo's generic properties
	// (D3, D4, D5): P from F, cv, w and s from Fluid.props, s on IAPWS's
	// reference by D6's constant, taken at Table 6's state.
	fl := fluid(-241814.0, 188.724)!
	_, ds := reference_offsets(fl)
	for row in table7 {
		rho := row.rho.f64() / molar_mass
		st := state_at_density(fl.model, row.t, rho)
		pr := fl.props(st)
		assert printed_ok(st.p * 1e-6, row.p), 'p at ${row.t} K, ${row.rho}: ${st.p * 1e-6} against ${row.p}'
		assert printed_ok(per_kg_k(pr.cv), row.cv), 'cv at ${row.t} K, ${row.rho}: ${per_kg_k(pr.cv)}'
		assert printed_ok(pr.w, row.w), 'w at ${row.t} K, ${row.rho}: ${pr.w} against ${row.w}'
		assert printed_ok(per_kg_k(pr.s - ds), row.s), 's at ${row.t} K, ${row.rho}: ${per_kg_k(pr.s - ds)}'
	}
}

fn test_the_references_differ_by_constants() {
	// D4 and D5 make the total Helmholtz energy IAPWS-95's but for a term
	// linear in T, so h and s differ from Table 3's by the same constants at
	// every state, liquid, vapour or supercritical, and cp, cv, w and P not
	// at all. Measured: h and s to 2.9e-14 and 2.4e-14; p, cp, cv and w to
	// 1.4e-13 (the liquid's p at 0.1 MPa, a difference of terms 1400 times
	// its size) but for cp at 647 K and 358 kg/m³, 4.3e-12, where (∂p/∂ρ)_T,
	// which cp divides by, is 0.4% of ρRT and the two routes round it
	// differently.
	fl := fluid(-241814.0, 188.724)!
	dh0, ds0 := reference_offsets(fl)
	mut worst_h := 0.0
	mut worst_s := 0.0
	mut worst_x := 0.0
	for row in table7 {
		rho := row.rho.f64() / molar_mass
		st := state_at_density(fl.model, row.t, rho)
		pr := fl.props(st)
		ia := props95(row.t, rho)
		worst_h = math.max(worst_h, math.abs(pr.h - ia.h - dh0) / math.abs(ia.h))
		worst_s = math.max(worst_s, math.abs(pr.s - ia.s - ds0) / math.abs(ia.s))
		for pair in [[st.p, ia.p], [pr.cp, ia.cp], [pr.cv, ia.cv], [pr.w, ia.w]] {
			worst_x = math.max(worst_x, math.abs(pair[0] / pair[1] - 1.0))
		}
	}
	assert worst_h < 1e-13, 'h: ${worst_h}'
	assert worst_s < 1e-13, 's: ${worst_s}'
	assert worst_x < 2e-11, 'p, cp, cv, w: ${worst_x}'
}

fn test_the_ideal_part_is_iapws_95_s() {
	// reflux's ideal gas with ideal_cp (D5) against Eq. (5)'s: Cp is R + R₉₅
	// times −τ²φ°_ττ... (D4's R − R₉₅ in it), and ∫Cp dT is the change in
	// R₉₅T(1 + τφ°_τ) plus (R − R₉₅)ΔT.
	cp := ideal_cp()
	r := thermo.r_gas
	for t in [130.0, 273.16, 500.0, 1273.0, 5000.0] {
		o := phi_ideal(critical_t / t, 1.0)
		tau := critical_t / t
		want := r95 * (1.0 - tau * tau * o.tt) + (r - r95)
		assert math.tolerance(cp.cp(t), want, 1e-14), 'Cp at ${t} K: ${cp.cp(t)} against ${want}'
		h := fn (t f64) f64 {
			tau := critical_t / t
			return r95 * t * (1.0 + tau * phi_ideal(tau, 1.0).t)
		}
		dh := cp.int_cp(298.15, t)
		want_dh := h(t) - h(298.15) + (r - r95) * (t - 298.15)
		assert math.abs(dh - want_dh) <= 1e-9 * math.abs(h(t)), '∫Cp at ${t} K: ${dh} against ${want_dh}'
	}
}

fn test_the_critical_point() {
	// D8: at (Tc, ρc) itself φʳ and its first derivatives are finite, and so
	// are φʳ_δδ and φʳ_δτ; φʳ_ττ is −∞. Nothing is NaN.
	r := phi_residual(1.0, 1.0)
	for x in [r.f, r.d, r.dd, r.t, r.dt] {
		assert math.is_finite(x)
	}
	assert math.is_inf(r.tt, -1)
	// δ = 1 at any other τ: Table 5's 0/0, finite here.
	for tau in [0.5, 0.99, 0.999999, 1.000001, 1.5] {
		x := phi_residual(tau, 1.0)
		for y in [x.f, x.d, x.dd, x.t, x.tt, x.dt] {
			assert math.is_finite(y), 'τ = ${tau}, δ = 1'
		}
	}
}

fn test_the_critical_point_s_properties() {
	// At (Tc, ρc) P is Pc, and Cv and Cp are +∞, through Table 3 and through
	// thermo's Fluid.props alike. Table 3's w is 0 to rounding (7.9e-5 m/s);
	// Fluid.props' w is NaN there, its Cp/Cv being ∞/∞ (D8, as amended in
	// step 1): thermo's formula is the cubics', and every cubic's w would
	// move with it.
	fl := fluid(-241814.0, 188.724)!
	st := state_at_density(fl.model, critical_t, critical_rho)
	assert math.abs(st.p / critical_p - 1.0) < 1e-12, 'P at the critical point: ${st.p}'
	pr := fl.props(st)
	assert math.is_inf(pr.cv, 1) && math.is_inf(pr.cp, 1)
	assert math.is_nan(pr.w)
	ia := props95(critical_t, critical_rho)
	assert math.is_inf(ia.cv, 1) && math.is_inf(ia.cp, 1)
	assert ia.w < 1e-3, 'w at the critical point: ${ia.w}'
}
