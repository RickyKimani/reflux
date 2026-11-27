// Internal test (module water): water's transport properties against their
// releases' verification tables (M5a's criterion 7): IAPWS R12-08's Tables 4
// and 5, R15-11's Tables 4, 5, 7, 8 and 9, and R1-76's Table 1, every value
// within half a unit of its last printed figure; and against the iapws
// package's values over a grid.
module water

import json2
import math
import os
import water.if97

fn half_unit_t(printed string) f64 {
	mant := printed.all_before('e')
	exp := if printed.contains('e') { printed.all_after('e').int() } else { 0 }
	decimals := mant.all_after('.').len
	return 0.5 * math.pow(10.0, f64(exp - decimals))
}

struct Row {
	got   f64
	want  string
	where string
}

fn check_rows(rows []Row) {
	mut bad := []string{}
	for r in rows {
		if !(math.abs(r.got - r.want.f64()) <= half_unit_t(r.want)) {
			bad << '${r.where}: ${r.got} against ${r.want}'
		}
	}
	assert bad.len == 0, bad.join('\n')
}

// R12-08's Table 4: μ with μ₂ = 1, in µPa·s.
fn test_viscosity_table_4() {
	pts := [
		['298.15', '998', '889.735100'],
		['298.15', '1200', '1437.649467'],
		['373.15', '1000', '307.883622'],
		['433.15', '1', '14.538324'],
		['433.15', '1000', '217.685358'],
		['873.15', '1', '32.619287'],
		['873.15', '100', '35.802262'],
		['873.15', '600', '77.430195'],
		['1173.15', '1', '44.217245'],
		['1173.15', '100', '47.640433'],
		['1173.15', '400', '64.154608'],
	]
	mut rows := []Row{}
	for pt in pts {
		mu := viscosity_industrial(pt[0].f64(), pt[1].f64())!
		rows << Row{mu * 1e6, pt[2], 'R12-08 Table 4 ${pt[0]} K, ${pt[1]} kg/m3'}
	}
	check_rows(rows)
}

// R12-08's Table 5: near the critical point, ξ (nm), μ₂ and μ (µPa·s),
// with IAPWS-95's ζ.
fn test_viscosity_table_5() {
	pts := [
		['122', '0.309247', '1.00000289', '25.520677'],
		['222', '1.571405', '1.00375120', '31.337589'],
		['272', '5.266522', '1.03416789', '36.228143'],
		['322', '16.590209', '1.09190440', '42.961579'],
		['372', '5.603768', '1.03665871', '45.688204'],
		['422', '1.876244', '1.00596332', '49.436256'],
	]
	mut rows := []Row{}
	for pt in pts {
		v := viscosity_parts(647.35, pt[0].f64())
		w := 'R12-08 Table 5 647.35 K, ${pt[0]} kg/m3'
		rows << Row{v.xi, pt[1], '${w} xi'}
		rows << Row{v.mu2, pt[2], '${w} mu2'}
		rows << Row{v.mu * 1e6, pt[3], '${w} mu'}
	}
	check_rows(rows)
	// and through the range-checked function
	assert math.abs(viscosity(647.35, 322.0)! * 1e6 - 42.961579) <= 5e-7
}

// R15-11's Tables 4 and 5: λ (mW/(m K)), and in Table 5 λ₁ and λ₂, with
// λ₀(647.35 K) = 51.5764797.
fn test_conductivity_tables_4_5() {
	mut rows := []Row{}
	for pt in [['298.15', '0', '18.4341883'], ['298.15', '998', '607.712868'],
		['298.15', '1200', '799.038144'], ['873.15', '0', '79.1034659']] {
		c := conductivity_parts(pt[0].f64(), pt[1].f64())
		rows << Row{c.lambda * 1e3, pt[2], 'R15-11 Table 4 ${pt[0]} K, ${pt[1]} kg/m3'}
		rows << Row{c.lambda2, '0', 'R15-11 Table 4 ${pt[0]} K, ${pt[1]} kg/m3 lambda2'}
	}
	pts := [
		['1', '1.0068497', '0.0001300', '51.9298924'],
		['122', '2.1445173', '20.3162320', '130.922885'],
		['222', '3.4840736', '188.091206', '367.787459'],
		['272', '4.2233708', '540.133176', '757.959776'],
		['322', '4.9681953', '1187.51354', '1443.75556'],
		['372', '5.6961250', '356.53333', '650.319402'],
		['422', '6.3973429', '118.931062', '448.883487'],
		['750', '11.5870532', '3.3419303', '600.961346'],
	]
	for pt in pts {
		c := conductivity_parts(647.35, pt[0].f64())
		w := 'R15-11 Table 5 647.35 K, ${pt[0]} kg/m3'
		rows << Row{c.lambda0, '51.5764797', '${w} lambda0'}
		rows << Row{c.lambda1, pt[1], '${w} lambda1'}
		rows << Row{c.lambda2, pt[2], '${w} lambda2'}
		rows << Row{c.lambda * 1e3, pt[3], '${w} lambda'}
	}
	check_rows(rows)
}

// R15-11's Tables 7 to 9: the industrial λ on IF97, regions 1 and 2 at
// (p, T) and region 3 at (ρ, T), with every quantity on the way: λ, λ₀, λ₁,
// λ₂, ρ, (∂ρ/∂p)_T at T and at T_R (kg m⁻³ MPa⁻¹), ξ, cp, cv (kJ/(kg K)),
// Z(y) and μ (µPa·s).
fn industrial_rows(st if97.State, want []string, w string) ![]Row {
	c := conductivity_industrial(st)!
	rho := 1.0 / st.v
	dr := rho / rho_star
	ztr := zeta_tr_industrial(dr)
	mu := viscosity_industrial(st.t, rho)!
	got := [c.lambda * 1e3, c.lambda0, c.lambda1, c.lambda2, rho, st.drhodp * 1e6,
		ztr * rho_star / p_star * 1e6, c.xi, st.cp / 1e3, st.cv / 1e3,
		z_of(cond_qd * c.xi, st.cp / st.cv, dr), mu * 1e6]
	names := ['lambda', 'lambda0', 'lambda1', 'lambda2', 'rho', 'drho/dp', 'drho/dp at TR', 'xi',
		'cp', 'cv', 'Z', 'mu']
	mut rows := []Row{}
	for k in 0 .. got.len {
		rows << Row{got[k], want[k], '${w} ${names[k]}'}
	}
	return rows
}

fn test_conductivity_industrial_tables_7_to_9() {
	mut rows := []Row{}
	rows << industrial_rows(if97.state_pt(20e6, 620.0)!, ['0.481485195e3', '0.484911627e2',
		'0.966869008e1', '0.126391714e2', '0.613227777e3', '0.520937820e1', '0.935037951', '0.377694973',
		'0.763433705e1', '0.303793441e1', '0.166942638', '0.709051068e2'],
		'R15-11 Table 7 20 MPa, 620 K')!
	rows << industrial_rows(if97.state_pt(50e6, 620.0)!, ['0.545038940e3', '0.484911627e2',
		'0.111212177e2', '0.575816285e1', '0.699226043e3', '0.184869007e1', '0.639306277', '0.189692422',
		'0.532047725e1', '0.291692653e1', '0.113592223', '0.841527945e2'],
		'R15-11 Table 7 50 MPa, 620 K')!
	rows << industrial_rows(if97.state_pt(0.3e6, 650.0)!, ['0.522311024e2', '0.518787461e2',
		'0.100678943e1', '0.129246457e-3', '0.100452141e1', '0.336351419e1', '0.223819386e1',
		'0.104305448e-2', '0.207010035e1', '0.159675313e1', '0.121437275e-2', '0.234877453e2'],
		'R15-11 Table 8 0.3 MPa, 650 K')!
	rows << industrial_rows(if97.state_pt(50e6, 800.0)!, ['0.177709914e3', '0.698329394e2',
		'0.244965343e1', '0.664341394e1', '0.218030012e3', '0.661484493e1', '0.312182530e1',
		'0.193491903', '0.590718707e1', '0.252343426e1', '0.137263826', '0.393727534e2'],
		'R15-11 Table 8 50 MPa, 800 K')!
	rows << industrial_rows(if97.region3(222.0, 647.35), ['0.366879411e3', '0.515764797e2',
		'0.348407362e1', '0.187183159e3', '0.222e3', '0.177778595e3', '0.311832789e1', '0.158223683e1',
		'0.101054488e3', '0.437466458e1', '0.217577777', '0.312204749e2'],
		'R15-11 Table 9 222 kg/m3, 647.35 K')!
	rows << industrial_rows(if97.region3(322.0, 647.35), ['0.124182415e4', '0.515764797e2',
		'0.496819532e1', '0.985582122e3', '0.322e3', '0.692651138e4', '0.275192511e1', '0.124722016e2',
		'0.312090124e4', '0.452163449e1', '0.322306729e-1', '0.393455495e2'],
		'R15-11 Table 9 322 kg/m3, 647.35 K')!
	check_rows(rows)
}

// R1-76's Table 1, column 4: σ by its equation, mN/m, every 5 °C from
// 0.01 °C to 370 °C.
fn test_surface_tension_table_1() {
	calc := ['75.65', '74.94', '74.22', '73.49', '72.74', '71.97', '71.19', '70.40', '69.60', '68.78',
		'67.94', '67.10', '66.24', '65.37', '64.48', '63.58', '62.67', '61.75', '60.82', '59.87',
		'58.91', '57.94', '56.96', '55.97', '54.97', '53.96', '52.93', '51.90', '50.86', '49.80',
		'48.74', '47.67', '46.59', '45.50', '44.41', '43.30', '42.19', '41.07', '39.95', '38.81',
		'37.67', '36.53', '35.38', '34.23', '33.07', '31.90', '30.74', '29.57', '28.39', '27.22',
		'26.04', '24.87', '23.69', '22.51', '21.34', '20.16', '18.99', '17.83', '16.66', '15.51',
		'14.36', '13.22', '12.09', '10.97', '9.86', '8.77', '7.70', '6.65', '5.63', '4.63', '3.67',
		'2.74', '1.88', '1.08', '0.39']
	mut rows := []Row{}
	for k, want in calc {
		tc := if k == 0 { 0.01 } else { 5.0 * f64(k) }
		rows << Row{surface_tension(tc + 273.15)! * 1e3, want, 'R1-76 Table 1 ${tc} C'}
	}
	check_rows(rows)
	assert surface_tension(critical_t)! == 0.0
}

// Each correlation's range (D11, as read in step 6).
fn test_the_ranges() {
	// the viscosity and conductivity from 250 K, to 1273.15 K at any density
	assert (viscosity(249.9, 1000.0) or { -1.0 }) == -1.0
	assert viscosity(250.0, 1000.0)! > 0.0
	assert thermal_conductivity(1273.15, 500.0)! > 0.0
	// above it, a dilute gas alone: to 2500 K for μ, 1500 K for λ
	assert viscosity(2400.0, 1.0)! > 0.0
	assert (viscosity(2400.0, 10.0) or { -1.0 }) == -1.0
	assert (viscosity(2600.0, 1.0) or { -1.0 }) == -1.0
	assert thermal_conductivity(1450.0, 1.0)! > 0.0
	assert (thermal_conductivity(1600.0, 1.0) or { -1.0 }) == -1.0
	// σ from 248.15 K to Tc
	assert (surface_tension(248.0) or { -1.0 }) == -1.0
	assert (surface_tension(647.1) or { -1.0 }) == -1.0
	mut kind_ok := false
	_ := viscosity(2600.0, 1.0) or {
		kind_ok = err is TransportError
		0.0
	}
	assert kind_ok
}

fn transport_rows(f json2.Any, key string) [][]f64 {
	list := f.as_map()[key] or { return [][]f64{} }
	return list.as_array().map(fn (r json2.Any) []f64 {
		return r.as_array().map(it.f64())
	})
}

// Against the iapws package (testdata/m5_transport.json,
// tools/gen_m5_transport.py): IAPWS-95's μ and λ at the package's own
// densities, so that the correlations alone are compared, measured within
// 1.0e-14 and 6.0e-14 (near Tc too); IF97's industrial μ and λ at reflux's
// IF97 states, within 4.9e-13 and 2.9e-11, the worst in region 5, where
// R15-11 sets λ₂ to 0 and reflux's λ is the lower; σ to 3.3e-16.
fn test_against_the_iapws_package() {
	text := os.read_file(os.join_path(os.dir(@FILE), '..', 'testdata', 'm5_transport.json'))!
	f := json2.decode[json2.Any](text)!
	mut wm, mut wk := 0.0, 0.0
	a := transport_rows(f, 'iapws95')
	assert a.len == 791
	for r in a {
		wm = math.max(wm, math.abs(viscosity(r[0], r[2])! / r[3] - 1.0))
		wk = math.max(wk, math.abs(thermal_conductivity(r[0], r[2])! / r[4] - 1.0))
	}
	assert wm <= 1e-13 && wk <= 2e-13, '${wm} ${wk}'
	mut im, mut ik := 0.0, 0.0
	b := transport_rows(f, 'if97')
	assert b.len == 775
	for r in b {
		st := if97.state_pt(r[1] * 1e6, r[0])!
		im = math.max(im, math.abs(viscosity_industrial(r[0], 1.0 / st.v)! / r[3] - 1.0))
		ik = math.max(ik, math.abs(conductivity_industrial(st)!.lambda / r[4] - 1.0))
	}
	assert im <= 2e-12 && ik <= 1e-10, '${im} ${ik}'
	mut ws := 0.0
	for r in transport_rows(f, 'tension') {
		ws = math.max(ws, math.abs(surface_tension(r[0])! / r[1] - 1.0))
	}
	assert ws <= 2e-15, '${ws}'
}

// IF97's region 3 close to the critical point (the fixture's
// "if97_critical": 647 K to 650 K and 21.8 MPa to 22.6 MPa, and the critical
// isochore 1 mK to 0.3 K above Tc), at the package's own densities so that
// the correlations alone are compared: where ρ̄ crosses Table 6's edges with
// ξ > 0, so the column Eq. (25) takes matters, and where ζ passes 1e3.
//
// Measured: μ within 6.7e-16; λ within 4.2e-12 on the grid and 1.3e-9 on
// the isochore at 1 mK, where (∂p/∂ρ)_T cancels to about 1/ζ of its terms
// (ζ = 1.2e5 there), so ζ, and λ₂ with it, carries about ζ·ε of rounding:
// the bound is 4e-11 + 4e-14·ζ.
fn test_industrial_near_the_critical_point() {
	text := os.read_file(os.join_path(os.dir(@FILE), '..', 'testdata', 'm5_transport.json'))!
	rows := transport_rows(json2.decode[json2.Any](text)!, 'if97_critical')
	assert rows.len == 539
	mut large_zeta, mut at_edge := 0, 0
	for r in rows {
		st := if97.region3(r[2], r[0])
		c := conductivity_industrial(st)!
		zeta := p_star / rho_star * st.drhodp
		dm := math.abs(viscosity_industrial(r[0], r[2])! / r[3] - 1.0)
		dk := math.abs(c.lambda / r[4] - 1.0)
		assert dm <= 1e-14, '${r[0]} K, ${r[2]} kg/m³: μ off by ${dm}'
		assert dk <= 4e-11 + 4e-14 * zeta, '${r[0]} K, ${r[2]} kg/m³: λ off by ${dk}'
		if zeta > 1e3 {
			large_zeta++
		}
		if c.xi > 0.0 && cond_a_rho.any(math.abs(r[2] / rho_star - it) < 0.01) {
			at_edge++
		}
	}
	assert large_zeta == 5 && at_edge == 19, '${large_zeta} ${at_edge}'
}

// R15-11's footnote 3: no critical enhancement in region 5, though its ξ
// there isn't always 0 (the package keeps λ₂ at 1075 K and 10 kPa, 2.9e-11
// of λ).
fn test_industrial_region_5_has_no_enhancement() {
	st := if97.state_pt(1e4, 1075.0)!
	assert st.region == 5
	c := conductivity_industrial(st)!
	assert c.lambda2 == 0.0 && c.xi == 0.0
	assert c.lambda == c.lambda0 * c.lambda1 * 1e-3
}
