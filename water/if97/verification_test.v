// IF97 against its releases' verification tables (M5a's criterion 4): R7-97,
// SR3-03 and SR5-05, every value within half a unit of its last printed
// figure, each backward equation's values as the backward equation gives them.
module if97

import math

fn half_unit(printed string) f64 {
	mant := printed.all_before('e')
	exp := if printed.contains('e') { printed.all_after('e').int() } else { 0 }
	decimals := mant.all_after('.').len
	return 0.5 * math.pow(10.0, f64(exp - decimals))
}

// Printed is a value computed, in the table's units, and the release's
// printing of it, with where it's printed.
struct Printed {
	got   f64
	want  string
	where string
}

fn check(rows []Printed) {
	mut bad := []string{}
	for row in rows {
		if !(math.abs(row.got - row.want.f64()) <= half_unit(row.want)) {
			bad << '${row.where}: ${row.got} against ${row.want}'
		}
	}
	assert bad.len == 0, bad.join('\n')
}

// A state's printed properties in the tables' units: v in m³/kg, h and u in
// kJ/kg, s and cp in kJ/(kg·K), w in m/s.
fn props(st State, want []string, table string, col string) []Printed {
	got := [st.v, st.h / 1e3, st.u / 1e3, st.s / 1e3, st.cp / 1e3, st.w]
	names := ['v', 'h', 'u', 's', 'cp', 'w']
	mut out := []Printed{}
	for k in 0 .. 6 {
		out << Printed{got[k], want[k], '${table} ${col} ${names[k]}'}
	}
	return out
}

// Table 5: region 1 at (300 K, 3 MPa), (300 K, 80 MPa), (500 K, 3 MPa).
fn test_table_5() {
	mut rows := []Printed{}
	rows << props(region1(3e6, 300.0), ['0.100215168e-2', '0.115331273e3', '0.112324818e3',
		'0.392294792', '0.417301218e1', '0.150773921e4'], 'Table 5', '300 K, 3 MPa')
	rows << props(region1(80e6, 300.0), ['0.971180894e-3', '0.184142828e3', '0.106448356e3',
		'0.368563852', '0.401008987e1', '0.163469054e4'], 'Table 5', '300 K, 80 MPa')
	rows << props(region1(3e6, 500.0), ['0.120241800e-2', '0.975542239e3', '0.971934985e3',
		'0.258041912e1', '0.465580682e1', '0.124071337e4'], 'Table 5', '500 K, 3 MPa')
	check(rows)
}

// Table 15: region 2 at (300 K, 0.0035 MPa), (700 K, 0.0035 MPa), (700 K,
// 30 MPa).
fn test_table_15() {
	mut rows := []Printed{}
	rows << props(region2(3500.0, 300.0), ['0.394913866e2', '0.254991145e4', '0.241169160e4',
		'0.852238967e1', '0.191300162e1', '0.427920172e3'], 'Table 15', '300 K, 0.0035 MPa')
	rows << props(region2(3500.0, 700.0), ['0.923015898e2', '0.333568375e4', '0.301262819e4',
		'0.101749996e2', '0.208141274e1', '0.644289068e3'], 'Table 15', '700 K, 0.0035 MPa')
	rows << props(region2(30e6, 700.0), ['0.542946619e-2', '0.263149474e4', '0.246861076e4',
		'0.517540298e1', '0.103505092e2', '0.480386523e3'], 'Table 15', '700 K, 30 MPa')
	check(rows)
}

// Table 18: the metastable vapour at (450 K, 1 MPa), (440 K, 1 MPa),
// (450 K, 1.5 MPa).
fn test_table_18() {
	mut rows := []Printed{}
	rows << props(region2_metastable(1e6, 450.0), ['0.192516540', '0.276881115e4', '0.257629461e4',
		'0.656660377e1', '0.276349265e1', '0.498408101e3'], 'Table 18', '450 K, 1 MPa')
	rows << props(region2_metastable(1e6, 440.0), ['0.186212297', '0.274015123e4', '0.255393894e4',
		'0.650218759e1', '0.298166443e1', '0.489363295e3'], 'Table 18', '440 K, 1 MPa')
	rows << props(region2_metastable(1.5e6, 450.0), ['0.121685206', '0.272134539e4', '0.253881758e4',
		'0.629170440e1', '0.362795578e1', '0.481941819e3'], 'Table 18', '450 K, 1.5 MPa')
	check(rows)
}

// Table 33: region 3 at (650 K, 500 kg/m³), (650 K, 200 kg/m³), (750 K,
// 500 kg/m³); p in MPa, then h, u, s, cp and w.
fn test_table_33() {
	mut rows := []Printed{}
	cases := [[650.0, 500.0], [650.0, 200.0], [750.0, 500.0]]
	wants := [
		['0.255837018e2', '0.186343019e4', '0.181226279e4', '0.405427273e1', '0.138935717e2',
			'0.502005554e3'],
		['0.222930643e2', '0.237512401e4', '0.226365868e4', '0.485438792e1', '0.446579342e2',
			'0.383444594e3'],
		['0.783095639e2', '0.225868845e4', '0.210206932e4', '0.446971906e1', '0.634165359e1',
			'0.760696041e3'],
	]
	for k, c in cases {
		st := region3(c[1], c[0])
		got := [st.p / 1e6, st.h / 1e3, st.u / 1e3, st.s / 1e3, st.cp / 1e3, st.w]
		for m, name in ['p', 'h', 'u', 's', 'cp', 'w'] {
			rows << Printed{got[m], wants[k][m], 'Table 33 ${c[0]} K, ${c[1]} kg/m3 ${name}'}
		}
	}
	check(rows)
}

// Tables 35 and 36: region 4's ps(T) and Ts(p).
fn test_tables_35_36() {
	check([
		Printed{psat97(300.0) / 1e6, '0.353658941e-2', 'Table 35 300 K'},
		Printed{psat97(500.0) / 1e6, '0.263889776e1', 'Table 35 500 K'},
		Printed{psat97(600.0) / 1e6, '0.123443146e2', 'Table 35 600 K'},
		Printed{tsat97(0.1e6), '0.372755919e3', 'Table 36 0.1 MPa'},
		Printed{tsat97(1e6), '0.453035632e3', 'Table 36 1 MPa'},
		Printed{tsat97(10e6), '0.584149488e3', 'Table 36 10 MPa'},
	])
}

// Table 42: region 5 at (1500 K, 0.5 MPa), (1500 K, 30 MPa), (2000 K,
// 30 MPa).
fn test_table_42() {
	mut rows := []Printed{}
	rows << props(region5(0.5e6, 1500.0), ['0.138455090e1', '0.521976855e4', '0.452749310e4',
		'0.965408875e1', '0.261609445e1', '0.917068690e3'], 'Table 42', '1500 K, 0.5 MPa')
	rows << props(region5(30e6, 1500.0), ['0.230761299e-1', '0.516723514e4', '0.447495124e4',
		'0.772970133e1', '0.272724317e1', '0.928548002e3'], 'Table 42', '1500 K, 30 MPa')
	rows << props(region5(30e6, 2000.0), ['0.311385219e-1', '0.657122604e4', '0.563707038e4',
		'0.853640523e1', '0.288569882e1', '0.106736948e4'], 'Table 42', '2000 K, 30 MPa')
	check(rows)
}

// B23's point (R7-97, Section 4): T = 623.15 K, p = 16.5291643 MPa.
fn test_b23_point() {
	check([
		Printed{p23(623.15) / 1e6, '0.165291643e2', 'B23 p(623.15 K)'},
		Printed{t23(0.165291643e8), '0.623150000e3', 'B23 T(16.5291643 MPa)'},
	])
}

// Tables 7 and 9: region 1's backward T(p, h) and T(p, s).
fn test_tables_7_9() {
	check([
		Printed{t1_ph(3e6, 500e3), '0.391798509e3', 'Table 7 3 MPa, 500 kJ/kg'},
		Printed{t1_ph(80e6, 500e3), '0.378108626e3', 'Table 7 80 MPa, 500 kJ/kg'},
		Printed{t1_ph(80e6, 1500e3), '0.611041229e3', 'Table 7 80 MPa, 1500 kJ/kg'},
		Printed{t1_ps(3e6, 0.5e3), '0.307842258e3', 'Table 9 3 MPa, 0.5 kJ/(kg K)'},
		Printed{t1_ps(80e6, 0.5e3), '0.309979785e3', 'Table 9 80 MPa, 0.5 kJ/(kg K)'},
		Printed{t1_ps(80e6, 3e3), '0.565899909e3', 'Table 9 80 MPa, 3 kJ/(kg K)'},
	])
}

// B2bc's point (R7-97, Section 6.3.1): p = 100 MPa, h = 3516.004323 kJ/kg.
fn test_b2bc_point() {
	check([
		Printed{p2bc(0.3516004323e7) / 1e6, '0.100000000e3', 'B2bc p(3516.004323 kJ/kg)'},
		Printed{h2bc(100e6) / 1e3, '0.3516004323e4', 'B2bc h(100 MPa)'},
	])
}

// Tables 24 and 29: region 2's backward T(p, h) and T(p, s), three points
// in each subregion.
fn test_tables_24_29() {
	check([
		Printed{t2_ph(0.001e6, 3000e3), '0.534433241e3', 'Table 24 2a 0.001 MPa, 3000 kJ/kg'},
		Printed{t2_ph(3e6, 3000e3), '0.575373370e3', 'Table 24 2a 3 MPa, 3000 kJ/kg'},
		Printed{t2_ph(3e6, 4000e3), '0.101077577e4', 'Table 24 2a 3 MPa, 4000 kJ/kg'},
		Printed{t2_ph(5e6, 3500e3), '0.801299102e3', 'Table 24 2b 5 MPa, 3500 kJ/kg'},
		Printed{t2_ph(5e6, 4000e3), '0.101531583e4', 'Table 24 2b 5 MPa, 4000 kJ/kg'},
		Printed{t2_ph(25e6, 3500e3), '0.875279054e3', 'Table 24 2b 25 MPa, 3500 kJ/kg'},
		Printed{t2_ph(40e6, 2700e3), '0.743056411e3', 'Table 24 2c 40 MPa, 2700 kJ/kg'},
		Printed{t2_ph(60e6, 2700e3), '0.791137067e3', 'Table 24 2c 60 MPa, 2700 kJ/kg'},
		Printed{t2_ph(60e6, 3200e3), '0.882756860e3', 'Table 24 2c 60 MPa, 3200 kJ/kg'},
		Printed{t2_ps(0.1e6, 7.5e3), '0.399517097e3', 'Table 29 2a 0.1 MPa, 7.5 kJ/(kg K)'},
		Printed{t2_ps(0.1e6, 8e3), '0.514127081e3', 'Table 29 2a 0.1 MPa, 8 kJ/(kg K)'},
		Printed{t2_ps(2.5e6, 8e3), '0.103984917e4', 'Table 29 2a 2.5 MPa, 8 kJ/(kg K)'},
		Printed{t2_ps(8e6, 6e3), '0.600484040e3', 'Table 29 2b 8 MPa, 6 kJ/(kg K)'},
		Printed{t2_ps(8e6, 7.5e3), '0.106495556e4', 'Table 29 2b 8 MPa, 7.5 kJ/(kg K)'},
		Printed{t2_ps(90e6, 6e3), '0.103801126e4', 'Table 29 2b 90 MPa, 6 kJ/(kg K)'},
		Printed{t2_ps(20e6, 5.75e3), '0.697992849e3', 'Table 29 2c 20 MPa, 5.75 kJ/(kg K)'},
		Printed{t2_ps(80e6, 5.25e3), '0.854011484e3', 'Table 29 2c 80 MPa, 5.25 kJ/(kg K)'},
		Printed{t2_ps(80e6, 5.75e3), '0.949017998e3', 'Table 29 2c 80 MPa, 5.75 kJ/(kg K)'},
	])
}

// SR3-03: h3ab's point (Section 3.1), and Tables 5, 8, 12, 15, 18 and 20.
fn test_sr3_03() {
	check([
		Printed{h3ab(25e6) / 1e3, '2.095936454e3', 'SR3-03 h3ab(25 MPa)'},
		Printed{t3_ph(20e6, 1700e3), '6.293083892e2', 'Table 5 3a 20 MPa, 1700 kJ/kg'},
		Printed{t3_ph(50e6, 2000e3), '6.905718338e2', 'Table 5 3a 50 MPa, 2000 kJ/kg'},
		Printed{t3_ph(100e6, 2100e3), '7.336163014e2', 'Table 5 3a 100 MPa, 2100 kJ/kg'},
		Printed{t3_ph(20e6, 2500e3), '6.418418053e2', 'Table 5 3b 20 MPa, 2500 kJ/kg'},
		Printed{t3_ph(50e6, 2400e3), '7.351848618e2', 'Table 5 3b 50 MPa, 2400 kJ/kg'},
		Printed{t3_ph(100e6, 2700e3), '8.420460876e2', 'Table 5 3b 100 MPa, 2700 kJ/kg'},
		Printed{v3_ph(20e6, 1700e3), '1.749903962e-3', 'Table 8 3a 20 MPa, 1700 kJ/kg'},
		Printed{v3_ph(50e6, 2000e3), '1.908139035e-3', 'Table 8 3a 50 MPa, 2000 kJ/kg'},
		Printed{v3_ph(100e6, 2100e3), '1.676229776e-3', 'Table 8 3a 100 MPa, 2100 kJ/kg'},
		Printed{v3_ph(20e6, 2500e3), '6.670547043e-3', 'Table 8 3b 20 MPa, 2500 kJ/kg'},
		Printed{v3_ph(50e6, 2400e3), '2.801244590e-3', 'Table 8 3b 50 MPa, 2400 kJ/kg'},
		Printed{v3_ph(100e6, 2700e3), '2.404234998e-3', 'Table 8 3b 100 MPa, 2700 kJ/kg'},
		Printed{t3_ps(20e6, 3.8e3), '6.282959869e2', 'Table 12 3a 20 MPa, 3.8 kJ/(kg K)'},
		Printed{t3_ps(50e6, 3.6e3), '6.297158726e2', 'Table 12 3a 50 MPa, 3.6 kJ/(kg K)'},
		Printed{t3_ps(100e6, 4.0e3), '7.056880237e2', 'Table 12 3a 100 MPa, 4.0 kJ/(kg K)'},
		Printed{t3_ps(20e6, 5.0e3), '6.401176443e2', 'Table 12 3b 20 MPa, 5.0 kJ/(kg K)'},
		Printed{t3_ps(50e6, 4.5e3), '7.163687517e2', 'Table 12 3b 50 MPa, 4.5 kJ/(kg K)'},
		Printed{t3_ps(100e6, 5.0e3), '8.474332825e2', 'Table 12 3b 100 MPa, 5.0 kJ/(kg K)'},
		Printed{v3_ps(20e6, 3.8e3), '1.733791463e-3', 'Table 15 3a 20 MPa, 3.8 kJ/(kg K)'},
		Printed{v3_ps(50e6, 3.6e3), '1.469680170e-3', 'Table 15 3a 50 MPa, 3.6 kJ/(kg K)'},
		Printed{v3_ps(100e6, 4.0e3), '1.555893131e-3', 'Table 15 3a 100 MPa, 4.0 kJ/(kg K)'},
		Printed{v3_ps(20e6, 5.0e3), '6.262101987e-3', 'Table 15 3b 20 MPa, 5.0 kJ/(kg K)'},
		Printed{v3_ps(50e6, 4.5e3), '2.332634294e-3', 'Table 15 3b 50 MPa, 4.5 kJ/(kg K)'},
		Printed{v3_ps(100e6, 5.0e3), '2.449610757e-3', 'Table 15 3b 100 MPa, 5.0 kJ/(kg K)'},
		Printed{p3sat_h(1700e3) / 1e6, '1.724175718e1', 'Table 18 1700 kJ/kg'},
		Printed{p3sat_h(2000e3) / 1e6, '2.193442957e1', 'Table 18 2000 kJ/kg'},
		Printed{p3sat_h(2400e3) / 1e6, '2.018090839e1', 'Table 18 2400 kJ/kg'},
		Printed{p3sat_s(3.8e3) / 1e6, '1.687755057e1', 'Table 20 3.8 kJ/(kg K)'},
		Printed{p3sat_s(4.2e3) / 1e6, '2.164451789e1', 'Table 20 4.2 kJ/(kg K)'},
		Printed{p3sat_s(5.2e3) / 1e6, '1.668968482e1', 'Table 20 5.2 kJ/(kg K)'},
	])
}

// SR5-05's Tables 3 and 11: the subregion boundaries.
fn test_sr5_05_boundaries() {
	check([
		Printed{t3ab(40e6), '6.930341408e2', 'Table 3 T3ab(40 MPa)'},
		Printed{t3cd(25e6), '6.493659208e2', 'Table 3 T3cd(25 MPa)'},
		Printed{t3ef(40e6), '7.139593992e2', 'Table 3 T3ef(40 MPa)'},
		Printed{t3gh(23e6), '6.498873759e2', 'Table 3 T3gh(23 MPa)'},
		Printed{t3ij(23e6), '6.515778091e2', 'Table 3 T3ij(23 MPa)'},
		Printed{t3jk(23e6), '6.558338344e2', 'Table 3 T3jk(23 MPa)'},
		Printed{t3mn(22.8e6), '6.496054133e2', 'Table 3 T3mn(22.8 MPa)'},
		Printed{t3op(22.8e6), '6.500106943e2', 'Table 3 T3op(22.8 MPa)'},
		Printed{t3qu(22e6), '6.456355027e2', 'Table 3 T3qu(22 MPa)'},
		Printed{t3rx(22e6), '6.482622754e2', 'Table 3 T3rx(22 MPa)'},
		Printed{t3uv(22.3e6), '6.477996121e2', 'Table 11 T3uv(22.3 MPa)'},
		Printed{t3wx(22.3e6), '6.482049480e2', 'Table 11 T3wx(22.3 MPa)'},
	])
}

// SR5-05's Tables 5 and 13: two points in each subregion, each through the
// subregion's own equation and through the choice of subregion, which must
// choose it.
fn test_sr5_05_volumes() {
	rows := [
		['a', '50', '630', '1.470853100e-3'],
		['a', '80', '670', '1.503831359e-3'],
		['b', '50', '710', '2.204728587e-3'],
		['b', '80', '750', '1.973692940e-3'],
		['c', '20', '630', '1.761696406e-3'],
		['c', '30', '650', '1.819560617e-3'],
		['d', '26', '656', '2.245587720e-3'],
		['d', '30', '670', '2.506897702e-3'],
		['e', '26', '661', '2.970225962e-3'],
		['e', '30', '675', '3.004627086e-3'],
		['f', '26', '671', '5.019029401e-3'],
		['f', '30', '690', '4.656470142e-3'],
		['g', '23.6', '649', '2.163198378e-3'],
		['g', '24', '650', '2.166044161e-3'],
		['h', '23.6', '652', '2.651081407e-3'],
		['h', '24', '654', '2.967802335e-3'],
		['i', '23.6', '653', '3.273916816e-3'],
		['i', '24', '655', '3.550329864e-3'],
		['j', '23.5', '655', '4.545001142e-3'],
		['j', '24', '660', '5.100267704e-3'],
		['k', '23', '660', '6.109525997e-3'],
		['k', '24', '670', '6.427325645e-3'],
		['l', '22.6', '646', '2.117860851e-3'],
		['l', '23', '646', '2.062374674e-3'],
		['m', '22.6', '648.6', '2.533063780e-3'],
		['m', '22.8', '649.3', '2.572971781e-3'],
		['n', '22.6', '649.0', '2.923432711e-3'],
		['n', '22.8', '649.7', '2.913311494e-3'],
		['o', '22.6', '649.1', '3.131208996e-3'],
		['o', '22.8', '649.9', '3.221160278e-3'],
		['p', '22.6', '649.4', '3.715596186e-3'],
		['p', '22.8', '650.2', '3.664754790e-3'],
		['q', '21.1', '640', '1.970999272e-3'],
		['q', '21.8', '643', '2.043919161e-3'],
		['r', '21.1', '644', '5.251009921e-3'],
		['r', '21.8', '648', '5.256844741e-3'],
		['s', '19.1', '635', '1.932829079e-3'],
		['s', '20', '638', '1.985387227e-3'],
		['t', '17', '626', '8.483262001e-3'],
		['t', '20', '640', '6.227528101e-3'],
		['u', '21.5', '644.6', '2.268366647e-3'],
		['u', '22.0', '646.1', '2.296350553e-3'],
		['v', '22.5', '648.6', '2.832373260e-3'],
		['v', '22.3', '647.9', '2.811424405e-3'],
		['w', '22.15', '647.5', '3.694032281e-3'],
		['w', '22.3', '648.1', '3.622226305e-3'],
		['x', '22.11', '648.0', '4.528072649e-3'],
		['x', '22.3', '649.0', '4.556905799e-3'],
		['y', '22.0', '646.84', '2.698354719e-3'],
		['y', '22.064', '647.05', '2.717655648e-3'],
		['z', '22.0', '646.89', '3.798732962e-3'],
	]
	mut checks := []Printed{}
	mut wrong := []string{}
	for row in rows {
		code := row[0][0]
		p := row[1].f64() * 1e6
		t := row[2].f64()
		table := if code <= `t` { 'Table 5' } else { 'Table 13' }
		checks << Printed{v3_sub(sub_of(code), p, t), row[3], '${table} v3${row[0]}(${row[1]} MPa, ${row[2]} K)'}
		got := subregion3(p, t, .by_t)
		if got != code {
			wrong << '${row[1]} MPa, ${row[2]} K: 3${got.ascii_str()}, not 3${row[0]}'
		}
	}
	check(checks)
	assert wrong.len == 0, wrong.join('\n')
	// Table 13's last row, v3z(22.064 MPa, 647.15 K), is printed 3.701940010e-3,
	// but Eq. (4) with Table A2.6's coefficients is 3.7019400094668e-3 (mpmath at
	// 50 digits): the printing is 0.53 units of its last figure from its own
	// equation. Its terms sum to 1/2.8e4 of their magnitudes and e = 4, so
	// rounding allows 4 × 2.8e4 × 2.2e-16 = 2.5e-11 of v; reflux is 1.7e-12 off.
	z := v3_sub(sub_of(`z`), 22.064e6, 647.15)
	assert math.abs(z / 3.7019400094668e-3 - 1.0) < 2.5e-11, '${z}'
	assert subregion3(22.064e6, 647.15, .by_t) == `z`
}
