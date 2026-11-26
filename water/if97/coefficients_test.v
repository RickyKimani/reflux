// Internal test (module if97): every coefficient and exponent transcribed
// from the releases (coefficients*.v, D10) is the iapws package's, as its
// source has them (testdata/m5_if97.json, tools/gen_m5_if97.py): two
// transcriptions of the same tables, made apart.
module if97

import json2
import os

fn fixture() !json2.Any {
	path := os.join_path(os.dir(@FILE), '..', '..', 'testdata', 'm5_if97.json')
	text := os.read_file(path) or {
		return error('${path}: ${err}. Regenerate with `uv run tools/gen_m5_if97.py`.')
	}
	return json2.decode[json2.Any](text)!
}

fn at(a json2.Any, keys ...string) json2.Any {
	mut cur := a
	for k in keys {
		cur = cur.as_map()[k] or { return json2.Any(json2.null) }
	}
	return cur
}

fn nums(a json2.Any) []f64 {
	return a.as_array().map(it.f64())
}

fn fl(a []int) []f64 {
	return a.map(f64(it))
}

fn abs_all(a []f64) []f64 {
	return a.map(if it < 0 { -it } else { it })
}

struct Same {
	what   string
	ours   []f64
	theirs []f64
}

// The R7-97 and SR3-03 tables kept as one array each in the iapws
// package's constants module: (reflux's name, the package's prefix).
const by_prefix = [
	['r1', 'Region1_L', 'Region1_n'],
	['t1ph', 'Backward1_T_Ph_L', 'Backward1_T_Ph_n'],
	['t1ps', 'Backward1_T_Ps_L', 'Backward1_T_Ps_n'],
	['r2', 'Region2_', 'Region2_nr'],
	['t2aph', 'Backward2a_T_Ph_L', 'Backward2a_T_Ph_n'],
	['t2bph', 'Backward2b_T_Ph_L', 'Backward2b_T_Ph_n'],
	['t2cph', 'Backward2c_T_Ph_L', 'Backward2c_T_Ph_n'],
	['t2aps', 'Backward2a_T_Ps_L', 'Backward2a_T_Ps_n'],
	['t2bps', 'Backward2b_T_Ps_L', 'Backward2b_T_Ps_n'],
	['t2cps', 'Backward2c_T_Ps_L', 'Backward2c_T_Ps_n'],
	['r3', 'Region3_L', 'Region3_n'],
	['r5', 'Region5_', 'Region5_nr'],
	['t3aph', 'Backward3a_T_Ph_L', 'Backward3a_T_Ph_n'],
	['t3bph', 'Backward3b_T_Ph_L', 'Backward3b_T_Ph_n'],
	['v3aph', 'Backward3a_v_Ph_L', 'Backward3a_v_Ph_n'],
	['v3bph', 'Backward3b_v_Ph_L', 'Backward3b_v_Ph_n'],
	['t3aps', 'Backward3a_T_Ps_L', 'Backward3a_T_Ps_n'],
	['t3bps', 'Backward3b_T_Ps_L', 'Backward3b_T_Ps_n'],
	['v3aps', 'Backward3a_v_Ps_L', 'Backward3a_v_Ps_n'],
	['v3bps', 'Backward3b_v_Ps_L', 'Backward3b_v_Ps_n'],
	['p3sath', 'PSat_h_L', 'PSat_h_n'],
	['p3sats', 'PSat_s_L', 'PSat_s_n'],
]

fn ours(name string) ([]f64, []f64, []f64) {
	return match name {
		'r1' { fl(r1_i), fl(r1_j), r1_n }
		't1ph' { fl(t1ph_i), fl(t1ph_j), t1ph_n }
		't1ps' { fl(t1ps_i), fl(t1ps_j), t1ps_n }
		'r2' { fl(r2_i), fl(r2_j), r2_n }
		't2aph' { fl(t2aph_i), fl(t2aph_j), t2aph_n }
		't2bph' { fl(t2bph_i), fl(t2bph_j), t2bph_n }
		't2cph' { fl(t2cph_i), fl(t2cph_j), t2cph_n }
		't2aps' { t2aps_i, fl(t2aps_j), t2aps_n }
		't2bps' { fl(t2bps_i), fl(t2bps_j), t2bps_n }
		't2cps' { fl(t2cps_i), fl(t2cps_j), t2cps_n }
		'r3' { fl(r3_i), fl(r3_j), r3_n }
		'r5' { fl(r5_i), fl(r5_j), r5_n }
		't3aph' { fl(t3aph_i), fl(t3aph_j), t3aph_n }
		't3bph' { fl(t3bph_i), fl(t3bph_j), t3bph_n }
		'v3aph' { fl(v3aph_i), fl(v3aph_j), v3aph_n }
		'v3bph' { fl(v3bph_i), fl(v3bph_j), v3bph_n }
		't3aps' { fl(t3aps_i), fl(t3aps_j), t3aps_n }
		't3bps' { fl(t3bps_i), fl(t3bps_j), t3bps_n }
		'v3aps' { fl(v3aps_i), fl(v3aps_j), v3aps_n }
		'v3bps' { fl(v3bps_i), fl(v3bps_j), v3bps_n }
		'p3sath' { fl(p3sath_i), fl(p3sath_j), p3sath_n }
		'p3sats' { fl(p3sats_i), fl(p3sats_j), p3sats_n }
		else { []f64{}, []f64{}, []f64{} }
	}
}

fn test_the_coefficients_are_the_iapws_package_s() {
	f := fixture()!
	c := at(f, 'coefficients', 'const')
	local := at(f, 'coefficients', 'local')
	inline := at(f, 'coefficients', 'inline')
	mut pairs := []Same{}
	for row in by_prefix {
		i, j, n := ours(row[0])
		// Region 2's and 5's are Ir and Jr, the rest Li and Lj.
		ik := if row[1].ends_with('_') { '${row[1]}Ir' } else { '${row[1]}i' }
		jk := if row[1].ends_with('_') { '${row[1]}Jr' } else { '${row[1]}j' }
		pairs << Same{'${row[0]} I', i, nums(at(c, ik))}
		pairs << Same{'${row[0]} J', j, nums(at(c, jk))}
		pairs << Same{'${row[0]} n', n, nums(at(c, row[2]))}
	}
	// The ideal-gas parts.
	pairs << Same{'r2o J', fl(r2o_j), nums(at(c, 'Region2_cp0_Jo'))}
	pairs << Same{'r2o n', r2o_n, nums(at(c, 'Region2_cp0_no'))}
	pairs << Same{'r5o J', fl(r5o_j), nums(at(c, 'Region5_cp0_Jo'))}
	pairs << Same{'r5o n', r5o_n, nums(at(c, 'Region5_cp0_no'))}
	// Region 4 (its tuple has a 0 in front, to count from 1), and B23 in its
	// two functions.
	pairs << Same{'r4 n', r4_n, nums(at(local, '_PSat_T', 'n'))[1..]}
	pairs << Same{'r4 n (Ts)', r4_n, nums(at(local, '_TSat_P', 'n'))[1..]}
	pairs << Same{'b23 n1-n3', b23_n[..3], nums(at(local, '_P23_T', 'n'))}
	// The package's n₅ in T23(p) has an 8 dropped: 13.9188397787 where R7-97's
	// Table 1 prints 0.139 188 397 788 70 × 10², 13.918839778870. Ours is the
	// release's, which B23's point holds (verification_test.v); the package's
	// moves T23 by 1e-10 relative. To report upstream.
	t_p := nums(at(local, '_t_P', 'n'))
	assert t_p[2] == 13.9188397787 && b23_n[4] == 13.91883977887
	pairs << Same{'b23 n3-n4', [b23_n[2], b23_n[3]], t_p[..2]}
	// SR5-05's boundaries: T3ab, T3op and T3wx in the constants module, the
	// rest as ng[code] (their Iᵢ are 0, 1, 2, ...).
	pairs << Same{'t3ab I', fl(t3ab_i), nums(at(c, 'tab_P_Li'))}
	pairs << Same{'t3ab n', t3ab_n, nums(at(c, 'tab_P_n'))}
	pairs << Same{'t3op I', fl(t3op_i), nums(at(c, 'top_P_Li'))}
	pairs << Same{'t3op n', t3op_n, nums(at(c, 'top_P_n'))}
	pairs << Same{'t3wx I', fl(t3wx_i), nums(at(c, 'twx_P_Li'))}
	pairs << Same{'t3wx n', t3wx_n, nums(at(c, 'twx_P_n'))}
	bounds := {
		'cd': [fl(t3cd_i), t3cd_n]
		'gh': [fl(t3gh_i), t3gh_n]
		'ij': [fl(t3ij_i), t3ij_n]
		'jk': [fl(t3jk_i), t3jk_n]
		'mn': [fl(t3mn_i), t3mn_n]
		'qu': [fl(t3qu_i), t3qu_n]
		'rx': [fl(t3rx_i), t3rx_n]
		'uv': [fl(t3uv_i), t3uv_n]
	}
	for code, b in bounds {
		pairs << Same{'t3${code} I', b[0], []f64{len: b[0].len, init: f64(index)}}
		pairs << Same{'t3${code} n', b[1], nums(at(local, '_txx_P', 'ng', code))}
	}
	// SR5-05's subregions: terms, and Tables 4 and 12's reducing quantities
	// and parameters, [v*, p* (MPa), T*, a, b, c, d, e] in the package.
	for s in subs3 {
		code := s.code.ascii_str()
		pairs << Same{'v3${code} I', fl(s.i), nums(at(local, '_Backward3x_v_PT', 'Li', code))}
		pairs << Same{'v3${code} J', fl(s.j), nums(at(local, '_Backward3x_v_PT', 'Lj', code))}
		pairs << Same{'v3${code} n', s.n, nums(at(local, '_Backward3x_v_PT', 'n', code))}
		par := nums(at(local, '_Backward3x_v_PT', 'par', code))
		// 3n's exponential form, Eq. (5), has no c, d and e: Table 4 prints
		// dashes there, and the package 0s.
		mut want := [s.vs, s.ps / 1e6, s.ts, s.a, s.b, s.c, s.d, f64(s.e)]
		if s.e == 0 {
			want = [s.vs, s.ps / 1e6, s.ts, s.a, s.b, 0.0, 0.0, 0.0]
		}
		pairs << Same{'v3${code} reducing', want, par}
	}
	mut bad := []string{}
	for pr in pairs {
		if pr.ours != pr.theirs {
			bad << '${pr.what}: ${pr.ours} against ${pr.theirs}'
		}
	}
	// The functions that write theirs inline: each of ours is among their
	// numbers, as magnitudes (a minus sign there is an operator).
	inl := {
		'B2bc n1 to n3': [abs_all(b2bc_n[..3]), abs_all(nums(at(inline, '_P_2bc')))]
		'B2bc n3 to n5': [abs_all(b2bc_n[2..]), abs_all(nums(at(inline, '_hbc_P')))]
		'h3ab':          [abs_all(h3ab_n), abs_all(nums(at(inline, '_h_3ab')))]
		'region 3 n1':   [[r3_n1], nums(at(inline, '_Region3'))]
		'T3ef slope':    [[3.727888004], nums(at(inline, '_tef_P'))]
	}
	for what, pr in inl {
		for x in pr[0] {
			if x !in pr[1] {
				bad << '${what}: ${x} not among ${pr[1]}'
			}
		}
	}
	assert bad.len == 0, bad.join('\n')
	// Every table was compared: 22 of three arrays, 4 ideal-gas arrays, 4 for
	// region 4 and B23, 6 + 16 boundary arrays and 26 subregions of four.
	assert pairs.len == 22 * 3 + 4 + 4 + 6 + 16 + 26 * 4
}
