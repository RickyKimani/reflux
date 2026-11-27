module if97

import math

// Region 3's v(p, T), IAPWS SR5-05(2016): the subregion boundaries, the
// choice of subregion (its Tables 2 and 10), and each subregion's equation.

// A subregion's v(p, T), Eq. (4): v/v* = [Σ nᵢ ((π − a)^c)^Iᵢ ((θ − b)^d)^Jᵢ]^e,
// with π = p/p* and θ = T/T*, or for 3n, Eq. (5): ln(v/v*) = Σ nᵢ (π − a)^Iᵢ
// (θ − b)^Jᵢ. Reducing quantities and parameters from Tables 4 and 12.
struct Sub {
	code u8
	vs   f64 // v*, m³/kg
	ps   f64 // p*, Pa
	ts   f64 // T*, K
	a    f64
	b    f64
	c    f64 // 1 or 0.5
	d    f64 // 1 or 0.25
	e    int // 1 or 4; 0 for 3n's exponential form
	i    []int
	j    []int
	n    []f64
}

// Table 4: subregions 3a to 3t.
const sub3a = Sub{`a`, 0.0024, 100e6, 760.0, 0.085, 0.817, 1.0, 1.0, 1, v3a_i, v3a_j, v3a_n}
const sub3b = Sub{`b`, 0.0041, 100e6, 860.0, 0.280, 0.779, 1.0, 1.0, 1, v3b_i, v3b_j, v3b_n}
const sub3c = Sub{`c`, 0.0022, 40e6, 690.0, 0.259, 0.903, 1.0, 1.0, 1, v3c_i, v3c_j, v3c_n}
const sub3d = Sub{`d`, 0.0029, 40e6, 690.0, 0.559, 0.939, 1.0, 1.0, 4, v3d_i, v3d_j, v3d_n}
const sub3e = Sub{`e`, 0.0032, 40e6, 710.0, 0.587, 0.918, 1.0, 1.0, 1, v3e_i, v3e_j, v3e_n}
const sub3f = Sub{`f`, 0.0064, 40e6, 730.0, 0.587, 0.891, 0.5, 1.0, 4, v3f_i, v3f_j, v3f_n}
const sub3g = Sub{`g`, 0.0027, 25e6, 660.0, 0.872, 0.971, 1.0, 1.0, 4, v3g_i, v3g_j, v3g_n}
const sub3h = Sub{`h`, 0.0032, 25e6, 660.0, 0.898, 0.983, 1.0, 1.0, 4, v3h_i, v3h_j, v3h_n}
const sub3i = Sub{`i`, 0.0041, 25e6, 660.0, 0.910, 0.984, 0.5, 1.0, 4, v3i_i, v3i_j, v3i_n}
const sub3j = Sub{`j`, 0.0054, 25e6, 670.0, 0.875, 0.964, 0.5, 1.0, 4, v3j_i, v3j_j, v3j_n}
const sub3k = Sub{`k`, 0.0077, 25e6, 680.0, 0.802, 0.935, 1.0, 1.0, 1, v3k_i, v3k_j, v3k_n}
const sub3l = Sub{`l`, 0.0026, 24e6, 650.0, 0.908, 0.989, 1.0, 1.0, 4, v3l_i, v3l_j, v3l_n}
const sub3m = Sub{`m`, 0.0028, 23e6, 650.0, 1.000, 0.997, 1.0, 0.25, 1, v3m_i, v3m_j, v3m_n}
const sub3n = Sub{`n`, 0.0031, 23e6, 650.0, 0.976, 0.997, 1.0, 1.0, 0, v3n_i, v3n_j, v3n_n}
const sub3o = Sub{`o`, 0.0034, 23e6, 650.0, 0.974, 0.996, 0.5, 1.0, 1, v3o_i, v3o_j, v3o_n}
const sub3p = Sub{`p`, 0.0041, 23e6, 650.0, 0.972, 0.997, 0.5, 1.0, 1, v3p_i, v3p_j, v3p_n}
const sub3q = Sub{`q`, 0.0022, 23e6, 650.0, 0.848, 0.983, 1.0, 1.0, 4, v3q_i, v3q_j, v3q_n}
const sub3r = Sub{`r`, 0.0054, 23e6, 650.0, 0.874, 0.982, 1.0, 1.0, 1, v3r_i, v3r_j, v3r_n}
const sub3s = Sub{`s`, 0.0022, 21e6, 640.0, 0.886, 0.990, 1.0, 1.0, 4, v3s_i, v3s_j, v3s_n}
const sub3t = Sub{`t`, 0.0088, 20e6, 650.0, 0.803, 1.020, 1.0, 1.0, 1, v3t_i, v3t_j, v3t_n}

// Table 12: the auxiliary subregions 3u to 3z, near the critical point.
const sub3u = Sub{`u`, 0.0026, 23e6, 650.0, 0.902, 0.988, 1.0, 1.0, 1, v3u_i, v3u_j, v3u_n}
const sub3v = Sub{`v`, 0.0031, 23e6, 650.0, 0.960, 0.995, 1.0, 1.0, 1, v3v_i, v3v_j, v3v_n}
const sub3w = Sub{`w`, 0.0039, 23e6, 650.0, 0.959, 0.995, 1.0, 1.0, 4, v3w_i, v3w_j, v3w_n}
const sub3x = Sub{`x`, 0.0049, 23e6, 650.0, 0.910, 0.988, 1.0, 1.0, 1, v3x_i, v3x_j, v3x_n}
const sub3y = Sub{`y`, 0.0031, 22e6, 650.0, 0.996, 0.994, 1.0, 1.0, 4, v3y_i, v3y_j, v3y_n}
const sub3z = Sub{`z`, 0.0038, 22e6, 650.0, 0.993, 0.994, 1.0, 1.0, 4, v3z_i, v3z_j, v3z_n}

const subs3 = [sub3a, sub3b, sub3c, sub3d, sub3e, sub3f, sub3g, sub3h, sub3i, sub3j, sub3k, sub3l,
	sub3m, sub3n, sub3o, sub3p, sub3q, sub3r, sub3s, sub3t, sub3u, sub3v, sub3w, sub3x, sub3y,
	sub3z]

// Footnotes to Tables 2 and 10.
const p3cd = 19.00881189173929e6 // Pa
const p3_v264 = 21.93161551e6 // psat at v = 0.00264 m³/kg, Pa
const p3_v385 = 21.90096265e6 // psat at v = 0.00385 m³/kg, Pa

fn sub_of(code u8) Sub {
	return subs3[int(code - `a`)]
}

// v3_sub is a subregion's v(p, T), Eqs. (4) and (5), as published.
fn v3_sub(s Sub, p f64, t f64) f64 {
	x0 := p / s.ps - s.a
	y0 := t / s.ts - s.b
	if s.e == 0 {
		return s.vs * math.exp(sum_xy(s.i, s.j, s.n, x0, y0))
	}
	x := if s.c == 0.5 { math.sqrt(x0) } else { x0 }
	y := if s.d == 0.25 { math.sqrt(math.sqrt(y0)) } else { y0 }
	w := sum_xy(s.i, s.j, s.n, x, y)
	return if s.e == 4 { s.vs * w * w * w * w } else { s.vs * w }
}

// Eq. (1): T/T* = Σ nᵢ π^Iᵢ, T* = 1 K, p* = 1 MPa.
fn t3_poly(ii []int, nn []f64, p f64) f64 {
	pi := p / 1e6
	mut s := 0.0
	for k in 0 .. nn.len {
		s += nn[k] * ipow(pi, ii[k])
	}
	return s
}

// Eq. (2): T/T* = Σ nᵢ (ln π)^Iᵢ.
fn t3_log(ii []int, nn []f64, p f64) f64 {
	l := math.log(p / 1e6)
	mut s := 0.0
	for k in 0 .. nn.len {
		s += nn[k] * ipow(l, ii[k])
	}
	return s
}

// The subregion boundaries, T in K from p in Pa (Tables 1 and 9).
pub fn t3ab(p f64) f64 {
	return t3_log(t3ab_i, t3ab_n, p)
}

pub fn t3cd(p f64) f64 {
	return t3_poly(t3cd_i, t3cd_n, p)
}

// t3ef is Eq. (3): the saturation line's slope at the critical point
// through it, (∂θsat/∂π)c = 3.727888004.
pub fn t3ef(p f64) f64 {
	return 3.727888004 * (p / 1e6 - 22.064) + 647.096
}

pub fn t3gh(p f64) f64 {
	return t3_poly(t3gh_i, t3gh_n, p)
}

pub fn t3ij(p f64) f64 {
	return t3_poly(t3ij_i, t3ij_n, p)
}

pub fn t3jk(p f64) f64 {
	return t3_poly(t3jk_i, t3jk_n, p)
}

pub fn t3mn(p f64) f64 {
	return t3_poly(t3mn_i, t3mn_n, p)
}

pub fn t3op(p f64) f64 {
	return t3_log(t3op_i, t3op_n, p)
}

pub fn t3qu(p f64) f64 {
	return t3_poly(t3qu_i, t3qu_n, p)
}

pub fn t3rx(p f64) f64 {
	return t3_poly(t3rx_i, t3rx_n, p)
}

pub fn t3uv(p f64) f64 {
	return t3_poly(t3uv_i, t3uv_n, p)
}

pub fn t3wx(p f64) f64 {
	return t3_log(t3wx_i, t3wx_n, p)
}

// Which side of the saturation line a state at T = Tsat(p) is on: the
// liquid's, the vapour's, or as the release's tables have it by T.
enum Side {
	by_t
	liquid
	vapour
}

// subregion3 is SR5-05's Table 2 and, between T3qu and T3rx near the
// critical point, its Table 10: the subregion of region 3 holding (p, T).
fn subregion3(p f64, t f64, side Side) u8 {
	if p > 40e6 {
		return if t <= t3ab(p) { `a` } else { `b` }
	}
	if p > 25e6 {
		return if t <= t3cd(p) {
			`c`
		} else if t <= t3ab(p) {
			`d`
		} else if t <= t3ef(p) {
			`e`
		} else {
			`f`
		}
	}
	if p > 23.5e6 {
		return if t <= t3cd(p) {
			`c`
		} else if t <= t3gh(p) {
			`g`
		} else if t <= t3ef(p) {
			`h`
		} else if t <= t3ij(p) {
			`i`
		} else if t <= t3jk(p) {
			`j`
		} else {
			`k`
		}
	}
	if p > 23e6 {
		return if t <= t3cd(p) {
			`c`
		} else if t <= t3gh(p) {
			`l`
		} else if t <= t3ef(p) {
			`h`
		} else if t <= t3ij(p) {
			`i`
		} else if t <= t3jk(p) {
			`j`
		} else {
			`k`
		}
	}
	if p > 22.5e6 {
		return if t <= t3cd(p) {
			`c`
		} else if t <= t3gh(p) {
			`l`
		} else if t <= t3mn(p) {
			`m`
		} else if t <= t3ef(p) {
			`n`
		} else if t <= t3op(p) {
			`o`
		} else if t <= t3ij(p) {
			`p`
		} else if t <= t3jk(p) {
			`j`
		} else {
			`k`
		}
	}
	if p > psat97(643.15) {
		if t <= t3cd(p) {
			return `c`
		}
		if t <= t3qu(p) {
			return `q`
		}
		if t > t3rx(p) {
			return if t <= t3jk(p) { `r` } else { `k` }
		}
		return auxiliary3(p, t, side)
	}
	ts := tsat97(p)
	liquid := match side {
		.liquid { true }
		.vapour { false }
		.by_t {
			if p > p3cd { t < ts } else { t <= ts }
		}
	}
	if p > 20.5e6 {
		return if t <= t3cd(p) {
			`c`
		} else if liquid {
			`s`
		} else if t <= t3jk(p) {
			`r`
		} else {
			`k`
		}
	}
	if p > p3cd {
		return if t <= t3cd(p) {
			`c`
		} else if liquid {
			`s`
		} else {
			`t`
		}
	}
	return if liquid { `c` } else { `t` }
}

// auxiliary3 is Table 10: the subregions 3u to 3z, between T3qu and T3rx.
fn auxiliary3(p f64, t f64, side Side) u8 {
	if p > 22.11e6 {
		return if t <= t3uv(p) {
			`u`
		} else if t <= t3ef(p) {
			`v`
		} else if t <= t3wx(p) {
			`w`
		} else {
			`x`
		}
	}
	if p > pc {
		return if t <= t3uv(p) {
			`u`
		} else if t <= t3ef(p) {
			`y`
		} else if t <= t3wx(p) {
			`z`
		} else {
			`x`
		}
	}
	ts := tsat97(p)
	liquid := match side {
		.liquid { true }
		.vapour { false }
		.by_t { t <= ts }
	}
	if liquid {
		if p > p3_v264 {
			return if t <= t3uv(p) { `u` } else { `y` }
		}
		return `u`
	}
	if p > p3_v385 {
		return if t <= t3wx(p) { `z` } else { `x` }
	}
	return `x`
}

// v3_backward is SR5-05's v(p, T) in region 3, as published: the backward
// equation of the subregion holding (p, T).
pub fn v3_backward(p f64, t f64) f64 {
	return v3_sub(sub_of(subregion3(p, t, .by_t)), p, t)
}
