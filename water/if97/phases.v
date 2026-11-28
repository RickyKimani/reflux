module if97

// The states a flash asks of IF97 as a phase model (M5b's D12): each side's
// equation at (p, T) across the saturation line, as far as it has a state,
// the saturated phases at a (p, T) of the line, and region 2's ideal-gas
// part, which the phase model's residuals are against (D14).

// ideal2_j and ideal2_n are region 2's ideal-gas part's exponents and
// coefficients, R7-97's Table 10.
pub const ideal2_j = r2o_j
pub const ideal2_n = r2o_n

// region2_ideal is region 2's ideal-gas part alone at (p, T), Eq. (16):
// γ° = ln π + Σ nᵢ° τ^Jᵢ°, by Table 12's relations with no residual.
pub fn region2_ideal(p f64, t f64) State {
	pi := p / 1e6
	tau := 540.0 / t
	o := gamma_ideal(tau, r2o_j, r2o_n, r2o_n[0], r2o_n[1])
	return gas(2, p, t, pi, tau, o, Gibbs{})
}

fn in_range(p f64, t f64) ! {
	if !(t >= t_min && t <= t_max) || !(p > 0.0 && p <= p_max) || (t > t5 && p > p5_max) {
		return out_of_range(range_text(p, '${t} K'))
	}
}

// liquid_pt is the liquid at (p, T): region 1 to 623.15 K, region 3's
// liquid side to Tc, at a p where the vapour is stable too, as far as its
// equation goes; at and above Tc, and at or above pc, the one state there.
pub fn liquid_pt(p f64, t f64) !State {
	in_range(p, t)!
	if t <= t13 {
		return region1(p, t)
	}
	if t < tc && p < pc {
		return at_p(region3(rho3(p, t, .liquid)!, t), p)
	}
	return state_pt(p, t)
}

// vapour_pt is the vapour at (p, T): region 2 to 623.15 K, and above it to
// Tc region 2 up to B23 and region 3's vapour side beyond, at a p where the
// liquid is stable too, as far as its equation goes; at and above Tc, and at
// or above pc, the one state there.
pub fn vapour_pt(p f64, t f64) !State {
	in_range(p, t)!
	if t <= t13 {
		return region2(p, t)
	}
	if t < tc && p < pc {
		if p <= p23(t) {
			return region2(p, t)
		}
		return at_p(region3(rho3(p, t, .vapour)!, t), p)
	}
	return state_pt(p, t)
}

// saturated_at is the liquid and the vapour at (p, T) on the saturation
// line, below the critical point: regions 1 and 2 to 623.15 K, region 3's
// two sides above it.
pub fn saturated_at(p f64, t f64) !(State, State) {
	if !(t >= t_min && t < tc) || !(p > 0.0 && p < pc) {
		return out_of_range('IF97 has saturated phases from ${t_min} K to its critical point, ${tc} K: not ${t} K and ${p} Pa')
	}
	return saturated(p, t)
}
