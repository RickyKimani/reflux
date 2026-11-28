module transport

import math

// A liquid (D6, D7, D8): the Grunberg–Nissan rule with Lucas's (1981)
// compressed-liquid factor, Li's rule with Missenard's (1970), and
// Winterfeld, Scriven and Davis' surface tension. Each mixing rule here
// takes compounds that all have their values; the phase's own rules (D9)
// choose which.

// lucas_liquid is Lucas's compressed-liquid factor μ(P)/μ(Pˢᵃᵗ) for a
// compound with critical constants tc, pc, its acentric factor and its
// vapour pressure at T, as chemicals' Lucas has it,
//
//     μ/μˢᵃᵗ = [1 + D(ΔPr/2.118)^A] / (1 + C·ω·ΔPr),  ΔPr = (P − Pˢᵃᵗ)/Pc,
//
// with Tr taken as at most 1, and ΔPr below 0, a liquid below its vapour
// pressure, as 0: 1 there, where DWSIM takes |P − Pˢᵃᵗ|.
pub fn lucas_liquid(t f64, p f64, tc f64, pc f64, omega f64, psat f64) f64 {
	tr := math.min(t / tc, 1.0)
	c := tr * (tr * (tr * (tr * (tr * (tr * (15.6719 * tr - 59.8127) + 96.1209) - 84.8291) + 44.1706) - 13.404) +
		2.1616) - 0.07921
	d := 0.3257 * math.pow(1.0039 - math.pow(tr, 2.573), -0.2906) - 0.2086
	a := 0.9991 - 4.674e-4 / (1.0523 * math.pow(tr, -0.03877) - 1.0513)
	mut dpr := (p - psat) / pc
	if dpr < 0.0 {
		dpr = 0.0
	}
	return (1.0 + d * math.pow(dpr * (1.0 / 2.118), a)) / (1.0 + c * omega * dpr)
}

// Missenard's Q, by Tr (0.5 to 0.8) and Pr (1 to 200), as chemicals' table
// has it (from Poling et al., §10-10).
const missenard_tr = [0.5, 0.6, 0.7, 0.8]!
const missenard_pr = [1.0, 5.0, 10.0, 50.0, 100.0, 200.0]!
const missenard_q = [
	[0.012, 0.0165, 0.017, 0.019, 0.020, 0.020]!, // Tr 0.5, by Pr
	[0.015, 0.020, 0.022, 0.024, 0.025, 0.025]!,
	[0.018, 0.025, 0.027, 0.031, 0.032, 0.032]!,
	[0.036, 0.038, 0.038, 0.038, 0.038, 0.038]!,
]!

// bracket is the index i of the interval [xs[i], xs[i+1]] holding x, and x's
// fraction along it, x clamped to the table's ends.
fn bracket(xs []f64, x f64) (int, f64) {
	n := xs.len
	if !(x > xs[0]) {
		return 0, 0.0
	}
	if !(x < xs[n - 1]) {
		return n - 2, 1.0
	}
	mut i := 0
	for i < n - 2 && x >= xs[i + 1] {
		i++
	}
	return i, (x - xs[i]) / (xs[i + 1] - xs[i])
}

// missenard is Missenard's compressed-liquid factor λ(P)/λ for a compound
// with critical constants tc and pc, 1 + Q·Pr^0.7, Q bilinear in Tr and Pr
// on Missenard's table and held at its edges beyond them, as chemicals'
// spline of degree 1 does. It's 1 at zero pressure.
pub fn missenard(t f64, p f64, tc f64, pc f64) f64 {
	tr := t / tc
	pr := p / pc
	i, u := bracket(missenard_tr[..], tr)
	j, w := bracket(missenard_pr[..], pr)
	q0 := missenard_q[i][j] * (1.0 - w) + missenard_q[i][j + 1] * w
	q1 := missenard_q[i + 1][j] * (1.0 - w) + missenard_q[i + 1][j + 1] * w
	q := q0 * (1.0 - u) + q1 * u
	return 1.0 + q * math.pow(pr, 0.7)
}

// Parameter is a pair's Grunberg–Nissan Gᵢⱼ at T (D6): an expression in T
// the file gives (flowsheet's), or a constant.
pub interface Parameter {
	value(t f64) !f64
}

// Constant is a Parameter that doesn't move with T.
pub struct Constant {
pub:
	g f64
}

pub fn (c Constant) value(_ f64) !f64 {
	return c.g
}

// grunberg_nissan is a liquid's viscosity, Pa·s, by Grunberg and Nissan's
// rule (D6) over compounds at mole fractions x (summing to 1),
//
//     ln μ = Σᵢ xᵢ ln μᵢ + ½ΣᵢΣⱼ xᵢxⱼGᵢⱼ,
//
// g n×n, symmetric with a zero diagonal: with every Gᵢⱼ = 0, the
// logarithmic rule.
pub fn grunberg_nissan(x []f64, mu []f64, g []f64) f64 {
	n := x.len
	mut s := 0.0
	for i in 0 .. n {
		if x[i] == 0 {
			continue
		}
		s += x[i] * math.log(mu[i])
	}
	mut e := 0.0
	for i in 0 .. n {
		for j in 0 .. n {
			e += x[i] * x[j] * g[i * n + j]
		}
	}
	return math.exp(s + 0.5 * e)
}

// li is a liquid's thermal conductivity, W/(m·K), by Li's rule (DIPPR 9I,
// D7) over compounds at mole fractions x, each with its molar volume (m³/mol)
// and conductivity:
//
//     φᵢ = xᵢVᵢ/ΣxⱼVⱼ,  λ = ΣᵢΣⱼ φᵢφⱼ·2/(1/λᵢ + 1/λⱼ).
pub fn li(x []f64, v []f64, lambda []f64) f64 {
	n := x.len
	mut total := 0.0
	for i in 0 .. n {
		total += x[i] * v[i]
	}
	mut out := 0.0
	for i in 0 .. n {
		if x[i] == 0 {
			continue
		}
		pi := x[i] * v[i] / total
		for j in 0 .. n {
			if x[j] == 0 {
				continue
			}
			pj := x[j] * v[j] / total
			out += pi * pj * 2.0 / (1.0 / lambda[i] + 1.0 / lambda[j])
		}
	}
	return out
}

// winterfeld is a liquid's surface tension, N/m, by Winterfeld, Scriven and
// Davis (D8) over compounds at mole fractions x, each with its molar volume
// (m³/mol) and surface tension:
//
//     σ = ΣᵢΣⱼ xᵢxⱼ(σᵢσⱼ)^½VᵢVⱼ / (ΣₖxₖVₖ)².
pub fn winterfeld(x []f64, v []f64, sigma []f64) f64 {
	n := x.len
	mut total := 0.0
	for i in 0 .. n {
		total += x[i] * v[i]
	}
	mut out := 0.0
	for i in 0 .. n {
		if x[i] == 0 {
			continue
		}
		for j in 0 .. n {
			if x[j] == 0 {
				continue
			}
			out += x[i] * x[j] * math.sqrt(sigma[i] * sigma[j]) * v[i] * v[j]
		}
	}
	return out / (total * total)
}
