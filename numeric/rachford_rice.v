module numeric

import math

// RRSolution is the split that rachford_rice finds for fixed K-values: the
// vapour fraction, the liquid fraction, and the two phases' mole fractions.
pub struct RRSolution {
pub:
	// beta is the vapour fraction V/F. It lies outside [0, 1] for a negative
	// flash: the K-values then describe a single-phase feed.
	beta f64
	// liquid is L/F = 1 − β, computed directly rather than by subtracting,
	// so it keeps its digits when the liquid is a trace.
	liquid     f64
	x          []f64 // liquid mole fractions
	y          []f64 // vapour mole fractions, yᵢ = Kᵢxᵢ
	iterations int
}

// The most Newton-or-bisection steps rachford_rice takes. Newton needs
// about 5–10, and bisection alone would need at most about 1,100 (halving
// down to the smallest f64).
const rr_max_iter = 2500

// rachford_rice solves the Rachford–Rice equation
//
//     Σᵢ zᵢ(Kᵢ − 1) / (1 + β(Kᵢ − 1)) = 0
//
// for the vapour fraction β of a feed z split into a liquid x and a vapour
// y = Kx. z should sum to 1; x and y then sum to 1 too.
//
// The equation has a pole wherever 1 + β(Kᵢ − 1) = 0. The root that gives
// positive compositions lies between the two poles nearest zero:
// c₁ = 1/(1 − K_max) < 0 and cₙ = 1/(1 − K_min) > 1, counting only
// components present in z. Across that window the sum falls from +∞ to −∞,
// so there is exactly one root there. It may lie outside [0, 1], which is
// Whitson and Michelsen's negative flash. rachford_rice fails only when no
// such window exists, because every K is on one side of 1, or when its
// inputs aren't valid.
//
// To keep extreme K-values accurate, β is written as a position in the
// window, β = c₁ + t(cₙ − c₁), with s = 1 − t. Then
//
//     1 + β(Kᵢ − 1) = s·uᵢ + t·vᵢ,
//     uᵢ = (K_max − Kᵢ)/(K_max − 1) ≥ 0,  vᵢ = (Kᵢ − K_min)/(1 − K_min) ≥ 0.
//
// Neither term is ever negative, so each denominator has full relative
// precision anywhere in the window: a heavy component's 1e-20, say, as the
// liquid fraction goes to zero. Subtracting β from a pole would lose it.
// The solver stores whichever of t and s is smaller, so a root 1e-30 from
// one end of the window still has all its digits.
//
// The function solved is Leibovici and Neoschil's (1992): the sum times the
// distances to the two poles, Φ(t) = t·s·Σᵢ zᵢ(Kᵢ − 1)/(s·uᵢ + t·vᵢ). It
// stays finite at both ends and is close to linear near each, so Newton's
// method converges fast even when the root is near an end. It's safeguarded
// by bisection on a bracket that always contains the root.
pub fn rachford_rice(z []f64, k []f64) !RRSolution {
	n := z.len
	if n == 0 || k.len != n {
		return error('rachford_rice: z has ${n} entries and K has ${k.len}')
	}
	mut kmax, mut kmin := 0.0, math.inf(1)
	mut total := 0.0
	for i in 0 .. n {
		if !math.is_finite(z[i]) || z[i] < 0.0 {
			return error('rachford_rice: z[${i}] = ${z[i]} is not a mole fraction')
		}
		if !math.is_finite(k[i]) || k[i] <= 0.0 {
			return error('rachford_rice: K[${i}] = ${k[i]} is not a positive number')
		}
		if z[i] > 0.0 {
			kmax = math.max(kmax, k[i])
			kmin = math.min(kmin, k[i])
			total += z[i]
		}
	}
	if total == 0.0 {
		return error('rachford_rice: z is all zeros')
	}
	if kmax <= 1.0 || kmin >= 1.0 {
		side := if kmax <= 1.0 { 'at most 1: all liquid' } else { 'at least 1: all vapour' }
		return error('rachford_rice: every K is ${side}, so there is no two-phase split')
	}
	mut u := []f64{len: n}
	mut v := []f64{len: n}
	for i in 0 .. n {
		u[i] = (kmax - k[i]) / (kmax - 1.0)
		v[i] = (k[i] - kmin) / (1.0 - kmin)
	}
	// Φ at the two ends of the window: only the components that set K_max
	// (u = 0) survive at t = 0, and only those that set K_min (v = 0) at
	// t = 1. Interpolating between them gives the first guess.
	mut phi0, mut phi1 := 0.0, 0.0
	for i in 0 .. n {
		if z[i] > 0.0 {
			if u[i] == 0.0 {
				phi0 += z[i] * (k[i] - 1.0) / v[i]
			}
			if v[i] == 0.0 {
				phi1 += z[i] * (k[i] - 1.0) / u[i]
			}
		}
	}
	mut p := WindowPoint{
		t: phi0 / (phi0 - phi1)
		s: -phi1 / (phi0 - phi1)
	}
	// The bracket: Φ > 0 at lo and Φ < 0 at hi.
	mut lo := WindowPoint{0.0, 1.0}
	mut hi := WindowPoint{1.0, 0.0}
	mut iterations := 0
	for {
		iterations++
		if iterations > rr_max_iter {
			return error('rachford_rice: no convergence in ${rr_max_iter} iterations')
		}
		phi, dphi := rr_phi(z, k, u, v, p)
		if phi == 0.0 {
			break
		}
		if phi > 0.0 {
			lo = p
		} else {
			hi = p
		}
		mut next := p
		mut newton_ok := false
		if dphi < 0.0 {
			step := -phi / dphi
			// Converged when the Newton step would move the stored
			// coordinate by no more than a few rounding errors. Such a step
			// may not change p at all, so it's tested before the bracket.
			if math.abs(step) <= 4.0 * f64_epsilon * math.min(p.t, p.s) {
				break
			}
			next = p.moved(step)
			newton_ok = lo.before(next) && next.before(hi)
		}
		if !newton_ok {
			next = midpoint(lo, hi)
			if !lo.before(next) || !next.before(hi) {
				break // lo and hi are adjacent: p is as close as it gets
			}
		}
		p = next
	}
	mut x := []f64{len: n}
	mut y := []f64{len: n}
	for i in 0 .. n {
		if z[i] > 0.0 {
			x[i] = z[i] / (p.s * u[i] + p.t * v[i])
			y[i] = k[i] * x[i]
		}
	}
	// β from the nearer end of the window, so its distance from that end is
	// exact. width = cₙ − c₁, written to avoid cancellation.
	width := (kmax - kmin) / ((kmax - 1.0) * (1.0 - kmin))
	mut beta, mut liquid := 0.0, 0.0
	if p.t <= p.s {
		// c₁ = −1/(K_max − 1). It's tiny when K_max is huge, so a trace of
		// vapour keeps its digits. L isn't small on this side for any reason
		// the K-values preserve, so 1 − β loses nothing.
		beta = -1.0 / (kmax - 1.0) + width * p.t
		liquid = 1.0 - beta
	} else {
		// cₙ = 1/(1 − K_min), and 1 − cₙ = −K_min/(1 − K_min).
		beta = 1.0 / (1.0 - kmin) - width * p.s
		liquid = -kmin / (1.0 - kmin) + width * p.s
	}
	return RRSolution{
		beta:       beta
		liquid:     liquid
		x:          x
		y:          y
		iterations: iterations
	}
}

// rr_phi returns Φ(t) = t·s·Σᵢ zᵢ(Kᵢ − 1)/(s·uᵢ + t·vᵢ) and dΦ/dt at p.
// Each term is written with the ratio t/s, so no product or quotient
// overflows or underflows near either end of the window.
fn rr_phi(z []f64, k []f64, u []f64, v []f64, p WindowPoint) (f64, f64) {
	ratio := p.t / p.s
	mut phi, mut dphi := 0.0, 0.0
	for i in 0 .. z.len {
		if z[i] == 0.0 || k[i] == 1.0 {
			continue
		}
		// t·s/(s·u + t·v) = t/(u + v·t/s) = s/(u·s/t + v).
		qt := 1.0 / (u[i] + v[i] * ratio)
		qs := 1.0 / (u[i] / ratio + v[i])
		a := z[i] * (k[i] - 1.0)
		phi += a * (if p.t <= p.s { p.t * qt } else { p.s * qs })
		// d/dt [t·s/(s·u + t·v)] = (s²u − t²v)/(s·u + t·v)² = u·qt² − v·qs².
		dphi += a * (u[i] * qt * qt - v[i] * qs * qs)
	}
	return phi, dphi
}

// WindowPoint is a position in the Rachford–Rice window, held as its
// distances t and s = 1 − t from the two ends. The smaller one is stored
// exactly and the other is 1 minus it, so a point 1e-30 from an end keeps
// all its digits.
struct WindowPoint {
	t f64
	s f64
}

// moved returns the point d further along in t (and d less in s), keeping
// whichever distance is now the smaller one exact.
fn (p WindowPoint) moved(d f64) WindowPoint {
	t, s := p.t + d, p.s - d
	return if t <= s { WindowPoint{t, 1.0 - t} } else { WindowPoint{1.0 - s, s} }
}

// before reports whether p lies strictly nearer the t = 0 end than q. Near
// the far end it compares the exact s values.
fn (p WindowPoint) before(q WindowPoint) bool {
	if p.t <= 0.5 || q.t <= 0.5 {
		return p.t < q.t
	}
	return p.s > q.s
}

// midpoint splits a bracket halfway, in s when both ends store s exactly
// and in t otherwise. It's the fallback when a Newton step would leave the
// bracket, which is rare.
fn midpoint(a WindowPoint, b WindowPoint) WindowPoint {
	if a.s <= 0.5 && b.s <= 0.5 {
		s := 0.5 * (a.s + b.s)
		return WindowPoint{1.0 - s, s}
	}
	t := 0.5 * (a.t + b.t)
	return WindowPoint{t, 1.0 - t}
}
