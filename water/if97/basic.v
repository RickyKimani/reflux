module if97

import math

// The basic equations, as published (R7-97): each gives a state from its
// own variables with no range check. The steam-table functions
// (functions.v) choose the region and refuse what's outside it.

// A dimensionless Gibbs energy γ(π, τ) and its derivatives.
struct Gibbs {
mut:
	g   f64
	gp  f64
	gpp f64
	gt  f64
	gtt f64
	gpt f64
}

// A dimensionless Helmholtz energy φ(δ, τ) and its derivatives.
struct Helmholtz {
mut:
	f   f64
	fd  f64
	fdd f64
	ft  f64
	ftt f64
	fdt f64
}

// gamma1 is region 1's γ = Σ nᵢ (7.1 − π)^Iᵢ (τ − 1.222)^Jᵢ, Eq. (7), with
// Table 4's derivatives.
fn gamma1(pi f64, tau f64) Gibbs {
	x := 7.1 - pi
	y := tau - 1.222
	mut g := Gibbs{}
	for k in 0 .. r1_n.len {
		i := f64(r1_i[k])
		j := f64(r1_j[k])
		a := r1_n[k] * ipow(x, r1_i[k]) * ipow(y, r1_j[k])
		g.g += a
		g.gp -= i * a / x
		g.gpp += i * (i - 1.0) * a / (x * x)
		g.gt += j * a / y
		g.gtt += j * (j - 1.0) * a / (y * y)
		g.gpt -= i * j * a / (x * y)
	}
	return g
}

// region1 is the liquid's basic equation, Eq. (7), with Table 3's
// relations: p* = 16.53 MPa, T* = 1386 K.
pub fn region1(p f64, t f64) State {
	pi := p / 16.53e6
	tau := 1386.0 / t
	g := gamma1(pi, tau)
	rt := r * t
	d := g.gp - tau * g.gpt
	return State{
		region: 1
		t:      t
		p:      p
		v:      pi * g.gp * rt / p
		h:      rt * tau * g.gt
		s:      r * (tau * g.gt - g.g)
		u:      rt * (tau * g.gt - pi * g.gp)
		cp:     -r * tau * tau * g.gtt
		cv:     r * (-tau * tau * g.gtt + d * d / g.gpp)
		w:      math.sqrt(rt * g.gp * g.gp / (d * d / (tau * tau * g.gtt) - g.gpp))
		// v = (RT/p*)γ_π, so (∂v/∂p)_T = (RT/p*²)γ_ππ
		drhodp: -rt * g.gpp / (16.53e6 * 16.53e6 * (pi * g.gp * rt / p) * (pi * g.gp * rt / p))
		// and (∂v/∂T)_p = (R/p*)(γ_π − τγ_πτ)
		dvdt:   r / 16.53e6 * d
	}
}

// gamma_res is a residual part Σ nᵢ π^Iᵢ y^Jᵢ, y = τ − shift: region 2's
// (Eq. 17, shift 0.5), the metastable vapour's (Eq. 19) and region 5's
// (Eq. 34, shift 0), with their derivatives in π and τ.
fn gamma_res(pi f64, tau f64, shift f64, ii []int, jj []int, nn []f64) Gibbs {
	y := tau - shift
	mut g := Gibbs{}
	for k in 0 .. nn.len {
		i := f64(ii[k])
		j := f64(jj[k])
		a := nn[k] * ipow(pi, ii[k]) * ipow(y, jj[k])
		g.g += a
		g.gp += i * a / pi
		g.gpp += i * (i - 1.0) * a / (pi * pi)
		g.gt += j * a / y
		g.gtt += j * (j - 1.0) * a / (y * y)
		g.gpt += i * j * a / (pi * y)
	}
	return g
}

// An ideal-gas part γ° = ln π + Σ nᵢ° τ^Jᵢ°: its τ-sum and the sum's two
// derivatives. n01 and n02 replace n₁° and n₂° (the metastable vapour's).
struct Ideal {
	g   f64 // Σ nᵢ° τ^Jᵢ°, without ln π
	gt  f64
	gtt f64
}

fn gamma_ideal(tau f64, jj []int, nn []f64, n01 f64, n02 f64) Ideal {
	mut g, mut gt, mut gtt := 0.0, 0.0, 0.0
	for k in 0 .. nn.len {
		n := if k == 0 {
			n01
		} else if k == 1 {
			n02
		} else {
			nn[k]
		}
		j := f64(jj[k])
		a := n * ipow(tau, jj[k])
		g += a
		gt += j * a / tau
		gtt += j * (j - 1.0) * a / (tau * tau)
	}
	return Ideal{g, gt, gtt}
}

// gas is Table 12's (and Table 39's) relations for a Gibbs energy split
// into an ideal-gas part and a residual part, so that the low pressure's
// properties keep their digits.
fn gas(region int, p f64, t f64, pi f64, tau f64, o Ideal, g Gibbs) State {
	rt := r * t
	gt := o.gt + g.gt
	gtt := o.gtt + g.gtt
	a := 1.0 + pi * g.gp - tau * pi * g.gpt
	b := 1.0 - pi * pi * g.gpp
	v := rt / p * (1.0 + pi * g.gp)
	ps := p / pi // p*
	return State{
		region: region
		t:      t
		p:      p
		v:      v
		h:      rt * tau * gt
		s:      r * (tau * gt - (math.log(pi) + o.g + g.g))
		u:      rt * (tau * gt - 1.0 - pi * g.gp)
		cp:     -r * tau * tau * gtt
		cv:     r * (-tau * tau * gtt - a * a / b)
		w:      math.sqrt(rt * (1.0 + 2.0 * pi * g.gp + pi * pi * g.gp * g.gp) / (b +
			a * a / (tau * tau * gtt)))
		// γ_ππ = −1/π² + γʳ_ππ, and (∂v/∂p)_T = (RT/p*²)γ_ππ
		drhodp: -rt * (g.gpp - 1.0 / (pi * pi)) / (ps * ps * v * v)
		// v = (RT/p)(1 + πγʳ_π), so (∂v/∂T)_p = (R/p)(1 + πγʳ_π − τπγʳ_πτ)
		dvdt:   r / p * a
	}
}

// region2 is the vapour's basic equation, Eqs. (15) to (17): p* = 1 MPa,
// T* = 540 K.
pub fn region2(p f64, t f64) State {
	pi := p / 1e6
	tau := 540.0 / t
	o := gamma_ideal(tau, r2o_j, r2o_n, r2o_n[0], r2o_n[1])
	g := gamma_res(pi, tau, 0.5, r2_i, r2_j, r2_n)
	return gas(2, p, t, pi, tau, o, g)
}

// region2_metastable is the supplementary equation for the metastable
// vapour, Eqs. (18) and (19), valid from the saturated vapour to the 5%
// equilibrium-moisture line, up to 10 MPa (R7-97, Section 6.2). Its states
// carry region 2.
pub fn region2_metastable(p f64, t f64) State {
	pi := p / 1e6
	tau := 540.0 / t
	o := gamma_ideal(tau, r2o_j, r2o_n, r2m_n01, r2m_n02)
	g := gamma_res(pi, tau, 0.5, r2m_i, r2m_j, r2m_n)
	return gas(2, p, t, pi, tau, o, g)
}

// region5 is the high temperature's basic equation, Eqs. (32) to (34):
// p* = 1 MPa, T* = 1000 K.
pub fn region5(p f64, t f64) State {
	pi := p / 1e6
	tau := 1000.0 / t
	o := gamma_ideal(tau, r5o_j, r5o_n, r5o_n[0], r5o_n[1])
	g := gamma_res(pi, tau, 0.0, r5_i, r5_j, r5_n)
	return gas(5, p, t, pi, tau, o, g)
}

// phi3 is region 3's φ = n₁ ln δ + Σ nᵢ δ^Iᵢ τ^Jᵢ, Eq. (28), with Table 32's
// derivatives.
fn phi3(delta f64, tau f64) Helmholtz {
	mut f := Helmholtz{
		f:   r3_n1 * math.log(delta)
		fd:  r3_n1 / delta
		fdd: -r3_n1 / (delta * delta)
	}
	for k in 0 .. r3_n.len {
		i := f64(r3_i[k])
		j := f64(r3_j[k])
		a := r3_n[k] * ipow(delta, r3_i[k]) * ipow(tau, r3_j[k])
		f.f += a
		f.fd += i * a / delta
		f.fdd += i * (i - 1.0) * a / (delta * delta)
		f.ft += j * a / tau
		f.ftt += j * (j - 1.0) * a / (tau * tau)
		f.fdt += i * j * a / (delta * tau)
	}
	return f
}

// region3 is the near-critical region's basic equation, Eq. (28), in ρ and
// T, with Table 31's relations: ρ* = ρc, T* = Tc.
pub fn region3(rho f64, t f64) State {
	delta := rho / rhoc
	tau := tc / t
	f := phi3(delta, tau)
	rt := r * t
	a := delta * f.fd - delta * tau * f.fdt
	b := 2.0 * delta * f.fd + delta * delta * f.fdd
	return State{
		region: 3
		t:      t
		p:      rho * rt * delta * f.fd
		v:      1.0 / rho
		h:      rt * (tau * f.ft + delta * f.fd)
		s:      r * (tau * f.ft - f.f)
		u:      rt * tau * f.ft
		cp:     r * (-tau * tau * f.ftt + a * a / b)
		cv:     -r * tau * tau * f.ftt
		w:      math.sqrt(rt * (b - a * a / (tau * tau * f.ftt)))
		drhodp: 1.0 / (rt * b) // (∂p/∂ρ)_T = RT(2δφ_δ + δ²φ_δδ)
		// (∂v/∂T)_p = (∂p/∂T)_ρ/(ρ²(∂p/∂ρ)_T), (∂p/∂T)_ρ = ρR(δφ_δ − δτφ_δτ)
		dvdt:   a / (rho * t * b)
	}
}
