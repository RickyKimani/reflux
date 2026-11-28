module transport

import math
import quantity

// A gas at pressure (D5): Lucas's factor on the low-pressure viscosity, and
// Stiel and Thodos' residual on the low-pressure conductivity. The formulas
// are as Poling, Prausnitz and O'Connell give them (The Properties of Gases
// and Liquids, 5th ed., 2001, §9-4, §9-6 and §10-5), Lucas's as polykin and
// Cantera have them too; docs/design.md, "Transport", writes them out.

const r_gas = quantity.molar_gas_constant

// debye is 1 D in C·m, the unit of the component files' dipole moments.
pub const debye = 3.33564e-30

// lucas_fp0 is Lucas's polar factor Fp° for a compound at its own Tr, from
// its reduced dipole moment μr = 52.46 μ²Pc/Tc² (μ in debye, Pc in bar).
fn lucas_fp0(c Pure, tr f64) !f64 {
	if math.is_nan(c.dipole) {
		return error("Lucas's polar term needs ${c.name}'s dipole moment, which its component file hasn't")
	}
	mur := 52.46 * c.dipole * c.dipole * (c.pc * 1e-5) / (c.tc * c.tc)
	if mur < 0.022 {
		return 1.0
	}
	if !(c.zc < 0.292) {
		return error("Lucas's polar term needs Zc below 0.292, and ${c.name}'s is ${c.zc}")
	}
	f := 30.55 * math.pow(0.292 - c.zc, 1.72)
	if mur < 0.075 {
		return 1.0 + f
	}
	return 1.0 + f * math.abs(0.96 + 0.1 * (tr - 0.7))
}

// lucas_fq0 is Lucas's quantum factor Fq° for He, H₂ and D₂ (Q = 1.38, 0.76,
// 0.52), 1 for any other compound, M in g/mol.
fn lucas_fq0(c Pure, tr f64) f64 {
	if c.quantum == 0 {
		return 1.0
	}
	x := tr - 12.0
	sign := if x > 0 { 1.0 } else { -1.0 }
	return 1.22 * math.pow(c.quantum, 0.15) * (1.0 + 0.00385 * math.pow(x * x, 1.0 / (c.mw * 1e3)) * sign)
}

// lucas_z1 is Lucas's low-pressure Z₁ = η°ξ at Tr, with the polar and
// quantum factors given.
fn lucas_z1(tr f64, fp0 f64, fq0 f64) f64 {
	return (0.807 * math.pow(tr, 0.618) - 0.357 * math.exp(-0.449 * tr) + 0.340 * math.exp(-4.058 * tr) +
		0.018) * fp0 * fq0
}

// lucas_dilute is Lucas's own low-pressure viscosity of a compound, Pa·s,
// η° = Z₁/ξ with ξ = 0.176(Tc/(M³Pc⁴))^⅙ (Tc in K, M in g/mol, Pc in bar,
// η in μP). Transport takes η° from the correlations; this is for the test
// of Fp° and Fq° against chemicals' Lucas_gas.
pub fn lucas_dilute(c Pure, t f64) !f64 {
	tr := t / c.tc
	m := c.mw * 1e3
	pc := c.pc * 1e-5
	xi := 0.176 * math.pow(c.tc / (m * m * m * pc * pc * pc * pc), 1.0 / 6.0)
	return lucas_z1(tr, lucas_fp0(c, tr)!, lucas_fq0(c, tr)) / xi * 1e-7
}

// lucas_z2_low is Z₂ below Tr = 1, Lucas's form for a vapour below its
// saturation pressure. At Pr = 0 it's 0.6·Tr, from the same arithmetic.
fn lucas_z2_low(tr f64, pr f64) f64 {
	alpha := 3.262 + 14.98 * math.pow(pr, 5.508)
	beta := 1.390 + 5.746 * pr
	return 0.600 + 0.760 * math.pow(pr, alpha) + (6.990 * math.pow(pr, beta) - 0.6) * (1.0 - tr)
}

// LucasMixture is a gas's constants by Lucas's mixture rules (Poling et
// al., Eqs. 9-5.18 to 9-5.23): Tcm (K) and Mm (kg/mol) by mole fractions,
// Pcm = R·Tcm·ΣyZc/ΣyVc (Pa), and Fp°m and Fq°m the mole-fraction means of
// each compound's at its own Tr, Fq°m times A.
pub struct LucasMixture {
pub:
	tc  f64
	pc  f64
	mw  f64
	fp0 f64
	fq0 f64
}

// lucas_mixture is the constants of compounds cs at mole fractions y and
// T by Lucas's mixture rules; A = 1 − 0.01(M_H/M_L)^0.87 where the heaviest
// compound is more than nine times the lightest and 0.05 < y_H < 0.7, else 1.
pub fn lucas_mixture(cs []Pure, y []f64, t f64) !LucasMixture {
	mut tcm, mut zs, mut vs, mut mm := 0.0, 0.0, 0.0, 0.0
	mut fp0, mut fq0 := 0.0, 0.0
	mut heavy, mut light := -1, -1
	for i, c in cs {
		if y[i] == 0 {
			continue
		}
		tcm += y[i] * c.tc
		zs += y[i] * c.zc
		vs += y[i] * c.vc
		mm += y[i] * c.mw
		tri := t / c.tc
		fp0 += y[i] * lucas_fp0(c, tri)!
		fq0 += y[i] * lucas_fq0(c, tri)
		if heavy < 0 || c.mw > cs[heavy].mw {
			heavy = i
		}
		if light < 0 || c.mw < cs[light].mw {
			light = i
		}
	}
	if !(vs > 0 && zs > 0) {
		return error("Lucas's mixture rules need each compound's Zc and Vc")
	}
	ratio := cs[heavy].mw / cs[light].mw
	if ratio > 9 && y[heavy] > 0.05 && y[heavy] < 0.7 {
		fq0 *= 1.0 - 0.01 * math.pow(ratio, 0.87)
	}
	return LucasMixture{
		tc:  tcm
		pc:  r_gas * tcm * zs / vs
		mw:  mm
		fp0: fp0
		fq0: fq0
	}
}

// dilute is the mixture's low-pressure viscosity by Lucas's own method,
// Pa·s, η°m = Z₁(Trm)/ξm (Poling et al., Eq. 9-4.16 on the mixture's
// constants). Transport takes η° from the correlations by Wilke's rule
// (D4); this is for the test of the mixture rules against the book's
// Example 9-7.
pub fn (m LucasMixture) dilute(t f64) f64 {
	mw := m.mw * 1e3
	pc := m.pc * 1e-5
	xi := 0.176 * math.pow(m.tc / (mw * mw * mw * pc * pc * pc * pc), 1.0 / 6.0)
	return lucas_z1(t / m.tc, m.fp0, m.fq0) / xi * 1e-7
}

// lucas_factor is μ/μ° = Y·Fp·Fq, Lucas's high-pressure factor for a gas
// of compounds cs at mole fractions y (summing to 1), T and P, by Lucas's
// mixture rules: Tcm and Mm by mole fractions, Pcm = R·Tcm·ΣyZc/ΣyVc, Fp°m
// and Fq°m the mole-fraction means of each compound's at its own Tr, Fq°m
// times A = 1 − 0.01(M_H/M_L)^0.87 where the heaviest compound is more than
// nine times the lightest and 0.05 < y_H < 0.7.
//
//     Fp = [1 + (Fp° − 1)Y⁻³]/Fp°,  Fq = [1 + (Fq° − 1)(Y⁻¹ − 0.007(ln Y)⁴)]/Fq°
//
// Y is taken against its own value at zero pressure (D5): above Tr = 1
// that's Lucas's own, Z₂/Z₁; below it Z₂(P)/Z₂(0), since Lucas's Z₂ there
// tends to 0.6·Tr, not Z₁, and the published factor isn't 1 at zero
// pressure. So the factor is exactly 1 at P = 0. Below Tr = 1 the form
// holds for a vapour below its saturation pressure, so a gas above its
// pseudo-critical pressure there is refused; so is one past Lucas's range,
// Tr to 40 and Pr to 100.
pub fn lucas_factor(cs []Pure, y []f64, t f64, p f64) !f64 {
	m := lucas_mixture(cs, y, t)!
	tcm, pcm, fp0, fq0 := m.tc, m.pc, m.fp0, m.fq0
	tr := t / tcm
	pr := p / pcm
	mut yy := 0.0
	if tr <= 1.0 {
		if !(pr < 1.0) {
			return error("Lucas's factor below its pseudo-critical temperature (${shown(tcm)} K) is for a vapour below its saturation pressure, and ${shown(p)} Pa is above its pseudo-critical pressure, ${shown(pcm)} Pa")
		}
		yy = lucas_z2_low(tr, pr) / lucas_z2_low(tr, 0.0)
	} else {
		if !(tr < 40.0 && pr <= 100.0) {
			return error("Lucas's factor is for Tr below 40 and Pr to 100, not ${shown(tr)} and ${shown(pr)}")
		}
		z1 := lucas_z1(tr, fp0, fq0)
		a := 1.245e-3 / tr * math.exp(5.1726 * math.pow(tr, -0.3286))
		b := a * (1.6553 * tr - 1.2723)
		c := 0.4489 / tr * math.exp(3.0578 * math.pow(tr, -37.7332))
		d := 1.7368 / tr * math.exp(2.2310 * math.pow(tr, -7.6351))
		f := 0.9425 * math.exp(-0.1853 * math.pow(tr, 0.4489))
		z2 := z1 * (1.0 + a * math.pow(pr, 1.3088) / (b * math.pow(pr, f) + 1.0 / (1.0 + c * math.pow(pr,
			d))))
		yy = z2 / z1
	}
	// Fp and Fq as 1 + (F° − 1)(…)/F°, the published forms rearranged, so
	// that at Y = 1 each is 1 to the bit and the factor is exactly 1.
	fp := 1.0 + (fp0 - 1.0) * (math.pow(yy, -3.0) - 1.0) / fp0
	ln := math.log(yy)
	fq := 1.0 + (fq0 - 1.0) * (1.0 / yy - 0.007 * ln * ln * ln * ln - 1.0) / fq0
	return yy * fp * fq
}

// Pseudo is a gas's pseudo-critical constants by Yorizane's rules (Poling
// et al., §10-7), for Stiel and Thodos' residual:
//
//     Vcᵢⱼ = (Vcᵢ^⅓ + Vcⱼ^⅓)³/8,  Vcm = ΣΣ yᵢyⱼVcᵢⱼ,
//     Tcm = ΣΣ yᵢyⱼVcᵢⱼ(TcᵢTcⱼ)^½ / Vcm,
//     ωm = Σ yᵢωᵢ,  Zcm = 0.291 − 0.08ωm,  Pcm = Zcm·R·Tcm/Vcm,  Mm = Σ yᵢMᵢ.
pub struct Pseudo {
pub:
	tc f64 // K
	pc f64 // Pa
	vc f64 // m³/mol
	zc f64
	mw f64 // kg/mol
}

// yorizane is the pseudo-critical constants of compounds cs at y.
pub fn yorizane(cs []Pure, y []f64) !Pseudo {
	n := cs.len
	mut vcm, mut tv, mut wm, mut mm := 0.0, 0.0, 0.0, 0.0
	for i in 0 .. n {
		if y[i] == 0 {
			continue
		}
		if math.is_nan(cs[i].vc) {
			return error("Stiel and Thodos' residual needs ${cs[i].name}'s critical volume, which its component file hasn't")
		}
		wm += y[i] * cs[i].omega
		mm += y[i] * cs[i].mw
		for j in 0 .. n {
			if y[j] == 0 {
				continue
			}
			s := math.cbrt(cs[i].vc) + math.cbrt(cs[j].vc)
			vij := s * s * s / 8.0
			vcm += y[i] * y[j] * vij
			tv += y[i] * y[j] * vij * math.sqrt(cs[i].tc * cs[j].tc)
		}
	}
	tcm := tv / vcm
	zcm := 0.291 - 0.08 * wm
	return Pseudo{
		tc: tcm
		pc: zcm * r_gas * tcm / vcm
		vc: vcm
		zc: zcm
		mw: mm
	}
}

// stiel_thodos is λ − λ°, W/(m·K), Stiel and Thodos' (1964) residual
// conductivity of a gas at molar volume v (m³/mol) with pseudo-critical
// constants pc, in its three published pieces in ρr = Vc/V,
//
//     (λ − λ°)ΓZc⁵ = 1.22e-2(e^(0.535ρr) − 1)      ρr < 0.5
//                  = 1.14e-2(e^(0.67ρr) − 1.069)    ρr < 2
//                  = 2.60e-3(e^(1.155ρr) + 2.016)   ρr ≤ 2.8
//
// Γ = 210(Tc·M³/Pc⁴)^⅙ (M in g/mol, Pc in bar). It's 0 at zero density;
// past ρr = 2.8, its range, it's refused.
pub fn stiel_thodos(pc Pseudo, v f64) !f64 {
	rr := pc.vc / v
	if !(rr >= 0 && rr <= 2.8) {
		return error("Stiel and Thodos' residual is for ρr up to 2.8, and the gas's is ${shown(rr)}")
	}
	m := pc.mw * 1e3
	gamma := 210.0 * math.pow(pc.tc * m * m * m * math.pow(pc.pc * 1e-5, -4.0), 1.0 / 6.0)
	term := if rr < 0.5 {
		1.22e-2 * (math.exp(0.535 * rr) - 1.0)
	} else if rr < 2.0 {
		1.14e-2 * (math.exp(0.67 * rr) - 1.069)
	} else {
		2.60e-3 * (math.exp(1.155 * rr) + 2.016)
	}
	return term * math.pow(pc.zc, -5.0) / gamma
}
