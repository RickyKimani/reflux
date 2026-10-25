module thermo

import math
import numeric

// Translated is a volume-shifted version of another model (Péneloux, Rauzy &
// Fréze, 1982). At volume V it behaves exactly as the inner model does at
// W = V + C, with C = Σ nᵢcᵢ:
//
//     P′(T, V, n) = P(T, V + C, n)
//     F′(T, V, n) = F(T, W, n) + n·ln(V/W)
//
// Every volume moves by −C, and each ln φᵢ by −cᵢP/RT, so vapour pressures
// and phase splits don't change. It exists to fix the cubics' liquid
// densities. A positive c shrinks volumes.
//
// It wraps any Model, including another Translated: a struct holding an
// interface that itself satisfies the interface. Nothing else in thermo knows
// it exists. The derivative formulas are in docs/design.md, "Volume shift".
pub struct Translated {
pub:
	inner Model
	c     []f64 // volume shift per component, m³/mol
}

// translate wraps a model with a constant volume shift per component.
pub fn translate(inner Model, c []f64) !Translated {
	comps := inner.components()
	nc := comps.len
	if c.len != nc {
		return error('got ${c.len} volume shifts for ${nc} components')
	}
	// Each pure component's packing limit must stay positive after the shift,
	// or the volume solver would be asked to start at V ≤ 0.
	for i in 0 .. nc {
		mut e := []f64{len: nc}
		e[i] = 1.0
		b := inner.min_volume(e)
		if c[i] >= b {
			return error('the volume shift for ${comps[i].name} (${c[i]} m³/mol) must be smaller than its covolume (${b} m³/mol)')
		}
	}
	return Translated{
		inner: inner
		c:     c.clone()
	}
}

// components implements Model.
pub fn (m Translated) components() []Component {
	return m.inner.components()
}

// min_volume implements Model: the inner limit, moved by −C.
pub fn (m Translated) min_volume(n []f64) f64 {
	return m.inner.min_volume(n) - numeric.dot(n, m.c)
}

// volume_derivs implements Model.
pub fn (m Translated) volume_derivs(t f64, v f64, n []f64) VolumeDerivs {
	ntot := numeric.sum(n)
	w := v + numeric.dot(n, m.c)
	d := m.inner.volume_derivs(t, w, n)
	return VolumeDerivs{
		f:    d.f + ntot * math.log(v / w)
		f_v:  d.f_v + ntot / v - ntot / w
		f_vv: d.f_vv - ntot / (v * v) + ntot / (w * w)
	}
}

// residual implements Model.
pub fn (m Translated) residual(t f64, v f64, n []f64) Residual {
	ntot := numeric.sum(n)
	w := v + numeric.dot(n, m.c)
	r := m.inner.residual(t, w, n)
	ln_vw := math.log(v / w)
	// ∂/∂nᵢ picks up the inner F_V through ∂W/∂nᵢ = cᵢ.
	f_n := []f64{len: n.len, init: r.f_n[index] + m.c[index] * r.f_v + ln_vw - ntot * m.c[index] / w}
	return Residual{
		f:    r.f + ntot * ln_vw
		f_v:  r.f_v + ntot / v - ntot / w
		f_vv: r.f_vv - ntot / (v * v) + ntot / (w * w)
		f_t:  r.f_t
		f_n:  f_n
	}
}

// residual2 implements Model.
pub fn (m Translated) residual2(t f64, v f64, n []f64) Residual2 {
	nc := n.len
	ntot := numeric.sum(n)
	w := v + numeric.dot(n, m.c)
	r := m.inner.residual2(t, w, n)
	ln_vw := math.log(v / w)
	mut f_n := []f64{len: nc}
	mut f_vn := []f64{len: nc}
	mut f_tn := []f64{len: nc}
	mut f_nn := []f64{len: nc * nc}
	for i in 0 .. nc {
		ci := m.c[i]
		f_n[i] = r.f_n[i] + ci * r.f_v + ln_vw - ntot * ci / w
		f_vn[i] = r.f_vn[i] + ci * r.f_vv + 1.0 / v - 1.0 / w + ntot * ci / (w * w)
		f_tn[i] = r.f_tn[i] + ci * r.f_tv
		for j in 0 .. nc {
			cj := m.c[j]
			f_nn[i * nc + j] = r.f_nn[i * nc + j] + cj * r.f_vn[i] + ci * r.f_vn[j] +
				ci * cj * r.f_vv - (ci + cj) / w + ntot * ci * cj / (w * w)
		}
	}
	return Residual2{
		Residual: Residual{
			f:    r.f + ntot * ln_vw
			f_v:  r.f_v + ntot / v - ntot / w
			f_vv: r.f_vv - ntot / (v * v) + ntot / (w * w)
			f_t:  r.f_t
			f_n:  f_n
		}
		f_tt:     r.f_tt
		f_tv:     r.f_tv
		f_vn:     f_vn
		f_tn:     f_tn
		f_nn:     f_nn
	}
}
