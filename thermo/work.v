module thermo

import math
import numeric

// CubicWork is a cubic's solve_tp for many states, allocating nothing once
// made (M8 step 7): a reactor reads its phase some twenty times an
// integration step, and solve_tp's 944 bytes a call were nearly all a study
// allocated, which the collector's lock made threads wait on. It keeps each
// component's √a at the last T, one volume function it moves to each
// composition, and ln φ, and does solve_tp's arithmetic in solve_tp's
// order, so that it gives solve_tp's bits. Its volume function is asked
// only for F and its V-derivatives, so ∂A/∂T isn't computed: it's NaN, which
// would show in anything that asked. One per thread: it's mutable.
@[heap]
pub struct CubicWork {
	c Cubic
mut:
	t   f64   // the T that s is at; 0 before the first call
	s   []f64 // √aᵢ at t
	a_n []f64 // the volume function's Aᵢ, which it shares
	at  &CubicAtTn
	f   VolumeFunction
pub mut:
	ln_phi []f64 // ln φᵢ at the last state solved
}

// WorkState is a state CubicWork solved: the volume and Z, its ln φ being
// the work's.
pub struct WorkState {
pub:
	v f64 // total volume, m³
	z f64
}

// work makes a cubic's workspace for solve_tp.
pub fn (c Cubic) work() &CubicWork {
	nc := c.comps.len
	a_n := []f64{len: nc}
	at := &CubicAtTn{
		c:   c
		att: Attraction{
			a_n: a_n
		}
	}
	return &CubicWork{
		c:      c
		s:      []f64{len: nc}
		a_n:    a_n
		at:     at
		f:      VolumeFunction(at)
		ln_phi: []f64{len: nc}
	}
}

// solve_tp is thermo.solve_tp for this cubic: the requested volume root at
// (T, P, n), with ln φ in w.ln_phi.
pub fn (mut w CubicWork) solve_tp(t f64, p f64, n []f64, root Root) !WorkState {
	c := w.c
	nc := c.comps.len
	// check_tpn's checks, without boxing the cubic as a Model.
	if t <= 0.0 {
		return InputError{
			what:  'temperature (K)'
			value: t
		}
	}
	if p <= 0.0 {
		return InputError{
			what:  'pressure (Pa)'
			value: p
		}
	}
	if n.len != nc {
		return error('composition has ${n.len} entries but the model has ${nc} components')
	}
	if n.any(it < 0.0) || numeric.sum(n) <= 0.0 {
		return error('mole numbers must be non-negative with a positive total, got ${n}')
	}
	// sqrt_a's √aᵢ at T, kept until T moves.
	if t != w.t {
		for i in 0 .. nc {
			tc := c.comps[i].tc
			r, _, _ := c.alphas[i].root(t / tc)
			w.s[i] = math.sqrt(c.ac[i]) * r
		}
		w.t = t
	}
	// attraction's sums, into the volume function's own a_n.
	mut a := 0.0
	for i in 0 .. nc {
		mut q := 0.0
		for j in 0 .. nc {
			q += (1.0 - c.kij[i * nc + j]) * w.s[j] * n[j]
		}
		a += n[i] * w.s[i] * q
		w.a_n[i] = 2.0 * w.s[i] * q
	}
	// volume_function's fields, and fix's.
	rt := r_gas * t
	ntot := numeric.sum(n)
	bb := c.min_volume(n)
	w.at.t = t
	w.at.rt = rt
	w.at.ntot = ntot
	w.at.bb = bb
	w.at.att = Attraction{
		a:   a
		a_t: math.nan()
		a_n: w.a_n
	}
	w.at.a_rt = a / rt
	fx := Fixed{
		f:    w.f
		t:    t
		ntot: numeric.sum(n)
		vmin: c.min_volume(n)
	}
	v := match root {
		.liquid {
			liquid_volume(fx, p) or { vapor_volume(fx, p) or { return no_root(t, p) } }
		}
		.vapor {
			vapor_volume(fx, p) or { liquid_volume(fx, p) or { return no_root(t, p) } }
		}
		.stable {
			stable_volume(fx, p)!
		}
	}
	// state_at's Z and ln φ = ∂F/∂nᵢ − ln Z, residual's ∂F/∂nᵢ term by term.
	q := c.vb_terms(v, bb)
	z := p * v / (numeric.sum(n) * r_gas * t)
	ln_z := math.log(z)
	a_rt := w.at.a_rt
	for i in 0 .. nc {
		f_n := -q.g - ntot * q.g_b * c.b[i] - q.f * w.a_n[i] / rt - a_rt * q.f_b * c.b[i]
		w.ln_phi[i] = f_n - ln_z
	}
	return WorkState{
		v: v
		z: z
	}
}
