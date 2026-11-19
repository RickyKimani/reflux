module kinetics

import expr
import math

// State is where rates are evaluated: T, P, each species' mole fraction,
// and, for a fugacity or concentration basis, each fugacity coefficient
// φᵢ and the compressibility Z (1 for the ideal gas, the default).
pub struct State {
pub:
	t   f64   @[required] // K
	p   f64   @[required] // Pa
	y   []f64 @[required]
	phi []f64 // 1 for each species if empty
	z   f64 = 1.0
}

// Point is a state's inputs to the tape, with what the chain rule needs:
// each reaction's d ln K/dT.
pub struct Point {
pub mut:
	x        []f64
	dln_k_dt []f64
	// kept says that its K and d ln K/dT are in, at k_t, where fill keeps
	// them until T moves (M8's D11).
	kept bool
	k_t  f64
	// Each species' ideal-gas H° and S° at k_t, and whether it's in: the
	// reactions that share a species share its two integrals.
	h    []f64
	s    []f64
	have []bool
	// integrals counts the species whose H° and S° a fill took, each two
	// integrals of its Cp, since the point was made.
	integrals int
}

// point makes a point's workspace, once, for fill to fill at each state.
pub fn (k &Kinetics) point() Point {
	return Point{
		x:        []f64{len: k.ctx.slots()}
		dln_k_dt: []f64{len: k.reactions.len}
		h:        []f64{len: k.species.len}
		s:        []f64{len: k.species.len}
		have:     []bool{len: k.species.len}
	}
}

// fill fills a point with a state's inputs: T, P, then each species'
// partial pressure pᵢ = yᵢP, fugacity fᵢ = φᵢyᵢP, concentration
// cᵢ = yᵢP/(ZRT) and mole fraction yᵢ, then K at 1 bar of each reaction
// some rate uses, from its standard properties at T (D7).
//
// K is computed only when T has moved since the point's last fill (M8's
// D11). Each species' Cp is integrated from 298.15 K for it, by quadrature
// for most compounds, which was 93% of a reforming PFR's solve (M9 step 11),
// though an isothermal reactor's T never moves. The K kept is the K a new
// point would compute, to the bit.
//
// At a new T, each species' two integrals are taken once, however many
// reactions it's in, and each reaction sums them as Reaction.standard
// does: its K, to the bit, in 10 integrals for the reforming model's five
// species where its three reactions took 24.
pub fn (k &Kinetics) fill(s State, mut pt Point) ! {
	n := k.species.len
	if s.y.len != n || (s.phi.len != 0 && s.phi.len != n) {
		return refuse('a state', 'mole fractions and fugacity coefficients for ${n} species')
	}
	pt.x[0] = s.t
	pt.x[1] = s.p
	ct := s.p / (s.z * gas_constant * s.t)
	for i in 0 .. n {
		phi := if s.phi.len == n { s.phi[i] } else { 1.0 }
		pt.x[k.ctx.slot(0, i)] = s.y[i] * s.p
		pt.x[k.ctx.slot(1, i)] = phi * s.y[i] * s.p
		pt.x[k.ctx.slot(2, i)] = s.y[i] * ct
		pt.x[k.ctx.slot(3, i)] = s.y[i]
	}
	if pt.kept && s.t == pt.k_t {
		return
	}
	// Not kept until every K is in: a reaction that fails leaves the point
	// to be filled again.
	pt.kept = false
	if pt.have.len != n {
		// A point made by hand, without point's workspace.
		pt.h = []f64{len: n}
		pt.s = []f64{len: n}
		pt.have = []bool{len: n}
	}
	for i in 0 .. n {
		pt.have[i] = false
	}
	for r in k.k_used {
		rx := k.reactions[r]
		rx.has_standard(s.t)!
		mut dh := 0.0
		mut ds := 0.0
		for j, q in rx.species {
			if !pt.have[q] {
				h, e := rx.species_standard(k.species[q], j, s.t)!
				pt.h[q] = h
				pt.s[q] = e
				pt.have[q] = true
				pt.integrals++
			}
			dh += rx.nu[j] * pt.h[q]
			ds += rx.nu[j] * pt.s[q]
		}
		st := standard_of(s.t, dh, ds)
		pt.x[k.ctx.slot(4, r)] = math.exp(st.ln_k)
		pt.dln_k_dt[r] = st.dln_k_dt
	}
	pt.k_t = s.t
	pt.kept = true
}

// buffer makes an evaluation's workspace.
pub fn (k &Kinetics) buffer() expr.Buffer {
	return k.tape.buffer()
}

// rates evaluates every rate at a point into out, in the order of rated,
// allocating nothing.
pub fn (k &Kinetics) rates(pt &Point, mut buf expr.Buffer, mut out []f64) {
	k.tape.eval(pt.x, mut buf)
	for i in 0 .. k.rates.len {
		out[i] = k.tape.output(buf, i)
	}
}

// jacobian is every rate's derivative with respect to T, P and each mole
// fraction, at a state whose rates were just evaluated into buf: d[i*(2+n)
// + 0] is dr_i/dT, + 1 dr_i/dP, and + 2 + j dr_i/dy_j. The chain rule runs
// through each family's dependence on the state and each K's on T,
// d ln K/dT = ΔH°/RT², with φ and Z held fixed: exact for the ideal gas,
// and a reactor that knows their derivatives adds their terms itself. g is
// the tape's gradient, inputs long for each rate.
pub fn (k &Kinetics) jacobian(s State, pt &Point, mut buf expr.Buffer, mut g []f64, mut d []f64) {
	k.tape.gradient(mut buf, mut g)
	n := k.species.len
	w := 2 + n
	m := k.tape.inputs
	ct := s.p / (s.z * gas_constant * s.t)
	for r in 0 .. k.rates.len {
		row := r * m
		mut dt := g[row]
		mut dp := g[row + 1]
		for i in 0 .. n {
			phi := if s.phi.len == n { s.phi[i] } else { 1.0 }
			gp := g[row + k.ctx.slot(0, i)]
			gf := g[row + k.ctx.slot(1, i)]
			gc := g[row + k.ctx.slot(2, i)]
			gy := g[row + k.ctx.slot(3, i)]
			ci := pt.x[k.ctx.slot(2, i)]
			dt -= gc * ci / s.t
			dp += (gp + gf * phi) * s.y[i] + gc * s.y[i] * ct / s.p
			d[r * w + 2 + i] = (gp + gf * phi) * s.p + gc * ct + gy
		}
		for q in k.k_used {
			slot := k.ctx.slot(4, q)
			dt += g[row + slot] * pt.x[slot] * pt.dln_k_dt[q]
		}
		d[r * w] = dt
		d[r * w + 1] = dp
	}
}
