module reactors

import expr
import kinetics
import math
import thermo

// Source is a reactor's rates as a source of each species, mol/(m³·s):
//
//     Sᵢ = Σⱼ sⱼ νᵢⱼ rⱼ(T, P, y),
//
// over the rates the reactor runs, sⱼ the catalyst loading (kg/m³) for a
// rate per mass of catalyst and 1 for one per volume. With its derivatives
// by the flows, T and P, from the rates' own by the chain rule (D5), φ and Z
// held fixed. Its buffers are made once, and evaluating it allocates
// nothing but the phase's volume root.
struct Source {
	k     kinetics.Kinetics
	rates []int         // the kinetics' rates the reactor runs
	scale []f64         // sⱼ, for each of those
	nu    [][]f64       // νᵢⱼ over the species, for each of those
	model ?thermo.Model // the phase's, for φ and Z; none for the ideal gas
mut:
	pt     kinetics.Point
	buf    expr.Buffer
	out    []f64 // every rate
	g      []f64 // the tape's gradients
	d      []f64 // every rate's derivatives by T, P and each y
	y      []f64
	phi    []f64
	z      f64
	volume f64          // the phase's volume, m³ per second of flow: m³/s
	st     thermo.State // the phase's, from the model, at the state last read
}

// per is a rate law's basis: per mass of catalyst or per volume.
fn per(law kinetics.RateLaw) kinetics.Per {
	return match law {
		kinetics.PowerLaw, kinetics.Lhhw, kinetics.EleyRideal, kinetics.Expression { law.per }
	}
}

// new_source runs the rates of the reactions given, by index into the
// kinetics' reactions: each must have a rate. A rate per mass needs a
// catalyst loading above zero.
fn new_source(unit string, k kinetics.Kinetics, reactions []int, loading f64, model ?thermo.Model) !Source {
	n := k.species.len
	mut rates, mut scale := []int{}, []f64{}
	mut nu := [][]f64{}
	for i in reactions {
		j := k.rated.index(i)
		if j < 0 {
			return spec_error(unit, '${k.reactions[i].name} has no rate law')
		}
		s := if per(k.laws[j]) == .mass { loading } else { 1.0 }
		if !(s > 0) || math.is_inf(s, 1) {
			return spec_error(unit, "${k.reactions[i].name}'s rate is per mass of catalyst, and the catalyst loading is ${loading} kg/m³")
		}
		mut row := []f64{len: n}
		r := k.reactions[i]
		for q, sp in r.species {
			row[sp] = r.nu[q]
		}
		rates << j
		scale << s
		nu << row
	}
	return Source{
		k:     k
		rates: rates
		scale: scale
		nu:    nu
		model: model
		pt:    k.point()
		buf:   k.buffer()
		out:   []f64{len: k.rates.len}
		g:     []f64{len: k.rates.len * k.tape.inputs}
		d:     []f64{len: k.rates.len * (2 + n)}
		y:     []f64{len: n}
		phi:   []f64{len: n, init: 1.0}
		z:     1.0
	}
}

// state reads the phase at T, P and flows: its mole fractions, a flow that
// rounding takes below zero read as zero, as DWSIM reads it, and its φ, Z
// and volume from the model.
fn (mut s Source) state(t f64, p f64, flows []f64) !kinetics.State {
	mut total := 0.0
	for f in flows {
		total += math.max(f, 0)
	}
	if !(total > 0) || !math.is_finite(total) {
		return error('no flow')
	}
	for i, f in flows {
		s.y[i] = math.max(f, 0) / total
	}
	if model := s.model {
		n := flows.map(math.max(it, 0))
		st := thermo.solve_tp(model, t, p, n, .vapor)!
		for i, l in st.ln_phi {
			s.phi[i] = math.exp(l)
		}
		s.z = st.z
		s.volume = st.v
		s.st = st
	} else {
		s.volume = total * thermo.r_gas * t / p
	}
	return kinetics.State{
		t:   t
		p:   p
		y:   s.y
		phi: s.phi
		z:   s.z
	}
}

// at writes each species' source into src, mol/(m³·s), at a state just read.
fn (mut s Source) at(st kinetics.State, mut src []f64) ! {
	s.k.fill(st, mut s.pt)!
	s.k.rates(&s.pt, mut s.buf, mut s.out)
	for i in 0 .. src.len {
		src[i] = 0
	}
	for q, j in s.rates {
		r := s.scale[q] * s.out[j]
		if !math.is_finite(r) {
			return error("${s.k.reactions[s.k.rated[j]].name}'s rate isn't finite")
		}
		for i, v in s.nu[q] {
			if v != 0 {
				src[i] += v * r
			}
		}
	}
}

// derivs writes the sources' derivatives, after at at the same state: by
// each flow into dflows (n×n row-major, ∂Sᵢ/∂Fₖ), through yₗ = Fₗ/ΣF; and by
// T and P.
fn (mut s Source) derivs(st kinetics.State, total f64, mut dflows []f64, mut dt []f64, mut dp []f64) {
	n := s.y.len
	s.k.jacobian(st, &s.pt, mut s.buf, mut s.g, mut s.d)
	w := 2 + n
	for i in 0 .. n * n {
		dflows[i] = 0
	}
	for i in 0 .. n {
		dt[i] = 0
		dp[i] = 0
	}
	for q, j in s.rates {
		row := j * w
		for i, v in s.nu[q] {
			if v == 0 {
				continue
			}
			c := v * s.scale[q]
			dt[i] += c * s.d[row]
			dp[i] += c * s.d[row + 1]
			// ∂r/∂Fₖ = Σₗ ∂r/∂yₗ·(δₗₖ − yₗ)/ΣF.
			mut mean := 0.0
			for l in 0 .. n {
				mean += s.d[row + 2 + l] * s.y[l]
			}
			for kk in 0 .. n {
				dflows[i * n + kk] += c * (s.d[row + 2 + kk] - mean) / total
			}
		}
	}
}
