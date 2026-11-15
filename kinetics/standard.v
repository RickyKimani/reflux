module kinetics

import math
import thermo

// gas_constant is R = N_A·k_B, the 2019 SI's exact value, as expr's
// built-in constant R is: K(T) appears in rates that use that R.
pub const gas_constant = 8.31446261815324

// Standard is a reaction's standard properties at T, its species ideal
// gases at 1 bar (D7): ΔH°, ΔS° and ΔG°, ln K = −ΔG°/RT, and d ln K/dT =
// ΔH°/RT², van 't Hoff's equation, for a rate's derivative with respect to
// T.
pub struct Standard {
pub:
	t        f64 // K
	dh       f64 // J/mol
	ds       f64 // J/(mol·K)
	dg       f64 // J/mol
	ln_k     f64
	dln_k_dt f64 // 1/K
}

// standard is the reaction's standard properties at T, from each species'
// ideal-gas H and S along its Cp, as thermo's Fluid takes them:
//
//     H°(T) = ΔH°f + ∫ Cp dT,   S°(T) = S° + ∫ Cp/T dT,   from 298.15 K.
//
// A species without ΔH°f, S° or ideal-gas data is refused, by name.
pub fn (r Reaction) standard(species []Species, t f64) !Standard {
	if !(t > 0) || math.is_inf(t, 1) {
		return fail(r.name, r.text, 0, 'no standard properties at T = ${t} K')
	}
	mut dh := 0.0
	mut ds := 0.0
	for j, k in r.species {
		gas, dhf, s0 := formation(r, species[k], r.written[j])!
		h := dhf + gas.cp.int_cp(thermo.t_ref, t)
		s := s0 + gas.cp.int_cp_over_t(thermo.t_ref, t)
		if !math.is_finite(h) || !math.is_finite(s) {
			return fail(r.name, r.text, 0, "${r.written[j]}'s Cp doesn't integrate from ${thermo.t_ref} to ${t} K")
		}
		dh += r.nu[j] * h
		ds += r.nu[j] * s
	}
	dg := dh - t * ds
	return Standard{
		t:        t
		dh:       dh
		ds:       ds
		dg:       dg
		ln_k:     -dg / (gas_constant * t)
		dln_k_dt: dh / (gas_constant * t * t)
	}
}

fn formation(r Reaction, s Species, written string) !(thermo.IdealGas, f64, f64) {
	gas := s.gas or { return fail(r.name, r.text, 0, '${written} has no ideal-gas data') }
	dhf := gas.dhf or { return fail(r.name, r.text, 0, '${written} has no ΔH°f') }
	s0 := gas.s0 or { return fail(r.name, r.text, 0, '${written} has no S°') }
	return gas, dhf, s0
}

// formation_gap is Σν ΔG°f less ΔG° at 298.15 K from ΔH°f and S°, J/mol:
// zero where a data set's Gibbs energies of formation agree with its
// enthalpies and entropies, the elements' entropies cancelling in a
// balanced equation.
pub fn (r Reaction) formation_gap(species []Species) !f64 {
	mut dgf := 0.0
	for j, k in r.species {
		g := species[k].dgf or { return fail(r.name, r.text, 0, '${r.written[j]} has no ΔG°f') }
		dgf += r.nu[j] * g
	}
	st := r.standard(species, thermo.t_ref)!
	return dgf - st.dg
}
