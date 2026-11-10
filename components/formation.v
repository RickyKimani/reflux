module components

import chem
import thermo

// element_entropies are the elements' standard molar entropies in their
// reference states at 298.15 K and 1 bar, J/(mol·K) per mole of atoms, by
// atomic number: CODATA's key values (Cox, Wagman and Medvedev, CODATA Key
// Values for Thermodynamics, 1989), for the elements in the data. The
// diatomic gases' are half the molecule's, and bromine's (a liquid) and
// iodine's (a solid) half theirs; carbon's is graphite's, sulfur's
// rhombic sulfur's, mercury's the liquid's.
const element_entropies = {
	1:  130.680 / 2
	2:  126.153
	6:  5.74
	7:  191.609 / 2
	8:  205.152 / 2
	9:  202.791 / 2
	10: 146.328
	14: 18.81
	16: 32.054
	17: 223.081 / 2
	18: 154.846
	35: 152.210 / 2
	36: 164.085
	53: 116.14 / 2
	54: 169.685
	80: 75.90
}

// formation_gap is how far a compound's ΔG°f is from its ΔH°f and S°, J/mol:
//
//     ΔG°f − [ΔH°f − T°(S° − Σ nₑ S°ₑ)],
//
// with the elements' entropies S°ₑ in their reference states. A data set
// whose three agree has no gap, and a reaction's K is the same from ΔG°f
// as from ΔH°f and S°. Refused for a compound without the three, a
// formula, or an element whose entropy isn't held.
pub fn (r Record) formation_gap() !f64 {
	what := '${r.name} (${r.ids.cas})'
	dhf := r.constants.dhf or { return error('${what} has no ΔH°f') }
	dgf := r.constants.dgf or { return error('${what} has no ΔG°f') }
	s0 := r.constants.s0 or { return error('${what} has no S°') }
	f := r.formula()!
	mut elements := 0.0
	for c in f.counts {
		s := element_entropies[c.z] or {
			return error('${what}: no reference entropy for ${chem.symbol(c.z)}')
		}
		elements += f64(c.n) * s
	}
	return dgf - (dhf - thermo.t_ref * (s0 - elements))
}
