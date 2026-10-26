// Module components holds reflux's component data: what a component file
// in data/components/ contains, and (from M2 step 5) the database that
// looks compounds up and builds thermo models from them.
//
// The format is described in docs/data-format.md. The field names of the
// structs below are the file's keys, in the file's order.
module components

// Record is one compound as a component file holds it: identifiers,
// constants in SI per mol, model parameters, correlations with their
// coefficients verbatim from the source, and group counts.
pub struct Record {
pub mut:
	name         string
	source       string // where every value comes from, unless `sources` says otherwise
	ids          Ids
	constants    Constants
	parameters   Parameters
	correlations Correlations
	groups       Groups
	// sources holds the source of each value that doesn't come from
	// `source`, keyed like 'constants.tc'.
	sources map[string]string
}

// Ids are a compound's identifiers. Only the CAS number is required.
pub struct Ids {
pub mut:
	cas           string
	inchikey      ?string
	pubchem       ?int // PubChem CID
	smiles        ?string
	structure     ?string // ChemSep's structural formula, e.g. CH3(CH2)2CH3
	chemsep_index ?int
	family        ?int // ChemSep's family number
	aliases       []string
}

// Constants are scalar properties in SI per mol. The thermochemical ones
// (dhf, dgf, s0) are for the ideal gas at 298.15 K, s0 also at 1 bar.
pub struct Constants {
pub mut:
	tc                   ?f64 // K
	pc                   ?f64 // Pa
	vc                   ?f64 // m³/mol
	zc                   ?f64
	omega                ?f64
	mw                   ?f64 // kg/mol
	tb                   ?f64 // K, normal boiling point
	tm                   ?f64 // K, melting point
	t_triple             ?f64 // K
	p_triple             ?f64 // Pa
	v_liquid_tb          ?f64 // m³/mol, liquid at the normal boiling point
	dhf                  ?f64 // J/mol
	dgf                  ?f64 // J/mol
	s0                   ?f64 // J/(mol·K)
	dh_fusion            ?f64 // J/mol, at the melting point
	dh_combustion        ?f64 // J/mol, net (LHV)
	dipole               ?f64 // C·m
	radius_of_gyration   ?f64 // m
	solubility_parameter ?f64 // J^0.5/m^1.5
	charge               ?int
}

// Parameters belong to particular models, in SI where they have units.
pub struct Parameters {
pub mut:
	rackett          ?f64
	costald_volume   ?f64 // m³/mol, COSTALD's V*
	costald_omega    ?f64 // COSTALD's acentric factor (ChemSep calls it the SRK one)
	wilson_volume    ?f64 // m³/mol
	uniquac_r        ?f64
	uniquac_q        ?f64
	uniquac_q_prime  ?f64
	vdw_volume       ?f64 // m³/mol, van der Waals (Bondi)
	vdw_area         ?f64 // m²/mol
	lj_sigma         ?f64 // m, Lennard-Jones diameter
	lj_epsilon_k     ?f64 // K, Lennard-Jones energy over k
	fuller_volume    ?f64 // Fuller's diffusion volume
	parachor         ?f64 // kg^0.25·m³/(s^0.5·mol)
	specific_gravity ?f64
	// Mathias–Copeman α constants. ChemSep doesn't say which EOS they were
	// fitted for, so reflux doesn't use them by default. They're separate
	// because five ChemSep compounds have only some of them.
	mathias_copeman_c1               ?f64
	mathias_copeman_c2               ?f64
	mathias_copeman_c3               ?f64
	api_srk_s1                       ?f64
	api_srk_s2                       ?f64
	chao_seader_omega                ?f64
	chao_seader_solubility_parameter ?f64 // J^0.5/m^1.5
	chao_seader_liquid_volume        ?f64 // m³/mol
}

// CorrelationData is a temperature-dependent correlation as stored: its
// equation number and coefficients exactly as the source gives them, in the
// source's units. tdep evaluates it.
pub struct CorrelationData {
pub mut:
	eq     int
	coeffs []f64  // A, B, C, …; never converted
	tmin   f64    // K
	tmax   f64    // K
	units  string // the output's unit, as the source writes it, e.g. 'J/kmol/K'
	tc     ?f64   // K; for the reduced-temperature forms, when it isn't constants.tc
}

// Correlations are the temperature-dependent properties.
pub struct Correlations {
pub mut:
	vapor_pressure              ?CorrelationData
	antoine_vapor_pressure      ?CorrelationData
	ideal_gas_cp                ?CorrelationData
	rpp_cp                      ?CorrelationData // ideal-gas Cp from Reid, Prausnitz and Poling
	liquid_cp                   ?CorrelationData
	solid_cp                    ?CorrelationData
	liquid_density              ?CorrelationData
	solid_density               ?CorrelationData
	heat_of_vaporization        ?CorrelationData
	second_virial               ?CorrelationData
	liquid_viscosity            ?CorrelationData
	liquid_viscosity_rps        ?CorrelationData
	vapor_viscosity             ?CorrelationData
	liquid_thermal_conductivity ?CorrelationData
	vapor_thermal_conductivity  ?CorrelationData
	surface_tension             ?CorrelationData
	relative_permittivity       ?CorrelationData
}

// Groups are group counts for group-contribution methods, keyed by the
// method's group number. An empty map means no groups for that method.
pub struct Groups {
pub mut:
	unifac_vle      map[int]int
	unifac_lle      map[int]int
	modified_unifac map[int]int // Dortmund
	asog            map[int]int
	umr             map[int]int
	gc_method       map[int]int // PPR78
}

// Setting a value by its key, for the importer and (step 4) the reader.
// Each struct's fields are visited at compile time, so a key is exactly a
// field name, and an unknown key is an error rather than a silent no-op.

// set sets constant `key`.
pub fn (mut c Constants) set(key string, v f64) ! {
	set_field(mut c, key, v)!
}

// set sets parameter `key`.
pub fn (mut p Parameters) set(key string, v f64) ! {
	set_field(mut p, key, v)!
}

// set sets correlation `key`.
pub fn (mut c Correlations) set(key string, d CorrelationData) ! {
	$for f in Correlations.fields {
		$if f.typ is ?CorrelationData {
			if f.name == key {
				c.$(f.name) = d
				return
			}
		}
	}
	return error('no correlation `${key}`')
}

// correlation_keys returns the names of every correlation in the format, in
// file order.
pub fn correlation_keys() []string {
	mut keys := []string{}
	$for f in Correlations.fields {
		keys << f.name
	}
	return keys
}

// get returns correlation `key`. It fails if the record has none, or if
// the format has no correlation of that name.
pub fn (c Correlations) get(key string) !CorrelationData {
	$for f in Correlations.fields {
		$if f.typ is ?CorrelationData {
			if f.name == key {
				// Unwrapped in a typed function: inside `$for`, `if d := …` on
				// an Option struct field leaves d an Option (a V 0.5.2 bug).
				return present(c.$(f.name), key)
			}
		}
	}
	return error('no correlation called `${key}` in the component format')
}

fn present(opt ?CorrelationData, key string) !CorrelationData {
	return opt or { return error('no `${key}` data') }
}

// set sets the group counts for method `key`.
pub fn (mut g Groups) set(key string, counts map[int]int) ! {
	$for f in Groups.fields {
		$if f.typ is map[int]int {
			if f.name == key {
				g.$(f.name) = counts.clone()
				return
			}
		}
	}
	return error('no group method `${key}`')
}

// set_field sets the ?f64 field of s named key.
fn set_field[T](mut s T, key string, v f64) ! {
	$for f in T.fields {
		$if f.typ is ?f64 {
			if f.name == key {
				s.$(f.name) = v
				return
			}
		}
	}
	return error('no ${typeof(s).name} key `${key}`')
}
