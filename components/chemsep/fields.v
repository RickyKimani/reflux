module chemsep

// The table that drives the import: what each ChemSep element becomes.

// What a ChemSep element becomes.
enum Kind {
	name         // the compound's name
	text_id      // a string identifier
	int_id       // an integer identifier
	constant     // a components.Constants value
	int_constant // the one integer constant, charge
	parameter    // a components.Parameters value
	correlation  // a components.Correlations entry
	groups       // a components.Groups entry
}

struct Field {
	kind     Kind
	key      string // the component-file key
	units    string // the units attribute the element must carry; '' for none
	per_kmol bool   // per kmol in ChemSep, so divided by 1000, exactly
}

// Every element a ChemSep 8.50 compound record can contain. The units are
// ChemSep's; '_' is its mark for dimensionless.
const fields = {
	'CompoundID':                       Field{
		kind: .name
	}
	'LibraryIndex':                     Field{
		kind: .int_id
		key:  'chemsep_index'
	}
	'StructureFormula':                 Field{
		kind: .text_id
		key:  'structure'
	}
	'Family':                           Field{
		kind: .int_id
		key:  'family'
	}
	'CAS':                              Field{
		kind: .text_id
		key:  'cas'
	}
	'SMILES':                           Field{
		kind: .text_id
		key:  'smiles'
	}
	'InChIKey':                         Field{
		kind: .text_id
		key:  'inchikey'
	}
	'PubChemID':                        Field{
		kind: .int_id
		key:  'pubchem'
	}
	'CriticalTemperature':              Field{
		kind:  .constant
		key:   'tc'
		units: 'K'
	}
	'CriticalPressure':                 Field{
		kind:  .constant
		key:   'pc'
		units: 'Pa'
	}
	'CriticalVolume':                   Field{
		kind:     .constant
		key:      'vc'
		units:    'm3/kmol'
		per_kmol: true
	}
	'CriticalCompressibility':          Field{
		kind:  .constant
		key:   'zc'
		units: '_'
	}
	'AcentricityFactor':                Field{
		kind:  .constant
		key:   'omega'
		units: '_'
	}
	'MolecularWeight':                  Field{
		kind:     .constant
		key:      'mw'
		units:    'kg/kmol'
		per_kmol: true
	}
	'NormalBoilingPointTemperature':    Field{
		kind:  .constant
		key:   'tb'
		units: 'K'
	}
	'NormalMeltingPointTemperature':    Field{
		kind:  .constant
		key:   'tm'
		units: 'K'
	}
	'TriplePointTemperature':           Field{
		kind:  .constant
		key:   't_triple'
		units: 'K'
	}
	'TriplePointPressure':              Field{
		kind:  .constant
		key:   'p_triple'
		units: 'Pa'
	}
	'LiquidVolumeAtNormalBoilingPoint': Field{
		kind:     .constant
		key:      'v_liquid_tb'
		units:    'm3/kmol'
		per_kmol: true
	}
	'HeatOfFormation':                  Field{
		kind:     .constant
		key:      'dhf'
		units:    'J/kmol'
		per_kmol: true
	}
	'GibbsEnergyOfFormation':           Field{
		kind:     .constant
		key:      'dgf'
		units:    'J/kmol'
		per_kmol: true
	}
	'AbsEntropy':                       Field{
		kind:     .constant
		key:      's0'
		units:    'J/kmol/K'
		per_kmol: true
	}
	'HeatOfFusionAtMeltingPoint':       Field{
		kind:     .constant
		key:      'dh_fusion'
		units:    'J/kmol'
		per_kmol: true
	}
	'HeatOfCombustion':                 Field{
		kind:     .constant
		key:      'dh_combustion'
		units:    'J/kmol'
		per_kmol: true
	}
	'DipoleMoment':                     Field{
		kind:  .constant
		key:   'dipole'
		units: 'Coulomb.m'
	}
	'RadiusOfGyration':                 Field{
		kind:  .constant
		key:   'radius_of_gyration'
		units: 'm'
	}
	'SolubilityParameter':              Field{
		kind:  .constant
		key:   'solubility_parameter'
		units: 'J0.5/m1.5'
	}
	'Charge':                           Field{
		kind:  .int_constant
		key:   'charge'
		units: '_'
	}
	'RacketParameter':                  Field{
		kind:  .parameter
		key:   'rackett'
		units: '_'
	}
	'COSTALDVolume':                    Field{
		kind:     .parameter
		key:      'costald_volume'
		units:    'm3/kmol'
		per_kmol: true
	}
	'CostaldAcentricFactor':            Field{
		kind:  .parameter
		key:   'costald_omega'
		units: '_'
	}
	'WilsonVolume':                     Field{
		kind:     .parameter
		key:      'wilson_volume'
		units:    'm3/kmol'
		per_kmol: true
	}
	'UniquacR':                         Field{
		kind:  .parameter
		key:   'uniquac_r'
		units: '_'
	}
	'UniquacQ':                         Field{
		kind:  .parameter
		key:   'uniquac_q'
		units: '_'
	}
	'UniquacQP':                        Field{
		kind:  .parameter
		key:   'uniquac_q_prime'
		units: '_'
	}
	'VanDerWaalsVolume':                Field{
		kind:     .parameter
		key:      'vdw_volume'
		units:    'm3/kmol'
		per_kmol: true
	}
	'VanDerWaalsArea':                  Field{
		kind:     .parameter
		key:      'vdw_area'
		units:    'm2/kmol'
		per_kmol: true
	}
	'DiameterLJ':                       Field{
		kind:  .parameter
		key:   'lj_sigma'
		units: 'm'
	}
	'EnergyLJ':                         Field{
		kind:  .parameter
		key:   'lj_epsilon_k'
		units: 'K'
	}
	'FullerVolume':                     Field{
		kind:  .parameter
		key:   'fuller_volume'
		units: '_'
	}
	'Parachor':                         Field{
		kind:     .parameter
		key:      'parachor'
		units:    'kg0.25.m3/s0.5/kmol'
		per_kmol: true
	}
	'SpecificGravity':                  Field{
		kind:  .parameter
		key:   'specific_gravity'
		units: '_'
	}
	'MatthiasCopemanC1':                Field{
		kind:  .parameter
		key:   'mathias_copeman_c1'
		units: '_'
	}
	'MatthiasCopemanC2':                Field{
		kind:  .parameter
		key:   'mathias_copeman_c2'
		units: '_'
	}
	'MatthiasCopemanC3':                Field{
		kind:  .parameter
		key:   'mathias_copeman_c3'
		units: '_'
	}
	'ApiSrkS1':                         Field{
		kind:  .parameter
		key:   'api_srk_s1'
		units: '_'
	}
	'ApiSrkS2':                         Field{
		kind:  .parameter
		key:   'api_srk_s2'
		units: '_'
	}
	'ChaoSeaderAcentricFactor':         Field{
		kind:  .parameter
		key:   'chao_seader_omega'
		units: '_'
	}
	'ChaoSeaderSolubilityParameter':    Field{
		kind:  .parameter
		key:   'chao_seader_solubility_parameter'
		units: 'J0.5/m1.5'
	}
	'ChaoSeaderLiquidVolume':           Field{
		kind:     .parameter
		key:      'chao_seader_liquid_volume'
		units:    'm3/kmol'
		per_kmol: true
	}
	'VaporPressure':                    Field{
		kind:  .correlation
		key:   'vapor_pressure'
		units: 'Pa'
	}
	'AntoineVaporPressure':             Field{
		kind:  .correlation
		key:   'antoine_vapor_pressure'
		units: 'Pa'
	}
	'IdealGasHeatCapacityCp':           Field{
		kind:  .correlation
		key:   'ideal_gas_cp'
		units: 'J/kmol/K'
	}
	'RPPHeatCapacityCp':                Field{
		kind:  .correlation
		key:   'rpp_cp'
		units: 'J/kmol/K'
	}
	'LiquidHeatCapacityCp':             Field{
		kind:  .correlation
		key:   'liquid_cp'
		units: 'J/kmol/K'
	}
	'SolidHeatCapacityCp':              Field{
		kind:  .correlation
		key:   'solid_cp'
		units: 'J/kmol/K'
	}
	'LiquidDensity':                    Field{
		kind:  .correlation
		key:   'liquid_density'
		units: 'kmol/m3'
	}
	'SolidDensity':                     Field{
		kind:  .correlation
		key:   'solid_density'
		units: 'kmol/m3'
	}
	'HeatOfVaporization':               Field{
		kind:  .correlation
		key:   'heat_of_vaporization'
		units: 'J/kmol'
	}
	'SecondVirialCoefficient':          Field{
		kind:  .correlation
		key:   'second_virial'
		units: 'm3/kmol'
	}
	'LiquidViscosity':                  Field{
		kind:  .correlation
		key:   'liquid_viscosity'
		units: 'Pa.s'
	}
	'LiquidViscosityRPS':               Field{
		kind:  .correlation
		key:   'liquid_viscosity_rps'
		units: 'Pa.s'
	}
	'VaporViscosity':                   Field{
		kind:  .correlation
		key:   'vapor_viscosity'
		units: 'Pa.s'
	}
	'LiquidThermalConductivity':        Field{
		kind:  .correlation
		key:   'liquid_thermal_conductivity'
		units: 'W/m/K'
	}
	'VaporThermalConductivity':         Field{
		kind:  .correlation
		key:   'vapor_thermal_conductivity'
		units: 'W/m/K'
	}
	'SurfaceTension':                   Field{
		kind:  .correlation
		key:   'surface_tension'
		units: 'N/m'
	}
	'RelativeStaticPermittivity':       Field{
		kind:  .correlation
		key:   'relative_permittivity'
		units: '_'
	}
	'UnifacVLE':                        Field{
		kind: .groups
		key:  'unifac_vle'
	}
	'UnifacLLE':                        Field{
		kind: .groups
		key:  'unifac_lle'
	}
	'ModifiedUnifac':                   Field{
		kind: .groups
		key:  'modified_unifac'
	}
	'Asog':                             Field{
		kind: .groups
		key:  'asog'
	}
	'Umr':                              Field{
		kind: .groups
		key:  'umr'
	}
	'GCmethod':                         Field{
		kind: .groups
		key:  'gc_method'
	}
}
