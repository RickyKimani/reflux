// Module transport is M5b's (docs/m5b-plan.md): a phase's viscosity,
// thermal conductivity and surface tension, from each compound's
// correlations by mixing rules, with pressure corrections for a dense gas
// and a compressed liquid. It knows nothing of component files: the
// components module builds its compounds (D2), as it builds activity's
// models.
module transport

import math
import tdep

// Quantity is one of a compound's correlations transport uses.
pub enum Quantity {
	liquid_viscosity
	vapour_viscosity
	liquid_conductivity
	vapour_conductivity
	surface_tension
	liquid_density
	vapour_pressure
}

// quantities are every Quantity, in its order, the order of Pure.sources.
pub const quantities = [Quantity.liquid_viscosity, .vapour_viscosity, .liquid_conductivity,
	.vapour_conductivity, .surface_tension, .liquid_density, .vapour_pressure]

// key is the quantity's name in a component file.
pub fn (q Quantity) key() string {
	return match q {
		.liquid_viscosity { 'liquid_viscosity' }
		.vapour_viscosity { 'vapor_viscosity' }
		.liquid_conductivity { 'liquid_thermal_conductivity' }
		.vapour_conductivity { 'vapor_thermal_conductivity' }
		.surface_tension { 'surface_tension' }
		.liquid_density { 'liquid_density' }
		.vapour_pressure { 'vapor_pressure' }
	}
}

// saturated says whether q is the saturated liquid's, which ends at Tc
// (D3): its correlations are fitted to it, and 190 of ChemSep's surface
// tensions go negative by Tc.
pub fn (q Quantity) saturated() bool {
	return q in [.liquid_viscosity, .liquid_conductivity, .surface_tension, .liquid_density,
		.vapour_pressure]
}

// Correlated is one of a compound's correlations in SI: value =
// scale × c.eval(T), the arithmetic of components.Property, so that a value
// is the same number either way.
pub struct Correlated {
pub:
	key   string // the correlation used, e.g. 'liquid_viscosity'
	c     tdep.Correlation
	scale f64
}

// at is the value at T, in SI.
pub fn (c Correlated) at(t f64) f64 {
	return c.scale * c.c.eval(t)
}

// Source is a compound's correlation for a quantity, or why it has none.
pub struct Source {
pub:
	has bool
	c   Correlated
	why string // when it has none: the component file's lack, in its words
}

// Pure is a compound as transport needs it: its constants in SI per mol,
// and its correlations, each by Quantity. A constant the component file
// lacks is NaN, and a method that needs it says so.
pub struct Pure {
pub:
	name    string
	cas     string
	tc      f64 // K
	pc      f64 // Pa
	vc      f64 // m³/mol
	zc      f64
	omega   f64
	mw      f64      // kg/mol
	dipole  f64      // debye; NaN where the file has none
	quantum f64      // Lucas's Q: 1.38 for He, 0.76 for H₂, 0.52 for D₂, 0 for the rest
	sources []Source // by Quantity, in quantities' order
}

// Reading is a quantity's value at T and whether T was inside its
// correlation's fitted range: outside it the value is the correlation's
// extrapolation (M2's D7), which a phase's notes name.
pub struct Reading {
pub:
	value    f64
	in_range bool
}

// value is q at T, or why there's none: the compound has no such
// correlation, a saturated liquid's quantity at or above the compound's Tc
// (D3), or a value that isn't a positive finite number.
pub fn (p Pure) value(q Quantity, t f64) !Reading {
	idx := int(q)
	if idx >= p.sources.len {
		return error('${p.name} has no ${q.key()} data')
	}
	s := p.sources[idx]
	if !s.has {
		return error(s.why)
	}
	if q.saturated() && !(t < p.tc) {
		return error("${p.name}'s ${q.key()} is the saturated liquid's, and ${shown(t)} K isn't below its Tc, ${p.tc} K")
	}
	v := s.c.at(t)
	if !(v > 0) || math.is_inf(v, 1) {
		return error("${p.name}'s ${q.key()} at ${shown(t)} K is ${v}, which isn't one")
	}
	return Reading{
		value:    v
		in_range: s.c.c.in_range(t)
	}
}

// outside is a note that q was taken at T beyond its fitted range.
fn (p Pure) outside(q Quantity, t f64) string {
	c := p.sources[int(q)].c.c
	return "${p.name}'s ${q.key()} at ${shown(t)} K, outside its fitted range, ${c.tmin} K to ${c.tmax} K"
}

// shown is x to 6 significant digits, for a message (M5b's D16): a report's
// note then reads the same from tcc's build and gcc's, whose last bits of a
// state can differ, and reads as a number a person would write. Rounded by
// math.round on an exact power of ten, never printf (whose ties differ).
pub fn shown(x f64) f64 {
	if x == 0.0 || !math.is_finite(x) {
		return x
	}
	k := 5 - int(math.floor(math.log10(math.abs(x))))
	if k > 0 && k <= 22 {
		f := math.pow10(k)
		return math.round(x * f) / f
	}
	if k <= 0 && k >= -22 {
		f := math.pow10(-k)
		return math.round(x / f) * f
	}
	return x
}
