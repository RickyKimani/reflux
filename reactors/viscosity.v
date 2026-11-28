// Module reactors is M9's reactors (D10): conversion, equilibrium over
// chosen reactions, Gibbs, CSTR and PFR, each a flowsheet unit with a maker
// for the file (D1).
module reactors

// Viscosity is a gas mixture's viscosity, Pa·s, at T (K), P (Pa) and mole
// fractions y, for Ergun's pressure drop (D10, as amended in step 0). The
// PFR takes it as an interface, so a test can give it DWSIM's. Its own is
// Wilke's rule at low pressure, transport.Wilke (moved there by M5b's D2),
// which components.wilke builds from the component files, and which the
// PFR keeps at an unmoved T (M8's D11).
pub interface Viscosity {
	viscosity(t f64, p f64, y []f64) !f64
}
