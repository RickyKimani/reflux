module thermo

// Gas constant, J/(mol·K). Exact since the 2019 SI redefinition.
pub const r_gas = 8.314462618

// Component holds the constants every model can rely on. Parameters that
// belong to one model (PRSV's κ₁, CPA's association energy, PC-SAFT's
// segment number) live in that model, not here.
//
// Everything inside the simulator is SI (K, Pa, m³, mol, J); unit
// conversion happens only at the edges (input files, printing).
//
// @[required] makes leaving a field out of a struct literal a compile
// error. Go would silently give you Tc = 0 and a NaN three calls later.
pub struct Component {
pub:
	name  string @[required]
	tc    f64    @[required] // critical temperature, K
	pc    f64    @[required] // critical pressure, Pa
	omega f64    @[required] // acentric factor, dimensionless
	mw    f64    @[required] // molar mass, kg/mol
}

// InputError is a typed error, the V counterpart of a Go error struct you
// would match with errors.As. Embedding `Error` gives it a default code()
// so only msg() has to be written.
pub struct InputError {
	Error
pub:
	what  string
	value f64
}

// msg makes InputError satisfy the built-in IError interface.
pub fn (e InputError) msg() string {
	return 'invalid ${e.what}: ${e.value}'
}

// check returns an error if the constants are not physical. A `!` with
// no type after it is the V spelling of a Go func that returns only `error`.
fn (c Component) check() ! {
	if c.tc <= 0.0 {
		return InputError{
			what:  '${c.name} critical temperature (K)'
			value: c.tc
		}
	}
	if c.pc <= 0.0 {
		return InputError{
			what:  '${c.name} critical pressure (Pa)'
			value: c.pc
		}
	}
	if c.mw <= 0.0 {
		return InputError{
			what:  '${c.name} molar mass (kg/mol)'
			value: c.mw
		}
	}
}
