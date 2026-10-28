module numeric

// The solvers take what they work on as an interface: a small struct of the
// caller's that carries its own data, with one method. A Go solver would
// take a function value, and callers would pass closures over their state.
// But V 0.5.2 never frees a closure, and a flash that made one on every
// call would leak a few hundred bytes each time, gigabytes over a
// sensitivity study (see docs/v-for-go.md).

// Function is f(x), for integrate and bisect.
pub interface Function {
	at(x f64) f64
}

// FallibleFunction is f(x) where evaluating f can fail, as a flash can: for
// brent.
pub interface FallibleFunction {
	at(x f64) !f64
}

// System is a square nonlinear system F(x) = 0, for newton. eval writes
// F(x) into f and the Jacobian ∂Fᵢ/∂xⱼ into jac, row-major (jac[i*n + j]),
// the layout lu_factor takes. It returns an error when x is outside the
// system's domain, a negative mole number say; newton then takes a
// shorter step.
//
// An interface method can mark parameters `mut`, as here: the solver owns
// f and jac and lets the system write into them, with no allocation per
// call.
pub interface System {
	eval(x []f64, mut f []f64, mut jac []f64) !
}

// Objective is a smooth function to minimise, for minimize. eval returns
// f(x) and writes the gradient into g and the Hessian into h, n×n
// row-major. It returns an error when x is outside the function's domain;
// minimize then takes a shorter step.
pub interface Objective {
	eval(x []f64, mut g []f64, mut h []f64) !f64
}
