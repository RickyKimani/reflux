module fit

import flowsheet
import math
import numeric

// Fitted is a fit's outcome (M8's D14, D15): each parameter's value, its
// standard error and 95% interval by linearisation, their correlations,
// and each residual; or why there's none. Every number is SI, each
// parameter in the file's own form.
pub struct Fitted {
pub:
	plan            flowsheet.FitPlan
	ok              bool
	error           string
	values          []f64
	standard_errors []f64
	lower           []f64
	upper           []f64
	correlation     [][]f64
	ssr             f64 // Σ of the weighted residuals squared
	dof             int // the residuals less the parameters
	t               f64 // Student's t at 0.975 on dof, the intervals' half-width over a standard error
	centred         []Centred
	residuals       []Residual
	iterations      int
	evaluations     int
	gradient        f64 // the gradient test's value at the answer
}

// Centred is a constant fitted centred on a temperature (D14), its A
// fitted as its value there, k_ref, and reported back as A.
pub struct Centred {
pub:
	constant string
	t_ref    f64
}

// Residual is one measurement against the model: its row in the data, from
// 0, its column, the measured and the model's value in SI, and the
// residual over its weight.
pub struct Residual {
pub:
	row      int
	column   string
	measured f64
	model    f64
	weighted f64
}

// statistics turns a least-squares solution into the file's parameters
// and their errors. to_file maps θ, as fitted, to the file's form, and
// jacobian is its derivative, n×n row-major: the covariance s²·(JᵀJ)⁻¹ in θ
// goes to M·C·Mᵀ in the file's form, M = ∂(file)/∂θ, the delta method, which
// is exact for the linearised covariance.
fn statistics(s numeric.LeastSquaresSolution, m int, values []f64, jacobian []f64) !([]f64, []f64, []f64, [][]f64, int, f64) {
	n := values.len
	dof := m - n
	if dof < 1 {
		return error('${m} residuals for ${n} parameters leave no degrees of freedom for their errors')
	}
	s2 := s.ssr / f64(dof)
	mut cov := []f64{len: n * n}
	for p in 0 .. n {
		for q in 0 .. n {
			mut v := 0.0
			for a in 0 .. n {
				for b in 0 .. n {
					v += jacobian[p * n + a] * s.inverse[a * n + b] * jacobian[q * n + b]
				}
			}
			cov[p * n + q] = s2 * v
		}
	}
	t := numeric.student_t_quantile(0.975, f64(dof))!
	mut se := []f64{len: n}
	mut lower := []f64{len: n}
	mut upper := []f64{len: n}
	for p in 0 .. n {
		se[p] = math.sqrt(cov[p * n + p])
		lower[p] = values[p] - t * se[p]
		upper[p] = values[p] + t * se[p]
	}
	mut corr := [][]f64{len: n, init: []f64{len: n}}
	for p in 0 .. n {
		for q in 0 .. n {
			corr[p][q] = if p == q { 1.0 } else { cov[p * n + q] / (se[p] * se[q]) }
		}
	}
	return se, lower, upper, corr, dof, t
}

// failed is a fit with no answer, and why.
fn failed(plan flowsheet.FitPlan, reason string) Fitted {
	return Fitted{
		plan:  plan
		error: reason
	}
}

// run runs one of a file's fits: of rate constants to rates (D14), or of
// kᵢⱼ or an activity model's pairs to bubble points (D15, M6's D16).
pub fn run(file flowsheet.FitFile, plan flowsheet.FitPlan) Fitted {
	data := read_data(plan.path) or { return failed(plan, err.msg()) }
	return match plan.kind {
		.rates { fit_rates(file, plan, data) }
		.vle { fit_vle(file, plan, data) }
	}
}

// weight_for is a column's weight, if the plan gives one, checked against
// the column's dimension; a Weight with no column if it gives none.
fn weight_for(plan flowsheet.FitPlan, data &Data, column string) !flowsheet.Weight {
	for w in plan.weights {
		if w.column == column {
			if !w.relative {
				k := data.index(column)
				d := data.columns[k].dim()!
				if !w.dim.eq(d) {
					return error('the weight of ${column} is ${w.dim.name()}, but the column is ${d.name()}')
				}
			}
			return w
		}
	}
	return flowsheet.Weight{}
}

// unsolved is least squares' failure in the fit's words: a parameter the
// data don't determine, by its name.
fn unsolved(plan flowsheet.FitPlan, err IError) string {
	if err is numeric.LeastSquaresError {
		e := err as numeric.LeastSquaresError
		if e.undetermined >= 0 && e.undetermined < plan.parameters.len {
			return "${plan.parameters[e.undetermined].name} is not determined by the data: its effect on the residuals is a combination of the other parameters'"
		}
	}
	return err.msg()
}
