module fit

import expr
import flowsheet
import kinetics
import math
import numeric
import quantity

// R as expr has it, N_A·k_B, the constant the rates' exponents use.
const gas_r = 8.31446261815324

// RateColumn is a measured rate: its column in the data, its rate among
// the kinetics', and its weight.
struct RateColumn {
	column   int
	rate     int
	name     string
	relative bool
	sigma    f64 // a fraction of each measurement, or SI
}

// RateProblem is rate data as least squares' residuals (D14): each
// measurement less the model, over its weight, at each row's state, on the
// ideal gas. The fitted constants are the tape's inputs after T and P, so
// that ∂r/∂θ is the tape's own reverse sweep, exact.
struct RateProblem {
	k      kinetics.Kinetics
	states []kinetics.State
	cols   []RateColumn
	data   Data
	n      int
mut:
	pt  kinetics.Point
	buf expr.Buffer
	out []f64
	g   []f64
}

fn (mut p RateProblem) eval(theta []f64, mut r []f64, mut j []f64) ! {
	inputs := p.k.tape.inputs
	mut i := 0
	for row, st in p.states {
		p.k.fill(st, mut p.pt)!
		for q in 0 .. p.n {
			p.pt.x[2 + q] = theta[q]
		}
		p.k.rates(&p.pt, mut p.buf, mut p.out)
		p.k.tape.gradient(mut p.buf, mut p.g)
		for c in p.cols {
			measured := p.data.rows[row][c.column]
			if math.is_nan(measured) {
				continue
			}
			sigma := if c.relative { c.sigma * math.abs(measured) } else { c.sigma }
			r[i] = (p.out[c.rate] - measured) / sigma
			for q in 0 .. p.n {
				j[i * p.n + q] = p.g[c.rate * inputs + 2 + q] / sigma
			}
			i++
		}
	}
}

// fit_rates fits a file's named constants to rate data (D14): each
// measured rate's residual weighted by its error, the constants the tape's
// inputs, an Arrhenius or van 't Hoff constant whose A and energy are both
// fitted centred on the data's mean T, and every answer reported in the
// file's own form.
fn fit_rates(file flowsheet.FitFile, plan flowsheet.FitPlan, data Data) Fitted {
	kin := file.kinetics or { return failed(plan, 'the file has no reactions') }
	t_col, p_col := data.index('T'), data.index('P')
	if t_col < 0 || p_col < 0 {
		return failed(plan, '${data.file}: rate data need columns T and P, with their units')
	}
	mut y_cols := map[int]int{} // a species, to its column
	mut cols := []RateColumn{}
	for k, c in data.columns {
		if k == t_col || k == p_col {
			continue
		}
		if c.name.starts_with('y[') && c.name.ends_with(']') {
			sp := c.name[2..c.name.len - 1]
			i := kinetics.find(kin.species, sp) or {
				return failed(plan, '${data.file}: column ${c.name}: ${err.msg()}')
			}
			if c.unit != '' && c.unit != '%' {
				return failed(plan, '${data.file}: column ${c.name} is a mole fraction, with no unit or %')
			}
			y_cols[i] = k
			continue
		}
		if c.name.starts_with('r[') && c.name.ends_with(']') {
			name := c.name[2..c.name.len - 1]
			mut reaction := -1
			for i, rx in kin.reactions {
				if rx.name == name {
					reaction = i
				}
			}
			rate := kin.rated.index(reaction)
			if reaction < 0 || rate < 0 {
				mut rated := []string{}
				for i in kin.rated {
					rated << kin.reactions[i].name
				}
				hint := quantity.closest(name, rated)
				more := if hint.len > 0 { ' (did you mean ${hint.join(', ')}?)' } else { '' }
				return failed(plan, '${data.file}: column ${c.name}: no reaction ${name} with a rate${more}')
			}
			d := c.dim() or { return failed(plan, '${data.file}: column ${c.name}: ${err.msg()}') }
			if !d.eq(kin.rates[rate].dim) {
				return failed(plan, '${data.file}: column ${c.name} is ${d.name()}, but its rate is ${kin.rates[rate].dim.name()}')
			}
			wt := weight_for(plan, &data, c.name) or { return failed(plan, err.msg()) }
			cols << if wt.column != '' {
				RateColumn{
					column:   k
					rate:     rate
					name:     c.name
					relative: wt.relative
					sigma:    wt.value
				}
			} else {
				RateColumn{
					column: k
					rate:   rate
					name:   c.name
					sigma:  1.0
				}
			}
			continue
		}
		return failed(plan, '${data.file}: no column ${c.name} in rate data: T, P, y[species] and r[reaction]')
	}
	if cols.len == 0 {
		return failed(plan, '${data.file}: no measured rate, a column r[reaction]')
	}
	for w in plan.weights {
		if !cols.any(it.name == w.column) {
			return failed(plan, 'a weight for ${w.column}, which is not a measured rate in ${data.file}')
		}
	}
	// The states, and the measurements' count.
	mut states := []kinetics.State{}
	mut m := 0
	mut tsum := 0.0
	for r, row in data.rows {
		line := r + 2
		if math.is_nan(row[t_col]) || math.is_nan(row[p_col]) {
			return failed(plan, '${data.file}:${line}: T and P are needed in every row')
		}
		mut y := []f64{len: kin.species.len}
		mut total := 0.0
		for i, k in y_cols {
			v := row[k]
			if math.is_nan(v) || v < 0.0 {
				return failed(plan, '${data.file}:${line}: ${data.columns[k].name} must be a fraction, at least 0')
			}
			y[i] = v
			total += v
		}
		if total > 1.0 + 1e-9 {
			return failed(plan, '${data.file}:${line}: the mole fractions add up to ${total}, above 1')
		}
		for c in cols {
			v := row[c.column]
			if math.is_nan(v) {
				continue
			}
			if c.relative && v == 0.0 {
				return failed(plan, '${data.file}:${line}: ${c.name} is 0, which a relative weight makes an error of 0')
			}
			m++
		}
		states << kinetics.State{
			t: row[t_col]
			p: row[p_col]
			y: y
		}
		tsum += row[t_col]
	}
	tbar := tsum / f64(states.len)
	// The constants, written with the fitted numbers as variables, and
	// centred where both of a pair are fitted.
	n := plan.parameters.len
	mut variables := []expr.Variable{}
	mut replace := map[string]kinetics.TempFn{}
	mut centred := []Centred{}
	mut partner := []int{len: n, init: -1} // a centred A's energy, by index
	for par in plan.parameters {
		variables << expr.Variable{var_name(par), par.dim}
		if par.constant in replace {
			continue
		}
		mut c := ?kinetics.NamedConstant(none)
		for nc in kin.constants() {
			if nc.name == par.constant {
				c = nc
			}
		}
		constant := c or { return failed(plan, 'no constant ${par.constant}') }
		mine := plan.parameters.filter(it.constant == par.constant)
		match constant.value {
			kinetics.Arrhenius {
				a := constant.value as kinetics.Arrhenius
				both := has(mine, .a) && has(mine, .ea)
				centre := both && a.t_ref.trim_space() == ''
				replace[par.constant] = kinetics.Arrhenius{
					a:     if has(mine, .a) { 'fit_${par.constant}_A' } else { a.a }
					ea:    if has(mine, .ea) { 'fit_${par.constant}_Ea' } else { a.ea }
					t_ref: if centre { '${tbar} K' } else { a.t_ref }
				}
				if centre {
					centred << Centred{par.constant, tbar}
				}
			}
			kinetics.VantHoff {
				v := constant.value as kinetics.VantHoff
				both := has(mine, .a) && has(mine, .dh)
				centre := both && v.t_ref.trim_space() == ''
				replace[par.constant] = kinetics.VantHoff{
					a:     if has(mine, .a) { 'fit_${par.constant}_A' } else { v.a }
					dh:    if has(mine, .dh) { 'fit_${par.constant}_dH' } else { v.dh }
					t_ref: if centre { '${tbar} K' } else { v.t_ref }
				}
				if centre {
					centred << Centred{par.constant, tbar}
				}
			}
			kinetics.Fixed, kinetics.Written {
				// A number with its unit, as the file's reader checked.
				replace[par.constant] = kinetics.Fixed{'fit_${par.constant}'}
			}
		}
	}
	for q, par in plan.parameters {
		if par.part == .a && centred.any(it.constant == par.constant) {
			for e, other in plan.parameters {
				if other.constant == par.constant && other.part in [.ea, .dh] {
					partner[q] = e
				}
			}
		}
	}
	fk := kin.with_constants(replace, variables) or { return failed(plan, err.msg()) }
	mut theta0 := plan.parameters.map(it.start)
	for q in 0 .. n {
		if partner[q] >= 0 {
			theta0[q] = plan.parameters[q].start * math.exp(-plan.parameters[partner[q]].start / (gas_r * tbar))
		}
	}
	mut problem := RateProblem{
		k:      fk
		states: states
		cols:   cols
		data:   data
		n:      n
		pt:     fk.point()
		buf:    fk.buffer()
		out:    []f64{len: fk.rates.len}
		g:      []f64{len: fk.rates.len * fk.tape.inputs}
	}
	s := numeric.least_squares(mut problem, theta0, m) or { return failed(plan, unsolved(plan, err)) }
	// θ back to the file's form, and its derivative for the delta method.
	mut values := s.theta.clone()
	mut jac := []f64{len: n * n}
	for q in 0 .. n {
		jac[q * n + q] = 1.0
	}
	for q in 0 .. n {
		e := partner[q]
		if e >= 0 {
			grow := math.exp(s.theta[e] / (gas_r * tbar))
			values[q] = s.theta[q] * grow
			jac[q * n + q] = grow
			jac[q * n + e] = values[q] / (gas_r * tbar)
		}
	}
	se, lower, upper, corr, dof, t := statistics(s, m, values, jac) or {
		return failed(plan, err.msg())
	}
	// Each measurement against the model at the answer.
	mut residuals := []Residual{}
	mut i := 0
	for row, st in states {
		fk.fill(st, mut problem.pt) or { return failed(plan, err.msg()) }
		for q in 0 .. n {
			problem.pt.x[2 + q] = s.theta[q]
		}
		fk.rates(&problem.pt, mut problem.buf, mut problem.out)
		for c in cols {
			measured := data.rows[row][c.column]
			if math.is_nan(measured) {
				continue
			}
			residuals << Residual{
				row:      row
				column:   c.name
				measured: measured
				model:    problem.out[c.rate]
				weighted: s.r[i]
			}
			i++
		}
	}
	return Fitted{
		plan:            plan
		ok:              true
		values:          values
		standard_errors: se
		lower:           lower
		upper:           upper
		correlation:     corr
		ssr:             s.ssr
		dof:             dof
		t:               t
		centred:         centred
		residuals:       residuals
		iterations:      s.iterations
		evaluations:     s.evaluations
		gradient:        s.gradient
	}
}

// has is whether a constant's fitted parameters include this part.
fn has(mine []flowsheet.FitParameter, part flowsheet.Part) bool {
	return mine.any(it.part == part)
}

// var_name is a parameter's variable in the rates' context.
fn var_name(p flowsheet.FitParameter) string {
	return match p.part {
		.a { 'fit_${p.constant}_A' }
		.ea { 'fit_${p.constant}_Ea' }
		.dh { 'fit_${p.constant}_dH' }
		.value { 'fit_${p.constant}' }
		.kij { 'fit_kij_${p.i}_${p.j}' }
	}
}
