module fit

import flash
import flowsheet
import math
import numeric
import thermo

// The step of a kᵢⱼ's central difference: kᵢⱼ is of order 0.01, its bubble
// points' rounding about 1e-12, so the difference is good to about 1e-7.
const kij_step = 1e-6

// VleRow is a measured bubble point: T and the liquid's composition, and
// the measured P and vapour fractions, NaN where one isn't measured.
struct VleRow {
	t f64
	x []f64
	p f64
	y []f64
}

// VleProblem is VLE data as least squares' residuals (D15): each bubble
// point's P and each measured vapour fraction less the model's, over their
// weights, the model the file's cubic with the fitted kᵢⱼ. Its Jacobian is
// by central differences, each kᵢⱼ ± 1e-6.
struct VleProblem {
	base       thermo.Cubic
	parameters []flowsheet.FitParameter
	rows       []VleRow
	p_weight   flowsheet.Weight
	y_weights  []?flowsheet.Weight // one a component, none where y isn't measured
	y_names    []string            // each component's y column, '' for none
mut:
	flashes int
}

// residuals are the weighted residuals at θ, each row's P first and then
// each measured y.
fn (mut p VleProblem) residuals(theta []f64) ![]f64 {
	nc := p.base.comps.len
	mut kij := p.base.kij.clone()
	for q, par in p.parameters {
		kij[par.i * nc + par.j] = theta[q]
		kij[par.j * nc + par.i] = theta[q]
	}
	model := thermo.Cubic{
		...p.base
		kij: kij
	}
	mut out := []f64{}
	for row in p.rows {
		p.flashes++
		b := flash.bubble_p(model, row.t, row.x) or {
			return error('no bubble point at T = ${row.t} K and x = ${row.x}: ${err.msg()}')
		}
		if b.phases.len != 2 {
			return error('the bubble point at T = ${row.t} K and x = ${row.x} has ${b.phases.len} phases')
		}
		y := b.phases[0].state.n
		sigma_p := if p.p_weight.relative { p.p_weight.value * row.p } else { p.p_weight.value }
		out << (b.p - row.p) / sigma_p
		for i in 0 .. nc {
			if w := p.y_weights[i] {
				if !math.is_nan(row.y[i]) {
					out << (y[i] - row.y[i]) / (if w.relative {
						w.value * row.y[i]
					} else {
						w.value
					})
				}
			}
		}
	}
	return out
}

fn (mut p VleProblem) eval(theta []f64, mut r []f64, mut j []f64) ! {
	r0 := p.residuals(theta)!
	for i, v in r0 {
		r[i] = v
	}
	n := theta.len
	for q in 0 .. n {
		mut up := theta.clone()
		mut down := theta.clone()
		up[q] += kij_step
		down[q] -= kij_step
		ru := p.residuals(up)!
		rd := p.residuals(down)!
		for i in 0 .. r0.len {
			j[i * n + q] = (ru[i] - rd[i]) / (up[q] - down[q])
		}
	}
}

// fit_vle fits a file's kᵢⱼ to bubble points (D15): T, P, the liquid's
// composition and the vapour's, P and each y weighted by its error, which
// the file must give, since they're of different kinds.
fn fit_vle(file flowsheet.FitFile, plan flowsheet.FitPlan, data Data) Fitted {
	fl := file.fluid or { return failed(plan, 'the file has no [model]') }
	if fl.model !is thermo.Cubic {
		return failed(plan, "kij is a cubic's, and the model is not one")
	}
	base := fl.model as thermo.Cubic
	nc := file.names.len
	t_col, p_col := data.index('T'), data.index('P')
	if t_col < 0 || p_col < 0 {
		return failed(plan, '${data.file}: VLE data need columns T and P, with their units')
	}
	mut x_cols := []int{len: nc, init: -1}
	mut y_cols := []int{len: nc, init: -1}
	for k, c in data.columns {
		if k == t_col || k == p_col {
			continue
		}
		family := c.name.all_before('[')
		if (family == 'x' || family == 'y') && c.name.ends_with(']') {
			name := c.name.all_after('[').all_before_last(']')
			i := file.names.index(name)
			if i < 0 {
				return failed(plan, '${data.file}: column ${c.name}: ${name} is not one of the components')
			}
			if c.unit != '' && c.unit != '%' {
				return failed(plan, '${data.file}: column ${c.name} is a mole fraction, with no unit or %')
			}
			if family == 'x' {
				x_cols[i] = k
			} else {
				y_cols[i] = k
			}
			continue
		}
		return failed(plan, '${data.file}: no column ${c.name} in VLE data: T, P, x[component] and y[component]')
	}
	missing := x_cols.filter(it < 0).len
	if missing > 1 {
		return failed(plan, '${data.file}: the liquid needs every component but one, x[component]')
	}
	p_weight := weight_for(plan, &data, 'P') or { return failed(plan, err.msg()) }
	if p_weight.column == '' {
		return failed(plan, 'a VLE fit needs the error of P and of each measured y, as weights = { P = "0.5 %", "y[${file.names[0]}]" = "0.003" }')
	}

	mut y_weights := []?flowsheet.Weight{len: nc}
	for i in 0 .. nc {
		if y_cols[i] < 0 {
			continue
		}
		column := data.columns[y_cols[i]].name
		w := weight_for(plan, &data, column) or { return failed(plan, err.msg()) }
		if w.column == '' {
			return failed(plan, 'a VLE fit needs the error of P and of each measured y: give ${column} its weight')
		}
		y_weights[i] = w
	}
	for w in plan.weights {
		if w.column != 'P' && !y_cols.any(it >= 0 && data.columns[it].name == w.column) {
			return failed(plan, 'a weight for ${w.column}, which is not a measured column in ${data.file}')
		}
	}
	mut rows := []VleRow{}
	mut m := 0
	for r, row in data.rows {
		line := r + 2
		if math.is_nan(row[t_col]) || math.is_nan(row[p_col]) {
			return failed(plan, '${data.file}:${line}: T and P are needed in every row')
		}
		mut x := []f64{len: nc}
		mut total := 0.0
		mut rest := -1
		for i in 0 .. nc {
			if x_cols[i] < 0 {
				rest = i
				continue
			}
			v := row[x_cols[i]]
			if math.is_nan(v) || v < 0.0 || v > 1.0 {
				return failed(plan, '${data.file}:${line}: ${data.columns[x_cols[i]].name} must be a fraction')
			}
			x[i] = v
			total += v
		}
		if rest >= 0 {
			x[rest] = 1.0 - total
			if x[rest] < -1e-12 {
				return failed(plan, '${data.file}:${line}: the liquid fractions add up to ${total}, above 1')
			}
			x[rest] = math.max(x[rest], 0.0)
		}
		mut y := []f64{len: nc, init: math.nan()}
		for i in 0 .. nc {
			if y_cols[i] >= 0 {
				y[i] = row[y_cols[i]]
				if !math.is_nan(y[i]) {
					m++
				}
			}
		}
		m++
		rows << VleRow{
			t: row[t_col]
			x: x
			p: row[p_col]
			y: y
		}
	}
	mut problem := VleProblem{
		base:       base
		parameters: plan.parameters
		rows:       rows
		p_weight:   p_weight
		y_weights:  y_weights
		y_names:    y_cols.map(if it >= 0 { data.columns[it].name } else { '' })
	}
	n := plan.parameters.len
	s := numeric.least_squares(mut problem, plan.parameters.map(it.start), m) or {
		return failed(plan, unsolved(plan, err))
	}
	mut jac := []f64{len: n * n}
	for q in 0 .. n {
		jac[q * n + q] = 1.0
	}
	se, lower, upper, corr, dof, t := statistics(s, m, s.theta, jac) or {
		return failed(plan, err.msg())
	}
	return Fitted{
		plan:            plan
		ok:              true
		values:          s.theta
		standard_errors: se
		lower:           lower
		upper:           upper
		correlation:     corr
		ssr:             s.ssr
		dof:             dof
		t:               t
		residuals:       vle_residuals(mut problem, s.theta, s.r) or {
			return failed(plan, err.msg())
		}
		iterations:      s.iterations
		evaluations:     s.evaluations
		gradient:        s.gradient
	}
}

// vle_residuals is each measurement against the model at the answer.
fn vle_residuals(mut p VleProblem, theta []f64, weighted []f64) ![]Residual {
	nc := p.base.comps.len
	mut kij := p.base.kij.clone()
	for q, par in p.parameters {
		kij[par.i * nc + par.j] = theta[q]
		kij[par.j * nc + par.i] = theta[q]
	}
	model := thermo.Cubic{
		...p.base
		kij: kij
	}
	mut out := []Residual{}
	mut k := 0
	for r, row in p.rows {
		b := flash.bubble_p(model, row.t, row.x)!
		out << Residual{
			row:      r
			column:   'P'
			measured: row.p
			model:    b.p
			weighted: weighted[k]
		}
		k++
		for i in 0 .. nc {
			if _ := p.y_weights[i] {
				if !math.is_nan(row.y[i]) {
					out << Residual{
						row:      r
						column:   p.y_names[i]
						measured: row.y[i]
						model:    b.phases[0].state.n[i]
						weighted: weighted[k]
					}
					k++
				}
			}
		}
	}
	return out
}
