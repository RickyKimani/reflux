module fit

import activity
import components
import flash
import flowsheet
import math
import numeric
import thermo

// The step of a kᵢⱼ's central difference: kᵢⱼ is of order 0.01, its bubble
// points' rounding about 1e-12, so the difference is good to about 1e-7.
const kij_step = 1e-6

// An activity model's pair is stepped further (M6's D16): a pair energy by
// 1e-4 of its value, or of 2500 J/mol, about RT, where that's larger, and
// NRTL's α by 1e-4. Its bubble points' rounding, 1e-12, over the step is
// then about 1e-8 of the Jacobian, and a central difference's bias about
// as much, which moves the answer by as little of its standard error but
// stops no convergence. With 1e-6 the rounding was 1e-6 of J, and kept the
// published methanol–water set's fit of NRTL's α from its gradient test.
const energy_step = 1e-4

const energy_floor = 2500.0

const alpha_step = 1e-4

// step_of is a parameter's central-difference step at its value v.
fn step_of(par flowsheet.FitParameter, v f64) f64 {
	return match par.part {
		.pair_energy { energy_step * math.max(math.abs(v), energy_floor) }
		.alpha { alpha_step }
		else { kij_step }
	}
}

// GammaPairs is a γ–φ fluid whose activity model's pairs a fit varies (M6's
// D16): the fluid as built, its NRTL or UNIQUAC, the solvents' names, each
// solvent pair's expressions, and UNIQUAC's r and q.
struct GammaPairs {
	base   thermo.GammaPhi
	liquid string
	names  []string
	uses   []components.PairUse
	r      []f64
	q      []f64
}

// pair_text is τ as an expression of a pair's energy A, J/mol: NRTL's A/RT,
// UNIQUAC's exp(−A/RT). V writes a double as the shortest decimal that
// reads back to it, so the expression's A is the fit's.
pub fn pair_text(liquid string, a f64) string {
	return if liquid == 'NRTL' {
		'${a} J/mol / (R*T)'
	} else {
		'exp(${-a} J/mol / (R*T))'
	}
}

// model is the fluid with θ's pairs: each fitted energy's τ and α in place
// of the expressions it was built with, the rest as they were.
fn (g &GammaPairs) model(parameters []flowsheet.FitParameter, theta []f64) !thermo.PhaseModel {
	mut texts := g.uses.map(it.texts.clone())
	for q, par in parameters {
		ia, ib := g.base.solvents.index(par.i), g.base.solvents.index(par.j)
		for k, u in g.uses {
			forward := u.i == ia && u.j == ib
			if !forward && !(u.i == ib && u.j == ia) {
				continue
			}
			if par.part == .alpha {
				texts[k][2] = theta[q].str()
			} else {
				texts[k][if forward { 0 } else { 1 }] = pair_text(g.liquid, theta[q])
			}
		}
	}
	liquid := if g.liquid == 'NRTL' {
		mut pairs := []activity.NrtlPair{}
		for k, u in g.uses {
			pairs << activity.NrtlPair{
				i:      u.i
				j:      u.j
				tau_ij: activity.plain(texts[k][0], 'τ')!
				tau_ji: activity.plain(texts[k][1], 'τ')!
				alpha:  activity.plain(texts[k][2], 'α')!
			}
		}
		activity.Model(activity.new_nrtl(g.names, pairs)!)
	} else {
		mut pairs := []activity.UniquacPair{}
		for k, u in g.uses {
			pairs << activity.UniquacPair{
				i:      u.i
				j:      u.j
				tau_ij: activity.plain(texts[k][0], 'τ')!
				tau_ji: activity.plain(texts[k][1], 'τ')!
			}
		}
		activity.Model(activity.new_uniquac(g.names, g.r, g.q, pairs)!)
	}
	return thermo.PhaseModel(thermo.new_gamma_phi(thermo.GammaPhi{
		...g.base
		liquid: liquid
	})!)
}

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
// weights, the model the file's cubic with the fitted kᵢⱼ, or its γ–φ fluid
// with the fitted pairs (M6's D16). Its Jacobian is by central differences,
// each kᵢⱼ and α ± 1e-6, each energy ± 1e-6 of itself (step_of).
struct VleProblem {
	nc         int
	cubic      ?thermo.Cubic
	gamma      ?GammaPairs
	parameters []flowsheet.FitParameter
	rows       []VleRow
	p_weight   flowsheet.Weight
	y_weights  []?flowsheet.Weight // one a component, none where y isn't measured
	y_names    []string            // each component's y column, '' for none
mut:
	flashes int
}

// model_at is the model at θ: the cubic with its fitted kᵢⱼ, or the γ–φ
// fluid with its fitted pairs.
fn (p &VleProblem) model_at(theta []f64) !thermo.PhaseModel {
	if base := p.cubic {
		nc := p.nc
		mut kij := base.kij.clone()
		for q, par in p.parameters {
			kij[par.i * nc + par.j] = theta[q]
			kij[par.j * nc + par.i] = theta[q]
		}
		return thermo.by_roots(thermo.Cubic{
			...base
			kij: kij
		})
	}
	g := p.gamma or { return error('the fit has no model') }
	return g.model(p.parameters, theta)!
}

// residuals are the weighted residuals at θ, each row's P first and then
// each measured y.
fn (mut p VleProblem) residuals(theta []f64) ![]f64 {
	nc := p.nc
	model := p.model_at(theta)!
	mut out := []f64{}
	for row in p.rows {
		p.flashes++
		b := flash.bubble_p(model, row.t, row.x) or {
			return error('no bubble point at T = ${row.t} K and x = ${row.x}: ${err.msg()}')
		}
		if b.phases.len != 2 {
			return error('the bubble point at T = ${row.t} K and x = ${row.x} has ${b.phases.len} phases')
		}
		y := b.phases.filter(it.label == .vapor)[0].state.n
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
		h := step_of(p.parameters[q], theta[q])
		up[q] += h
		down[q] -= h
		ru := p.residuals(up)!
		rd := p.residuals(down)!
		for i in 0 .. r0.len {
			j[i * n + q] = (ru[i] - rd[i]) / (up[q] - down[q])
		}
	}
}

// fit_vle fits a file's kᵢⱼ to bubble points (D15), or its activity model's
// pairs (M6's D16): T, P, the liquid's composition and the vapour's, P and
// each y weighted by its error, which the file must give, since they're of
// different kinds.
fn fit_vle(file flowsheet.FitFile, plan flowsheet.FitPlan, data Data) Fitted {
	fl := file.fluid or { return failed(plan, 'the file has no [model]') }
	mut cubic := ?thermo.Cubic(none)
	mut gamma := ?GammaPairs(none)
	if plan.parameters[0].part == .kij {
		eos := fl.helmholtz() or { return failed(plan, "kij is a cubic's, and the model is not one") }
		if eos !is thermo.Cubic {
			return failed(plan, "kij is a cubic's, and the model is not one")
		}
		cubic = eos as thermo.Cubic
	} else {
		g := file.gamma_phi or {
			return failed(plan, "A and alpha are an activity model's pairs, and the model is not one")
		}
		mut r, mut q := []f64{}, []f64{}
		if file.liquid == 'UNIQUAC' {
			for i in g.model.solvents {
				r << g.records[i].parameters.uniquac_r or { 0.0 }
				q << g.records[i].parameters.uniquac_q or { 0.0 }
			}
		}
		gamma = GammaPairs{
			base:   g.model
			liquid: file.liquid
			names:  g.model.solvents.map(g.records[it].name)
			uses:   g.pairs
			r:      r
			q:      q
		}
	}
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
		nc:         nc
		cubic:      cubic
		gamma:      gamma
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
	nc := p.nc
	model := p.model_at(theta)!
	mut out := []Residual{}
	mut k := 0
	for r, row in p.rows {
		b := flash.bubble_p(model, row.t, row.x)!
		vapour := b.phases.filter(it.label == .vapor)[0]
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
						model:    vapour.state.n[i]
						weighted: weighted[k]
					}
					k++
				}
			}
		}
	}
	return out
}
