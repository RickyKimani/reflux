module numeric

import math

// Gauss–Kronrod 7–15 nodes and weights on [−1, 1], from QUADPACK (Piessens
// et al., 1983). xgk holds the 15 Kronrod nodes' non-negative half,
// largest first. The 7-point Gauss rule uses every other one (indices 1,
// 3, 5 and the centre).
const xgk = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
	0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
	0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
	0.207784955007898467600689403773245, 0.000000000000000000000000000000000]!
const wgk = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
	0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
	0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
	0.204432940075298892414161999234649, 0.209482141084727828012999174891714]!
const wg = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
	0.381830050505118944950369775488975, 0.417959183673469387755102040816327]!

// The most panels integrate will split [a, b] into before giving up: about
// 30 000 evaluations of f. A smooth correlation needs one or two.
const max_panels = 2000

// Panel is one piece of [a, b] with its Gauss–Kronrod estimates.
struct Panel {
	a      f64
	b      f64
	value  f64 // the Kronrod estimate of ∫ f over [a, b]
	err    f64 // its estimated error
	resabs f64 // ∫ |f|, for the roundoff limit
}

// integrate returns ∫ f(x) dx from a to b, to a relative accuracy of about
// rtol, by globally adaptive Gauss–Kronrod 7–15 quadrature (QUADPACK's QAG).
// It keeps a list of panels and repeatedly halves the one with the largest
// error estimate, until the estimates add up to less than rtol·|result|.
//
// Floating-point roundoff limits the accuracy to about 1e-14 relative to
// ∫ |f|; a smaller rtol gets that instead. integrate fails, with an error,
// if f isn't finite somewhere, or if max_panels aren't enough.
pub fn integrate(f fn (f64) f64, a f64, b f64, rtol f64) !f64 {
	if a == b {
		return 0.0
	}
	mut panels := [gk15(f, a, b)]
	for {
		mut total := 0.0
		mut err := 0.0
		mut resabs := 0.0
		mut worst := 0
		for i, p in panels {
			total += p.value
			err += p.err
			resabs += p.resabs
			if p.err > panels[worst].err {
				worst = i
			}
		}
		if !math.is_finite(total) {
			return error('integrate: integrand is not finite on [${a}, ${b}]')
		}
		// Each panel's estimate is at least 50ε·∫|f| (see gk15), so 100ε
		// leaves room for the panels to be roundoff-limited.
		if err <= math.max(rtol * math.abs(total), 100.0 * f64_epsilon * resabs) {
			return total
		}
		if panels.len >= max_panels {
			return error('integrate: no convergence on [${a}, ${b}] with ${max_panels} panels')
		}
		p := panels[worst]
		m := 0.5 * (p.a + p.b)
		if m == p.a || m == p.b {
			// The panel is too narrow to halve: f is singular near here.
			return error('integrate: no convergence near x = ${p.a}')
		}
		panels[worst] = gk15(f, p.a, m)
		panels << gk15(f, m, p.b)
	}
	return error('unreachable')
}

// gk15 applies both rules to [a, b]. The error estimate is QUADPACK's: the
// difference between the rules, scaled down because the Kronrod estimate
// is far more accurate than the Gauss one it's compared with, and never
// below what roundoff allows.
fn gk15(f fn (f64) f64, a f64, b f64) Panel {
	c := 0.5 * (a + b)
	h := 0.5 * (b - a)
	mut fv1 := [7]f64{}
	mut fv2 := [7]f64{}
	fc := f(c)
	mut kronrod := fc * wgk[7]
	mut gauss := fc * wg[3]
	mut resabs := math.abs(kronrod)
	for j in 0 .. 7 {
		dx := h * xgk[j]
		fv1[j] = f(c - dx)
		fv2[j] = f(c + dx)
		kronrod += wgk[j] * (fv1[j] + fv2[j])
		resabs += wgk[j] * (math.abs(fv1[j]) + math.abs(fv2[j]))
		if j % 2 == 1 {
			gauss += wg[j / 2] * (fv1[j] + fv2[j])
		}
	}
	// resasc approximates ∫ |f − mean f|, the scale the error is judged on.
	mean := 0.5 * kronrod
	mut resasc := wgk[7] * math.abs(fc - mean)
	for j in 0 .. 7 {
		resasc += wgk[j] * (math.abs(fv1[j] - mean) + math.abs(fv2[j] - mean))
	}
	ah := math.abs(h)
	resasc *= ah
	resabs *= ah
	mut err := math.abs((kronrod - gauss) * h)
	if resasc != 0.0 && err != 0.0 {
		err = resasc * math.min(1.0, math.pow(200.0 * err / resasc, 1.5))
	}
	if resabs > math.smallest_non_zero_f64 / (50.0 * f64_epsilon) {
		err = math.max(50.0 * f64_epsilon * resabs, err)
	}
	return Panel{
		a:      a
		b:      b
		value:  kronrod * h
		err:    err
		resabs: resabs
	}
}
