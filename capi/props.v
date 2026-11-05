module capi

// Many points at once (M10a step 3): a phase's properties, and each
// compound's vapour pressure. A point that fails gets its status and NaN,
// and the others are as they'd be alone (D7). Each point is the V API's
// answer, divided by the moles where it's extensive, so that a V program
// calling the V API gets the same bits (criterion 3).
import math
import numeric
import thermo

// A point's status (reflux.h). The flash's kinds, from bad input to
// outside the envelope, are flash.Failure's in its order.
const point_ok = i32(0)
const point_bad_input = i32(1)
const point_no_convergence = i32(2)
const point_no_volume_root = i32(3)
const point_more_phases = i32(4)
const point_phases_swapped = i32(5)
const point_outside_envelope = i32(6)
const point_supercritical = i32(7)
const point_not_finite = i32(8)

// A phase's label (reflux.h).
const no_phase = i32(-1)
const liquid = i32(0)
const vapor = i32(1)

// Volume roots, as reflux.h numbers them.
const roots = {
	i32(0): thermo.Root.liquid
	1:      .vapor
	2:      .stable
}

// The properties reflux_props writes, in reflux.h's order: REFLUX_PROP_Z
// and so on.
const prop_names = ['z', 'v', 'h', 's', 'g', 'u', 'a', 'cp', 'cv', 'w', 'mu_jt']
const n_props = prop_names.len

// PointError is why a point has no answer.
struct PointError {
	Error
	status i32
	reason string
}

fn (e PointError) msg() string {
	return e.reason
}

fn point_error(status i32, reason string) IError {
	return PointError{
		status: status
		reason: reason
	}
}

// Points is the shape of a many-point call: n points, and inputs with 1
// entry or n, read through row.
struct Points {
	n i64
}

// check says whether an input's count is one the call can use.
fn (pts Points) check(len i64, what string) ! {
	if len != 1 && len != pts.n {
		return bad('${what} has ${len} entries: 1, or n = ${pts.n}')
	}
}

// row is the entry of an input of len entries that point i uses.
fn row(len i64, i i64) int {
	return if len == 1 { 0 } else { int(i) }
}

fn not_null(p voidptr, what string) ! {
	if p == unsafe { nil } {
		return bad('${what} is NULL')
	}
}

// Reasons keeps the first failed point's reason, for err. The status array
// says which point that is, and the same call for that point alone gives its
// reason, which is how a caller asks for one point's (D7).
struct Reasons {
mut:
	first string
	any   bool
}

fn (mut r Reasons) add(e IError) {
	if !r.any {
		r.first = e.msg()
		r.any = true
	}
}

fn status_of(e IError) i32 {
	return if e is PointError { e.status } else { point_bad_input }
}

// check_tp refuses a temperature or pressure that isn't finite and above
// zero. thermo's own check, `t <= 0.0`, lets NaN and infinity through.
fn check_tp(t f64, p f64) ! {
	check_t(t)!
	check_p(p)!
}

fn check_t(t f64) ! {
	if !math.is_finite(t) || t <= 0.0 {
		return point_error(point_bad_input, 'T = ${t} K: a temperature must be finite and above 0')
	}
}

fn check_p(p f64) ! {
	if !math.is_finite(p) || p <= 0.0 {
		return point_error(point_bad_input, 'P = ${p} Pa: a pressure must be finite and above 0')
	}
}

fn check_z(z []f64) ! {
	for j, x in z {
		if !math.is_finite(x) || x < 0.0 {
			return point_error(point_bad_input, 'z[${j}] = ${x}: amounts must be finite and not negative')
		}
	}
	total := numeric.sum(z)
	if total <= 0.0 || !math.is_finite(total) {
		return point_error(point_bad_input, 'z adds up to ${total}: it must be finite and above 0')
	}
}

// PhaseProps are one point's answers, as reflux_props writes them.
struct PhaseProps {
	values []f64 // the n_props properties, in reflux.h's order
	ln_phi []f64
	phase  i32
}

// phase_props is one point of reflux_props: the V API's state and
// properties, per mole.
fn phase_props(f thermo.Fluid, root thermo.Root, t f64, p f64, z []f64) !PhaseProps {
	check_tp(t, p)!
	check_z(z)!
	st := thermo.solve_tp(f.model, t, p, z, root) or {
		return point_error(point_no_volume_root, err.msg())
	}
	pr := f.props(st)
	n := numeric.sum(z)
	values := [st.z, st.v / n, pr.h / n, pr.s / n, pr.g / n, pr.u / n, pr.a / n, pr.cp / n, pr.cv / n,
		pr.w, pr.mu_jt]
	// An answer that isn't finite isn't one: far beyond the correlations'
	// range (above about 1e10 K their ideal-gas Cp overflows), or the
	// arithmetic's.
	for k, x in values {
		finite(x, prop_names[k], t, p)!
	}
	for j, x in st.ln_phi {
		finite(x, 'ln φ[${j}]', t, p)!
	}
	return PhaseProps{
		values: values
		ln_phi: st.ln_phi
		phase:  if thermo.phase_label(f.model, st) == .liquid { liquid } else { vapor }
	}
}

fn finite(x f64, what string, t f64, p f64) ! {
	if !math.is_finite(x) {
		return point_error(point_not_finite, '${what} is ${x} at T = ${t} K, P = ${p} Pa: the model gives no finite answer there')
	}
}

// reflux_props evaluates one phase at many points.
@[export: 'reflux_props']
pub fn reflux_props(fluid i64, root i32, n i64, t &f64, t_len i64, p &f64, p_len i64, z &f64, z_rows i64, out &f64, ln_phi &f64, phase &i32, status &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	rt := roots[root] or {
		return fail(bad_argument, buf, buf_len, 'no root ${root}: REFLUX_ROOT_LIQUID (0), REFLUX_ROOT_VAPOR (1) or REFLUX_ROOT_STABLE (2)')
	}
	pts := check_points(n, [voidptr(t), p, z, out, ln_phi, phase, status], [
		't',
		'p',
		'z',
		'out',
		'ln_phi',
		'phase',
		'status',
	]) or { return failed(err, buf, buf_len) }
	pts.check(t_len, 't') or { return failed(err, buf, buf_len) }
	pts.check(p_len, 'p') or { return failed(err, buf, buf_len) }
	pts.check(z_rows, 'z') or { return failed(err, buf, buf_len) }
	nc := s.fluid.model.components().len
	mut reasons := Reasons{}
	for i in i64(0) .. n {
		zr := row(z_rows, i) * nc
		zi := []f64{len: nc, init: unsafe { z[zr + index] }}
		ti := unsafe { t[row(t_len, i)] }
		pi := unsafe { p[row(p_len, i)] }
		k := int(i)
		res := phase_props(s.fluid, rt, ti, pi, zi) or {
			reasons.add(err)
			unsafe {
				for m in 0 .. n_props {
					out[m * int(n) + k] = math.nan()
				}
				for j in 0 .. nc {
					ln_phi[k * nc + j] = math.nan()
				}
				phase[k] = no_phase
				status[k] = status_of(err)
			}
			continue
		}
		unsafe {
			for m in 0 .. n_props {
				out[m * int(n) + k] = res.values[m]
			}
			for j in 0 .. nc {
				ln_phi[k * nc + j] = res.ln_phi[j]
			}
			phase[k] = res.phase
			status[k] = point_ok
		}
	}
	write_message(buf, buf_len, reasons.first)
	return ok
}

// check_points checks a many-point call's count, and its arrays when there
// are points to write.
fn check_points(n i64, arrays []voidptr, names []string) !Points {
	if n < 0 {
		return bad('n is ${n}')
	}
	if n > 0 {
		for i, a in arrays {
			not_null(a, names[i])!
		}
	}
	return Points{n}
}

// pure_model is a cubic's j'th compound alone, with its own α: what its
// vapour pressure is computed from.
fn pure_model(m thermo.Model, j int) !thermo.Cubic {
	if m is thermo.Cubic {
		return thermo.new_cubic(m.kind, [m.comps[j]], alphas: [m.alphas[j]])!
	}
	return error('vapour pressure is for cubic equations of state')
}

// saturation is one point of reflux_psat.
fn saturation(pure thermo.Cubic, t f64) !f64 {
	c := pure.comps[0]
	check_t(t)!
	if t >= c.tc {
		return point_error(point_supercritical, '${c.name}: T = ${t} K is at or above its Tc, ${c.tc} K, so it has no vapour pressure')
	}
	// saturation_pressure returns only once it has converged on finite
	// volumes, so its answer is always finite.
	return thermo.saturation_pressure(pure, t) or {
		return point_error(point_no_convergence, err.msg())
	}
}

// reflux_psat writes one compound's vapour pressure at many temperatures.
@[export: 'reflux_psat']
pub fn reflux_psat(fluid i64, component i32, n i64, t &f64, out &f64, status &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	nc := s.records.len
	if component < 0 || component >= nc {
		return fail(bad_argument, buf, buf_len, 'component ${component}: the fluid has ${nc}, from 0')
	}
	check_points(n, [voidptr(t), out, status], ['t', 'out', 'status']) or {
		return failed(err, buf, buf_len)
	}
	pure := pure_model(s.fluid.model, int(component)) or {
		return fail(data_error, buf, buf_len, err.msg())
	}
	mut reasons := Reasons{}
	for i in 0 .. int(n) {
		ps := saturation(pure, unsafe { t[i] }) or {
			reasons.add(err)
			unsafe {
				out[i] = math.nan()
				status[i] = status_of(err)
			}
			continue
		}
		unsafe {
			out[i] = ps
			status[i] = point_ok
		}
	}
	write_message(buf, buf_len, reasons.first)
	return ok
}
