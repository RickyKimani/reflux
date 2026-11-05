module capi

// Stability and flashes at many points (M10a step 4), in one layout for every
// flash (D6): the whole's β, T, P, H, S and V; each phase's fraction, V, H, S
// and composition, in the place of its label; and the number of phases. Each
// point is the V API's answer, with H, S and V summed over the phases in the
// order the flash gives them, as flowsheet.Stream does.
import flash
import math
import numeric
import thermo

// Flashes, as reflux.h numbers them: what a point's two values are.
const flash_pt = i32(0)

// A flash's values and a phase's, in reflux.h's order.
const flash_names = ['beta', 't', 'p', 'h', 's', 'v']
const phase_names = ['fraction', 'v', 'h', 's']
const n_flash = flash_names.len
const n_phase = phase_names.len

// flash_status is the point status of a flash's failure: its kind's, in
// flash.Failure's order, which reflux.h's numbers follow from 1.
fn flash_status(e IError) i32 {
	if e is flash.FlashError {
		return match e.kind {
			.bad_input { point_bad_input }
			.no_convergence { point_no_convergence }
			.no_volume_root { point_no_volume_root }
			.more_phases { point_more_phases }
			.phases_swapped { point_phases_swapped }
			.outside_envelope { point_outside_envelope }
		}
	}
	return status_of(e)
}

// StabilityPoint is one point of reflux_stability.
struct StabilityPoint {
	stable bool
	tm     f64
	w      []f64
}

// stability_point tests one feed with the V API's defaults: every trial, so
// that tm is the lowest any reached.
fn stability_point(m thermo.Model, t f64, p f64, z []f64) !StabilityPoint {
	check_tp(t, p)!
	check_z(z)!
	st := flash.stability(m, t, p, z)!
	mut best := -1
	for i, tr in st.trials {
		if !tr.trivial && (best < 0 || tr.tm < st.trials[best].tm) {
			best = i
		}
	}
	if best < 0 {
		total := numeric.sum(z)
		return StabilityPoint{st.stable, 0.0, z.map(it / total)}
	}
	tr := st.trials[best]
	return StabilityPoint{st.stable, tr.tm, tr.w}
}

// reflux_stability tests many feeds' stability.
@[export: 'reflux_stability']
pub fn reflux_stability(fluid i64, n i64, t &f64, t_len i64, p &f64, p_len i64, z &f64, z_rows i64, stable &i32, tm &f64, w &f64, status &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	pts := check_points(n, [voidptr(t), p, z, stable, tm, w, status], ['t', 'p', 'z', 'stable',
		'tm', 'w', 'status']) or { return failed(err, buf, buf_len) }
	pts.check(t_len, 't') or { return failed(err, buf, buf_len) }
	pts.check(p_len, 'p') or { return failed(err, buf, buf_len) }
	pts.check(z_rows, 'z') or { return failed(err, buf, buf_len) }
	nc := s.records.len
	mut reasons := Reasons{}
	for i in 0 .. int(n) {
		zr := row(z_rows, i) * nc
		zi := []f64{len: nc, init: unsafe { z[zr + index] }}
		res := stability_point(s.fluid.model, unsafe { t[row(t_len, i)] }, unsafe { p[row(p_len, i)] },
			zi) or {
			reasons.add(err)
			unsafe {
				stable[i] = -1
				tm[i] = math.nan()
				for j in 0 .. nc {
					w[i * nc + j] = math.nan()
				}
				status[i] = flash_status(err)
			}
			continue
		}
		unsafe {
			stable[i] = if res.stable { 1 } else { 0 }
			tm[i] = res.tm
			for j in 0 .. nc {
				w[i * nc + j] = res.w[j]
			}
			status[i] = point_ok
		}
	}
	write_message(buf, buf_len, reasons.first)
	return ok
}

// FlashPoint is one point of reflux_flash, as it writes it.
struct FlashPoint {
	values   []f64   // n_flash, in reflux.h's order
	phase    [][]f64 // by label: n_phase values each, NaN for an absent phase
	x        [][]f64 // by label: the composition, NaN for an absent phase
	n_phases int
}

// flash_at runs one point's flash.
fn flash_at(spec i32, f thermo.Fluid, a f64, b f64, z []f64) !flash.FlashResult {
	return match spec {
		flash_pt {
			check_tp(a, b)!
			check_z(z)!
			flash.pt(f.model, a, b, z)!
		}
		else {
			error('no flash ${spec}')
		}
	}
}

// flash_point is one point of reflux_flash: the flash, then each phase's
// properties and the whole's.
fn flash_point(spec i32, f thermo.Fluid, a f64, b f64, z []f64) !FlashPoint {
	res := flash_at(spec, f, a, b, z)!
	nc := z.len
	mut phase := [][]f64{len: 2, init: []f64{len: n_phase, init: math.nan()}}
	mut x := [][]f64{len: 2, init: []f64{len: nc, init: math.nan()}}
	mut h, mut s, mut v := 0.0, 0.0, 0.0
	for ph in res.phases {
		pr := f.props(ph.state)
		// A phase's state is one mole of it, so its properties are molar.
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
		v += ph.fraction * ph.state.v
		j := if ph.label == .liquid { 0 } else { 1 }
		phase[j] = [ph.fraction, ph.state.v, pr.h, pr.s]
		x[j] = ph.state.n.clone()
	}
	values := [res.beta, res.t, res.p, h, s, v]
	for k, val in values {
		finite(val, flash_names[k], res.t, res.p)!
	}
	// A phase's values need no check of their own: each enters the whole's
	// (its fraction, β or 1 − β), so one that isn't finite makes it so.
	return FlashPoint{values, phase, x, res.phases.len}
}

// reflux_flash flashes many feeds.
@[export: 'reflux_flash']
pub fn reflux_flash(fluid i64, spec i32, n i64, a &f64, a_len i64, b &f64, b_len i64, z &f64, z_rows i64, out &f64, phase &f64, x &f64, n_phases &i32, status &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	if spec != flash_pt {
		return fail(bad_argument, buf, buf_len, 'no flash ${spec}: REFLUX_FLASH_PT (0)')
	}
	pts := check_points(n, [voidptr(a), b, z, out, phase, x, n_phases, status], [
		'a',
		'b',
		'z',
		'out',
		'phase',
		'x',
		'n_phases',
		'status',
	]) or { return failed(err, buf, buf_len) }
	pts.check(a_len, 'a') or { return failed(err, buf, buf_len) }
	pts.check(b_len, 'b') or { return failed(err, buf, buf_len) }
	pts.check(z_rows, 'z') or { return failed(err, buf, buf_len) }
	nc := s.records.len
	nn := int(n)
	mut reasons := Reasons{}
	for i in 0 .. nn {
		zr := row(z_rows, i) * nc
		zi := []f64{len: nc, init: unsafe { z[zr + index] }}
		res := flash_point(spec, s.fluid, unsafe { a[row(a_len, i)] }, unsafe { b[row(b_len, i)] },
			zi) or {
			reasons.add(err)
			unsafe {
				for k in 0 .. n_flash {
					out[k * nn + i] = math.nan()
				}
				for k in 0 .. n_phase {
					for j in 0 .. 2 {
						phase[(k * nn + i) * 2 + j] = math.nan()
					}
				}
				for c in 0 .. 2 * nc {
					x[i * 2 * nc + c] = math.nan()
				}
				n_phases[i] = 0
				status[i] = flash_status(err)
			}
			continue
		}
		unsafe {
			for k in 0 .. n_flash {
				out[k * nn + i] = res.values[k]
			}
			for k in 0 .. n_phase {
				for j in 0 .. 2 {
					phase[(k * nn + i) * 2 + j] = res.phase[j][k]
				}
			}
			for j in 0 .. 2 {
				for c in 0 .. nc {
					x[(i * 2 + j) * nc + c] = res.x[j][c]
				}
			}
			n_phases[i] = i32(res.n_phases)
			status[i] = point_ok
		}
	}
	write_message(buf, buf_len, reasons.first)
	return ok
}
