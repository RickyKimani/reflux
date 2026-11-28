module capi

// A phase's transport at many points (M5b's D17): its viscosity, thermal
// conductivity and surface tension, each with its own status, the phase's
// volume from the fluid's model at its label's root. Each value is the V
// API's, transport.Model.phase at T, P, the mole fractions and the molar
// volume, so a V program calling it gets the same bits (criterion 11).
import math
import numeric
import thermo
import transport

// REFLUX_POINT_NO_VALUE: a property with no value at a point that has a
// phase, its reason in err.
const point_no_value = i32(9)

// reflux_transport's values, in reflux.h's order: REFLUX_TRANSPORT_VISCOSITY
// and so on.
const transport_names = ['viscosity', 'conductivity', 'surface_tension']
const n_transport = transport_names.len

// Labels reflux_transport takes, and the volume root each is had at.
const transport_labels = {
	i32(0): thermo.PhaseLabel.liquid
	1:      .vapor
}

// TransportPoint is one point of reflux_transport: each property's value
// and status, and why any has none (but a vapour's surface tension).
struct TransportPoint {
	values   []f64
	statuses []i32
	reasons  string
}

// transport_point is one point of reflux_transport. A point with no phase
// fails as a whole, every status its reason's.
fn transport_point(f thermo.Fluid, tm transport.Model, label thermo.PhaseLabel, t f64, p f64, z []f64) !TransportPoint {
	check_tp(t, p)!
	check_z(z)!
	root := if label == .vapor { thermo.Root.vapor } else { thermo.Root.liquid }
	st := f.model.phase(t, p, z, root) or { return point_error(point_no_volume_root, err.msg()) }
	n := numeric.sum(z)
	tr := tm.phase(label, t, p, z.map(it / n), st.v / n)
	mut values := []f64{len: n_transport, init: math.nan()}
	mut statuses := []i32{len: n_transport, init: point_no_value}
	for k, r in [tr.viscosity, tr.conductivity, tr.surface_tension] {
		if r.ok() {
			values[k] = r.value
			statuses[k] = point_ok
		}
	}
	return TransportPoint{
		values:   values
		statuses: statuses
		reasons:  tr.reasons(label).join('; ')
	}
}

// reflux_transport evaluates a phase's transport at many points.
@[export: 'reflux_transport']
pub fn reflux_transport(fluid i64, label i32, n i64, t &f64, t_len i64, p &f64, p_len i64, z &f64, z_rows i64, out &f64, status &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	lb := transport_labels[label] or {
		return fail(bad_argument, buf, buf_len, 'no label ${label}: REFLUX_LIQUID (0) or REFLUX_VAPOR (1)')
	}
	pts := check_points(n, [voidptr(t), p, z, out, status], ['t', 'p', 'z', 'out', 'status']) or {
		return failed(err, buf, buf_len)
	}
	pts.check(t_len, 't') or { return failed(err, buf, buf_len) }
	pts.check(p_len, 'p') or { return failed(err, buf, buf_len) }
	pts.check(z_rows, 'z') or { return failed(err, buf, buf_len) }
	tm := s.transport() or { return fail(data_error, buf, buf_len, err.msg()) }
	nc := s.fluid.model.components().len
	mut reasons := Reasons{}
	for i in i64(0) .. n {
		zr := row(z_rows, i) * nc
		zi := []f64{len: nc, init: unsafe { z[zr + index] }}
		ti := unsafe { t[row(t_len, i)] }
		pi := unsafe { p[row(p_len, i)] }
		k := int(i)
		res := transport_point(s.fluid, tm, lb, ti, pi, zi) or {
			reasons.add(err)
			unsafe {
				for m in 0 .. n_transport {
					out[m * int(n) + k] = math.nan()
					status[m * int(n) + k] = status_of(err)
				}
			}
			continue
		}
		if res.reasons != '' {
			reasons.add(error(res.reasons))
		}
		unsafe {
			for m in 0 .. n_transport {
				out[m * int(n) + k] = res.values[m]
				status[m * int(n) + k] = res.statuses[m]
			}
		}
	}
	write_message(buf, buf_len, reasons.first)
	return ok
}
