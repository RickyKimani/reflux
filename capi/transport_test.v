module capi

// reflux_transport (M5b's D17), called in process: every value is the V
// API's, transport.Model.phase at the phase's volume root, bit for bit,
// with its own status; a point without a phase fails as a whole; and err
// holds the first point's reasons.
import components
import math
import thermo
import water

fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

struct TransportOut {
	code   i32
	out    []f64
	status []i32
	reason string
}

fn call_transport(fluid i64, label i32, nc int, t []f64, p []f64, z []f64) TransportOut {
	n := math.max(t.len, math.max(p.len, z.len / nc))
	mut out := []f64{len: n_transport * n, init: 7.0}
	mut status := []i32{len: n_transport * n, init: 99}
	mut buf := []u8{len: 1024, init: 0xff}
	code := reflux_transport(fluid, label, n, unsafe { &f64(t.data) }, t.len, unsafe { &f64(p.data) },
		p.len, unsafe { &f64(z.data) }, z.len / nc, unsafe { &f64(out.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 1024)
	return TransportOut{code, out, status, c_string(buf)}
}

// want is the V API's transport for a point: the phase at its label's root,
// then the setup's transport model there.
fn want(s components.Setup, label thermo.PhaseLabel, t f64, p f64, z []f64) ![]f64 {
	root := if label == .vapor { thermo.Root.vapor } else { thermo.Root.liquid }
	st := s.fluid.model.phase(t, p, z, root)!
	mut n := 0.0
	for x in z {
		n += x
	}
	tr := s.transport()!.phase(label, t, p, z.map(it / n), st.v / n)
	return [tr.viscosity.value, tr.conductivity.value, tr.surface_tension.value]
}

fn same_bits(got f64, w f64) bool {
	return math.f64_bits(got) == math.f64_bits(w) || (math.is_nan(got) && math.is_nan(w))
}

fn test_each_value_is_the_v_apis() {
	db := components.load_default()!
	s := db.setup(.pr, ['n-butane', 'n-pentane'])!
	fl := add(FluidItem{s})
	defer {
		remove(fl) or {}
	}
	// A liquid at two states and two compositions, three of them per mole
	// (amounts summing to 3), then the vapour at one composition.
	t := [300.0, 320.0]
	p := [5e5, 8e5]
	z := [0.3, 0.7, 1.2, 1.8]
	r := call_transport(fl, liquid, 2, t, p, z)
	assert r.code == ok && r.reason == '', r.reason
	for i in 0 .. 2 {
		w := want(s, .liquid, t[i], p[i], z[i * 2..i * 2 + 2])!
		for k in 0 .. n_transport {
			assert r.status[k * 2 + i] == point_ok, '${i} ${k}'
			assert same_bits(r.out[k * 2 + i], w[k]), '${i} ${k}: ${r.out[k * 2 + i]} against ${w[k]}'
		}
	}
	v := call_transport(fl, vapor, 2, [400.0], [1e5], [0.5, 0.5])
	assert v.code == ok
	// A vapour's surface tension has no value, and isn't a reason.
	assert v.reason == '', v.reason
	w := want(s, .vapor, 400.0, 1e5, [0.5, 0.5])!
	for k in 0 .. 2 {
		assert v.status[k] == point_ok && same_bits(v.out[k], w[k])
	}
	assert v.status[2] == point_no_value && math.is_nan(v.out[2])
	// Per mole: 3 mol of the vapour has 1 mol's values, its volume divided by
	// its moles (a gas's conductivity takes Stiel and Thodos' residual at it).
	v3 := call_transport(fl, vapor, 2, [400.0], [1e5], [1.5, 1.5])
	w3 := want(s, .vapor, 400.0, 1e5, [1.5, 1.5])!
	for k in 0 .. 2 {
		assert same_bits(v3.out[k], w3[k]) && math.abs(v3.out[k] - w[k]) <= 1e-12 * w[k], '${k}'
	}
	// At 300 K and 1 bar the mixture has a liquid root and a vapour root, and
	// each label takes its own.
	for label in [liquid, vapor] {
		lb := if label == vapor { thermo.PhaseLabel.vapor } else { thermo.PhaseLabel.liquid }
		two := call_transport(fl, label, 2, [300.0], [1e5], [0.5, 0.5])
		w2 := want(s, lb, 300.0, 1e5, [0.5, 0.5])!
		assert same_bits(two.out[0], w2[0]) && same_bits(two.out[1], w2[1]), '${label}'
	}
	gas := want(s, .vapor, 300.0, 1e5, [0.5, 0.5])!
	liq := want(s, .liquid, 300.0, 1e5, [0.5, 0.5])!
	assert liq[0] > 10.0 * gas[0], 'the liquid root: ${liq[0]} against ${gas[0]}'
}

fn test_iapws95_water_is_m5as_functions() {
	db := components.load_default()!
	s := db.iapws95('water')!
	fl := add(FluidItem{s})
	defer {
		remove(fl) or {}
	}
	// One mole, and three, whose density is the same.
	for amount in [1.0, 3.0] {
		r := call_transport(fl, liquid, 1, [300.0], [1e5], [amount])
		assert r.code == ok && r.status == [point_ok, point_ok, point_ok], r.reason
		st := s.fluid.model.phase(300.0, 1e5, [amount], .liquid)!
		rho := water.molar_mass * amount / st.v
		assert r.out[0] == water.viscosity(300.0, rho)!, '${amount}'
		assert r.out[1] == water.thermal_conductivity(300.0, rho)!, '${amount}'
		assert r.out[2] == water.surface_tension(300.0)!
		assert math.abs(rho - 996.5) <= 0.1, '${rho} kg/m³'
	}
}

fn test_a_point_without_a_phase_or_values_says_why() {
	db := components.load_default()!
	s := db.setup(.pr, ['methane', 'n-butane'])!
	fl := add(FluidItem{s})
	defer {
		remove(fl) or {}
	}
	// A bad T fails as a whole; a liquid at 250 K with 30% methane, above
	// its Tc, has no value (D9); the next is a liquid of n-butane alone.
	r := call_transport(fl, liquid, 2, [-1.0, 250.0, 250.0], [20e5], [0.3, 0.7, 0.3, 0.7, 0.0,
		1.0])
	assert r.code == ok
	for k in 0 .. n_transport {
		assert r.status[k * 3] == point_bad_input && math.is_nan(r.out[k * 3])
		assert r.status[k * 3 + 1] == point_no_value && math.is_nan(r.out[k * 3 + 1])
		assert r.status[k * 3 + 2] == point_ok
	}
	assert r.reason.starts_with('T = -1.0 K'), r.reason
	// The second point alone gives its own reasons.
	one := call_transport(fl, liquid, 2, [250.0], [20e5], [0.3, 0.7])
	assert one.reason == "no viscosity, thermal conductivity or surface tension: Methane, above its Tc, is 0.3 of the liquid's moles, more than 0.1 (D9)", one.reason
}

fn test_the_call_fails_only_for_what_is_wrong_with_it() {
	db := components.load_default()!
	s := db.setup(.pr, ['methane'])!
	fl := add(FluidItem{s})
	mut out := []f64{len: 9}
	mut status := []i32{len: 9}
	mut buf := []u8{len: 256}
	t, p, z := [300.0], [1e5], [1.0]
	assert reflux_transport(fl, 2, 1, unsafe { &f64(t.data) }, 1, unsafe { &f64(p.data) },
		1, unsafe { &f64(z.data) }, 1, unsafe { &f64(out.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 256) == bad_argument
	assert c_string(buf).contains('no label 2')
	assert reflux_transport(fl, liquid, 3, unsafe { &f64(t.data) }, 2, unsafe { &f64(p.data) },
		1, unsafe { &f64(z.data) }, 1, unsafe { &f64(out.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 256) == bad_argument
	remove(fl)!
	assert reflux_transport(fl, liquid, 1, unsafe { &f64(t.data) }, 1, unsafe { &f64(p.data) },
		1, unsafe { &f64(z.data) }, 1, unsafe { &f64(out.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 256) == bad_handle
	// No points: nothing read or written, NULLs and all.
	assert reflux_transport(add(FluidItem{s}), liquid, 0, unsafe { nil }, 0, unsafe { nil },
		0, unsafe { nil }, 0, unsafe { nil }, unsafe { nil }, &char(buf.data), 256) == ok
}
