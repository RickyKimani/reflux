module capi

// Many points at once (M10a step 3), called in process: every point is the
// V API's answer, bit for bit; a point that fails gets its status, NaN and
// no phase, and leaves the others as they'd be alone; and the call itself
// fails only for what's wrong with it as a whole.
import components
import math
import thermo

// c_string reads a NUL-terminated string out of a buffer. (Each _test.v file
// is a program of its own.)
fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

// fluid_handle registers a fluid from the database, as reflux_fluid_new
// would, and returns its handle and its setup.
fn fluid_handle(kind thermo.CubicKind, names []string, opts components.FluidOptions) !(i64, components.Setup) {
	db := components.load_default()!
	s := db.setup(kind, names, opts)!
	return add(FluidItem{s}), s
}

// PropsCall is reflux_props' inputs, with counts of 1 or n.
struct PropsCall {
	root i32
	n    i64
	t    []f64
	p    []f64
	z    []f64 // z_rows × nc
}

struct PropsOut {
	code   i32
	out    []f64
	ln_phi []f64
	phase  []i32
	status []i32
	reason string
}

fn call_props(fluid i64, nc int, c PropsCall) PropsOut {
	n := int(c.n)
	mut out := []f64{len: n_props * n, init: 7.0}
	mut ln_phi := []f64{len: n * nc, init: 7.0}
	mut phase := []i32{len: n, init: 9}
	mut status := []i32{len: n, init: 9}
	mut buf := []u8{len: 1024, init: 0xff}
	code := reflux_props(fluid, c.root, c.n, unsafe { &f64(c.t.data) }, c.t.len, unsafe { &f64(c.p.data) },
		c.p.len, unsafe { &f64(c.z.data) }, c.z.len / nc, unsafe { &f64(out.data) }, unsafe { &f64(ln_phi.data) },
		unsafe { &i32(phase.data) }, unsafe { &i32(status.data) }, &char(buf.data), 1024)
	return PropsOut{code, out, ln_phi, phase, status, c_string(buf)}
}

// expect checks point k of a call against the V API for (t, p, z).
fn expect(f thermo.Fluid, root thermo.Root, t f64, p f64, z []f64, got PropsOut, k int, n int) ! {
	st := thermo.solve_tp(f.model, t, p, z, root)!
	pr := f.props(st)
	tot := z[0] + z[1]
	want := {
		'z':     st.z
		'v':     st.v / tot
		'h':     pr.h / tot
		's':     pr.s / tot
		'g':     pr.g / tot
		'u':     pr.u / tot
		'a':     pr.a / tot
		'cp':    pr.cp / tot
		'cv':    pr.cv / tot
		'w':     pr.w
		'mu_jt': pr.mu_jt
	}
	for m, name in prop_names {
		assert got.out[m * n + k] == want[name], '${name} at point ${k}'
	}
	assert got.ln_phi[k * 2..k * 2 + 2] == st.ln_phi
	label := if thermo.phase_label(f.model, st) == .liquid { liquid } else { vapor }
	assert got.phase[k] == label && got.status[k] == point_ok
}

fn test_every_point_is_the_v_apis() {
	fl, s := fluid_handle(.pr, ['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	// A liquid, a vapour, and a two-phase state's two roots; z as amounts
	// that don't add up to 1 at the last point.
	t := [250.0, 400.0, 300.0, 300.0]
	p := [50e5, 10e5, 20e5, 20e5]
	z := [0.2, 0.8, 0.9, 0.1, 0.5, 0.5, 1.5, 0.5]
	for root, want in roots {
		got := call_props(fl, 2, PropsCall{root, 4, t, p, z})
		assert got.code == ok && got.reason == ''
		for k in 0 .. 4 {
			expect(s.fluid, want, t[k], p[k], z[k * 2..k * 2 + 2], got, k, 4)!
		}
	}
	// The two-phase state's liquid and vapour roots differ, and stable is one.
	liq := call_props(fl, 2, PropsCall{0, 1, [300.0], [20e5], [0.5, 0.5]})
	vap := call_props(fl, 2, PropsCall{1, 1, [300.0], [20e5], [0.5, 0.5]})
	stable := call_props(fl, 2, PropsCall{2, 1, [300.0], [20e5], [0.5, 0.5]})
	assert liq.phase == [liquid] && vap.phase == [vapor]
	assert liq.out[0] < 0.3 && vap.out[0] > 0.7
	assert stable.out == liq.out || stable.out == vap.out
}

fn test_one_value_serves_every_point() {
	fl, _ := fluid_handle(.pr, ['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	full := call_props(fl, 2, PropsCall{2, 3, [300.0, 300.0, 300.0], [1e5, 10e5, 50e5], [
		0.3,
		0.7,
		0.3,
		0.7,
		0.3,
		0.7,
	]})
	short := call_props(fl, 2, PropsCall{2, 3, [300.0], [1e5, 10e5, 50e5], [0.3, 0.7]})
	assert full.code == ok && short.code == ok
	assert short.out == full.out && short.ln_phi == full.ln_phi && short.phase == full.phase
	assert short.status == [point_ok, point_ok, point_ok]
}

fn test_a_bad_point_leaves_the_others_alone() {
	fl, s := fluid_handle(.pr, ['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	inf := math.inf(1)
	cases := [
		[math.nan(), 10e5, 0.5, 0.5],
		[300.0, 0.0, 0.5, 0.5],
		[-inf, 10e5, 0.5, 0.5],
		[300.0, inf, 0.5, 0.5],
		[300.0, 10e5, -0.1, 1.1],
		[300.0, 10e5, 0.0, 0.0],
		[300.0, 10e5, math.nan(), 1.0],
		[300.0, 10e5, 1e308, 1e308],
	]
	reasons := [
		'T = nan K: a temperature must be finite and above 0',
		'P = 0.0 Pa: a pressure must be finite and above 0',
		'T = -inf K: a temperature must be finite and above 0',
		'P = +inf Pa: a pressure must be finite and above 0',
		'z[0] = -0.1: amounts must be finite and not negative',
		'z adds up to 0.0: it must be finite and above 0',
		'z[0] = nan: amounts must be finite and not negative',
		'z adds up to +inf: it must be finite and above 0',
	]
	for c, bad_point in cases {
		// The bad point between two good ones.
		t := [300.0, bad_point[0], 350.0]
		p := [10e5, bad_point[1], 5e5]
		z := [0.4, 0.6, bad_point[2], bad_point[3], 0.7, 0.3]
		got := call_props(fl, 2, PropsCall{2, 3, t, p, z})
		assert got.code == ok && got.status == [point_ok, point_bad_input, point_ok]
		assert got.reason == reasons[c]
		for m in 0 .. n_props {
			assert math.is_nan(got.out[m * 3 + 1])
		}
		assert got.ln_phi[2..4].all(math.is_nan(it)) && got.phase[1] == no_phase
		expect(s.fluid, .stable, 300.0, 10e5, [0.4, 0.6], got, 0, 3)!
		expect(s.fluid, .stable, 350.0, 5e5, [0.7, 0.3], got, 2, 3)!
		// Alone, the bad point says the same.
		alone := call_props(fl, 2, PropsCall{2, 1, [t[1]], [p[1]], z[2..4]})
		assert alone.status == [point_bad_input] && alone.reason == reasons[c]
	}
}

fn test_beyond_the_model_is_no_answer() {
	// Above about 1e10 K the ideal-gas Cp correlations overflow, so H and S
	// are NaN: not an answer. At 1e-300 Pa the vapour root is past what a
	// float holds. Points like these get their status, and the first one's
	// reason is err's.
	fl, _ := fluid_handle(.pr, ['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	got := call_props(fl, 2, PropsCall{2, 4, [300.0, 1e10, 300.0, 1e300], [1e5, 1e5, 1e-300, 1e5], [
		0.5,
		0.5,
	]})
	assert got.code == ok
	assert got.status == [point_ok, point_not_finite, point_no_volume_root, point_not_finite]
	assert got.reason == 'h is nan at T = 1e+10 K, P = 100000.0 Pa: the model gives no finite answer there'
	for k in 1 .. 4 {
		assert math.is_nan(got.out[k]) && got.phase[k] == no_phase
	}
	alone := call_props(fl, 2, PropsCall{2, 1, [300.0], [1e-300], [0.5, 0.5]})
	assert alone.status == [point_no_volume_root]
	assert alone.reason.starts_with('no volume root at T = 300.0 K')
}

fn test_a_call_that_cant_be_made() {
	fl, _ := fluid_handle(.pr, ['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	t, p, z := [300.0, 310.0], [1e5, 1e5], [0.5, 0.5]
	for c, why in {
		3: 'no root 3: REFLUX_ROOT_LIQUID (0), REFLUX_ROOT_VAPOR (1) or REFLUX_ROOT_STABLE (2)'
	} {
		got := call_props(fl, 2, PropsCall{i32(c), 2, t, p, z})
		assert got.code == bad_argument && got.reason == why
	}
	mut got := call_props(fl, 2, PropsCall{2, 3, t, p, z})
	assert got.code == bad_argument && got.reason == 't has 2 entries: 1, or n = 3'
	got = call_props(fl, 2, PropsCall{2, 2, t, [1e5, 1e5, 1e5], z})
	assert got.code == bad_argument && got.reason == 'p has 3 entries: 1, or n = 2'
	got = call_props(fl, 2, PropsCall{2, 3, [300.0], [1e5], [0.5, 0.5, 0.5, 0.5]})
	assert got.code == bad_argument && got.reason == 'z has 2 entries: 1, or n = 3'
	// Nothing out of place was written.
	assert got.out.all(it == 7.0) && got.status.all(it == 9)
	// No points: nothing to do, and nothing needed.
	mut buf := []u8{len: 256}
	assert reflux_props(fl, 2, 0, unsafe { nil }, 0, unsafe { nil }, 0, unsafe { nil }, 0,
		unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil }, &char(buf.data), 256) == ok
	assert reflux_props(fl, 2, -1, unsafe { nil }, 0, unsafe { nil }, 0, unsafe { nil }, 0,
		unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil }, &char(buf.data), 256) == bad_argument
	assert c_string(buf) == 'n is -1'
	assert reflux_props(fl, 2, 1, unsafe { &f64(t.data) }, 1, unsafe { &f64(p.data) }, 1,
		unsafe { &f64(z.data) }, 1, unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil },
		&char(buf.data), 256) == bad_argument
	assert c_string(buf) == 'out is NULL'
	// A database's handle isn't a fluid's.
	db := add(DatabaseItem{components.load_default()!})
	defer {
		remove(db) or {}
	}
	got = call_props(db, 2, PropsCall{2, 2, t, p, z})
	assert got.code == bad_handle && got.reason == 'handle ${db} is a database, not a fluid'
}

fn psat_of(fl i64, j i32, t []f64) (i32, []f64, []i32, string) {
	mut out := []f64{len: t.len, init: 7.0}
	mut status := []i32{len: t.len, init: 9}
	mut buf := []u8{len: 1024}
	code := reflux_psat(fl, j, t.len, unsafe { &f64(t.data) }, unsafe { &f64(out.data) },
		unsafe { &i32(status.data) }, &char(buf.data), 1024)
	return code, out, status, c_string(buf)
}

fn test_vapour_pressure_is_each_compounds_own() {
	// Mathias–Copeman α, so that the pure models must take the fluid's α.
	fl, s := fluid_handle(.pr, ['n-butane', 'methane'], alpha: .mathias_copeman)!
	defer {
		remove(fl) or {}
	}
	m := s.fluid.model as thermo.Cubic
	t := [150.0, 300.0, math.nan(), 425.0]
	// n-butane's Tc is 425.12 K, methane's 190.56 K.
	for j, want in [
		[point_ok, point_ok, point_bad_input, point_ok],
		[point_ok, point_supercritical, point_bad_input, point_supercritical],
	] {
		code, out, status, _ := psat_of(fl, i32(j), t)
		assert code == ok && status == want
		pure := thermo.new_cubic(.pr, [m.comps[j]], alphas: [m.alphas[j]])!
		for k, tk in t {
			if status[k] == point_ok {
				assert out[k] == thermo.saturation_pressure(pure, tk)!, '${k}, ${j}'
			} else {
				assert math.is_nan(out[k])
			}
		}
	}
	_, _, _, why := psat_of(fl, 1, t)
	assert why == 'Methane: T = 300.0 K is at or above its Tc, 190.56 K, so it has no vapour pressure'
	_, _, _, why2 := psat_of(fl, 1, [math.nan()])
	assert why2 == 'T = nan K: a temperature must be finite and above 0'
	// At Tc itself there's none; and thermo's successive substitution doesn't
	// converge for n-butane at 10 K, where Psat is around 1e-116 Pa.
	_, _, at_tc, _ := psat_of(fl, 0, [m.comps[0].tc, 10.0])
	assert at_tc == [point_supercritical, point_no_convergence]
	// n-butane's α is Mathias–Copeman's, so its Psat isn't Soave's.
	_, mc, _, _ := psat_of(fl, 0, [300.0])
	assert mc[0] != thermo.saturation_pressure(thermo.new_cubic(.pr, [m.comps[0]])!, 300.0)!
	// Calls that can't be made.
	for j in [i32(-1), 2] {
		code, out, _, why3 := psat_of(fl, j, t)
		assert code == bad_argument && why3 == 'component ${j}: the fluid has 2, from 0'
		assert out.all(it == 7.0)
	}
	mut buf := []u8{len: 1024}
	assert reflux_psat(fl, 0, -1, unsafe { nil }, unsafe { nil }, unsafe { nil }, &char(buf.data),
		1024) == bad_argument
	assert c_string(buf) == 'n is -1'
	assert reflux_psat(fl, 0, 0, unsafe { nil }, unsafe { nil }, unsafe { nil }, &char(buf.data),
		1024) == ok
}
