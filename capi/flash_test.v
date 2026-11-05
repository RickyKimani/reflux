module capi

// Stability and the PT flash at many points (M10a step 4), called in process:
// every point is the V API's answer, bit for bit; each phase is in the place
// of its label, and a split's balances close; and a point that fails gets the
// status of its kind.
import components
import math
import flash
import thermo

// c_string reads a NUL-terminated string out of a buffer. (Each _test.v file
// is a program of its own.)
fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

fn fluid_of(names []string, opts components.FluidOptions) !(i64, components.Setup) {
	db := components.load_default()!
	s := db.setup(.pr, names, opts)!
	return add(FluidItem{s}), s
}

// FlashOut is reflux_flash's outputs for n points of nc compounds.
struct FlashOut {
	code     i32
	out      []f64
	phase    []f64
	x        []f64
	n_phases []i32
	status   []i32
	reason   string
}

fn call_flash(fl i64, spec i32, nc int, a []f64, b []f64, z []f64, n int) FlashOut {
	mut out := []f64{len: n_flash * n, init: 7.0}
	mut phase := []f64{len: n_phase * n * 2, init: 7.0}
	mut x := []f64{len: n * 2 * nc, init: 7.0}
	mut n_phases := []i32{len: n, init: 9}
	mut status := []i32{len: n, init: 9}
	mut buf := []u8{len: 1024, init: 0xff}
	code := reflux_flash(fl, spec, n, unsafe { &f64(a.data) }, a.len, unsafe { &f64(b.data) },
		b.len, unsafe { &f64(z.data) }, z.len / nc, unsafe { &f64(out.data) }, unsafe { &f64(phase.data) },
		unsafe { &f64(x.data) }, unsafe { &i32(n_phases.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 1024)
	return FlashOut{code, out, phase, x, n_phases, status, c_string(buf)}
}

// value is out's value k at point i of n.
fn (o FlashOut) value(k int, i int, n int) f64 {
	return o.out[k * n + i]
}

// of_phase is phase value k of the phase labelled j at point i of n.
fn (o FlashOut) of_phase(k int, i int, j int, n int) f64 {
	return o.phase[(k * n + i) * 2 + j]
}

fn (o FlashOut) comp(i int, j int, nc int) []f64 {
	return o.x[(i * 2 + j) * nc..(i * 2 + j + 1) * nc]
}

// expect checks point i against the V API's PT flash.
fn expect(f thermo.Fluid, t f64, p f64, z []f64, got FlashOut, i int, n int) ! {
	expect_result(f, flash.pt(f.model, t, p, z)!, z.len, got, i, n)
}

// expect_result checks point i against a V API flash's result.
fn expect_result(f thermo.Fluid, r flash.FlashResult, nc int, got FlashOut, i int, n int) {
	mut h, mut s, mut v := 0.0, 0.0, 0.0
	mut seen := [false, false]
	for ph in r.phases {
		pr := f.props(ph.state)
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
		v += ph.fraction * ph.state.v
		j := if ph.label == .liquid { liquid } else { vapor }
		seen[j] = true
		assert got.of_phase(0, i, j, n) == ph.fraction && got.of_phase(1, i, j, n) == ph.state.v
		assert got.of_phase(2, i, j, n) == pr.h && got.of_phase(3, i, j, n) == pr.s
		assert got.comp(i, j, nc) == ph.state.n
	}
	for j in [liquid, vapor] {
		if !seen[j] {
			for k in 0 .. n_phase {
				assert math.is_nan(got.of_phase(k, i, j, n))
			}
			assert got.comp(i, j, nc).all(math.is_nan(it))
		}
	}
	assert [got.value(0, i, n), got.value(1, i, n), got.value(2, i, n)] == [r.beta, r.t, r.p]
	assert [got.value(3, i, n), got.value(4, i, n), got.value(5, i, n)] == [h, s, v]
	assert got.n_phases[i] == r.phases.len && got.status[i] == point_ok
}

fn test_every_point_is_the_v_apis() {
	fl, s := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	// A vapour, a liquid, two splits, and amounts that don't add up to 1.
	t := [400.0, 250.0, 300.0, 300.0, 300.0]
	p := [5e5, 50e5, 20e5, 40e5, 20e5]
	z := [0.5, 0.5, 0.2, 0.8, 0.5, 0.5, 0.6, 0.4, 1.0, 1.0]
	got := call_flash(fl, flash_pt, 2, t, p, z, 5)
	assert got.code == ok && got.reason == ''
	assert got.n_phases == [i32(1), 1, 2, 2, 2]
	for i in 0 .. 5 {
		expect(s.fluid, t[i], p[i], z[i * 2..i * 2 + 2], got, i, 5)!
	}
	// The vapour is in the vapour's place, with β = 1, and the liquid in the
	// liquid's, with β = 0.
	assert got.value(0, 0, 5) == 1.0 && got.of_phase(0, 0, vapor, 5) == 1.0
	assert math.is_nan(got.of_phase(0, 0, liquid, 5))
	assert got.value(0, 1, 5) == 0.0 && got.of_phase(0, 1, liquid, 5) == 1.0
	// One value for every point is the same as one per point.
	short := call_flash(fl, flash_pt, 2, [300.0], [20e5, 40e5], [0.5, 0.5], 2)
	long := call_flash(fl, flash_pt, 2, [300.0, 300.0], [20e5, 40e5], [0.5, 0.5, 0.5, 0.5], 2)
	assert short.out == long.out && short.phase == long.phase && short.x == long.x
}

fn test_a_splits_balances_close() {
	// Checked from the arrays alone: fractions add to 1, β is the vapour's,
	// each compound's moles balance, and its fugacity is the same in both
	// phases, from states solved here.
	fl, s := fluid_of(['methane', 'ethane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	z := [0.3, 0.3, 0.4]
	t := [250.0, 280.0, 300.0]
	got := call_flash(fl, flash_pt, 3, t, [20e5], z, 3)
	for i in 0 .. 3 {
		assert got.n_phases[i] == 2
		bv, bl := got.of_phase(0, i, vapor, 3), got.of_phase(0, i, liquid, 3)
		assert math.abs(bv + bl - 1.0) < 1e-15 && got.value(0, i, 3) == bv
		y, x := got.comp(i, vapor, 3), got.comp(i, liquid, 3)
		for c in 0 .. 3 {
			assert math.abs(bv * y[c] + bl * x[c] - z[c]) < 1e-14
		}
		sv := thermo.solve_tp(s.fluid.model, t[i], 20e5, y, .stable)!
		sl := thermo.solve_tp(s.fluid.model, t[i], 20e5, x, .stable)!
		for c in 0 .. 3 {
			fv := sv.ln_phi[c] + math.log(y[c])
			fl_ := sl.ln_phi[c] + math.log(x[c])
			assert math.abs(fv - fl_) < 1e-11, '${c}: ${fv} vs ${fl_}'
		}
		// The vapour is the less dense: per mole, a larger volume and a
		// smaller mass.
		assert got.of_phase(1, i, vapor, 3) > got.of_phase(1, i, liquid, 3)
	}
}

fn test_each_failure_kind_has_its_own_status() {
	// Every flash.Failure kind gets the REFLUX_POINT_* of its own name, checked
	// for each kind here, since a state that reaches it isn't always to be had.
	// Swapped phases above an envelope come at isolated points, among
	// "outside the envelope", and rounding moves them from one platform to
	// another (M10a's exit review). capi_test.v holds the numbers to the
	// header's.
	by_name := {
		'bad_input':        point_bad_input
		'no_convergence':   point_no_convergence
		'no_volume_root':   point_no_volume_root
		'more_phases':      point_more_phases
		'phases_swapped':   point_phases_swapped
		'outside_envelope': point_outside_envelope
	}
	kinds := [flash.Failure.bad_input, .no_convergence, .no_volume_root, .more_phases, .phases_swapped,
		.outside_envelope]
	assert kinds.len == by_name.len
	mut seen := map[i32]bool{}
	for kind in kinds {
		got := flash_status(flash.FlashError{ kind: kind, reason: 'why' })
		assert got == by_name[kind.str()], '${kind}: ${got}'
		seen[got] = true
	}
	// No two kinds share one, and none is OK.
	assert seen.len == kinds.len && point_ok !in seen
}

fn test_a_flash_that_fails_says_why() {
	fl, _ := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	got := call_flash(fl, flash_pt, 2, [300.0, math.nan(), 300.0], [20e5, 20e5, 0.0], [
		0.5,
		0.5,
	], 3)
	assert got.status == [point_ok, point_bad_input, point_bad_input]
	assert got.reason == 'T = nan K: a temperature must be finite and above 0'
	for i in 1 .. 3 {
		assert got.n_phases[i] == 0
		for k in 0 .. n_flash {
			assert math.is_nan(got.value(k, i, 3))
		}
		for k in 0 .. n_phase {
			assert math.is_nan(got.of_phase(k, i, liquid, 3))
			assert math.is_nan(got.of_phase(k, i, vapor, 3))
		}
		assert got.x[i * 4..i * 4 + 4].all(math.is_nan(it))
	}
	// Water, n-decane and methane with a large water–hydrocarbon kᵢⱼ have a
	// third phase at 350 K and 5 bar: a two-phase flash refuses it.
	three, _ := fluid_of(['water', 'n-decane', 'methane'],
		kij:       false
		overrides: [
			components.Kij{
				cas1: '7732-18-5'
				cas2: '124-18-5'
				k:    0.5
			},
			components.Kij{
				cas1: '7732-18-5'
				cas2: '74-82-8'
				k:    0.5
			},
			components.Kij{
				cas1: '124-18-5'
				cas2: '74-82-8'
				k:    0.04
			},
		]
	)!
	defer {
		remove(three) or {}
	}
	got3 := call_flash(three, flash_pt, 3, [350.0], [5e5], [0.5, 0.3, 0.2], 1)
	assert got3.status == [point_more_phases]
	assert got3.reason.starts_with('flash: a phase of the split at 350.0 K and 500000.0 Pa is unstable')
	// Beyond the model: H is NaN above 1e10 K, and at 1e-300 Pa the vapour's
	// volume is past what a float holds.
	hot := call_flash(fl, flash_pt, 2, [1e10], [1e5], [0.5, 0.5], 1)
	assert hot.status == [point_not_finite] && hot.reason.starts_with('h is nan')
	thin := call_flash(fl, flash_pt, 2, [300.0], [1e-300], [0.5, 0.5], 1)
	assert thin.status == [point_no_volume_root], thin.reason
}

fn test_a_flash_call_that_cant_be_made() {
	fl, _ := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	got := call_flash(fl, 7, 2, [300.0], [1e5], [0.5, 0.5], 1)
	assert got.code == bad_argument
	assert got.reason == 'no flash 7: REFLUX_FLASH_PT (0), REFLUX_FLASH_PH (1), REFLUX_FLASH_PS (2), REFLUX_FLASH_TV (3), REFLUX_FLASH_UV (4), REFLUX_FLASH_TVF (5), REFLUX_FLASH_PVF (6)'
	assert got.out.all(it == 7.0)
	bad := call_flash(fl, flash_pt, 2, [300.0, 310.0], [1e5], [0.5, 0.5], 3)
	assert bad.code == bad_argument && bad.reason == 'a has 2 entries: 1, or n = 3'
	mut buf := []u8{len: 256}
	assert reflux_flash(fl, flash_pt, 0, unsafe { nil }, 0, unsafe { nil }, 0, unsafe { nil },
		0, unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil }, &char(buf.data),
		256) == ok
	t := [300.0]
	assert reflux_flash(fl, flash_pt, 1, unsafe { &f64(t.data) }, 1, unsafe { &f64(t.data) },
		1, unsafe { &f64(t.data) }, 1, unsafe { nil }, unsafe { nil }, unsafe { nil }, unsafe { nil },
		unsafe { nil }, &char(buf.data), 256) == bad_argument
	assert c_string(buf) == 'out is NULL'
}

// StabilityOut is reflux_stability's outputs.
struct StabilityOut {
	code   i32
	stable []i32
	tm     []f64
	w      []f64
	status []i32
	reason string
}

fn call_stability(fl i64, nc int, t []f64, p []f64, z []f64, n int) StabilityOut {
	mut stable := []i32{len: n, init: 9}
	mut tm := []f64{len: n, init: 7.0}
	mut w := []f64{len: n * nc, init: 7.0}
	mut status := []i32{len: n, init: 9}
	mut buf := []u8{len: 1024}
	code := reflux_stability(fl, n, unsafe { &f64(t.data) }, t.len, unsafe { &f64(p.data) },
		p.len, unsafe { &f64(z.data) }, z.len / nc, unsafe { &i32(stable.data) }, unsafe { &f64(tm.data) },
		unsafe { &f64(w.data) }, unsafe { &i32(status.data) }, &char(buf.data), 1024)
	return StabilityOut{code, stable, tm, w, status, c_string(buf)}
}

fn test_stability_is_the_v_apis() {
	fl, s := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	// Stable vapour, stable liquid, and a feed that splits.
	t := [400.0, 250.0, 300.0]
	p := [5e5, 50e5, 20e5]
	z := [0.5, 0.5, 0.2, 0.8, 0.5, 0.5]
	got := call_stability(fl, 2, t, p, z, 3)
	assert got.code == ok && got.status == [point_ok, point_ok, point_ok]
	assert got.stable == [i32(1), 1, 0]
	for i in 0 .. 3 {
		st := flash.stability(s.fluid.model, t[i], p[i], z[i * 2..i * 2 + 2])!
		mut best := -1
		for k, tr in st.trials {
			if !tr.trivial && (best < 0 || tr.tm < st.trials[best].tm) {
				best = k
			}
		}
		if best < 0 {
			// Every trial came back to the feed: tm is 0 and w the feed.
			assert got.tm[i] == 0.0 && got.w[i * 2..i * 2 + 2] == z[i * 2..i * 2 + 2]
		} else {
			assert got.tm[i] == st.trials[best].tm && got.w[i * 2..i * 2 + 2] == st.trials[best].w
		}
	}
	assert got.tm[2] < -1e-3 && got.tm[0] > -1e-10 && got.tm[1] > -1e-10
	// A pure compound is always stable; its every trial is the feed.
	pure, _ := fluid_of(['n-butane'])!
	defer {
		remove(pure) or {}
	}
	one := call_stability(pure, 1, [300.0], [2.5e5], [2.0], 1)
	assert one.stable == [i32(1)] && one.tm == [0.0] && one.w == [1.0]
	bad := call_stability(fl, 2, [300.0, -1.0], [1e5], [0.5, 0.5], 2)
	assert bad.status == [point_ok, point_bad_input] && bad.stable[1] == -1
	assert math.is_nan(bad.tm[1]) && bad.w[2..].all(math.is_nan(it))
	assert bad.reason == 'T = -1.0 K: a temperature must be finite and above 0'
}

// v_api is the V API's flash for a spec, as reflux_flash runs it.
fn v_api(spec i32, f thermo.Fluid, a f64, b f64, z []f64) !flash.FlashResult {
	return match spec {
		flash_pt { flash.pt(f.model, a, b, z)! }
		flash_ph { flash.ph(f, a, b, z)! }
		flash_ps { flash.ps(f, a, b, z)! }
		flash_tv { flash.tv(f.model, a, b, z)! }
		flash_uv { flash.uv(f, a, b, z)! }
		flash_tvf { flash.tvf(f.model, a, b, z)! }
		else { flash.pvf(f.model, a, b, z)! }
	}
}

fn test_every_flash_is_the_v_apis() {
	// From a split and a vapour, each flash gets its own pair of values back:
	// its answer is the V API's, bit for bit, and lands on the PT state.
	fl, s := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	f := s.fluid
	z := [0.5, 0.5]
	for tp in [[300.0, 20e5], [400.0, 5e5]] {
		t, p := tp[0], tp[1]
		r := flash.pt(f.model, t, p, z)!
		mut h, mut en, mut v, mut u := 0.0, 0.0, 0.0, 0.0
		for ph in r.phases {
			pr := f.props(ph.state)
			h += ph.fraction * pr.h
			en += ph.fraction * pr.s
			v += ph.fraction * ph.state.v
			u += ph.fraction * pr.u
		}
		mut cases := [[f64(flash_ph), p, h], [f64(flash_ps), p, en], [f64(flash_tv), t, v],
			[f64(flash_uv), u, v]]
		if r.phases.len == 2 {
			cases << [f64(flash_tvf), t, r.beta]
			cases << [f64(flash_pvf), p, r.beta]
		}
		for c in cases {
			spec := i32(c[0])
			got := call_flash(fl, spec, 2, [c[1]], [c[2]], z, 1)
			assert got.code == ok && got.status == [point_ok], '${spec}: ${got.reason}'
			expect_result(f, v_api(spec, f, c[1], c[2], z)!, 2, got, 0, 1)
			assert math.abs(got.value(1, 0, 1) - t) <= 1e-8 * t, '${spec}: T'
			assert math.abs(got.value(2, 0, 1) - p) <= 1e-8 * p, '${spec}: P'
			assert got.n_phases[0] == r.phases.len, '${spec}'
		}
	}
}

fn test_a_bubble_point_has_its_incipient_vapour() {
	// A vapour-fraction flash gives both phases, even one with no moles: at a
	// bubble point the vapour's fraction is 0 and its composition the first
	// bubble's; at a dew point the liquid's fraction is 0.
	fl, _ := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	z := [0.2, 0.8]
	got := call_flash(fl, flash_tvf, 2, [300.0], [0.0, 1.0], z, 2)
	assert got.status == [point_ok, point_ok] && got.n_phases == [i32(2), 2]
	assert got.value(0, 0, 2) == 0.0 && got.of_phase(0, 0, vapor, 2) == 0.0
	assert got.of_phase(0, 0, liquid, 2) == 1.0 && got.comp(0, liquid, 2) == z
	y := got.comp(0, vapor, 2)
	assert math.abs(y[0] + y[1] - 1.0) < 1e-15 && y[0] > z[0]
	assert got.value(0, 1, 2) == 1.0 && got.of_phase(0, 1, liquid, 2) == 0.0
	assert got.comp(1, vapor, 2) == z
	// The bubble pressure is above the dew pressure.
	assert got.value(2, 0, 2) > got.value(2, 1, 2)
}

fn test_each_flash_checks_its_values() {
	fl, _ := fluid_of(['methane', 'n-butane'])!
	defer {
		remove(fl) or {}
	}
	z := [0.5, 0.5]
	for c in [
		[f64(flash_ph), 0.0, -1e4],
		[f64(flash_ph), 1e5, math.nan()],
		[f64(flash_ps), 1e5, math.inf(1)],
		[f64(flash_tv), 300.0, 0.0],
		[f64(flash_uv), math.nan(), 1e-3],
		[f64(flash_tvf), 300.0, 1.5],
		[f64(flash_pvf), 1e5, math.nan()],
	] {
		got := call_flash(fl, i32(c[0]), 2, [c[1]], [c[2]], z, 1)
		assert got.status == [point_bad_input], '${c}'
	}
	why := fn [fl, z] (spec i32, a f64, b f64) string {
		return call_flash(fl, spec, 2, [a], [b], z, 1).reason
	}
	assert why(flash_ph, 0.0, -1e4) == 'P = 0.0 Pa: a pressure must be finite and above 0'
	assert why(flash_ph, 1e5, math.nan()) == 'H = nan J/mol: it must be finite'
	assert why(flash_ps, 1e5, math.inf(1)) == 'S = +inf J/(mol K): it must be finite'
	assert why(flash_tv, 300.0, 0.0) == 'V = 0.0 m3/mol: a molar volume must be finite and above 0'
	assert why(flash_uv, math.nan(), 1e-3) == 'U = nan J/mol: it must be finite'
	assert why(flash_tvf, 300.0, 1.5) == 'β = 1.5: a vapour fraction must be between 0 and 1'
	assert why(flash_pvf, 1e5, math.nan()) == 'β = nan: a vapour fraction must be between 0 and 1'
	// Every flash checks the feed too.
	badz := call_flash(fl, flash_ph, 2, [1e5], [-1e4], [-0.1, 1.1], 1)
	assert badz.status == [point_bad_input]
	assert badz.reason == 'z[0] = -0.1: amounts must be finite and not negative'
	// An enthalpy no temperature from 10 K to 10⁴ K reaches is the flash's
	// own bad input.
	far := call_flash(fl, flash_ph, 2, [1e5], [1e10], z, 1)
	assert far.status == [point_bad_input] && far.reason.starts_with('flash: '), far.reason
}
