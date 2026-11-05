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
	r := flash.pt(f.model, t, p, z)!
	nc := z.len
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
	got := call_flash(fl, 1, 2, [300.0], [1e5], [0.5, 0.5], 1)
	assert got.code == bad_argument && got.reason == 'no flash 1: REFLUX_FLASH_PT (0)'
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
