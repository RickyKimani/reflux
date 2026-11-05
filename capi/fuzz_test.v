module capi

// Seeded fuzzing of every exported function (M10a step 6, exit criterion 2):
// handles of every kind, counts right and wrong, codes out of range, strings
// empty, junk or not UTF-8, and values from NaN and ±∞ to denormals. Nothing
// may crash, and every call must keep its contract: a known code, a known
// status for every point, NaN for a point without an answer and finite
// values for one with. Buffers are always the size the call is told: a C
// caller that lies about a buffer's size can't be caught.
import components
import math
import os
import rand

// c_string reads a NUL-terminated string out of a buffer. (Each _test.v file
// is a program of its own.)
fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

const call_codes = [ok, bad_argument, bad_handle, not_found, data_error, too_small]
const point_statuses = [point_ok, point_bad_input, point_no_convergence, point_no_volume_root,
	point_more_phases, point_phases_swapped, point_outside_envelope, point_supercritical,
	point_not_finite]

// pick is one of xs, chosen by the seeded generator.
fn pick[T](xs []T) T {
	return xs[rand.intn(xs.len) or { 0 }]
}

fn chance(p f64) bool {
	return rand.f64() < p
}

// wild is a value anything might be handed: the edges of a double, or
// something ordinary on a log scale, of either sign.
fn wild() f64 {
	return match rand.intn(12) or { 0 } {
		0 { math.nan() }
		1 { math.inf(1) }
		2 { math.inf(-1) }
		3 { 0.0 }
		4 { -0.0 }
		5 { pick([5e-324, 1e-308, 1e308, -1e308, -1.0, 1.0]) }
		6, 7 { -math.pow(10.0, rand.f64() * 24.0 - 12.0) }
		else { math.pow(10.0, rand.f64() * 24.0 - 12.0) }
	}
}

// plausible is a value near where a flash's inputs live, a T, P, H, V or β,
// so that the flashes get past their checks often enough to run.
fn plausible() f64 {
	return match rand.intn(6) or { 0 } {
		0 { 50.0 + rand.f64() * 600.0 }
		1 { math.pow(10.0, 2.0 + rand.f64() * 6.0) }
		2 { -1e5 + rand.f64() * 1e5 }
		3 { math.pow(10.0, -5.0 + rand.f64() * 3.0) }
		4 { rand.f64() }
		else { wild() }
	}
}

// feed is the amounts of n points of nc compounds: mostly not negative, so
// that points get past the checks, but now and then anything.
fn feed(n int) []f64 {
	z := inputs(n, true)
	return if chance(0.2) { z } else { z.map(math.abs(it)) }
}

// sane says whether a point's T, P (or the like) and feed are what the calls
// accept, by the rules reflux.h gives: finite and above zero, and amounts
// finite, not negative and not all zero. A point that breaks them must be
// REFLUX_POINT_BAD_INPUT.
fn sane(t f64, p f64, z []f64) bool {
	if !math.is_finite(t) || t <= 0.0 || !math.is_finite(p) || p <= 0.0 {
		return false
	}
	mut total := 0.0
	for x in z {
		if !math.is_finite(x) || x < 0.0 {
			return false
		}
		total += x
	}
	return total > 0.0 && math.is_finite(total)
}

// row_of is point i's composition, from z_rows rows.
fn row_of(z []f64, rows int, nc int, i int) []f64 {
	r := if rows == 1 { 0 } else { i }
	return z[r * nc..(r + 1) * nc]
}

fn inputs(n int, nice bool) []f64 {
	return []f64{len: n, init: if nice && chance(0.8) { plausible() } else { wild() }}
}

// count is how many entries an input has: 1, n, or a wrong count now and then.
fn count(n int) int {
	return match rand.intn(10) or { 0 } {
		0 { pick([0, 2, n + 1, -1]) }
		1, 2, 3 { 1 }
		else { n }
	}
}

const names = ['methane', 'n-butane', 'nitrogen', 'water', 'CO2', 'carbon dioxide', '74-82-8',
	'', ' ', 'unobtanium', 'methan', 'METHANE', '3967-55-3', 'VNWKTOKETHGBQD-UHFFFAOYSA-N',
	'x'.repeat(5000), 'méthane', '\t', 'n-decane', 'ethane']

// c_text is a name to hand over as C text: from the list, or bytes that aren't
// UTF-8. Every one ends in a NUL.
fn c_text() []u8 {
	if chance(0.1) {
		return [u8(0xff), 0xfe, 0x80, 0]
	}
	mut b := pick(names).bytes()
	b << 0
	return b
}

// Pool is what the fuzzing calls on: live handles of each kind, and ids that
// aren't, or are no longer, handles.
struct Pool {
	db     i64
	fluids []i64
	nc     []int
	dead   []i64
}

fn pool() !Pool {
	db := components.load_default()!
	mut fluids := []i64{}
	mut nc := []int{}
	for s in [
		db.setup(.pr, ['methane', 'n-butane'])!,
		db.setup(.srk, ['methane', 'ethane', 'n-butane'])!,
		db.setup(.pr, ['water'], alpha: .mathias_copeman)!,
		db.setup(.pr, ['water', 'n-decane', 'methane'])!,
	] {
		fluids << add(FluidItem{s})
		nc << s.records.len
	}
	freed := add(FluidItem{db.setup(.pr, ['methane'])!})
	remove(freed)!
	return Pool{add(DatabaseItem{db}), fluids, nc, [freed, 0, -1, 1 << 62, -(i64(1) << 63)]}
}

// handle is a live fluid's (with its number of compounds), the database's,
// or one that isn't a handle.
fn (p Pool) handle() (i64, int) {
	r := rand.intn(10) or { 0 }
	if r < 7 {
		k := rand.intn(p.fluids.len) or { 0 }
		return p.fluids[k], p.nc[k]
	}
	if r == 7 {
		return p.db, 3
	}
	return pick(p.dead), 3
}

fn check_code(code i32, what string) {
	assert code in call_codes, '${what}: code ${code}'
}

// points checks a many-point call's statuses and values: a point with an
// answer has finite values, and one without has NaN. It returns the
// statuses, for the tally, or none if the call failed as a whole.
fn points(code i32, status []i32, vals [][]f64, what string) []i32 {
	check_code(code, what)
	if code != ok {
		return []i32{}
	}
	for i, st in status {
		assert st in point_statuses, '${what}: status ${st} at ${i}'
		for v in vals {
			x := v[i]
			if st == point_ok {
				assert math.is_finite(x), '${what}: ${x} at ${i}, an answer'
			} else {
				assert math.is_nan(x), '${what}: ${x} at ${i}, status ${st}'
			}
		}
	}
	return status
}

fn fuzz_props(p Pool) []i32 {
	fl, nc := p.handle()
	n := pick([0, 1, 2, 5, -1])
	nn := if n < 0 { 0 } else { n }
	tl, pl, zr := count(nn), count(nn), count(nn)
	t := inputs(if tl > 0 { tl } else { 1 }, true)
	pp := inputs(if pl > 0 { pl } else { 1 }, true)
	z := feed(nc * (if zr > 0 { zr } else { 1 }))
	mut out := []f64{len: n_props * nn}
	mut lnphi := []f64{len: nn * nc}
	mut phase := []i32{len: nn}
	mut status := []i32{len: nn}
	mut buf := []u8{len: 256}
	root := i32(rand.intn(5) or { 0 }) - 1
	code := reflux_props(fl, root, n, unsafe { &f64(t.data) }, tl, unsafe { &f64(pp.data) },
		pl, unsafe { &f64(z.data) }, zr, unsafe { &f64(out.data) }, unsafe { &f64(lnphi.data) },
		unsafe { &i32(phase.data) }, unsafe { &i32(status.data) }, &char(buf.data), 256)
	mut cols := [][]f64{}
	for k in 0 .. n_props {
		cols << out[k * nn..(k + 1) * nn]
	}
	got := points(code, status, cols, 'props')
	if code == ok {
		for i, ph in phase {
			assert ph == if status[i] == point_ok { ph } else { no_phase }
			assert ph in [no_phase, liquid, vapor]
			// Bad input exactly where the rules are broken.
			ok_in := sane(t[if tl == 1 { 0 } else { i }], pp[if pl == 1 { 0 } else { i }],
				row_of(z, zr, nc, i))
			assert (status[i] == point_bad_input) == !ok_in, 'props: status ${status[i]} for ${t} ${pp} ${z}'
		}
	}
	return got
}

fn fuzz_psat(p Pool) []i32 {
	fl, nc := p.handle()
	n := pick([0, 1, 3, -1])
	nn := if n < 0 { 0 } else { n }
	t := inputs(if nn > 0 { nn } else { 1 }, true)
	mut out := []f64{len: nn}
	mut status := []i32{len: nn}
	mut buf := []u8{len: 256}
	j := i32(rand.intn(nc + 3) or { 0 }) - 1
	code := reflux_psat(fl, j, n, unsafe { &f64(t.data) }, unsafe { &f64(out.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 256)
	got := points(code, status, [out], 'psat')
	if code == ok {
		for i, st in status {
			if !sane(t[i], 1.0, [1.0]) {
				assert st == point_bad_input, 'psat: status ${st} at T = ${t[i]}'
			}
		}
	}
	return got
}

fn fuzz_stability(p Pool) []i32 {
	fl, nc := p.handle()
	n := pick([0, 1, 2, -1])
	nn := if n < 0 { 0 } else { n }
	tl, pl, zr := count(nn), count(nn), count(nn)
	t := inputs(if tl > 0 { tl } else { 1 }, true)
	pp := inputs(if pl > 0 { pl } else { 1 }, true)
	z := feed(nc * (if zr > 0 { zr } else { 1 }))
	mut stable := []i32{len: nn}
	mut tm := []f64{len: nn}
	mut w := []f64{len: nn * nc}
	mut status := []i32{len: nn}
	mut buf := []u8{len: 256}
	code := reflux_stability(fl, n, unsafe { &f64(t.data) }, tl, unsafe { &f64(pp.data) }, pl,
		unsafe { &f64(z.data) }, zr, unsafe { &i32(stable.data) }, unsafe { &f64(tm.data) },
		unsafe { &f64(w.data) }, unsafe { &i32(status.data) }, &char(buf.data), 256)
	got := points(code, status, [tm], 'stability')
	if code == ok {
		for i, s in stable {
			assert s == if status[i] == point_ok { s } else { -1 } && s in [i32(-1), 0, 1]
			ok_in := sane(t[if tl == 1 { 0 } else { i }], pp[if pl == 1 { 0 } else { i }],
				row_of(z, zr, nc, i))
			assert (status[i] == point_bad_input) == !ok_in, 'stability: status ${status[i]}'
		}
	}
	return got
}

fn fuzz_flash(p Pool) []i32 {
	fl, nc := p.handle()
	n := pick([0, 1, 2, -1])
	nn := if n < 0 { 0 } else { n }
	al, bl, zr := count(nn), count(nn), count(nn)
	a := inputs(if al > 0 { al } else { 1 }, true)
	b := inputs(if bl > 0 { bl } else { 1 }, true)
	z := feed(nc * (if zr > 0 { zr } else { 1 }))
	mut out := []f64{len: n_flash * nn}
	mut phase := []f64{len: n_phase * nn * 2}
	mut x := []f64{len: nn * 2 * nc}
	mut n_phases := []i32{len: nn}
	mut status := []i32{len: nn}
	mut buf := []u8{len: 256}
	spec := i32(rand.intn(9) or { 0 }) - 1
	code := reflux_flash(fl, spec, n, unsafe { &f64(a.data) }, al, unsafe { &f64(b.data) },
		bl, unsafe { &f64(z.data) }, zr, unsafe { &f64(out.data) }, unsafe { &f64(phase.data) },
		unsafe { &f64(x.data) }, unsafe { &i32(n_phases.data) }, unsafe { &i32(status.data) },
		&char(buf.data), 256)
	mut cols := [][]f64{}
	for k in 0 .. n_flash {
		cols << out[k * nn..(k + 1) * nn]
	}
	got := points(code, status, cols, 'flash ${spec}')
	if code == ok {
		for i, st in status {
			assert n_phases[i] == if st == point_ok { n_phases[i] } else { 0 }
			assert n_phases[i] in [i32(0), 1, 2]
			// A PT flash's T and P follow props' rules; a feed that breaks
			// them is bad input for every flash.
			ai, bi := a[if al == 1 { 0 } else { i }], b[if bl == 1 { 0 } else { i }]
			zi := row_of(z, zr, nc, i)
			if !sane(1.0, 1.0, zi) || (spec == flash_pt && !sane(ai, bi, zi)) {
				assert st == point_bad_input, 'flash ${spec}: status ${st} for ${ai} ${bi} ${zi}'
			}
		}
	}
	return got
}

fn fuzz_fluid_new(p Pool) {
	db := if chance(0.9) { p.db } else { pick(p.dead) }
	nc := pick([0, 1, 2, 3, -1])
	keep := [][]u8{len: if nc > 0 { nc } else { 1 }, init: c_text()}
	ptrs := keep.map(&char(it.data))
	no := pick([0, 0, 1, 2, -1])
	k := inputs(if no > 0 { no } else { 1 }, false)
	ka := [][]u8{len: if no > 0 { no } else { 1 }, init: c_text()}
	kb := [][]u8{len: if no > 0 { no } else { 1 }, init: c_text()}
	pa := ka.map(&char(it.data))
	pb := kb.map(&char(it.data))
	mut id := i64(0)
	mut buf := []u8{len: 256}
	code := reflux_fluid_new(db, i32(rand.intn(6) or { 0 }) - 1, unsafe { &&char(ptrs.data) },
		nc, i32(rand.intn(4) or { 0 }) - 1, i32(rand.intn(4) or { 0 }) - 1, unsafe { &&char(pa.data) },
		unsafe { &&char(pb.data) }, unsafe { &f64(k.data) }, no, &id, &char(buf.data), 256)
	check_code(code, 'fluid_new')
	if code == ok {
		assert id > 0
		remove(id) or { assert false, 'a new fluid that was no handle' }
	} else {
		assert c_string(buf) != '', 'fluid_new: code ${code} without a reason'
	}
}

// text_call is one of the functions that write text, with its buffer.
fn text_call(which int, h i64, q []u8, out &char, len int, needed &i32, buf &char) i32 {
	return match which {
		0 { reflux_compound(h, &char(q.data), out, len, needed, buf, 256) }
		1 { reflux_database_index(h, out, len, needed, buf, 256) }
		2 { reflux_fluid_setup(h, out, len, needed, buf, 256) }
		else { reflux_data_root(out, len, needed) }
	}
}

fn fuzz_text_out(p Pool) {
	h, _ := p.handle()
	which := rand.intn(4) or { 0 }
	q := c_text()
	mut buf := []u8{len: 256}
	// First the room the text takes, then a buffer at the edge of it: a byte
	// short, exact, a byte over, or anything.
	mut room := i32(0)
	text_call(which, h, q, unsafe { nil }, 0, &room, &char(buf.data))
	len := if room > 0 && chance(0.7) {
		int(room) + pick([-2, -1, 0, 1])
	} else {
		pick([0, 1, 5, 64, 1 << 16])
	}
	if len < 0 {
		return
	}
	mut out := []u8{len: len + 1, init: 0xab}
	mut needed := i32(0)
	code := text_call(which, h, q, &char(out.data), len, &needed, &char(buf.data))
	check_code(code, 'text')
	// Nothing written past the length given.
	assert out[len] == 0xab
	if code == ok {
		assert int(needed) <= len && out[needed - 1] == 0
	}
}

fn fuzz_open() {
	dirs := ['', 'nowhere', components.default_root, components.default_dir, @FILE, 'x'.repeat(3000)]
	mut root := pick(dirs).bytes()
	root << 0
	if chance(0.1) {
		root = [u8(0xff), 0]
	}
	extra := [][]u8{len: 1, init: c_text()}
	pe := extra.map(&char(it.data))
	mut id := i64(0)
	mut buf := []u8{len: 256}
	code := reflux_database_open(&char(root.data), unsafe { &&char(pe.data) }, pick([
		0,
		0,
		1,
		-1,
	]),
		&id, &char(buf.data), 256)
	check_code(code, 'open')
	if code == ok {
		remove(id) or { assert false, 'an open database that was no handle' }
	}
}

fn fuzz_small(p Pool) {
	mut buf := []u8{len: 256}
	match rand.intn(3) or { 0 } {
		0 {
			// Only a dead id may be freed here: the pool's live ones are needed.
			check_code(reflux_free(pick(p.dead), &char(buf.data), 256), 'free')
		}
		1 {
			// Around the version's own length, as well as anywhere.
			k := version().len
			len := pick([0, 1, 3, k - 1, k, k + 1, 64])
			mut v := []u8{len: len + 1, init: 0xab}
			code := reflux_version(&char(v.data), len)
			check_code(code, 'version')
			assert v[len] == 0xab
			// An answer is the whole version, never a cut one.
			if code == ok {
				assert c_string(v) == version()
			}
		}
		else {
			mut heap := i64(0)
			assert reflux_collect(&heap) == ok && heap > 0
		}
	}
}

fn test_every_function_keeps_its_contract_whatever_its_arguments() {
	// CI's seed and count; REFLUX_FUZZ_SEED and REFLUX_FUZZ_CALLS set others,
	// for a longer search by hand.
	seed := os.getenv_opt('REFLUX_FUZZ_SEED') or { '20261105' }.u32()
	n := os.getenv_opt('REFLUX_FUZZ_CALLS') or { '2000' }.int()
	rand.seed([seed, 7])
	p := pool()!
	mut calls := map[string]int{}
	// How many points of each status each function gave, to show the fuzzing
	// gets past the checks into the flashes.
	mut seen := map[string]int{}
	for _ in 0 .. n {
		match rand.intn(9) or { 0 } {
			0, 1 {
				for st in fuzz_props(p) {
					seen['props ${st}']++
				}
				calls['props']++
			}
			2 {
				for st in fuzz_psat(p) {
					seen['psat ${st}']++
				}
				calls['psat']++
			}
			3 {
				for st in fuzz_stability(p) {
					seen['stability ${st}']++
				}
				calls['stability']++
			}
			4, 5 {
				for st in fuzz_flash(p) {
					seen['flash ${st}']++
				}
				calls['flash']++
			}
			6 {
				fuzz_fluid_new(p)
				calls['fluid_new']++
			}
			7 {
				fuzz_text_out(p)
				fuzz_open()
				calls['text and open']++
			}
			else {
				fuzz_small(p)
				calls['small']++
			}
		}
	}
	println(calls)
	mut keys := seen.keys()
	keys.sort()
	println(keys.map('${it}: ${seen[it]}').join(', '))
}
