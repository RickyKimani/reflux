module capi

// The database and fluid functions (M10a step 2), called in process. What
// they write is compared with what the V API gives for the same request, and
// every argument a caller could get wrong gets its code and its reason.
import components
import math
import os

const db_api = components.load_default() or { panic(err) }

// c_string reads a NUL-terminated string out of a buffer. (capi_test.v has
// one too: each _test.v file is a program of its own.)
fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

// cs is a V string as C text. V keeps a NUL after a string's bytes.
fn cs(s string) &char {
	return &char(s.str)
}

fn cs_array(xs []string) []&char {
	return xs.map(&char(it.str))
}

// TextCall is a function that writes text: (out, out_len, needed, err,
// err_len).
type TextCall = fn (&char, i32, &i32, &char, i32) i32

// text_of calls f twice, as reflux.h says to: first for the room the text
// takes, then with that room.
fn text_of(f TextCall) !string {
	mut needed := i32(-1)
	mut reason := []u8{len: 1024}
	code := f(unsafe { nil }, 0, &needed, &char(reason.data), 1024)
	if code != too_small {
		return error('code ${code}: ${c_string(reason)}')
	}
	mut out := []u8{len: int(needed)}
	again := f(&char(out.data), needed, unsafe { nil }, &char(reason.data), 1024)
	if again != ok {
		return error('code ${again}: ${c_string(reason)}')
	}
	return c_string(out)
}

fn open_db(root string) (i32, i64, string) {
	mut id := i64(0)
	mut reason := []u8{len: 1024}
	code := reflux_database_open(cs(root), unsafe { nil }, 0, &id, &char(reason.data), 1024)
	return code, id, c_string(reason)
}

fn compound_of(db i64, query string) !string {
	return text_of(fn [db, query] (o &char, l i32, n &i32, b &char, bl i32) i32 {
		return reflux_compound(db, cs(query), o, l, n, b, bl)
	})
}

fn setup_text(fluid i64) !string {
	return text_of(fn [fluid] (o &char, l i32, n &i32, b &char, bl i32) i32 {
		return reflux_fluid_setup(fluid, o, l, n, b, bl)
	})
}

// Mistake is a fluid's arguments with one thing wrong, and the reason.
struct Mistake {
	args FluidArgs
	why  string
}

@[params]
struct FluidArgs {
	equation i32      = 3
	names    []string = ['methane', 'n-butane']
	alpha    i32
	kij      i32 = 1
	pair_a   []string
	pair_b   []string
	pair_k   []f64
}

fn new_fluid(db i64, a FluidArgs) (i32, i64, string) {
	names := cs_array(a.names)
	firsts := cs_array(a.pair_a)
	seconds := cs_array(a.pair_b)
	mut id := i64(0)
	mut reason := []u8{len: 1024}
	code := reflux_fluid_new(db, a.equation, unsafe { &&char(names.data) }, i32(names.len),
		a.alpha, a.kij, unsafe { &&char(firsts.data) }, unsafe { &&char(seconds.data) },
		unsafe { &f64(a.pair_k.data) }, i32(a.pair_k.len), &id, &char(reason.data), 1024)
	return code, id, c_string(reason)
}

fn test_a_database_opens_on_a_data_root() {
	mut buf := []u8{len: 1024}
	mut needed := i32(0)
	assert reflux_data_root(&char(buf.data), 1024, &needed) == ok
	assert c_string(buf) == components.data_root() && needed == components.data_root().len + 1
	// NULL and "" are the default root.
	mut a := i64(0)
	assert reflux_database_open(unsafe { nil }, unsafe { nil }, 0, &a, unsafe { nil }, 0) == ok
	code, b, _ := open_db('')
	assert code == ok && b > a
	assert database_of(a)!.root == components.data_root()
	assert database_of(b)!.root == components.data_root()
	remove(a)!
	remove(b)!
	// A directory that isn't a data root can't be read.
	nowhere := os.join_path(os.temp_dir(), 'reflux_capi_nowhere')
	code2, _, why := open_db(nowhere)
	assert code2 == data_error
	assert why == '${nowhere} is not a data root: it has no components directory'
	// Arguments it can't use, each checked as it's made: a map keyed by the
	// calls' codes would keep only the last.
	mut id := i64(0)
	bad_utf8 := [u8(0xff), 0xfe, 0]
	assert reflux_database_open(cs(''), unsafe { nil }, 0, unsafe { nil }, &char(buf.data),
		1024) == bad_argument
	assert c_string(buf) == 'db is NULL'
	assert reflux_database_open(cs(''), unsafe { nil }, 1, &id, &char(buf.data), 1024) == bad_argument
	assert c_string(buf) == 'extra is NULL'
	assert reflux_database_open(cs(''), unsafe { nil }, -1, &id, &char(buf.data), 1024) == bad_argument
	assert c_string(buf) == 'extra has -1 entries'
	assert reflux_database_open(&char(bad_utf8.data), unsafe { nil }, 0, &id, &char(buf.data),
		1024) == bad_argument
	assert c_string(buf) == "root isn't UTF-8"
	assert id == 0
}

fn test_a_users_directory_follows_the_roots() {
	dir := os.join_path(os.temp_dir(), 'reflux_capi_mine_${os.getpid()}')
	os.rmdir_all(dir) or {}
	os.mkdir_all(dir)!
	defer {
		os.rmdir_all(dir) or {}
	}
	os.write_file(os.join_path(dir, '12345-67-4.toml'), 'format = 1\nname = "Unobtainium"\nsource = "a test"\n\n[ids]\ncas = "12345-67-4"\n\n[constants]\ntc = 500.0\n')!
	extra := cs_array([dir])
	mut id := i64(0)
	assert reflux_database_open(unsafe { nil }, unsafe { &&char(extra.data) }, 1, &id,
		unsafe { nil }, 0) == ok
	defer {
		remove(id) or {}
	}
	assert compound_of(id, 'unobtainium')!.contains('tc = 500.0')
	assert database_of(id)!.dirs == [components.default_dir, dir]
}

fn test_text_needs_room_for_its_nul() {
	code, db, _ := open_db('')
	assert code == ok
	defer {
		remove(db) or {}
	}
	want := components.write_toml(db_api.find('methane')!)!
	n := want.len + 1
	mut reason := []u8{len: 1024}
	for len in [0, 1, n - 1] {
		mut out := []u8{len: n, init: 0xff}
		mut needed := i32(0)
		assert reflux_compound(db, cs('methane'), &char(out.data), i32(len), &needed,
			&char(reason.data), 1024) == too_small
		assert needed == n && out.all(it == 0xff)
		assert c_string(reason) == 'the text takes ${n} bytes, and the buffer has ${len}'
	}
	mut out := []u8{len: n + 1, init: 0xff}
	assert reflux_compound(db, cs('methane'), &char(out.data), i32(n), unsafe { nil },
		unsafe { nil }, 0) == ok
	assert c_string(out) == want && out[n] == 0xff
}

fn test_a_compound_by_any_name() {
	code, db, _ := open_db('')
	assert code == ok
	defer {
		remove(db) or {}
	}
	want := components.write_toml(db_api.find('74-82-8')!)!
	for q in ['methane', 'METHANE', ' Methane ', '74-82-8', 'VNWKTOKETHGBQD-UHFFFAOYSA-N'] {
		assert compound_of(db, q)! == want
	}
	index := text_of(fn [db] (o &char, l i32, n &i32, b &char, bl i32) i32 {
		return reflux_database_index(db, o, l, n, b, bl)
	})!
	assert index == components.write_index(db_api.list())!
	// Not found, with the closest named; and a query it can't read.
	mut out := []u8{len: 64}
	mut reason := []u8{len: 1024}
	for q, why in {
		'methan': 'no compound `methan` in the database. Did you mean Methane (74-82-8)'
		'':       'no compound `` in the database'
	} {
		assert reflux_compound(db, cs(q), &char(out.data), 64, unsafe { nil }, &char(reason.data),
			1024) == not_found
		assert c_string(reason).starts_with(why), c_string(reason)
	}
	assert reflux_compound(db, unsafe { nil }, &char(out.data), 64, unsafe { nil }, &char(reason.data),
		1024) == bad_argument
	assert c_string(reason) == 'query is NULL'
}

fn test_a_fluid_is_the_v_apis() {
	code, db, _ := open_db('')
	assert code == ok
	defer {
		remove(db) or {}
	}
	c, fl, _ := new_fluid(db,
		names:  ['Methane', '106-97-8', 'nitrogen']
		alpha:  1
		pair_a: ['7727-37-9']
		pair_b: ['methane']
		pair_k: [0.0311]
	)
	assert c == ok
	defer {
		remove(fl) or {}
	}
	want := db_api.setup(.pr, ['74-82-8', '106-97-8', '7727-37-9'],
		alpha:     .mathias_copeman
		overrides: [
			components.Kij{
				cas1:    '7727-37-9'
				cas2:    '74-82-8'
				k:       0.0311
				comment: 'set by the caller'
			},
		]
	)!
	assert setup_text(fl)! == want.to_toml()!
	got := setup_of(fl)!
	assert got.kij.k == want.kij.k && got.alphas == want.alphas
	assert got.fluid.model.components() == want.fluid.model.components()
	// Without kᵢⱼ, every pair is missing; SRK takes srk.ipd's.
	c2, off, _ := new_fluid(db, kij: 0)
	assert c2 == ok && setup_of(off)!.kij.k == [0.0, 0.0, 0.0, 0.0]
	c3, srk, _ := new_fluid(db, equation: 2)
	assert c3 == ok && setup_of(srk)!.kij.k == db_api.setup(.srk, ['methane', 'n-butane'])!.kij.k
	remove(off)!
	remove(srk)!
}

fn test_a_fluids_mistakes_are_refused() {
	code, db, _ := open_db('')
	assert code == ok
	defer {
		remove(db) or {}
	}
	for m in [
		Mistake{FluidArgs{
			equation: 4
		}, 'no equation 4: REFLUX_VDW (0), REFLUX_RK (1), REFLUX_SRK (2) or REFLUX_PR (3)'},
		Mistake{FluidArgs{
			alpha: 2
		}, 'no alpha 2: REFLUX_ALPHA_EQUATION (0) or REFLUX_ALPHA_MATHIAS_COPEMAN (1)'},
		Mistake{FluidArgs{
			equation: 2
			alpha:    1
		}, "ChemSep's Mathias–Copeman constants (mc-PR.lib) are fitted for PR, not SRK"},
		Mistake{FluidArgs{
			kij: 2
		}, "kij is 2: 1 for ChemSep's, 0 for none"},
		Mistake{FluidArgs{
			names: []
		}, 'a fluid needs at least one compound, and nc is 0'},
		Mistake{FluidArgs{
			names: ['methane', '74-82-8']
		}, '`74-82-8` is Methane (74-82-8), which is in the list already'},
		Mistake{FluidArgs{
			pair_a: ['methane']
			pair_b: ['n-butane']
			pair_k: [math.nan()]
		}, 'pair_k[0] is nan'},
		Mistake{FluidArgs{
			pair_a: ['methane']
			pair_b: ['nitrogen']
			pair_k: [0.1]
		}, "pair_b[0] is Nitrogen (7727-37-9), which isn't in the fluid"},
		Mistake{FluidArgs{
			pair_a: ['methane']
			pair_b: ['74-82-8']
			pair_k: [0.1]
		}, 'override 0 pairs `methane` with itself'},
		Mistake{FluidArgs{
			pair_a: ['methane', 'n-butane']
			pair_b: ['n-butane', 'methane']
			pair_k: [0.1, 0.2]
		}, 'override 1 sets the pair `n-butane`, `methane` again'},
	] {
		c, id, why := new_fluid(db, m.args)
		assert c == bad_argument && why == m.why && id == 0
	}
	// A name the database lacks is not found; a compound it can't model is
	// the data's.
	c, _, why := new_fluid(db, names: ['methane', 'unobtanium'])
	assert c == not_found && why.starts_with('no compound `unobtanium` in the database')
	c2, _, why2 := new_fluid(db, names: ['3967-55-3'])
	assert c2 == data_error
	assert why2 == '1,2-Dichloroethylene Carbonate (3967-55-3) has no Tc, Pc, ω, so no equation of state can use it'
	// NULLs where something must be.
	names := cs_array(['methane'])
	mut reason := []u8{len: 1024}
	mut id := i64(0)
	k := [0.1]
	assert reflux_fluid_new(db, 3, unsafe { &&char(names.data) }, 1, 0, 1, unsafe { nil },
		unsafe { nil }, unsafe { nil }, 0, unsafe { nil }, &char(reason.data), 1024) == bad_argument
	assert c_string(reason) == 'fluid is NULL'
	assert reflux_fluid_new(db, 3, unsafe { nil }, 1, 0, 1, unsafe { nil }, unsafe { nil },
		unsafe { nil }, 0, &id, &char(reason.data), 1024) == bad_argument
	assert c_string(reason) == 'names is NULL'
	assert reflux_fluid_new(db, 3, unsafe { &&char(names.data) }, 1, 0, 1, unsafe { nil },
		unsafe { nil }, unsafe { nil }, 1, &id, &char(reason.data), 1024) == bad_argument
	assert c_string(reason) == 'pair_k is NULL'
	assert reflux_fluid_new(db, 3, unsafe { &&char(names.data) }, 1, 0, 1, unsafe { nil },
		unsafe { &&char(names.data) }, unsafe { &f64(k.data) }, 1, &id, &char(reason.data), 1024) == bad_argument
	assert c_string(reason) == 'pair_a is NULL'
	assert id == 0
}

fn test_a_handle_is_of_its_kind() {
	code, db, _ := open_db('')
	assert code == ok
	c, fl, _ := new_fluid(db)
	assert c == ok
	mut out := []u8{len: 64}
	mut reason := []u8{len: 1024}
	assert reflux_compound(fl, cs('methane'), &char(out.data), 64, unsafe { nil }, &char(reason.data),
		1024) == bad_handle
	assert c_string(reason) == 'handle ${fl} is a fluid, not a database'
	assert reflux_fluid_setup(db, &char(out.data), 64, unsafe { nil }, &char(reason.data), 1024) == bad_handle
	assert c_string(reason) == 'handle ${db} is a database, not a fluid'
	// A fluid outlives its database, and a freed one is gone.
	assert reflux_free(db, unsafe { nil }, 0) == ok
	assert setup_text(fl)!.contains('name = "Methane"')
	c2, _, why := new_fluid(db)
	assert c2 == bad_handle && why == 'handle ${db} was freed already'
	assert reflux_free(fl, unsafe { nil }, 0) == ok
	assert reflux_fluid_setup(fl, &char(out.data), 64, unsafe { nil }, &char(reason.data), 1024) == bad_handle
	assert c_string(reason) == 'handle ${fl} was freed already'
}
