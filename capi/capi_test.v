module capi

// The C API's own functions, called in process (D11), so that every CI job
// tests them, tcc's included, though only gcc can build the library: the
// header and the code agree, a reason is cut to fit at a character, a handle
// once freed is an error, and the registry holds up under threads. The Python
// tests (python/tests/) call the built library.
import components
import encoding.utf8
import os
import v.vmod

const header = os.read_file(os.join_path(@VMODROOT, 'capi', 'reflux.h')) or { panic(err) }

// Thing is something to keep in the registry.
struct Thing {
	n int
}

// c_string reads a NUL-terminated string out of a buffer.
fn c_string(buf []u8) string {
	return unsafe { cstring_to_vstring(&char(buf.data)) }
}

// defines are the header's #define REFLUX_* values.
fn defines() map[string]int {
	mut out := map[string]int{}
	for line in header.split_into_lines() {
		f := line.split(' ').filter(it != '')
		if f.len >= 3 && f[0] == '#define' && f[1].starts_with('REFLUX_') {
			out[f[1]] = f[2].trim('()').int()
		}
	}
	return out
}

// declared are the functions the header declares.
fn declared() []string {
	mut out := []string{}
	for line in header.split_into_lines() {
		t := line.trim_space()
		if t.starts_with('/*') || t.starts_with('*') || t.starts_with('#') || !t.contains('(') {
			continue
		}
		name := t.all_before('(').all_after_last(' ').trim_left('*')
		if name.starts_with('reflux_') {
			out << name
		}
	}
	return out
}

fn test_the_header_and_the_code_agree() {
	// Every #define is the code's own value, and the code has no other.
	mut want := {
		'REFLUX_ABI_VERSION':  int(abi_version)
		'REFLUX_OK':           ok
		'REFLUX_BAD_ARGUMENT': bad_argument
		'REFLUX_BAD_HANDLE':   bad_handle
		'REFLUX_NOT_FOUND':    not_found
		'REFLUX_DATA_ERROR':   data_error
		'REFLUX_TOO_SMALL':    too_small
	}
	for code, kind in equations {
		want['REFLUX_' + components.equation_name(kind).to_upper()] = code
	}
	for code, a in alphas {
		want['REFLUX_ALPHA_' + a.str().to_upper()] = code
	}
	for name, code in {
		'OK':               point_ok
		'BAD_INPUT':        point_bad_input
		'NO_CONVERGENCE':   point_no_convergence
		'NO_VOLUME_ROOT':   point_no_volume_root
		'MORE_PHASES':      point_more_phases
		'PHASES_SWAPPED':   point_phases_swapped
		'OUTSIDE_ENVELOPE': point_outside_envelope
		'SUPERCRITICAL':    point_supercritical
		'NOT_FINITE':       point_not_finite
	} {
		want['REFLUX_POINT_' + name] = code
	}
	want['REFLUX_NO_PHASE'] = no_phase
	want['REFLUX_LIQUID'] = liquid
	want['REFLUX_VAPOR'] = vapor
	for code, root in roots {
		want['REFLUX_ROOT_' + root.str().to_upper()] = code
	}
	for k, name in prop_names {
		want['REFLUX_PROP_' + name.to_upper()] = k
	}
	want['REFLUX_PROPS'] = n_props
	assert defines() == want
	assert reflux_abi_version() == abi_version
	// Every @[export] in this module is declared, and every declaration
	// exported. The C compiler checks their types on Windows; the Python
	// tests check the built library's symbol table.
	mut exports := []string{}
	for f in os.ls(os.join_path(@VMODROOT, 'capi'))! {
		if f.ends_with('.v') && !f.ends_with('_test.v') {
			for line in os.read_lines(os.join_path(@VMODROOT, 'capi', f))! {
				if line.starts_with("@[export: '") {
					exports << line.all_after("@[export: '").all_before("'")
				}
			}
		}
	}
	exports.sort()
	mut names := declared()
	names.sort()
	assert exports == names
}

fn test_the_version_is_v_mods() {
	want := vmod.decode(os.read_file(os.join_path(@VMODROOT, 'v.mod'))!)!.version
	mut buf := []u8{len: 64, init: 0xff}
	assert reflux_version(&char(buf.data), 64) == ok
	assert c_string(buf) == want
	// With no room for its NUL, or nowhere to write, nothing is written.
	mut small := []u8{len: want.len, init: 0xff}
	assert reflux_version(&char(small.data), i32(want.len)) == bad_argument
	assert small.all(it == 0xff)
	assert reflux_version(unsafe { nil }, 64) == bad_argument
	mut exact := []u8{len: want.len + 1, init: 0xff}
	assert reflux_version(&char(exact.data), i32(want.len + 1)) == ok && c_string(exact) == want
}

fn test_a_reason_is_cut_to_fit_at_a_character() {
	// ᵢ and ⱼ take three bytes each in UTF-8, so most cuts fall inside one.
	msg := 'kᵢⱼ at 300 K'
	for n in 0 .. msg.len + 3 {
		// One guard byte past the buffer the call is given.
		mut buf := []u8{len: n + 1, init: 0xff}
		write_message(&char(buf.data), i32(n), msg)
		assert buf[n] == 0xff, 'wrote past a buffer of ${n}'
		if n == 0 {
			continue
		}
		got := c_string(buf)
		assert msg.starts_with(got) && utf8.validate_str(got), '${n}: ${got}'
		// As much as fits: the next character, if any, wouldn't have.
		if got.len < msg.len {
			assert got.len + utf8_char_len(msg[got.len]) > n - 1, '${n}: ${got}'
		}
	}
	write_message(unsafe { nil }, 10, msg)
	write_message(&char(unsafe { nil }), -1, msg)
}

fn test_a_handle_once_freed_is_an_error() {
	a := add(Thing{1})
	b := add(Thing{2})
	assert b == a + 1
	item := get(a)!
	assert item is Thing && (item as Thing).n == 1
	mut buf := []u8{len: 128}
	assert reflux_free(a, &char(buf.data), 128) == ok
	for _ in 0 .. 2 {
		assert reflux_free(a, &char(buf.data), 128) == bad_handle
		assert c_string(buf) == 'handle ${a} was freed already'
	}
	if _ := get(a) {
		assert false, 'a freed handle was found'
	}
	for id in [i64(0), -1, b + 1000] {
		assert reflux_free(id, &char(buf.data), 128) == bad_handle
		assert c_string(buf) == '${id} is not a handle'
	}
	// Nowhere to write the reason: the code all the same.
	assert reflux_free(a, unsafe { nil }, 0) == bad_handle
	// The last id given out, freed, is freed already: the counter is at it.
	assert reflux_free(b, &char(buf.data), 128) == ok
	assert reflux_free(b, &char(buf.data), 128) == bad_handle
	assert c_string(buf) == 'handle ${b} was freed already'
	// Ids aren't reused.
	c := add(Thing{3})
	assert c > b
	remove(c)!
}

fn test_the_registry_holds_up_under_threads() ! {
	// It's the one state the library shares: 8 threads adding, reading and
	// removing at once must each get ids of their own, read back what they
	// put, and leave the count right. A read without the lock goes wrong only
	// when it lands in a write's resize of the map, in one round of three or
	// more (M10a step 1's planted mistakes), so there are 20 rounds, each on a
	// map that starts empty and must grow.
	registry.mu.lock()
	kept := registry.items.clone()
	registry.mu.unlock()
	for _ in 0 .. 20 {
		registry.mu.lock()
		registry.items = map[i64]Item{}
		registry.mu.unlock()
		registry_round()!
	}
	registry.mu.lock()
	registry.items = kept.clone()
	registry.mu.unlock()
}

// registry_round is one round of 8 threads at once on the registry.
fn registry_round() ! {
	before := registry.items.len
	mut ts := []thread []i64{}
	for k in 0 .. 8 {
		ts << spawn fn (k int) []i64 {
			mut ids := []i64{}
			for i in 0 .. 500 {
				id := add(Thing{k * 1000 + i})
				item := get(id) or { panic(err) }
				if item !is Thing || (item as Thing).n != k * 1000 + i {
					panic('handle ${id} has the wrong item')
				}
				ids << id
			}
			for id in ids[..250] {
				remove(id) or { panic(err) }
				if _ := get(id) {
					panic('handle ${id} is still there')
				}
			}
			return ids
		}(k)
	}
	all := ts.wait()
	mut seen := map[i64]bool{}
	for k, ids in all {
		for i, id in ids {
			assert id !in seen
			seen[id] = true
			if i >= 250 {
				item := get(id)!
				assert item is Thing && (item as Thing).n == k * 1000 + i
			}
		}
	}
	assert registry.items.len == before + 8 * 250
	for ids in all {
		for id in ids[250..] {
			remove(id)!
		}
	}
	assert registry.items.len == before
}

fn test_collect_gives_the_heap() {
	mut heap := i64(0)
	assert reflux_collect(&heap) == ok && heap > 0
	assert reflux_collect(unsafe { nil }) == ok
}
