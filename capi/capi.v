module capi

// Package capi is reflux's C API (M10a): the functions capi/reflux.h declares,
// built into a shared library by `v run tools/python.vsh build` and called
// from Python through ctypes, or from any language that loads a C library.
//
// Every exported function follows the header's conventions (D2): integers are
// i32 or i64, never `int`, which is 64-bit in V and 32-bit in C; arrays are a
// pointer and a length; the result is ok or a negative code, with the reason
// written into the caller's buffer. Each one registers the calling thread with
// the collector for the length of the call (D5), and none may panic: a panic
// would take the caller's process with it.
import v.vmod

// The header is included where V's integer types are C's: V's i64 is `long
// long`, which is int64_t on Windows and macOS but not on Linux, where
// int64_t is `long`. The C compiler then checks every exported function
// against its declaration, and refuses a header that disagrees.
#flag -I @VMODROOT/capi

$if windows {
	#include "reflux.h"
}

// gc_threads.h registers a thread the collector doesn't know (Python's).
// @VEXEROOT is V's install directory, where the collector's headers are.
#flag -I @VEXEROOT/thirdparty/libgc/include
#include "gc_threads.h"

fn C.capi_gc_enter() int
fn C.capi_gc_leave(registered int)

// abi_version is reflux.h's REFLUX_ABI_VERSION.
const abi_version = i32(1)

// A call's result (reflux.h).
const ok = i32(0)
const bad_argument = i32(-1)
const bad_handle = i32(-2)

// enter registers the calling thread with the collector if it doesn't know
// it, and says whether it did. Every exported function calls it first, and
// leave on its way out.
fn enter() int {
	return C.capi_gc_enter()
}

fn leave(registered int) {
	C.capi_gc_leave(registered)
}

// write_message copies msg into the caller's buffer, NUL-terminated and cut
// to fit at the start of a character, so that a reason in UTF-8 stays valid.
// A nil buffer, or one of no length, gets nothing.
fn write_message(buf &char, len i32, msg string) {
	if buf == unsafe { nil } || len <= 0 {
		return
	}
	mut n := msg.len
	if n > int(len) - 1 {
		n = int(len) - 1
		// A UTF-8 continuation byte is 10xxxxxx: back off to a first byte.
		for n > 0 && msg[n] & 0xC0 == 0x80 {
			n--
		}
	}
	unsafe {
		vmemcpy(buf, msg.str, n)
		buf[n] = char(0)
	}
}

// fail writes a reason and returns its code, for a call's error paths.
fn fail(code i32, buf &char, len i32, msg string) i32 {
	write_message(buf, len, msg)
	return code
}

// version is v.mod's, as `reflux version` prints it.
fn version() string {
	m := vmod.decode(@VMOD_FILE) or { return 'unknown' }
	return m.version
}

// reflux_abi_version returns REFLUX_ABI_VERSION (reflux.h).
@[export: 'reflux_abi_version']
pub fn reflux_abi_version() i32 {
	r := enter()
	defer {
		leave(r)
	}
	return abi_version
}

// reflux_version writes v.mod's version into buf.
@[export: 'reflux_version']
pub fn reflux_version(buf &char, len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	text := version()
	if buf == unsafe { nil } || int(len) <= text.len {
		return bad_argument
	}
	write_message(buf, len, text)
	return ok
}

// reflux_free drops a handle.
@[export: 'reflux_free']
pub fn reflux_free(id i64, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	remove(id) or { return fail(bad_handle, buf, buf_len, err.msg()) }
	return ok
}

// reflux_collect runs a full collection, and gives the heap's size.
@[export: 'reflux_collect']
pub fn reflux_collect(heap_bytes &i64) i32 {
	r := enter()
	defer {
		leave(r)
	}
	gc_collect()
	if heap_bytes != unsafe { nil } {
		unsafe {
			*heap_bytes = i64(gc_heap_usage().heap_size)
		}
	}
	return ok
}
