// Module ffi is the M1 step 9 spike: a few C-ABI functions over reflux's
// thermo, built as a shared library and called from Python. It is not the
// M10a API, only enough to find out what breaks. Build and run it with
// `v run spike/python_ffi/spike.vsh`.
//
// Every integer at the boundary is i32, never `int`. V's `int` is 64-bit
// on 64-bit platforms (like Go's), while C's `int`, ctypes.c_int and
// NumPy's intc are 32-bit. The first version of this spike used `int`,
// wrote 8-byte values into a 4-byte status array, and corrupted the heap.
module ffi

import math
import thermo

// A handle Python keeps between calls. The question the spike asks: does
// V's garbage collector know that Python still holds it?
struct ModelBox {
	eos thermo.Cubic
}

// pinned_models keeps handles reachable from V's own memory, so the
// collector can't free them while Python holds them. Globals need
// -enable-globals.
__global pinned_models = []voidptr{}

fn pure_pr(tc f64, pc f64, omega f64) !thermo.Cubic {
	return thermo.new_cubic(.pr, [
		thermo.Component{
			name:  'ffi'
			tc:    tc
			pc:    pc
			omega: omega
			mw:    0.1
		},
	])
}

fn root_of(r i32) thermo.Root {
	return match r {
		0 { thermo.Root.liquid }
		1 { thermo.Root.vapor }
		else { thermo.Root.stable }
	}
}

// reflux_scale multiplies n doubles in place: the zero-copy check, where
// Python passes a NumPy array's own buffer.
@[export: 'reflux_scale']
pub fn reflux_scale(x &f64, n i32, k f64) {
	for i in 0 .. n {
		unsafe {
			x[i] *= k
		}
	}
}

// reflux_pr_z_many solves PR for n (T, P) points in one call, writing Z and
// a status per point (0 ok, 1 no root). It returns the number of failures,
// or −1 if the model can't be built.
@[export: 'reflux_pr_z_many']
pub fn reflux_pr_z_many(n i32, t &f64, p &f64, root i32, tc f64, pc f64, omega f64, z_out &f64, status &i32) i32 {
	eos := pure_pr(tc, pc, omega) or { return -1 }
	mut failures := i32(0)
	for i in 0 .. n {
		ti, pi := unsafe { t[i] }, unsafe { p[i] }
		st := thermo.solve_tp(eos, ti, pi, [1.0], root_of(root)) or {
			unsafe {
				z_out[i] = math.nan()
				status[i] = 1
			}
			failures++
			continue
		}
		unsafe {
			z_out[i] = st.z
			status[i] = 0
		}
	}
	return failures
}

// reflux_model_new returns a handle that only Python refers to.
@[export: 'reflux_model_new']
pub fn reflux_model_new(tc f64, pc f64, omega f64) voidptr {
	eos := pure_pr(tc, pc, omega) or { return unsafe { nil } }
	return voidptr(&ModelBox{
		eos: eos
	})
}

// reflux_model_new_pinned returns a handle that V also keeps a reference
// to, until reflux_model_free.
@[export: 'reflux_model_new_pinned']
pub fn reflux_model_new_pinned(tc f64, pc f64, omega f64) voidptr {
	h := reflux_model_new(tc, pc, omega)
	if h != unsafe { nil } {
		pinned_models << h
	}
	return h
}

// reflux_model_free drops V's reference to a pinned handle.
@[export: 'reflux_model_free']
pub fn reflux_model_free(h voidptr) {
	i := pinned_models.index(h)
	if i >= 0 {
		pinned_models.delete(i)
	}
}

// reflux_model_z solves one point with a model handle; NaN on failure.
@[export: 'reflux_model_z']
pub fn reflux_model_z(h voidptr, t f64, p f64, root i32) f64 {
	m := unsafe { &ModelBox(h) }
	st := thermo.solve_tp(m.eos, t, p, [1.0], root_of(root)) or { return math.nan() }
	return st.z
}

// reflux_gc_collect runs a full collection and returns the heap size, bytes.
@[export: 'reflux_gc_collect']
pub fn reflux_gc_collect() u64 {
	gc_collect()
	return u64(gc_heap_usage().heap_size)
}

// reflux_churn allocates and drops n small arrays, so that memory freed by
// a collection gets reused and overwritten.
@[export: 'reflux_churn']
pub fn reflux_churn(n i32) f64 {
	mut s := 0.0
	for i in 0 .. n {
		a := []f64{len: 16, init: f64(i + index)}
		s += a[15]
	}
	return s
}
