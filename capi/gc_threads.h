// Thread registration for V code called from threads V didn't create.
//
// Boehm GC must know every thread that runs V code, to scan its stack and
// pause it during a collection. A thread created by Python is unknown to it,
// and a collection started from one aborts with "Collecting from unknown
// thread": on Linux (found by M1's spike in CI) and on Windows (found by
// M10a's probes; V 0.5.2 builds the collector with GC_THREADS on both). V's
// start-up calls GC_allow_register_threads(), so a thread may then register
// itself. Without GC_THREADS these are no-ops.
#ifndef REFLUX_GC_THREADS_H
#define REFLUX_GC_THREADS_H

#include <gc.h>

// The library's first call starts the collector in the calling thread (V
// starts a shared library then), and the collector registers that thread as
// its main one, for good. If the thread then ends, as a Python thread that
// imports the package may, it stays registered:
// - on Windows that's harmless: before suspending a thread for a collection,
//   the collector asks whether it has exited, and drops it if so;
// - on Linux it signals the thread and waits for an answer that never comes,
//   and aborts: "Signals delivery fails constantly" (CI, M10a step 3).
// So on Linux, in the library (CAPI_LIBRARY, defined by tools/python.vsh),
// the first call's thread is unregistered when the call ends, as Boehm
// allows once for the main thread, unless it's the process's first thread,
// which lives as long as the process: Python's main thread, usually. It then
// registers on each call, like any other thread. A V program calling these
// functions in process, as capi's tests do, owns its main thread, so there
// it's left alone. Elsewhere than Linux and Windows, the thread is always
// unregistered: correct, if slower.
#if defined(CAPI_LIBRARY) && !defined(_WIN32)
#if defined(__linux__)
#include <sys/syscall.h>
#include <unistd.h>

// capi_may_end says whether the calling thread could end before the
// process: any but the process's first, whose thread id is the process id.
static int capi_may_end(void) {
	return syscall(SYS_gettid) != getpid();
}
#else
static int capi_may_end(void) {
	return 1;
}
#endif
#endif

// capi_gc_enter registers the calling thread if the collector doesn't know
// it, and returns 1 if it did, so that capi_gc_leave unregisters it. On the
// library's first call it may return 1 for the thread the collector started
// in: see above.
static int capi_gc_enter(void) {
#if defined(GC_THREADS)
	struct GC_stack_base sb;
	if (GC_thread_is_registered()) {
#if defined(CAPI_LIBRARY) && !defined(_WIN32)
		// The first call doesn't race (reflux.h), but the flag is atomic all
		// the same: it's read on every call from a registered thread.
		static int started = 0;
		if (!__atomic_exchange_n(&started, 1, __ATOMIC_SEQ_CST) && capi_may_end()) {
			return 1;
		}
#endif
		return 0;
	}
	if (GC_get_stack_base(&sb) != GC_SUCCESS) {
		return 0;
	}
	GC_register_my_thread(&sb);
	return 1;
#else
	return 0;
#endif
}

// capi_gc_leave undoes capi_gc_enter, so a thread that later exits leaves no
// stale registration behind.
static void capi_gc_leave(int registered) {
#if defined(GC_THREADS)
	if (registered) {
		GC_unregister_my_thread();
	}
#else
	(void)registered;
#endif
}

#endif
