// Thread registration for V code called from threads V didn't create.
//
// Boehm GC must know every thread that runs V code, to scan its stack and
// pause it during a collection. A thread created by Python is unknown to
// it, and on Linux a collection started from one aborts with "Collecting
// from unknown thread". This is only needed, and only possible, when the
// collector is built with thread support (GC_THREADS), as it is on Linux.
// V's own startup code already calls GC_allow_register_threads() in that
// case. Without GC_THREADS, as in V's Windows build, these are no-ops.
#ifndef REFLUX_GC_THREADS_H
#define REFLUX_GC_THREADS_H

#include <gc.h>

// reflux_enter registers the calling thread if the collector doesn't know
// it, and returns 1 if it did.
static int reflux_enter(void) {
#if defined(GC_THREADS)
	struct GC_stack_base sb;
	if (GC_thread_is_registered()) {
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

// reflux_leave undoes reflux_enter, so a thread that later exits leaves no
// stale registration behind.
static void reflux_leave(int registered) {
#if defined(GC_THREADS)
	if (registered) {
		GC_unregister_my_thread();
	}
#else
	(void)registered;
#endif
}

#endif
