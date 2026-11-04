/*
 * reflux.h: reflux's C API (M10a), ABI version 1.
 *
 * The library is built from the V module `capi` by `v run tools/python.vsh
 * build`, with gcc, and exports exactly the functions declared here.
 *
 * Conventions (docs/m10a-plan.md, D2):
 * - Every function returns an int32_t: REFLUX_OK, or a negative code for a
 *   failure of the whole call. Where a function takes a buffer (err, err_len),
 *   it writes the reason there, NUL-terminated and cut to fit at a character
 *   boundary. A NULL buffer, or one of no length, gets nothing.
 * - Integers are int32_t or int64_t, never int. An array is a pointer and a
 *   length, and the caller allocates every array a function fills.
 * - A handle is an int64_t id that is never used twice, so a stale or freed
 *   one is REFLUX_BAD_HANDLE rather than a crash.
 * - Any thread may call, and several at once. But the library starts itself
 *   on its first call, and that start must not race: make one call, such as
 *   reflux_abi_version(), from one thread before calling from several.
 */
#ifndef REFLUX_H
#define REFLUX_H

#include <stdint.h>

#define REFLUX_ABI_VERSION 1

/* A call's result. */
#define REFLUX_OK 0
#define REFLUX_BAD_ARGUMENT (-1) /* an argument it can't use */
#define REFLUX_BAD_HANDLE (-2)   /* an id that isn't a live handle */

/* reflux_abi_version returns REFLUX_ABI_VERSION, for a caller to check before
   it calls anything else. */
int32_t reflux_abi_version(void);

/* reflux_version writes reflux's version, such as "0.1.0", into buf. A buffer
   too small for it is REFLUX_BAD_ARGUMENT. */
int32_t reflux_version(char *buf, int32_t len);

/* reflux_free drops a handle. */
int32_t reflux_free(int64_t id, char *err, int32_t err_len);

/* reflux_collect runs a full garbage collection, and writes the collector's
   heap size in bytes to heap_bytes unless it's NULL: for tests of memory. */
int32_t reflux_collect(int64_t *heap_bytes);

#endif
