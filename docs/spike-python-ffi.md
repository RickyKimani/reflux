# Spike: calling reflux from Python (M1 step 9)

**Verdict: go**, with the conditions listed under "What M10a must do".
A V shared library works inside a Python process on Windows. NumPy arrays
cross without copying, and the results agree with the `thermo` fixtures to
5e-15. Getting there turned up five problems, and any one of them would
have broken M10a if it had been found late.

The spike lives in [spike/python_ffi/](../spike/python_ffi):

| File | Role |
|---|---|
| `ffi/ffi.v` | A V module with a few `@[export]` C-ABI functions over `thermo` |
| `run.py` | Python + NumPy scenarios, each in its own process so a crash is reported, not hidden (uv script) |
| `spike.vsh` | Builds the library and runs the scenarios |

```bash
v run spike/python_ffi/spike.vsh            # gcc
```

```bash
v run spike/python_ffi/spike.vsh -cc tcc    # shows what tcc does
```

## Results (Windows 11, V eab77e7, MSYS2 gcc 16.2, Python 3.14.9, NumPy 2.5.4)

| Scenario | gcc | gcc `-prod` | tcc |
|---|---|---|---|
| NumPy buffer modified in place (zero-copy) | ok | — | no exports |
| 7 PR points vs the `thermo` fixture | ok, 4.6e-15 | — | no exports |
| 10,000 points in one call | ok, 8.6–11.5 µs/point | 5.7 µs/point | no exports |
| the same points from a Python loop | ok, 12–15 µs/point | 9.1 µs/point | no exports |
| 500 calls × 1,000 points, heap after each 100 calls | stable at 4.3 MB after one step up | same | no exports |
| model handle held only by Python | **freed by the GC** (next call returned NaN) | — | no exports |
| model handle pinned in a V global | survives | — | no exports |
| call and collect from a non-main Python thread | ok | — | no exports |

## What went wrong, and why

1. **V's `int` is 64-bit; C's `int` is 32-bit.** The first version declared
   `n int` and `status &int`. Python passed 32-bit values (`ctypes.c_int`,
   NumPy `intc`). V wrote 8-byte statuses into a 4-byte array, so any call
   with about 1,000 points or more corrupted the heap (Windows exit code
   0xC0000374) or crashed with an access violation. Smaller calls survived
   by luck. A build with the garbage collector turned off crashed the same
   way, which ruled the collector out. `sizeof(int)` prints 8. The fix is
   `i32` at every boundary, matched by `ctypes.c_int32` and `np.int32`.
2. **tcc can't build a usable Windows DLL.** V marks exports with
   `__attribute__((visibility("default")))`. That controls exports on ELF
   (Linux) but means nothing to a Windows DLL, which needs
   `__declspec(dllexport)`. gcc falls back to exporting *every* symbol, so
   the functions are there, along with 1,600 others including the GC's
   internals. tcc exports none.
3. **The gcc DLL needs `libgcc_s_seh-1.dll`, and Python won't find it.**
   Since 3.8, Python on Windows doesn't search PATH for a DLL's
   dependencies. Linking with `-cflags -static-libgcc` removes the
   dependency; the DLL then needs only the Windows UCRT.
4. **The garbage collector frees handles that only Python holds.** A
   `&ModelBox` returned to Python as a `void*` was collected. The next call
   through it returned NaN; it could as easily have crashed. The collector
   can't see Python's memory. A handle stored in a V global array survives.
5. **A V panic inside the library takes Python down with it.** V's crash
   handler reported "print_backtrace_skipping_top_frames is not
   implemented", and the process exited. There's no way to catch a panic
   across the boundary, so the library must never panic on bad input.

What works without any help:
- **Initialisation:** V's module initialisation runs when the DLL loads, so
  there's no explicit init call.
- **The garbage collector:** it runs normally inside the Python process, and
  its heap levels off rather than growing with the number of calls.
- **Globals:** they keep their values between calls.
- **Threads:** a non-main Python thread could call in and trigger a
  collection.

## What M10a must do

- **Explicit widths at the ABI:** `i32` or `i64`, `f64`, pointers. Never
  `int`.
- **Windows builds with gcc and `-static-libgcc`**, plus a `.def` file or
  linker flags so that only the `reflux_*` functions are exported, not all
  1,600 symbols.
- **Handles are entries in a V-side registry**, with an explicit `free`.
  Better still, give Python an integer id rather than a raw pointer, so that
  a stale or double-freed handle is an error rather than a crash. (The
  registry needs `-enable-globals`, or a registry struct reached from a
  single global.)
- **No panics reachable from the API.** Every failure becomes a status code.
  `thermo` already panics only on internal bugs (`CubicKind.consts`); the
  API layer has to keep it that way.
- **Vectorised calls.** Each Python→V call costs about 3.6 µs of `ctypes`
  overhead, compared with a 5.5 µs PR solve. That's a 1.7× penalty for the
  solve, and far more for cheap functions such as evaluating a property at
  a known volume.
- **Check arrays on the Python side.** C-contiguous `float64` (and `int32`
  for statuses) via `np.ascontiguousarray`; otherwise V reads the wrong
  bytes.

## Still open

- **Concurrent calls.** `ctypes` releases the GIL during a call, so two
  Python threads can be inside V at once. The single-thread test passed, but
  Boehm GC needs every thread that allocates to be registered. Concurrency
  needs its own test, and possibly a lock or explicit thread registration.
- **Linux.** It isn't tested here. CI now runs the spike on Linux with gcc
  and with tcc, where ELF visibility may make the tcc build usable.
