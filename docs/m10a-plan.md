# M10a plan — thermo and the flash from Python

M10a makes reflux's thermodynamics and flash callable from Python, with NumPy
arrays in and out, through a C API that any language can load. It answers one
of the README's DWSIM complaints directly: you can't bring NumPy, pandas or
Jupyter to your own analysis. It's planned like M1 to M4: research first, then
exit criteria, decisions and steps. The roadmap entry is in
[roadmap.md](roadmap.md), and M1's spike, which retired the biggest risks, in
[spike-python-ffi.md](spike-python-ffi.md).

Status: **approved** on 2026-11-04. Steps 1 to 9 are done, and step 10 is
next.

## Goal

- **A C API** in a shared library built with `v -shared`: plain functions
  with explicit integer widths, `f64` arrays passed as pointer and length, so
  that NumPy arrays cross without copying (roadmap I1).
- **Vectorised calls:** one call evaluates or flashes a whole array of
  states, with a status for each, so a failure at one point doesn't throw away
  the rest.
- **What it covers:** the component database and fluids built from it; a
  phase's properties at (T, P, z); pure-component vapour pressure; and every
  flash M3 has: stability, PT, PH, PS, TV, UV, and the vapour-fraction flashes
  with bubble and dew points.
- **A `reflux` Python package** over the C API with `ctypes`, so users need no
  compiler. Results come back as NumPy arrays.
- **Safe from Python threads,** which `ctypes` lets into the library at the
  same time.

Flowsheets from Python are M10b. Wheels that bundle the library come later
(I2): M10a's package installs from the source tree.

## What the research found

**What M1's spike settled** ([spike-python-ffi.md](spike-python-ffi.md)), and
CI has run since on Windows and Linux:
- V's runtime and its garbage collector work inside a Python process. Module
  initialisation runs when the library loads, so there's no init call. (Not
  quite: step 1 found that it runs on the first call, which mustn't race; see
  D5.)
- **Integers at the boundary must be `i32` or `i64`.** V's `int` is 64-bit
  and C's 32-bit; the first version wrote 8-byte statuses into a 4-byte array
  and corrupted the heap.
- **Shared libraries need gcc** (or clang): tcc's Windows DLLs export nothing,
  and its Linux `.so` asks for an executable stack that glibc refuses. On
  Windows, `-static-libgcc` removes a DLL Python wouldn't find, and gcc
  exports all 1,600 symbols unless told otherwise.
- **A handle only Python holds is freed by the collector,** which can't see
  Python's memory. Handles must live in a V-side registry.
- **A V panic takes Python down with it,** so no input may reach one.
- **A Python call costs about 3.6 µs** in `ctypes` overhead, against 5.5 µs
  for a PR volume root. Calls must be vectorised.

**What this plan's probes found** (V 0.5.2, gcc 16.2, Windows 11, Python 3.13,
NumPy 2):
- **The collector has thread support on Windows too.** V 0.5.2 defines
  `GC_THREADS` and `GC_WIN32_THREADS` for gcc builds on every OS
  (`vlib/builtin/builtin_d_gcboehm.c.v`). A probe DLL reports the thread that
  loaded it as registered, and a thread Python starts as not. So the spike's
  note that "V's Windows build has none" is wrong. Registering the calling
  thread at every entry point, as the spike's `gc_threads.h` does, is needed on
  Windows as much as on Linux. The spike's helper already does it, since it
  tests `GC_THREADS`.
- **Concurrent calls are correct.** 8 Python threads, each solving PR for
  10,000 points 30 times through the spike's library while a ninth thread
  churned allocations and forced collections: every result equal to the serial
  one bit for bit, and no crash.
- **They don't scale far.** Without the forced collections, on 8 logical
  CPUs: 1.6× the serial rate with 2 threads, 2.0× with 4 and 2.1× with 8. Each
  point allocates, and the collector's lock and its collections are shared.
  Parallel speed is M8's question (its sweeps use `spawn`), not M10a's.
- **reflux has no global state:** no `__global` and no mutable module-level
  value in `thermo`, `flash`, `numeric`, `tdep` or `components`. The only state
  the library will share between threads is its registry of handles.

**Where the data comes from.** `components.load_default()` reads the compound
files from `$REFLUX_DATA`, or else from the source tree's `data/components/`,
whose path is baked in when V compiles (`@VMODROOT`). But the kᵢⱼ banks
(`pr.ipd`, `srk.ipd`) and Mathias–Copeman's `mc-PR.lib` always come from the
source tree's `data/vendor/chemsep/`, whatever `REFLUX_DATA` says. A library
moved away from its source tree would lose its kᵢⱼ. The CLI has the same flaw,
which v0.1.0's "Known limitations" describes only in part.

**Bit for bit means one build.** M4 found tcc's and gcc's builds differing in
the last bits (a drum's duty of 0 under tcc, −9.3e-13 kW under gcc), and
Linux's libm and Windows' UCRT needn't agree in the last bit of `exp`, `log`
and `pow`. So "the same numbers as the V API" can only be checked against a V
program built by the same compiler, with the same flags, on the same machine.

**What other libraries do.** CoolProp's `PropsSI` accepts arrays and loops in
C++, one property at a time. Caleb Bell's `thermo` flashes one state per call,
as a Python object. reflux's flashes return a whole state, so a vectorised
call fills several arrays at once: β, the phases' compositions and
properties, and a status.

## Exit criteria

All checked by tests in CI: V tests of the C API under both compilers, and
Python tests against the library in the gcc jobs, on Windows and Linux.

1. **The C API is what its header says.** `capi/reflux.h` declares every
   function, and the library exports exactly those, checked against its symbol
   table on both OSes. Every integer is `i32` or `i64`, every array a pointer
   and a length, and `reflux_abi_version()` gives the ABI's version.
2. **Nothing panics.** Every function, called with bad input (NaN, ±∞, zero or
   negative sizes, mismatched lengths, a mole fraction below 0, an unknown,
   stale or freed handle, an empty name) returns an error, never a crash. Each
   is a test, and seeded fuzzing calls every function with mutated arguments.
3. **The same numbers as the V API, bit for bit.** For a reference set of
   states covering every function, Python's results equal a V program's that
   calls the V API directly, built with the same compiler and flags on the same
   machine.
4. **Vectorised, without copies.** A call takes arrays of N states and fills
   arrays the caller allocated. NumPy inputs that are already contiguous
   `float64` aren't copied. Python broadcasts T, P and z by NumPy's rules, and
   z can be one composition for every point.
5. **Each point has its own status.** A failure at one point leaves the others'
   answers as they would be alone. Each status names the flash's `Failure` kind
   (`no_convergence`, `more_phases`, …) or bad input, and the reason for a
   point can be had in words.
6. **Handles are safe and memory is flat.** A handle is an id that is never
   reused, so a stale or freed one is an error. Python frees its handles when
   its objects go, and the collector's heap levels off over 100,000 calls.
7. **Threads are safe.** 8 Python threads calling every vectorised function at
   once get the serial answers bit for bit, on Windows and Linux. The
   throughput is measured and recorded, with no target.
8. **The exit script.** A script computes a T–x–y diagram and a phase-fraction
   map into NumPy arrays and plots them with matplotlib. It runs in CI, and its
   numbers match the V API's (criterion 3).
9. **It installs, and finds its data.** On Windows and Linux, one command
   builds the library and another installs the package from the source tree.
   The package finds the component data, the kᵢⱼ banks and `mc-PR.lib` with it,
   from wherever it's installed.
10. **Speed is measured.** Per-call and per-point overhead against V's own
    loop, for properties and a PT flash at 1, 100 and 10,000 points, with gcc
    `-prod`, in `docs/benchmarks.md`. There's no target, as in M3 and M4.
11. CI is green, and the docs are updated: README, roadmap, design notes, V
    notes, the C API's header and a new `docs/python.md`, whose every example
    runs in CI.

## Decisions

| # | Decision | Choice | Why |
|---|---|---|---|
| D1 | Where it lives | A V module `capi` with the exported functions, the registry and the thread registration, and `capi/reflux.h`. A Python package in `python/` (`python/reflux/`, its `pyproject.toml` and `tests/`). `tools/python.vsh` builds the library into the package and runs the Python tests. The spike's scenarios move into those tests, and its code and CI step go; its write-up stays as the record, corrected on Windows threads. | One module owns the boundary, so nothing else in reflux needs `unsafe` or C types. The package beside it keeps the two in step. |
| D2 | ABI conventions | Every function is `reflux_*`, takes `i32`, `i64`, `f64` and pointers only, and returns an `i32`: 0, or a negative code for a failure of the whole call (a bad handle, wrong sizes, a lookup that fails), whose message it writes into a buffer the caller passes (pointer and length). Outputs go into arrays the caller allocated. No V memory crosses the boundary, only handles. | The caller owns every byte it sees, so nothing needs freeing on the Python side, and a message needs no shared state between threads. |
| D3 | Exports | Only the header's functions: on Linux by `-fvisibility=hidden`, since V marks its exports visible; on Windows by a `.def` file generated from the header. A test compares the symbol table with the header on both OSes. **Amended in step 1:** on Linux by a linker version script, not `-fvisibility=hidden`, which can't hide the collector's symbols: they come from an object V compiled before (V's own macOS code says as much). And on Windows the header is compiled with the module, so the C compiler checks every declaration against the code. | Not the collector's internals, which gcc exports today on Windows, and one list of truth. |
| D4 | Handles | An `i64` id into a registry: a map from id to fluid (or database) behind a `sync.Mutex`, ids from a counter that never repeats. `reflux_free(id)` drops one. In Python each handle belongs to an object with `close()`, a context manager and a `weakref.finalize`. | The collector can't see Python's references (the spike). A fresh id per object turns a stale or doubly freed handle into an error, not a crash. |
| D5 | Threads | Every exported function registers the calling thread with the collector for the length of the call, on every OS (the spike's `gc_threads.h`). Concurrent calls are allowed; the registry's lock is the only shared state. **Amended in step 1:** V starts a shared library on its first call, not when it loads, behind a flag that isn't thread-safe, so the first call mustn't race. The package makes it at import, which Python's import lock serialises, and `reflux.h` tells other callers to make one call from one thread first. **Amended in step 3:** that first call's thread is where the collector starts, and it keeps it registered. On Linux a registered thread that has ended aborts the next collection, so the library unregisters it after its first call, unless it's the process's first thread. | The probes: needed on Windows too, and correct under load. Their speed isn't M10a's target. |
| D6 | Vectorised results | Arrays of N points, C-contiguous `float64` (and `int32` for statuses), shaped for two phases: β, T, P, the number of phases, each phase's label, fraction, composition (N × 2 × c), molar V, H and S, and the whole's H, S and V. An absent phase is NaN. M6's third phase will be a new ABI version. **Amended in step 4:** the two phases are in the places of their labels, the liquid's (0) then the vapour's (1), as `REFLUX_LIQUID` and `REFLUX_VAPOR` number them. So a phase's label is its place, with no array of labels, and the vapour's composition is always the same column. A single phase takes its label's place, the other's values NaN. | One layout for every flash, which is what a plot or a DataFrame wants. Two phases is all M3's flash returns. |
| D7 | Statuses | Per point, an `int32`: 0 for an answer, a stable number for each `FlashError` kind, and one for bad input. Python maps them to an `IntEnum`. A point's reason in words comes from the single-point form of the same call, which writes the message into the caller's buffer. | Numbers fill an array cheaply; words for every point would cost a string each. |
| D8 | The Python API | `reflux.Database()` (or a data directory), `db.fluid('PR', ['methane', 'n-butane'], alpha=..., kij=...)`, then on the fluid `props`, `psat`, `stability`, `flash_pt`, `flash_ph`, `flash_ps`, `flash_tv`, `flash_uv`, `flash_tvf`, `flash_pvf`, `bubble_t`, `bubble_p`, `dew_t` and `dew_p`, each taking arrays or scalars. A flash returns a `FlashResult` of arrays with a `status` array and `ok` mask. SI throughout, as the V API. Python 3.11 or later, NumPy 2. | The V API's names and units, so that the V docs and design notes serve Python too. |
| D9 | Data | A database opens on a data directory: the caller's, or `$REFLUX_DATA`, or the source tree's. The kᵢⱼ banks and `mc-PR.lib` come from the same root (`components` changes so they do), which fixes the CLI too. The package records the source tree it was built from. | A library away from its source tree must still find its kᵢⱼ. |
| D10 | Build and install | `v run tools/python.vsh build` makes `python/reflux/_lib/reflux.dll` (or `.so`) with gcc `-prod`, then `pip install -e python` (or `uv pip install -e python`). The package's version is v.mod's. | No compiler for users of a built tree, and wheels can come later without changing the package (I2). |
| D11 | Tests | V tests of `capi` call the exported functions in process, under tcc and gcc. Python tests with pytest, pinned, run by `tools/python.vsh test` in the gcc jobs on both OSes. A scenario that could crash runs in a process of its own, as the spike's did, so a crash is reported. | Most of the boundary's logic is tested in every job; what only Python can show is tested where the library can be built. |
| D12 | The V reference | A V program built by `tools/python.vsh` with the library's compiler and flags writes the V API's answers for the reference states to a temporary JSON file, and the Python tests compare bit for bit. It's never committed. | The research: bits differ between compilers and between OSes' maths libraries. |
| D13 | Versions | The ABI has its own number (`reflux_abi_version()`), as the file and report formats do. M10a ends with v0.2.0: a minor release before 1.0 for an addition, by the CHANGELOG's rules. | The CHANGELOG's policy, and a version Python can check before it calls anything. |
| D14 | Who writes what | As before, I write it all unless you claim a step. | |

## Steps

```mermaid
flowchart LR
  S1[1 Library skeleton] --> S2[2 Database, fluids, data]
  S2 --> S3[3 Properties, Psat]
  S3 --> S4[4 Stability, PT flash]
  S4 --> S5[5 The other flashes]
  S5 --> S6[6 Threads, memory, fuzzing]
  S5 --> S7[7 V reference, speed]
  S6 --> S8[8 Exit script, docs/python.md]
  S7 --> S8
  S8 --> S9[9 Release v0.2.0]
  S9 --> S10[10 Docs and exit review]
```

Each step from 2 on adds C functions with their V tests, their Python wrappers
and their Python tests together.

### 1. The library skeleton — done
- **The C API,** in [capi/](../capi/):
  - [reflux.h](../capi/reflux.h) declares four functions:
    - `reflux_abi_version`, which returns 1;
    - `reflux_version`;
    - `reflux_free`;
    - `reflux_collect`: a full collection, and the heap's size, for step 6's
      memory test.
  - A call returns `REFLUX_OK` (0), `REFLUX_BAD_ARGUMENT` (−1) or
    `REFLUX_BAD_HANDLE` (−2), with the reason in the caller's buffer.
  - [capi.v](../capi/capi.v) has the exported functions. Each registers its
    thread with the collector for the call (D5), through the spike's helper,
    now [gc_threads.h](../capi/gc_threads.h).
  - A reason is cut to fit the caller's buffer at the start of a character, so
    it stays valid UTF-8, and it always ends in a NUL. A nil buffer, or one of
    no length, gets nothing, and the call returns its code all the same.
  - [registry.v](../capi/registry.v) is the registry of handles (D4):
    - a map from id to item behind a `sync.Mutex`, with ids from a counter;
    - an id freed twice is "freed already", and one never given out "not a
      handle";
    - it's empty until step 2 puts databases and fluids in it.
  - On Windows the header is compiled with the module, so the C compiler checks
    each declaration against its function. A probe header with an `int32_t` id
    failed the build with "conflicting types". This can't be done on Linux
    (D3, as amended): there `int64_t` is `long` and V's `i64` is `long long`,
    so the same declarations would conflict.
- **Exports:** [tools/python.vsh](../tools/python.vsh) reads the header's
  declarations and gives the linker a `.def` file on Windows, or a version
  script on Linux.
  - Without it the DLL exported about 800 symbols; with it, the four.
  - On Linux, `-fvisibility=hidden` wouldn't have done (D3, amended). The
    collector's symbols come from an object V compiled earlier, which the flag
    can't reach.
  - Linux is checked only in CI, since this machine has no WSL.
- **Building and testing:** `tools/python.vsh build` puts the library in
  `python/reflux/_lib/`, built with gcc, with `-static-libgcc` on Windows.
  - The build is `-prod` by default (about 40 s here), or `-debug` for a quicker
    one.
  - `test` builds the library, then runs the package's tests with uv, pytest
    9.2.0 and NumPy 2.5.4 pinned. `--python` names the interpreter.
- **The package,** in [python/](../python/): `pyproject.toml`, `reflux/_native.py`
  and `reflux/__init__.py`.
  - `_native.py` loads the library (or `$REFLUX_LIB`) and declares each function
    with fixed-width ctypes.
  - It makes the first call, checks the ABI version, and maps the codes to
    `RefluxError`, `BadArgument` and `BadHandle`.
  - `reflux` itself gives `abi_version()`, `version()`, `collect()` and
    `__version__`.
  - Installing it (a `[build-system]`, and the test that installs) is step 2's,
    with D9's data.
- **The first call mustn't race** (D5, amended). The spike said V's runtime
  starts when the library loads. In fact it starts on the first call:
  - each exported function's wrapper calls `_vno_main_init_caller()`;
  - that runs `_vinit`, which starts the collector and calls
    `GC_allow_register_threads()`, once, behind a plain `static bool`;
  - the flag is set before `_vinit` runs, so a second thread goes straight on
    while the first is still starting the runtime;
  - eight Python threads making their first calls at once aborted in Boehm:
    "Threads explicit registering is not previously enabled".
  - So `_native` makes the first call at import, which Python's import lock
    serialises, and `reflux.h` tells any other caller to make one call from
    one thread first.
  - One thread's first call is fine, even from a thread that then ends: three
    runs, and a test. (On Windows. On Linux it wasn't: CI found it after
    step 3, and step 3's write-up has the fix.)
- **CI:** the gcc jobs run `tools/python.vsh test` on both OSes.
  - The spike's step stays for now, because only its thread scenario has moved
    (into `test_threads.py`). The rest need later steps' functions:
    - its handle scenarios, step 2's fluids;
    - its zero-copy, vectorised and `thermo` scenarios, step 3's properties;
    - its collector stress, step 6;
    - its per-call cost, step 7.
  - Each job now stops after 30 minutes. It took 4 to 9, but a race can hang a
    test rather than fail it (planted mistakes, below), and GitHub's own limit
    is six hours.
- **Tests:**
  - [capi/capi_test.v](../capi/capi_test.v), 6, under tcc and gcc:
    - the header against the code: its defines against the consts, and its
      declarations against the `@[export]` names;
    - the version is v.mod's;
    - a reason cut to fit at a character, at every length, with guard bytes;
    - a handle freed, then freed again, the last one given out among them, and
      ids that were never handles;
    - 20 rounds of 8 threads adding, reading and removing at once, each round
      on a map that starts empty;
    - `reflux_collect`.
  - [python/tests/](../python/tests/), 10, with each scenario that could crash
    run in a process of its own (D11):
    - the versions agree: the ABI's in Python, the library and the header, and
      the package's in v.mod and `pyproject.toml`;
    - a library of another ABI version is refused;
    - the codes are the header's;
    - the library exports exactly the header's functions;
    - a reason cut to fit at every length, without writing past the buffer;
    - ids that aren't handles, at both ends of `i64`;
    - the version needs room for its NUL;
    - `collect`;
    - 8 threads making 200 calls each at once, with collections;
    - the import from a thread that ends.
- **Found:**
  - The spike's write-up was wrong twice about Windows. V 0.5.2's collector has
    thread support there, so a Python thread must be registered, and V starts
    the library on its first call. [spike-python-ffi.md](spike-python-ffi.md)
    is corrected.
  - In V, inside `or {}`, `err` is the error even where a parameter is called
    `err`.
  - `@[has_globals]` allows `__global` only in the file it heads.
  - `v fmt -w` rewrote a file that didn't parse, dropping lines.
  - A V map read without its lock goes wrong only while another thread's
    write resizes it, so a test must make it grow (planted mistakes, below).
  - These are in [v-for-go.md](v-for-go.md), §1, §4, §11 and §12.
- **Planted mistakes:** 36, in a scratch copy of `capi/`, `python/` and
  `tools/python.vsh`. 31 were caught at once.
  - Two of those were caught by a hang, not a failure:
    - `add` without its lock: with 8 threads writing the map at once, the V
      test never ended;
    - a thread never unregistered: once Python threads had ended while still
      registered, the Python tests never ended.
    - Hence the CI jobs' timeout. The mutation script's own timeout killed `v`
      but not the test program `v` had started. That program kept the pipes
      open, so the run waited two hours on the first hang. The script now kills
      the whole tree.
  - Four survivors got tests:
    - `get` without its lock: reads in the threads, between the writes. At
      first that failed only some runs: 11 of 30 under gcc `-prod`, 18 of 30
      under tcc. It failed in the first round alone, since a read goes wrong
      only while a write resizes the map.
    - So the test runs 20 rounds, each on a map that starts empty. The mistake
      then failed 30 runs of 30 under both compilers, and the code as it is
      passed 30 of 30.
    - The "freed already" bound one short: freeing the last handle given out,
      twice.
    - No ABI check, and no first call, which skipped the check along with the
      call: a library of another ABI version is refused.
  - One is equivalent: dropping `-static-libgcc`.
    - With `-shared`, gcc links `-lgcc_s` before `-lgcc` (its specs). So a DLL
      imports `libgcc_s_seh-1.dll` only if it calls one of libgcc's helpers.
    - Neither `capi` nor the spike's library does with gcc 16.2, with `-g` or
      `-prod`. M1's spike did.
    - The flag stays for later steps' code. Should the DLL ever need that DLL,
      every Python test fails at import, since Python doesn't search PATH for a
      DLL's dependencies.
  - Planted again with the new tests, the four survivors were caught, and so
    were the registry's other five: 35 of 36, all but the equivalent one.
    `add` without its lock now fails rather than hangs.

### 2. The database, fluids and data — done
- **The data root (D9),** in `components`:
  - A data root is a directory holding `components/` and `vendor/chemsep/`.
    `components.open(root, extra)` opens one, with directories of the
    caller's own compounds after its own. The `Database` keeps its `root` and
    its `vendor` directory.
  - A fluid's kᵢⱼ banks and `mc-PR.lib` come from that `vendor` directory.
    `load_kij` and `load_mathias_copeman_pr` now take the directory, and
    `setup` passes its database's.
  - `data_root()` is `$REFLUX_DATA`, or else the source tree's `data/`, and
    `load_default()` opens it.
  - A directory that isn't a data root is named as such. One that looks like a
    components directory, as `REFLUX_DATA` named before, gets a hint to give
    its parent.
  - `load(dirs)` is gone, so nothing can open compounds without saying where
    their kᵢⱼ come from. The flowsheet file's `data` directories now follow
    the caller's database's root, not the environment's.
- **A fluid says how it was built:** `Setup` gains its equation and
  `to_toml()`. That gives the equation, the kᵢⱼ bank, each compound's CAS
  number, name, Tc, Pc, ω, molar mass and α, and each pair's kᵢⱼ with its
  source, line, comment and unused listings. Floats are V's shortest
  round-trip form, so they read back exactly.
- **The C API,** [capi/database.v](../capi/database.v), with six functions:
  - `reflux_data_root`;
  - `reflux_database_open`;
  - `reflux_database_index`, the compounds in `index.toml`'s format;
  - `reflux_compound`, a compound's file;
  - `reflux_fluid_new`;
  - `reflux_fluid_setup`, the setup's TOML.
- **Text out:** text goes into the caller's buffer, with the room it takes in
  `*needed`, or `REFLUX_TOO_SMALL` with nothing written.
- **Text in:** NUL-terminated UTF-8, checked. A list of texts is a pointer
  array and a count.
- **New codes:** `REFLUX_NOT_FOUND` (a lookup, whose reason names the
  closest), `REFLUX_DATA_ERROR` and `REFLUX_TOO_SMALL`.
- **Arguments first:** everything a caller can get wrong is checked before the
  V API is called, so what `components` refuses after is the data's.
  - The checks: the equation, α and kᵢⱼ numbers; Mathias–Copeman with
    anything but PR; no compounds; a compound named twice, two ways; and an
    override that's NaN, pairs a compound with itself, repeats a pair or
    names a compound outside the fluid.
  - A handle of the wrong kind is `REFLUX_BAD_HANDLE`, "handle 7 is a fluid,
    not a database". A fluid outlives its database.
- **No `const` in the header:** V's definitions take `char *`, and gcc refused
  a probe header's `const char *` with "conflicting types". So the header
  says in words that the library never writes through text going in.
- **The package:**
  - `reflux.Database(root=None, extra=())` has `compound(query)`, a dict read
    from the TOML with `tomllib`, `compounds()`, and `fluid(equation, names,
    alpha=, kij=, overrides=)`.
  - `reflux.Fluid` has `names`, `cas`, `nc`, `tc`, `pc`, `omega`,
    `molar_mass` and `kij`, read-only NumPy arrays, and `setup`.
  - `reflux.data_root()`, and the errors `NotFound`, a `LookupError`, and
    `DataError`.
  - Each object frees its handle on `close()`, at the end of a `with`, or by
    a `weakref.finalize` when it's garbage (D4).
  - Text comes back in a 64 KiB buffer at first, and again with the room the
    library asks for.
- **Installing (D10):** `pyproject.toml` has a setuptools build with the
  library as package data.
  - Installed from a copy into a fresh environment, the package uses its own
    library and finds the data in the source tree the library was built from:
    kᵢⱼ included.
  - `uv pip install -e python` gives the source tree's own package.
- **Tests:**
  - [capi/database_test.v](../capi/database_test.v), 7 tests: every function
    against the V API's answer for the same request; the room text needs;
    each mistake with its code and reason; and handles of the wrong kind.
  - In `components`:
    - a root needs both its directories, with the hint;
    - a moved copy is used whole: its `pr.ipd` lists 0.0313 for
      methane/n-butane where the source tree's lists 0.0244, and its
      `mc-PR.lib` has no methane;
    - `REFLUX_DATA` moves the default;
    - the setup's TOML reads back exactly;
    - an empty query finds nothing.
  - In `flowsheet`: a file's `data` follows its database's root.
  - In Python, `test_database.py` (19 tests):
    - a compound is its file in the repository, read with `tomllib`;
    - the fluid's constants are its compounds' files', and its kᵢⱼ the lines
      of `pr.ipd`;
    - the arrays are read-only;
    - every mistake raises its error;
    - closing is closing;
    - a fluid only Python held is freed;
    - your own directory replaces the root's methane;
    - text longer than the first buffer comes back whole;
    - the moved copy, opened directly and through `REFLUX_DATA`, set before
      the process and by Python after the import.
  - `test_handles.py`, the spike's handle scenarios: a fluid only Python holds
    survives 40 rounds of databases and fluids made, freed and collected
    around it.
  - `test_install.py`; and `test_threads.py`, whose 8 threads now also build
    fluids from one shared database and read its files.
  - The Python suite is 31 tests, in 43 s here: the install takes 14, and each
    thread test 11.
- **Found:**
  - An empty query matched the 35 compounds without an InChIKey, so it was
    "ambiguous" and listed them all: `find` compared the query with each
    compound's InChIKey, or `''`. An empty or blank query now finds nothing.
    It came from M2.
  - On Windows, Python's `os.environ` and the library share the UCRT's
    environment, so `REFLUX_DATA` set after the import is seen. The test
    checks it.
  - V: a map literal keyed by calls keeps one entry per distinct key, so a
    table of cases checked only its last until it became an array of structs;
    a struct can't be a map key; and `_test.v` files don't share helpers
    ([v-for-go.md](v-for-go.md), §9). V has no `const`, and a C array of
    strings is `&&char` (§11).
- **For 0.2.0's CHANGELOG entry** (step 9):
  - `REFLUX_DATA` now names the data root, not its `components/`, and the
    kᵢⱼ and Mathias–Copeman constants come from it too;
  - `components.load(dirs)` is now `components.open(root, extra)`;
  - `load_kij` and `load_mathias_copeman_pr` take a directory;
  - `Database` has `root` and `vendor`, and `Setup` has `kind` and
    `to_toml()`;
  - an empty name finds no compound.
- **Planted mistakes:** 49, in a scratch copy of the repository, across
  `components`, `flowsheet/file.v`, `capi/database.v` and the package. All
  were caught, none by a hang.
  - Three were caught by tests written for them before the run, for gaps seen
    while writing the list: a file's `data` following its database's root,
    your own directory from Python, and text longer than its first buffer.
  - The catches came from the module's own tests: `components` caught its 18,
    the flowsheet file test its one, `capi`'s V tests their 17, and the Python
    tests the package's 13.
  - **Corrected in step 3:** three of these didn't compile, so no test had
    caught them. Two left an import unused and one left a `mut` never
    changed, each an error under `-W`, and the script counted the failed
    build as caught. Rewritten to compile and rerun with the rest in step 3,
    all 49 were caught.

### 3. Properties and vapour pressure — done
- **Many points in a call** ([capi/props.v](../capi/props.v)), the
  conventions every vectorised function from here on keeps (D6, D7):
  - n points, and each per-point input with a count of 1 or n, so a scalar or
    one composition isn't repeated;
  - outputs property-major, `out[k·n + i]`, so each property is a contiguous
    array of the points;
  - a status per point, `REFLUX_POINT_*` from 0 to 8. The flash's kinds
    follow `flash.Failure`'s order, and 4 to 6 wait for step 4. A point
    without an answer gets NaN and `REFLUX_NO_PHASE`, and the others are as
    they'd be alone;
  - `err` holds the first failed point's reason, so the reason for a point is
    the same call made for it alone. Python's `reason(index)` does that.
- **`reflux_props`:**
  - the phase at each (T, P, z), by the root asked for (liquid, vapour or
    stable), with Z, V, H, S, G, U, A, Cp, Cv, w and μ_JT per mole;
  - ln φ, and the phase's label by its phase identification parameter;
  - each value is the V API's (`solve_tp`, then `Fluid.props`), divided by
    Σz where it's extensive, so a V program doing the same gets the same bits
    (criterion 3, step 7).
- **A point's inputs are checked first:** T and P finite and above zero, and z
  finite, not negative and not all zero. thermo's own check, `t <= 0.0`, lets
  NaN and ∞ through.
- **An answer that isn't finite isn't one:** a new status,
  `REFLUX_POINT_NOT_FINITE` (8). Probing the extremes found points returned
  as answers that weren't:
  - above about 1e10 K, ChemSep's ideal-gas Cp correlations overflow, so H
    and S were NaN and Cp infinite, with the status OK;
  - at 1e300 K the volume solver returned a spurious root, Z = 6e-301.
  - Both now get `NOT_FINITE`. A volume root past what a float holds (1e-300
    Pa, 1e-300 K, 1e100 Pa) is `REFLUX_POINT_NO_VOLUME_ROOT`, as thermo says.
- **`reflux_psat`:** one compound's vapour pressure at n temperatures.
  - It's the fluid's own cubic for that compound alone, with its own α, so
    Mathias–Copeman's when the fluid has it.
  - At or above Tc it's `REFLUX_POINT_SUPERCRITICAL`, and when thermo's
    iteration fails it's `REFLUX_POINT_NO_CONVERGENCE`.
  - It began as every compound in one call. Then a single-point call's reason
    was the first failed compound's, which isn't always the one asked about,
    so it became one compound per call.
- **ΔH°f and S° in the setup's TOML,** where a compound has them, so the
  package can say whether H, S, G, U and A have a formation basis
  (`Fluid.absolute`).
- **The package:**
  - `Fluid.props(T, P, z=None, root="stable")` gives a `Properties`, and
    `Fluid.psat(T, compound=None)` gives a `Saturation`: arrays of the
    points' shape, with `status`, `ok` and `reason(index)`.
  - `reflux.Status` and `reflux.Phase` are IntEnums of the header's numbers,
    and `Fluid.index` finds a compound by position, name or CAS number.
  - T, P and z broadcast by NumPy's rules, z's last axis the compounds'. An
    input already of the points' shape, C-contiguous `float64`, isn't copied:
    a test checks it shares memory.
  - A pure fluid needs no z.
- **The spike's scenarios** for zero copies, many points in one call and
  agreement with `thermo` now live in the package's tests. Its CI step stays
  for the collector stress (step 6) and the per-call cost (step 7).
- **Against `thermo`:** `test_thermo_reference.py` builds M1's fixture
  through the C API. M1's constants go in a directory of the caller's own
  compounds, then all 14 cases are compared, PR and SRK, pure and ternary,
  with their kᵢⱼ as overrides. The spike compared Z for pure n-butane alone.
  - Hʳ and Sʳ are taken as the difference from the same state at 1 mPa, where
    the residual part is below 1e-11 of the answer.
  - Worst, against `thermo`:

    | Quantity | Worst, relative |
    |---|---|
    | Z | 5.3e-15 |
    | ln φ | 3.3e-15 |
    | V | 1.8e-11 |
    | Hʳ | 7.2e-10 of RT |
    | Sʳ | 4.6e-10 of R |

    V's is about what the two programs' gas constants' last digits give; Hʳ
    and Sʳ are within M1's 1e-9, as in thermo's own fixture test.
- **Tests:**
  - [capi/props_test.v](../capi/props_test.v), 6: every point against the V
    API, bit for bit, for every root; one value for every point; each bad
    input between two good points, and alone; beyond the model; calls that
    can't be made; and Psat, each compound's own, at Tc and where thermo
    can't converge.
  - `test_props.py`, 11 tests:
    - a grid is its points one by one, bit for bit, and twice the moles is
      the same per mole to 1e-13;
    - inputs already in shape aren't copied;
    - G, A and U keep their definitions to 1e-12, Cp > Cv, and Z = Pv/RT;
    - the stable root is the one of lower G, and where there's one root every
      root is it;
    - at Psat the liquid's ln φ is the vapour's to 1e-9;
    - Psat by compound, and a pure fluid's compound axis;
    - failures, beyond the model, mistakes;
    - `absolute`, with methane of your own without ΔH°f.
  - `test_thermo_reference.py`, 15.
  - `test_library.py` now checks every number in the package against the
    header's.
  - The Python suite is 57 tests.
- **Found:**
  - The non-finite answers above.
  - **In CI after the step was pushed,** the Linux gcc job aborted in
    `test_an_import_from_a_thread_that_ends`: "Signals delivery fails
    constantly at GC #2". This had gone unseen since step 1, since this
    machine has no Linux.
    - The library's first call starts the collector in its thread, which the
      collector keeps registered as its main thread. When a Python thread
      made that call, through its import, and ended, the next collection
      signalled a dead thread.
    - Windows asks whether a thread has exited before suspending it, which is
      why the test passed here.
    - The fix, in [capi/gc_threads.h](../capi/gc_threads.h): on Linux, in the
      library, the first call's thread is unregistered when the call ends,
      as Boehm allows once for the main thread. That isn't done if it's the
      process's first thread (thread id = process id), which lives as long
      as the process.
    - Unregistering it always would have been simpler. But then Python's
      main thread registers on every call, and an empty call from it took
      2.29 µs instead of 0.37.
    - The library is told it's the library by `-DCAPI_LIBRARY`, from
      `tools/python.vsh`. A probe header that refused to compile with it
      showed that it reaches the header there, and not in capi's own tests,
      whose main thread is V's. D5 is amended.
  - thermo's `saturation_pressure`, successive substitution from Edmister's
    start, gives up at Tc(1−1e-6) for n-butane though it converges at
    Tc(1−1e-9), and doesn't converge at or below 10 K. The status says so;
    a better method is carried forward.
  - `Saturation.reason` read a one-compound fluid's `psat(T)` as one
    compound's, without its last axis. Fixed, with a test, before the commit.
  - Step 2's planted mistakes weren't all proven: three didn't compile, an
    unused import twice and a `mut` never changed, and a failed build counted
    as caught. The scripts now report a mutation that doesn't compile as NO
    BUILD. Step 2's are rerun below, those three rewritten to compile.
  - V writes 1e10 as `1e+10` and +∞ as `+inf`
    ([v-for-go.md](v-for-go.md), §8).
- **Planted mistakes:** 42 for this step, in a scratch copy, in one run with
  step 2's 49. The script now reports a mutation that doesn't compile as NO
  BUILD; none didn't.
  - 40 were caught: 28 by `capi`'s V tests and 12 by the Python tests.
  - "psat reason's compound" survived: with every compound's Psat, a reason
    was always read for the first compound, and no test asked for another's.
    A test now asks for n-butane's, and the mistake is caught.
  - "ln φ not checked finite" survived, and is kept. For ln φ alone to be
    infinite, Z would have to underflow to exactly 0, a finite value whose
    log is −∞, with every property finite. No state found reaches it: at
    1e300 K, where Z fell to 6e-301, H was already NaN. It's the same rule
    every value keeps, that one that isn't finite is no answer.

### 4. Stability and the PT flash — done
- **`reflux_stability`** ([capi/flash.v](../capi/flash.v)) tests n feeds with
  the V API's defaults: Wilson's two trial phases and a nearly pure one of
  each compound. For each feed it gives:
  - the verdict (1, 0, or −1 without an answer);
  - tm, the lowest tangent-plane distance any trial reached away from the
    feed, negative where the feed is unstable;
  - w, that trial's composition, the phase that would appear.
  - If every trial came back to the feed, as for a vapour well away from its
    dew point or a pure compound, tm is 0 and w the feed.
- **`reflux_flash(fluid, spec, …)`** is one function for every flash.
  - Each point has two values, a and b, whose meaning `spec` gives:
    `REFLUX_FLASH_PT`, (T, P), is this step's. Step 5's flashes are new specs
    with the same layout.
  - It writes the whole's β, T, P, H, S and V per mole of feed, and each
    phase's fraction, V, H and S per mole of the phase.
  - It also writes each phase's composition, and how many phases there are.
- **Each phase is in the place of its label** (D6, as amended):
  - the liquid's (0) then the vapour's (1), so `x[..., VAPOR, :]` is always
    the vapour's composition;
  - a split's less dense phase is the vapour (D6 of M3), and a single phase
    takes its label's place, with NaN in the other;
  - the whole's H, S and V are its phases', summed in the order the flash
    gives them, as `flowsheet.Stream` sums them.
- **Statuses:** a flash's failure gets the status of its kind, in
  `flash.Failure`'s order: no convergence, no volume root, or more phases
  than two. An answer that isn't finite is `NOT_FINITE`, checked on the whole's
  values. Each phase's values enter the whole's, so a phase's needs no check
  of its own.
- **The package:**
  - `Fluid.stability(T, P, z)` gives a `Stability` (`stable`, `tm`, `w`).
  - `Fluid.flash_pt(T, P, z)` gives a `FlashResult`: `beta`, `t`, `p`, `h`,
    `s`, `v` and `phases`, with `fraction` and `x` on an axis of two phases.
  - Its `liquid` and `vapor` are `PhaseView`s, each with `fraction`, `x`,
    `v`, `h`, `s` and `present`.
  - Broadcasting is `props`', through one function now shared by props,
    stability and the flashes.
- **Against `thermo`:** `test_flash.py` flashes all 83 of M3's PT states
  (`testdata/m3_thermo.json`) through the package, with M3's constants as
  the caller's own compounds.
  - The splits agree to 2.0e-12 and the phases' volumes to 1.9e-11, the same
    worst deviations `flash/thermo_test.v` measures for the V API.
  - The fixture's ideal-gas Cp isn't ChemSep's, so H and S aren't compared
    here.
  - Two of the 39 splits are ones `thermo` labels liquid–liquid. The C API
    calls the less dense phase of each the vapour (D6), and the fixture lists
    its phases by density, so they're matched by order.
- **Tests:**
  - [capi/flash_test.v](../capi/flash_test.v), 5 tests:
    - every point against the V API's flash, bit for bit;
    - a three-compound split's balances and equal fugacities, from the arrays
      and from states solved again;
    - a single phase in its label's place;
    - a failure's status: bad input, three phases (water, n-decane and
      methane with a water–hydrocarbon kᵢⱼ of 0.5), not finite, no volume
      root;
    - calls that can't be made;
    - stability against the V API, a pure compound, and a bad point.
  - `test_flash.py`, 7 tests:
    - M3's splits;
    - a grid is its points one by one, bit for bit;
    - a split's balances close and its fugacities agree, checked through
      `props`;
    - the whole is its phases, and a single phase is `props`' stable root;
    - a single phase takes its label's place;
    - stability agrees with the flash on a 9 × 9 grid;
    - failures.
  - The Python suite is 64 tests.
- **Found:** I'd written a closure (`reduce` with an anonymous function) to sum
  a composition in code that runs per point, which CLAUDE.md forbids. It's
  `numeric.sum` now. A per-phase finite check that could never fire on its own
  is deleted.
- **Planted mistakes:** 40, in `capi/flash.v` and the package, in a scratch
  copy. 36 were caught at once: 28 by `capi`'s V tests, 8 by the Python
  tests.
  - One didn't compile: V refuses `if (…)`'s parentheses. Rewritten to
    compile, it was caught.
  - "failed phase left" survived: a failed point's phase values were never
    checked to be NaN. Both test files now check them, and it's caught.
  - Two survive, and are explained:
    - **The flash's own bad input**, mapped to another status, is never
      seen. The C API checks a point's inputs before the flash does, so the
      flash's `.bad_input` never comes back. The arm stays, since V's
      `match` on an enum must cover every value.
    - **No convergence**, mapped to another status, is never seen either.
      With the defaults the C API uses, no state is known that makes the PT
      flash or stability fail to converge: M3's 10,000-point sweep found
      none, and the V tests reach it only with `tol: 1e-16` or
      `max_steps: 1`. The vapour-fraction flashes do reach it at their
      defaults (`flash/vf_test.v`), so step 5's tests can check this
      mapping.
    - **Caught in step 5:** both are now reached through the public API (an
      H no temperature reaches, and a dew point that doesn't converge), and
      both mistakes are caught.

### 5. The other flashes — done
- **Six more specs for `reflux_flash`,** each a pair of values a and b, in
  the one layout of step 4:

  | Spec | a | b |
  |---|---|---|
  | PH | P (Pa) | H (J/mol) |
  | PS | P (Pa) | S (J/(mol·K)) |
  | TV | T (K) | V (m³/mol) |
  | UV | U (J/mol) | V (m³/mol) |
  | TVF | T (K) | β |
  | PVF | P (Pa) | β |

  - Each spec checks its own two values first:
    - T and P finite and above zero;
    - H, S and U finite;
    - V finite and above zero;
    - β between 0 and 1, NaN refused.
  - Then the feed is checked, and the V API's flash called with its
    defaults.
- **Bubble and dew points** are the vapour-fraction flashes at β = 0 and 1.
  - Both phases come back, even the one with no moles. Its composition is
    the incipient phase's: the first bubble, or the first drop.
  - The package names them: `bubble_p`, `dew_p`, `bubble_t`, `dew_t`.
- **Every failure kind can now be reached,** and is tested through the
  public API. With M3's methane/n-butane, 40/60:
  - `dew_t` at 90 bar fails to converge;
  - both points at 100 bar are outside the envelope;
  - at 110 bar the dew-point iteration reaches a boundary at 26 K the wrong
    way round, so its phases are swapped.
  - An H that no temperature from 10 K to 10⁴ K reaches is the PH flash's
    own bad input.
  - Step 4 couldn't test the first and last of these mappings, and explained
    them; step 5 tests them.
- **The package:** `Fluid.flash_ph(P, H, z)`, `flash_ps(P, S, z)`,
  `flash_tv(T, V, z)`, `flash_uv(U, V, z)`, `flash_tvf(T, beta, z)`,
  `flash_pvf(P, beta, z)`, `bubble_p(T, z)`, `dew_p(T, z)`, `bubble_t(P, z)`
  and `dew_t(P, z)`. Each gives a `FlashResult`, broadcast as `props`.
- **Against `thermo`,** through the package with M3's constants:
  - its 40 bubble and dew points to 9.0e-13 in T or P and 1.3e-13 in the
    incipient phase;
  - its 36 pure saturation states to 9.1e-13.
  - These are the V API's own measured deviations (`flash/vf_test.v`).
- **Round trips:** a grid of PT states, splits and single phases, given back
  to each flash (its P and H, its T and V, and so on) lands on T, P and β to
  1e-8.
  - One of ten splits came back from TVF at 80.6 bar where the PT state was
    at 50. Near the critical region two pressures have the same β, and
    Newton finds the one its start leads to (documented in `vf.v`).
  - So the vapour-fraction flashes are held to the claim: a PT flash at
    their answer has the β asked for, to 1e-8. Nine of ten land on the state
    they came from.
- **Tests:**
  - [capi/flash_test.v](../capi/flash_test.v) now has 8 tests:
    - every spec against the V API, bit for bit, landing on the PT state it
      came from;
    - a bubble and a dew point's incipient phases;
    - each spec's checks of its values, with their reasons, a bad feed, and
      the PH flash's own bad input.
  - `test_flashes.py`, 5 tests: the round trips on arrays; M3's boundaries
    and saturation states; bubble and dew points; every failure kind.
  - `test_flash.py` and `test_flashes.py` share M3's database through
    `support.py`.
  - The Python suite is 69 tests.
- **Planted mistakes:** 23, in `capi/flash.v` and the package, in a scratch
  copy. All were caught, every one compiling: 15 by `capi`'s V tests and 8
  by the Python tests.
  - Among them are step 4's two explained survivors, the flash's own bad
    input and no convergence each mapped to another status, and the
    swapped-phases and outside-the-envelope mappings. Step 5's tests catch
    all four.
  - One gap was seen while writing the list and filled first: no test gave a
    flash other than PT a bad feed.

### 6. Threads, memory and fuzzing — done
- **Seeded fuzzing** (criterion 2): [capi/fuzz_test.v](../capi/fuzz_test.v)
  calls every exported function with mutated arguments, in process, so it
  runs in every CI job under both compilers.
  - Handles: live ones of each kind, freed ones, 0, −1, ±2⁶³.
  - Counts: right, wrong or negative.
  - Codes out of range: root, spec, equation, α, kᵢⱼ and compound.
  - Names: empty, junk, 5,000 characters, or not UTF-8.
  - Values: NaN, ±∞, ±0, denormals, ±1e308, and ordinary ones on log scales,
    mixed with plausible T, P, H, V and β, so that the flashes run too.
  - Buffers are always the size the call is told: a C caller that lies about
    a buffer can't be caught.
  - Every call must keep its contract:
    - a known code, and a known status for every point;
    - NaN for a point without an answer, and finite values for one with;
    - bad input exactly where a point's T, P (or the like) or feed break the
      header's rules, an oracle the fuzzer computes itself;
    - nothing written past a buffer, at lengths around the room each text
      needs;
    - the whole version or none of it, and a reason for a refused fluid.
  - CI runs 2,000 calls (seed 20261105) in about 10 s under tcc.
    `REFLUX_FUZZ_SEED` and `REFLUX_FUZZ_CALLS` set others: five more seeds of
    20,000 calls each, 100,000 in all with gcc `-prod`, ran clean here, each
    in under a minute.
  - It tallies each function's statuses, which shows it gets past the checks.
    Under CI's seed:
    - the flashes answered 16 points and reached no convergence, no volume
      root, more phases, outside the envelope and not finite;
    - properties answered 69 points and reached not finite and no volume root;
    - stability answered 21 and failed to converge at 11 extreme points;
    - Psat reached supercritical and no convergence.
  - Nothing panicked or hung, so the risk table's fallback (each panic
    becoming a typed error in its module) wasn't needed.
- **Threads** (criterion 7): `test_every_vectorised_function_under_threads`.
  - 8 Python threads call every vectorised function three times each: props
    at two roots, Psat, stability, and all seven flashes with bubble and dew
    points, on a 32-point grid.
  - Every array equals the serial run's bit for bit, NaN included.
  - It runs on Windows and on Linux in CI.
  - Throughput, with the `-prod` library on 8 logical CPUs, and no target:

    | Work | Serial | 8 threads | Ratio |
    |---|---|---|---|
    | props, 1,000 points a call | 25,763 points/s | 42,195 points/s | 1.6× |
    | PT flash | 626 points/s | 1,454 points/s | 2.3× |
    | PH flash | 241 points/s | 307 points/s | 1.3× |

    Each point allocates, and the collector's lock is shared, so threads add
    little. The PH flash, which allocates most, gains least. Parallel speed
    is M8's question, as the plan's probes said.
- **Memory** (criterion 6): `test_memory.py` makes 100,000 calls.
  - Most are a single point of properties, where a per-call leak would show.
    One call in a hundred each makes and frees a fluid, flashes PT, tests
    stability, looks up a compound, or flashes PH.
  - The heap, measured every 10,000 calls, stayed at 4,580 KiB throughout. A
    leak of 100 bytes a call would have added 10 MB.
  - It runs in a process of its own, in about 35 s. This replaces the
    spike's collector stress, which made 500 calls.
- **The Python suite** is 71 tests, in about two minutes here. The memory,
  install and thread tests take most of that.
  - Stale and freed handles were step 2's (`test_closed_is_closed`,
    `test_a_fluid_only_python_held_is_freed`), and step 1's ids are never
    reused.
- **The spike's CI step** keeps one scenario now, the per-call cost, which
  step 7 moves.
- **Planted mistakes:** step 6 adds tests, not library code, so its six
  mistakes are ones each new test is there to catch, run against those tests
  alone:
  - an answer not checked finite;
  - the feed unchecked;
  - text written a byte past its buffer;
  - the version a byte past its buffer;
  - a fluid that leaks;
  - threads left unregistered.
  The first run caught three. The other three showed the fuzzer was weaker
  than it should be:
  - Its contract was loose. A NaN feed that skipped its check came back as
    "not finite", which keeps NaN with a failure status, so the fuzzer
    couldn't tell. It now computes which points break the rules and holds
    them to bad input, and the feed is negative now and then.
  - Its buffer lengths never met the edge: 0, 1, 5, 64 and 65,536 are never a
    text's length exactly. It now asks each text's room first and tries a
    byte short, exact and a byte over.
  - The version's check, made one short, only cut the version, since
    `write_message` bounds itself. So nothing was written past, but the
    answer was wrong. The fuzzer now checks that an answer is the whole
    version.
  Rerun, all six were caught.

### 7. The V reference and speed — done
- **The V reference** (D12, criterion 3):
  [tools/capi_reference.v](../tools/capi_reference.v) writes the V API's
  answers for a reference set of states to JSON, each float as the hex of
  its bits.
  - It calls `components`, `thermo` and `flash`, never `capi`, and computes
    each point as `reflux.h` documents it: properties per mole as the
    extensive ones over Σz, a phase in its label's place, the whole summed in
    the flash's order, and an answer that isn't finite as none.
  - `tools/python.vsh test` builds it with the library's compiler and flags
    (gcc, with or without `-prod`, as the library), runs it, and points
    `REFLUX_REFERENCE` at what it wrote. The file is never committed.
  - `python/tests/test_reference.py` computes the same through the package
    and the library, and every bit must agree, NaN and −0 included:
    - each fluid's setup, and three compounds' files, as text;
    - 240 property points: four fluids (PR binary, SRK ternary, PR n-butane
      with Mathias–Copeman α, and a five-compound PR with a kᵢⱼ override),
      three roots, 18 states and two beyond the model;
    - 44 vapour pressures, 36 stability tests;
    - every flash: PT at nine states each, with PH, PS, TV and UV from its
      answer, and TVF and PVF at three β, failures included.
  - pytest run alone skips these, since there's no reference without
    `python.vsh`.
- **Speed** (criterion 10), in [benchmarks.md](benchmarks.md), "From
  Python".
  - The work: properties and a PT flash of methane, ethane and n-butane at
    20 bar, at 1, 100 and 10,000 points. Each is timed in V's own loop
    (`tools/bench_case.v --capi`, the work the C API does at each point) and
    from Python through the package.
  - `tools/gen_m10a_bench.py` writes `testdata/m10a_bench.json`, and a test
    checks that the committed report shows every timing.
  - The overheads measured:

    | Work | Per call (at 1 point) | Per point (100) | Per point (10,000) |
    |---|---|---|---|
    | Properties | 109 µs | 4.6 µs | 6.2 µs |
    | PT flash | 250 µs | 188 µs | 117 µs |

  - **What they are,** found by a scratch library that ran V's own loop
    inside itself, from Python:
    - A call costs 0.1 to 0.2 ms more from Python at any size: the package's
      NumPy setup and result objects.
    - The C API adds a few µs a point to properties (the checks, the copies
      out). For the PT flash it adds nothing measurable: V's loop inside the
      library took as long as the call.
    - V's code runs slower inside the library, in a Python process, than in a
      V program of its own: the flash at 100 points took 41 to 49 ms against
      27 to 31. With a 2 GB initial heap, so that the collector rarely runs,
      it was still 41 against 31 to 34. So it's more than the collections, and
      its cause isn't known: carried forward to M8, with parallel speed.
  - **Timing on this laptop moves by up to 2×** within minutes, with its
    power state, as M3 found. So each pair is timed in eight rounds, V then
    Python, about four seconds in all, which keeps the comparison fair even
    when the absolute times move.
  - The first version of V's loop skipped the phase label that the C API
    computes for each point, four extra pressure derivatives, and so made
    the C API look 2.5 µs a point dearer than it is.
- **The spike is retired** (D1): its code and its CI step are gone, and its
  write-up stays as the record.
  - [spike-python-ffi.md](spike-python-ffi.md) now lists where each scenario
    went: zero copies, `thermo`, the vectorised call, the collector, the
    handles and the thread all to the package's tests, and the per-call cost
    to the benchmarks.
  - CI's gcc jobs run `tools/python.vsh test` alone.
- **Tests:**
  - `test_reference.py`, 5 tests: text; properties; Psat; stability; every
    flash, each with answers.
  - In `benchmarks/flowsheet_test.v`, the committed report has every Python
    timing.
  - The Python suite is 76 tests.
- **Planted mistakes:** four changes of an ulp or so, which only a
  bit-for-bit comparison can see, run against `test_reference.py` alone, the
  library and the reference rebuilt each time:
  - a property scaled by (1 + 1e-15);
  - a phase's S scaled the same way;
  - the compositions normalised again;
  - V computed as V·(1/n).
  - Three were caught. "V·(1/n)" survived, because every feed's amounts
    added to 1 or 2. Dividing by a power of two is exact, so x/n and x·(1/n)
    agree to the bit.
  - The reference's second feeds now add to 3, and all four are caught.

### 8. The exit script and the Python docs — done
- **The exit script** (criterion 8):
  [python/examples/phase_diagrams.py](../python/examples/phase_diagrams.py)
  computes two phase diagrams into NumPy arrays and draws them with
  matplotlib. Both are PR with ChemSep's kᵢⱼ.
  - A T–x–y of methane and n-butane at 20 bar, below the binary's critical
    locus: 101 bubble and 101 dew points, with tie lines from bubble points to
    their first bubbles.
  - A vapour-fraction map of a natural gas (80% methane, then ethane to
    n-pentane) from 6,480 PT flashes, 150 to 350 K and 1.5 to 120 bar. That's
    past its cricondenbar, about 108 bar, and the two-phase region is outlined.
  - Its grids are exact in binary, or quotients both sides round alike: x =
    i/100, T = 150 + 2.5i K, P = 1.5j bar. So
    [tools/capi_reference.v](../tools/capi_reference.v)'s new `examples`
    section flashes the same states through the V API.
  - Its figure, [images/phase-diagrams.png](images/phase-diagrams.png), is in
    the guide. Drawn again, it's the same byte for byte (matplotlib 3.11.2).
- **Its tests,** `python/tests/test_examples.py`, 4:
  - It runs as a program, says every point has an answer, and draws a PNG of
    1,210 × 506 pixels.
  - The T–x–y is a phase diagram, by the claims alone:
    - at the pure ends, bubble and dew are the same T, the compound's
      saturation at 20 bar (its Psat there misses by 5.2e-11 for n-butane and
      3.1e-11 for methane);
    - each tie line closes: a first bubble's dew point is the same T to
      2.0e-13, and its first drop the liquid to 4.2e-13;
    - inside, the dew point is above the bubble point, by 5.1 K at least.
  - The map splits exactly where the stability test calls the feed unstable,
    at all 6,480 points (2,445 splits).
  - Its numbers are the V API's bit for bit, its inputs included: 202
    vapour-fraction flashes and 6,480 PT flashes.
- **[python.md](python.md),** the guide (criterion 11):
  - What it covers: installing, a first flash, many points at once, what a
    flash gives, properties, vapour pressure, stability, the other flashes,
    statuses and errors, fluids, closing, threads, speed, the C API from
    other languages, and the exit script.
  - Its 15 examples run in order as one session, in
    `python/tests/test_docs.py`, and what each prints must be exactly the text
    block under it. So the page can't drift from the package.
  - Its outputs were read as claims before they were kept. Two examples
    changed for it: a Psat at 150 K printed n-butane's 11 Pa as `0.` bar, and
    a stable feed's trial isn't "the phase that would appear".
- **A bug the exit script found, in M3's PT flash, fixed.** On a first grid,
  the gas failed at 265.0 and 265.1 K at 81 bar: "the split did not converge
  in 200 steps". That's well inside its two-phase region, where K runs from 1.7
  to 0.07, and 264.9 K split in 23 steps.
  - The cause: from the first start, substitution was converging, its
    residual down to 0.005 by step 5. Then the extrapolation at step 6 threw β
    from 0.14 to 0.87, and the residual up to 1.5. The split slid off towards
    the feed, β passing 1 and growing to 42,000 in the negative flash, until it
    collapsed onto it. The other starts began near the trivial solution and
    crept away too slowly to finish in 200 steps. At 264.9 K the same jump
    happened, and substitution happened to recover.
  - The split extrapolated as stability does, keeping every extrapolation
    Rachford–Rice accepts. Michelsen and Mollerup keep one only if it lowers G,
    a test stability skips on the evidence of 6,000 states
    ([design.md](design.md)).
  - The fix: an extrapolation is kept only if it lowers the residual, max
    |ln fᵢⱽ − ln fᵢᴸ|, below that of the split it was made from. Otherwise the
    plain step is taken, and the rejected one counts as a step. The residual,
    not G, because G means nothing at a β outside [0, 1], which the negative
    flash passes through.
  - What it changed:
    - on 16,160 flashes of the gas, 150 to 350 K and 1 to 160 bar, it undoes
      392 extrapolations, and the failures go from 1 (262 K, 105 bar) to 0;
    - the steps total 67,754 against 67,751, the point that failed included;
    - M3's 166 fixture splits take 2,248 steps against 2,246, and 44 Newton
      finishes as before;
    - the robustness sweep's 10,000 flashes change at one point. Point 5752
      is still refused as three phases, but by a tm of −6.643 rather than
      −4.047, since its split converged to another one before the
      verification;
    - the golden flowsheets solve to their reports, to 1e-9, as before.
  - The fix goes in 0.2.0's CHANGELOG entry. `test_an_extrapolation_that_overshoots_is_undone`
    in `flash/pt_test.v` flashes the failing states, each to criterion 3's
    equilibrium checks and a recorded 21 steps. The map's own test asks every
    point for an answer.
- **Smaller:**
  - `tools/python.vsh test` pins matplotlib 3.11.2.
  - `tools/gen_m10a_bench.py` writes LF line endings, as the other generators
    do. Git had normalised the committed file, but a regenerated one was CRLF.
- **Tests:** `test_examples.py`, 4, and `test_docs.py`, 2 (the page's
  examples, and its block finder). In V, the PT flash's new test. The Python
  suite is 82 tests.
- **Planted mistakes:** 11, all caught:
  - in the flash, no extrapolation judged, and a rejected one kept;
  - in the script, the compositions over 101, bubble and dew swapped, the
    figure at 100 dpi, and the gas reversed;
  - in the reference, the gas's composition, and its pressures one step off;
  - on the page, one printed digit, and an example's UV flash given H for U;
  - in the package, `bubble_t` at β = 1.

### 9. Release v0.2.0 — done
- **The version:** v.mod's and the package's are 0.2.0. `CHANGELOG.md` has its
  entry, undated until the tag. As for v0.1.0, the tag comes on `main` after
  step 10 and the merge, so 0.2.0 is the reviewed commit. Then a commit dates
  the entry, and `tools/release.vsh` runs there.
- **The entry:**
  - Added: the C API, the package, the same numbers as the V API, the exit
    script, and the speed from Python.
  - Changed: `REFLUX_DATA` as a data root, with the kᵢⱼ from it;
    `components`' API, as step 2 listed it; the golden reports written again.
  - Fixed: step 8's PT flash, and an empty name.
  - Known limitations: 0.1.0's, with no wheels, threads that don't scale far,
    and V's code slower inside the library.
- **What a version promises from 0.2.0 on** (D13), in the CHANGELOG:
  - a change to the C API that would break a caller is a new ABI version, and
    the package refuses a library of another;
  - the package's API changes as the V modules' does, in a minor release
    before 1.0 and never in a patch.
- **The golden reports were written again** (`tools/golden.vsh`), for the
  version each carries, and the entry names all five. Compared with v0.1.0's
  number by number, none moved, so step 8's flash fix doesn't touch their
  flashes.
- **A version's golden reports are its own now.** `flowsheet/golden_test.v`
  checks that each report's `reflux` is v.mod's version. Before, a new version
  could leave them saying the old one, and the gate wouldn't see a change to
  name.
- **The gate on this branch** refuses for exactly the reasons it should: the
  branch isn't `main`, and the entry is undated. While the step was
  uncommitted, it also refused on the tree's changes.
- **The guide's speed section** now says what a flash costs from Python: V's
  code runs 1.3 to 1.8 times slower inside the library. The entry's Known
  limitations say so too.
- **Tests:** the golden test's version check.
- **Planted mistakes:** 4, in a local clone with the step committed in it,
  since the gate reads git. All were caught:
  - a golden report left at 0.1.0, by the golden test;
  - the package at 0.1.0, by `test_library.py`;
  - a golden report not named in the entry, by the gate, which names it;
  - v.mod at 0.1.0, by the gate.
  - The third first survived, in a clone with the step uncommitted. The gate
    reads the golden reports' changes from git, between the last tag and
    `HEAD`, and they weren't in `HEAD`. A release's tree is clean, so `HEAD` is
    the tree, and committed, the gate named it.

### 10. Docs and exit review
- README, roadmap, design notes, the V notes, and a criterion-by-criterion
  review, as in M1 to M4.

## Risks

| Risk | Effect | Fallback |
|---|---|---|
| Restricting exports on Windows breaks V's own initialisation, which runs from the DLL's entry | a library that loads but doesn't start | export by `__declspec(dllexport)` through a C wrapper instead of a `.def` file; or keep gcc's export-all and test that the header's functions are present |
| Concurrent calls fail on Linux, where the probes didn't run | crashes under threads | a lock around every call (correct, if slower), and a report to V upstream |
| A flash path panics on input the V tests never gave it | Python dies | the fuzzing of step 6; each panic found becomes a typed error in the module it came from |
| A V update changes the collector or `-shared` | the library breaks | V is pinned (`.v-version`); the Python tests run in CI on every change |
| `ctypes` overhead dominates small calls | slow scalar use from Python | vectorised calls are the way in; the benchmarks say how small is too small |

## Not in M10a

- Flowsheets from Python, and units written in Python: M10b.
- Wheels with the library inside, and publishing to PyPI (I2).
- Parallel calls inside the library (`spawn`): M8's sweeps.
- A third phase in the results: M6, as a new ABI version.
- Units of measure in Python: SI only, as the V API. `quantity` could be
  exposed later.
