# V for Go programmers, as used in this repo

Each section shows Go first, then V, and names the file in this repo where
the V version appears. Written against **V 0.5.2**. V is pre-1.0 and
changes between releases, so if something here stops compiling, check
`v version` first.

---

## 1. Toolchain

| You want to…             | Go                     | V                           |
|--------------------------|------------------------|-----------------------------|
| run                      | `go run .`             | `v run .`                   |
| build                    | `go build`             | `v .` (debug) / `v -prod .` |
| test                     | `go test ./...`        | `v -W test .`               |
| format                   | `gofmt -w .`           | `v fmt -w .`                |
| lint                     | `go vet ./...`         | `v vet .`                   |
| docs                     | `go doc thermo`        | `v doc thermo`              |
| add a dependency         | `go get …`             | `v install …`               |
| script a task            | Makefile / shell       | `v run tools/check.vsh`     |

A `.vsh` file is a V script. Top-level statements run in order, and the
`os` module is imported implicitly, so `system_args`, `chdir` and friends
need no prefix. It's how reflux avoids a Makefile that wouldn't work in
nushell or PowerShell. See [tools/check.vsh](../tools/check.vsh).

**The big difference:** V compiles to C and then calls a C compiler. On your
machine the only C compiler is the **tcc** that ships inside V (`v doctor`
shows it). tcc compiles very fast but does almost no optimisation. That's
fine while you learn. Before you benchmark anything against DWSIM, install
GCC, Clang or the MSVC Build Tools and build with `v -prod -cc gcc .` (or
`-cc clang`, `-cc msvc`). A `-prod` build is usually much faster.

To see the C that V generates for a file, run `v -o out.c file.v`.
It is a good way to find out what a feature costs.

`-W` turns warnings into errors. Use it for tests and CI (see §4 for why).

A `.vsh` script inside the project can import the project's own modules,
just as a program can. [tools/import_chemsep.vsh](../tools/import_chemsep.vsh)
is a 30-line script around `import components.chemsep`. Go would need a
`cmd/` package to do the same.

**Windows gotcha: keep the project path short.** tcc puts temporary files
with long names next to the output. In a deeply nested folder it fails with
`could not write … .o: No such file or directory`, then tries to fall back to
gcc, which you don't have. A path like `C:\Users\Ricky\dev\reflux` is fine.

---

## 2. Modules and visibility

```go
// Go: package name ≠ path; capital letter = exported
package thermo
func Solve() {}   // exported
func check() {}   // private
```

```v
// V: directory name = module name; `pub` = exported
module thermo
pub fn solve_tp() {}  // exported
fn check_tpn() {}     // private
```

- Every file in `thermo/` starts with `module thermo`. There is no
  `github.com/...` import path. `import thermo` is resolved relative to the
  folder that holds `v.mod`.
- Naming is **enforced**: functions and variables are `snake_case`, types are
  `PascalCase`. A function called `Solve` is a compile error.
- **Selective import:** `import thermo { Component, CubicKind }` lets you
  write `Component{…}` instead of `thermo.Component{…}`. See [main.v](../main.v).
- Struct fields have visibility sections instead of per-field capitals:

```v
pub struct Cubic {
pub:          // readable outside the module, not writable
	kind CubicKind
mut:          // writable inside the module only
	cache f64
pub mut:      // readable and writable from anywhere
	note string
}
```

---

## 3. Immutability by default

```go
x := 1.0
x = 2.0          // fine in Go
```

```v
x := 1.0
x = 2.0          // compile error: `x` is immutable
mut y := 1.0
y = 2.0          // fine
```

The same applies to function arguments. To mutate an argument, **both** the
callee and the caller have to say so, which makes the call site honest:

```v
fn scale(mut v []f64, k f64) { for mut x in v { x *= k } }
scale(mut volumes, 1e6)
```

`LU.solve_into(b, mut x)` in [numeric/lu.v](../numeric/lu.v) is the pattern
for iterations: the caller owns `x`, says `mut` at the call, and the solve
writes into it instead of allocating a new slice on every Newton step.

`for mut x in arr` changes the elements in place. In Go, `for _, x := range arr`
gives you a **copy**, which is a classic bug. V's `for mut` does what it
looks like it does. See `cubic_real_roots` in [numeric/cubic.v](../numeric/cubic.v).

V also has no shadowing. Redeclaring `x` in an inner block is a compile
error ("redefinition of `x`"), which removes another family of Go bugs.
Two relatives of that rule:
- A variable can't take the name of the module it's in: `components`
  inside module `components` is an error.
- A local named like a module-level const is a warning, and so an error
  under `-W`. A test const `db` clashed with the `db` inside
  `components.load` this way.

**Closures capture explicitly, and by value.** A Go closure silently
captures variables by reference, so a loop variable or a later assignment
can change what it sees. A V closure must list what it uses, and it gets a
copy taken when the closure is created:

```v
mut x := 1.0
f := fn [x] () f64 { return x }   // without [x]: compile error
x = 5.0
f()                               // 1.0, not 5.0
```

See `critical_constants` in [thermo/cubic.v](../thermo/cubic.v), which
passes `fn [s, m] (b f64) f64 { … }` to `numeric.bisect`.

`mut` in the capture list doesn't change this: `fn [mut n] () { n++ }`
increments the closure's *own* copy, and the outer `n` stays where it was.
To share state the way a Go closure does, capture a reference:

```v
struct Counter {
mut:
	n int
}

mut calls := &Counter{}
f := fn [mut calls] (x f64) f64 {
	calls.n++ // the pointer is copied; the Counter it points to is shared
	return x * x
}
```

[numeric/integrate_test.v](../numeric/integrate_test.v) counts function
evaluations this way. My first version captured an `int` and counted zero.

---

## 4. Errors and options: `!T`, `?T`, and `or {}`

This is the feature you'll use most, and it's where V differs most from Go.

```go
// Go
func SaturationPressure(m Model, t float64) (float64, error) {
	if t >= tc {
		return 0, fmt.Errorf("T above Tc")
	}
	...
}
p, err := SaturationPressure(eos, 350)
if err != nil {
	return 0, err
}
```

```v
// V: `!f64` means "an f64 or an error"
pub fn saturation_pressure(m Model, t f64) !f64 {
	if t >= c.tc {
		return error('T above Tc')
	}
	...
}
p := thermo.saturation_pressure(eos, 350.0)!           // propagate (Go's `return err`)
p := thermo.saturation_pressure(eos, 350.0) or { 0.0 }  // fall back to a default
p := thermo.saturation_pressure(eos, 350.0) or {       // handle it; `err` is in scope
	eprintln(err)
	return
}
```

- **You can't silently drop a Result.** Using the value of a `!` function
  without `!` or `or {}` is a compile error. Calling one as a bare statement,
  `f()`, is only a *warning* in V 0.5.2. Build and test with `-W`, which
  turns warnings into errors. Go lets you write `p, _ := …` without any
  complaint.
- `!` with no type (`fn check() !`) is Go's `func check() error`. See
  `Component.check` in [thermo/component.v](../thermo/component.v).
- **`?T` is an Option**: a value or `none`, with no error message. Use it
  when "nothing here" is a normal answer, not a failure. `liquid_volume` in
  [thermo/model.v](../thermo/model.v) returns `?f64`, and `none` there
  means "the liquid branch doesn't reach this pressure". `or {}` chains read
  naturally:

```v
v := liquid_volume(m, t, p, n) or { vapor_volume(m, t, p, n) or { return no_root(t, p) } }
```

- Map lookups use the same form, in place of Go's `v, ok := m[k]`:
  `want := book[kind] or { panic('missing') }`.
- **Option fields.** A struct field can be `?f64`. Leave it out of the
  literal and it's `none`; assign a plain `f64` and it's set. Read it with
  `or`. Go would use a `*float64` (and risk a nil dereference) or a
  separate `HasDhf bool`. See `IdealGas` in
  [thermo/ideal.v](../thermo/ideal.v):

```v
pub struct IdealGas {
	cp  CpForm @[required]
	dhf ?f64   // ΔH°f, J/mol, if known
}

ig := thermo.IdealGas{ cp: form }        // dhf is none
dhf := ig.dhf or { 0.0 }                 // a value either way
```
- **Typed errors** (Go's `errors.As`): make a struct that embeds `Error`
  and implements `msg()`, then test it with `is`:

```v
pub struct InputError {
	Error            // embedded: supplies default msg() and code()
pub:
	what  string
	value f64
}

pub fn (e InputError) msg() string {
	return 'invalid ${e.what}: ${e.value}'
}

thermo.new_cubic(.pr, [broken]) or {
	if err is thermo.InputError { /* … */ }
}
```

  Inside the `if err is …` branch, `err` *is* that type, so its own fields
  can be used without a cast. A failed lookup's `LookupError` carries the
  closest names, so a caller can offer them, not just print the message
  ([components/database.v](../components/database.v)):

```v
rec := db.find('watr') or {
	if err is components.LookupError {
		println(err.suggestions[0]) // Water (7732-18-5)
	}
	return
}
```

- **Wrapping** (Go's `fmt.Errorf("…: %w", err)`):
  `f() or { return error('loading ${path}: ${err}') }`. See `arg_or` in
  [main.v](../main.v). V doesn't keep an error chain the way `%w` does.
  You get the message, not something you can unwrap.
- `or {}` blocks that have to produce a value can end in `exit()`,
  `panic()` or any function marked `@[noreturn]` (see `die` in main.v).

---

## 5. Interfaces vs sum types

V has both, and reflux uses each where it fits.

**Interfaces** work as in Go: a type satisfies one implicitly by having the
methods. `Model` in [thermo/model.v](../thermo/model.v) is an interface
because the set of models is *open*: CPA, PC-SAFT or a wrapper someone
writes later must plug in without editing `thermo`.

```v
pub interface Model {
	components() []Component
	min_volume(n []f64) f64
	volume_derivs(t f64, v f64, n []f64) VolumeDerivs
	residual(t f64, v f64, n []f64) Residual
}

pub fn solve_tp(m Model, t f64, p f64, n []f64, root Root) !State
```

`Cubic` never says `implements Model`. Passing a `Cubic` to `solve_tp`
compiles because it has the five methods. Inside a function, `if m is Cubic`
checks the concrete type and smart-casts `m`, like a Go type assertion.

**An interface can hold an interface.** `Translated` in
[thermo/translated.v](../thermo/translated.v) has a field `inner Model` and
implements `Model` itself, by evaluating `inner` at a shifted volume. So
`thermo.translate(thermo.translate(eos, c1)!, c2)!` is a valid model, and
`solve_tp` can't tell it from a cubic. It's the same decorator pattern as
wrapping an `io.Reader` in Go.

**Sum types** are what Go lacks: a *closed* set of variants, checked by the
compiler. `Alpha` in [thermo/alpha.v](../thermo/alpha.v) is one, because the
module owns the full list of α functions:

```v
pub type Alpha = InvSqrtTr | Prsv | Soave | Unity

pub fn (a Alpha) eval(tr f64) (f64, f64) {
	match a {
		Unity { return 1.0, 0.0 }
		Soave { /* here `a` is a Soave, so a.m works */ }
		...
	}
}
```

- `match` must cover every variant (or have an `else`). Add `Twu{…}` and the
  compiler lists every `match` that still needs a branch. A Go type switch
  would compile and silently fall through.
- Variants carry their own data (`Prsv{kappa0, kappa1}`), and you can store
  a mix of them in one array: `alphas []Alpha`, one per component.
- You can define methods on the sum type itself. Go can't put methods on
  an interface.
- `match` and `if` are **expressions**, so they can return a value. Go has
  no ternary, but V's `x := if c { a } else { b }` does the same job.

Rule of thumb: interface when others will add implementations, sum type
when you own the list and want the compiler to check you've handled every
case.

---

## 6. Enums

```go
type CubicKind int
const ( VdW CubicKind = iota; RK; SRK; PR )
cubic.New(cubic.PR, comps)
```

```v
pub enum CubicKind { vdw rk srk pr }
thermo.new_cubic(.pr, comps)!       // the type is inferred, so `.pr` is enough
kinds := [CubicKind.vdw, .rk, .srk, .pr]  // name it once, then shorthand
```

A `match` on an enum is exhaustive too. See `CubicKind.consts()` in
[thermo/cubic.v](../thermo/cubic.v).

---

## 7. Structs

- **`@[required]` fields.** Leaving one out of a literal is a compile error.
  Go would give you `Tc = 0` without complaint. Every field of `Component`
  is required. See [thermo/component.v](../thermo/component.v).
- **Defaults:** `tol f64 = 1e-10` in the struct definition.
- **`@[params]` structs** give you named, optional arguments, which Go
  doesn't have. `new_cubic` takes a `CubicOptions`:

```v
thermo.new_cubic(.pr, comps)!                                   // no options
thermo.new_cubic(.pr, comps, kij: [0.0, 0.02, 0.02, 0.0])!       // one by name
thermo.new_cubic(.pr, [hexane], alphas: [thermo.prsv(0.2975, 0.05104)])!
```

- **Consts can be structs and arrays,** which Go doesn't allow:
  `const butane = Component{…}` and `const kinds = [CubicKind.vdw, .rk, .srk, .pr]`.
- Receivers: `fn (c Cubic) residual(…)` is a value receiver, and
  `fn (mut s Solver) step()` is the equivalent of Go's `func (s *Solver) Step()`.
- **Read-only exported fields.** In Go an exported field can always be
  written by anyone who has the struct. V separates the two: fields under
  `pub:` can be read outside the module but not assigned, and only
  `pub mut:` fields can be assigned from anywhere. The importer's `Source`,
  `Repair` and `Duplicate`
  ([components/chemsep/chemsep.v](../components/chemsep/chemsep.v)) are
  `pub:`, so a report can't be edited after the fact. (I first assumed
  `pub mut:` meant "mutable only inside the module", as older V docs said.
  A probe showed it doesn't.)
- **Struct update syntax**, which Go lacks: `Repair{...it, cas: rec.ids.cas}`
  copies `it` and changes one field.

---

## 8. Arrays, maps, strings

| Go                                | V                                   |
|-----------------------------------|-------------------------------------|
| `make([]float64, 0, 3)`           | `[]f64{cap: 3}`                     |
| `make([]float64, n)`              | `[]f64{len: n}`                     |
| make + loop filling `s[i] = f(i)` | `[]f64{len: n, init: f(index)}`     |
| `s = append(s, x)`                | `s << x`                            |
| `slices.Sort(s)`                  | `s.sort()`                          |
| `slices.Clone(s)`                 | `s.clone()`                         |
| loop + `if`                       | `s.filter(it > big_b)`              |
| loop + append                     | `s.map(it * 1e6)`                   |
| loop + early return               | `s.any(it < 0.0)`                   |
| `s[0]`, `s[len(s)-1]`             | `s.first()`, `s.last()`             |
| `[3]float64`                      | `[3]f64` (fixed size, stack)        |
| `map[CubicKind][]float64{…}`      | `{CubicKind.pr: [112.6, 2486.0]}`   |
| `fmt.Sprintf("%8.4f", x)`         | `'${x:8.4f}'`                       |
| `fmt.Sprintf("%-20s", s)`         | `'${s:-20}'`                        |
| `strconv.ParseFloat(s, 64)`       | `strconv.atof64(s)!`                |

`it` is the implicit element name inside `filter`, `map`, `any` and `all`.

**Gotcha:** arrays are **shared** when assigned, like Go slices. Use
`.clone()` when you keep a copy of an argument. `new_cubic` and `State`
clone the arrays they store, so the caller can reuse theirs.

**Gotcha:** `'abc'.f64()` returns `0.0` with no error. Always parse user input
with `strconv.atof64`, which returns `!f64`.

**Gotcha: `toml.Any` converts without complaint.** `doc.value('name').f64()`
on a string returns `0.0`, and a missing key gives a `Null` that also
converts to zero. Check the type first (`if v is f64 || v is i64`), or use
`value_opt`, which returns an error for a missing key. This is why the
component reader in M2 checks every value's type.

**Gotcha (a V bug): don't interpolate an inline `.map()`.**
`println('${m.values().map(it.int())}')` prints the array's internal
struct instead of `[2, 3]`. Assign it to a variable first and interpolate
that. Found in the M2 data-format spike, with V `eab77e7`.

**Gotcha: `int` is 64-bit, but C's `int` is 32-bit.** Like Go's `int`, V's
`int` is 64-bit on 64-bit platforms (`sizeof(int)` is 8 in this V). In any
function that C, or Python's `ctypes`, calls, use explicit widths: `i32`
or `i64`. The Python spike found this the hard way. A `status &int` array
filled from a 4-byte NumPy buffer corrupted the heap
([docs/spike-python-ffi.md](spike-python-ffi.md)). Go programmers who've
used cgo know this as `C.int` versus `int`.

**Gotcha:** the string format spec has no `>`. `'${s:>8}'` prints `s` without
padding. Right-align is `'${s:8}'` and left-align is `'${s:-8}'`.

---

## 9. Tests

- Test files end in `_test.v` and test functions are named `test_*`. They
  run with `v -W test .`.
- `assert cond, 'message'` replaces `t.Errorf`. A failing assert prints both
  sides of the comparison for you.
- `!` works inside test functions, and an error fails the test.
- **Internal vs external** works like Go's `package x` and `package x_test`:
  - [numeric/cubic_test.v](../numeric/cubic_test.v) starts with
    `module numeric`, so it can call the private `eval`.
  - [thermo/cubic_test.v](../thermo/cubic_test.v) has no `module` line and
    imports `thermo`, so it can only see `pub` names. That keeps the public
    API honest.
- `math.tolerance(a, b, tol)` is the relative-or-absolute float compare.
- **An expensive fixture, computed once.** A `const` whose value is a
  function call is evaluated when the program starts, so every test function
  in the file shares it. Go would need `TestMain` or a `sync.Once`.
  [components/chemsep/chemsep_test.v](../components/chemsep/chemsep_test.v)
  imports the whole ChemSep database this way, once (about 4 s), for 16
  tests:

```v
const imported = import_files(@VMODROOT, vendored) or { panic(err) }
```
- The thermo tests are mostly **identities**, not stored numbers: derivatives
  against finite differences, Euler homogeneity, Σ nᵢ ln φᵢ = Gʳ/RT, and the
  Newton solver against the polynomial roots. Identities catch bugs that
  matching a single textbook value would miss.

---

## 10. Memory and performance

- By default V uses a garbage collector (Boehm), so your mental model from Go
  carries over. `-autofree` and `-gc none` exist, but you don't need them.
- Dynamic arrays `[]f64` live on the heap. `[3]f64` lives on the stack. In
  hot loops such as flash iterations or column stages, allocate once and reuse.
  This is the same advice as in Go.
- Passing a struct as an interface (`solve_tp(eos, …)`) boxes it, as in Go.
  That's cheap compared with a volume solve, but worth knowing.
- **Measure with `-prod` and a real C compiler** (see §1). tcc numbers
  don't tell you anything about speed.

---

## 11. C interop, both directions

This is V's biggest practical advantage over Go: there's no cgo tax and no
separate build mode. The Python spike ([spike/python_ffi/](../spike/python_ffi))
does it both ways.

**V calling C.** `#flag` adds compiler flags, `#include` pulls in a header,
and `fn C.name(...)` declares what you'll call. `@VEXEROOT` is V's install
directory and `@VMODROOT` the folder holding `v.mod`. From
[ffi/ffi.v](../spike/python_ffi/ffi/ffi.v):

```v
#flag -I @VEXEROOT/thirdparty/libgc/include
#flag -I @VMODROOT/spike/python_ffi/ffi
#include "gc_threads.h"

fn C.reflux_enter() int

reg := C.reflux_enter()
```

When the C side needs the preprocessor (`#if defined(GC_THREADS)`), write a
small C header and call it, as `gc_threads.h` does. V code can't see C
macros.

**C (or Python) calling V.** `@[export: 'name']` gives a function a fixed C
symbol, and `v -shared` builds a `.dll` or `.so`:

```v
@[export: 'reflux_scale']
pub fn reflux_scale(x &f64, n i32, k f64) {
	for i in 0 .. n {
		unsafe {
			x[i] *= k // pointer indexing is allowed only inside unsafe
		}
	}
}
```

```bash
v -cc gcc -shared -enable-globals -o reflux_ffi spike/python_ffi/ffi
```

What the spike learned, in V terms:
- **`int` is 64-bit.** Use `i32` or `i64` at the boundary (§8).
- **Globals** (`__global pinned = []voidptr{}`) need `-enable-globals`. V
  discourages them, but a registry of handles that the garbage collector
  must see is a legitimate use.
- **A module's `fn init()`** runs when the library loads, like a Go package's
  `init()`.
- **`defer { … }`** runs at function exit, as in Go, except that it takes a
  block rather than a call.
- **tcc can't build usable shared libraries** on Windows or Linux; use gcc.

For dense linear algebra beyond the built-in LU, the
[`vsl`](https://github.com/vlang/vsl) package wraps BLAS and LAPACK, with
pure-V fallbacks.

---

## 12. Concurrency, briefly

`spawn f()` is `go f()`, and `chan f64` and `select` exist. `shared` and
`lock` are how you guard shared state. This part of V is less mature than
Go's. A single flowsheet solve is mostly sequential, but the cases of a
sensitivity analysis are independent, which makes them the first place
`spawn` will pay off.

---

## 13. Generics, JSON and reflection

Generics use square brackets, as in Go 1.18+, and are monomorphised (one
copy per type). V's `json2` decodes straight into a struct by field name.
Together they make the fixture loader in
[thermo/fixtures_test.v](../thermo/fixtures_test.v) one function for every
fixture type:

```v
fn read_fixture[T](name string) !T {
	text := os.read_file(os.join_path(os.dir(@FILE), '..', 'testdata', name))!
	return json2.decode[T](text)!
}

fx := read_fixture[ThermoFixture]('m1_thermo.json')!
```

`@FILE` is the path of the source file, filled in at compile time, so the
test finds `testdata/` from wherever it's run. (Note: V 0.5's old `json`
module is gone; it's `json2` now.)

A key that may be missing from the JSON fits an Option field. Where Go
would use a `*float64` that's nil when absent, V's `json2` leaves a `?f64`
as `none`, and reading it forces you to handle that case. From
[tdep/fixtures_test.v](../tdep/fixtures_test.v):

```v
struct Point {
	t            f64
	chemicals_y  ?f64 // only for the DIPPR forms
}

if want := p.chemicals_y {
	assert agree(corr.eval(p.t), want, 1e-12, 0.0)
}
```

A JSON key can't be a V keyword such as `int`. The fixture generator writes
`integral` instead.

**Reflection happens at compile time.** Go's `reflect` inspects values at
run time, through `interface{}`, and a mistake is a run-time panic. V's
`$for` loops over a type's fields while compiling, and generates ordinary
typed code for each one. [components/record.v](../components/record.v)
sets a struct's float field by its name this way:

```v
fn set_field[T](mut s T, key string, v f64) ! {
	$for f in T.fields {
		$if f.typ is ?f64 {
			if f.name == key {
				s.$(f.name) = v // becomes s.tc = v, s.pc = v, … : one branch per field
				return
			}
		}
	}
	return error('no ${typeof(s).name} key `${key}`')
}
```

`$if f.typ is ?f64` is decided at compile time, so fields of other types
never reach the assignment, and the assignment is type-checked like any
other. The component file's keys are simply the struct's field names. The
importer, and the TOML writer and reader still to come, share one
definition and can't drift apart.

Two smaller gotchas from the importer:
- `strconv.atof64` rejects `abc` and `1.2.3` but accepts `nan` and `inf`,
  so check `math.is_finite` yourself.
- `arr.filter(…)` and `arr.map(…)` give you `it` but not `index`, which
  exists only in array initialisers (`[]int{len: n, init: index}`). Use a
  loop.

Two limits of `$for`, both met in the component-file writer
([components/file.v](../components/file.v)):
- `s.$(f.name)` can't go inside a string interpolation. Assign it to a
  variable first.
- In V 0.5.2, `if d := s.$(f.name) { … }` unwraps an Option *float* field,
  but leaves an Option *struct* field (`?CorrelationData`) still an Option
  inside the block. Pass the Option to an ordinary typed function and
  unwrap it there.

**Struct equality is deep.** In Go, `==` on a struct that holds a map or a
slice doesn't compile, and you reach for `reflect.DeepEqual`. V's `==`
compares structs field by field, through arrays, maps and Options. The
round-trip test is just `assert read_toml(path)! == record`, for each of
642 compounds.

V's `toml` module has quirks to code around:
- `nan` reads back as the *string* "nan";
- `inf` reads back as the integer 9218868437227405312;
- integers come back as `i64`, `int` or `u64` depending on their size.

[components/file.v](../components/file.v) checks every value's type, and
takes an integer as a float only when it's exact (|i| ≤ 2⁵³).

---

## 14. A V limit that shaped this code

**Operator overloading needs both operands to be the same type.** `2.0 * x`
doesn't compile when `x` is a user-defined `Dual`. That's why the models use
hand-written derivatives rather than automatic differentiation (see
[docs/design.md](design.md)).
