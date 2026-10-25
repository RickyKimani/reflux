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

**The big difference:** V compiles to C and then calls a C compiler. On your
machine the only C compiler is the **tcc** that ships inside V (`v doctor`
shows it). tcc compiles very fast but does almost no optimisation. That's
fine while you learn. Before you benchmark anything against DWSIM, install
GCC, Clang or the MSVC Build Tools and build with `v -prod -cc gcc .` (or
`-cc clang`, `-cc msvc`). A `-prod` build is usually much faster.

To see the C that V generates for a file, run `v -o out.c file.v`.
It is a good way to find out what a feature costs.

`-W` turns warnings into errors. Use it for tests and CI (see §4 for why).

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

`for mut x in arr` changes the elements in place. In Go, `for _, x := range arr`
gives you a **copy**, which is a classic bug. V's `for mut` does what it
looks like it does. See `cubic_real_roots` in [numeric/cubic.v](../numeric/cubic.v).

V also has no shadowing. Redeclaring `x` in an inner block is a compile
error ("redefinition of `x`"), which removes another family of Go bugs.

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

---

## 8. Arrays, maps, strings

| Go                                | V                                   |
|-----------------------------------|-------------------------------------|
| `make([]float64, 0, 3)`           | `[]f64{cap: 3}`                     |
| `make([]float64, n)`              | `[]f64{len: n}`                     |
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

## 11. C interop (for LAPACK, SUNDIALS, CoolProp later)

This is V's biggest practical advantage over Go: there's no cgo tax and no
separate build mode.

```v
#flag -lCoolProp
#include "CoolPropLib.h"
fn C.PropsSI(output &char, n1 &char, v1 f64, n2 &char, v2 f64, fluid &char) f64

p := C.PropsSI(c'P', c'T', 350.0, c'Q', 0.0, c'n-Butane')
```

You won't need any of this for a while. The [`vsl`](https://github.com/vlang/vsl)
package wraps BLAS and LAPACK, and has pure-V fallbacks, if a dense solver is
ever the bottleneck.

---

## 12. Concurrency, briefly

`spawn f()` is `go f()`, and `chan f64` and `select` exist. `shared` and
`lock` are how you guard shared state. This part of V is less mature than
Go's. A single flowsheet solve is mostly sequential, but the cases of a
sensitivity analysis are independent, which makes them the first place
`spawn` will pay off.

---

## 13. A V limit that shaped this code

**Operator overloading needs both operands to be the same type.** `2.0 * x`
doesn't compile when `x` is a user-defined `Dual`. That's why the models use
hand-written derivatives rather than automatic differentiation (see
[docs/design.md](design.md)).
