module capi

// Component data and fluids through the C API (M10a step 2): a database
// opened on a data root (D9), its compounds as TOML, and fluids built from
// them with their setup as TOML. The text formats are components' own, so
// that the docs of the data serve every caller.
import components
import encoding.utf8
import math
import thermo

// The codes beyond step 1's (reflux.h).
const not_found = i32(-3)
const data_error = i32(-4)
const too_small = i32(-5)

// Equations of state and α functions, as reflux.h numbers them.
const equations = {
	i32(0): thermo.CubicKind.vdw
	1:      .rk
	2:      .srk
	3:      .pr
}
const alphas = {
	i32(0): components.AlphaChoice.equation
	1:      .mathias_copeman
}

// DatabaseItem is what a database's handle refers to.
struct DatabaseItem {
	db components.Database
}

// FluidItem is what a fluid's handle refers to: the fluid, and how it was
// built.
struct FluidItem {
	setup components.Setup
}

// CallError is a failure whose code isn't the default for where it happens.
struct CallError {
	Error
	result i32
	reason string
}

fn (e CallError) msg() string {
	return e.reason
}

fn bad(reason string) IError {
	return CallError{
		result: bad_argument
		reason: reason
	}
}

// failed writes an error's reason and returns its code: a CallError's own, a
// lookup's REFLUX_NOT_FOUND, and anything else REFLUX_DATA_ERROR, which is
// what components' other failures are (a file it can't read, or a compound it
// can't build a fluid from).
fn failed(e IError, buf &char, buf_len i32) i32 {
	code := if e is CallError {
		e.result
	} else if e is components.LookupError {
		not_found
	} else {
		data_error
	}
	return fail(code, buf, buf_len, e.msg())
}

// write_text puts text in the caller's buffer, if it fits with its NUL, and
// the room it takes in *needed (reflux.h). If it doesn't fit, the reason says
// so and nothing is written.
fn write_text(text string, out &char, out_len i32, needed &i32, buf &char, buf_len i32) i32 {
	n := text.len + 1
	if needed != unsafe { nil } {
		unsafe {
			*needed = i32(n)
		}
	}
	if out == unsafe { nil } || int(out_len) < n {
		return fail(too_small, buf, buf_len, 'the text takes ${n} bytes, and the buffer has ${out_len}')
	}
	unsafe {
		vmemcpy(out, text.str, text.len)
		out[text.len] = char(0)
	}
	return ok
}

// text_in reads a caller's NUL-terminated UTF-8 text.
fn text_in(p &char, what string) !string {
	if p == unsafe { nil } {
		return bad('${what} is NULL')
	}
	s := unsafe { cstring_to_vstring(p) }
	if !utf8.validate_str(s) {
		return bad("${what} isn't UTF-8")
	}
	return s
}

// texts_in reads a caller's array of n texts.
fn texts_in(p &&char, n i32, what string) ![]string {
	if n < 0 {
		return bad('${what} has ${n} entries')
	}
	if n == 0 {
		return []string{}
	}
	if p == unsafe { nil } {
		return bad('${what} is NULL')
	}
	mut out := []string{cap: int(n)}
	for i in 0 .. int(n) {
		out << text_in(unsafe { p[i] }, '${what}[${i}]')!
	}
	return out
}

// what names an item's kind, for a handle of the wrong one.
fn what(item Item) string {
	return if item is DatabaseItem {
		'a database'
	} else if item is FluidItem {
		'a fluid'
	} else {
		'something else'
	}
}

fn database_of(id i64) !components.Database {
	item := get(id) or {
		return CallError{
			result: bad_handle
			reason: err.msg()
		}
	}
	if item is DatabaseItem {
		return item.db
	}
	return CallError{
		result: bad_handle
		reason: 'handle ${id} is ${what(item)}, not a database'
	}
}

fn setup_of(id i64) !components.Setup {
	item := get(id) or {
		return CallError{
			result: bad_handle
			reason: err.msg()
		}
	}
	if item is FluidItem {
		return item.setup
	}
	return CallError{
		result: bad_handle
		reason: 'handle ${id} is ${what(item)}, not a fluid'
	}
}

// reflux_data_root writes the data root a database opens when given none.
@[export: 'reflux_data_root']
pub fn reflux_data_root(out &char, out_len i32, needed &i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	return write_text(components.data_root(), out, out_len, needed, unsafe { nil }, 0)
}

// reflux_database_open opens a data root, with directories of the caller's
// own compounds after its own.
@[export: 'reflux_database_open']
pub fn reflux_database_open(root &char, extra &&char, n_extra i32, db &i64, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	if db == unsafe { nil } {
		return fail(bad_argument, buf, buf_len, 'db is NULL')
	}
	mut dir := ''
	if root != unsafe { nil } {
		dir = text_in(root, 'root') or { return failed(err, buf, buf_len) }
	}
	if dir == '' {
		dir = components.data_root()
	}
	dirs := texts_in(extra, n_extra, 'extra') or { return failed(err, buf, buf_len) }
	d := components.open(dir, dirs) or { return failed(err, buf, buf_len) }
	unsafe {
		*db = add(DatabaseItem{d})
	}
	return ok
}

// reflux_database_index writes every compound's index entry.
@[export: 'reflux_database_index']
pub fn reflux_database_index(db i64, out &char, out_len i32, needed &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	d := database_of(db) or { return failed(err, buf, buf_len) }
	text := components.write_index(d.list()) or { return failed(err, buf, buf_len) }
	return write_text(text, out, out_len, needed, buf, buf_len)
}

// reflux_compound writes a compound's component file.
@[export: 'reflux_compound']
pub fn reflux_compound(db i64, query &char, out &char, out_len i32, needed &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	d := database_of(db) or { return failed(err, buf, buf_len) }
	q := text_in(query, 'query') or { return failed(err, buf, buf_len) }
	rec := d.find(q) or { return failed(err, buf, buf_len) }
	text := components.write_toml(rec) or { return failed(err, buf, buf_len) }
	return write_text(text, out, out_len, needed, buf, buf_len)
}

// reflux_fluid_new builds a fluid from the database.
@[export: 'reflux_fluid_new']
pub fn reflux_fluid_new(db i64, equation i32, names &&char, nc i32, alpha i32, kij i32, pair_a &&char, pair_b &&char, pair_k &f64, n_overrides i32, fluid &i64, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := new_setup(db, equation, names, nc, alpha, kij, pair_a, pair_b, pair_k, n_overrides,
		fluid) or { return failed(err, buf, buf_len) }
	unsafe {
		*fluid = add(FluidItem{s})
	}
	return ok
}

// new_setup checks a fluid's arguments and builds it. Everything the caller
// could have got wrong is checked here, so that what setup refuses is the
// data's doing.
fn new_setup(db i64, equation i32, names &&char, nc i32, alpha i32, kij i32, pair_a &&char, pair_b &&char, pair_k &f64, n_overrides i32, fluid &i64) !components.Setup {
	if fluid == unsafe { nil } {
		return bad('fluid is NULL')
	}
	d := database_of(db)!
	kind := equations[equation] or {
		return bad('no equation ${equation}: REFLUX_VDW (0), REFLUX_RK (1), REFLUX_SRK (2) or REFLUX_PR (3)')
	}
	a := alphas[alpha] or {
		return bad('no alpha ${alpha}: REFLUX_ALPHA_EQUATION (0) or REFLUX_ALPHA_MATHIAS_COPEMAN (1)')
	}
	if a == .mathias_copeman && kind != .pr {
		return bad("ChemSep's Mathias–Copeman constants (mc-PR.lib) are fitted for PR, not ${components.equation_name(kind)}")
	}
	if kij !in [i32(0), 1] {
		return bad("kij is ${kij}: 1 for ChemSep's, 0 for none")
	}
	if nc < 1 {
		return bad('a fluid needs at least one compound, and nc is ${nc}')
	}
	// Each compound is found here, so that a compound named twice, two ways,
	// is caught, and an override's names can be matched to the list.
	mut cas := []string{}
	for q in texts_in(names, nc, 'names')! {
		rec := d.find(q)!
		if rec.ids.cas in cas {
			return bad('`${q}` is ${rec.name} (${rec.ids.cas}), which is in the list already')
		}
		cas << rec.ids.cas
	}
	if n_overrides > 0 && pair_k == unsafe { nil } {
		return bad('pair_k is NULL')
	}
	firsts := texts_in(pair_a, n_overrides, 'pair_a')!
	seconds := texts_in(pair_b, n_overrides, 'pair_b')!
	mut overrides := []components.Kij{}
	mut pairs := []string{}
	for i in 0 .. int(n_overrides) {
		k := unsafe { pair_k[i] }
		if !math.is_finite(k) {
			return bad('pair_k[${i}] is ${k}')
		}
		a1 := in_fluid(d, cas, firsts[i], 'pair_a[${i}]')!
		a2 := in_fluid(d, cas, seconds[i], 'pair_b[${i}]')!
		if a1 == a2 {
			return bad('override ${i} pairs `${firsts[i]}` with itself')
		}
		key := if a1 < a2 { '${a1} ${a2}' } else { '${a2} ${a1}' }
		if key in pairs {
			return bad('override ${i} sets the pair `${firsts[i]}`, `${seconds[i]}` again')
		}
		pairs << key
		overrides << components.Kij{
			cas1:    a1
			cas2:    a2
			k:       k
			comment: 'set by the caller'
		}
	}
	return d.setup(kind, cas, alpha: a, kij: kij == 1, overrides: overrides)!
}

// in_fluid finds an override's compound, which must be one of the fluid's.
fn in_fluid(d components.Database, cas []string, q string, what string) !string {
	rec := d.find(q)!
	if rec.ids.cas !in cas {
		return bad("${what} is ${rec.name} (${rec.ids.cas}), which isn't in the fluid")
	}
	return rec.ids.cas
}

// reflux_fluid_setup writes how a fluid was built.
@[export: 'reflux_fluid_setup']
pub fn reflux_fluid_setup(fluid i64, out &char, out_len i32, needed &i32, buf &char, buf_len i32) i32 {
	r := enter()
	defer {
		leave(r)
	}
	s := setup_of(fluid) or { return failed(err, buf, buf_len) }
	text := s.to_toml() or { return failed(err, buf, buf_len) }
	return write_text(text, out, out_len, needed, buf, buf_len)
}
