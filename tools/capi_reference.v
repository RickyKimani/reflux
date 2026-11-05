// capi_reference writes the V API's answers for a reference set of states,
// for the Python tests to compare with the C API's bit for bit (M10a, D12 and
// exit criterion 3). It calls components, thermo and flash directly, never
// capi, and computes each point as reflux.h documents: properties per mole as
// the extensive ones over Σz, a phase in the place of its label, the whole
// summed over the phases in the flash's order, and an answer that isn't
// finite as none. tools/python.vsh builds it with the library's compiler and
// flags, runs it, and points the tests at what it wrote:
//
//     v -cc gcc -prod -o capi_reference tools/capi_reference.v
//     capi_reference reference.json
//
// Floats are written as the hex of their bits, so that NaN and -0 compare too.
// It's never committed: bits differ between compilers and maths libraries.
module main

import components
import flash
import math
import numeric
import os
import strings
import thermo

// A fluid of the reference set, as the Python side builds it again.
struct FluidCase {
	equation  string
	names     []string
	alpha     string
	kij       bool
	overrides [][]string // pairs of names
	ks        []f64      // their kij
}

const fluids = [
	FluidCase{
		equation: 'PR'
		names:    ['methane', 'n-butane']
		alpha:    'equation'
		kij:      true
	},
	FluidCase{
		equation: 'SRK'
		names:    ['methane', 'ethane', 'n-butane']
		alpha:    'equation'
		kij:      true
	},
	FluidCase{
		equation: 'PR'
		names:    ['n-butane']
		alpha:    'mathias_copeman'
		kij:      true
	},
	FluidCase{
		equation:  'PR'
		names:     ['methane', 'ethane', 'propane', 'n-butane', 'nitrogen']
		alpha:     'equation'
		kij:       true
		overrides: [['methane', 'nitrogen']]
		ks:        [0.03]
	},
]

const equations = {
	'PR':  thermo.CubicKind.pr
	'SRK': .srk
}

const ts = [200.0, 300.0, 450.0]
const ps = [1e5, 20e5, 80e5]

// Statuses and labels, as reflux.h numbers them.
const ok = 0
const supercritical = 7
const not_finite = 8
const no_volume_root = 3
const no_convergence = 2

fn status_of(e IError) int {
	if e is flash.FlashError {
		return int(e.kind) + 1
	}
	return no_convergence
}

fn bits(x f64) string {
	return '"0x${math.f64_bits(x).hex()}"'
}

fn bits_list(xs []f64) string {
	return '[' + xs.map(bits(it)).join(', ') + ']'
}

fn nans(n int) []f64 {
	return []f64{len: n, init: math.nan()}
}

fn all_finite(xs []f64) bool {
	return xs.all(math.is_finite(it))
}

// quote writes s as a JSON string.
fn quote(s string) string {
	mut b := strings.new_builder(s.len + 2)
	b.write_u8(`"`)
	for c in s {
		match c {
			`"` { b.write_string('\\"') }
			`\\` { b.write_string('\\\\') }
			`\n` { b.write_string('\\n') }
			`\r` { b.write_string('\\r') }
			`\t` { b.write_string('\\t') }
			else {
				if c < 0x20 {
					b.write_string('\\u00${c.hex()}')
				} else {
					b.write_u8(c)
				}
			}
		}
	}
	b.write_u8(`"`)
	return b.str()
}

// feeds are the compositions a fluid of nc compounds is tried at: one that
// adds up to 1, and three times as much of another, so that per mole means
// per mole. A total of 2 wouldn't do: dividing by a power of two is exact, so
// x/n and x·(1/n) would agree to the bit (step 7's planted mistakes).
fn feeds(nc int) [][]f64 {
	return match nc {
		1 { [[1.0], [3.0]] }
		2 { [[0.5, 0.5], [2.25, 0.75]] }
		3 { [[0.3, 0.3, 0.4], [0.9, 0.6, 1.5]] }
		else { [[0.4, 0.2, 0.15, 0.1, 0.15], [1.2, 0.3, 0.6, 0.45, 0.45]] }
	}
}

// props_point is one point of reflux_props.
fn props_point(f thermo.Fluid, root thermo.Root, t f64, p f64, z []f64) string {
	nc := z.len
	mut status := ok
	mut values := nans(11)
	mut ln_phi := nans(nc)
	mut phase := -1
	if st := thermo.solve_tp(f.model, t, p, z, root) {
		pr := f.props(st)
		n := numeric.sum(z)
		v := [st.z, st.v / n, pr.h / n, pr.s / n, pr.g / n, pr.u / n, pr.a / n, pr.cp / n, pr.cv / n,
			pr.w, pr.mu_jt]
		if all_finite(v) && all_finite(st.ln_phi) {
			values, ln_phi = v.clone(), st.ln_phi.clone()
			phase = if thermo.phase_label(f.model, st) == .liquid { 0 } else { 1 }
		} else {
			status = not_finite
		}
	} else {
		status = no_volume_root
	}
	return '{"root": "${root}", "t": ${bits(t)}, "p": ${bits(p)}, "z": ${bits_list(z)}, "status": ${status}, "values": ${bits_list(values)}, "ln_phi": ${bits_list(ln_phi)}, "phase": ${phase}}'
}

fn psat_point(m thermo.Cubic, j int, t f64) string {
	pure := thermo.new_cubic(m.kind, [m.comps[j]], alphas: [m.alphas[j]]) or { panic(err) }
	mut status, mut p := ok, math.nan()
	if t >= m.comps[j].tc {
		status = supercritical
	} else {
		p = thermo.saturation_pressure(pure, t) or {
			status = no_convergence
			math.nan()
		}
	}
	return '{"component": ${j}, "t": ${bits(t)}, "status": ${status}, "p": ${bits(p)}}'
}

fn stability_point(m thermo.Model, t f64, p f64, z []f64) string {
	st := flash.stability(m, t, p, z) or {
		return '{"t": ${bits(t)}, "p": ${bits(p)}, "z": ${bits_list(z)}, "status": ${status_of(err)}, "stable": -1, "tm": ${bits(math.nan())}, "w": ${bits_list(nans(z.len))}}'
	}
	mut best := -1
	for i, tr in st.trials {
		if !tr.trivial && (best < 0 || tr.tm < st.trials[best].tm) {
			best = i
		}
	}
	total := numeric.sum(z)
	tm, w := if best < 0 { 0.0, z.map(it / total) } else { st.trials[best].tm, st.trials[best].w }
	stable := if st.stable { 1 } else { 0 }
	return '{"t": ${bits(t)}, "p": ${bits(p)}, "z": ${bits_list(z)}, "status": ${ok}, "stable": ${stable}, "tm": ${bits(tm)}, "w": ${bits_list(w)}}'
}

// FlashPoint is one point of reflux_flash, and with an answer the PT values
// the other flashes are given.
struct FlashPoint {
	json string
	ok   bool
	h    f64
	s    f64
	v    f64
	p    f64
}

fn run_flash(spec string, f thermo.Fluid, a f64, b f64, z []f64) !flash.FlashResult {
	return match spec {
		'pt' { flash.pt(f.model, a, b, z)! }
		'ph' { flash.ph(f, a, b, z)! }
		'ps' { flash.ps(f, a, b, z)! }
		'tv' { flash.tv(f.model, a, b, z)! }
		'uv' { flash.uv(f, a, b, z)! }
		'tvf' { flash.tvf(f.model, a, b, z)! }
		else { flash.pvf(f.model, a, b, z)! }
	}
}

fn flash_point(spec string, f thermo.Fluid, a f64, b f64, z []f64) FlashPoint {
	nc := z.len
	head := '"spec": "${spec}", "a": ${bits(a)}, "b": ${bits(b)}, "z": ${bits_list(z)}'
	failed := fn [head, nc] (status int) FlashPoint {
		return FlashPoint{
			json: '{${head}, "status": ${status}, "out": ${bits_list(nans(6))}, "liquid": ${bits_list(nans(4))}, "vapor": ${bits_list(nans(4))}, "x_liquid": ${bits_list(nans(nc))}, "x_vapor": ${bits_list(nans(nc))}, "n_phases": 0}'
		}
	}
	r := run_flash(spec, f, a, b, z) or { return failed(status_of(err)) }
	mut phase := [nans(4), nans(4)]
	mut x := [nans(nc), nans(nc)]
	mut h, mut s, mut v := 0.0, 0.0, 0.0
	for ph in r.phases {
		pr := f.props(ph.state)
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
		v += ph.fraction * ph.state.v
		j := if ph.label == .liquid { 0 } else { 1 }
		phase[j] = [ph.fraction, ph.state.v, pr.h, pr.s]
		x[j] = ph.state.n.clone()
	}
	out := [r.beta, r.t, r.p, h, s, v]
	if !all_finite(out) {
		return failed(not_finite)
	}
	return FlashPoint{
		json: '{${head}, "status": ${ok}, "out": ${bits_list(out)}, "liquid": ${bits_list(phase[0])}, "vapor": ${bits_list(phase[1])}, "x_liquid": ${bits_list(x[0])}, "x_vapor": ${bits_list(x[1])}, "n_phases": ${r.phases.len}}'
		ok:   true
		h:    h
		s:    s
		v:    v
		p:    r.p
	}
}

fn fluid_json(db components.Database, fc FluidCase) !string {
	mut cas := []string{}
	for name in fc.names {
		cas << db.find(name)!.ids.cas
	}
	mut overrides := []components.Kij{}
	for i, pair in fc.overrides {
		overrides << components.Kij{
			cas1:    db.find(pair[0])!.ids.cas
			cas2:    db.find(pair[1])!.ids.cas
			k:       fc.ks[i]
			comment: 'set by the caller'
		}
	}
	alpha := if fc.alpha == 'mathias_copeman' {
		components.AlphaChoice.mathias_copeman
	} else {
		components.AlphaChoice.equation
	}
	s := db.setup(equations[fc.equation] or { thermo.CubicKind.pr }, cas,
		alpha:     alpha
		kij:       fc.kij
		overrides: overrides
	)!
	f := s.fluid
	m := f.model as thermo.Cubic
	nc := cas.len
	mut props := []string{}
	mut stab := []string{}
	mut flashes := []string{}
	for root in [thermo.Root.liquid, .vapor, .stable] {
		for z in feeds(nc) {
			for t in ts {
				for p in ps {
					props << props_point(f, root, t, p, z)
				}
			}
		}
		// Beyond the model: no finite answer, and no volume root.
		props << props_point(f, root, 1e10, 1e5, feeds(nc)[0])
		props << props_point(f, root, 300.0, 1e-300, feeds(nc)[0])
	}
	mut psat := []string{}
	for j in 0 .. nc {
		for t in [150.0, 250.0, 350.0, 450.0] {
			psat << psat_point(m, j, t)
		}
	}
	z := feeds(nc)[0]
	for t in ts {
		for p in ps {
			stab << stability_point(f.model, t, p, z)
			pt := flash_point('pt', f, t, p, z)
			flashes << pt.json
			if pt.ok {
				flashes << flash_point('ph', f, p, pt.h, z).json
				flashes << flash_point('ps', f, p, pt.s, z).json
				flashes << flash_point('tv', f, t, pt.v, z).json
				flashes << flash_point('uv', f, pt.h - pt.p * pt.v, pt.v, z).json
			}
		}
	}
	for beta in [0.0, 0.5, 1.0] {
		for t in [250.0, 300.0] {
			flashes << flash_point('tvf', f, t, beta, z).json
		}
		for p in [5e5, 20e5] {
			flashes << flash_point('pvf', f, p, beta, z).json
		}
	}
	flashes << flash_point('pt', f, 300.0, 1e-300, z).json
	mut fields := []string{}
	fields << '"equation": ${quote(fc.equation)}'
	fields << '"names": [${fc.names.map(quote(it)).join(', ')}]'
	fields << '"alpha": ${quote(fc.alpha)}'
	fields << '"kij": ${fc.kij}'
	fields << '"overrides": [${fc.overrides.map('[${quote(it[0])}, ${quote(it[1])}]').join(', ')}]'
	fields << '"ks": ${bits_list(fc.ks)}'
	fields << '"setup": ${quote(s.to_toml()!)}'
	fields << '"props": [\n    ${props.join(',\n    ')}]'
	fields << '"psat": [\n    ${psat.join(',\n    ')}]'
	fields << '"stability": [\n    ${stab.join(',\n    ')}]'
	fields << '"flash": [\n    ${flashes.join(',\n    ')}]'
	return '{\n  ${fields.join(',\n  ')}\n}'
}

fn main() {
	if os.args.len != 2 {
		eprintln('usage: capi_reference <file to write>')
		exit(2)
	}
	db := components.load_default() or {
		eprintln(err.msg())
		exit(1)
	}
	mut parts := []string{}
	for fc in fluids {
		parts << fluid_json(db, fc) or {
			eprintln('${fc.names}: ${err.msg()}')
			exit(1)
		}
	}
	mut compounds := []string{}
	for name in ['methane', 'water', 'n-butane'] {
		rec := db.find(name) or { panic(err) }
		compounds << '${quote(name)}: ${quote(components.write_toml(rec) or { panic(err) })}'
	}
	how := $if prod { 'prod' } $else { 'debug' }
	text := '{"built": ${quote(how)}, "fluids": [\n${parts.join(',\n')}\n], "compounds": {${compounds.join(', ')}}}\n'
	os.write_file(os.args[1], text) or {
		eprintln(err.msg())
		exit(1)
	}
}
