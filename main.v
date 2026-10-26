// A small demo of the thermo layer, on component data from the database.
//
//     v run .                 # n-butane at 350 K and 9.4573 bar (SVA Example 3.9)
//     v run . 300 2.5         # T in K, P in bar
module main

import math
import os
import strconv
import components
import thermo { CubicKind }

fn main() {
	t := arg_or(1, 350.0) or { die(err) }
	p := (arg_or(2, 9.4573) or { die(err) }) * 1e5
	// The built-in database: ChemSep 8.50, converted (data/components/).
	db := components.load_default() or { die(err) }
	pure_table(db, t, p) or { die(err) }
	println('')
	alpha_example(db) or { die(err) }
	println('')
	mixture_example(db) or { die(err) }
	println('')
	fluid_example(db) or { die(err) }
	println('')
	lookup_example(db)
}

fn pure_table(db components.Database, t f64, p f64) ! {
	r := db.find('n-butane')!
	butane := r.component()!
	println('${butane.name} (${r.ids.cas}) at T = ${t} K, P = ${p / 1e5} bar, constants from ${r.source}\n')
	println('${'EOS':-20} ${'Z liq':8} ${'Z vap':8} ${'V liq':9} ${'V vap':9} ${'Psat':8}')
	println('${'':-20} ${'':8} ${'':8} ${'cm3/mol':9} ${'cm3/mol':9} ${'bar':8}')
	for kind in [CubicKind.vdw, .rk, .srk, .pr] {
		eos := thermo.new_cubic(kind, [butane])!
		liq := thermo.solve_tp(eos, t, p, [1.0], .liquid)!
		vap := thermo.solve_tp(eos, t, p, [1.0], .vapor)!
		// `if` is an expression: the last line of the taken branch is its value.
		psat := if t < butane.tc {
			'${thermo.saturation_pressure(eos, t)! / 1e5:.4f}'
		} else {
			'-'
		}
		vl, vv := liq.molar_volume() * 1e6, vap.molar_volume() * 1e6
		if math.abs(vl - vv) <= 1e-9 * vv {
			println('${eos.name():-20} ${liq.z:8.5f} ${'':8} ${vl:9.1f} ${'':9} ${psat:8}  (single root)')
		} else {
			println('${eos.name():-20} ${liq.z:8.5f} ${vap.z:8.5f} ${vl:9.1f} ${vv:9.1f} ${psat:8}')
		}
	}
	vp := r.property('vapor_pressure')!
	v := vp.at(t)
	println("ChemSep's own vapour-pressure correlation: ${v.value / 1e5:.4f} bar${if v.in_range {
		''
	} else {
		' (extrapolated)'
	}}")
}

// alpha_example swaps PR's Soave α for ChemSep's Mathias–Copeman constants,
// fitted for PR, which takes one option. Soave's α is weakest for polar
// compounds, so water at its normal boiling point shows the difference.
fn alpha_example(db components.Database) ! {
	tb := 373.124 // K, water's normal boiling point (IAPWS-95)
	soave := db.fluid(.pr, ['water'])!
	mc := db.setup(.pr, ['water'], alpha: .mathias_copeman)!
	println('Water at its normal boiling point, ${tb} K (measured Psat: 1.01325 bar)')
	println('  PR, Soave α             ${thermo.saturation_pressure(soave.model, tb)! / 1e5:.4f} bar')
	println('  PR, Mathias–Copeman α   ${thermo.saturation_pressure(mc.fluid.model, tb)! / 1e5:.4f} bar  (${mc.alphas[0]})')
}

fn mixture_example(db components.Database) ! {
	t, p := 300.0, 10e5
	n := [0.3, 0.7]
	s := db.setup(.pr, ['methane', 'n-butane'])!
	pair := s.kij.pairs[0]
	println('30% methane + 70% n-butane, PR at ${t} K and ${p / 1e5} bar')
	if u := pair.used {
		println('  kij = ${pair.k}, from ${s.kij.bank} line ${u.line} (${u.comment}; ${pair.others.len} other listings)')
	}
	for root in [thermo.Root.liquid, .vapor] {
		st := thermo.solve_tp(s.fluid.model, t, p, n, root)!
		phi := st.ln_phi.map(math.exp(it))
		println('  ${root:-7} Z = ${st.z:.5f}  phi(methane) = ${phi[0]:.4f}  phi(n-butane) = ${phi[1]:.4f}')
	}
}

// fluid_example gets total properties from PR plus the database's ideal gas:
// ChemSep's Cp, and its ΔH°f and S°, so H and S are absolute, on the
// elements at 298.15 K and 1 bar.
fn fluid_example(db components.Database) ! {
	fl := db.fluid(.pr, ['n-butane'])!
	println("n-butane, PR + ChemSep's ideal-gas data, per mole")
	println('${'':-20} ${'H':9} ${'S':8} ${'Cp':8} ${'Cv':8} ${'w':8} ${'mu_JT':8}')
	println('${'':-20} ${'kJ/mol':9} ${'J/mol K':8} ${'J/mol K':8} ${'J/mol K':8} ${'m/s':8} ${'K/bar':8}')
	for label, pt in {
		'vapor  400 K  5 bar': [400.0, 5e5]
		'liquid 300 K 30 bar': [300.0, 30e5]
	} {
		root := if label.starts_with('vapor') { thermo.Root.vapor } else { thermo.Root.liquid }
		pr := fl.props(thermo.solve_tp(fl.model, pt[0], pt[1], [1.0], root)!)
		println('${label:-20} ${pr.h / 1e3:9.3f} ${pr.s:8.2f} ${pr.cp:8.2f} ${pr.cv:8.2f} ${pr.w:8.1f} ${pr.mu_jt * 1e5:8.3f}')
	}
}

// lookup_example shows what an unknown name gets: the closest matches.
fn lookup_example(db components.Database) {
	db.find('butan') or {
		if err is components.LookupError {
			println('Looking up "butan" finds nothing. Closest: ${err.suggestions[..3].join(', ')}')
		}
	}
}

// arg_or parses os.args[i] as a number, or returns `fallback` if there is
// no such argument. Beware: os.args[i].f64() would quietly turn "abc"
// into 0.0; strconv.atof64 returns an error instead.
fn arg_or(i int, fallback f64) !f64 {
	if os.args.len <= i {
		return fallback
	}
	// Wrapping an error with context, like Go's fmt.Errorf("...: %w", err).
	return strconv.atof64(os.args[i]) or { error('argument ${i} (${os.args[i]}): ${err}') }
}

// @[noreturn] lets die() end an `or { }` block that must otherwise
// produce a value, the same way exit() and panic() can.
@[noreturn]
fn die(err IError) {
	eprintln('error: ${err}')
	exit(1)
}
