// reflux demo: the thermo and flash layers, on component data from the
// database.
//
//     reflux demo            # n-butane at 350 K and 9.4573 bar (SVA Example 3.9)
//     reflux demo 300 2.5    # T in K, P in bar
module main

import math
import strconv
import components
import flash
import thermo { CubicKind }

fn demo(args []string, mut c Console) int {
	if args.len > 2 {
		c.complain('reflux demo: it takes T and P, not ${args.join(' ')}\n\n${usage}')
		return 2
	}
	t := number(args, 0, 350.0) or {
		c.complain('reflux demo: ${err}\n\n${usage}')
		return 2
	}
	p := (number(args, 1, 9.4573) or {
		c.complain('reflux demo: ${err}\n\n${usage}')
		return 2
	}) * 1e5
	examples(t, p, mut c) or {
		c.complain('reflux demo: ${err}')
		return 1
	}
	return 0
}

fn examples(t f64, p f64, mut c Console) ! {
	// The built-in database: ChemSep 8.50, converted (data/components/).
	db := components.load_default()!
	pure_table(db, t, p, mut c)!
	c.say('')
	alpha_example(db, mut c)!
	c.say('')
	mixture_example(db, mut c)!
	c.say('')
	fluid_example(db, mut c)!
	c.say('')
	flash_example(db, mut c)!
	c.say('')
	lookup_example(db, mut c)
}

fn pure_table(db components.Database, t f64, p f64, mut c Console) ! {
	r := db.find('n-butane')!
	butane := r.component()!
	c.say('${butane.name} (${r.ids.cas}) at T = ${t} K, P = ${p / 1e5} bar, constants from ${r.source}\n')
	c.say('${'EOS':-20} ${'Z liq':8} ${'Z vap':8} ${'V liq':9} ${'V vap':9} ${'Psat':8}')
	c.say('${'':-20} ${'':8} ${'':8} ${'cm3/mol':9} ${'cm3/mol':9} ${'bar':8}')
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
			c.say('${eos.name():-20} ${liq.z:8.5f} ${'':8} ${vl:9.1f} ${'':9} ${psat:8}  (single root)')
		} else {
			c.say('${eos.name():-20} ${liq.z:8.5f} ${vap.z:8.5f} ${vl:9.1f} ${vv:9.1f} ${psat:8}')
		}
	}
	vp := r.property('vapor_pressure')!
	v := vp.at(t)
	c.say("ChemSep's own vapour-pressure correlation: ${v.value / 1e5:.4f} bar${if v.in_range {
		''
	} else {
		' (extrapolated)'
	}}")
}

// alpha_example swaps PR's Soave α for ChemSep's Mathias–Copeman constants,
// fitted for PR, which takes one option. Soave's α is weakest for polar
// compounds, so water at its normal boiling point shows the difference.
fn alpha_example(db components.Database, mut c Console) ! {
	tb := 373.124 // K, water's normal boiling point (IAPWS-95)
	soave := db.fluid(.pr, ['water'])!
	mc := db.setup(.pr, ['water'], alpha: .mathias_copeman)!
	c.say('Water at its normal boiling point, ${tb} K (measured Psat: 1.01325 bar)')
	c.say('  PR, Soave α             ${thermo.saturation_pressure(soave.model, tb)! / 1e5:.4f} bar')
	c.say('  PR, Mathias–Copeman α   ${thermo.saturation_pressure(mc.fluid.model, tb)! / 1e5:.4f} bar  (${mc.alphas[0]})')
}

fn mixture_example(db components.Database, mut c Console) ! {
	t, p := 300.0, 10e5
	n := [0.3, 0.7]
	s := db.setup(.pr, ['methane', 'n-butane'])!
	pair := s.kij.pairs[0]
	c.say('30% methane + 70% n-butane, PR at ${t} K and ${p / 1e5} bar')
	if u := pair.used {
		c.say('  kij = ${pair.k}, from ${s.kij.bank} line ${u.line} (${u.comment}; ${pair.others.len} other listings)')
	}
	for root in [thermo.Root.liquid, .vapor] {
		st := thermo.solve_tp(s.fluid.model, t, p, n, root)!
		phi := st.ln_phi.map(math.exp(it))
		c.say('  ${root:-7} Z = ${st.z:.5f}  phi(methane) = ${phi[0]:.4f}  phi(n-butane) = ${phi[1]:.4f}')
	}
}

// fluid_example gets total properties from PR plus the database's ideal gas:
// ChemSep's Cp, and its ΔH°f and S°, so H and S are absolute, on the
// elements at 298.15 K and 1 bar.
fn fluid_example(db components.Database, mut c Console) ! {
	fl := db.fluid(.pr, ['n-butane'])!
	c.say("n-butane, PR + ChemSep's ideal-gas data, per mole")
	c.say('${'':-20} ${'H':9} ${'S':8} ${'Cp':8} ${'Cv':8} ${'w':8} ${'mu_JT':8}')
	c.say('${'':-20} ${'kJ/mol':9} ${'J/mol K':8} ${'J/mol K':8} ${'J/mol K':8} ${'m/s':8} ${'K/bar':8}')
	for label, pt in {
		'vapor  400 K  5 bar': [400.0, 5e5]
		'liquid 300 K 30 bar': [300.0, 30e5]
	} {
		root := if label.starts_with('vapor') { thermo.Root.vapor } else { thermo.Root.liquid }
		pr := fl.props(thermo.solve_tp(fl.model, pt[0], pt[1], [1.0], root)!)
		c.say('${label:-20} ${pr.h / 1e3:9.3f} ${pr.s:8.2f} ${pr.cp:8.2f} ${pr.cv:8.2f} ${pr.w:8.1f} ${pr.mu_jt * 1e5:8.3f}')
	}
}

// flash_example splits a mixture into its two phases at (T, P), finds its
// bubble point at that P, and flashes back to T from the split's H.
fn flash_example(db components.Database, mut c Console) ! {
	fl := db.fluid(.pr, ['methane', 'n-butane'])!
	z := [0.4, 0.6]
	t, p := 250.0, 30e5
	c.say('40% methane + 60% n-butane, PR, flashed at ${t} K and ${p / 1e5} bar')
	r := flash.pt(fl.model, t, p, z)!
	if r.phases.len == 2 {
		y, x := r.phases[0].state.n, r.phases[1].state.n
		c.say('  PT: vapour fraction ${r.beta:.4f}, methane ${y[0]:.4f} in the vapour and ${x[0]:.4f} in the liquid')
	}
	bubble := flash.pvf(fl.model, p, 0.0, z)!
	c.say('  bubble point at ${p / 1e5} bar: ${bubble.t:.2f} K')
	mut h := 0.0
	for ph in r.phases {
		h += ph.fraction * fl.props(ph.state).h
	}
	back := flash.ph(fl, p, h, z)!
	c.say("  PH from the split's H (${h / 1e3:.3f} kJ/mol): back to ${back.t:.6f} K")
}

// lookup_example shows what an unknown name gets: the closest matches.
fn lookup_example(db components.Database, mut c Console) {
	db.find('butan') or {
		if err is components.LookupError {
			c.say('Looking up "butan" finds nothing. Closest: ${err.suggestions[..3].join(', ')}')
		}
	}
}

// number parses args[i] as a number above zero, or returns `fallback` if
// there is no such argument. Beware: args[i].f64() would quietly turn "abc"
// into 0.0; strconv.atof64 returns an error instead, though it takes "nan"
// and "inf".
fn number(args []string, i int, fallback f64) !f64 {
	if args.len <= i {
		return fallback
	}
	// Wrapping an error with context, like Go's fmt.Errorf("...: %w", err).
	x := strconv.atof64(args[i]) or { return error('argument ${i + 1} (${args[i]}): ${err}') }
	if !(x > 0.0) || !math.is_finite(x) {
		return error('argument ${i + 1} (${args[i]}): T and P must be above zero')
	}
	return x
}
