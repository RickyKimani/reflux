// A small demo of the thermo layer.
//
//     v run .                 # n-butane at 350 K and 9.4573 bar (SVA Example 3.9)
//     v run . 300 2.5         # T in K, P in bar
module main

import math
import os
import strconv
import thermo { Component, CubicKind }

// Smith, Van Ness & Abbott, Table B.1.
const butane = Component{
	name:  'n-butane'
	tc:    425.1
	pc:    37.96e5
	omega: 0.200
}

const methane = Component{
	name:  'methane'
	tc:    190.6
	pc:    45.99e5
	omega: 0.012
}

// The PRSV example in Caleb Bell's `thermo` library. κ₁ is fitted together
// with a particular ω, so the two have to come from the same source.
const hexane = Component{
	name:  'n-hexane'
	tc:    507.6
	pc:    30.25e5
	omega: 0.2975
}

fn main() {
	t := arg_or(1, 350.0) or { die(err) }
	p := (arg_or(2, 9.4573) or { die(err) }) * 1e5
	pure_table(t, p) or { die(err) }
	println('')
	prsv_example() or { die(err) }
	println('')
	mixture_example() or { die(err) }
}

fn pure_table(t f64, p f64) ! {
	println('${butane.name} at T = ${t} K, P = ${p / 1e5} bar\n')
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
}

// prsv_example swaps PR's α for PRSV's on one component, which takes one
// option. At Tr ≈ 0.67 the two barely differ; PRSV earns its keep at low
// Tr and for polar compounds.
fn prsv_example() ! {
	tb := 341.88 // normal boiling point of n-hexane, K
	pr := thermo.new_cubic(.pr, [hexane])!
	prsv := thermo.new_cubic(.pr, [hexane], alphas: [thermo.prsv(hexane.omega, 0.05104)])!
	p_pr := thermo.saturation_pressure(pr, tb)!
	p_prsv := thermo.saturation_pressure(prsv, tb)!
	println('${hexane.name} at its normal boiling point, ${tb} K (measured Psat: 1.01325 bar)')
	println('  PR    ${p_pr / 1e5:.4f} bar')
	println('  PRSV  ${p_prsv / 1e5:.4f} bar')
}

fn mixture_example() ! {
	t, p := 300.0, 10e5
	n := [0.3, 0.7]
	eos := thermo.new_cubic(.pr, [methane, butane], kij: [0.0, 0.02, 0.02, 0.0])!
	println('30% methane + 70% n-butane, PR with kij = 0.02, at ${t} K and ${p / 1e5} bar')
	for root in [thermo.Root.liquid, .vapor] {
		s := thermo.solve_tp(eos, t, p, n, root)!
		phi := s.ln_phi.map(math.exp(it))
		println('  ${root:-7} Z = ${s.z:.5f}  phi(methane) = ${phi[0]:.4f}  phi(n-butane) = ${phi[1]:.4f}')
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
