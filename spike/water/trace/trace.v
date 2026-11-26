// M5a step 2's trace: every cubic answer, bit for bit, before and after D7's
// changes to thermo's volume solvers and saturation_pressure.
//
//     v run spike/water/trace > before.txt     (on the code before)
//     v run spike/water/trace > after.txt      (on the code after)
//     cmp before.txt after.txt
//
// It prints every float of every answer as its bits, and every refusal by
// its message:
// - the sweep's seeded mixtures (M3), PR: the stability test and the PT
//   flash at each point; at one point in five, SRK and a PR with Péneloux's
//   shift too, the three volume roots and their properties, the
//   vapour-fraction flashes at β 0, 0.4 and 1 by T and by P, PH and PS back
//   to the PT answer, and TV and UV from its state;
// - the first 160 pure compounds of the sweep's pool: saturation_pressure
//   and the pure vapour-fraction flash at five temperatures up to 0.9999 Tc,
//   and PT flashes either side of the saturation line.
//
// It's run from the repository root. Not part of the simulator.
module main

import components
import flash
import math
import os
import sweep
import thermo

struct Out {
mut:
	lines []string
}

fn b(x f64) string {
	return math.f64_bits(x).hex()
}

fn bs(xs []f64) string {
	return xs.map(b(it)).join(',')
}

fn (mut o Out) say(s string) {
	o.lines << s
}

fn (mut o Out) state(tag string, st thermo.State, fl thermo.Fluid) {
	pr := fl.props(st)
	o.say('${tag} v=${b(st.v)} z=${b(st.z)} lnphi=${bs(st.ln_phi)} h=${b(pr.h)} s=${b(pr.s)} cp=${b(pr.cp)} cv=${b(pr.cv)} w=${b(pr.w)} jt=${b(pr.mu_jt)} label=${thermo.phase_label(fl.model, st)}')
}

fn (mut o Out) result(tag string, r flash.FlashResult, fl thermo.Fluid) {
	o.say('${tag} t=${b(r.t)} p=${b(r.p)} beta=${b(r.beta)} path=${r.path} it=${r.iterations} stab=${r.stability.stable} trials=${r.stability.trials.len}')
	for i, ph in r.phases {
		o.say('${tag} phase${i} ${ph.label} frac=${b(ph.fraction)} x=${bs(ph.state.n)}')
		o.state('${tag} phase${i}', ph.state, fl)
	}
}

fn (mut o Out) fail(tag string, err IError) {
	o.say('${tag} ERR ${err.msg()}')
}

fn (mut o Out) flash_all(tag string, fl thermo.Fluid, pt sweep.Point, full bool) {
	m := fl.model
	st := flash.stability(m, pt.t, pt.p, pt.z) or {
		o.fail('${tag} stab', err)
		flash.Stability{}
	}
	o.say('${tag} stab stable=${st.stable} trials=${st.trials.len} tm=${bs(st.trials.map(it.tm))}')
	r := flash.pt(m, pt.t, pt.p, pt.z) or {
		o.fail('${tag} pt', err)
		return
	}
	o.result('${tag} pt', r, fl)
	if !full {
		return
	}
	for root in [thermo.Root.liquid, .vapor, .stable] {
		s := thermo.solve_tp(m, pt.t, pt.p, pt.z, root) or {
			o.fail('${tag} root ${root}', err)
			continue
		}
		o.state('${tag} root ${root}', s, fl)
	}
	for beta in [0.0, 0.4, 1.0] {
		if a := flash.tvf(m, pt.t, beta, pt.z) {
			o.result('${tag} tvf ${beta}', a, fl)
		} else {
			o.fail('${tag} tvf ${beta}', err)
		}
		if a := flash.pvf(m, pt.p, beta, pt.z) {
			o.result('${tag} pvf ${beta}', a, fl)
		} else {
			o.fail('${tag} pvf ${beta}', err)
		}
	}
	mut h := 0.0
	mut s := 0.0
	mut u := 0.0
	mut v := 0.0
	for ph in r.phases {
		pr := fl.props(ph.state)
		h += ph.fraction * pr.h
		s += ph.fraction * pr.s
		u += ph.fraction * pr.u
		v += ph.fraction * ph.state.v
	}
	if a := flash.ph(fl, pt.p, h, pt.z) {
		o.result('${tag} ph', a, fl)
	} else {
		o.fail('${tag} ph', err)
	}
	if a := flash.ps(fl, pt.p, s, pt.z) {
		o.result('${tag} ps', a, fl)
	} else {
		o.fail('${tag} ps', err)
	}
	if a := flash.tv(m, pt.t, v, pt.z) {
		o.result('${tag} tv', a, fl)
	} else {
		o.fail('${tag} tv', err)
	}
	if a := flash.uv(fl, u, v, pt.z) {
		o.result('${tag} uv', a, fl)
	} else {
		o.fail('${tag} uv', err)
	}
}

fn main() {
	count := if os.args.len > 1 { os.args[1].int() } else { 3000 }
	db := components.load_default() or { panic(err) }
	pool := sweep.candidates(db) or { panic(err) }
	mut o := Out{}
	for pt in sweep.points(pool, count) {
		full := pt.index % 5 == 0
		pr := db.fluid(.pr, pt.cas) or {
			o.fail('${pt.index} fluid', err)
			continue
		}
		o.flash_all('${pt.index} pr', pr, pt, full)
		if full {
			srk := db.fluid(.srk, pt.cas) or { panic(err) }
			o.flash_all('${pt.index} srk', srk, pt, true)
			shift := pr.model.components().map(0.1 * it.tc * thermo.r_gas / it.pc * 0.05)
			tr := thermo.translate(pr.model, shift) or { panic(err) }
			trf := thermo.new_fluid(tr, pr.ideal) or { panic(err) }
			o.flash_all('${pt.index} shifted', trf, pt, true)
		}
		if o.lines.len > 20000 {
			println(o.lines.join('\n'))
			o.lines.clear()
		}
	}
	for k, e in pool[..160] {
		fl := db.fluid(.pr, [e.cas]) or {
			o.fail('pure ${k} fluid', err)
			continue
		}
		c := fl.model.components()[0]
		for f in [0.5, 0.7, 0.9, 0.99, 0.9999] {
			t := f * c.tc
			tag := 'pure ${k} ${f}'
			ps := thermo.saturation_pressure(fl.model, t) or {
				o.fail('${tag} psat', err)
				0.0
			}
			o.say('${tag} psat ${b(ps)}')
			if a := flash.tvf(fl.model, t, 0.0, [1.0]) {
				o.result('${tag} tvf', a, fl)
			} else {
				o.fail('${tag} tvf', err)
			}
			if ps > 0.0 {
				for g in [0.9, 1.1] {
					if a := flash.pt(fl.model, t, g * ps, [1.0]) {
						o.result('${tag} pt ${g}', a, fl)
					} else {
						o.fail('${tag} pt ${g}', err)
					}
				}
			}
		}
	}
	println(o.lines.join('\n'))
}
