// M8 step 0 spike: where a PFR's solve allocates. It reads the reactors'
// own parts, so it's a test of that module, kept here: copy it there to run
// it, and remove it after.
//
//     cp spike/study/probe/reactors_alloc.v reactors/alloc_probe_test.v
//     v -cc gcc -prod -stats test reactors/alloc_probe_test.v
//     rm reactors/alloc_probe_test.v
//
// For each reforming flowsheet it counts what the collector hands out, by
// its own total, for each thing a step of the integration calls: the
// phase's state, the rates, their derivatives and the viscosity. Then for
// a whole solve, so that what the parts don't explain is the integrator's
// and the unit's own.
module reactors

import flowsheet
import os

fn handed_out() u64 {
	return u64(gc_heap_usage().total_bytes)
}

fn probe(file string) ! {
	l := flowsheet.load(os.join_path(@VMODROOT, 'testdata', 'flowsheets', file), makers: makers())!
	inlet := l.sheet.feeds[0].stream
	for p in l.sheet.units {
		u := p.unit
		if u is Pfr {
			n := u.k.species.len
			mut src := new_source(u.name, u.k, u.reactions, u.catalyst.loading, l.fluid.model)!
			mut s := []f64{len: n}
			mut dflows := []f64{len: n * n}
			mut dt := []f64{len: n}
			mut dp := []f64{len: n}
			calls := 1000
			mut st := src.state(inlet.t, inlet.p, inlet.flows)!
			mut b := handed_out()
			for _ in 0 .. calls {
				st = src.state(inlet.t, inlet.p, inlet.flows)!
			}
			state_bytes := f64(handed_out() - b) / calls
			b = handed_out()
			for _ in 0 .. calls {
				src.at(st, mut s)!
			}
			rate_bytes := f64(handed_out() - b) / calls
			b = handed_out()
			for _ in 0 .. calls {
				src.derivs(st, inlet.total(), mut dflows, mut dt, mut dp)
			}
			deriv_bytes := f64(handed_out() - b) / calls
			mut visc_bytes := 0.0
			if visc := u.viscosity {
				b = handed_out()
				for _ in 0 .. calls {
					visc.viscosity(inlet.t, inlet.p, src.y)!
				}
				visc_bytes = f64(handed_out() - b) / calls
			}
			u.calculate(l.fluid, [inlet])!
			b = handed_out()
			out := u.calculate(l.fluid, [inlet])!
			whole := f64(handed_out() - b)
			steps := (out.profile or { flowsheet.Profile{} }).points.len - 1
			println('${file}: a call allocates, in bytes: the state ${state_bytes:.0f}, the rates ${rate_bytes:.0f}, their derivatives ${deriv_bytes:.0f}, the viscosity ${visc_bytes:.0f}; the whole solve ${whole / 1e6:.2f} MB over ${steps} steps, ${whole / steps / 1e3:.1f} kB a step')
		}
	}
}

fn test_where_a_pfr_allocates() ! {
	probe('reformer_dwsim.toml')!
	probe('reformer.toml')!
}
