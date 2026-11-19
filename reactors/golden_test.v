// The golden flowsheets with reactors (M9 step 10): the published reformer
// on Xu and Froment's kinetics, and the owner's Tutorial 2 as its tutorial
// prints its conversion, whose report warns of an outlet past equilibrium.
// Each report in testdata/flowsheets/golden/reactors/ is tools/golden.vsh's,
// and solving it again gives the same to 1e-9 of each value's scale, as
// flowsheet/golden_test.v holds M4's: T and P of themselves, a flow of the
// feed, H of RT, S of R and a duty of F·RT; and the same warnings, and the
// same profile, its points' volumes to 1e-6 (below).
import flowsheet { Report }
import json2
import math
import os
import reactors
import thermo

const golden_dir = os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'golden', 'reactors')

fn close(got f64, want f64, scale f64, at string) {
	assert math.abs(got - want) <= 1e-9 * scale, '${at}: ${got} against ${want}'
}

fn test_every_golden_flowsheet_with_reactors_solves_to_its_results() {
	files := os.ls(golden_dir)!.filter(it.ends_with('.json'))
	assert files.len == 2, '${files}'
	for f in files {
		want := json2.decode[Report](os.read_file(os.join_path(golden_dir, f))!)!
		assert want.file == 'testdata/flowsheets/${f.all_before('.json')}.toml'
		assert want.reflux == flowsheet.version(), '${f} was written by ${want.reflux}: run tools/golden.vsh for ${flowsheet.version()}'
		l := flowsheet.load(os.join_path(@VMODROOT, want.file), makers: reactors.makers())!
		got := l.report(l.solve()!)
		assert got.streams.map(it.name) == want.streams.map(it.name), f
		assert got.units.map(it.name) == want.units.map(it.name), f
		mut feed := 0.0
		for w in want.streams {
			if w.from == none {
				feed += w.molar_flow
			}
		}
		rt := thermo.r_gas * got.streams[0].t
		for k, w in want.streams {
			g := got.streams[k]
			at := '${f} ${w.name}'
			close(g.t, w.t, w.t, '${at} T')
			close(g.p, w.p, w.p, '${at} P')
			close(g.h, w.h, rt, '${at} H')
			close(g.s, w.s, thermo.r_gas, '${at} S')
			for i, x in w.molar_flows {
				close(g.molar_flows[i], x, feed, '${at} flow ${i}')
			}
		}
		for k, w in want.units {
			g := got.units[k]
			assert g.warnings == w.warnings, '${f} ${w.name}: ${g.warnings}'
			for i, e in w.energy {
				close(g.energy[i].rate, e.rate, feed * rt, '${f} ${w.name} ${e.name}')
			}
			wp := w.profile or { []flowsheet.ProfileReport{} }
			gp := g.profile or { []flowsheet.ProfileReport{} }
			assert gp.len == wp.len, '${f} ${w.name}: ${gp.len} profile points, not ${wp.len}'
			for i, pt in wp {
				// A profile's points are where the integrator stepped, and its
				// step control is moved by rounding: gcc puts the reformer's
				// within 1.1e-9 of tcc's, which wrote the report, and their P
				// within 1.0e-9 with them. So each point's volume is held to
				// 1e-6 of the reactor's, and its T and P to 1e-6 of themselves;
				// the outlet, in the streams, is held to 1e-9 above.
				assert math.abs(gp[i].volume - pt.volume) <= 1e-6 * wp.last().volume, '${f} ${w.name} V ${i}'
				assert math.abs(gp[i].t - pt.t) <= 1e-6 * pt.t, '${f} ${w.name} T ${i}'
				assert math.abs(gp[i].p - pt.p) <= 1e-6 * pt.p, '${f} ${w.name} P ${i}'
			}
		}
	}
}
