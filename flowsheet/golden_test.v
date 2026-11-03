// The golden flowsheets (D13): whole flowsheets whose results can't drift
// unnoticed. Each file in testdata/flowsheets/golden/ is the report that
// tools/golden.vsh wrote for one, and solving it again gives the same, to 1e-9
// of each value's scale: T and P of themselves, a flow of the feed, H of RT and
// S of R, a duty of F·RT, and the same passes. They're regenerated only on
// purpose, with the change explained.
import flowsheet { Report }
import json2
import math
import os
import thermo

const golden_dir = os.join_path(@VMODROOT, 'testdata', 'flowsheets', 'golden')

fn close(got f64, want f64, scale f64, at string) {
	assert math.abs(got - want) <= 1e-9 * scale, '${at}: ${got} against ${want}'
}

fn test_every_golden_flowsheet_solves_to_its_results() {
	files := os.ls(golden_dir)!.filter(it.ends_with('.json'))
	assert files.len == 5, '${files}'
	for f in files {
		want := json2.decode[Report](os.read_file(os.join_path(golden_dir, f))!)!
		assert want.file == 'testdata/flowsheets/${f.all_before('.json')}.toml'
		l := flowsheet.load(os.join_path(@VMODROOT, want.file))!
		same(l.report(l.solve()!), want, f)
	}
}

fn same(got Report, want Report, at string) {
	assert got.status == 'solved' && want.status == 'solved', at
	assert got.streams.map(it.name) == want.streams.map(it.name), at
	assert got.units.map(it.name) == want.units.map(it.name), at
	mut feed := 0.0
	mut streams := map[string]flowsheet.StreamReport{}
	for w in want.streams {
		streams[w.name] = w
		if w.from == none {
			feed += w.molar_flow
		}
	}
	for k, w in want.streams {
		g := got.streams[k]
		where := '${at} ${w.name}'
		assert g.from == w.from && g.to == w.to, where
		close(g.t, w.t, w.t, '${where} T')
		close(g.p, w.p, w.p, '${where} P')
		close(g.vapour_fraction, w.vapour_fraction, 1.0, '${where} beta')
		close(g.h, w.h, thermo.r_gas * w.t, '${where} H')
		close(g.s, w.s, thermo.r_gas, '${where} S')
		close(g.v, w.v, w.v, '${where} V')
		for i, n in w.molar_flows {
			close(g.molar_flows[i], n, feed, '${where} flow ${i}')
			close(g.mole_fractions[i], w.mole_fractions[i], 1.0, '${where} x ${i}')
		}
	}
	for k, w in want.units {
		g := got.units[k]
		assert g.kind == w.kind && g.inlets == w.inlets && g.outlets == w.outlets, '${at} ${w.name}'
		mut flow, mut t := 0.0, 0.0
		for name in w.inlets {
			flow += streams[name].molar_flow
			t = math.max(t, streams[name].t)
		}
		assert g.energy.len == w.energy.len, '${at} ${w.name}'
		for i, e in w.energy {
			assert g.energy[i].name == e.name
			close(g.energy[i].rate, e.rate, math.max(flow, feed) * thermo.r_gas * t, '${at} ${w.name} ${e.name}')
		}
	}
	assert got.blocks.len == want.blocks.len, at
	for k, w in want.blocks {
		g := got.blocks[k]
		assert g.units == w.units && g.tears == w.tears && g.cycles == w.cycles, at
		// The same passes, none more.
		assert g.residuals.len == w.residuals.len, '${at}: ${g.residuals.len} passes, not ${w.residuals.len}'
	}
}
