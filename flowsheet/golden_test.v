// The golden flowsheets (D13): whole flowsheets whose results can't drift
// unnoticed. Each file in testdata/flowsheets/golden/ is the report that
// tools/golden.vsh wrote for one, and solving it again gives the same, to 1e-9
// of each value's scale: T and P of themselves, a flow of the feed, H of RT and
// S of R, a duty of F·RT, and the same passes. They're regenerated only on
// purpose, with the change explained.
import flowsheet { Report, Stream }
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
		// A version's golden reports are its own, written again when v.mod's
		// version changes, so the release gate sees each and its entry names
		// it (M10a step 9).
		assert want.reflux == flowsheet.version(), '${f} was written by ${want.reflux}: run tools/golden.vsh for ${flowsheet.version()}'
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

fn test_every_golden_flowsheet_balances_from_scratch() {
	// Exit criterion 3, which solve checks for itself: here each unit's
	// component and energy balances are worked out again on the inlets it used,
	// with H from T, P and each phase's composition, to 1e-10 of its flow and of
	// F·RT; and the whole, feeds against products, to what the tears leave open.
	// solve_test.v does the same for its flowsheets and report_test.v from a
	// report's numbers, but no_recycle and saturated_loop are checked only here.
	for f in os.ls(golden_dir)!.filter(it.ends_with('.json')) {
		l := flowsheet.load(os.join_path(@VMODROOT, 'testdata', 'flowsheets', f.replace('.json',
			'.toml')))!
		sol := l.solve()!
		n := l.names.len
		for r in sol.units {
			mut c := []f64{len: n}
			mut q, mut flow, mut hot := 0.0, 0.0, 0.0
			for s in r.inlets {
				q += tally(mut c, l.fluid, s, 1.0)!
				flow += s.total()
				hot = math.max(hot, s.t)
			}
			for s in r.outlets {
				q += tally(mut c, l.fluid, s, -1.0)!
			}
			for e in r.energy {
				q += e.rate
			}
			if flow == 0.0 {
				continue
			}
			for i, x in c {
				assert math.abs(x) <= 1e-10 * flow, '${f} ${r.name}: component ${i} is off by ${x}'
			}
			assert math.abs(q) <= 1e-10 * flow * thermo.r_gas * hot, '${f} ${r.name}: energy is off by ${q}'
		}
		mut c := []f64{len: n}
		mut q, mut feed, mut feed_t := 0.0, 0.0, 0.0
		for x in l.sheet.feeds {
			q += tally(mut c, l.fluid, x.stream, 1.0)!
			feed += x.stream.total()
			feed_t = math.max(feed_t, x.stream.t)
		}
		for name in sol.products {
			q += tally(mut c, l.fluid, sol.streams[name], -1.0)!
		}
		for r in sol.units {
			for e in r.energy {
				q += e.rate
			}
		}
		// What the last pass can leave open in a tear: tol of its flow in each
		// component, which carries the stream's H, and tol of RT in its H.
		mut open_flow, mut open_energy := 0.0, 0.0
		for b in sol.plan.blocks {
			for name in b.tears {
				s := sol.streams[name]
				open_flow += l.params.tol * s.total()
				open_energy += l.params.tol * s.total() * (f64(n) * math.abs(s.h) +
					thermo.r_gas * s.t)
			}
		}
		for i, x in c {
			assert math.abs(x) <= open_flow + 1e-12 * feed, '${f}: component ${i} overall is off by ${x}'
		}
		assert math.abs(q) <= open_energy + 1e-10 * feed * thermo.r_gas * feed_t, '${f}: energy overall is off by ${q}'
	}
}

// tally adds a stream's flows to c with its sign, and returns its enthalpy
// flow, F·H, with H from scratch: each phase's volume root solved again at the
// stream's T and P.
fn tally(mut c []f64, fl thermo.Fluid, s Stream, sign f64) !f64 {
	for i, x in s.flows {
		c[i] += sign * x
	}
	mut h := 0.0
	for phase in s.result.phases {
		if phase.fraction == 0.0 {
			continue
		}
		root := if phase.label == .vapor { thermo.Root.vapor } else { thermo.Root.liquid }
		h += phase.fraction * fl.props(thermo.solve_tp(fl.model, s.t, s.p, phase.state.n, root)!).h
	}
	return sign * s.total() * h
}
