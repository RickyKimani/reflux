// The graph and its tears (exit criterion 1), against testdata/m4_graphs.json:
// the nine classical graphs of Gundersen and Hertzberg (1983), which its
// generator holds to the paper's tables, graphs worked by hand, and seeded
// random ones. Each plan is checked from scratch: every stream is calculated
// before it's used, once its tears are cut, and the tears are as few as the
// exact answer, with as few second breaks. Then a malformed flowsheet is
// refused, each for its own reason.
import flowsheet { Flowsheet, FlowsheetError, Outcome, Stream, Unit }
import json2
import os
import thermo

struct GraphCase {
	name          string
	vertices      int
	edges         [][]int
	nets          [][]int
	cycles        int
	cycle_list    [][]int
	min_tears     int
	second_breaks int
}

struct Fixture {
	graphs []GraphCase
}

// Stub is a unit that's only its ports: the graph never calculates.
struct Stub {
	name    string
	inlets  []string
	outlets []string
}

fn (s Stub) inlet_ports() []string {
	return s.inlets
}

fn (s Stub) outlet_ports() []string {
	return s.outlets
}

fn (s Stub) calculate(_ thermo.Fluid, _ []Stream) !Outcome {
	return error('a stub')
}

fn load() !Fixture {
	path := os.join_path(@VMODROOT, 'testdata', 'm4_graphs.json')
	return json2.decode[Fixture](os.read_file(path)!)!
}

// sheet builds a flowsheet with a stub per node, a stream si for edge i, and
// a feed to every node, so that no loop goes without one.
fn sheet(n int, edges [][]int) !Flowsheet {
	return sheet_in_order(n, edges, []int{len: edges.len, init: index})
}

// sheet_in_order lists each stub's ports in the order the edges are given.
fn sheet_in_order(n int, edges [][]int, order []int) !Flowsheet {
	mut f := Flowsheet{}
	for node in 1 .. n + 1 {
		mut ins := ['feed']
		mut outs := []string{}
		mut ports := {
			'feed': 'feed${node}'
		}
		for i in order {
			e := edges[i]
			if e[0] == node {
				port := 'out${i}'
				outs << port
				ports[port] = 's${i}'
			}
			if e[1] == node {
				port := 'in${i}'
				ins << port
				ports[port] = 's${i}'
			}
		}
		f.add_feed('feed${node}', Stream{})!
		f.add_unit(Stub{'u${node}', ins, outs}, ports)!
	}
	return f
}

// tears_of is the edge indices of all a plan's tears.
fn tears_of(p flowsheet.Plan) []int {
	mut out := []int{}
	for b in p.blocks {
		for name in b.tears {
			out << name[1..].int()
		}
	}
	out.sort()
	return out
}

// check_order checks a plan from scratch: every unit appears once, and each
// stream that isn't torn leaves its unit before it enters the next.
fn check_order(n int, edges [][]int, p flowsheet.Plan, at string) {
	mut when := []int{len: n + 1, init: -1}
	mut step := 0
	for b in p.blocks {
		for u in b.units {
			assert when[u + 1] == -1, '${at}: unit ${u} twice'
			when[u + 1] = step
			step++
		}
	}
	assert step == n, at
	torn := tears_of(p)
	for i, e in edges {
		if i !in torn {
			assert when[e[0]] < when[e[1]], '${at}: stream ${i} is used before it is calculated'
		}
	}
}

// acyclic is whether a graph has no cycle left, by repeatedly removing nodes
// with nothing coming in. It knows nothing of the code under test.
fn acyclic(n int, edges [][]int, cut []int) bool {
	mut gone := []bool{len: n + 1}
	for _ in 0 .. n + 1 {
		for v in 1 .. n + 1 {
			mut blocked := false
			for i, e in edges {
				if i !in cut && !gone[e[0]] && e[1] == v {
					blocked = true
				}
			}
			if !blocked {
				gone[v] = true
			}
		}
	}
	return gone[1..].all(it)
}

fn test_every_graph_matches_its_reference() {
	// Each graph's nets, cycles and tears against the exact answers, and each
	// plan against the checks above. Problem 1 has 415 cycles, and its 21
	// tears take the most search.
	f := load()!
	assert f.graphs.len == 59
	for g in f.graphs {
		s := sheet(g.vertices, g.edges)!
		p := s.plan()!
		assert p.notes.len == 0, '${g.name}: ${p.notes}'
		check_order(g.vertices, g.edges, p, g.name)

		// The nets are the blocks with cycles in them.
		mut nets := p.blocks.filter(it.cycles > 0).map(it.units.map(it + 1))
		for mut net in nets {
			net.sort()
		}
		nets.sort_with_compare(fn (a &[]int, b &[]int) int {
			return (*a)[0] - (*b)[0]
		})
		assert nets == g.nets, '${g.name}: ${nets}'
		mut cycles := 0
		for b in p.blocks {
			cycles += b.cycles
		}
		assert cycles == g.cycles, g.name

		// As few tears as there can be, with as few second breaks, and each one
		// needed.
		torn := tears_of(p)
		assert torn.len == g.min_tears, '${g.name}: ${torn}'
		assert acyclic(g.vertices, g.edges, torn), '${g.name}: a cycle is left'
		mut breaks := 0
		for c in g.cycle_list {
			hits := c.filter(it in torn).len
			assert hits > 0, "${g.name}: cycle ${c} isn't torn"
			breaks += hits - 1
		}
		assert breaks == g.second_breaks, '${g.name}: ${breaks} vs ${g.second_breaks}'
		for t in torn {
			assert !acyclic(g.vertices, g.edges, torn.filter(it != t)), "${g.name}: ${t} isn't needed"
		}
	}
}

fn test_a_plan_is_reproducible() {
	// The same flowsheet gives the same plan, whichever order its units'
	// ports are listed in.
	f := load()!
	for g in f.graphs {
		backwards := []int{len: g.edges.len, init: g.edges.len - 1 - index}
		a := sheet(g.vertices, g.edges)!.plan()!
		b := sheet_in_order(g.vertices, g.edges, backwards)!.plan()!
		assert a == b, g.name
	}
}

fn test_a_named_tear_is_honoured() {
	f := load()!
	p5 := f.graphs[4]
	assert p5.name.starts_with('problem 5')
	// Whichever loop stream is named, it's torn, and the plan still works:
	// more are added if the loops need them, none if they don't.
	for i, _ in p5.edges {
		s := sheet(p5.vertices, p5.edges)!
		in_a_loop := p5.cycle_list.any(i in it)
		p := s.plan(tears: ['s${i}']) or {
			assert !in_a_loop, '${i}: ${err}'
			assert err.msg().contains('in no loop'), err.msg()
			continue
		}
		assert in_a_loop
		torn := tears_of(p)
		assert i in torn, '${i}: ${torn}'
		check_order(p5.vertices, p5.edges, p, 'tear ${i}')
		assert torn.len >= p5.min_tears
	}
	// Tears that break everything are the whole set, though not the fewest.
	cycle := [[1, 2], [2, 3], [3, 1], [3, 2]]
	s := sheet(3, cycle)!
	assert tears_of(s.plan()!).len == 1
	both := s.plan(tears: ['s3', 's2'])!
	assert both.blocks.map(it.tears).filter(it.len > 0) == [['s2', 's3']]
	// Named twice is named once.
	assert tears_of(s.plan(tears: ['s2', 's3', 's2'])!) == [2, 3]
	// One of the two loops named: the other still gets a tear.
	two := sheet(4, [[1, 2], [2, 1], [2, 3], [3, 4], [4, 3]])!
	mut named := two.plan(tears: ['s1'])!
	assert tears_of(named).len == 2 && 1 in tears_of(named)
	assert named.blocks.len == 2 && named.blocks[0].tears == ['s1']
}

fn test_too_many_cycles_falls_back_to_a_greedy_cover() {
	// The complete graph on 9 nodes has far too many cycles to list. It takes
	// at least 45 tears: its 9 self-loops and one of each of 36 pairs. Greedy
	// isn't promised to reach that, but the plan must still be right and say
	// so.
	n := 9
	mut edges := [][]int{}
	for a in 1 .. n + 1 {
		for b in 1 .. n + 1 {
			edges << [a, b]
		}
	}
	p := sheet(n, edges)!.plan()!
	check_order(n, edges, p, 'K9')
	torn := tears_of(p)
	assert torn.len >= 45 && torn.len <= 50, '${torn.len}'
	assert acyclic(n, edges, torn)
	assert p.notes.len == 1 && p.notes[0].contains('greedy'), '${p.notes}'
	// The listing stopped at its limit: it has 125,673 to find.
	assert p.blocks[0].cycles == 50000
	// A named tear goes in too, once: it's the one greedy would take first.
	q := sheet(n, edges)!.plan(tears: ['s1'])!
	named := tears_of(q)
	assert named.filter(it == 1).len == 1, '${named}'
	check_order(n, edges, q, 'K9 with a tear')
}

fn test_ties_go_to_the_earliest_streams() {
	// a redundant pair to avoid: {s0, s3} and {s1, s2} are both 2 tears with no
	// second break, and {s0, s2} has one. The earliest of the first two wins.
	edges := [[1, 2], [2, 1], [2, 3], [3, 2], [3, 1]]
	assert tears_of(sheet(3, edges)!.plan()!) == [0, 3]
}

fn test_a_search_out_of_work_says_so() {
	// With no work to spend it keeps the first set it finds, which still
	// breaks every cycle, and the report says it may not be the fewest.
	f := load()!
	g := f.graphs[0]
	p := sheet(g.vertices, g.edges)!.plan(work: 0)!
	assert p.notes.len == 1 && p.notes[0].contains('work limit'), '${p.notes}'
	check_order(g.vertices, g.edges, p, g.name)
	assert acyclic(g.vertices, g.edges, tears_of(p))
}

// refused asserts that planning fails as a bad flowsheet, for this reason.
fn refused(f Flowsheet, unit string, want string, tears []string) {
	f.plan(tears: tears) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .bad_flowsheet && e.unit == unit, err.msg()
		assert err.msg().contains(want), err.msg()
		return
	}
	assert false, 'accepted: ${want}'
}

// add_refused asserts that adding a unit fails, for this reason.
fn add_refused(u Unit, ports map[string]string, want string) {
	mut f := Flowsheet{}
	f.add_unit(u, ports) or {
		assert err is FlowsheetError, err.msg()
		e := err as FlowsheetError
		assert e.kind == .bad_flowsheet && e.unit == u.name, err.msg()
		assert err.msg().contains(want), err.msg()
		return
	}
	assert false, 'accepted: ${want}'
}

fn stub(name string, ins []string, outs []string) Stub {
	return Stub{name, ins, outs}
}

fn test_a_malformed_flowsheet_is_refused_for_its_reason() {
	// A port left open, in or out.
	u := stub('u', ['in'], ['out'])
	add_refused(u, {
		'out': 's'
	}, 'port in is open')
	add_refused(u, {
		'in': 's'
	}, 'port out is open')
	add_refused(u, {
		'in':  ''
		'out': 's'
	}, 'port in is open')
	// A port the unit doesn't have: the wrong number of ports for a mixer.
	mixer := flowsheet.new_mixer('M', 2)!
	add_refused(mixer, {
		'in1': 'a'
		'in2': 'b'
		'in3': 'c'
		'out': 'd'
	}, 'no port in3, only in1, in2, out')
	add_refused(stub('', [], []), {}, 'a unit has no name')

	// Two units of one name, and two feeds.
	mut twice := Flowsheet{}
	twice.add_unit(u, {
		'in':  'a'
		'out': 'b'
	})!
	if _ := twice.add_unit(u, {
		'in':  'c'
		'out': 'd'
	})
	{
		assert false, 'two units of one name'
	} else {
		assert err.msg().contains('two units'), err.msg()
	}
	twice.add_feed('a', Stream{})!
	if _ := twice.add_feed('a', Stream{}) {
		assert false, 'two feeds of one name'
	} else {
		assert err.msg().contains('two feeds'), err.msg()
	}
	if _ := twice.add_feed('', Stream{}) {
		assert false, 'a feed with no name'
	} else {
		assert err.msg().contains('no name'), err.msg()
	}

	// A flowsheet with nothing in it.
	refused(Flowsheet{}, '', 'needs a unit', [])

	// A stream with two sources: two outlets, or a feed and an outlet.
	mut a := Flowsheet{}
	a.add_feed('f', Stream{})!
	a.add_unit(stub('p', ['in'], ['out']), {
		'in':  'f'
		'out': 'x'
	})!
	a.add_unit(stub('q', [], ['out']), {
		'out': 'x'
	})!
	refused(a, 'q', 'stream x has two sources, this and p', [])
	mut b := Flowsheet{}
	b.add_feed('f', Stream{})!
	b.add_unit(stub('p', [], ['out']), {
		'out': 'f'
	})!
	refused(b, 'p', 'stream f has two sources, this and a feed', [])

	// A stream to two inlets.
	mut c := Flowsheet{}
	c.add_feed('f', Stream{})!
	c.add_unit(stub('p', ['in'], []), {
		'in': 'f'
	})!
	c.add_unit(stub('q', ['in'], []), {
		'in': 'f'
	})!
	refused(c, 'q', 'stream f goes to two inlets, here and p', [])

	// A stream from nowhere: a misspelt name, usually.
	mut d := Flowsheet{}
	d.add_unit(stub('p', ['in'], []), {
		'in': 'nothing'
	})!
	refused(d, 'p', 'stream nothing comes from no feed or outlet', [])

	// A feed that goes nowhere.
	mut e := Flowsheet{}
	e.add_feed('lost', Stream{})!
	e.add_unit(stub('p', [], ['out']), {
		'out': 'x'
	})!
	refused(e, '', 'feed lost goes to no inlet', [])

	// A loop fed by a unit is fed, and a unit with no inlets, which makes what
	// it sends out, needs nothing to feed it.
	mut fed := Flowsheet{}
	fed.add_unit(stub('source', [], ['out']), {
		'out': 'a'
	})!
	fed.add_unit(stub('p', ['a', 'back'], ['out']), {
		'a':    'a'
		'back': 'y'
		'out':  'x'
	})!
	fed.add_unit(stub('q', ['in'], ['out']), {
		'in':  'x'
		'out': 'y'
	})!
	assert fed.plan()!.blocks.len == 2

	// A loop with nothing feeding it.
	mut g := Flowsheet{}
	g.add_unit(stub('p', ['in'], ['out']), {
		'in':  'y'
		'out': 'x'
	})!
	g.add_unit(stub('q', ['in'], ['out']), {
		'in':  'x'
		'out': 'y'
	})!
	refused(g, '', 'nothing feeds the loop of p, q', [])
	mut h := Flowsheet{}
	h.add_unit(stub('p', ['in'], ['out']), {
		'in':  'x'
		'out': 'x'
	})!
	refused(h, '', 'nothing feeds the loop of p', [])

	// A tear that isn't a stream between units, or isn't in a loop.
	mut k := sheet(3, [[1, 2], [2, 3], [3, 2]])!
	refused(k, '', 'tear nothing is not a stream between two units', ['nothing'])
	refused(k, '', 'tear feed1 is not a stream between two units', ['feed1'])
	refused(k, '', 'tear s0 is in no loop', ['s0'])
}

fn test_a_product_is_a_stream_with_no_inlet_after_it() {
	// A real flowsheet's ports: the reference flowsheet's shape, with stream
	// names, and blocks in order.
	mixer := flowsheet.new_mixer('mixer', 2)!
	cooler := flowsheet.new_heater('cooler', mode: .outlet_t, value: 230.0)!
	drum := flowsheet.new_flash_drum('drum', mode: .duty, value: 0.0)!
	split := flowsheet.new_splitter('split', [0.5, 0.5])!
	pump := flowsheet.new_pump('pump', mode: .delta_p, value: 1e5, efficiency: 0.75)!
	valve := flowsheet.new_valve('valve', mode: .outlet_p, value: 5e5)!
	mut f := Flowsheet{}
	f.add_feed('feed', Stream{})!
	// Added out of order: a stream can be named before its other end exists.
	f.add_unit(valve, {
		'in':  'on'
		'out': 'product'
	})!
	f.add_unit(drum, {
		'in1':    'cold'
		'vapour': 'gas'
		'liquid': 'liquid'
	})!
	f.add_unit(pump, {
		'in':  'back'
		'out': 'recycle'
	})!
	f.add_unit(mixer, {
		'in1': 'feed'
		'in2': 'recycle'
		'out': 'mixed'
	})!
	f.add_unit(cooler, {
		'in':  'mixed'
		'out': 'cold'
	})!
	f.add_unit(split, {
		'in':   'liquid'
		'out1': 'back'
		'out2': 'on'
	})!
	p := f.plan()!
	assert p.notes.len == 0
	names := p.blocks.map(it.units.map(f.units[it].unit.name))
	assert names.len == 2
	assert names[1] == ['valve']
	assert names[0].len == 5
	assert p.blocks[0].cycles == 1 && p.blocks[0].tears.len == 1
	assert p.blocks[1].cycles == 0 && p.blocks[1].tears.len == 0
	// Any stream of a single loop would do, so the choice falls to the one
	// coming back to where the feed enters, though it's far from the earliest.
	assert p.blocks[0].tears == ['recycle']
	assert names[0] == ['mixer', 'cooler', 'drum', 'split', 'pump']
	// A named tear still wins.
	q := f.plan(tears: ['liquid'])!
	assert q.blocks[0].tears == ['liquid']
	assert q.blocks[0].units.map(f.units[it].unit.name) == ['split', 'pump', 'mixer', 'cooler',
		'drum']
}
