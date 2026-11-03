module flowsheet

// Placed is a unit in a flowsheet, with the stream on each of its ports.
pub struct Placed {
pub:
	unit    Unit
	inlets  []string // stream names, in the order of inlet_ports
	outlets []string // stream names, in the order of outlet_ports
}

// Feed is a stream that enters the flowsheet from outside.
pub struct Feed {
pub:
	name   string
	stream Stream
}

// Flowsheet is units joined by named streams (D5). A stream joins one outlet
// to one inlet, a feed has no outlet behind it, and a product has no inlet
// after it. Units are added in any order, so a stream can be named before
// the unit at its other end exists. That's why the whole is only checked by
// plan.
pub struct Flowsheet {
pub mut:
	units []Placed
	feeds []Feed
}

// add_feed adds a stream from outside.
pub fn (mut f Flowsheet) add_feed(name string, s Stream) ! {
	if name == '' {
		return bad_sheet('', 'a feed has no name')
	}
	if f.feeds.any(it.name == name) {
		return bad_sheet('', 'there are two feeds named ${name}')
	}
	f.feeds << Feed{name, s}
}

// add_unit adds a unit, with ports mapping each of its port names to a
// stream's name. Every port needs a stream, and no others may be given.
pub fn (mut f Flowsheet) add_unit(u Unit, ports map[string]string) ! {
	if u.name == '' {
		return bad_sheet('', 'a unit has no name')
	}
	if f.units.any(it.unit.name == u.name) {
		return bad_sheet(u.name, 'there are two units of this name')
	}
	inlets := u.inlet_ports()
	outlets := u.outlet_ports()
	for port, _ in ports {
		if port !in inlets && port !in outlets {
			mut all := inlets.clone()
			all << outlets
			return bad_sheet(u.name, 'it has no port ${port}, only ${all.join(', ')}')
		}
	}
	f.units << Placed{
		unit:    u
		inlets:  streams_on(u, inlets, ports)!
		outlets: streams_on(u, outlets, ports)!
	}
}

// streams_on is the stream on each of these ports.
fn streams_on(u Unit, names []string, ports map[string]string) ![]string {
	mut out := []string{cap: names.len}
	for port in names {
		name := ports[port] or { '' }
		if name == '' {
			return bad_sheet(u.name, 'port ${port} is open')
		}
		out << name
	}
	return out
}

// graph checks that the streams join up, and returns the units as a graph.
fn (f Flowsheet) graph() !Graph {
	if f.units.len == 0 {
		return bad_sheet('', 'a flowsheet needs a unit')
	}
	mut source := map[string]int{}
	for feed in f.feeds {
		source[feed.name] = -1
	}
	for i, p in f.units {
		for name in p.outlets {
			if name in source {
				from := source[name]
				other := if from < 0 { 'a feed' } else { f.units[from].unit.name }
				return bad_sheet(p.unit.name, 'stream ${name} has two sources, this and ${other}')
			}
			source[name] = i
		}
	}
	mut target := map[string]int{}
	mut edges := []Edge{}
	mut fed := []bool{len: f.units.len}
	for i, p in f.units {
		for name in p.inlets {
			if name in target {
				other := f.units[target[name]].unit.name
				return bad_sheet(p.unit.name, 'stream ${name} goes to two inlets, here and ${other}')
			}
			target[name] = i
			from := source[name] or {
				return bad_sheet(p.unit.name, 'stream ${name} comes from no feed or outlet')
			}
			if from < 0 {
				fed[i] = true
			} else {
				edges << Edge{from, i, name}
			}
		}
	}
	for feed in f.feeds {
		if feed.name !in target {
			return bad_sheet('', 'feed ${feed.name} goes to no inlet')
		}
	}
	// In the order the streams leave their units, so that tears don't depend
	// on the order of the inlets.
	edges.sort_with_compare(fn (a &Edge, b &Edge) int {
		return if a.from != b.from { a.from - b.from } else { compare_strings(a.name, b.name) }
	})
	return new_graph(f.units.len, edges, fed)
}

// Block is a part of a flowsheet to solve in one go: one unit, or a net of
// units joined in loops, whose tears must be converged together.
pub struct Block {
pub:
	units  []int    // indices into Flowsheet.units, in calculation order
	tears  []string // the streams torn, whose first values the block starts from
	cycles int      // elementary cycles in the block, listed up to 50,000
}

// Plan is how to solve a flowsheet: blocks in order, each needing only
// feeds, the blocks before it, and its own tears.
pub struct Plan {
pub:
	blocks []Block
	notes  []string // anything the report should say about the choices
}

// PlanParams are plan's options.
@[params]
pub struct PlanParams {
pub:
	tears []string // streams to tear, whatever the search would choose
	// work is how long to search for the best tears, roughly in cycles looked
	// at. It's a count and not a time, so that a plan is reproducible. At this
	// default, a graph like the classical test problems takes a moment.
	work i64 = 200_000_000
}

// plan checks a flowsheet and decides the order to solve it in (D5).
//
// Tarjan's strongly connected components give the blocks, in an order where
// each needs only the ones before it. Within a block of more than one unit,
// Johnson's algorithm lists the elementary cycles, and tear_set picks the
// fewest streams that break them all. A tear in params is honoured, and
// more are chosen if the loops need them.
pub fn (f Flowsheet) plan(params PlanParams) !Plan {
	g := f.graph()!
	named := f.named_tears(g, params.tears)!
	mut blocks := []Block{}
	mut notes := []string{}
	for group in g.components() {
		inside := g.mask(group)
		if g.is_net(group) && !g.has_input(group, inside) {
			names := group.map(f.units[it].unit.name)
			return bad_sheet('', 'nothing feeds the loop of ${names.join(', ')}')
		}
		chosen := g.tear_set(group, named.filter(inside[g.edges[it].from]), params.work)
		notes << chosen.notes
		blocks << Block{
			units:  g.order(group, chosen.tears)
			tears:  chosen.tears.map(g.edges[it].name)
			cycles: chosen.cycles
		}
	}
	return Plan{
		blocks: blocks
		notes:  notes
	}
}

// named_tears finds the edges of the tears asked for. A tear has to be in a
// loop: tearing anything else would do nothing but cost an iteration.
fn (f Flowsheet) named_tears(g Graph, names []string) ![]int {
	mut out := []int{cap: names.len}
	mut group_of := []int{len: g.n}
	for i, group in g.components() {
		for u in group {
			group_of[u] = i
		}
	}
	for name in names {
		e := g.edges.filter(it.name == name)
		if e.len == 0 {
			return bad_sheet('', 'tear ${name} is not a stream between two units')
		}
		if group_of[e[0].from] != group_of[e[0].to] {
			return bad_sheet('', 'tear ${name} is in no loop')
		}
		i := g.edges.index(e[0])
		if i !in out {
			out << i
		}
	}
	return out
}
