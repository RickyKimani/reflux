module flowsheet

import math

// Edge is a stream from one unit to another.
struct Edge {
	from int
	to   int
	name string
}

// Graph is a flowsheet's units as nodes and the streams between them as
// edges. Two units can be joined by parallel edges, which are separate
// cycles, and an edge can leave and enter one unit. Everything here is by
// index, and in a fixed order, so that a plan is reproducible.
struct Graph {
	n     int
	edges []Edge
	out   [][]int // the edges leaving each node
	fed   []bool  // whether each node takes a feed
}

fn new_graph(n int, edges []Edge, fed []bool) Graph {
	mut out := [][]int{len: n}
	for i, e in edges {
		out[e.from] << i
	}
	return Graph{
		n:     n
		edges: edges
		out:   out
		fed:   fed
	}
}

// mask marks the nodes of a group.
fn (g Graph) mask(group []int) []bool {
	mut inside := []bool{len: g.n}
	for u in group {
		inside[u] = true
	}
	return inside
}

// is_net is whether a strongly connected group has a cycle: more than one
// node, or a stream from its one node to itself.
fn (g Graph) is_net(group []int) bool {
	return group.len > 1 || g.out[group[0]].any(g.edges[it].to == group[0])
}

// entries marks the streams that enter a unit where the group takes something
// from outside: a recycle, as it were, coming back to its mixer.
fn (g Graph) entries(group []int) []bool {
	inside := g.mask(group)
	mut open := []bool{len: g.n}
	for u in group {
		open[u] = g.fed[u]
	}
	for e in g.edges {
		if inside[e.to] && !inside[e.from] {
			open[e.to] = true
		}
	}
	return g.edges.map(inside[it.from] && inside[it.to] && open[it.to])
}

// has_input is whether anything comes into a group from outside it.
fn (g Graph) has_input(group []int, inside []bool) bool {
	return group.any(g.fed[it]) || g.edges.any(inside[it.to] && !inside[it.from])
}

// Tarjan is Tarjan's strongly connected components (SIAM J. Comput. 1, 1972).
struct Tarjan {
	g Graph
mut:
	count    int
	index    []int // 0 until visited
	low      []int
	on_stack []bool
	stack    []int
	groups   [][]int
}

fn (mut t Tarjan) visit(v int) {
	t.count++
	t.index[v] = t.count
	t.low[v] = t.count
	t.stack << v
	t.on_stack[v] = true
	for e in t.g.out[v] {
		w := t.g.edges[e].to
		if t.index[w] == 0 {
			t.visit(w)
			t.low[v] = math.min(t.low[v], t.low[w])
		} else if t.on_stack[w] {
			t.low[v] = math.min(t.low[v], t.index[w])
		}
	}
	if t.low[v] == t.index[v] {
		mut group := []int{}
		for {
			w := t.stack.pop()
			t.on_stack[w] = false
			group << w
			if w == v {
				break
			}
		}
		group.sort()
		t.groups << group
	}
}

// components is the strongly connected components, with each one after the
// ones it needs.
fn (g Graph) components() [][]int {
	mut t := Tarjan{
		g:        g
		index:    []int{len: g.n}
		low:      []int{len: g.n}
		on_stack: []bool{len: g.n}
	}
	for v in 0 .. g.n {
		if t.index[v] == 0 {
			t.visit(v)
		}
	}
	return t.groups.reverse()
}

// Circuits is Johnson's enumeration of elementary circuits (SIAM J. Comput.
// 4, 1975), by edges, so parallel streams are separate circuits. Johnson
// looks only within the component of the nodes from start up. Here a node
// that can't get back to start stays blocked instead, which costs a little
// time and finds the same circuits.
struct Circuits {
	g     Graph
	cut   []bool
	limit int
mut:
	start   int
	blocked []bool
	waiting [][]int // Johnson's B lists
	path    []int
	found   [][]int
}

// leads_on is the node an edge reaches, if the search may use the edge.
fn (c Circuits) leads_on(e int) ?int {
	w := c.g.edges[e].to
	return if c.cut[e] || w < c.start { none } else { w }
}

fn (mut c Circuits) unblock(v int) {
	c.blocked[v] = false
	waiting := c.waiting[v].clone()
	c.waiting[v] = []
	for w in waiting {
		if c.blocked[w] {
			c.unblock(w)
		}
	}
}

// circuit extends the path from v, and says if it got back to start.
fn (mut c Circuits) circuit(v int) bool {
	mut closed := false
	c.blocked[v] = true
	for e in c.g.out[v] {
		if c.found.len >= c.limit {
			return true
		}
		w := c.leads_on(e) or { continue }
		c.path << e
		if w == c.start {
			c.found << c.path.clone()
			closed = true
		} else if !c.blocked[w] && c.circuit(w) {
			closed = true
		}
		c.path.delete_last()
	}
	if closed {
		c.unblock(v)
	} else {
		for e in c.g.out[v] {
			w := c.leads_on(e) or { continue }
			if v !in c.waiting[w] {
				c.waiting[w] << v
			}
		}
	}
	return closed
}

// cycles lists a group's elementary cycles, each as its edges, leaving out
// the cut edges. It stops at limit.
fn (g Graph) cycles(group []int, cut []bool, limit int) [][]int {
	mut c := Circuits{
		g:     g
		cut:   cut
		limit: limit
	}
	for s in group {
		c.start = s
		c.blocked = []bool{len: g.n}
		c.waiting = [][]int{len: g.n}
		c.circuit(s)
	}
	return c.found
}

// order is a group's nodes in calculation order, once its tears are cut:
// each after the ones it needs from inside the group, and otherwise in
// flowsheet order.
fn (g Graph) order(group []int, tears []int) []int {
	inside := g.mask(group)
	mut waits := []int{len: g.n}
	for i, e in g.edges {
		if inside[e.from] && inside[e.to] && i !in tears {
			waits[e.to]++
		}
	}
	mut done := []bool{len: g.n}
	mut out := []int{cap: group.len}
	for out.len < group.len {
		v := group.filter(!done[it] && waits[it] == 0)[0]
		done[v] = true
		out << v
		for e in g.out[v] {
			if inside[g.edges[e].to] && e !in tears {
				waits[g.edges[e].to]--
			}
		}
	}
	return out
}
