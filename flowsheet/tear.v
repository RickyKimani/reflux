module flowsheet

// The most cycles listed for a tear search. Past this a block is cut by
// tear_greedy, which only ever needs a sample of them.
const cycle_limit = 50000

// TearChoice is what tear_set decided for one block.
struct TearChoice {
	tears  []int // edges
	cycles int   // cycles in the block
	notes  []string
}

// tear_set picks the block's tears (D5), starting from any named ones: the
// fewest streams that break every cycle, and among those the set that tears
// the fewest cycles twice (Upadhye and Grens's nonredundant sets), and then
// the earliest streams in the flowsheet. Every M4 stream has the same tear
// variables, so D5's fewest torn variables never decides.
fn (g Graph) tear_set(group []int, named []int, work i64) TearChoice {
	listed := g.cycles(group, []bool{len: g.edges.len}, cycle_limit)
	mut tears := named.clone()
	mut notes := []string{}
	if listed.len >= cycle_limit {
		tears = g.tear_greedy(group, named)
		notes << 'block of ${group.len} units has ${cycle_limit} or more cycles, too many to list, so its ${tears.len} tears came from a greedy cover and may not be the fewest'
	} else {
		mut search := new_tear_search(listed.filter(!it.any(it in named)), work)
		search.search()
		tears << search.best
		if search.stopped {
			notes << 'the search for tears stopped at its work limit, so the ${tears.len} tears of the block of ${group.len} units may not be the fewest'
		}
	}
	tears.sort()
	return TearChoice{
		tears:  tears
		cycles: listed.len
		notes:  notes
	}
}

// tear_greedy cuts the stream in the most cycles, again and again, until no
// cycle is left.
fn (g Graph) tear_greedy(group []int, named []int) []int {
	mut cut := []bool{len: g.edges.len}
	mut tears := named.clone()
	for e in named {
		cut[e] = true
	}
	for {
		listed := g.cycles(group, cut, cycle_limit)
		if listed.len == 0 {
			return tears
		}
		mut count := []int{len: g.edges.len}
		for cycle in listed {
			for e in cycle {
				count[e]++
			}
		}
		mut best := listed[0][0]
		for e, c in count {
			if c > count[best] {
				best = e
			}
		}
		cut[best] = true
		tears << best
	}
	return tears
}

// TearSearch finds the best set of streams hitting every cycle, by branch
// and bound. It's a minimum hitting set, NP-hard in general, but flowsheets'
// are small. Streams in exactly the same cycles, as those in series are, make
// one class, since the choice among them changes nothing here: Upadhye and
// Grens's families.
struct TearSearch {
	cycles       [][]int // each cycle's classes, shortest first
	class_cycles [][]int // each class's cycles
	class_stream []int   // each class's earliest stream
	limit        i64
mut:
	hits       []int  // times each cycle is torn so far
	chosen     []int  // classes
	excluded   []bool // classes a sibling branch has already tried
	extra      int    // Σ (hits − 1) over the torn cycles
	found      bool
	best       []int // the streams of the best set so far, sorted
	best_extra int
	marks      []int // for bound: the last call to use each class
	stamp      int
	work       i64
	stopped    bool
}

fn new_tear_search(cycles [][]int, limit i64) TearSearch {
	// Classes by the cycles each stream is in. Earlier streams come first.
	mut memberships := map[int][]int{}
	mut streams := []int{}
	for c, cycle in cycles {
		for e in cycle {
			if e !in memberships {
				streams << e
			}
			memberships[e] << c
		}
	}
	streams.sort()
	mut class_of := map[string]int{}
	mut stream_class := map[int]int{}
	mut class_cycles := [][]int{}
	mut class_stream := []int{}
	for e in streams {
		key := memberships[e].map(it.str()).join(',')
		if key !in class_of {
			class_of[key] = class_cycles.len
			class_cycles << memberships[e]
			class_stream << e
		}
		stream_class[e] = class_of[key]
	}
	// Shortest cycles first: they branch least, and self-loops decide at once.
	mut order := []int{len: cycles.len, init: index}
	order.sort_with_compare(fn [cycles] (a &int, b &int) int {
		x, y := cycles[*a].len, cycles[*b].len
		return if x != y { x - y } else { *a - *b }
	})
	mut renumber := []int{len: cycles.len}
	for new_index, old in order {
		renumber[old] = new_index
	}
	mut classes := [][]int{cap: cycles.len}
	for old in order {
		classes << cycles[old].map(stream_class[it])
	}
	for mut list in class_cycles {
		list = list.map(renumber[it])
	}
	return TearSearch{
		cycles:       classes
		class_cycles: class_cycles
		class_stream: class_stream
		limit:        limit
		hits:         []int{len: cycles.len}
		excluded:     []bool{len: class_cycles.len}
		marks:        []int{len: class_cycles.len}
	}
}

// take adds a class to the tears.
fn (mut s TearSearch) take(class int) {
	s.chosen << class
	for c in s.class_cycles[class] {
		if s.hits[c] > 0 {
			s.extra++
		}
		s.hits[c]++
	}
}

// drop takes a class back out.
fn (mut s TearSearch) drop(class int) {
	s.chosen.delete_last()
	for c in s.class_cycles[class] {
		s.hits[c]--
		if s.hits[c] > 0 {
			s.extra--
		}
	}
}

// bound is a lower bound on the tears still needed: cycles untouched so far
// that share no class need one tear each.
fn (mut s TearSearch) bound() int {
	s.stamp++
	mut n := 0
	for c, cycle in s.cycles {
		if s.hits[c] == 0 && !cycle.any(s.marks[it] == s.stamp) {
			n++
			for class in cycle {
				s.marks[class] = s.stamp
			}
		}
	}
	s.work += s.cycles.len
	return n
}

// earlier is whether a comes before b, element by element. V has == on
// arrays but no <.
fn earlier(a []int, b []int) bool {
	for i, x in a {
		if x != b[i] {
			return x < b[i]
		}
	}
	return false
}

// consider keeps the current set if it beats the best so far.
fn (mut s TearSearch) consider() {
	mut set := s.chosen.map(s.class_stream[it])
	set.sort()
	better := !s.found || set.len < s.best.len
		|| (set.len == s.best.len && (s.extra < s.best_extra
			|| (s.extra == s.best_extra && earlier(set, s.best))))
	if better {
		s.found = true
		s.best = set
		s.best_extra = s.extra
	}
}

// search branches on the first cycle not yet torn, one branch for each of
// its classes. A class tried in one branch is excluded from the next, so no
// set is reached twice.
fn (mut s TearSearch) search() {
	if s.work > s.limit && s.found {
		s.stopped = true
		return
	}
	mut open := -1
	for c, h in s.hits {
		if h == 0 {
			open = c
			break
		}
	}
	if open < 0 {
		s.consider()
		return
	}
	if s.found {
		least := s.chosen.len + s.bound()
		if least > s.best.len || (least == s.best.len && s.extra > s.best_extra) {
			return
		}
	}
	mut tried := []int{}
	for class in s.cycles[open] {
		if s.excluded[class] {
			continue
		}
		s.take(class)
		s.search()
		s.drop(class)
		s.excluded[class] = true
		tried << class
	}
	for class in tried {
		s.excluded[class] = false
	}
}
