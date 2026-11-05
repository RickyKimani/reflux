@[has_globals]
module capi

import sync

// Item is anything a handle refers to: a database or a fluid, from step 2.
interface Item {}

// Registry keeps what Python holds by handle (D4). The collector can't see
// Python's references, so an object only Python held would be freed under it
// (M1's spike); here it stays reachable from V until it's freed. Ids come from
// a counter and are never used twice, so a stale or doubly freed id is an
// error, not a crash. It's the one state the library shares between threads,
// so a lock guards it.
struct Registry {
mut:
	mu    &sync.Mutex = sync.new_mutex()
	next  i64
	items map[i64]Item
}

// registry needs @[has_globals] at the top of this file, which allows globals
// in this file alone, where -enable-globals would allow them everywhere.
__global registry = Registry{}

// add keeps an item, and returns its id.
fn add(item Item) i64 {
	registry.mu.lock()
	defer {
		registry.mu.unlock()
	}
	registry.next++
	registry.items[registry.next] = item
	return registry.next
}

// get returns the item with an id, or says why there's none.
fn get(id i64) !Item {
	registry.mu.lock()
	defer {
		registry.mu.unlock()
	}
	return registry.items[id] or { return error(missing(id)) }
}

// remove drops an id, or says why it can't.
fn remove(id i64) ! {
	registry.mu.lock()
	defer {
		registry.mu.unlock()
	}
	if id !in registry.items {
		return error(missing(id))
	}
	registry.items.delete(id)
}

// missing says why an id has no item: freed already, if it was ever given
// out, since ids aren't reused. It's called with the lock held.
fn missing(id i64) string {
	if id >= 1 && id <= registry.next {
		return 'handle ${id} was freed already'
	}
	return '${id} is not a handle'
}
