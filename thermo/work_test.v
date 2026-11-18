// CubicWork (M8 step 7): solve_tp's bits for any cubic, root and state,
// from one workspace reused across T and n, with nothing allocated once
// it's made.
import math
import rand
import thermo { Component, CubicKind, Root }

const methane = Component{
	name:  'methane'
	tc:    190.6
	pc:    45.99e5
	omega: 0.012
	mw:    0.016043
}

const butane = Component{
	name:  'n-butane'
	tc:    425.1
	pc:    37.96e5
	omega: 0.200
	mw:    0.058123
}

const pentane = Component{
	name:  'n-pentane'
	tc:    469.7
	pc:    33.70e5
	omega: 0.252
	mw:    0.072150
}

fn mixture(kind CubicKind) thermo.Cubic {
	return thermo.new_cubic(kind, [methane, butane, pentane],
		kij: [0.0, 0.02, 0.03, 0.02, 0.0, 0.004, 0.03, 0.004, 0.0]
	) or { panic(err) }
}

fn test_a_workspace_gives_solve_tps_bits() {
	rand.seed([u32(7), 17])
	mut compared := 0
	for kind in [CubicKind.vdw, .rk, .srk, .pr] {
		c := mixture(kind)
		mut w := c.work()
		for _ in 0 .. 400 {
			// A few temperatures, each for several compositions, as a
			// reactor's states come.
			t := 150.0 + 450.0 * rand.f64()
			for _ in 0 .. 5 {
				p := math.pow(10.0, 4.0 + 3.0 * rand.f64())
				n := [rand.f64(), rand.f64(), rand.f64()]
				root := [Root.liquid, .vapor, .stable][rand.intn(3) or { 0 }]
				want := thermo.solve_tp(c, t, p, n, root) or {
					w.solve_tp(t, p, n, root) or {
						assert err.msg().contains('no volume root')
						continue
					}
					assert false, 'the workspace solved what solve_tp refused'
					continue
				}
				got := w.solve_tp(t, p, n, root)!
				assert math.f64_bits(got.v) == math.f64_bits(want.v), '${kind} ${t} ${p} ${n}'
				assert math.f64_bits(got.z) == math.f64_bits(want.z)
				for i, x in want.ln_phi {
					assert math.f64_bits(w.ln_phi[i]) == math.f64_bits(x), '${kind} ${i}'
				}
				compared++
			}
		}
	}
	assert compared > 7000, '${compared}'
}

fn test_a_workspace_refuses_what_solve_tp_does() {
	c := mixture(.pr)
	mut w := c.work()
	for args in [[0.0, 1e5, 1.0], [300.0, -1.0, 1.0], [300.0, 1e5, -1.0]] {
		n := [args[2], 1.0, 1.0]
		mut got := ''
		w.solve_tp(args[0], args[1], n, .vapor) or { got = err.msg() }
		mut want := ''
		thermo.solve_tp(c, args[0], args[1], n, .vapor) or { want = err.msg() }
		assert got != '' && got == want, '${args}: ${got} | ${want}'
	}
	w.solve_tp(300.0, 1e5, [1.0, 1.0], .vapor) or {
		assert err.msg().contains('composition has 2 entries but the model has 3 components')
		return
	}
	assert false
}

fn test_a_workspace_allocates_nothing() {
	c := mixture(.pr)
	mut w := c.work()
	mut n := [0.5, 0.3, 0.2]
	w.solve_tp(800.0, 2e5, n, .vapor)!
	before := u64(gc_heap_usage().total_bytes)
	for k in 0 .. 1000 {
		n[0] = 0.5 + 1e-4 * f64(k)
		w.solve_tp(800.0, 2e5, n, .vapor)!
	}
	// A string or two of the collector's own, not a byte a call.
	assert u64(gc_heap_usage().total_bytes) - before < 1000
}
