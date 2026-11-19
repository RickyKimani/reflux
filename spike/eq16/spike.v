// M8 step 1 spike: what K(T) costs at a T that moves, and whether ChemSep's
// equation 16, Cp = A + exp(B/T + C + D·T + E·T²), which has no closed-form
// integral, could be integrated faster and still to tdep's 1e-13.
//
//     v -prod -cc gcc run spike/eq16
//
// 1. What a fill costs the published reforming model at a T that moves, and
//    what its parts cost: one reaction's standard properties, one species'
//    two integrals, and the integrand's evaluations each takes.
// 2. A fixed Gauss–Legendre rule on [298.15 K, T] in place of the adaptive
//    Gauss–Kronrod: for each number of nodes, how many of the database's
//    equation-16 compounds it holds to 1e-13 at every T tried, and the worst.
//    One pass of it gives ∫Cp dT and ∫Cp/T dT both.
//
// It's run from the repository root. Not part of the simulator.
module main

import components
import flowsheet
import kinetics
import math
import numeric
import reactors
import tdep
import thermo
import time

// Counting is a correlation that counts its evaluations.
struct Counting {
	c tdep.Correlation
	n &int
}

fn (f Counting) at(t f64) f64 {
	unsafe {
		(*f.n)++
	}
	return f.c.eval(t)
}

// OverT is Y(T)/T, as tdep integrates it for entropies.
struct OverT {
	f Counting
}

fn (f OverT) at(t f64) f64 {
	return f.f.at(t) / t
}

// legendre is Gauss–Legendre's nodes and weights on [−1, 1], by Newton's
// method on Pₙ from Chebyshev's first guesses.
fn legendre(n int) ([]f64, []f64) {
	mut x := []f64{len: n}
	mut w := []f64{len: n}
	for i in 0 .. n {
		mut z := math.cos(math.pi * (f64(i) + 0.75) / (f64(n) + 0.5))
		mut dp := 0.0
		for _ in 0 .. 100 {
			mut p0 := 1.0
			mut p1 := z
			for k in 2 .. n + 1 {
				p2 := ((2 * f64(k) - 1) * z * p1 - (f64(k) - 1) * p0) / f64(k)
				p0 = p1
				p1 = p2
			}
			dp = f64(n) * (z * p1 - p0) / (z * z - 1)
			dz := p1 / dp
			z -= dz
			if math.abs(dz) < 1e-16 {
				break
			}
		}
		x[i] = z
		w[i] = 2 / ((1 - z * z) * dp * dp)
	}
	return x, w
}

// fixed is ∫ Y dT and ∫ Y/T dT from a to b by one rule, in one pass.
fn fixed(c tdep.Correlation, a f64, b f64, x []f64, w []f64) (f64, f64) {
	half := 0.5 * (b - a)
	mid := 0.5 * (a + b)
	mut h := 0.0
	mut s := 0.0
	for i, xi in x {
		t := mid + half * xi
		y := c.eval(t)
		h += w[i] * y
		s += w[i] * y / t
	}
	return h * half, s * half
}

// composite is the same over m equal pieces: the reference.
fn composite(c tdep.Correlation, a f64, b f64, x []f64, w []f64, m int) (f64, f64) {
	mut h := 0.0
	mut s := 0.0
	for k in 0 .. m {
		lo := a + (b - a) * f64(k) / f64(m)
		hi := a + (b - a) * f64(k + 1) / f64(m)
		dh, ds := fixed(c, lo, hi, x, w)
		h += dh
		s += ds
	}
	return h, s
}

fn ns_since(t0 u64, calls int) f64 {
	return f64(time.sys_mono_now() - t0) / f64(calls)
}

fn eq16(gas thermo.IdealGas) ?tdep.Correlation {
	cp := gas.cp
	if cp is thermo.Correlated {
		if cp.c.eq == 16 {
			return cp.c
		}
	}
	return none
}

fn costs() ! {
	println('1. What K(T) costs the published model (testdata/flowsheets/reformer.toml) at a T that moves')
	l := flowsheet.load('testdata/flowsheets/reformer.toml', makers: reactors.makers())!
	for p in l.sheet.units {
		u := p.unit
		if u is reactors.Pfr {
			k := u.k
			n := k.species.len
			y := []f64{len: n, init: 1.0 / n}
			mut pt := k.point()
			calls := 2000
			mut t0 := time.sys_mono_now()
			mut before := gc_heap_usage().total_bytes
			for i in 0 .. calls {
				k.fill(kinetics.State{ t: 900.0 + 0.01 * i, p: 1e5, y: y }, mut pt)!
			}
			moving := ns_since(t0, calls)
			bytes := f64(gc_heap_usage().total_bytes - before) / calls
			t0 = time.sys_mono_now()
			for i in 0 .. calls {
				k.fill(kinetics.State{ t: 900.0, p: 1e5 + i, y: y }, mut pt)!
			}
			still := ns_since(t0, calls)
			println('  a fill at a new T: ${moving / 1e3:.1f} µs and ${bytes:.0f} bytes allocated, for ${k.k_used.len} reactions over ${n} species; at an unmoved T, ${still / 1e3:.3f} µs')
			t0 = time.sys_mono_now()
			for i in 0 .. calls {
				k.reactions[k.k_used[0]].standard(k.species, 900.0 + 0.01 * i)!
			}
			one := ns_since(t0, calls)
			mut integrals := 0
			for r in k.k_used {
				integrals += 2 * k.reactions[r].species.len
			}
			println("  one reaction's standard properties: ${one / 1e3:.1f} µs; the fill integrates ${integrals} times, two for each species of each reaction, where ${2 * n} would do, two for each species")
			for s in k.species {
				gas := s.gas or { continue }
				c := eq16(gas) or {
					t0 = time.sys_mono_now()
					mut sum := 0.0
					for i in 0 .. calls {
						sum += gas.cp.int_cp(thermo.t_ref, 900.0 + 0.01 * i)
						sum += gas.cp.int_cp_over_t(thermo.t_ref, 900.0 + 0.01 * i)
					}
					println('  ${s.name}: a closed form, ${ns_since(t0, calls):.0f} ns for both integrals (${sum:.3e})')
					continue
				}
				mut count := 0
				f := Counting{c, &count}
				for t in [600.0, 900.0, 1300.0] {
					count = 0
					numeric.integrate(f, thermo.t_ref, t, 1e-13)!
					h_n := count
					count = 0
					numeric.integrate(OverT{f}, thermo.t_ref, t, 1e-13)!
					println('  ${s.name} to ${t:.0f} K: ${h_n} evaluations for ∫Cp dT, ${count} for ∫Cp/T dT')
				}
				t0 = time.sys_mono_now()
				before = gc_heap_usage().total_bytes
				mut sum := 0.0
				for i in 0 .. calls {
					sum += gas.cp.int_cp(thermo.t_ref, 900.0 + 0.01 * i)
					sum += gas.cp.int_cp_over_t(thermo.t_ref, 900.0 + 0.01 * i)
				}
				both := ns_since(t0, calls)
				alloc := f64(gc_heap_usage().total_bytes - before) / calls
				x, w := legendre(16)
				t0 = time.sys_mono_now()
				for i in 0 .. calls {
					dh, ds := fixed(c, thermo.t_ref, 900.0 + 0.01 * i, x, w)
					sum += dh + ds
				}
				println('  ${s.name}: both integrals ${both / 1e3:.2f} µs and ${alloc:.0f} bytes; by one 16-node pass ${ns_since(t0, calls) / 1e3:.2f} µs (${sum:.3e})')
			}
		}
	}
}

fn accuracy() ! {
	println("2. A fixed Gauss–Legendre rule from 298.15 K, over the database's equation-16 compounds")
	db := components.load_default()!
	ts := [200.0, 250.0, 350.0, 500.0, 700.0, 900.0, 1100.0, 1300.0, 1500.0]
	rx, rw := legendre(32)
	mut cs := []tdep.Correlation{}
	mut names := []string{}
	for e in db.list() {
		r := db.find(e.cas)!
		gas := r.ideal_gas() or { continue }
		c := eq16(gas) or { continue }
		cs << c
		names << e.name
	}
	println('  ${cs.len} compounds, each to ${ts.map(int(it))} K; the reference is 32 nodes on each of 16 pieces')
	// The adaptive rule's own error against the reference, first.
	mut worst_now := 0.0
	mut failed := 0
	for c in cs {
		for t in ts {
			h, s := composite(c, thermo.t_ref, t, rx, rw, 16)
			a := c.integral(thermo.t_ref, t) or {
				failed++
				continue
			}
			b := c.integral_over_t(thermo.t_ref, t) or {
				failed++
				continue
			}
			worst_now = math.max(worst_now, math.max(math.abs(a / h - 1), math.abs(b / s - 1)))
		}
	}
	println("  tdep's adaptive integrals: within ${worst_now:.1e} of the reference at worst, ${failed} failed")
	for n in [8, 12, 16, 20, 24, 32, 48] {
		x, w := legendre(n)
		mut held := 0
		mut worst := 0.0
		mut worst_name := ''
		mut held_1000 := 0
		for q, c in cs {
			mut mine := 0.0
			mut mine_1000 := 0.0
			for t in ts {
				h, s := composite(c, thermo.t_ref, t, rx, rw, 16)
				a, b := fixed(c, thermo.t_ref, t, x, w)
				e := math.max(math.abs(a / h - 1), math.abs(b / s - 1))
				mine = math.max(mine, e)
				if t <= 1100 {
					mine_1000 = math.max(mine_1000, e)
				}
			}
			if mine <= 1e-13 {
				held++
			}
			if mine_1000 <= 1e-13 {
				held_1000++
			}
			if mine > worst || math.is_nan(mine) {
				worst = mine
				worst_name = names[q]
			}
		}
		println('  ${n:2} nodes: ${held} of ${cs.len} compounds within 1e-13 at every T (${held_1000} up to 1100 K); the worst ${worst:.1e}, ${worst_name}')
	}
}

fn main() {
	costs() or {
		eprintln(err.msg())
		exit(1)
	}
	accuracy() or {
		eprintln(err.msg())
		exit(1)
	}
}
