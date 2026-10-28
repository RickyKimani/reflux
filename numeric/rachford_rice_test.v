// Checks rachford_rice against testdata/m3_rachford_rice.json, written by
// tools/gen_m3_fixtures.py (run by hand with uv). The reference is the
// exact root for each case's floats, found by bisection at 100 digits and
// rounded once. The fixture also records what chemicals' solvers return.
import json2
import math
import numeric
import os

struct Source {
	generated string
	generator string
	python    string
	chemicals string
	mpmath    string
	seed      i64
}

struct Split {
	beta   f64
	liquid f64
	x      []f64
	y      []f64
}

struct Other {
	method string
	beta   ?f64
	x      []f64
	y      []f64
	error  string
}

struct Case {
	name      string
	z         []f64
	k         []f64
	reference Split
	chemicals []Other
}

struct Fixture {
	source Source
	note   string
	cases  []Case
}

fn read_fixture() !Fixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm3_rachford_rice.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m3_fixtures.py`.') }
	return json2.decode[Fixture](text)!
}

// Errors is how far one split is from the reference: β and L absolutely,
// and the compositions relatively, component by component.
struct Errors {
	beta        f64
	liquid      f64
	composition f64
}

fn errors(beta f64, liquid f64, x []f64, y []f64, ref Split) Errors {
	mut worst := 0.0
	for i in 0 .. x.len {
		for pair in [[x[i], ref.x[i]], [y[i], ref.y[i]]] {
			worst = math.max(worst, if pair[1] != 0.0 {
				math.abs(pair[0] - pair[1]) / math.abs(pair[1])
			} else {
				math.abs(pair[0])
			})
		}
	}
	return Errors{math.abs(beta - ref.beta), math.abs(liquid - ref.liquid), worst}
}

// window_scale is how far β and L can be trusted in absolute terms,
// beyond their own size. Each 1 + β(Kᵢ − 1) is accurate to a few ulps,
// and β = (eᵢ − 1)/(Kᵢ − 1) is best read from the K furthest from 1, which
// carries an error of about ε·|β + 1/(Kᵢ − 1)|: ε·(|β| + min(|c₁|, cₙ)),
// with c₁ and cₙ the window's ends. min(|c₁|, cₙ) is at most 1 for an
// ordinary mixture, and 1e8 when every K is within 1e-8 of 1. L = 1 − β
// inherits the same absolute error, plus ε·|L| of its own.
fn window_scale(c Case) f64 {
	mut kmax, mut kmin := 0.0, math.inf(1)
	for i, z in c.z {
		if z > 0.0 {
			kmax = math.max(kmax, c.k[i])
			kmin = math.min(kmin, c.k[i])
		}
	}
	return math.min(1.0 / (kmax - 1.0), 1.0 / (1.0 - kmin))
}

const eps = 2.220446049250313e-16

fn test_rachford_rice_matches_the_exact_root() {
	f := read_fixture()!
	assert f.cases.len == 87
	for c in f.cases {
		r := numeric.rachford_rice(c.z, c.k)!
		e := errors(r.beta, r.liquid, r.x, r.y, c.reference)
		// Every composition to a few ulps, extreme K-values included.
		assert e.composition <= 1e-15, '${c.name}: compositions off by ${e.composition}'
		// β and L to a few ulps of what the K-values allow.
		scale := window_scale(c)
		beta_bound := 4.0 * eps * (math.abs(c.reference.beta) + scale)
		liquid_bound := 4.0 * eps * (math.abs(c.reference.liquid) + scale)
		assert e.beta <= beta_bound, '${c.name}: β off by ${e.beta}, bound ${beta_bound}'
		assert e.liquid <= liquid_bound, '${c.name}: L off by ${e.liquid}, bound ${liquid_bound}'
		// A trace phase keeps its digits: L or β of 1e-15 to a few ulps of
		// itself, not of 1.
		if c.name.starts_with('trace liquid') {
			assert e.liquid <= 4.0 * eps * c.reference.liquid, '${c.name}: L off by ${e.liquid}'
		}
		if c.name.starts_with('trace vapour') {
			assert e.beta <= 4.0 * eps * c.reference.beta, '${c.name}: β off by ${e.beta}'
		}
		assert r.iterations <= 8, '${c.name}: ${r.iterations} iterations'
	}
}

fn test_rachford_rice_is_at_least_as_accurate_as_chemicals() {
	// Over the whole fixture, reflux's largest error in β and in the
	// compositions is no larger than that of any of chemicals' float
	// solvers. (Their mpmath solver is the cross-check on the reference.)
	f := read_fixture()!
	mut ours := Errors{}
	mut theirs := map[string]Errors{}
	for c in f.cases {
		r := numeric.rachford_rice(c.z, c.k)!
		e := errors(r.beta, r.liquid, r.x, r.y, c.reference)
		ours = Errors{math.max(ours.beta, e.beta), 0.0, math.max(ours.composition, e.composition)}
		for o in c.chemicals {
			if o.method == 'mpmath' {
				continue
			}
			if b := o.beta {
				oe := errors(b, 0.0, o.x, o.y, c.reference)
				w := theirs[o.method]
				theirs[o.method] = Errors{math.max(w.beta, oe.beta), 0.0, math.max(w.composition,
					oe.composition)}
			}
		}
	}
	assert theirs.len == 4
	for method, w in theirs {
		assert ours.beta <= w.beta, '${method}: β ${w.beta}, reflux ${ours.beta}'
		assert ours.composition <= w.composition, '${method}: compositions ${w.composition}, reflux ${ours.composition}'
	}
}

fn test_the_fixture_reference_agrees_with_chemicals_mpmath() {
	// The generator stops if they disagree; this checks that the committed
	// file is one it let through.
	f := read_fixture()!
	for c in f.cases {
		m := c.chemicals.filter(it.method == 'mpmath')
		assert m.len == 1, c.name
		b := m[0].beta or { panic('${c.name}: chemicals mpmath failed: ${m[0].error}') }
		assert math.abs(b - c.reference.beta) <= 1e-15 * math.max(1.0, math.abs(b)), c.name
	}
}

// Splitmix64, for seeded test cases that come out the same everywhere.
struct Rng {
mut:
	state u64
}

fn (mut r Rng) next() f64 {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	z = z ^ (z >> 31)
	return f64(z >> 11) / f64(u64(1) << 53)
}

fn test_rachford_rice_balances_mass() {
	// For seeded random feeds with K over twenty decades: zᵢ = L·xᵢ + β·yᵢ
	// component by component, Σx = Σy = 1, and y = K·x exactly. In a
	// negative flash L·xᵢ and β·yᵢ have opposite signs and can be much
	// larger than zᵢ, so the sum is checked to the rounding it carries,
	// ε·(|L|·xᵢ + |β|·yᵢ). For 0 ≤ β ≤ 1 that's ε·zᵢ.
	mut rng := Rng{20261027}
	mut done := 0
	for done < 2000 {
		n := 2 + int(rng.next() * 15)
		mut z := []f64{len: n, init: 0.01 + rng.next()}
		total := numeric.sum(z)
		for i in 0 .. n {
			z[i] /= total
		}
		k := []f64{len: n, init: math.pow(10.0, -10.0 + 20.0 * rng.next())}
		r := numeric.rachford_rice(z, k) or {
			assert err.msg().contains('no two-phase split')
			continue
		}
		done++
		for i in 0 .. n {
			assert r.y[i] == k[i] * r.x[i]
			back := r.liquid * r.x[i] + r.beta * r.y[i]
			size := math.abs(r.liquid) * r.x[i] + math.abs(r.beta) * r.y[i]
			assert math.abs(back - z[i]) <= 8.0 * eps * size, 'component ${i} of ${z} with K ${k}'
		}
		assert math.abs(numeric.sum(r.x) - 1.0) <= 1e-14
		assert math.abs(numeric.sum(r.y) - 1.0) <= 1e-14
	}
}

fn test_a_binary_takes_one_step() {
	// For two components the Leibovici–Neoschil function is linear in the
	// window position, so the first guess, interpolated between its values
	// at the two ends, is already the root.
	mut rng := Rng{7}
	for _ in 0 .. 200 {
		z0 := 0.001 + 0.998 * rng.next()
		k := [math.pow(10.0, 8.0 * rng.next()), math.pow(10.0, -8.0 * rng.next())]
		r := numeric.rachford_rice([z0, 1.0 - z0], k)!
		assert r.iterations == 1, '${z0} ${k}'
	}
}

fn test_absent_components_have_no_composition() {
	r := numeric.rachford_rice([0.0, 0.5, 0.5, 0.0], [100.0, 2.0, 0.5, 0.001])!
	assert r.x[0] == 0.0 && r.y[0] == 0.0 && r.x[3] == 0.0 && r.y[3] == 0.0
	// The absent K-values don't set the window: this is the binary
	// 0.5/0.5 with K = 2 and 0.5, whose root is β = 1/2.
	assert math.abs(r.beta - 0.5) <= 2.0 * eps
}

fn test_rachford_rice_refuses_what_it_cannot_split() {
	for args in [
		[[0.5, 0.5], [2.0, 3.0]], // every K above 1
		[[0.5, 0.5], [0.2, 1.0]], // every K at most 1
		[[1.0], [2.0]], // one component
		[[0.0, 0.0], [2.0, 0.5]], // nothing in the feed
		[[0.5, -0.5], [2.0, 0.5]], // a negative mole fraction
		[[0.5, 0.5], [2.0, 0.0]], // K = 0
		[[0.5, 0.5], [math.nan(), 0.5]],
		[[0.5, 0.5], [2.0]], // lengths differ
	] {
		if r := numeric.rachford_rice(args[0], args[1]) {
			assert false, 'z ${args[0]}, K ${args[1]} gave β = ${r.beta}'
		}
	}
	// The absent component's K = 0.1 doesn't count: both present K are above 1.
	if _ := numeric.rachford_rice([0.5, 0.5, 0.0], [2.0, 3.0, 0.1]) {
		assert false, 'only absent components below 1'
	} else {
		assert err.msg() == 'rachford_rice: every K is at least 1: all vapour, so there is no two-phase split'
	}
}
