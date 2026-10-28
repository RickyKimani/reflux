// Package sweep is M3's robustness sweep (decision D8, exit criterion 5):
// PT flashes of seeded random PR mixtures drawn from the component
// database, each answer checked from scratch and each failure classified.
// tools/flash_sweep.vsh runs 10,000 points and writes
// docs/flash-robustness.md; sweep_test.v runs the first 1,000 in CI.
module sweep

import components
import flash
import math
import numeric
import thermo

pub const seed = u64(20261028)

// The ranges: T and P log-uniform, 2 to 10 components, z uniform over the
// simplex.
pub const t_min = 100.0
pub const t_max = 1000.0
pub const p_min = 1e4
pub const p_max = 1e8
pub const nc_min = 2
pub const nc_max = 10

// Rng is splitmix64 (Steele, Lea & Flood 2014): small, fast, and the same
// on every platform and V version, so a seed always gives the same points.
pub struct Rng {
mut:
	state u64
}

// u64 is the next 64 random bits.
pub fn (mut r Rng) u64() u64 {
	r.state += 0x9e3779b97f4a7c15
	mut z := r.state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

// f64 is uniform in (0, 1): 53 random bits, never 0.
pub fn (mut r Rng) f64() f64 {
	return (f64(r.u64() >> 11) + 0.5) / 9007199254740992.0
}

// below is uniform in 0 .. n.
pub fn (mut r Rng) below(n int) int {
	return int(r.u64() % u64(n))
}

// Point is one state of the sweep: the compounds by CAS number, the feed,
// T in K and P in Pa.
pub struct Point {
pub:
	index int
	cas   []string
	names []string
	z     []f64
	t     f64
	p     f64
}

// candidates are the compounds the sweep draws from: those with Tc, Pc, ω,
// M and an ideal-gas Cp, by CAS number.
pub fn candidates(db components.Database) ![]components.IndexEntry {
	mut out := []components.IndexEntry{}
	for e in db.list() {
		r := db.find(e.cas)!
		r.component() or { continue }
		r.ideal_gas() or { continue }
		out << e
	}
	return out
}

// points draws count points from seed.
pub fn points(pool []components.IndexEntry, count int) []Point {
	mut rng := Rng{seed}
	mut out := []Point{}
	for k in 0 .. count {
		nc := nc_min + rng.below(nc_max - nc_min + 1)
		mut picked := []int{}
		for picked.len < nc {
			i := rng.below(pool.len)
			if i !in picked {
				picked << i
			}
		}
		picked.sort()
		mut z := picked.map(-math.log(rng.f64()))
		total := numeric.sum(z)
		z = z.map(it / total)
		t := math.exp(math.log(t_min) + rng.f64() * math.log(t_max / t_min))
		p := math.exp(math.log(p_min) + rng.f64() * math.log(p_max / p_min))
		out << Point{
			index: k
			cas:   picked.map(pool[it].cas)
			names: picked.map(pool[it].name)
			z:     z
			t:     t
			p:     p
		}
	}
	return out
}

// Class is what became of a point.
pub enum Class {
	single // one phase, checked stable
	split  // two phases, every invariant of exit criterion 3 checked
	// The explained failures (D8):
	three_phase   // the flash found a third phase (more_phases): M6's job
	near_critical // the split collapsed onto the feed: a critical point
	// Anything else is unexplained: no convergence, a broken invariant, or a
	// contradiction.
	unexplained
}

pub struct Outcome {
pub:
	point  Point
	class  Class
	detail string
}

// run flashes a point and classifies the answer.
pub fn run(db components.Database, pt Point) Outcome {
	fl := db.fluid(.pr, pt.cas) or {
		return Outcome{pt, .unexplained, 'no fluid: ${err}'}
	}
	m := fl.model
	r := flash.pt(m, pt.t, pt.p, pt.z) or {
		if err is flash.FlashError {
			e := err as flash.FlashError
			if e.kind == .more_phases {
				proof := prove_more_phases(m, pt) or {
					return Outcome{pt, .unexplained, "a third phase the sweep can't confirm: ${err}"}
				}
				return Outcome{pt, .three_phase, proof}
			}
			if e.kind == .no_convergence && e.reason.contains('collapsed onto the feed') {
				return Outcome{pt, .near_critical, e.reason}
			}
		}
		return Outcome{pt, .unexplained, err.msg()}
	}
	if r.phases.len == 1 {
		check_single(m, pt, r) or { return Outcome{pt, .unexplained, err.msg()} }
		return Outcome{pt, .single, ''}
	}
	check_split(m, pt, r) or { return Outcome{pt, .unexplained, err.msg()} }
	return Outcome{pt, .split, ''}
}

// trials are five seeded random compositions, an independent probe of a
// verdict on top of stability's own trials.
fn trials(pt Point) [][]f64 {
	mut rng := Rng{seed ^ u64(pt.index + 1) * 0x2545f4914f6cdd1d}
	mut out := [][]f64{}
	for _ in 0 .. 5 {
		mut w := pt.z.map(-math.log(rng.f64()))
		total := numeric.sum(w)
		out << w.map(it / total)
	}
	return out
}

// check_single re-tests a single-phase answer for stability, with extra
// random trials and to the flash's own tolerance, 1e-10: a composition
// below the feed's tangent plane contradicts the flash.
fn check_single(m thermo.Model, pt Point, r flash.FlashResult) ! {
	lowest := lowest_tpd(m, pt.t, pt.p, r.z, pt, 1e-10)!
	if lowest < -1e-10 {
		return error('a single-phase answer with a phase below its tangent plane (tm = ${lowest})')
	}
}

// tpd is the normalised tangent-plane distance of w against a phase x,
// Σ wᵢ(ln wᵢ + ln φᵢ(w) − ln xᵢ − ln φᵢ(x)), straight from the model and
// independent of flash.stability. Negative proves x unstable.
fn tpd(m thermo.Model, t f64, p f64, x []f64, w []f64) !f64 {
	sx := thermo.solve_tp(m, t, p, x, .stable)!
	sw := thermo.solve_tp(m, t, p, w, .stable)!
	mut sum := 0.0
	for i in 0 .. x.len {
		if w[i] > 0.0 {
			sum += w[i] * (math.log(w[i]) + sw.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
		}
	}
	return sum
}

// lowest_tpd is the lowest tangent-plane distance against x, by tpd, at
// every composition that flash.stability's trials reached, the seeded
// random ones among them. Stability is only the search: the verdict is
// tpd's, so a wrong verdict from stability can't pass. (A phase of a
// converged split is no use as a trial for the other: it sits on that
// phase's tangent plane, a stationary point with tm = 0.)
fn lowest_tpd(m thermo.Model, t f64, p f64, x []f64, pt Point, tol f64) !f64 {
	s := flash.stability(m, t, p, x, trials: trials(pt), tol: tol)!
	mut lowest := 0.0
	for tr in s.trials {
		if !tr.trivial {
			lowest = math.min(lowest, tpd(m, t, p, x, tr.w)!)
		}
	}
	return lowest
}

// prove_more_phases confirms a third-phase verdict from scratch, and says
// which phase is unstable and by how much. Without verification the flash
// returns the first split it finds, which is the one the verified flash
// found first and refused. One of its phases must then have a composition
// below its tangent plane by more than 1e-8.
fn prove_more_phases(m thermo.Model, pt Point) !string {
	r := flash.pt(m, pt.t, pt.p, pt.z, verify: false)!
	if r.phases.len != 2 {
		return error('unverified, the flash gives one phase')
	}
	for ph in r.phases {
		lowest := lowest_tpd(m, r.t, r.p, ph.state.n, pt, 1e-8)!
		if lowest < -1e-8 {
			name := if ph.label == .vapor { 'vapour' } else { 'liquid' }
			return 'the ${name}, tm = ${round_sig(lowest, 4)}'
		}
	}
	return error('neither phase of the split is unstable')
}

// check_split checks exit criterion 3 on a split, from scratch.
fn check_split(m thermo.Model, pt Point, r flash.FlashResult) ! {
	vap, liq := r.phases[0], r.phases[1]
	if !(r.beta > 0.0 && r.beta < 1.0) {
		return error('a split with β = ${r.beta}')
	}
	y, x := vap.state.n, liq.state.n
	// Two copies of the feed pass every other check, so the phases must
	// differ by more than the flash's own test for a collapsed split.
	mut distinct := false
	for i, z in r.z {
		if z > 0.0 && !(math.abs(math.log(y[i] / x[i])) < 1e-4) {
			distinct = true
		}
	}
	if !distinct {
		return error('the two phases are the same')
	}
	sy := thermo.solve_tp(m, r.t, r.p, y, .stable)!
	sx := thermo.solve_tp(m, r.t, r.p, x, .stable)!
	mut sum_x, mut sum_y := 0.0, 0.0
	for i, z in r.z {
		sum_x += x[i]
		sum_y += y[i]
		if math.abs(r.beta * y[i] + (1.0 - r.beta) * x[i] - z) > 1e-12 {
			return error('the mass balance of component ${i} is off')
		}
		// Below flash.trace_floor in a phase a component's amount there is
		// zero to double precision, and its fugacities can't agree.
		if z > 0.0 && x[i] >= flash.trace_floor && y[i] >= flash.trace_floor {
			res := math.abs(math.log(y[i]) + sy.ln_phi[i] - math.log(x[i]) - sx.ln_phi[i])
			if res > 1e-10 {
				return error('|ln fV − ln fL| of component ${i} is ${res}')
			}
		}
	}
	if math.abs(sum_x - 1.0) > 1e-14 || math.abs(sum_y - 1.0) > 1e-14 {
		return error('Σx − 1 = ${sum_x - 1.0}, Σy − 1 = ${sum_y - 1.0}')
	}
	for ph in [x, y] {
		lowest := lowest_tpd(m, r.t, r.p, ph, pt, 1e-8)!
		if lowest < -1e-8 {
			return error('a phase of the split is unstable (tm = ${lowest})')
		}
	}
}
