// Checks every supported correlation form against testdata/m2_tdep.json,
// written by tools/gen_m2_fixtures.py (run by hand with uv). The reference
// values are mpmath at 40 digits from formulas written independently of
// this module, on real ChemSep coefficient sets. For the DIPPR forms the
// fixture also carries chemicals' values.
import json2
import math
import os
import tdep

struct Source {
	generated string
	generator string
	python    string
	chemicals string
	mpmath    string
	chemsep   string
}

struct Point {
	t               f64
	y               f64
	dy              f64
	integral        f64
	integral_over_t f64
	// Only present for the DIPPR forms, and where chemicals didn't overflow.
	chemicals_y  ?f64
	chemicals_dy ?f64
}

struct Case {
	eq     int
	coeffs []f64
	tmin   f64
	tmax   f64
	tc     f64
	origin string
	t0     f64
	points []Point
}

struct Fixture {
	source                    Source
	note                      string
	unsupported_forms_in_data []int
	cases                     []Case
}

fn read_fixture() !Fixture {
	path := os.join_path(os.dir(@FILE), '..', 'testdata', 'm2_tdep.json')
	text := os.read_file(path) or { return error('${path}: ${err}. Regenerate with `uv run tools/gen_m2_fixtures.py`.') }
	return json2.decode[Fixture](text)!
}

// agree checks got against want to a relative tol, judged on the larger of
// |want| and scale, so that values near zero aren't held to a relative
// standard they can't meet.
fn agree(got f64, want f64, tol f64, scale f64) bool {
	return math.abs(got - want) <= tol * math.max(math.abs(want), math.max(scale, 1e-300))
}

fn test_every_supported_form_has_fixtures() {
	f := read_fixture()!
	mut seen := map[int]bool{}
	for c in f.cases {
		seen[c.eq] = true
	}
	for eq in [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 16, 17, 45, 100, 101, 102, 104, 105, 106,
		107, 114, 115, 116] {
		assert tdep.is_supported(eq)
		assert seen[eq], 'no fixture for equation ${eq}'
	}
	// The forms ChemSep uses that reflux can't evaluate yet. A new one
	// appearing means the vendored data changed: decide whether to support it.
	assert f.unsupported_forms_in_data == [118, 120, 121, 124, 200]
	for eq in f.unsupported_forms_in_data {
		assert !tdep.is_supported(eq)
	}
}

fn test_values_derivatives_and_integrals_match_mpmath() {
	f := read_fixture()!
	mut checked := 0
	for c in f.cases {
		corr := tdep.new(c.eq, c.coeffs, c.tmin, c.tmax, c.tc) or {
			assert false, '${c.origin}: ${err}'
			continue
		}
		for p in c.points {
			what := 'eq ${c.eq} ${c.origin} at ${p.t} K'
			y := corr.eval(p.t)
			assert agree(y, p.y, 1e-13, 0.0), '${what}: Y = ${y}, want ${p.y}'
			// The derivative's natural scale is Y/T; the integrals' are
			// Y·ΔT and Y·Δln T.
			dy := corr.deriv(p.t)
			assert agree(dy, p.dy, 1e-11, math.abs(p.y) / p.t), '${what}: dY/dT = ${dy}, want ${p.dy}'
			integral := corr.integral(c.t0, p.t)!
			assert agree(integral, p.integral, 1e-11, math.abs(p.y) * math.abs(p.t - c.t0)), '${what}: ∫Y dT = ${integral}, want ${p.integral}'
			over_t := corr.integral_over_t(c.t0, p.t)!
			assert agree(over_t, p.integral_over_t, 1e-11, math.abs(p.y) * math.abs(math.log(p.t / c.t0))), '${what}: ∫Y/T dT = ${over_t}, want ${p.integral_over_t}'
			checked++
		}
	}
	assert checked >= 250
}

fn test_dippr_forms_match_chemicals() {
	f := read_fixture()!
	mut checked := 0
	for c in f.cases {
		corr := tdep.new(c.eq, c.coeffs, c.tmin, c.tmax, c.tc)!
		for p in c.points {
			what := 'eq ${c.eq} ${c.origin} at ${p.t} K'
			if want := p.chemicals_y {
				assert agree(corr.eval(p.t), want, 1e-12, 0.0), '${what}: Y = ${corr.eval(p.t)}, chemicals ${want}'
				checked++
			}
			if want := p.chemicals_dy {
				assert agree(corr.deriv(p.t), want, 1e-9, math.abs(p.y) / p.t), '${what}: dY/dT = ${corr.deriv(p.t)}, chemicals ${want}'
			}
		}
	}
	assert checked >= 100
}
