// Module quality sweeps the component data for problems and writes
// docs/data-quality.md (M2 exit criteria 5 and 6).
//
// Nothing here fails the build: the data is ChemSep's, and its defects are
// facts to report, not bugs to fix. A test regenerates the report and fails
// if it differs from the committed one, so a new problem can't appear
// without someone seeing it in a diff.
module quality

import math
import os
import components
import components.chemsep
import thermo

// Where the report goes, relative to the repository root.
pub const report_path = 'docs/data-quality.md'

// Second virial coefficients may have either sign. Heats of vaporization
// and surface tensions vanish at Tc, so they may be zero. Every other
// property must be positive.
const may_be_negative = ['second_virial']
const may_be_zero = ['heat_of_vaporization', 'surface_tension']

// How many temperatures, evenly spaced over a correlation's range, are
// checked.
const sweep_points = 11

// 20 representative compounds for comparing PR's vapour pressure with
// ChemSep's correlation: light and heavy alkanes, an iso-alkane, a
// cycloalkane, aromatics, alkenes, inorganic gases, and polar compounds,
// where PR with the Soave α is expected to do worst.
const representative = ['74-82-8', '74-84-0', '74-98-6', '106-97-8', '75-28-5', '109-66-0', '110-54-3',
	'142-82-5', '124-18-5', '110-82-7', '71-43-2', '108-88-3', '74-85-1', '115-07-1', '124-38-9',
	'7727-37-9', '7783-06-4', '7732-18-5', '67-56-1', '67-64-1']
const reduced_temperatures = [0.5, 0.65, 0.8, 0.95]

// A Problem is one thing wrong with one compound's data.
pub struct Problem {
pub:
	cas  string
	name string
	what string
}

// BoilingPoint is a vapour pressure at the normal boiling point that isn't
// within 2% of 1 atm.
pub struct BoilingPoint {
pub:
	cas       string
	name      string
	key       string // the correlation used
	tb        f64
	p         f64
	deviation f64 // P/101325 − 1
	in_range  bool
}

// PrPsat compares PR's vapour pressure with ChemSep's correlation.
pub struct PrPsat {
pub:
	cas        string
	name       string
	key        string
	deviations []f64  // PR/ChemSep − 1, at each of reduced_temperatures
	in_range   []bool // whether ChemSep's correlation was in range there
}

// Sweep is everything the sweep found.
pub struct Sweep {
pub mut:
	compounds       int
	evaluated       int       // correlations evaluated
	properties      []Problem // not finite, or negative where they must be positive
	near_tc         int       // how many of those fail within 1% of Tc
	boiling_checked int
	boiling         []BoilingPoint
	pr              []PrPsat
	shared_inchikey []Problem
	past_tc         []Problem // eq.-106 ranges running past Tc
	no_ideal_gas    []Problem
	no_component    []Problem
	fallbacks       []Problem // vapour pressure from the Antoine record
}

// sweep checks every compound in the database.
pub fn sweep(db components.Database) !Sweep {
	mut s := Sweep{}
	keys := components.correlation_keys()
	mut by_inchikey := map[string][]components.IndexEntry{}
	for e in db.list() {
		r := db.find(e.cas)!
		s.compounds++
		if k := e.inchikey {
			by_inchikey[k] << e
		}
		for key in keys {
			p := r.property(key) or { continue } // unusable: the import report lists why
			if p.key != key {
				continue // the Antoine fallback, which is checked under its own key
			}
			s.evaluated++
			tc := r.constants.tc or { 0.0 }
			if f := check_property(key, p, tc) {
				s.properties << Problem{r.ids.cas, r.name, '${key}: ${f.why}'}
				if tc > 0.0 && math.abs(f.t / tc - 1.0) <= 0.01 {
					s.near_tc++
				}
			}
			d := r.correlations.get(key) or { continue }
			if d.eq == 106 && tc > 0.0 && d.tmax > tc {
				s.past_tc << Problem{r.ids.cas, r.name, '${key}: Tmax ${d.tmax} K, Tc ${tc} K'}
			}
		}
		r.ideal_gas() or { s.no_ideal_gas << Problem{r.ids.cas, r.name, err.msg()} }
		r.component() or { s.no_component << Problem{r.ids.cas, r.name, err.msg()} }
		vp := r.property('vapor_pressure') or { continue }
		if vp.key != 'vapor_pressure' {
			s.fallbacks << Problem{r.ids.cas, r.name, vp.note}
		}
		if tb := r.constants.tb {
			s.boiling_checked++
			v := vp.at(tb)
			dev := v.value / 101325.0 - 1.0
			if !(math.abs(dev) <= 0.02) {
				s.boiling << BoilingPoint{r.ids.cas, r.name, vp.key, tb, v.value, dev, v.in_range}
			}
		}
	}
	for k, list in by_inchikey {
		if list.len > 1 {
			s.shared_inchikey << Problem{list[0].cas, list[0].name, '${k}, also ${list[1..].map('${it.name} (${it.cas})').join(', ')}'}
		}
	}
	s.shared_inchikey.sort(a.cas < b.cas)
	for cas in representative {
		s.pr << pr_psat(db.find(cas)!)!
	}
	return s
}

// check_property evaluates a property across its range, and says what's
// wrong, or none. T = 0 K is skipped: some solid-phase ranges start there,
// and forms in 1/T aren't defined at it. tc, when known, puts the failing
// temperature in context: most failures sit at the critical point.
fn check_property(key string, p components.Property, tc f64) ?Failure {
	lo, hi := p.c.tmin, p.c.tmax
	n := if lo == hi { 1 } else { sweep_points }
	for i in 0 .. n {
		t := if n == 1 { lo } else { lo + (hi - lo) * f64(i) / f64(n - 1) }
		if t <= 0.0 {
			continue
		}
		v := p.at(t).value
		if !math.is_finite(v) {
			return Failure{'not finite at ${at(t, tc)}', t}
		}
		if key in may_be_negative {
			continue
		}
		if v < 0.0 {
			return Failure{'negative at ${at(t, tc)}: ${g4(v)} ${p.unit}', t}
		}
		if v == 0.0 && key !in may_be_zero {
			return Failure{'zero at ${at(t, tc)}', t}
		}
	}
	return none
}

// Failure is where a property first goes wrong.
struct Failure {
	why string
	t   f64
}

// at writes a temperature, with T/Tc when Tc is known.
fn at(t f64, tc f64) string {
	return if tc > 0.0 { '${g5(t)} K (T/Tc ${f3(t / tc)})' } else { '${g5(t)} K' }
}

// pr_psat compares PR's vapour pressure, from the compound's Tc, Pc and ω
// with the Soave α, against its vapour-pressure correlation.
fn pr_psat(r components.Record) !PrPsat {
	comp := r.component()!
	eos := thermo.new_cubic(.pr, [comp])!
	vp := r.property('vapor_pressure')!
	mut devs := []f64{}
	mut in_range := []bool{}
	for tr in reduced_temperatures {
		t := tr * comp.tc
		ps := thermo.saturation_pressure(eos, t)!
		v := vp.at(t)
		devs << ps / v.value - 1.0
		in_range << v.in_range
	}
	return PrPsat{r.ids.cas, r.name, vp.key, devs, in_range}
}

// generate renders the whole report, from the repository reflux was built
// from (@VMODROOT): the component files, the vendored XML and parameter
// files, and the CoolProp fixture.
pub fn generate() !string {
	db := components.load([components.default_dir])!
	s := sweep(db)!
	imp := chemsep.import_files(@VMODROOT, chemsep.vendored)!
	cool := read_coolprop(os.join_path(@VMODROOT, 'testdata', 'm2_coolprop.json'))!
	pr := components.load_kij(.pr)!
	srk := components.load_kij(.srk)!
	mc := components.load_mathias_copeman_pr()!
	return render(s, imp, cool, [pr, srk], mc, db)
}

// Number formatting that comes out the same whichever C compiler built
// reflux. The C library's printf breaks exact ties its own way: a tcc build
// printed 267.25 K, to four figures, as 267.3 and a gcc build as 267.2. So
// each value is rounded here first, half away from zero, with V's own
// math.round, and only then printed. A rounded value is never a tie.

// g4 and g5 print x to 4 and 5 significant figures.
fn g4(x f64) string {
	return '${round_sig(x, 4):.4g}'
}

fn g5(x f64) string {
	return '${round_sig(x, 5):.5g}'
}

// f2 and f3 print x to 2 and 3 decimal places.
fn f2(x f64) string {
	return '${round_to(x, 2):.2f}'
}

fn f3(x f64) string {
	return '${round_to(x, 3):.3f}'
}

// signed puts a + in front of a positive number, as printf's `+` flag does.
fn signed(s string) string {
	return if s.starts_with('-') { s } else { '+' + s }
}

fn round_sig(x f64, n int) f64 {
	if x == 0.0 || !math.is_finite(x) {
		return x
	}
	e := int(math.floor(math.log10(math.abs(x))))
	return round_to(x, n - 1 - e)
}

fn round_to(x f64, decimals int) f64 {
	s := math.pow(10.0, f64(decimals))
	return math.round(x * s) / s
}
