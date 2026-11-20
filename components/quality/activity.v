module quality

import activity
import components
import json2
import math
import os

// M6's activity-model data (step 1, criterion 2): the NRTL, UNIQUAC, Wilson
// and Henry's-law banks, and the three UNIFAC tables, each against thermo's
// independent copy of the published table (testdata/m6_unifac_tables.json,
// from tools/gen_m6_activity.py).

struct ThermoSubgroup {
	id   int
	name string
	main int
	r    f64 @[json: 'R']
	q    f64 @[json: 'Q']
}

struct ThermoInteraction {
	m   int
	n   int
	abc []f64
}

struct ThermoTable {
	subgroups    []ThermoSubgroup
	interactions []ThermoInteraction
}

struct ThermoTables {
	original      ThermoTable @[json: 'UNIFAC']
	lle           ThermoTable @[json: 'UNIFAC-LLE']
	dortmund      ThermoTable @[json: 'UNIFAC-Dortmund']
	dortmund_2006 ThermoTable @[json: 'UNIFAC-Dortmund-2006']
}

fn read_thermo_tables() !ThermoTables {
	return json2.decode[ThermoTables](os.read_file(os.join_path(@VMODROOT, 'testdata', 'm6_unifac_tables.json'))!)!
}

// thermo_number is ChemSep's subgroup id as thermo's table numbers it: in
// the original, ChemSep's C=C is 9, as in the 1987 Properties of Gases and
// Liquids, and thermo's 70, so ChemSep's 10 to 70 are thermo's 9 to 69
// (unifacrq.gct's first line). The LLE and Dortmund tables share numbers.
fn thermo_number(v activity.Variant, id int) int {
	if v != .original {
		return id
	}
	return if id == 9 {
		70
	} else if id >= 10 && id <= 70 {
		id - 1
	} else {
		id
	}
}

// sig_digits is how many significant digits x's shortest form has.
fn sig_digits(x f64) int {
	mut s := x.str().trim_left('-')
	if s.contains('e') {
		s = s.all_before('e')
	}
	digits := s.replace('.', '').trim_left('0')
	return math.max(digits.trim_right('0').len, 1)
}

// last_unit is one unit in the last significant digit of x as its
// shortest form prints it: 0.1 for 287.7, 1 for 727.0.
fn last_unit(x f64) f64 {
	if x == 0.0 {
		return 0.0
	}
	e := math.floor(math.log10(math.abs(x)))
	return math.pow(10.0, e - f64(sig_digits(x)) + 1.0)
}

// rounds_to is whether a and b are one number, printed to different
// precision: within a unit of the last digit of whichever is printed to
// fewer, which takes rounding (8484 for 8483.5) and cutting short (287.7 for
// 287.77) alike.
fn rounds_to(a f64, b f64) bool {
	if a == b {
		return true
	}
	unit := if sig_digits(a) <= sig_digits(b) { last_unit(a) } else { last_unit(b) }
	return math.abs(a - b) <= unit * (1.0 + 1e-9)
}

// all_round is whether parameters p and q are one, each of a, b and c
// within rounding of the other's, or too close to matter: moving ψ's
// exponent, (a + bT + cT²)/T, by under 1e-6 at 373.15 K (Dortmund's c of
// −7.4e-14 against −7.4e-16).
fn all_round(p []f64, q []f64) bool {
	t := 373.15
	for k in 0 .. p.len {
		d := math.abs(p[k] - q[k])
		effect := match k {
			0 { d / t }
			1 { d }
			else { d * t }
		}
		if !rounds_to(p[k], q[k]) && effect >= 1e-6 {
			return false
		}
	}
	return true
}

// Comparison is one of ChemSep's UNIFAC tables against thermo's.
pub struct Comparison {
pub mut:
	subgroups_same    int
	subgroups_rounded int
	subgroups_diff    []string
	subgroups_own     int // ChemSep's, thermo hasn't
	params_same       int
	params_rounded    int
	params_older      int // as thermo's older table has them: Dortmund's 2006
	params_diff       []string
	params_keys       []int // the differing parameters, by ChemSep's main groups
	params_own        int
}

// compare compares a table with thermo's, older being thermo's previous
// table where it has one. thermo's main groups are found through the
// subgroups each main group of ChemSep's has.
pub fn compare(v activity.Variant, t activity.Table, th ThermoTable, older ThermoTable) Comparison {
	mut c := Comparison{}
	mut tsub := map[int]ThermoSubgroup{}
	for s in th.subgroups {
		tsub[s.id] = s
	}
	mut tparam := map[int][]f64{}
	for p in th.interactions {
		tparam[activity.interaction_key(p.m, p.n)] = p.abc
	}
	mut oparam := map[int][]f64{}
	for p in older.interactions {
		oparam[activity.interaction_key(p.m, p.n)] = p.abc
	}
	mut mains := map[int]int{}
	mut ids := t.subgroups.keys()
	ids.sort()
	for id in ids {
		sg := t.subgroups[id] or { continue }
		ts := tsub[thermo_number(v, id)] or {
			c.subgroups_own++
			continue
		}
		mains[sg.main] = ts.main
		if ts.r == sg.r && ts.q == sg.q {
			c.subgroups_same++
		} else if rounds_to(ts.r, sg.r) && rounds_to(ts.q, sg.q) {
			c.subgroups_rounded++
		} else {
			c.subgroups_diff << '${sg.id} ${sg.name}: R ${sg.r}, Q ${sg.q}; thermo ${ts.r}, ${ts.q}'
		}
	}
	mut keys := t.interactions.keys()
	keys.sort()
	for key in keys {
		p := t.interactions[key] or { continue }
		m, n := key / 1000, key % 1000
		tm := mains[m] or {
			c.params_own++
			continue
		}
		tn := mains[n] or {
			c.params_own++
			continue
		}
		abc := tparam[activity.interaction_key(tm, tn)] or {
			c.params_own++
			continue
		}
		ours := [p.a, p.b, p.c]
		old := oparam[activity.interaction_key(tm, tn)] or { []f64{} }
		if ours == abc {
			c.params_same++
		} else if all_round(ours, abc) {
			c.params_rounded++
		} else if old.len == 3 && all_round(ours, old) {
			c.params_older++
		} else {
			c.params_diff << '${m} ${t.mains[m] or { '?' }} with ${n} ${t.mains[n] or { '?' }}: ${ours}; thermo ${abc}'
			c.params_keys << key
		}
	}
	return c
}

// ddbst_dortmund are the parameters of Dortmund's table, by ChemSep's main
// groups, that differ from thermo's 2016 copy and that DDBST's published
// page prints as ChemSep has them (https://www.ddbst.com/PublishedParametersUNIFACDO.html,
// read on 2026-11-20): so thermo's copy is the one that differs.
pub const ddbst_dortmund = [3027, 27003, 3029, 29003, 4012, 12004, 4027, 27004, 7011, 11007, 7013,
	13007]

// activity_lines is the report's section on the activity-model data.
fn activity_lines(db components.Database) ![]string {
	mut o := []string{}
	o << ''
	o << '## Activity-model parameters (M6)'
	o << ''
	o << "ChemSep's binary banks and Henry's-law constants, read as vendored:"
	o << ''
	o << '| File | Entries | Pairs | Listed more than once | Unreadable lines | Units |'
	o << '|---|---|---|---|---|---|'
	for kind in [components.ActivityBank.nrtl, .uniquac, .wilson] {
		b := components.load_pair_bank(db.vendor, kind)!
		mut pairs := map[string]int{}
		for e in b.entries {
			key := if e.cas1 < e.cas2 { '${e.cas1}/${e.cas2}' } else { '${e.cas2}/${e.cas1}' }
			pairs[key]++
		}
		o << '| `${b.file}` | ${b.entries.len} | ${pairs.len} | ${pairs.values().filter(it > 1).len} | ${b.defects.len} | ${b.units} |'
	}
	henry := components.load_henry(db.vendor)!
	o << '| `${henry.file}` | ${henry.entries.len} | ${henry.entries.len} gases, in water | 0 | ${henry.defects.len} | Pa |'
	o << ''
	o << 'Where a pair is listed more than once, reflux uses the first listing, as'
	o << "for kᵢⱼ. The banks' A are in cal/mol over DECHEMA's R, 1.98721 cal/(mol·K)"
	o << '(M6 step 0).'
	o << ''
	th := read_thermo_tables()!
	o << '### UNIFAC tables'
	o << ''
	o << "Each table is ChemSep's base file with its extension. Where a file"
	o << 'contradicts itself, or the published table shows a typo, the reader'
	o << 'repairs it and says so (M6 step 1):'
	o << ''
	for v in [activity.Variant.original, .lle, .dortmund] {
		load := components.load_unifac(db.vendor, v)!
		t := load.table
		o << '- **${v.name()}** (${load.files.map('`${it}`').join(', ')}): ${t.subgroups.len} subgroups, ${t.mains.len} main groups, ${t.interactions.len} parameters.'
		for r in load.repairs {
			place := if r.line > 0 { ' line ${r.line}' } else { '' }
			o << '  - Repaired, `${r.file}`${place}: ${md(r.before)} → ${md(r.after)}, ${md(r.why)}.'
		}
		for d in load.defects {
			o << '  - `${d.file}` line ${d.line}: ${md(d.why)}.'
		}
	}
	o << ''
	o << "Compounds whose groups aren't all in their table, so that the table"
	o << "can't take them:"
	o << ''
	for v in [activity.Variant.original, .lle, .dortmund] {
		t := components.load_unifac(db.vendor, v)!.table
		mut missing := []string{}
		for e in db.list() {
			r := db.find(e.cas)!
			g := match v {
				.original { r.groups.unifac_vle }
				.lle { r.groups.unifac_lle }
				.dortmund { r.groups.modified_unifac }
			}
			absent := g.keys().filter(it !in t.subgroups)
			if absent.len > 0 {
				missing << '${r.name} (${e.cas}), subgroups ${absent.map(it.str()).join(', ')}'
			}
		}
		o << '- ${v.name()}: ${if missing.len == 0 { 'none' } else { missing.join('; ') }}.'
	}
	o << ''
	o << "### Against thermo's copies"
	o << ''
	o << "thermo 0.6.1's tables are an independent copy of the published ones:"
	o << "the original (UFSG, UFIP), the LLE table and Dortmund's 2016, with its"
	o << '2006 parameters beside them. A value is *rounded* when the two are within'
	o << 'a unit of the last digit of whichever is printed to fewer: 8484 for'
	o << '8483.5, 287.7 for 287.77.'
	o << ''
	o << "| Table | Subgroups same | Rounded | Different | ChemSep's own | Parameters same | Rounded | As the 2006 table | Different | ChemSep's own |"
	o << '|---|---|---|---|---|---|---|---|---|---|'
	for v, pair in {
		activity.Variant.original: [th.original, ThermoTable{}]
		.lle:                      [th.lle, ThermoTable{}]
		.dortmund:                 [th.dortmund, th.dortmund_2006]
	} {
		c := compare(v, components.load_unifac(db.vendor, v)!.table, pair[0], pair[1])
		o << '| ${v.name()} | ${c.subgroups_same} | ${c.subgroups_rounded} | ${c.subgroups_diff.len} | ${c.subgroups_own} | ${c.params_same} | ${c.params_rounded} | ${c.params_older} | ${c.params_diff.len} | ${c.params_own} |'
	}
	o << ''
	o << 'Each difference, and what settles it:'
	o << ''
	for v, pair in {
		activity.Variant.original: [th.original, ThermoTable{}]
		.lle:                      [th.lle, ThermoTable{}]
		.dortmund:                 [th.dortmund, th.dortmund_2006]
	} {
		c := compare(v, components.load_unifac(db.vendor, v)!.table, pair[0], pair[1])
		for s in c.subgroups_diff {
			why := if v == .lle && s.starts_with('30 ') {
				"thermo's typo, 9183 for the original table's 0.9183; ChemSep's stands"
			} else {
				'unsettled: no published source for these was read'
			}
			o << '- ${v.name()}, subgroup ${md(s)}: ${why}.'
		}
		for k, s in c.params_diff {
			why := if v == .dortmund && c.params_keys[k] in ddbst_dortmund {
				"DDBST's published page prints ChemSep's values; thermo's copy is the one that differs"
			} else {
				'unsettled: no published source for these was read'
			}
			o << '- ${v.name()}, ${md(s)}: ${why}.'
		}
	}
	return o
}
