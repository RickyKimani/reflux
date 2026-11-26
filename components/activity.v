module components

import activity
import math
import os

// ChemSep's activity-model parameters (M6 step 1, D7): the binary banks
// nrtl.ipd, uniquac.ipd and wilson.ipd, Henry's-law constants in
// henrys_law.ipd, and the three UNIFAC tables, each with its extension.
// Read straight from data/vendor/chemsep/, as the kᵢⱼ banks are, every line
// that can't be read kept as a defect.
//
// The banks give A in cal/mol and state no formula and no R. DECHEMA's
// forms are NRTL's τ = A/RT, UNIQUAC's τ = exp(−A/RT) and Wilson's
// Λ = (Vⱼ/Vᵢ)·exp(−A/RT), with R = 1.98721 cal/(mol·K). That's confirmed
// for NRTL alone (step 0: DWSIM's, on the same bank, to 4e-10). For UNIQUAC
// and Wilson, whose bank is Andrew Mason's compilation, that R is inferred,
// and so is Wilson's form, until step 12. A bank's parameter is written as
// the expression with that R, the calorie cancelling (D7 as amended):
// `-189.0469 K / (1.98721*T)`.

// dechema_r is DECHEMA's R, in cal/(mol·K), as text.
pub const dechema_r = '1.98721'

// PairParams is one line of a binary bank: its two compounds, its numbers
// as values and as the file prints them, and the rest of the line.
pub struct PairParams {
pub:
	cas1    string
	cas2    string
	values  []f64
	texts   []string // each number as printed
	comment string   // e.g. 'Methanol/Water p61 1/1a' (the DECHEMA volume and page)
	line    int
}

// PairBank is a binary bank: NRTL's A₁₂, A₂₁ and α, or UNIQUAC's or
// Wilson's A₁₂ and A₂₁.
pub struct PairBank {
pub:
	file    string
	comment string
	units   string // its Units= line: cal/mol
	entries []PairParams
	defects []Defect
}

// ActivityBank names a binary bank.
pub enum ActivityBank {
	nrtl
	uniquac
	wilson
}

fn (k ActivityBank) file() string {
	return match k {
		.nrtl { 'nrtl.ipd' }
		.uniquac { 'uniquac.ipd' }
		.wilson { 'wilson.ipd' }
	}
}

fn (k ActivityBank) columns() int {
	return if k == .nrtl { 3 } else { 2 }
}

// load_pair_bank loads one of ChemSep's binary banks from a directory of its
// files (vendor_dir, or a database's vendor).
pub fn load_pair_bank(dir string, kind ActivityBank) !PairBank {
	path := os.join_path(dir, kind.file())
	return parse_pair_bank(os.read_file(path)!, kind.file(), kind.columns())!
}

// parse_pair_bank reads a binary bank's text: each data line is `CAS1 CAS2`,
// then columns numbers, then a comment.
pub fn parse_pair_bank(text string, file string, columns int) !PairBank {
	comment, units, lines := lib_section(text, 'IPD') or { return error('${file}: ${err}') }
	mut entries := []PairParams{}
	mut defects := []Defect{}
	for l in lines {
		f, rest := leading_fields(l.text, 2 + columns)
		if f.len < 2 + columns {
			defects << Defect{l.n, l.text, 'fewer than ${2 + columns} fields'}
			continue
		}
		check_cas(f[0]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		check_cas(f[1]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		if f[0] == f[1] {
			defects << Defect{l.n, l.text, 'a compound paired with itself'}
			continue
		}
		mut values := []f64{}
		for s in f[2..] {
			values << finite_number(s) or {
				defects << Defect{l.n, l.text, err.msg()}
				break
			}
		}
		if values.len < columns {
			continue
		}
		entries << PairParams{
			cas1:    f[0]
			cas2:    f[1]
			values:  values
			texts:   f[2..].map(number_text(it))
			comment: rest
			line:    l.n
		}
	}
	return PairBank{
		file:    file
		comment: comment
		units:   units
		entries: entries
		defects: defects
	}
}

// number_text is a number as printed, with a 0 before a leading point, which
// expr needs: `.187e-1` as `0.187e-1`, the same double.
fn number_text(s string) string {
	if s.starts_with('.') {
		return '0' + s
	}
	if s.starts_with('-.') {
		return '-0' + s[1..]
	}
	return s
}

// for_pair returns every listing of the pair a, b, in file order, whichever
// way round each line names it.
pub fn (b PairBank) for_pair(a string, c string) []PairParams {
	return b.entries.filter((it.cas1 == a && it.cas2 == c) || (it.cas1 == c && it.cas2 == a))
}

// HenryParams is one line of henrys_law.ipd: a gas's Henry's-law constant
// in a solvent, by ChemSep's equation eq over [tmin, tmax], in Pa.
pub struct HenryParams {
pub:
	gas     string
	solvent string
	b       int // the file's `b` column, 1 throughout, undocumented
	eq      int
	tmin    f64
	tmax    f64
	coeffs  []f64 // A to E
	comment string
	line    int
}

// HenryBank is henrys_law.ipd.
pub struct HenryBank {
pub:
	file    string
	comment string
	entries []HenryParams
	defects []Defect
}

// load_henry loads ChemSep's Henry's-law constants from a directory of its
// files.
pub fn load_henry(dir string) !HenryBank {
	return parse_henry(os.read_file(os.join_path(dir, 'henrys_law.ipd'))!, 'henrys_law.ipd')!
}

// parse_henry reads henrys_law.ipd's text: each data line is `gas solvent b
// eq tmin tmax A B C D E`, then the gas's name.
pub fn parse_henry(text string, file string) !HenryBank {
	comment, _, lines := lib_section(text, 'IPD') or { return error('${file}: ${err}') }
	mut entries := []HenryParams{}
	mut defects := []Defect{}
	for l in lines {
		f, rest := leading_fields(l.text, 11)
		if f.len < 11 {
			defects << Defect{l.n, l.text, 'fewer than eleven fields'}
			continue
		}
		check_cas(f[0]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		check_cas(f[1]) or {
			defects << Defect{l.n, l.text, err.msg()}
			continue
		}
		mut nums := []f64{}
		for s in f[2..] {
			nums << finite_number(s) or {
				defects << Defect{l.n, l.text, err.msg()}
				break
			}
		}
		if nums.len < 9 {
			continue
		}
		if nums[0] != f64(int(nums[0])) || nums[1] != f64(int(nums[1])) {
			defects << Defect{l.n, l.text, 'b and the equation must be whole numbers'}
			continue
		}
		if !(nums[2] < nums[3]) {
			defects << Defect{l.n, l.text, 'an empty temperature range'}
			continue
		}
		entries << HenryParams{
			gas:     f[0]
			solvent: f[1]
			b:       int(nums[0])
			eq:      int(nums[1])
			tmin:    nums[2]
			tmax:    nums[3]
			coeffs:  nums[4..9].clone()
			comment: rest
			line:    l.n
		}
	}
	return HenryBank{
		file:    file
		comment: comment
		entries: entries
		defects: defects
	}
}

// TableDefect is a line of a UNIFAC table that couldn't be read, or that
// repeats one before it.
pub struct TableDefect {
pub:
	file string
	line int // 1-based
	text string
	why  string
}

// TableRepair is a value of a UNIFAC table changed on reading, because the
// file contradicts itself there, and why.
pub struct TableRepair {
pub:
	file   string
	line   int
	text   string
	before string
	after  string
	why    string
}

// UnifacLoad is a UNIFAC table as read, with what couldn't be, and what was
// repaired.
pub struct UnifacLoad {
pub:
	table   activity.Table
	files   []string
	defects []TableDefect
	repairs []TableRepair
}

// UnifacFiles are a variant's files: its main groups' and parameters' (one
// for each of a, b and c in Dortmund's), then its subgroups', each followed
// by its extension.
struct UnifacFiles {
	params [][]string
	subs   []string
}

fn unifac_files(variant activity.Variant) UnifacFiles {
	return match variant {
		.original {
			UnifacFiles{[['unifacvl.gct', 'unifacvl_ext.gct']], ['unifacrq.gct', 'unifacrq_ext.gct']}
		}
		.lle {
			UnifacFiles{[['unifacll.gct', 'unifacll_ext.gct']], ['llefacrq.gct', 'llefacrq_ext.gct']}
		}
		.dortmund {
			UnifacFiles{[['modfac_a.gct', 'modfac_a_ext.gct'], ['modfac_b.gct', 'modfac_b_ext.gct'],
				['modfac_c.gct', 'modfac_c_ext.gct']], [
				'modfacrq.gct',
				'modfacrq_ext.gct',
			]}
		}
	}
}

// load_unifac loads a UNIFAC table from a directory of ChemSep's files, its
// extension's groups and parameters added to the base's. A group or a
// parameter listed twice keeps its first listing, the other reported.
//
// The files contradict themselves in places (M6 step 1): a subgroup or a
// parameter row names a main group whose number, in the list, is another's.
// Where the name matches another main group's in the list exactly, case
// aside, and neither name begins the other (spellings such as Morph and
// Morpholine, or PYRIDIN and Pyridine, aren't contradictions), the number is
// taken from the name, and the repair reported. A number a row uses that the
// list doesn't have is matched to the list by its name likewise, or else is
// a new main group, named by its rows.
pub fn load_unifac(dir string, variant activity.Variant) !UnifacLoad {
	uf := unifac_files(variant)
	mut mains := map[int]string{}
	mut defects := []TableDefect{}
	mut repairs := []TableRepair{}
	mut files := []string{}
	mut all_rows := [][]TableRow{len: uf.params.len}
	mut row_files := [][]string{len: uf.params.len}
	for coef, pair in uf.params {
		for name in pair {
			files << name
			text := os.read_file(os.join_path(dir, name))!
			ms, rows, ds := parse_main_groups(text, name)
			defects << ds
			for id, mname in ms {
				if id !in mains {
					mains[id] = mname
				}
			}
			for row in rows {
				all_rows[coef] << row
				row_files[coef] << name
			}
		}
	}
	// The a rows' pairs by their numbers, to tell a row's own pair from one
	// another row lists; then each row's pair resolved, the a rows first,
	// whose resolved pairs decide the b and c rows'.
	mut numbered := map[int]int{}
	for row in all_rows[0] {
		numbered[mains_key(row.m, row.n)]++
	}
	mut resolved_a := map[int]bool{}
	mut inter := map[int]activity.Interaction{}
	for coef, rows in all_rows {
		for k, row in rows {
			file := row_files[coef][k]
			m, n := resolve_row(mut mains, mut repairs, mut defects, row, file, coef, numbered,
				resolved_a) or { continue }
			// A pair listed before, in this coefficient's rows, keeps its
			// first listing, the row reported once.
			if coef == 0 {
				if resolved_a[mains_key(m, n)] {
					defects << TableDefect{file, row.line, row.text, 'a pair listed again'}
					continue
				}
				resolved_a[mains_key(m, n)] = true
			} else {
				old := inter[activity.interaction_key(m, n)] or { activity.Interaction{} }
				if (coef == 1 && old.b != 0.0) || (coef == 2 && old.c != 0.0) {
					defects << TableDefect{file, row.line, row.text, 'a pair listed again'}
					continue
				}
			}
			for key, v in {
				activity.interaction_key(m, n): row.mn
				activity.interaction_key(n, m): row.nm
			} {
				old := inter[key] or { activity.Interaction{} }
				inter[key] = match coef {
					0 {
						activity.Interaction{
							...old
							a: v
						}
					}
					1 {
						activity.Interaction{
							...old
							b: v
						}
					}
					else {
						activity.Interaction{
							...old
							c: v
						}
					}
				}
			}
		}
	}
	mut subgroups := map[int]activity.Subgroup{}
	for name in uf.subs {
		files << name
		text := os.read_file(os.join_path(dir, name))!
		lines, ds, rs := parse_subgroups(text, name)
		defects << ds
		repairs << rs
		for l in lines {
			mut sg := l.sg
			if sg.id in subgroups {
				defects << TableDefect{name, l.line, l.text, 'a subgroup listed again'}
				continue
			}
			own := mains[sg.main] or { '' }
			if contradicts(own, l.main_name) {
				named := named_main(mains, l.main_name)
				if named > 0 {
					repairs << TableRepair{name, l.line, l.text, 'main group ${sg.main}', 'main group ${named}', 'its main group is named ${l.main_name}, which is ${named} in the list; ${sg.main} is ${own}'}
					sg = activity.Subgroup{
						...sg
						main: named
					}
				}
			}
			subgroups[sg.id] = sg
		}
	}
	for c in unifac_corrections {
		if c.variant != variant {
			continue
		}
		if c.subgroup > 0 {
			sg := subgroups[c.subgroup] or { continue }
			shipped := if c.field == 'R' { sg.r } else { sg.q }
			if shipped != c.shipped {
				defects << TableDefect{uf.subs[0], 0, '${sg.id} ${sg.name}', "a correction's shipped ${c.field}, ${c.shipped}, isn't the file's ${shipped}: not applied"}
				continue
			}
			subgroups[c.subgroup] = if c.field == 'R' {
				activity.Subgroup{
					...sg
					r: c.published
				}
			} else {
				activity.Subgroup{
					...sg
					q: c.published
				}
			}
			repairs << TableRepair{uf.subs[0], 0, '${sg.id} ${sg.name}', '${c.field} ${c.shipped}', '${c.field} ${c.published}', 'the published table: ${c.source}'}
			continue
		}
		key := activity.interaction_key(c.m, c.n)
		p := inter[key] or { continue }
		if p.a != c.shipped {
			defects << TableDefect{uf.params[0][0], 0, '${c.m} with ${c.n}', "a correction's shipped a, ${c.shipped}, isn't the file's ${p.a}: not applied"}
			continue
		}
		inter[key] = activity.Interaction{
			...p
			a: c.published
		}
		repairs << TableRepair{uf.params[0][0], 0, '${c.m} ${mains[c.m] or { '?' }} with ${c.n} ${mains[c.n] or {
			'?'
		}}', 'a ${c.shipped}', 'a ${c.published}', 'the published table: ${c.source}'}
	}
	return UnifacLoad{
		table:   activity.Table{
			variant:      variant
			subgroups:    subgroups
			mains:        mains
			interactions: inter
		}
		files:   files
		defects: defects
		repairs: repairs
	}
}

// UnifacCorrection is a value of ChemSep's UNIFAC tables that the published
// table contradicts (M6 step 1, D7 as amended): a subgroup's R or Q, or a
// parameter aₘₙ, by ChemSep's numbers, with the value as shipped and as
// published. A correction applies only where the file still holds the
// shipped value, and each is reported.
struct UnifacCorrection {
	variant   activity.Variant
	subgroup  int    // a subgroup's id, or 0 for a parameter
	field     string // R or Q
	m         int
	n         int
	shipped   f64
	published f64
	source    string
}

const ddbst_unifac = "DDBST's published UNIFAC parameters, https://www.ddbst.com/published-parameters-unifac.html, read on 2026-11-20, which thermo 0.6.1's copy agrees with"

const ddbst_groups = "R and Q are a group's own, the same in the LLE table: DDBST's published UNIFAC parameters and thermo 0.6.1's copy of Magnussen et al.'s table"

const ddbst_dortmund = "DDBST's published modified UNIFAC (Dortmund) parameters, https://www.ddbst.com/PublishedParametersUNIFACDO.html, read on 2026-11-20, which thermo 0.6.1's copy agrees with"

// The original table's typos (a digit dropped or changed) and its OH–Br
// pair, which isn't the published one at all; the LLE table's copies of the
// two subgroups; and Dortmund's quaternary C, whose published Q is 0, and
// six ionic-liquid subgroups, whose R and Q aren't the published ones.
// ChemSep's Dortmund parameters that differ from thermo's are the published
// ones, so none is corrected. The differences of rounding (ChemSep printing
// 8484 for 8483.5) are left as shipped.
const unifac_corrections = [
	UnifacCorrection{.original, 47, 'R', 0, 0, 1.006, 1.0106, ddbst_unifac},
	UnifacCorrection{.original, 62, 'Q', 0, 0, 2.481, 2.484, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 5, 33, 133.4, 202.1, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 33, 5, 721.9, 524.9, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 17, 24, 822.0, 882.0, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 33, 27, 4228.0, 4288.0, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 5, 34, 727.0, 727.8, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 4, 10, 586.6, 586.8, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 11, 24, 54.47, 54.57, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 32, 28, -27.545, -27.45, ddbst_unifac},
	UnifacCorrection{.original, 0, '', 38, 13, 167.1, 167.3, ddbst_unifac},
	UnifacCorrection{.lle, 22, 'Q', 0, 0, 2.481, 2.484, ddbst_groups},
	UnifacCorrection{.lle, 33, 'R', 0, 0, 1.006, 1.0106, ddbst_groups},
	UnifacCorrection{.dortmund, 4, 'Q', 0, 0, 0.4582, 0.0, ddbst_dortmund},
	UnifacCorrection{.dortmund, 178, 'R', 0, 0, 1.7989, 1.3662, ddbst_dortmund},
	UnifacCorrection{.dortmund, 178, 'Q', 0, 0, 0.64, 0.6797, ddbst_dortmund},
	UnifacCorrection{.dortmund, 184, 'R', 0, 0, 2.411, 1.843, ddbst_dortmund},
	UnifacCorrection{.dortmund, 184, 'Q', 0, 0, 2.409, 1.6997, ddbst_dortmund},
	UnifacCorrection{.dortmund, 189, 'R', 0, 0, 2.7986, 2.7867, ddbst_dortmund},
	UnifacCorrection{.dortmund, 189, 'Q', 0, 0, 2.7744, 2.7723, ddbst_dortmund},
	UnifacCorrection{.dortmund, 195, 'R', 0, 0, 4.62, 3.9628, ddbst_dortmund},
	UnifacCorrection{.dortmund, 195, 'Q', 0, 0, 1.1707, 0.6214, ddbst_dortmund},
	UnifacCorrection{.dortmund, 196, 'R', 0, 0, 2.4878, 2.1094, ddbst_dortmund},
	UnifacCorrection{.dortmund, 196, 'Q', 0, 0, 2.474, 2.5106, ddbst_dortmund},
	UnifacCorrection{.dortmund, 197, 'R', 0, 0, 3.3854, 3.371, ddbst_dortmund},
	UnifacCorrection{.dortmund, 197, 'Q', 0, 0, 2.009, 2.0001, ddbst_dortmund},
]

// contradicts is whether a name a line gives a main group can't be the
// list's name for its number: they differ, case aside, and neither begins
// the other. An empty list name (no such number) contradicts nothing here.
fn contradicts(listed string, stated string) bool {
	a, b := listed.to_lower(), stated.to_lower()
	return a != '' && a != b && !a.starts_with(b) && !b.starts_with(a)
}

// named_main is the one main group the list names name, case aside, or 0.
fn named_main(mains map[int]string, name string) int {
	ids := mains.keys().filter((mains[it] or { '' }).to_lower() == name.to_lower())
	return if ids.len == 1 { ids[0] } else { 0 }
}

// mains_key is an unordered pair of main groups' key.
fn mains_key(m int, n int) int {
	return if m < n { m * 1000 + n } else { n * 1000 + m }
}

// resolve_row is the pair of main groups a parameter row means. By its
// numbers, unless a name contradicts its number's in the list, or names a
// number the list hasn't; then by its names where the evidence says so:
// - an a row (or the only table's row) takes its names' pair when its
//   numbers' pair is listed by another row too, or a number isn't in the
//   list at all;
// - a b or c row takes whichever pair the a rows have, resolved: Dortmund's
//   c rows have names copied from the row above, and the a and b rows the
//   right ones.
// A number neither the list nor a name accounts for is a new main group,
// named by its rows. A contradiction the evidence doesn't settle is
// reported, and the row skipped.
fn resolve_row(mut mains map[int]string, mut repairs []TableRepair, mut defects []TableDefect, row TableRow, file string, coef int, numbered map[int]int, resolved_a map[int]bool) ?(int, int) {
	mut ends := [row.m, row.n]
	mut by_name := [row.m, row.n]
	mut odd := false
	for k, name in [row.m_name, row.n_name] {
		listed := mains[ends[k]] or { '' }
		if listed == '' || contradicts(listed, name) {
			named := named_main(mains, name)
			if named > 0 {
				by_name[k] = named
				odd = true
			}
		}
	}
	if !odd {
		for k, name in [row.m_name, row.n_name] {
			if ends[k] !in mains {
				mains[ends[k]] = name
				repairs << TableRepair{file, row.line, row.text, 'no main group ${ends[k]}', 'main group ${ends[k]}, ${name}', "a main group the list doesn't number, named by its parameters' rows"}
			}
		}
		if ends[0] == ends[1] {
			defects << TableDefect{file, row.line, row.text, 'not a pair of main groups'}
			return none
		}
		return ends[0], ends[1]
	}
	numbers := mains_key(ends[0], ends[1])
	names := mains_key(by_name[0], by_name[1])
	missing := ends[0] !in mains || ends[1] !in mains
	take_names := if coef == 0 {
		missing || (numbered[numbers] or { 0 }) > 1
	} else {
		!resolved_a[numbers] && resolved_a[names]
	}
	if take_names {
		repairs << TableRepair{file, row.line, row.text, 'main groups ${ends[0]} and ${ends[1]}', 'main groups ${by_name[0]} and ${by_name[1]}', 'its names, ${row.m_name} and ${row.n_name}, are those numbers in the list'}
		return by_name[0], by_name[1]
	}
	if !missing && (coef == 0 || resolved_a[numbers]) {
		defects << TableDefect{file, row.line, row.text, "a name that isn't its number's in the list; the numbers stand, as ${if coef == 0 {
			'no other row lists their pair'
		} else {
			'the a rows pair them'
		}}"}
		return ends[0], ends[1]
	}
	defects << TableDefect{file, row.line, row.text, "a name that isn't its number's, and nothing to settle which is right"}
	return none
}

// TableRow is one interaction line: main groups m and n, aₘₙ and aₙₘ.
pub struct TableRow {
pub:
	m      int
	m_name string
	n      int
	n_name string
	mn     f64
	nm     f64
	line   int
	text   string
}

// parse_main_groups reads a main-group table: a count, that many `id name`
// lines, then `m name n name aₘₙ aₙₘ` lines. A `#` starts a comment.
pub fn parse_main_groups(text string, file string) (map[int]string, []TableRow, []TableDefect) {
	mut mains := map[int]string{}
	mut rows := []TableRow{}
	mut defects := []TableDefect{}
	for i, raw in text.split_into_lines() {
		if i == 0 {
			continue // the count, and what the file is
		}
		l := raw.all_before('#').trim_space()
		if l == '' {
			continue
		}
		f, _ := leading_fields(l, 6)
		match f.len {
			2 {
				id := f[0].int()
				if id <= 0 || f[0] != id.str() {
					defects << TableDefect{file, i + 1, raw, 'a main group with no number'}
					continue
				}
				mains[id] = f[1]
			}
			6 {
				m, n := f[0].int(), f[2].int()
				mn := finite_number(f[4]) or {
					defects << TableDefect{file, i + 1, raw, err.msg()}
					continue
				}
				nm := finite_number(f[5]) or {
					defects << TableDefect{file, i + 1, raw, err.msg()}
					continue
				}
				if m <= 0 || n <= 0 || m == n {
					defects << TableDefect{file, i + 1, raw, 'not a pair of main groups'}
					continue
				}
				rows << TableRow{m, f[1], n, f[3], mn, nm, i + 1, raw}
			}
			else {
				defects << TableDefect{file, i + 1, raw, 'neither a main group nor a pair'}
			}
		}
	}
	return mains, rows, defects
}

// SubgroupLine is a subgroup as a table's line gives it, with the name it
// gives its main group.
pub struct SubgroupLine {
pub:
	sg        activity.Subgroup
	main_name string
	line      int
	text      string
}

// parse_subgroups reads a subgroup table: a count, then `id name main
// mainname R Q M` lines.
pub fn parse_subgroups(text string, file string) ([]SubgroupLine, []TableDefect, []TableRepair) {
	mut out := []SubgroupLine{}
	mut defects := []TableDefect{}
	mut repairs := []TableRepair{}
	for i, raw in text.split_into_lines() {
		if i == 0 {
			continue
		}
		l := raw.all_before('#').trim_space()
		if l == '' {
			continue
		}
		f, _ := leading_fields(l, 7)
		if f.len < 7 {
			defects << TableDefect{file, i + 1, raw, 'fewer than seven fields'}
			continue
		}
		mut nums := []f64{}
		for k, s in f[4..7] {
			// The molar mass, which nothing uses, with a decimal comma
			// (unifacrq_ext.gct's `57,028`) is read as the point it means.
			if k == 2 && s.count(',') == 1 && !s.contains('.') {
				repairs << TableRepair{file, i + 1, raw, s, s.replace(',', '.'), 'a molar mass with a decimal comma'}
				nums << finite_number(s.replace(',', '.')) or { math.nan() }
				continue
			}
			nums << finite_number(s) or {
				defects << TableDefect{file, i + 1, raw, err.msg()}
				break
			}
		}
		id, main := f[0].int(), f[2].int()
		if nums.len < 3 || id <= 0 || main <= 0 {
			if nums.len == 3 {
				defects << TableDefect{file, i + 1, raw, 'no subgroup or main group number'}
			}
			continue
		}
		out << SubgroupLine{
			sg:        activity.Subgroup{
				id:   id
				name: f[1]
				main: main
				r:    nums[0]
				q:    nums[1]
				mw:   nums[2]
			}
			main_name: f[3]
			line:      i + 1
			text:      raw
		}
	}
	return out, defects, repairs
}

// PairSource says where a binary model's pair came from.
pub enum PairSource {
	bank        // its bank's first listing of the pair
	file        // the caller's, a file's [[model.pairs]] (D12)
	ideal       // none listed, by the caller's leave: NRTL's τ = 0 or Wilson's Λ = 1, an ideal mixture of the two
	no_residual // none listed, by the caller's leave: UNIQUAC's τ = 1, no residual interaction. r and q still count, so it isn't ideal
}

// PairUse records how a binary model got one pair, i before j.
pub struct PairUse {
pub:
	i       int
	j       int
	source  PairSource
	listing ?PairParams  // the bank's line used
	others  []PairParams // its other listings of the pair, unused
	texts   []string     // the parameters as expressions, i–j first
}

// PairOverride is a pair's parameters as a file gives them (D12), in place
// of the bank's: a and b by CAS number, and each parameter an expression in
// T, a's with b first: NRTL's τ_ab, τ_ba and α, UNIQUAC's τ_ab and τ_ba, or
// Wilson's Λ_ab and Λ_ba.
pub struct PairOverride {
pub:
	a     string
	b     string
	texts []string
}

@[params]
pub struct ActivityOptions {
pub:
	// ideal lets a pair the bank doesn't list have no energy of
	// interaction, where by default it's refused (D7). That's an ideal pair
	// in NRTL and Wilson (PairSource.ideal). UNIQUAC's keeps its
	// combinatorial term (PairSource.no_residual).
	ideal bool
	// overrides are pairs the caller gives, each in place of the bank's.
	overrides []PairOverride
}

// BinaryModel is a binary model built from a bank, with how it got each
// pair.
pub struct BinaryModel[M] {
pub:
	model M
	pairs []PairUse
}

// names are the compounds' names, and their records, by CAS number.
fn (db Database) records(cas []string) ![]Record {
	mut out := []Record{}
	for c in cas {
		out << db.find(c)!
	}
	return out
}

// pair_uses gives each pair's expressions, i–j orientation first: the
// caller's override as given, or the bank's listing by its formula. An
// unlisted pair is refused unless opts.ideal, when its source is unlisted
// and its texts none_listed. swapped turns an override's texts from b–a to
// a–b, the parameters' count being the model's.
fn pair_uses(bank PairBank, cas []string, names []string, opts ActivityOptions, texts fn (int, int, []string) []string, unlisted PairSource, none_listed []string, swapped fn ([]string) []string) ![]PairUse {
	mut out := []PairUse{}
	mut missing := []string{}
	for o in opts.overrides {
		ia, ib := cas.index(o.a), cas.index(o.b)
		if ia < 0 || ib < 0 || ia == ib {
			return error("a pair ${o.a} with ${o.b}, which isn't two of the compounds")
		}
		if o.texts.len != none_listed.len {
			return error('a pair of ${names[ia]} with ${names[ib]} has ${o.texts.len} parameters, where the model takes ${none_listed.len}')
		}
		if opts.overrides.filter((it.a == o.a && it.b == o.b) || (it.a == o.b && it.b == o.a)).len > 1 {
			return error('the pair ${names[ia]} with ${names[ib]} is given twice')
		}
	}
	for i in 0 .. cas.len {
		for j in i + 1 .. cas.len {
			listed := bank.for_pair(cas[i], cas[j])
			given := opts.overrides.filter((it.a == cas[i] && it.b == cas[j])
				|| (it.a == cas[j] && it.b == cas[i]))
			if given.len > 0 {
				o := given[0]
				out << PairUse{
					i:      i
					j:      j
					source: .file
					others: listed
					texts:  if o.a == cas[i] { o.texts.clone() } else { swapped(o.texts) }
				}
				continue
			}
			if listed.len == 0 {
				if !opts.ideal {
					missing << '${names[i]} with ${names[j]}'
				}
				out << PairUse{
					i:      i
					j:      j
					source: unlisted
					texts:  none_listed
				}
				continue
			}
			first := listed[0]
			// The line's own orientation, turned to i–j: A₁₂ is i's with j.
			mut t := first.texts.clone()
			if first.cas1 != cas[i] {
				t[0], t[1] = first.texts[1], first.texts[0]
			}
			out << PairUse{
				i:       i
				j:       j
				source:  .bank
				listing: first
				others:  listed[1..]
				texts:   texts(i, j, t)
			}
		}
	}
	if missing.len > 0 {
		return error("${bank.file} doesn't list ${missing.join(', ')}")
	}
	return out
}

// swap_two is a pair's two parameters, b–a turned to a–b.
fn swap_two(t []string) []string {
	return [t[1], t[0]]
}

// nrtl is NRTL for compounds given by CAS number, on ChemSep's bank:
// τᵢⱼ = Aᵢⱼ/(1.98721·T), A in cal/mol.
pub fn (db Database) nrtl(cas []string, opts ActivityOptions) !BinaryModel[activity.Nrtl] {
	recs := db.records(cas)!
	names := recs.map(it.name)
	bank := load_pair_bank(db.vendor, .nrtl)!
	uses := pair_uses(bank, cas, names, opts, fn (i int, j int, t []string) []string {
		return ['${t[0]} K / (${dechema_r}*T)', '${t[1]} K / (${dechema_r}*T)', t[2]]
	}, .ideal, ['0', '0', '0.3'], fn (t []string) []string {
		return [t[1], t[0], t[2]]
	})!
	mut pairs := []activity.NrtlPair{}
	for u in uses {
		pairs << activity.NrtlPair{
			i:      u.i
			j:      u.j
			tau_ij: activity.plain(u.texts[0], 'τ')!
			tau_ji: activity.plain(u.texts[1], 'τ')!
			alpha:  activity.plain(u.texts[2], 'α')!
		}
	}
	return BinaryModel[activity.Nrtl]{
		model: activity.new_nrtl(names, pairs)!
		pairs: uses
	}
}

// uniquac is UNIQUAC for compounds given by CAS number, on ChemSep's bank,
// τᵢⱼ = exp(−Aᵢⱼ/(1.98721·T)), with each compound's r and q. That the
// bank's R is DECHEMA's is inferred, until step 12.
pub fn (db Database) uniquac(cas []string, opts ActivityOptions) !BinaryModel[activity.Uniquac] {
	recs := db.records(cas)!
	names := recs.map(it.name)
	mut r := []f64{}
	mut q := []f64{}
	for rec in recs {
		r << rec.parameters.uniquac_r or { return error('${rec.name} has no UNIQUAC r') }
		q << rec.parameters.uniquac_q or { return error('${rec.name} has no UNIQUAC q') }
	}
	bank := load_pair_bank(db.vendor, .uniquac)!
	uses := pair_uses(bank, cas, names, opts, fn (i int, j int, t []string) []string {
		return ['exp(${t[0]} K / (-${dechema_r}*T))', 'exp(${t[1]} K / (-${dechema_r}*T))']
	}, .no_residual, ['1', '1'], swap_two)!
	mut pairs := []activity.UniquacPair{}
	for u in uses {
		pairs << activity.UniquacPair{
			i:      u.i
			j:      u.j
			tau_ij: activity.plain(u.texts[0], 'τ')!
			tau_ji: activity.plain(u.texts[1], 'τ')!
		}
	}
	return BinaryModel[activity.Uniquac]{
		model: activity.new_uniquac(names, r, q, pairs)!
		pairs: uses
	}
}

// wilson is Wilson's equation for compounds given by CAS number, on
// ChemSep's bank, Λᵢⱼ = (Vⱼ/Vᵢ)·exp(−Aᵢⱼ/(1.98721·T)), V each compound's
// Wilson volume. That form and R are inferred for this bank, until step 12.
pub fn (db Database) wilson(cas []string, opts ActivityOptions) !BinaryModel[activity.Wilson] {
	recs := db.records(cas)!
	names := recs.map(it.name)
	mut v := []string{}
	for rec in recs {
		vol := rec.parameters.wilson_volume or { return error('${rec.name} has no Wilson volume') }
		v << vol.str()
	}
	bank := load_pair_bank(db.vendor, .wilson)!
	uses := pair_uses(bank, cas, names, opts, fn [v] (i int, j int, t []string) []string {
		return ['${v[j]} / ${v[i]} * exp(${t[0]} K / (-${dechema_r}*T))',
			'${v[i]} / ${v[j]} * exp(${t[1]} K / (-${dechema_r}*T))']
	}, .ideal, ['1', '1'], swap_two)!
	mut pairs := []activity.WilsonPair{}
	for u in uses {
		pairs << activity.WilsonPair{
			i:         u.i
			j:         u.j
			lambda_ij: activity.plain(u.texts[0], 'Λ')!
			lambda_ji: activity.plain(u.texts[1], 'Λ')!
		}
	}
	return BinaryModel[activity.Wilson]{
		model: activity.new_wilson(names, pairs)!
		pairs: uses
	}
}

// unifac is a UNIFAC for compounds given by CAS number, each with its
// groups for the variant from its component file.
pub fn (db Database) unifac(variant activity.Variant, cas []string) !activity.Unifac {
	recs := db.records(cas)!
	load := load_unifac(db.vendor, variant)!
	groups := recs.map(match variant {
		.original { it.groups.unifac_vle }
		.lle { it.groups.unifac_lle }
		.dortmund { it.groups.modified_unifac }
	})
	return activity.new_unifac(load.table, recs.map(it.name), groups)!
}
