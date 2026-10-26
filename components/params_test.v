module components

import thermo

// ChemSep's parameter files, read once. The expected counts were checked
// independently, with a Python whitespace-split over the same files.
const pr = load_kij(.pr) or { panic(err) }
const srk = load_kij(.srk) or { panic(err) }
const mc = load_mathias_copeman_pr() or { panic(err) }

const h2 = '1333-74-0'
const n2 = '7727-37-9'
const methane = '74-82-8'
const co2 = '124-38-9'
const water = '7732-18-5'

fn test_the_kij_files_read_completely() {
	assert pr.file == 'pr.ipd' && pr.comment == 'DECHEMA Peng-Robinson EOS data'
	assert pr.entries.len == 302 && pr.defects.len == 0
	assert srk.file == 'srk.ipd' && srk.comment == 'DECHEMA Soave-Redlich-Kwong EOS data'
	assert srk.entries.len == 306 && srk.defects.len == 0
	// Lines checked by hand against pr.ipd.
	assert pr.entries[0] == Kij{
		cas1:    '7440-59-7'
		cas2:    '630-08-0'
		k:       0.5463
		comment: 'Helium-4/Carbon Monoxide p203'
		line:    19
	}
	assert pr.for_pair(n2, water) == [
		Kij{
			cas1:    n2
			cas2:    water
			k:       0.08
			comment: 'Nitrogen/Water Pedersen Table 16.1'
			line:    319
		},
	]
	// A pair is found whichever way round it's asked for.
	assert pr.for_pair(water, n2) == pr.for_pair(n2, water)
	assert pr.for_pair(h2, methane).map(it.k) == [0.0263, 0.0151, 0.0596, -0.0044, -0.0144]
	assert pr.for_pair(h2, n2).map(it.k) == [0.0711, 0.1196, 0.12]
	assert srk.for_pair(h2, n2).map(it.k) == [0.0563, 0.1141, 0.1133]
	assert pr.for_pair(h2, water).len == 0
}

fn test_unknown_compounds_in_the_kij_files() {
	// Six CAS numbers in the files aren't in ChemSep LITE's database. Their
	// pairs can never be asked for; step 8's report lists them.
	db := load_default()!
	for bank in [pr, srk] {
		mut unknown := map[string]bool{}
		for e in bank.entries {
			for cas in [e.cas1, e.cas2] {
				if cas !in db.entries {
					unknown[cas] = true
				}
			}
		}
		mut list := unknown.keys()
		list.sort()
		assert list == ['101-81-5', '75-10-5', '75-45-6', '75-72-9', '76-15-3', '92-51-3'], bank.file
	}
}

fn test_a_matrix_uses_the_first_listing_and_reports_the_rest() {
	m := pr.matrix([h2, n2, methane, water])!
	nc := 4
	assert m.k.len == nc * nc
	for i in 0 .. nc {
		assert m.k[i * nc + i] == 0.0
		for j in 0 .. nc {
			assert m.k[i * nc + j] == m.k[j * nc + i]
		}
	}
	// H2/N2: three DECHEMA fits; the first is used.
	assert m.k[0 * nc + 1] == 0.0711
	p01 := m.pairs[0]
	assert p01.i == 0 && p01.j == 1 && p01.source == .file
	assert (p01.used or { Kij{} }).line == 20
	assert p01.others.map(it.k) == [0.1196, 0.12]
	assert m.k[0 * nc + 2] == 0.0263
	assert m.k[1 * nc + 3] == 0.08
	assert m.k[2 * nc + 3] == 0.4850
	// N2/methane: six fits, lines 62–67; the first is 0.0289.
	assert m.k[1 * nc + 2] == 0.0289 && m.pairs[3].others.len == 5
	// Nothing for H2/water in pr.ipd: 0, and listed.
	missing := m.missing().map('${m.cas[it.i]}/${m.cas[it.j]}')
	assert missing == ['${h2}/${water}']
	assert m.k[0 * nc + 3] == 0.0
	report := m.report()
	assert report.starts_with('kij from pr.ipd for 4 components:\n')
	assert report.contains('  ${h2}/${n2}: 0.0711, line 20 (Hydrogen/Nitrogen T=90-113K p210); also listed: 0.1196 on line 21, 0.12 on line 22\n')
	assert report.contains('  ${h2}/${water}: 0, no data\n')
	// thermo accepts it as it is.
	db := load_default()!
	mut comps := []thermo.Component{}
	for cas in m.cas {
		comps << db.find(cas)!.component()!
	}
	thermo.new_cubic(.pr, comps, kij: m.k)!
}

fn test_overrides() {
	// First listed for CO2/water is DECHEMA's 523–623 K fit; Whitson's
	// widely used 0.1896 is the third line. An override picks it.
	assert pr.for_pair(co2, water).map(it.k) == [0.0063, 0.0952, 0.1896]
	plain := pr.matrix([co2, water])!
	assert plain.k[1] == 0.0063
	m := pr.matrix([co2, water],
		overrides: [Kij{
			cas1: water
			cas2: co2
			k:    0.1896
		}]
	)!
	assert m.k == [0.0, 0.1896, 0.1896, 0.0]
	assert m.pairs[0].source == .override
	assert m.pairs[0].others.len == 3 // all three file listings stay in the record
	assert m.report().contains('${co2}/${water}: 0.1896, given by hand; also listed: 0.0063 on line 232')
	// Mistakes in the request are errors.
	bad := {
		'not in the list': KijOptions{
			overrides: [Kij{
				cas1: co2
				cas2: methane
				k:    0.1
			}]
		}
		'with itself':     KijOptions{
			overrides: [Kij{
				cas1: co2
				cas2: co2
				k:    0.1
			}]
		}
		'two overrides':   KijOptions{
			overrides: [Kij{
				cas1: co2
				cas2: water
				k:    0.1
			}, Kij{
				cas1: water
				cas2: co2
				k:    0.2
			}]
		}
	}
	for what, opts in bad {
		if _ := pr.matrix([co2, water], opts) {
			assert false, what
		}
	}
	if _ := pr.matrix([co2, water, co2]) {
		assert false, 'a repeated component'
	}
}

fn test_no_kij_for_the_older_cubics() {
	for kind in [thermo.CubicKind.vdw, .rk] {
		if _ := load_kij(kind) {
			assert false, '${kind}'
		}
	}
}

fn test_mathias_copeman_for_pr() {
	assert mc.file == 'mc-PR.lib'
	assert mc.entries.len == 668
	mut seen := map[string]int{}
	for e in mc.entries {
		seen[e.cas] += 1
	}
	assert seen.len == 597
	assert seen.values().filter(it > 1).len == 60
	// The one unreadable line: propane, with no CAS number.
	assert mc.defects.len == 1
	assert mc.defects[0].line == 73 && mc.defects[0].text.starts_with('0.603312\t0\t0\tPropane')
	// n-butane, by hand from line 99.
	b := mc.find('106-97-8') or { panic('no n-butane') }
	assert b.line == 99 && b.comment == 'N-butane'
	assert b.alpha()! == thermo.Alpha(thermo.MathiasCopeman{
		c1: 0.652753
		c2: 0.0998805
		c3: 0.0
	})
	// Nonanal is listed four times, two different ways; the first wins.
	nonanal := mc.find('124-19-6') or { panic('no nonanal') }
	assert nonanal.line == 488 && nonanal.c1 or { 0.0 } == 1.14765
	// A tab-separated line with spaces between its numbers, as line 67 is.
	assert (mc.find('115-07-1') or { panic('no propylene') }).comment == 'Propylene, 0.1% on vapor pressure'
	// `*` is an unknown value, so no α can be made from it.
	t := mc.find('100-21-0') or { panic('no terephthalic acid') }
	assert t.c1 == none
	if _ := t.alpha() {
		assert false, 'an α from unknown constants'
	}
	// Every compound in it is in the database.
	db := load_default()!
	for cas, _ in seen {
		assert cas in db.entries, cas
	}
}

fn test_the_file_format_rules() {
	text := 'free header text\n[IPD]\nComment=test\n# a comment\n\n1333-74-0 7727-37-9 .5 a pair\n1333-74-0 7727-37-9\n1333-74-1 7727-37-9 .5\n1333-74-0 1333-74-0 .5\n1333-74-0 7727-37-9 x\n1333-74-0 7727-37-9 nan\n1333-74-0\t \t7727-37-9\t-.44e-2\tmixed separators\n[NEXT]\n1333-74-0 7727-37-9 .9\n'
	b := parse_ipd(text, 't.ipd')!
	assert b.comment == 'test'
	assert b.entries.map(it.k) == [0.5, -0.0044]
	assert b.entries[0].comment == 'a pair' && b.entries[0].line == 6
	assert b.entries[1].comment == 'mixed separators'
	assert b.defects.map(it.why) == ['fewer than three fields',
		'CAS number `1333-74-1` fails its check digit', 'a compound paired with itself',
		'`x` is not a number', '`nan` is not finite']
	assert b.defects.map(it.line) == [7, 8, 9, 10, 11]
	if _ := parse_ipd('no section here', 't.ipd') {
		assert false
	} else {
		assert err.msg() == 't.ipd: no [IPD] section'
	}
	lib := parse_lib('[LIB]\n106-97-8 0.1 * 0.3 # name\n106-97-8 0.1\n106-97-8 0.1 y 0.3\n',
		't.lib')!
	assert lib.entries.len == 1 && lib.entries[0].c2 == none && lib.entries[0].c3 or { 0.0 } == 0.3
	assert lib.entries[0].comment == 'name'
	assert lib.defects.map(it.why) == ['fewer than four fields', '`y` is not a number']
}
