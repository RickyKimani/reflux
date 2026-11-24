// M6 step 1, criterion 2: ChemSep's activity banks and UNIFAC tables, read
// as vendored; every defect, repair and correction; and the compounds each
// table covers. Each table against thermo's copy of the published one is
// components/quality's (docs/data-quality.md).
module components

import activity
import os

fn test_the_banks_read() {
	db := load_default()!
	mut counts := map[string]int{}
	mut repeated := map[string]int{}
	for kind in [ActivityBank.nrtl, .uniquac, .wilson] {
		b := load_pair_bank(db.vendor, kind)!
		assert b.units == 'cal/mol'
		assert b.defects.len == 0, '${b.file}: ${b.defects}'
		counts[b.file] = b.entries.len
		mut seen := map[string]bool{}
		for e in b.entries {
			key := if e.cas1 < e.cas2 { '${e.cas1} ${e.cas2}' } else { '${e.cas2} ${e.cas1}' }
			if seen[key] {
				repeated[b.file] = (repeated[b.file] or { 0 }) + 1
			}
			seen[key] = true
		}
	}
	println('pairs ${counts}, listed again ${repeated}')
	// Each listing after a pair's first, which reflux doesn't use (as kᵢⱼ).
	assert repeated == {
		'nrtl.ipd':    23
		'uniquac.ipd': 18
		'wilson.ipd':  4
	}
	assert counts == {
		'nrtl.ipd':    352
		'uniquac.ipd': 291
		'wilson.ipd':  366
	}
	nrtl := load_pair_bank(db.vendor, .nrtl)!
	mw := nrtl.for_pair('7732-18-5', '67-56-1')
	assert mw.len == 1 && mw[0].cas1 == '67-56-1' && mw[0].values == [-189.0469, 792.8020, 0.2999]
	assert mw[0].comment == 'Methanol/Water p61 1/1a'
	// A number printed with a leading point gets a 0 for expr, the same
	// double: Halothane/Methanol's α, `.187e-1`.
	halo := nrtl.for_pair('151-67-7', '67-56-1')[0]
	assert halo.texts[2] == '0.187e-1' && halo.values[2] == 0.0187
	henry := load_henry(db.vendor)!
	assert henry.entries.len == 41 && henry.defects.len == 0
	assert henry.entries.all(it.solvent == '7732-18-5' && it.eq == 121 && it.b == 1)
	h2 := henry.entries.filter(it.gas == '1333-74-0')
	assert h2.len == 1 && h2[0].tmin == 273.162 && h2[0].coeffs[0] == 86.3574
}

fn test_the_unifac_tables_read() {
	db := load_default()!
	mut summary := map[string][]int{}
	for v in [activity.Variant.original, .lle, .dortmund] {
		load := load_unifac(db.vendor, v)!
		t := load.table
		summary[v.name()] = [t.subgroups.len, t.mains.len, t.interactions.len, load.defects.len,
			load.repairs.len]
	}
	println(summary)
	// Subgroups, main groups, parameters, defects and repairs (with the
	// published corrections), all listed in docs/data-quality.md.
	assert summary == {
		'UNIFAC':          [133, 61, 1364, 0, 15]
		'UNIFAC-LLE':      [59, 33, 516, 0, 3]
		'UNIFAC-Dortmund': [127, 64, 1528, 3, 17]
	}
	// The repairs that change a model: furfural's LLE subgroup on its own
	// main group, not DMFA's; benzonitrile's AC-CN on AC-CN's, not AC-Br's;
	// Dortmund's quaternary C with its published Q of 0.
	lle := load_unifac(db.vendor, .lle)!.table
	assert lle.subgroups[22] or { activity.Subgroup{} }.main == 12
	original := load_unifac(db.vendor, .original)!.table
	assert original.subgroups[128] or { activity.Subgroup{} }.main == 64
	assert original.subgroups[129] or { activity.Subgroup{} }.main == 57
	dortmund := load_unifac(db.vendor, .dortmund)!.table
	assert dortmund.subgroups[4] or { activity.Subgroup{} }.q == 0.0
	// And the c rows whose names were copied from the row above keep their
	// numbers' pairs: ACCH2 with PF6, CCl4 with HCONR.
	assert dortmund.interactions[activity.interaction_key(4, 99)] or { activity.Interaction{} }.c == 2.3928e-2
	assert dortmund.interactions[activity.interaction_key(24, 49)] or { activity.Interaction{} }.c == 1.8550e-3
}

fn test_every_compound_has_its_groups() {
	db := load_default()!
	mut counts := map[string][]int{}
	for v in [activity.Variant.original, .lle, .dortmund] {
		t := load_unifac(db.vendor, v)!.table
		mut with, mut usable := 0, 0
		for e in db.list() {
			r := db.find(e.cas)!
			g := match v {
				.original { r.groups.unifac_vle }
				.lle { r.groups.unifac_lle }
				.dortmund { r.groups.modified_unifac }
			}
			if g.len == 0 {
				continue
			}
			with++
			if g.keys().all(it in t.subgroups) {
				usable++
			}
		}
		counts[v.name()] = [with, usable]
	}
	// Dimethyl sulfoxide's LLE group is subgroup 68, which the LLE table
	// hasn't (its DMSO is 57): listed in docs/data-quality.md.
	assert counts == {
		'UNIFAC':          [584, 584]
		'UNIFAC-LLE':      [484, 483]
		'UNIFAC-Dortmund': [548, 548]
	}
}

// The reader's rules for a table that contradicts itself, on small files
// written for them, since today's files can't reach every rule: a pair
// listed twice keeps its first; a row whose name isn't its number's, where
// its numbers' pair is listed already, takes its name's; a number the list
// hasn't, named nothing the list names, is a new main group; a subgroup's
// main group by its name; a molar mass with a decimal comma.
fn test_the_rules_for_a_table_that_contradicts_itself() {
	dir := os.join_path(os.temp_dir(), 'reflux_unifac_rules_${os.getpid()}')
	os.mkdir_all(dir)!
	defer {
		os.rmdir_all(dir) or {}
	}
	files := {
		'unifacvl.gct':     '3 Total number of main groups\n1 CH2\n2 OH\n3 H2O\n1 CH2 2 OH 986.5 156.4\n1 CH2 3 H2O 1318.0 300.0\n2 OH 3 H2O 353.5 -229.1\n1 CH2 2 OH 1.0 2.0\n'
		'unifacvl_ext.gct': '1 Total number of new main groups\n4 ACCN\n4 ACCN 1 CH2 10.0 20.0\n4 ACCN 2 OH 5.0 6.0\n4 ACCN 2 H2O 30.0 40.0\n5 CHX 3 H2O 7.0 8.0\n'
		'unifacrq.gct':     '3 Total number of subgroups\n1 CH3 1 CH2 0.9011 0.848 15.035\n2 OH 2 OH 1.0 1.2 17.007\n3 H2O 3 H2O 0.92 1.4 18.015\n'
		'unifacrq_ext.gct': '2 Total number of new sub groups\n4 ACCN 2 ACCN 1.5 1.2 39,02\n5 CHX 5 CHX 1.0 1.0 13.0\n'
	}
	for name, text in files {
		os.write_file(os.join_path(dir, name), text)!
	}
	load := load_unifac(dir, .original)!
	t := load.table
	assert t.mains == {
		1: 'CH2'
		2: 'OH'
		3: 'H2O'
		4: 'ACCN'
		5: 'CHX'
	}
	got := fn [t] (m int, n int) []f64 {
		p := t.interactions[activity.interaction_key(m, n)] or { return [] }
		return [p.a, p.b, p.c]
	}
	assert got(1, 2) == [986.5, 0.0, 0.0] && got(2, 1) == [156.4, 0.0, 0.0] // the first listing
	assert got(4, 2) == [5.0, 0.0, 0.0] // the row whose names are its numbers'
	assert got(4, 3) == [30.0, 0.0, 0.0] && got(3, 4) == [40.0, 0.0, 0.0] // by its name
	assert got(5, 3) == [7.0, 0.0, 0.0]
	assert t.subgroups[4] or { activity.Subgroup{} }.main == 4
	assert t.subgroups[4] or { activity.Subgroup{} }.mw == 39.02
	assert load.defects.len == 1 && load.defects[0].line == 8 && load.defects[0].why == 'a pair listed again'
	whys := load.repairs.map(it.why)
	assert load.repairs.len == 4, '${whys}'
	assert whys.any(it.contains('its names, ACCN and H2O, are those numbers'))
	assert whys.any(it.contains("a main group the list doesn't number"))
	assert whys.any(it.contains('its main group is named ACCN, which is 4'))
	assert whys.any(it == 'a molar mass with a decimal comma')
}
