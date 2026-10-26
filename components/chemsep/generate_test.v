module chemsep

import crypto.sha256
import os
import components

// M2 exit criterion 2: data/components/ is exactly what the importer
// writes from the vendored XML, and the vendored XML is still what its
// NOTICE says it is.

const imported = import_files(@VMODROOT, vendored) or { panic(err) }
const generated = generate(imported) or { panic(err) }

const out_dir = os.join_path(@VMODROOT, output_dir)

fn test_the_committed_files_are_what_the_importer_writes() {
	on_disk := os.ls(out_dir)!
	for name in on_disk {
		assert name in generated, '${output_dir}/${name} is not generated any more. Run `v run tools/import_chemsep.vsh`.'
	}
	for name, text in generated {
		path := os.join_path(out_dir, name)
		assert os.exists(path), '${output_dir}/${name} is missing. Run `v run tools/import_chemsep.vsh`.'
		assert os.read_file(path)! == text, '${output_dir}/${name} is out of date. Run `v run tools/import_chemsep.vsh`.'
	}
	assert generated.len == 642 + 3
}

fn test_every_file_reads_back_as_the_record_it_was_written_from() {
	for r in imported.records {
		path := os.join_path(out_dir, '${r.ids.cas}.toml')
		back := components.read_toml(path)!
		assert back == r, '${r.ids.cas} changed in the round trip'
	}
}

fn test_the_index_lists_every_compound_file() {
	entries := components.parse_index(generated['index.toml'] or { '' })!
	assert entries.len == 642
	for e in entries {
		r := components.read_toml(os.join_path(out_dir, '${e.cas}.toml'))!
		assert components.IndexEntry.of(r) == e
	}
}

fn test_the_report_says_what_happened() {
	changes := generated['CHANGES.md'] or { '' }
	assert changes.contains('644 records became 642 compound files.')
	assert changes.contains('## Repairs (48)')
	assert changes.contains('## Duplicates dropped (2)')
	assert changes.contains('### Impossible correlation ranges (42)')
	assert changes.contains("### Equation forms reflux can't evaluate (334)")
	assert changes.contains('| 3178-22-1 | Tert-butylcyclohexane | `ids.smiles` | `CC(C)(C)C1CCCCC1<` | `CC(C)(C)C1CCCCC1` |')
}

fn test_the_vendored_files_match_their_notice() {
	dir := os.join_path(@VMODROOT, 'data', 'vendor', 'chemsep')
	notice := os.read_file(os.join_path(dir, 'NOTICE.md'))!
	mut checked := 0
	for line in notice.split_into_lines() {
		// | `chemsep1.xml` | `pcd/` | … | 3818631 | `453b0f…` |
		cells := line.split('|').map(it.trim_space())
		if cells.len != 7 || !cells[1].starts_with('`') || cells[4].len == 0
			|| !cells[4].bytes().all(it.is_digit()) {
			continue
		}
		name := cells[1].trim('`')
		bytes := os.read_bytes(os.join_path(dir, name))!
		assert bytes.len == cells[4].int(), '${name}: ${bytes.len} bytes, NOTICE says ${cells[4]}'
		assert sha256.sum(bytes).hex() == cells[5].trim('`'), '${name} no longer matches the SHA-256 in its NOTICE'
		checked++
	}
	assert checked == 5
}
