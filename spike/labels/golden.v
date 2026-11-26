module main

import flowsheet
import reactors
import thermo

// golden_labels counts, in every golden flowsheet's solution, the streams
// split in two whose phases D9's Π labels would rename.
fn golden_labels() ! {
	mut renamed := 0
	mut splits := 0
	for name in ['reference', 'every_unit', 'no_recycle', 'interlocking', 'saturated_loop', 'reformer',
		'methanol'] {
		path := 'testdata/flowsheets/${name}.toml'
		loaded := flowsheet.load(path, makers: reactors.makers())!
		sol := loaded.solve()!
		eos := loaded.fluid.helmholtz()!
		for sname, st in sol.streams {
			r := st.result
			if r.phases.len != 2 {
				continue
			}
			splits++
			ll := thermo.phase_label(eos, r.phases[0].state)
			lh := thermo.phase_label(eos, r.phases[1].state)
			if !(ll == .vapor && lh == .liquid) {
				renamed++
				println('golden ${name}, stream ${sname}: ${ll} and ${lh}')
			}
		}
	}
	println('golden flowsheets: ${splits} two-phase streams, ${renamed} renamed by Π')
}
