module quality

import components
import math
import transport

// transport_lines is the report's section on what transport can take of
// each compound (M5b's D3): each quantity in the files and usable, and the
// compounds whose constants a pressure correction needs and lacks.
fn transport_lines(db components.Database) ![]string {
	n := transport.quantities.len
	mut in_files := []int{len: n}
	mut usable := []int{len: n}
	mut all_files, mut all_usable := 0, 0
	mut no_dipole := []Problem{}
	mut no_critical := []Problem{}
	mut unmade := []Problem{}
	for e in db.list() {
		r := db.find(e.cas)!
		mut every_file := true
		for i, q in transport.quantities {
			if _ := r.correlations.get(q.key()) {
				in_files[i]++
			} else if i < 6 {
				every_file = false
			}
		}
		if every_file {
			all_files++
		}
		p := r.transport_pure() or {
			unmade << Problem{r.ids.cas, r.name, err.msg()}
			continue
		}
		mut every := true
		for i, s in p.sources {
			if s.has {
				usable[i]++
			} else if i < 6 {
				every = false
			}
		}
		if every {
			all_usable++
		}
		if math.is_nan(p.dipole) {
			no_dipole << Problem{r.ids.cas, r.name, 'no dipole moment'}
		}
		if math.is_nan(p.vc) || math.is_nan(p.zc) {
			no_critical << Problem{r.ids.cas, r.name, 'no critical volume or Zc'}
		}
	}
	mut o := []string{}
	o << ''
	o << '## Transport (M5b)'
	o << ''
	o << "Each compound's correlations as transport takes them (M5b's D3): one"
	o << "with an impossible range, a form reflux can't evaluate or a unit it"
	o << "can't convert is none, as are a saturated liquid's quantities at or"
	o << "above Tc and any value that isn't positive, which the correlations'"
	o << 'table above lists where the sweep finds one.'
	o << ''
	o << '| Quantity | In the files | Usable |'
	o << '|---|---|---|'
	for i, q in transport.quantities {
		o << '| `${q.key()}` | ${in_files[i]} | ${usable[i]} |'
	}
	o << '| all six but the vapour pressure | ${all_files} | ${all_usable} |'
	o << ''
	o << 'A vapour pressure falls back to the Antoine record, as everywhere.'
	o << ''
	o << '### No dipole moment (${no_dipole.len})'
	o << ''
	o << "Lucas's polar term needs it, so a gas holding one of these has no"
	o << "viscosity (M5b's D5): reflux doesn't take a missing dipole moment as 0."
	o << ''
	o << problems_table(no_dipole)
	o << '### No critical volume or Zc (${no_critical.len})'
	o << ''
	o << "Lucas's mixture rules and Stiel and Thodos' pseudo-critical constants"
	o << 'need both (D5).'
	o << ''
	o << problems_table(no_critical)
	o << '### No Pure (${unmade.len})'
	o << ''
	o << 'Without Tc, Pc, ω or M a compound has no transport, as no equation of'
	o << 'state can use it.'
	o << ''
	o << problems_table(unmade)
	return o
}
