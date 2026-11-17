module main

import flowsheet
import quantity
import reactors
import study
import time

// study_command runs `reflux study FILE [NAME] [--json OUT] [--csv OUT]
// [--workers N]` (M8's D17): a study of the file's, named or its only one,
// on N workers, the processors by default (D8), printed as a table, with
// its JSON and CSV if asked for. It exits 0 when the study ran, failed
// points and all, 1 when it couldn't run, and 2 for a wrong command.
fn study_command(args []string, mut c Console) int {
	mut file_args := []string{}
	mut name := ''
	mut workers := cpus()
	mut positional := 0
	mut i := 0
	for i < args.len {
		arg := args[i]
		i++
		if arg in ['--json', '--csv'] {
			file_args << arg
			if i < args.len {
				file_args << args[i]
				i++
			}
			continue
		}
		if arg == '--workers' || arg.starts_with('--workers=') {
			mut value := arg.all_after_first('=')
			if !arg.contains('=') {
				value = args[i] or { '' }
				i++
			}
			n := value.int()
			if n < 1 || n > 1024 || n.str() != value {
				c.complain('reflux study: --workers needs a whole number from 1 to 1024, not "${value}"\n\n${usage}')
				return 2
			}
			workers = n
			continue
		}
		if arg == '--trace' {
			c.complain("reflux study: --trace is solve's\n\n${usage}")
			return 2
		}
		if arg.starts_with('-') {
			file_args << arg
			continue
		}
		positional++
		if positional == 2 {
			name = arg
			continue
		}
		if positional > 2 {
			c.complain('reflux study: one study at a time, not ${name} and ${arg}\n\n${usage}')
			return 2
		}
		file_args << arg
	}
	o := options(file_args) or {
		c.complain('reflux study: ${err.msg()}\n\n${usage}')
		return 2
	}
	loaded := flowsheet.load(o.file, makers: reactors.makers()) or {
		r := study.load_failure(o.file, err)
		e := r.error or { flowsheet.ErrorReport{} }
		c.complain('reflux: ${located(o.file, e)}')
		if o.json != '' {
			write(o.json, r.json(), mut c)
		}
		return 1
	}
	names := loaded.studies.map(it.name)
	if names.len == 0 {
		c.complain('reflux: ${o.file} has no studies: add a [study.NAME] (docs/flowsheet-format.md)')
		return 1
	}
	if name == '' && names.len > 1 {
		c.complain('reflux study: ${o.file} has studies ${names.join(', ')}: name one\n\n${usage}')
		return 2
	}
	plans := if name == '' { loaded.studies } else { loaded.studies.filter(it.name == name) }
	if plans.len == 0 {
		close := quantity.closest(name, names)
		hint := if close.len > 0 { ' (did you mean ${close.join(', ')}?)' } else { '' }
		c.complain('reflux study: ${o.file} has no study ${name}${hint}\n\n${usage}')
		return 2
	}
	t0 := time.sys_mono_now()
	result := study.run(loaded, plans[0], workers: workers)
	seconds := f64(time.sys_mono_now() - t0) / 1e9
	r := result.report(o.file, loaded.specs.map(it.name))
	took := if seconds < 1 {
		'${flowsheet.figures(1e3 * seconds)} ms'
	} else {
		'${flowsheet.figures(seconds)} s'
	}
	on := if workers == 1 { '1 worker' } else { '${workers} workers' }
	c.say(r.table() + '\nIn ${took}, on ${on}.')
	if o.json != '' && !write(o.json, r.json(), mut c) {
		return 1
	}
	if o.csv != '' && !write(o.csv, r.csv(), mut c) {
		return 1
	}
	return 0
}
