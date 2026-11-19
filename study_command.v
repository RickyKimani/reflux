module main

import flowsheet
import quantity
import reactors
import study
import time

// Named is a command over one of a file's named sections, `reflux study`'s
// or `reflux optimize`'s: the file and its --json and --csv, the section's
// name if given, and the workers.
struct Named {
	file_args []string
	name      string
	workers   int
}

// parse_named reads `reflux COMMAND FILE [NAME] [--json OUT] [--csv OUT]
// [--workers N]`, or complains and returns none, for an exit status of 2.
fn parse_named(command string, args []string, mut c Console) ?Named {
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
				c.complain('reflux ${command}: --workers needs a whole number from 1 to 1024, not "${value}"\n\n${usage}')
				return none
			}
			workers = n
			continue
		}
		if arg == '--trace' {
			c.complain("reflux ${command}: --trace is solve's\n\n${usage}")
			return none
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
			c.complain('reflux ${command}: one ${runs(command)} at a time, not ${name} and ${arg}\n\n${usage}')
			return none
		}
		file_args << arg
	}
	return Named{
		file_args: file_args
		name:      name
		workers:   workers
	}
}

// runs is what a command runs, as its messages name it.
fn runs(command string) string {
	return match command {
		'study' { 'study' }
		'fit' { 'fit' }
		else { 'optimisation' }
	}
}

// pick is the section a command names, or the file's only one, or none
// after a complaint: the exit status, 1 for a file with none, 2 for a name
// it doesn't have or one left out among several.
fn pick(command string, file string, names []string, name string, mut c Console) !int {
	plural := match command {
		'study' { 'studies' }
		'fit' { 'fits' }
		else { 'optimisations' }
	}
	table := match command {
		'study' { 'study' }
		'fit' { 'fit' }
		else { 'optimize' }
	}
	if names.len == 0 {
		c.complain('reflux: ${file} has no ${plural}: add a [${table}.NAME] (docs/flowsheet-format.md)')
		return error('1')
	}
	if name == '' && names.len > 1 {
		c.complain('reflux ${command}: ${file} has ${plural} ${names.join(', ')}: name one\n\n${usage}')
		return error('2')
	}
	k := if name == '' { 0 } else { names.index(name) }
	if k < 0 {
		close := quantity.closest(name, names)
		hint := if close.len > 0 { ' (did you mean ${close.join(', ')}?)' } else { '' }
		c.complain('reflux ${command}: ${file} has no ${runs(command)} ${name}${hint}\n\n${usage}')
		return error('2')
	}
	return k
}

// took is a time in seconds as the commands print it.
fn took(seconds f64, workers int) string {
	t := if seconds < 1 {
		'${flowsheet.figures(1e3 * seconds)} ms'
	} else {
		'${flowsheet.figures(seconds)} s'
	}
	on := if workers == 1 { '1 worker' } else { '${workers} workers' }
	return 'In ${t}, on ${on}.'
}

// study_command runs `reflux study FILE [NAME] [--json OUT] [--csv OUT]
// [--workers N]` (M8's D17): a study of the file's, named or its only one,
// on N workers, the processors by default (D8), printed as a table, with
// its JSON and CSV if asked for. It exits 0 when the study ran, failed
// points and all, 1 when it couldn't run, and 2 for a wrong command.
fn study_command(args []string, mut c Console) int {
	a := parse_named('study', args, mut c) or { return 2 }
	o := options(a.file_args) or {
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
	k := pick('study', o.file, loaded.studies.map(it.name), a.name, mut c) or {
		return err.msg().int()
	}
	t0 := time.sys_mono_now()
	result := study.run(loaded, loaded.studies[k], workers: a.workers)
	seconds := f64(time.sys_mono_now() - t0) / 1e9
	r := result.report(o.file, loaded.specs.map(it.name))
	c.say(r.table() + '\n' + took(seconds, a.workers))
	if o.json != '' && !write(o.json, r.json(), mut c) {
		return 1
	}
	if o.csv != '' && !write(o.csv, r.csv(), mut c) {
		return 1
	}
	return 0
}

// optimize_command runs `reflux optimize FILE [NAME] [--json OUT] [--csv
// OUT] [--workers N]` (D13, D17): an optimisation of the file's, named or
// its only one, its gradient's inputs on N workers, printed as a table,
// with its JSON and CSV if asked for. It exits 0 with an optimum, 1 when it
// couldn't run or found none, and 2 for a wrong command.
fn optimize_command(args []string, mut c Console) int {
	a := parse_named('optimize', args, mut c) or { return 2 }
	o := options(a.file_args) or {
		c.complain('reflux optimize: ${err.msg()}\n\n${usage}')
		return 2
	}
	loaded := flowsheet.load(o.file, makers: reactors.makers()) or {
		r := study.optimize_load_failure(o.file, err)
		e := r.error or { flowsheet.ErrorReport{} }
		c.complain('reflux: ${located(o.file, e)}')
		if o.json != '' {
			write(o.json, r.json(), mut c)
		}
		return 1
	}
	k := pick('optimize', o.file, loaded.optimizations.map(it.name), a.name, mut c) or {
		return err.msg().int()
	}
	t0 := time.sys_mono_now()
	optimum := study.optimize(loaded, loaded.optimizations[k], workers: a.workers)
	seconds := f64(time.sys_mono_now() - t0) / 1e9
	r := optimum.report(o.file, loaded.specs.map(it.name))
	c.say(r.table() + '\n' + took(seconds, a.workers))
	if o.json != '' && !write(o.json, r.json(), mut c) {
		return 1
	}
	if o.csv != '' && !write(o.csv, r.csv(), mut c) {
		return 1
	}
	return if optimum.ok { 0 } else { 1 }
}
