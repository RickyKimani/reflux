// reflux's command line (M4, D11; M9, D12):
//
//     reflux solve plant.toml [--json results.json] [--csv streams.csv] [--trace]
//     reflux kinetics check FILE [T P NAME=y ...]
//     reflux kinetics explain FILE T P NAME=y ...
//     reflux demo [T P]
//
// From the source, `v run . solve testdata/flowsheets/reference.toml`.
module main

import flowsheet
import os
import quantity
import strings

const usage = 'usage:
  reflux solve FILE [--json OUT] [--csv OUT] [--trace]
      Solve a flowsheet file (docs/flowsheet-format.md) and print its stream
      table. --json writes the results (docs/report-format.md), --csv the
      streams, and --trace prints every pass of every loop.
  reflux kinetics check FILE [T P NAME=y ...]
  reflux kinetics explain FILE T P NAME=y ...
      Check the reactions in a file, or explain its rates at a state:
      reflux kinetics help says more.
  reflux demo [T P]
      The thermodynamics and the flash, on compounds from the database:
      n-butane at T in K and P in bar, by default 350 K and 9.4573 bar.
  reflux version
  reflux help

The exit status is 0 when it worked, 1 when the flowsheet failed or a report
could not be written, and 2 when the command is wrong.'

fn main() {
	mut c := Console{
		live: true
	}
	exit(run(os.args[1..], mut c))
}

// Console is where a command's output goes: printed as it comes when the
// program runs, and kept for main_test.v to read when it doesn't.
struct Console {
	live bool
mut:
	out strings.Builder = strings.new_builder(1024)
	err strings.Builder = strings.new_builder(256)
}

fn (mut c Console) say(s string) {
	if c.live {
		println(s)
	} else {
		c.out.writeln(s)
	}
}

fn (mut c Console) complain(s string) {
	if c.live {
		eprintln(s)
	} else {
		c.err.writeln(s)
	}
}

// run runs the command in args, and returns its exit status.
fn run(args []string, mut c Console) int {
	command := args[0] or { '' }
	match command {
		'solve' {
			return solve(args[1..], mut c)
		}
		'demo' {
			return demo(args[1..], mut c)
		}
		'kinetics' {
			return kinetics_command(args[1..], mut c)
		}
		'version', '--version' {
			c.say('reflux ${flowsheet.version()}')
			return 0
		}
		'help', '--help', '-h' {
			c.say(usage)
			return 0
		}
		'' {
			c.complain(usage)
			return 2
		}
		else {
			hint := did_you_mean(command, ['solve', 'kinetics', 'demo', 'version', 'help'])
			c.complain('reflux: no command ${command}${hint}\n\n${usage}')
			return 2
		}
	}
}

fn did_you_mean(name string, accepted []string) string {
	close := quantity.closest(name, accepted)
	return if close.len > 0 { ' (did you mean ${close.join(', ')}?)' } else { '' }
}

// Options are solve's.
struct Options {
mut:
	file  string
	json  string
	csv   string
	trace bool
}

// options reads solve's arguments. An option's file can follow it, or an `=`.
fn options(args []string) !Options {
	mut o := Options{}
	mut i := 0
	for i < args.len {
		arg := args[i]
		i++
		name := arg.all_before('=')
		if name in ['--json', '--csv'] {
			mut value := arg.all_after_first('=')
			if !arg.contains('=') {
				if i >= args.len || args[i].starts_with('-') {
					return error('${name} needs a file name')
				}
				value = args[i]
				i++
			}
			if value == '' {
				return error('${name} needs a file name')
			}
			if (name == '--json' && o.json != '') || (name == '--csv' && o.csv != '') {
				return error('${name} is given twice')
			}
			if name == '--json' {
				o.json = value
			} else {
				o.csv = value
			}
		} else if arg == '--trace' {
			o.trace = true
		} else if arg.starts_with('-') {
			return error('no option ${arg}${did_you_mean(name, ['--json', '--csv', '--trace'])}')
		} else if o.file != '' {
			return error('one flowsheet at a time, not ${o.file} and ${arg}')
		} else {
			o.file = arg
		}
	}
	if o.file == '' {
		return error('it needs a flowsheet file')
	}
	// A report mustn't be written over the flowsheet, or over the other report.
	for flag, path in {
		'--json': o.json
		'--csv':  o.csv
	} {
		if path != '' && same_file(path, o.file) {
			return error('${flag} ${path} would overwrite the flowsheet')
		}
	}
	// os.real_path('') is '', so no file is the same as no --csv.
	if o.json != '' && same_file(o.json, o.csv) {
		return error('--json and --csv name the same file, ${o.json}')
	}
	return o
}

// same_file is whether two paths name one file. On Windows names differ in
// case and still name one file.
fn same_file(a string, b string) bool {
	x, y := os.real_path(a), os.real_path(b)
	return $if windows { x.to_lower() == y.to_lower() } $else { x == y }
}

// solve solves a flowsheet file and reports it (D11): a stream table on
// stdout, and the JSON and CSV it's asked for.
fn solve(args []string, mut c Console) int {
	o := options(args) or {
		c.complain('reflux solve: ${err.msg()}\n\n${usage}')
		return 2
	}
	loaded := flowsheet.load(o.file) or { return failed(o, flowsheet.load_failure(err), mut c) }
	sol := loaded.solve() or { return failed(o, loaded.failure(err), mut c) }
	r := flowsheet.Report{
		...loaded.report(sol)
		file: o.file
	}
	c.say(r.table())
	if o.trace {
		c.say(r.trace())
	}
	if o.json != '' && !write(o.json, r.json(), mut c) {
		return 1
	}
	if o.csv != '' && !write(o.csv, r.csv(), mut c) {
		return 1
	}
	return 0
}

// failed reports a failure: on stderr, with the passes so far if asked for,
// and as JSON if that was asked for. A CSV left from an earlier run is removed,
// so that nothing on disk can be taken for this run's streams.
fn failed(o Options, report flowsheet.Report, mut c Console) int {
	r := flowsheet.Report{
		...report
		file: o.file
	}
	if o.trace {
		c.say(r.trace())
	}
	e := r.error or { flowsheet.ErrorReport{} }
	c.complain('reflux: ${located(o.file, e)}')
	if !o.trace && r.blocks.any(it.residuals.len > 0) {
		c.complain('(--trace shows every pass)')
	}
	if o.csv != '' && os.exists(o.csv) {
		os.rm(o.csv) or { c.complain('reflux: ${o.csv}: ${err.msg()}') }
	}
	if o.json != '' {
		write(o.json, r.json(), mut c)
	}
	return 1
}

// located is a failure as a compiler puts one: the file and line first, which
// an editor can jump to, then what and why.
fn located(file string, e flowsheet.ErrorReport) string {
	mut place := file
	if e.line > 0 {
		place += ':${e.line}'
	}
	// A file that can't be read is its own key.
	what := if e.key != '' && e.key != file {
		'${e.key}: '
	} else if e.unit != '' {
		'${e.unit}: '
	} else {
		''
	}
	hint := if e.suggestions.len > 0 { ' (did you mean ${e.suggestions.join(', ')}?)' } else { '' }
	return '${place}: ${what}${e.reason}${hint}'
}

// write writes a report's file, and says so if it can't.
fn write(path string, text string, mut c Console) bool {
	os.write_file(path, text) or {
		c.complain('reflux: ${path}: ${err.msg()}')
		return false
	}
	return true
}
