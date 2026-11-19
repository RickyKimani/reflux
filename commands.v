// reflux's command tree, which its help is rendered from (cmdline/milestones.md,
// C5). The commands still parse their own arguments: only help reads the
// tree so far.
module main

import cmdline
import flowsheet

// usage is reflux's help without colour, which a wrong command prints after
// what was wrong; kinetics_usage is the kinetics commands'.
const usage = plain_help([]string{})
const kinetics_usage = plain_help(['kinetics'])

const state_section = cmdline.Section{
	title: 'State'
	text:  "A temperature and a pressure with their units, then each species' mole
fraction by name, formula or CAS number, summing to 1."
}

// help_text is the help that args ask for, without its last newline, which
// Console adds.
fn help_text(args []string, style cmdline.Style) !string {
	return reflux_app().help(args, style)!.trim_right('\n')
}

fn plain_help(path []string) string {
	return help_text(path, cmdline.Style{}) or { err.msg() }
}

fn reflux_app() cmdline.App {
	return cmdline.App{
		root:    reflux_commands()
		version: 'reflux ${flowsheet.version()}'
	}
}

fn reflux_commands() cmdline.Command {
	return cmdline.Command{
		name:        'reflux'
		description: 'reflux, a steady-state chemical process simulator.'
		usage:       'reflux <command> [arguments] [options]'
		commands:    [
			cmdline.Command{
				name:        'solve'
				summary:     'Solve a flowsheet and print its stream table'
				description: "Solve a flowsheet file (docs/flowsheet-format.md) and print its stream
table. A failure goes to stderr with its file and line, and to the JSON if
--json names one; a CSV left from an earlier run is removed, so that nothing
on disk can be taken for this run's streams."
				usage:       'reflux solve FILE [--json OUT] [--csv OUT] [--trace]'
				arguments:   [
					cmdline.Argument{
						name:        'FILE'
						description: 'The flowsheet'
						required:    true
					},
				]
				flags:       [
					cmdline.Flag{
						name:        'json'
						description: 'Write the results as JSON (docs/report-format.md)'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'csv'
						description: 'Write the stream table as CSV'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'trace'
						description: 'Print every pass of every loop, after the table'
					},
				]
				examples:    [
					'reflux solve testdata/flowsheets/reference.toml',
					'reflux solve plant.toml --json results.json --csv streams.csv',
				]
			},
			cmdline.Command{
				name:        'study'
				summary:     'Run a study: the flowsheet at many points'
				description: 'Run a study in a flowsheet file, by NAME if it has more than one, and print
its points. A study that ran exits 0, even with failed points, which its
table lists.'
				usage:       'reflux study FILE [NAME] [--json OUT] [--csv OUT] [--workers N]'
				arguments:   [
					cmdline.Argument{
						name:        'FILE'
						description: 'The flowsheet, with its studies'
						required:    true
					},
					cmdline.Argument{
						name:        'NAME'
						description: 'The study, if the file has more than one'
					},
				]
				flags:       [
					cmdline.Flag{
						name:        'json'
						description: 'Write the points as JSON (docs/report-format.md)'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'csv'
						description: 'Write the points as CSV'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'workers'
						description: 'Run on N threads: the processors by default'
						takes_value: true
						value_name:  'N'
					},
				]
				examples:    [
					'reflux study plant.toml',
					'reflux study plant.toml conversion --workers 4 --csv points.csv',
				]
			},
			cmdline.Command{
				name:        'optimize'
				summary:     'Run an optimisation and print its optimum'
				description: 'Run an optimisation in a flowsheet file, by NAME if it has more than one,
and print its optimum and its iterations. It exits 0 only with an optimum.'
				usage:       'reflux optimize FILE [NAME] [--json OUT] [--csv OUT] [--workers N]'
				arguments:   [
					cmdline.Argument{
						name:        'FILE'
						description: 'The flowsheet, with its optimisations'
						required:    true
					},
					cmdline.Argument{
						name:        'NAME'
						description: 'The optimisation, if the file has more than one'
					},
				]
				flags:       [
					cmdline.Flag{
						name:        'json'
						description: 'Write the optimum and its iterations as JSON'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'csv'
						description: 'Write the iterations as CSV'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'workers'
						description: 'Differentiate by N inputs at once: the processors by default'
						takes_value: true
						value_name:  'N'
					},
				]
				examples:    [
					'reflux optimize plant.toml --json optimum.json',
				]
			},
			cmdline.Command{
				name:        'fit'
				summary:     'Fit rate constants, or kij, to data'
				description: 'Fit constants in a file to its data, by NAME if it has more than one fit:
rate constants to rates, or kij to bubble points. It prints each constant
with its 95% interval, and exits 0 only with a fit.'
				usage:       'reflux fit FILE [NAME] [--json OUT] [--csv OUT]'
				arguments:   [
					cmdline.Argument{
						name:        'FILE'
						description: 'The file with the fit and its data'
						required:    true
					},
					cmdline.Argument{
						name:        'NAME'
						description: 'The fit, if the file has more than one'
					},
				]
				flags:       [
					cmdline.Flag{
						name:        'json'
						description: 'Write the fit as JSON (docs/report-format.md)'
						takes_value: true
						value_name:  'OUT'
					},
					cmdline.Flag{
						name:        'csv'
						description: 'Write the residuals as CSV'
						takes_value: true
						value_name:  'OUT'
					},
				]
				examples:    [
					'reflux fit rates.toml --csv residuals.csv',
				]
			},
			kinetics_commands(),
			cmdline.Command{
				name:        'demo'
				summary:     'The thermodynamics and the flash, on n-butane'
				description: 'The thermodynamics and the flash, on compounds from the database:
n-butane at T in K and P in bar.'
				usage:       'reflux demo [T P]'
				arguments:   [
					cmdline.Argument{
						name:        'T'
						description: 'The temperature in K: 350 by default'
					},
					cmdline.Argument{
						name:        'P'
						description: 'The pressure in bar: 9.4573 by default'
					},
				]
				examples:    [
					'reflux demo',
					'reflux demo 300 2.5',
				]
			},
			cmdline.Command{
				name:    'version'
				aliases: ['--version']
				summary: 'Print the version'
				usage:   'reflux version'
			},
			cmdline.Command{
				name:      'help'
				aliases:   [
					'--help',
					'-h',
				]
				summary:   'Help for reflux, or for a command'
				usage:     'reflux help [COMMAND...]'
				arguments: [
					cmdline.Argument{
						name:        'COMMAND'
						description: 'The command, and its own: kinetics check'
						variadic:    true
					},
				]
				examples:  [
					'reflux help kinetics check',
					'reflux kinetics check --help',
				]
			},
		]
		examples:    [
			'# Solve a flowsheet, and write its results as JSON',
			'reflux solve plant.toml --json results.json',
			'',
			'# Check the reactions in a file',
			'reflux kinetics check rates.toml',
			'',
			'# More on one command',
			'reflux help solve',
		]
		sections:    [
			cmdline.Section{
				title: 'Exit status'
				text:  '0 when it worked, 1 when the flowsheet failed or a report could not be
written, and 2 when the command is wrong. A study that ran is 0, even with
failed points, which its table lists; an optimisation is 0 only with an
optimum, and a fit only with a fit.'
			},
		]
	}
}

fn kinetics_commands() cmdline.Command {
	file := cmdline.Argument{
		name:        'FILE'
		description: 'The file with the reactions'
		required:    true
	}
	return cmdline.Command{
		name:        'kinetics'
		summary:     "Check a file's reactions, or explain their rates"
		description: 'Check the reactions in a file, or explain their rates at a state.'
		usage:       'reflux kinetics <command> FILE [T P NAME=y ...]'
		commands:    [
			cmdline.Command{
				name:        'check'
				summary:     "Check a file's reactions, at a state if given"
				description: "Read a file's reactions (docs/flowsheet-format.md), which refuses an
equation that does not balance or a rate in the wrong units, and print what
every other check finds. Given a state, a fitted K is judged at its T, each
typed constant says how far it moves the rates there, and each reversible
rate is held to zero at the equilibrium the reactions reach from it, in the
ideal gas, and to running towards it at the state."
				usage:       'reflux kinetics check FILE [T P NAME=y ...]'
				arguments:   [
					file,
					cmdline.Argument{
						name:        'T'
						description: "The temperature, with its unit: '900 K'"
					},
					cmdline.Argument{
						name:        'P'
						description: "The pressure, with its unit: '15 bar'"
					},
					cmdline.Argument{
						name:        'NAME=y'
						description: "Each species' mole fraction"
						variadic:    true
					},
				]
				examples:    [
					'reflux kinetics check rates.toml',
					"reflux kinetics check rates.toml '900 K' '15 bar' CH4=0.25 H2O=0.75",
				]
				sections:    [state_section]
			},
			cmdline.Command{
				name:        'explain'
				summary:     "Explain a file's rates at a state"
				description: 'Print each rate as reflux read it, and its factor table at the state.'
				usage:       'reflux kinetics explain FILE T P NAME=y ...'
				arguments:   [
					file,
					cmdline.Argument{
						name:        'T'
						description: "The temperature, with its unit: '900 K'"
						required:    true
					},
					cmdline.Argument{
						name:        'P'
						description: "The pressure, with its unit: '15 bar'"
						required:    true
					},
					cmdline.Argument{
						name:        'NAME=y'
						description: "Each species' mole fraction"
						required:    true
						variadic:    true
					},
				]
				examples:    [
					"reflux kinetics explain plant.toml '900 K' '15 bar' CH4=0.2 H2O=0.45 H2=0.2 CO=0.1 CO2=0.05",
				]
				sections:    [state_section]
			},
		]
		examples:    [
			'reflux kinetics check rates.toml',
			"reflux kinetics explain plant.toml '900 K' '15 bar' CH4=0.2 H2O=0.45 H2=0.2 CO=0.1 CO2=0.05",
		]
		sections:    [
			state_section,
			cmdline.Section{
				title: 'Exit status'
				text:  '0 when the file was read, 1 when it was not or the state is wrong, and 2
when the command is wrong.'
			},
		]
	}
}
