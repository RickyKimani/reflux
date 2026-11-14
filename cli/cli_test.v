module cli

fn test_help_is_deterministic_and_aligned() {
	command := Command{
		name:        'reflux'
		usage:       'reflux <command>'
		description: 'A simulator.'
		commands:    [
			Command{
				name:    'solve'
				summary: 'Solve a flowsheet.'
			},
			Command{
				name:    'version'
				aliases: ['--version']
				summary: 'Print the version.'
			},
		]
	}
	expected := 'Usage:\n  reflux <command>\n\nA simulator.\n\nCommands:\n  solve                Solve a flowsheet.\n  version (--version)  Print the version.\n'
	assert help(command) == expected
	assert !expected.contains('\t')
}

fn test_nested_command_help_contains_all_sections() {
	app := App{
		root: Command{
			name:     'reflux'
			commands: [
				Command{
					name:     'kinetics'
					commands: [
						Command{
							name:        'check'
							usage:       'reflux kinetics check FILE'
							description: 'Check reaction definitions.'
							arguments:   [
								Argument{
									name:        'FILE'
									description: 'Kinetics file.'
									required:    true
								},
							]
							flags:       [
								Flag{
									name:        'strict'
									description: 'Treat warnings as errors.'
								},
							]
							examples:    ['reflux kinetics check rates.toml --strict']
						},
					]
				},
			]
		}
	}
	got := app.help(['kinetics', 'check'], Style{}) or {
		assert false
		return
	}
	expected := 'Usage:\n  reflux kinetics check FILE\n\nCheck reaction definitions.\n\nArguments:\n  FILE (required)  Kinetics file.\n\nOptions:\n  --strict  Treat warnings as errors.\n\nExamples:\n  reflux kinetics check rates.toml --strict\n'
	assert got == expected
	plain := app.help(['kinetics', 'check', '--help'], Style{}) or {
		assert false
		return
	}
	assert plain == expected
	colored := app.help(['kinetics', 'check', '-h'], Style{
		terminal: Terminal{ color: true }
	}) or {
		assert false
		return
	}
	assert colored.starts_with('\x1b[1;36mUsage:\x1b[0m\n')
	assert colored.contains('\x1b[1;36mArguments:\x1b[0m\n')
	assert colored.contains('\x1b[1;36mOptions:\x1b[0m\n')
	assert colored.contains('\x1b[1;36mExamples:\x1b[0m\n')
}

fn test_parser_routes_and_normalizes_flags() {
	command := Command{
		name:     'reflux'
		commands: [
			Command{
				name:      'solve'
				arguments: [Argument{ name: 'FILE', required: true }]
				flags:     [
					Flag{ name: 'json', alias: 'j', takes_value: true },
					Flag{ name: 'trace' },
				]
			},
		]
	}
	invocation := parse(command, ['solve', 'plant.toml', '-j', 'result.json', '--trace']) or {
		assert false
		return
	}
	assert invocation.command == ['solve']
	assert invocation.arguments == ['plant.toml']
	assert invocation.values['json'] == 'result.json'
	assert invocation.switches['trace']
}

fn test_parser_reports_suggestions() {
	command := Command{
		name:     'reflux'
		commands: [Command{ name: 'solve' }]
	}
	parse(command, ['solv']) or {
		assert err.msg() == 'unknown command solv (did you mean solve?)'
		return
	}
	assert false
}

fn test_parser_rejects_invalid_forms() {
	command := Command{
		name:     'reflux'
		commands: [
			Command{
				name:      'solve'
				arguments: [Argument{ name: 'FILE', required: true }]
				flags:     [
					Flag{ name: 'json', takes_value: true },
					Flag{ name: 'trace' },
				]
			},
		]
	}
	assert_parse_error(command, ['solve', '--json'], '--json needs a value')
	assert_parse_error(command, ['solve', 'input.toml', '--trace', '--trace'],
		'--trace is given twice')
	assert_parse_error(command, ['solve', 'input.toml', '--json=', '--json=next.json'],
		'--json is given twice')
	assert_parse_error(command, ['solve', 'input.toml', '--trac'],
		'unknown option --trac (did you mean --trace?)')
	assert_parse_error(command, ['solve'], 'solve needs at least 1 argument')
	assert_parse_error(command, ['solve', 'a', 'b'], 'solve takes 1 argument, not 2')
}

fn assert_parse_error(command Command, args []string, expected string) {
	parse(command, args) or {
		assert err.msg() == expected
		return
	}
	panic('expected parser error: ${args.join(' ')}')
}

fn test_parser_supports_command_and_flag_aliases() {
	command := Command{
		name:     'reflux'
		commands: [
			Command{
				name:      'solve'
				aliases:   ['run']
				arguments: [Argument{ name: 'FILE', required: true }]
				flags:     [Flag{ name: 'json', alias: 'j', takes_value: true }]
			},
		]
	}
	invocation := parse(command, ['run', 'input.toml', '-j=output.json']) or {
		assert false
		return
	}
	assert invocation.command == ['solve']
	assert invocation.values['json'] == 'output.json'
}

fn test_validation_rejects_ambiguous_metadata() {
	app := App{
		root: Command{
			name:  'reflux'
			flags: [
				Flag{ name: 'json' },
				Flag{ name: 'output', alias: 'json' },
			]
		}
	}
	app.validate() or {
		assert err.msg() == 'root: duplicate flag json'
		return
	}
	assert false
}

fn test_diagnostic_respects_style() {
	plain := Diagnostic{
		kind:    .error
		message: 'bad input'
		hint:    'check the file'
		usage:   'reflux solve FILE'
	}.render(Style{})
	assert plain == 'Error: bad input\nHint: check the file\nUsage: reflux solve FILE\n'

	colored := Diagnostic{
		kind:    .warning
		message: 'using a default'
	}.render(Style{
		terminal: Terminal{ color: true }
	})
	assert colored == '\x1b[33mWarning:\x1b[0m using a default\n'
}

struct TestHandler {}

fn (_ TestHandler) execute(context Context) !int {
	invocation := context.invocation
	if invocation.arguments.len == 0 {
		return 3
	}
	return 0
}

fn test_app_builtin_outcomes_and_dispatch() {
	app := App{
		root:     Command{
			name:     'reflux'
			commands: [
				Command{
					name:      'run'
					arguments: [Argument{ name: 'FILE', required: true }]
				},
				Command{ name: 'version', aliases: ['--version'] },
			]
		}
		version:  'reflux 0.3.0'
		commands: [
			RegisteredCommand{
				path:    ['run']
				handler: TestHandler{}
			},
		]
	}
	help := app.run([], Style{}) or {
		assert false
		return
	}
	assert help.kind == .help
	assert help.code == 0
	version := app.run(['--version'], Style{}) or {
		assert false
		return
	}
	assert version.kind == .version
	assert version.text == 'reflux 0.3.0'
	executed := app.run(['run', 'input.toml'], Style{}) or {
		assert false
		return
	}
	assert executed.kind == .executed
	assert executed.code == 0
}

fn test_app_rejects_unregistered_command() {
	app := App{
		root: Command{
			name:     'reflux'
			commands: [Command{ name: 'run' }]
		}
	}
	app.execute(['run']) or {
		assert err.msg() == 'no handler for run'
		return
	}
	assert false
}
