module cmdline

import os

const colored = Style{
	terminal: Terminal{
		color: true
	}
}

// sgr is text in one of reflux_theme's codes, as Style writes it.
fn sgr(code string, text string) string {
	return '\x1b[${code}m${text}\x1b[0m'
}

// kinetics_tree is a command tree with every part of help: a group, a
// command with required, optional and variadic arguments, flags with and
// without a short name, examples and a section.
fn kinetics_tree() App {
	return App{
		root: Command{
			name:     'reflux'
			usage:    'reflux <command>'
			commands: [
				Command{
					name:    'solve'
					summary: 'Solve a flowsheet.'
				},
				Command{
					name:     'kinetics'
					summary:  'Check reactions.'
					commands: [
						Command{
							name:        'check'
							usage:       'reflux kinetics check FILE [T P NAME=y ...]'
							description: 'Check reaction definitions.'
							arguments:   [
								Argument{
									name:        'FILE'
									description: 'Kinetics file.'
									required:    true
								},
								Argument{
									name:        'T'
									description: 'Temperature.'
								},
								Argument{
									name:        'NAME=y'
									description: 'Mole fractions.'
									variadic:    true
								},
							]
							flags:       [
								Flag{
									name:        'json'
									alias:       'j'
									description: 'Write JSON.'
									takes_value: true
									value_name:  'OUT'
								},
								Flag{
									name:        'strict'
									description: 'Treat warnings as errors.'
								},
							]
							examples:    [
								'# Check a file',
								'reflux kinetics check rates.toml --strict',
								'',
								"reflux kinetics check rates.toml '900 K' CH4=1",
							]
							sections:    [
								Section{
									title: 'Exit status'
									text:  '0 when it read.\n1 when not.'
								},
							]
						},
					]
				},
			]
		}
	}
}

const check_help = "Check reaction definitions.

USAGE

  reflux kinetics check FILE [T P NAME=y ...]

ARGUMENTS

  FILE         Kinetics file.
  [T]          Temperature.
  [NAME=y...]  Mole fractions.

OPTIONS

  -j, --json OUT  Write JSON.
      --strict    Treat warnings as errors.

EXAMPLES

  # Check a file
  reflux kinetics check rates.toml --strict

  reflux kinetics check rates.toml '900 K' CH4=1

EXIT STATUS

  0 when it read.
  1 when not.
"

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
	expected := 'A simulator.\n\nUSAGE\n\n  reflux <command>\n\nCOMMANDS\n\n  solve                Solve a flowsheet.\n  version (--version)  Print the version.\n'
	assert help(command) == expected
	assert !expected.contains('\t')
}

// Every way of asking for one command's help gives the same help, and none
// needs the command's arguments, or minds them.
fn test_nested_command_help_contains_all_sections() {
	app := kinetics_tree()
	for args in [
		['kinetics', 'check'],
		['kinetics', 'check', '--help'],
		['kinetics', 'check', '-h'],
		['kinetics', 'check', 'rates.toml', '--strict', '--help'],
		['help', 'kinetics', 'check'],
		['kinetics', 'help', 'check'],
	] {
		got := app.help(args, Style{}) or { panic('${args}: ${err}') }
		assert got == check_help, args.str()
	}
}

// The first `help` is the request, and a second is the command it is for.
fn test_help_on_a_group_and_on_help() {
	group := kinetics_tree().help(['kinetics', 'help'], Style{}) or { panic(err) }
	assert group == 'COMMANDS\n\n  check\n'
	app := App{
		root: Command{
			name:     'reflux'
			usage:    'reflux <command>'
			commands: [Command{
				name:  'help'
				usage: 'reflux help [COMMAND]'
			}]
		}
	}
	root := app.help(['help'], Style{}) or { panic(err) }
	assert root.starts_with('USAGE\n\n  reflux <command>\n')
	on_help := app.help(['help', 'help'], Style{}) or { panic(err) }
	assert on_help == 'USAGE\n\n  reflux help [COMMAND]\n'
}

fn test_help_reports_an_unknown_command_with_suggestions() {
	app := kinetics_tree()
	app.help(['solv', '--help'], Style{}) or {
		assert err is ParseError
		assert err.msg() == 'no command solv (did you mean solve?)'
		return
	}
	assert false
}

fn test_wants_help() {
	root := kinetics_tree().root
	for args in [['help'], ['--help'], ['-h'], ['solve', 'a.toml', '--help'], ['kinetics', 'help'],
		['solv', '--help'], ['kinetics', 'check', 'f', '-h']] {
		assert wants_help(root, args), args.str()
	}
	// `solve help` is a file called help; after `--` --help is a value.
	for args in [[]string{}, ['solve', 'help'], ['solve', 'a.toml', '--', '--help'],
		['kinetics', 'check', 'help'], ['--version']] {
		assert !wants_help(root, args), args.str()
	}
}

// The colours land on the parts they are for: the label of a heading, the
// program, commands, flags and placeholders in a usage line, and in an
// example a comment, while a file name, a quoted value and CH4=1 keep the
// terminal's colour.
fn test_colored_help_highlights_each_part() {
	got := kinetics_tree().help(['kinetics', 'check'], colored) or { panic(err) }
	t := reflux_theme
	usage := '${sgr(t.program, 'reflux')} ${sgr(t.command, 'kinetics')} ${sgr(t.command, 'check')} ${sgr(t.placeholder, 'FILE')} ${sgr(t.muted, '[')}${sgr(t.placeholder, 'T')} ${sgr(t.placeholder, 'P')} ${sgr(t.placeholder, 'NAME=y')} ${sgr(t.muted, '...')}${sgr(t.muted, ']')}'
	assert got.starts_with('Check reaction definitions.\n\n${sgr(t.heading, ' USAGE ')}\n\n  ${usage}\n'), got
	for title in ['ARGUMENTS', 'OPTIONS', 'EXAMPLES', 'EXIT STATUS'] {
		assert got.contains('\n${sgr(t.heading, ' ${title} ')}\n\n'), title
	}
	// Alignment is the plain text's: the escapes take no columns.
	assert got.contains('  ${sgr(t.muted, '[')}${sgr(t.placeholder, 'T')}${sgr(t.muted, ']')}          Temperature.\n')
	assert got.contains('  ${sgr(t.flag, '-j')}${sgr(t.muted, ',')} ${sgr(t.flag, '--json')} ${sgr(t.placeholder, 'OUT')}  Write JSON.\n')
	assert got.contains('      ${sgr(t.flag, '--strict')}    Treat warnings as errors.\n')
	assert got.contains('  ${sgr(t.comment, '# Check a file')}\n')
	program := '${sgr(t.program, 'reflux')} ${sgr(t.command, 'kinetics')} ${sgr(t.command, 'check')}'
	assert got.contains('  ${program} rates.toml ${sgr(t.flag, '--strict')}\n')
	assert got.contains("  ${program} rates.toml '900 K' CH4=1\n")
	// The plain help is the coloured one without its escapes and labels' spaces.
	assert strip_escapes(got).replace(' USAGE ', 'USAGE').replace(' ARGUMENTS ', 'ARGUMENTS').replace(' OPTIONS ',
		'OPTIONS').replace(' EXAMPLES ', 'EXAMPLES').replace(' EXIT STATUS ', 'EXIT STATUS') == check_help
}

// strip_escapes removes every `ESC [ ... m`.
fn strip_escapes(s string) string {
	mut out := []u8{}
	mut i := 0
	for i < s.len {
		if s[i] == 0x1b {
			for i < s.len && s[i] != `m` {
				i++
			}
		} else {
			out << s[i]
		}
		i++
	}
	return out.bytestr()
}

fn test_plain_styles_write_no_escapes() {
	app := kinetics_tree()
	for style in [Style{}, Style{
		terminal: Terminal{
			interactive: true
		}
	}] {
		for args in [['help'], ['kinetics'], ['kinetics', 'check']] {
			text := app.help(args, style) or { panic(err) }
			assert !text.contains('\x1b'), args.str()
		}
		assert !Diagnostic{
			kind:    .error
			message: 'm'
			hint:    'h'
			usage:   'reflux solve FILE'
		}.render(style).contains('\x1b')
	}
}

// A label's width is what it shows, so T₁ (two characters in four bytes)
// lines up with FILE.
fn test_rows_align_by_visible_width() {
	text := help(Command{
		name:      'x'
		arguments: [
			Argument{
				name:        'FILE'
				required:    true
				description: 'a'
			},
			Argument{
				name:        'T₁'
				required:    true
				description: 'b'
			},
		]
	})
	assert text == 'ARGUMENTS\n\n  FILE  a\n  T₁    b\n'
}

fn test_terminal_honours_no_color_and_clicolor_force() {
	old_no_color := os.getenv_opt('NO_COLOR')
	old_force := os.getenv_opt('CLICOLOR_FORCE')
	defer {
		restore_env('NO_COLOR', old_no_color)
		restore_env('CLICOLOR_FORCE', old_force)
	}
	os.unsetenv('NO_COLOR')
	os.setenv('CLICOLOR_FORCE', '1', true)
	assert terminal().color && stderr_terminal().color
	os.setenv('NO_COLOR', '1', true)
	assert !terminal().color && !stderr_terminal().color
}

fn restore_env(name string, value ?string) {
	if v := value {
		os.setenv(name, v, true)
	} else {
		os.unsetenv(name)
	}
}

fn solve_tree() Command {
	return Command{
		name:     'reflux'
		commands: [
			Command{
				name:      'solve'
				aliases:   ['run']
				usage:     'reflux solve FILE [--json OUT] [--trace]'
				arguments: [Argument{
					name:     'FILE'
					required: true
				}]
				flags:     [
					Flag{
						name:        'json'
						alias:       'j'
						takes_value: true
					},
					Flag{
						name:        'csv'
						takes_value: true
					},
					Flag{
						name: 'trace'
					},
				]
			},
			Command{
				name:      'demo'
				arguments: [Argument{
					name: 'T'
				}, Argument{
					name: 'P'
				}]
			},
			Command{
				name:    'version'
				aliases: ['--version']
			},
		]
	}
}

fn test_parser_routes_and_normalizes_flags() {
	invocation := parse(solve_tree(), ['solve', 'plant.toml', '-j', 'result.json', '--trace']) or {
		panic(err)
	}
	assert invocation.command == ['solve']
	assert invocation.arguments == ['plant.toml']
	assert invocation.values['json'] == 'result.json'
	assert invocation.switches['trace']
	inline := parse(solve_tree(), ['run', 'input.toml', '-j=output.json', '--csv=']) or {
		panic(err)
	}
	assert inline.command == ['solve']
	assert inline.values['json'] == 'output.json'
	// An empty value is a value, which the application can refuse.
	assert inline.values['csv'] == ''
	version := parse(solve_tree(), ['--version']) or { panic(err) }
	assert version.command == ['version']
}

// A negative number and a lone dash are values, positional or a flag's; a
// flag where a value should be is a value forgotten.
fn test_parser_takes_negative_numbers_and_dashes_as_values() {
	demo := parse(solve_tree(), ['demo', '300', '-1']) or { panic(err) }
	assert demo.arguments == ['300', '-1']
	half := parse(solve_tree(), ['demo', '-.5', '-']) or { panic(err) }
	assert half.arguments == ['-.5', '-']
	value := parse(solve_tree(), ['solve', 'a.toml', '--json', '-2']) or { panic(err) }
	assert value.values['json'] == '-2'
	assert_parse_error(solve_tree(), ['solve', 'a.toml', '--csv', '--trace'], '--csv needs a value')
	assert_parse_error(solve_tree(), ['solve', 'a.toml', '--csv', '--'], '--csv needs a value')
	assert_parse_error(solve_tree(), ['demo', '-x'], 'no option -x')
}

fn test_parser_rejects_invalid_forms() {
	tree := solve_tree()
	assert_parse_error(tree, ['solve', '--json'], '--json needs a value')
	assert_parse_error(tree, ['solve', 'input.toml', '--trace', '--trace'], '--trace is given twice')
	assert_parse_error(tree, ['solve', 'input.toml', '--json=', '--json=next.json'],
		'--json is given twice')
	assert_parse_error(tree, ['solve', 'input.toml', '--trac'], 'no option --trac (did you mean --trace?)')
	assert_parse_error(tree, ['solve', 'input.toml', '--jsn'], 'no option --jsn (did you mean --json, --csv?)')
	assert_parse_error(tree, ['solve', 'input.toml', '--trace='], '--trace does not take a value')
	assert_parse_error(tree, ['solve', 'input.toml', '--trace=yes'], '--trace does not take a value')
	// A long name is --name and a short one -alias, not the other way round.
	assert_parse_error(tree, ['solve', 'input.toml', '-trace'], 'no option -trace (did you mean --trace?)')
	assert_parse_error(tree, ['solve', 'input.toml', '--j', 'x'], 'no option --j')
	assert_parse_error(tree, ['solve'], 'solve needs at least 1 argument')
	assert_parse_error(tree, ['solve', 'a', 'b'], 'solve takes 1 argument, not 2')
	assert_parse_error(tree, ['version', 'x'], 'version takes no positional arguments')
}

// Every failure is a ParseError, an unknown option too, and carries its
// command's usage for a diagnostic to show.
fn test_parse_errors_are_typed_and_carry_the_usage() {
	parse(solve_tree(), ['solve', 'a.toml', '--zzz']) or {
		assert err is ParseError
		e := err as ParseError
		assert e.usage == 'reflux solve FILE [--json OUT] [--trace]'
		assert e.diagnostic().render(Style{}) == 'Error: no option --zzz\nUsage: reflux solve FILE [--json OUT] [--trace]\n'
		return
	}
	assert false
}

fn test_parser_needs_a_command_in_a_group() {
	tree := kinetics_tree().root
	assert_parse_error(tree, ['kinetics'], 'kinetics needs a command: one of check')
	assert_parse_error(tree, []string{}, 'reflux needs a command: one of solve, kinetics')
	assert_parse_error(tree, ['--', 'solve'], 'no command --')
	assert_parse_error(tree, ['solv'], 'no command solv (did you mean solve?)')
	// zzz is three edits from every command: too far to guess at.
	assert_parse_error(tree, ['zzz'], 'no command zzz')
}

// A command is named exactly: `--run` is not the alias `run`.
fn test_parser_matches_aliases_exactly() {
	assert_parse_error(solve_tree(), ['--run', 'a.toml'], 'no command --run (did you mean run?)')
	assert_parse_error(solve_tree(), ['-version'], 'no command -version (did you mean version, --version?)')
}

fn test_parser_takes_the_rest_for_a_variadic_argument() {
	tree := Command{
		name:     'kinetics'
		commands: [
			Command{
				name:      'explain'
				arguments: [
					Argument{
						name:     'FILE'
						required: true
					},
					Argument{
						name:     'T'
						required: true
					},
					Argument{
						name:     'NAME=y'
						required: true
						variadic: true
					},
				]
			},
		]
	}
	got := parse(tree, ['explain', 'f.toml', '900 K', 'CH4=0.5', 'H2O=0.5']) or { panic(err) }
	assert got.arguments == ['f.toml', '900 K', 'CH4=0.5', 'H2O=0.5']
	assert_parse_error(tree, ['explain', 'f.toml', '900 K'], 'explain needs at least 3 arguments')
}

fn assert_parse_error(command Command, args []string, expected string) {
	parse(command, args) or {
		assert err is ParseError, args.str()
		assert err.msg() == expected, args.str()
		return
	}
	panic('expected parser error: ${args.join(' ')}')
}

fn test_suggestions_are_near_and_few() {
	assert did_you_mean('zzz', ['solve', 'fit']) == ''
	assert did_you_mean('-q', ['--json', '--csv', '--trace']) == ''
	assert did_you_mean('--jsn', ['--json', '--csv', '--trace']) == ' (did you mean --json, --csv?)'
	assert did_you_mean('Solve', ['solve']) == ' (did you mean solve?)'
	assert did_you_mean('ab', ['a', 'b', 'abc', 'xb']) == ' (did you mean a, b, abc?)'
}

fn assert_invalid(app App, expected string) {
	app.validate() or {
		assert err.msg() == expected
		return
	}
	panic('expected invalid metadata: ${expected}')
}

fn test_validation_rejects_ambiguous_metadata() {
	assert_invalid(App{
		root: Command{
			name:  'reflux'
			flags: [Flag{
				name: 'json'
			}, Flag{
				name:  'output'
				alias: 'json'
			}]
		}
	}, 'reflux: duplicate flag json')
	assert_invalid(App{
		root: Command{
			name:     'reflux'
			commands: [Command{
				name:    'solve'
				aliases: ['run', '--run']
			}]
		}
	}, 'reflux solve: duplicate command name --run')
	assert_invalid(App{
		root: Command{
			name:     'reflux'
			commands: [Command{
				name: '--solve'
			}]
		}
	}, 'reflux --solve: command name --solve starts with a dash')
	assert_invalid(App{
		root: Command{
			name:     'reflux'
			commands: [
				Command{
					name:      'x'
					arguments: [Argument{
						name:     'A'
						variadic: true
					}, Argument{
						name: 'B'
					}]
				},
			]
		}
	}, 'reflux x: variadic argument A is not the last')
	// An error names the command by its path.
	assert_invalid(App{
		root: Command{
			name:     'reflux'
			commands: [
				Command{
					name:     'kinetics'
					commands: [
						Command{
							name:      'check'
							arguments: [Argument{
								name: 'FILE'
							}, Argument{
								name: 'FILE'
							}]
						},
					]
				},
			]
		}
	}, 'reflux kinetics check: duplicate argument FILE')
}

// A command's own name with dashes is an alias wherever it comes in the list.
fn test_validation_takes_a_dashed_name_anywhere() {
	App{
		root: Command{
			name:     'reflux'
			commands: [Command{
				name:    'version'
				aliases: ['v', '--version']
			}]
		}
	}.validate() or { panic(err) }
}

fn test_diagnostic_respects_style() {
	plain := Diagnostic{
		kind:    .error
		message: 'bad input'
		hint:    'check the file'
		usage:   'reflux solve FILE'
	}.render(Style{})
	assert plain == 'Error: bad input\nHint: check the file\nUsage: reflux solve FILE\n'

	t := reflux_theme
	warning := Diagnostic{
		kind:    .warning
		message: 'using a default'
	}.render(colored)
	assert warning == '${sgr(t.warning, ' WARNING ')} using a default\n'
	failure := Diagnostic{
		kind:    .error
		message: 'bad input'
		usage:   'reflux solve FILE'
	}.render(colored)
	assert failure == '${sgr(t.error, ' ERROR ')} bad input\n${sgr(t.muted, 'Usage:')} ${sgr(t.program, 'reflux')} ${sgr(t.command, 'solve')} ${sgr(t.placeholder, 'FILE')}\n'
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
					arguments: [Argument{
						name:     'FILE'
						required: true
					}]
				},
				Command{
					name:    'version'
					aliases: ['--version']
				},
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
	for args in [[]string{}, ['run', 'input.toml', '--help']] {
		outcome := app.run(args, Style{}) or { panic(err) }
		assert outcome.kind == .help && outcome.code == 0, args.str()
	}
	version := app.run(['--version'], Style{}) or { panic(err) }
	assert version.kind == .version
	assert version.text == 'reflux 0.3.0'
	executed := app.run(['run', 'input.toml'], Style{}) or { panic(err) }
	assert executed.kind == .executed
	assert executed.code == 0
}

fn test_app_rejects_unregistered_command() {
	app := App{
		root: Command{
			name:     'reflux'
			commands: [Command{
				name: 'run'
			}]
		}
	}
	app.execute(['run']) or {
		assert err.msg() == 'no handler for run'
		return
	}
	assert false
}

fn test_a_handler_cannot_be_for_a_group() {
	tree := kinetics_tree()
	assert_invalid(App{
		root:     tree.root
		commands: [RegisteredCommand{
			path:    ['kinetics']
			handler: TestHandler{}
		}]
	}, 'handler path kinetics names a group of commands, which parse never returns')
	tree.execute(['kinetics']) or {
		assert err.msg() == 'kinetics needs a command: one of check'
		return
	}
	assert false
}
