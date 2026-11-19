module cmdline

import math
import strings

// Invocation is the result of parsing one command line. It contains no
// application-specific values or side effects.
pub struct Invocation {
pub:
	command []string
pub mut:
	arguments []string
	values    map[string]string
	switches  map[string]bool
}

// parse parses a command line against the command tree: its commands first,
// then the last command's flags and positional arguments, in any order. A
// flag is --name or -alias, and its value follows it or an `=`. After `--`
// everything is positional, and so is a negative number, -1 or -.5, and a
// lone `-`. It returns data or a ParseError, and prints and opens nothing.
pub fn parse(app Command, args []string) !Invocation {
	mut command := app
	mut command_path := []string{}
	mut index := 0
	for command.commands.len > 0 {
		if index >= args.len {
			names := command.commands.map(it.name).join(', ')
			return parse_error('${display_name(app, command_path)} needs a command: one of ${names}',
				command.usage)
		}
		token := args[index]
		child := command_index(command.commands, token)
		if child < 0 {
			return parse_error('no command ${token}${command_suggestion(command.commands, token)}', command.usage)
		}
		command = command.commands[child]
		command_path << command.name
		index++
	}
	name := display_name(app, command_path)

	mut invocation := Invocation{
		command: command_path
	}
	mut positional := []string{}
	mut end_of_options := false
	for index < args.len {
		token := args[index]
		index++
		if end_of_options || !is_option(token) {
			positional << token
			continue
		}
		if token == '--' {
			end_of_options = true
			continue
		}
		flag := find_flag(command.flags, token) or {
			return parse_error('no option ${token}${flag_suggestion(command.flags, token)}',
				command.usage)
		}
		long := '--${flag.name}'
		inline := token.contains('=')
		if flag.takes_value {
			mut value := token.all_after_first('=')
			if !inline {
				// A flag after a flag is a value forgotten, not a value.
				if index >= args.len || is_option(args[index]) {
					return parse_error('${long} needs a value', command.usage)
				}
				value = args[index]
				index++
			}
			if flag.name in invocation.values {
				return parse_error('${long} is given twice', command.usage)
			}
			invocation.values[flag.name] = value
		} else {
			if inline {
				return parse_error('${long} does not take a value', command.usage)
			}
			if invocation.switches[flag.name] {
				return parse_error('${long} is given twice', command.usage)
			}
			invocation.switches[flag.name] = true
		}
	}

	required := command.arguments.filter(it.required).len
	variadic := command.arguments.len > 0 && command.arguments.last().variadic
	if command.arguments.len == 0 && positional.len > 0 {
		return parse_error('${name} takes no positional arguments', command.usage)
	}
	if positional.len < required {
		return parse_error('${name} needs at least ${required} ${arguments_word(required)}',
			command.usage)
	}
	if !variadic && positional.len > command.arguments.len {
		return parse_error('${name} takes ${command.arguments.len} ${arguments_word(command.arguments.len)}, not ${positional.len}',
			command.usage)
	}
	invocation.arguments = positional
	return invocation
}

// display_name is how an error names a command: its path below the root,
// `kinetics check`, or the root's name.
fn display_name(root Command, path []string) string {
	return if path.len > 0 { path.join(' ') } else { root.name }
}

fn arguments_word(n int) string {
	return if n == 1 { 'argument' } else { 'arguments' }
}

// is_option is whether a token is a flag or `--`, rather than a value.
fn is_option(token string) bool {
	return token.len > 1 && token[0] == `-` && !is_negative_number(token)
}

// find_flag finds the flag a token names, as --name or -alias, with or
// without an `=value`.
fn find_flag(flags []Flag, token string) ?Flag {
	name := token.all_before('=')
	for flag in flags {
		if name == '--${flag.name}' || (flag.alias != '' && name == '-${flag.alias}') {
			return flag
		}
	}
	return none
}

// command_index finds the command a token names exactly, by its name or an
// alias, so that `--run` is not the command `run`.
fn command_index(commands []Command, token string) int {
	for i, child in commands {
		if child.name == token || token in child.aliases {
			return i
		}
	}
	return -1
}

fn normalize_name(name string) string {
	return name.trim_left('-')
}

fn command_suggestion(commands []Command, token string) string {
	mut names := []string{}
	for command in commands {
		names << command.name
		names << command.aliases
	}
	return did_you_mean(token, names)
}

// flag_suggestion compares a flag with the long names only: a one-letter
// flag is one edit from every other, which says nothing.
fn flag_suggestion(flags []Flag, token string) string {
	return did_you_mean(token.all_before('='), flags.map('--${it.name}'))
}

// did_you_mean suggests up to three names within two edits of a token,
// nearest first, case aside: `zzz` is three from `fit`, too far to be a
// slip of the keyboard.
fn did_you_mean(token string, names []string) string {
	mut close := []string{}
	for distance in 0 .. 3 {
		for name in names {
			if name !in close
				&& strings.levenshtein_distance(token.to_lower(), name.to_lower()) == distance {
				close << name
			}
		}
	}
	if close.len == 0 {
		return ''
	}
	return ' (did you mean ${close[..math.min(3, close.len)].join(', ')}?)'
}
