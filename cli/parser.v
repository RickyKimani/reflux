module cli

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

// parse parses commands, positional arguments and flags from args.
pub fn parse(app Command, args []string) !Invocation {
	mut command := app
	mut command_path := []string{}
	mut index := 0
	for index < args.len {
		token := args[index]
		if token.starts_with('-') && command_index(command.commands, token) < 0 {
			break
		}
		child := command_index(command.commands, token)
		if child < 0 {
			if command.commands.len > 0 {
				return parse_error('unknown command ${token}${command_suggestion(command.commands, token)}')
			}
			break
		}
		command = command.commands[child]
		command_path << command.name
		index++
	}

	mut invocation := Invocation{
		command:  command_path
		values:   map[string]string{}
		switches: map[string]bool{}
	}
	mut seen_values := map[string]bool{}
	mut positional := []string{}
	mut end_of_options := false
	for index < args.len {
		token := args[index]
		if token == '--' && !end_of_options {
			end_of_options = true
			index++
			continue
		}
		if end_of_options || !token.starts_with('-') {
			positional << token
			index++
			continue
		}
		name, inline_value := split_flag(token)
		has_inline_value := token.contains('=')
		index_of_flag := find_flag_index(command.flags, name)
		if index_of_flag < 0 {
			return error('unknown option ${token}${flag_suggestion(command.flags, name)}')
		}
		flag := command.flags[index_of_flag]
		if flag.takes_value {
			mut value := inline_value
			if !has_inline_value {
				if index + 1 >= args.len {
					return parse_error('--${flag.name} needs a value')
				}
				value = args[index + 1]
				index++
			}
			if seen_values[flag.name] {
				return parse_error('--${flag.name} is given twice')
			}
			invocation.values[flag.name] = value
			seen_values[flag.name] = true
		} else {
			if inline_value != '' {
				return parse_error('--${flag.name} does not take a value')
			}
			if invocation.switches[flag.name] {
				return parse_error('--${flag.name} is given twice')
			}
			invocation.switches[flag.name] = true
		}
		index++
	}
	if command.commands.len > 0 && command_path.len == 0 && positional.len > 0 {
		return parse_error('unknown command ${positional[0]}${command_suggestion(command.commands, positional[0])}')
	}
	if command.arguments.len > 0 {
		required := command.arguments.filter(it.required).len
		if positional.len < required {
			return parse_error('${command.name} needs at least ${required} argument${if required == 1 {
				''
			} else {
				's'
			}}')
		}
		if positional.len > command.arguments.len {
			return parse_error('${command.name} takes ${command.arguments.len} argument${if command.arguments.len == 1 {
				''
			} else {
				's'
			}}, not ${positional.len}')
		}
	} else if positional.len > 0 {
		return parse_error('${command.name} takes no positional arguments')
	}
	invocation.arguments = positional
	return invocation
}

fn split_flag(token string) (string, string) {
	if token.starts_with('--') {
		name := token[2..]
		if name.contains('=') {
			return name.all_before('='), name.all_after_first('=')
		}
		return name, ''
	}
	name := token[1..]
	if name.contains('=') {
		return name.all_before('='), name.all_after_first('=')
	}
	return name, ''
}

fn command_index(commands []Command, token string) int {
	for i, child in commands {
		if child.name == token || child.aliases.any(normalize_name(it) == normalize_name(token)) {
			return i
		}
	}
	return -1
}

fn normalize_name(name string) string {
	return name.trim_left('-')
}

fn command_suggestion(commands []Command, token string) string {
	mut best := ''
	mut distance := 4
	for command in commands {
		mut candidates := [command.name]
		candidates << command.aliases
		for candidate in candidates {
			d := edit_distance(normalize_name(token), normalize_name(candidate))
			if d < distance {
				best = candidate
				distance = d
			}
		}
	}
	return if best == '' { '' } else { ' (did you mean ${best}?)' }
}

fn flag_suggestion(flags []Flag, token string) string {
	mut best := ''
	mut distance := 4
	for flag in flags {
		for candidate in [flag.name, flag.alias] {
			if candidate == '' {
				continue
			}
			d := edit_distance(token, candidate)
			if d < distance {
				best = flag.name
				distance = d
			}
		}
	}
	return if best == '' { '' } else { ' (did you mean --${best}?)' }
}

fn edit_distance(a string, b string) int {
	mut previous := []int{len: b.runes().len + 1}
	for k in 0 .. previous.len {
		previous[k] = k
	}
	for i, left in a.runes() {
		mut current := []int{len: previous.len}
		current[0] = i + 1
		for j, right in b.runes() {
			cost := if left == right { 0 } else { 1 }
			current[j + 1] = min3(current[j] + 1, previous[j + 1] + 1, previous[j] + cost)
		}
		previous = current.clone()
	}
	return previous.last()
}

fn min3(a int, b int, c int) int {
	return if a < b {
		if a < c { a } else { c }
	} else if b < c {
		b
	} else {
		c
	}
}

fn find_flag_index(flags []Flag, name string) int {
	for i, flag in flags {
		if flag.name == name || flag.alias == name {
			return i
		}
	}
	return -1
}
