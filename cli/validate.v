module cli

// validate checks command metadata before the application starts accepting
// user input.
pub fn (app App) validate() ! {
	validate_command(app.root, 'root')!
	mut paths := [][]string{}
	for registered in app.commands {
		if registered.path.len == 0 {
			return error('a handler must have a command path')
		}
		canonical := canonical_path(app.root, registered.path)!
		if canonical != registered.path {
			return error('handler path ${registered.path.join(' ')} is not canonical; use ${canonical.join(' ')}')
		}
		if paths.any(it == canonical) {
			return error('duplicate handler for ${canonical.join(' ')}')
		}
		paths << canonical
	}
}

fn canonical_path(root Command, path []string) ![]string {
	mut command := root
	mut canonical := []string{}
	for token in path {
		index := command_index(command.commands, token)
		if index < 0 {
			return error('handler path ${path.join(' ')} names no command')
		}
		command = command.commands[index]
		canonical << command.name
	}
	return canonical
}

fn validate_command(command Command, location string) ! {
	if command.name == '' {
		return error('${location}: command name is empty')
	}
	if command.commands.len > 0 && command.arguments.len > 0 {
		return error('${location}: command groups cannot define positional arguments')
	}
	if command.commands.len > 0 && command.flags.len > 0 {
		return error('${location}: command groups cannot define flags')
	}
	mut names := [command.name]
	names << command.aliases
	for i, name in names {
		if name == '' {
			return error('${location}: command alias ${i} is empty')
		}
		if name.contains(' ') {
			return error('${location}: command name ${name} contains a space')
		}
		for previous in names[..i] {
			if i == 1 && name.starts_with('-') && normalize_name(previous) == normalize_name(name) {
				continue
			}
			if normalize_name(previous) == normalize_name(name) {
				return error('${location}: duplicate command name ${name}')
			}
		}
	}
	for i, child in command.commands {
		validate_command(child, '${location} ${command.name}[${i}]')!
		for previous in command.commands[..i] {
			if command_conflicts(previous, child) {
				return error('${location}: duplicate child command ${child.name}')
			}
		}
	}
	validate_arguments(command, location)!
	validate_flags(command, location)!
}

fn validate_arguments(command Command, location string) ! {
	mut names := []string{}
	mut optional_seen := false
	for argument in command.arguments {
		if argument.name == '' {
			return error('${location}: argument name is empty')
		}
		if argument.name.contains(' ') || argument.name.starts_with('-') {
			return error('${location}: invalid argument name ${argument.name}')
		}
		if argument.name in names {
			return error('${location}: duplicate argument ${argument.name}')
		}
		if !argument.required {
			optional_seen = true
		} else if optional_seen {
			return error('${location}: required argument ${argument.name} follows an optional argument')
		}
		names << argument.name
	}
}

fn validate_flags(command Command, location string) ! {
	mut names := []string{}
	for flag in command.flags {
		if flag.name == '' {
			return error('${location}: flag name is empty')
		}
		if flag.name.starts_with('-') || flag.name.contains('=') || flag.name.contains(' ') {
			return error('${location}: invalid flag name ${flag.name}')
		}
		if flag.name in names {
			return error('${location}: duplicate flag ${flag.name}')
		}
		names << flag.name
		if flag.alias != '' {
			if flag.alias.starts_with('-') || flag.alias.contains('=') || flag.alias.contains(' ') {
				return error('${location}: invalid flag alias ${flag.alias}')
			}
			if flag.alias in names {
				return error('${location}: duplicate flag ${flag.alias}')
			}
			names << flag.alias
		}
	}
}

fn command_conflicts(left Command, right Command) bool {
	if normalize_name(left.name) == normalize_name(right.name) {
		return true
	}
	for name in left.aliases {
		if normalize_name(name) == normalize_name(right.name) || right.aliases.any(normalize_name(it) == normalize_name(name)) {
			return true
		}
	}
	for name in right.aliases {
		if normalize_name(left.name) == normalize_name(name) {
			return true
		}
	}
	return false
}
