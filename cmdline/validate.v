module cmdline

// validate checks command metadata before the application starts accepting
// user input. An error names the command by its path, `reflux kinetics`.
pub fn (app App) validate() ! {
	validate_command(app.root, if app.root.name == '' { 'root' } else { app.root.name })!
	mut paths := [][]string{}
	for registered in app.commands {
		if registered.path.len == 0 {
			return error('a handler must have a command path')
		}
		canonical, command := canonical_path(app.root, registered.path)!
		if canonical != registered.path {
			return error('handler path ${registered.path.join(' ')} is not canonical; use ${canonical.join(' ')}')
		}
		if command.commands.len > 0 {
			return error('handler path ${canonical.join(' ')} names a group of commands, which parse never returns')
		}
		if paths.any(it == canonical) {
			return error('duplicate handler for ${canonical.join(' ')}')
		}
		paths << canonical
	}
}

fn canonical_path(root Command, path []string) !([]string, Command) {
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
	return canonical, command
}

fn validate_command(command Command, location string) ! {
	if command.name == '' {
		return error('${location}: command name is empty')
	}
	if command.name.starts_with('-') {
		return error('${location}: command name ${command.name} starts with a dash')
	}
	if command.commands.len > 0 && command.arguments.len > 0 {
		return error('${location}: command groups cannot define positional arguments')
	}
	if command.commands.len > 0 && command.flags.len > 0 {
		return error('${location}: command groups cannot define flags')
	}
	// An alias can be the name with dashes, `--version` for `version`; any
	// other two names that differ only in dashes are one name twice.
	mut keys := [command.name]
	for i, alias in command.aliases {
		if alias == '' {
			return error('${location}: command alias ${i + 1} is empty')
		}
		if alias.contains(' ') {
			return error('${location}: command name ${alias} contains a space')
		}
		dashed_name := alias.starts_with('-') && normalize_name(alias) == command.name
		if alias == command.name || command.aliases[..i].contains(alias)
			|| (!dashed_name && normalize_name(alias) in keys) {
			return error('${location}: duplicate command name ${alias}')
		}
		if !dashed_name {
			keys << normalize_name(alias)
		}
	}
	for i, child in command.commands {
		child_location := if child.name == '' {
			'${location} [${i}]'
		} else {
			'${location} ${child.name}'
		}
		validate_command(child, child_location)!
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
	for i, argument in command.arguments {
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
		if argument.variadic && i != command.arguments.len - 1 {
			return error('${location}: variadic argument ${argument.name} is not the last')
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
	mut left_names := [normalize_name(left.name)]
	for alias in left.aliases {
		left_names << normalize_name(alias)
	}
	if normalize_name(right.name) in left_names {
		return true
	}
	return right.aliases.any(normalize_name(it) in left_names)
}
