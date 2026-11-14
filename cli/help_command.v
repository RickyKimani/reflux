module cli

// help renders the command selected by a help-style argument. It deliberately
// does not call parse because a help request must work without required args.
pub fn (app App) help(args []string, style Style) !string {
	app.validate()!
	mut path := []string{}
	if args.len > 0 && args[0] == 'help' {
		path = args[1..].clone()
	} else {
		for arg in args {
			if arg in ['--help', '-h'] {
				break
			}
			if arg.starts_with('-') {
				break
			}
			path << arg
		}
	}
	command := find_command(app.root, path)!
	return render_help(command, style)
}

fn find_command(root Command, path []string) !Command {
	mut command := root
	for token in path {
		index := command_index(command.commands, token)
		if index < 0 {
			return error('unknown command ${token}${command_suggestion(command.commands, token)}')
		}
		command = command.commands[index]
	}
	return command
}
