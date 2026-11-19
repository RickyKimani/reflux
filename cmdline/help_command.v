module cmdline

// wants_help is whether args ask for help: --help or -h anywhere before
// `--`, or the word help where a command is expected, as in `reflux help
// solve` or `reflux kinetics help`. Past a command that takes arguments,
// help is an argument like any other.
pub fn wants_help(root Command, args []string) bool {
	mut command := root
	mut in_group := root.commands.len > 0
	for token in args {
		if token == '--' {
			return false
		}
		if token == '--help' || token == '-h' {
			return true
		}
		if !in_group {
			continue
		}
		if token == 'help' {
			return true
		}
		index := command_index(command.commands, token)
		if index < 0 {
			// Not a command: App.help says so, if a --help follows.
			in_group = false
			continue
		}
		command = command.commands[index]
		in_group = command.commands.len > 0
	}
	return false
}

// help renders the help a help request asks for. It doesn't call parse,
// because help must work without a command's required arguments: it follows
// the commands the request names, and stops at the first word that isn't
// one, where a command takes arguments. Where it expects a command, a word
// that isn't one is an error, with the commands it is near.
pub fn (app App) help(args []string, style Style) !string {
	app.validate()!
	mut command := app.root
	mut skipped_help := false
	for token in args {
		if token.starts_with('-') || command.commands.len == 0 {
			break
		}
		// The first help is the request; a second is the command, for help on
		// help.
		if token == 'help' && !skipped_help {
			skipped_help = true
			continue
		}
		index := command_index(command.commands, token)
		if index < 0 {
			return parse_error('no command ${token}${command_suggestion(command.commands, token)}', command.usage)
		}
		command = command.commands[index]
	}
	return render_help(command, style)
}
