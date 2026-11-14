module cli

// App combines a command tree with handlers. It is the framework boundary
// used by an application's main function; it does not terminate the process.
pub struct App {
pub:
	root     Command
	commands []RegisteredCommand
	version  string
}

pub fn (app App) execute(args []string) !int {
	app.validate()!
	invocation := parse(app.root, args)!
	return dispatch(invocation, app.commands)
}

// run handles framework-owned help and returns presentation data without
// writing to a stream. Application code remains responsible for printing it.
pub fn (app App) run(args []string, style Style) !Outcome {
	app.validate()!
	if args.len == 1 && (args[0] == 'version' || args[0] == '--version') {
		return Outcome{
			kind: .version
			text: app.version
		}
	}
	if args.len == 0 || is_help_request(args) {
		return Outcome{
			kind: .help
			text: app.help(args, style)!
		}
	}

	invocation := parse(app.root, args)!
	return Outcome{
		kind: .executed
		code: dispatch(invocation, app.commands)!
	}
}

fn is_help_request(args []string) bool {
	if args.len > 0 && args[0] == 'help' {
		return true
	}
	for arg in args {
		if arg == '--help' || arg == '-h' {
			return true
		}
		if arg == '--' {
			return false
		}
	}
	return false
}
