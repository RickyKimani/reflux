module cli

// Context carries one already-parsed command to a handler. It can grow with
// application-owned services without exposing raw argv.
pub struct Context {
pub:
	invocation Invocation
}

// Handler executes one already-parsed command and does not own process
// termination.
pub interface Handler {
	execute(Context) !int
}

pub struct RegisteredCommand {
pub:
	path    []string
	handler Handler
}

// dispatch selects a handler by its canonical command path.
pub fn dispatch(invocation Invocation, commands []RegisteredCommand) !int {
	for command in commands {
		if command.path == invocation.command {
			return command.handler.execute(Context{
				invocation: invocation
			})
		}
	}
	return error('no handler for ${invocation.command.join(' ')}')
}
