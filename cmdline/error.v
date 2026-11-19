module cmdline

// ParseError is a command line that doesn't fit the command tree: the
// adapter's exit status 2. Its message is what was wrong, and usage the
// command's usage line, kept apart so that a diagnostic can show it as one.
pub struct ParseError {
	Error
pub:
	message string
	usage   string
}

pub fn (e ParseError) msg() string {
	return e.message
}

// diagnostic is the error as a diagnostic to render, with its usage line.
pub fn (e ParseError) diagnostic() Diagnostic {
	return Diagnostic{
		kind:    .error
		message: e.message
		usage:   e.usage
	}
}

fn parse_error(message string, usage string) IError {
	return ParseError{
		message: message
		usage:   usage
	}
}
