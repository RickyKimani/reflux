module cli

pub struct ParseError {
	Error
pub:
	message string
	usage   string
}

pub fn (e ParseError) msg() string {
	if e.usage == '' {
		return e.message
	}
	return '${e.message}\nUsage: ${e.usage}'
}

fn parse_error(message string) IError {
	return ParseError{
		message: message
	}
}
