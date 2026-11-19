module cmdline

// Command describes one node in the command tree. It holds metadata only:
// parse, help and dispatch read it.
pub struct Command {
pub:
	name        string
	aliases     []string
	summary     string   // one line, in its parent's list of commands
	description string   // paragraphs at the top of its own help
	usage       string   // such as `reflux solve FILE [--json OUT]`
	examples    []string // lines: one that starts with `#` is a comment
	arguments   []Argument
	flags       []Flag
	commands    []Command
	sections    []Section // more help, after the examples
}

pub struct Argument {
pub:
	name        string
	description string
	required    bool
	// variadic takes every positional argument left. Only the last can be.
	variadic bool
}

pub struct Flag {
pub:
	name        string // used as --name
	alias       string // used as -alias
	description string
	takes_value bool
	value_name  string // the value's name in help: VALUE if empty
}

// Section is a titled block of text at the end of a command's help, such as
// what its exit statuses mean.
pub struct Section {
pub:
	title string
	text  string
}

fn (f Flag) value_label() string {
	return if f.value_name == '' { 'VALUE' } else { f.value_name }
}
