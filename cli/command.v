module cli

// Command describes one node in the CLI command tree. It contains metadata
// only; parsing and execution are introduced in later milestones.
pub struct Command {
pub:
	name        string
	aliases     []string
	summary     string
	description string
	usage       string
	examples    []string
	arguments   []Argument
	flags       []Flag
	commands    []Command
}

pub struct Argument {
pub:
	name        string
	description string
	required    bool
}

pub struct Flag {
pub:
	name        string
	alias       string
	description string
	takes_value bool
}
