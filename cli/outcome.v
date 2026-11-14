module cli

pub enum OutcomeKind {
	executed
	help
	version
}

pub struct Outcome {
pub:
	kind OutcomeKind
	code int
	text string
}
