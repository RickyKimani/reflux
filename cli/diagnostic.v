module cli

pub enum DiagnosticKind {
	error
	warning
	info
}

pub struct Diagnostic {
pub:
	kind    DiagnosticKind
	message string
	hint    string
	usage   string
}

pub fn (d Diagnostic) render(style Style) string {
	label := match d.kind {
		.error { style.failure('Error:') }
		.warning { style.warning('Warning:') }
		.info { style.accent('Info:') }
	}
	mut output := '${label} ${d.message}\n'
	if d.hint != '' {
		output += '${style.muted('Hint:')} ${d.hint}\n'
	}
	if d.usage != '' {
		output += '${style.muted('Usage:')} ${d.usage}\n'
	}
	return output
}

pub fn render_error(message string, style Style) string {
	return Diagnostic{
		kind:    .error
		message: message
	}.render(style)
}
