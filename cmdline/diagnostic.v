module cmdline

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

// render is the diagnostic as text, and writes it nowhere. Its kind is a
// word as well as a colour, so it reads the same copied into an editor.
pub fn (d Diagnostic) render(style Style) string {
	mut output := '${style.label(d.kind)} ${d.message}\n'
	if d.hint != '' {
		output += '${style.muted('Hint:')} ${d.hint}\n'
	}
	if d.usage != '' {
		output += '${style.muted('Usage:')} ${highlight(d.usage, true, style)}\n'
	}
	return output
}

pub fn render_error(message string, style Style) string {
	return Diagnostic{
		kind:    .error
		message: message
	}.render(style)
}
