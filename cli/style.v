module cli

pub struct Style {
pub:
	terminal Terminal
}

pub fn (s Style) heading(text string) string {
	return s.apply(text, '1;36')
}

pub fn (s Style) muted(text string) string {
	return s.apply(text, '2')
}

pub fn (s Style) success(text string) string {
	return s.apply(text, '32')
}

pub fn (s Style) warning(text string) string {
	return s.apply(text, '33')
}

pub fn (s Style) failure(text string) string {
	return s.apply(text, '1;31')
}

pub fn (s Style) accent(text string) string {
	return s.apply(text, '36')
}

fn (s Style) apply(text string, code string) string {
	if !s.terminal.color {
		return text
	}
	return '\x1b[${code}m${text}\x1b[0m'
}
