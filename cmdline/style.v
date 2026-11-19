module cmdline

// Theme is the palette: an SGR code for each part of the output. It is the
// only place colours are named, so changing it changes how everything looks
// and nothing else.
pub struct Theme {
pub:
	heading     string // a section's title, drawn as a label
	program     string // the program's name in a usage line or an example
	command     string
	flag        string
	placeholder string // FILE, <command>, NAME=y in a usage line
	comment     string // an example's `# ...` line
	muted       string // brackets, `...`, aliases, Hint: and Usage:
	error       string // the label on an error
	warning     string // the label on a warning
	info        string // the label on a note
	success     string
}

// reflux_theme is reflux's palette, in the xterm 256 colours that every
// current terminal has, Windows Terminal and macOS's Terminal among them. It
// takes the two sides of a reflux condenser: copper for the heat, in the
// headings and the program's name, and teal for the cooling water, in the
// commands. Flags are a cornflower blue, and placeholders and comments two
// greys. Each is a mid tone, so that it reads on a dark background and a
// light one: the flags, placeholders and comments at 3:1 or better on both.
pub const reflux_theme = Theme{
	heading:     '1;38;5;230;48;5;166' // cream on copper, #ffffd7 on #d75f00
	program:     '1;38;5;173'          // light copper, #d7875f
	command:     '38;5;37'             // teal, #00afaf
	flag:        '38;5;68'             // cornflower blue, #5f87d7
	placeholder: '3;38;5;246'          // italic grey, #949494
	comment:     '38;5;243'            // dark grey, #767676
	muted:       '38;5;244'            // grey, #808080
	error:       '1;38;5;231;48;5;160' // white on red, #d70000
	warning:     '1;38;5;232;48;5;178' // black on gold, #d7af00
	info:        '1;38;5;232;48;5;37'  // black on teal
	success:     '38;5;71'             // green, #5faf5f
}

// Style turns text into what a terminal shows: coloured when it can show
// colour, and the text itself when it can't.
pub struct Style {
pub:
	terminal Terminal
	theme    Theme = reflux_theme
}

// heading is a section's title. In colour it is a label with a space each
// side; without, the title alone.
pub fn (s Style) heading(text string) string {
	if !s.terminal.color {
		return text
	}
	return s.apply(' ${text} ', s.theme.heading)
}

pub fn (s Style) program(text string) string {
	return s.apply(text, s.theme.program)
}

pub fn (s Style) command(text string) string {
	return s.apply(text, s.theme.command)
}

pub fn (s Style) flag(text string) string {
	return s.apply(text, s.theme.flag)
}

pub fn (s Style) placeholder(text string) string {
	return s.apply(text, s.theme.placeholder)
}

pub fn (s Style) comment(text string) string {
	return s.apply(text, s.theme.comment)
}

pub fn (s Style) muted(text string) string {
	return s.apply(text, s.theme.muted)
}

pub fn (s Style) success(text string) string {
	return s.apply(text, s.theme.success)
}

// label is a diagnostic's kind. In colour it is a label like a heading's,
// ` ERROR `; without, `Error:`.
pub fn (s Style) label(kind DiagnosticKind) string {
	word, code := match kind {
		.error { 'Error', s.theme.error }
		.warning { 'Warning', s.theme.warning }
		.info { 'Info', s.theme.info }
	}
	if !s.terminal.color {
		return '${word}:'
	}
	return s.apply(' ${word.to_upper()} ', code)
}

fn (s Style) apply(text string, code string) string {
	if !s.terminal.color || code == '' || text == '' {
		return text
	}
	return '\x1b[${code}m${text}\x1b[0m'
}
