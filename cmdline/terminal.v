module cmdline

import os
import term

const default_terminal_width = 80

pub struct Terminal {
pub:
	interactive bool
	color       bool
	width       int
}

// terminal detects what stdout can show. Tests construct a Terminal directly
// when they need output that doesn't depend on where it runs.
pub fn terminal() Terminal {
	return detect(1)
}

// stderr_terminal detects what stderr can show, for diagnostics: stdout can
// be a file while stderr is still the terminal.
pub fn stderr_terminal() Terminal {
	return detect(2)
}

// detect is colour when the stream is a terminal that understands escape
// sequences, unless NO_COLOR is set (no-color.org); CLICOLOR_FORCE asks for
// colour even into a pipe, for `reflux help | less -R`.
//
// On Windows the terminal must also have virtual-terminal processing on.
// V's runtime turns it on for stdout and stderr when stdout is a console,
// and term.can_show_color_on_stderr() reads that mode, so a stderr whose
// stdout is redirected gets no escapes it can't show.
fn detect(fd int) Terminal {
	interactive := os.is_atty(fd) > 0
	mut width := default_terminal_width
	if interactive {
		columns, _ := term.get_terminal_size()
		if columns > 0 {
			width = columns
		}
	}
	shows_escapes := if fd == 2 {
		term.can_show_color_on_stderr()
	} else {
		term.can_show_color_on_stdout()
	}
	forced := os.getenv('CLICOLOR_FORCE') !in ['', '0']
	return Terminal{
		interactive: interactive
		color:       os.getenv('NO_COLOR') == '' && (forced || (interactive && shows_escapes))
		width:       width
	}
}
