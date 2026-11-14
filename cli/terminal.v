module cli

import os
import term

const default_terminal_width = 80

pub struct Terminal {
pub:
	interactive bool
	color       bool
	width       int
}

// terminal detects presentation capabilities for stdout. Tests can construct
// Terminal directly when they need deterministic output.
pub fn terminal() Terminal {
	interactive := os.is_atty(1) > 0
	no_color := os.getenv('NO_COLOR') != ''
	mut width := default_terminal_width
	if interactive {
		width, _ = term.get_terminal_size()
		if width <= 0 {
			width = default_terminal_width
		}
	}
	return Terminal{
		interactive: interactive
		color:       interactive && !no_color && term.can_show_color_on_stdout()
		width:       width
	}
}
