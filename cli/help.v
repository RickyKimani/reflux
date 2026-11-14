module cli

// help renders deterministic plain-text help from command metadata.
pub fn help(command Command) string {
	return render_help(command, Style{
		terminal: Terminal{}
	})
}

// render_help renders help with presentation controlled by the supplied
// terminal capabilities. It does not parse arguments or execute commands.
pub fn render_help(command Command, style Style) string {
	mut lines := []string{}
	if command.usage != '' {
		lines << style.heading('Usage:')
		lines << '  ${command.usage}'
	}
	if command.description != '' {
		if lines.len > 0 {
			lines << ''
		}
		lines << command.description
	}
	if command.commands.len > 0 {
		if lines.len > 0 {
			lines << ''
		}
		lines << style.heading('Commands:')
		mut rows := []HelpRow{}
		for child in command.commands {
			name := if child.aliases.len > 0 {
				'${child.name} (${child.aliases.join(', ')})'
			} else {
				child.name
			}
			rows << HelpRow{
				label:       name
				description: child.summary
			}
		}
		for row in rows {
			lines << format_row(row, rows)
		}
	}
	if command.arguments.len > 0 {
		if lines.len > 0 {
			lines << ''
		}
		lines << style.heading('Arguments:')
		mut rows := []HelpRow{}
		for argument in command.arguments {
			required := if argument.required { ' (required)' } else { '' }
			rows << HelpRow{
				label:       '${argument.name}${required}'
				description: argument.description
			}
		}
		for row in rows {
			lines << format_row(row, rows)
		}
	}
	if command.flags.len > 0 {
		if lines.len > 0 {
			lines << ''
		}
		lines << style.heading('Options:')
		mut rows := []HelpRow{}
		for flag in command.flags {
			mut name := if flag.alias != '' {
				'--${flag.name}, -${flag.alias}'
			} else {
				'--${flag.name}'
			}
			if flag.takes_value {
				name += ' VALUE'
			}
			rows << HelpRow{
				label:       name
				description: flag.description
			}
		}
		for row in rows {
			lines << format_row(row, rows)
		}
	}
	if command.examples.len > 0 {
		if lines.len > 0 {
			lines << ''
		}
		lines << style.heading('Examples:')
		for example in command.examples {
			lines << '  ${example}'
		}
	}
	return lines.join('\n') + '\n'
}

struct HelpRow {
	label       string
	description string
}

fn format_row(row HelpRow, rows []HelpRow) string {
	mut width := 0
	for other in rows {
		if other.label.len > width {
			width = other.label.len
		}
	}
	return '  ${row.label}${' '.repeat(width - row.label.len + 2)}${row.description}'
}
