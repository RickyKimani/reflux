module cmdline

import math
import strings

// help renders a command's help without colour.
pub fn help(command Command) string {
	return render_help(command, Style{})
}

// render_help renders a command's help: its description, then its usage,
// commands, arguments, options and examples, then its own sections, each
// under a heading. It doesn't parse arguments or run anything.
pub fn render_help(command Command, style Style) string {
	mut blocks := []string{}
	if command.description != '' {
		blocks << command.description
	}
	if command.usage != '' {
		blocks << section('Usage', ['  ${highlight(command.usage, true, style)}'], style)
	}
	if command.commands.len > 0 {
		mut list := []Row{}
		for child in command.commands {
			list << command_row(child, style)
		}
		blocks << section('Commands', rows(list), style)
	}
	if command.arguments.len > 0 {
		mut list := []Row{}
		for argument in command.arguments {
			list << argument_row(argument, style)
		}
		blocks << section('Arguments', rows(list), style)
	}
	if command.flags.len > 0 {
		mut alias_width := 0
		for flag in command.flags {
			if flag.alias != '' {
				alias_width = math.max(alias_width, '-${flag.alias}, '.len)
			}
		}
		mut list := []Row{}
		for flag in command.flags {
			list << flag_row(flag, alias_width, style)
		}
		blocks << section('Options', rows(list), style)
	}
	if command.examples.len > 0 {
		mut lines := []string{}
		for example in command.examples {
			lines << if example == '' { '' } else { '  ${highlight(example, false, style)}' }
		}
		blocks << section('Examples', lines, style)
	}
	for s in command.sections {
		mut lines := []string{}
		for line in s.text.split('\n') {
			lines << if line == '' { '' } else { '  ${line}' }
		}
		blocks << section(s.title, lines, style)
	}
	return blocks.join('\n\n') + '\n'
}

fn section(title string, body []string, style Style) string {
	return '${style.heading(title.to_upper())}\n\n${body.join('\n')}'
}

// Row is a line of a two-column list: a label, and its text aligned after
// the longest label. The width is the plain label's, which the escapes in
// the styled one don't add to.
struct Row {
	label  string
	styled string
	text   string
}

fn rows(list []Row) []string {
	mut width := 0
	for row in list {
		width = math.max(width, utf8_str_visible_length(row.label))
	}
	mut lines := []string{}
	for row in list {
		text := row.text.split('\n')
		pad := ' '.repeat(width - utf8_str_visible_length(row.label) + 2)
		lines << '  ${row.styled}${pad}${text[0]}'.trim_right(' ')
		for more in text[1..] {
			lines << if more == '' { '' } else { '${' '.repeat(width + 4)}${more}' }
		}
	}
	return lines
}

fn command_row(command Command, style Style) Row {
	if command.aliases.len == 0 {
		return Row{
			label:  command.name
			styled: style.command(command.name)
			text:   command.summary
		}
	}
	aliases := '(${command.aliases.join(', ')})'
	return Row{
		label:  '${command.name} ${aliases}'
		styled: '${style.command(command.name)} ${style.muted(aliases)}'
		text:   command.summary
	}
}

// argument_row is FILE for a required argument, [NAME] for an optional one,
// and NAME=y... for one that takes the rest.
fn argument_row(argument Argument, style Style) Row {
	dots := if argument.variadic { '...' } else { '' }
	if argument.required {
		return Row{
			label:  '${argument.name}${dots}'
			styled: '${style.placeholder(argument.name)}${style.muted(dots)}'
			text:   argument.description
		}
	}
	return Row{
		label:  '[${argument.name}${dots}]'
		styled: '${style.muted('[')}${style.placeholder(argument.name)}${style.muted('${dots}]')}'
		text:   argument.description
	}
}

// flag_row is `-j, --json FILE`, with the long names aligned when some flags
// have a short one and others don't.
fn flag_row(flag Flag, alias_width int, style Style) Row {
	mut label := ''
	mut styled := ''
	if flag.alias != '' {
		short := '-${flag.alias}'
		pad := ' '.repeat(alias_width - short.len - 2)
		label = '${short}, ${pad}'
		styled = '${style.flag(short)}${style.muted(',')} ${pad}'
	} else {
		label = ' '.repeat(alias_width)
		styled = label
	}
	label += '--${flag.name}'
	styled += style.flag('--${flag.name}')
	if flag.takes_value {
		label += ' ${flag.value_label()}'
		styled += ' ${style.placeholder(flag.value_label())}'
	}
	return Row{
		label:  label
		styled: styled
		text:   flag.description
	}
}

enum PieceKind {
	space
	word
	punctuation
	placeholder
	quoted
}

struct Piece {
	kind PieceKind
	text string
}

// highlight colours a usage line or an example. The program is its first
// word, and its commands the lowercase words straight after. In a usage line
// a word in capitals is a placeholder, FILE or NAME=y; in an example it is a
// value, CH4=0.2, and keeps the terminal's colour, as file names and numbers
// do. An example's `# ...` line is a comment.
fn highlight(line string, usage bool, style Style) string {
	if !style.terminal.color {
		return line
	}
	if line.trim_space().starts_with('#') {
		return style.comment(line)
	}
	mut out := strings.new_builder(2 * line.len)
	mut words := 0
	mut commanding := true
	for piece in pieces(line) {
		text := piece.text
		match piece.kind {
			.space {
				out.write_string(text)
			}
			.punctuation {
				out.write_string(style.muted(text))
			}
			.placeholder {
				out.write_string(style.placeholder(text))
				commanding = false
			}
			.quoted {
				out.write_string(text)
				commanding = false
			}
			.word {
				if words == 0 {
					out.write_string(style.program(text))
				} else if commanding && is_command_word(text) {
					out.write_string(style.command(text))
				} else {
					commanding = false
					if is_flag_token(text) {
						out.write_string(style.flag(text))
					} else if usage && is_placeholder(text) {
						out.write_string(style.placeholder(text))
					} else {
						out.write_string(text)
					}
				}
				words++
			}
		}
	}
	return out.str()
}

// pieces splits a line into words, quoted strings, <placeholders>,
// punctuation ([, ], | and ...) and the spaces between, which together are
// the line again.
fn pieces(line string) []Piece {
	mut list := []Piece{}
	mut i := 0
	for i < line.len {
		c := line[i]
		mut j := i + 1
		mut kind := PieceKind.word
		if c == ` ` {
			for j < line.len && line[j] == ` ` {
				j++
			}
			kind = .space
		} else if c == `'` || c == `"` {
			for j < line.len && line[j] != c {
				j++
			}
			j = math.min(j + 1, line.len)
			kind = .quoted
		} else if c == `<` {
			for j < line.len && line[j] != `>` {
				j++
			}
			j = math.min(j + 1, line.len)
			kind = .placeholder
		} else if line[i..].starts_with('...') {
			j = i + 3
			kind = .punctuation
		} else if is_punctuation(c) {
			kind = .punctuation
		} else {
			for j < line.len && line[j] != ` ` && !is_punctuation(line[j])
				&& !line[j..].starts_with('...') {
				j++
			}
		}
		list << Piece{
			kind: kind
			text: line[i..j]
		}
		i = j
	}
	return list
}

fn is_punctuation(c u8) bool {
	return c == `[` || c == `]` || c == `|`
}

// is_command_word is whether a word could name a command: lowercase letters,
// digits and dashes, with no dot or `=` to make it a file or a value.
fn is_command_word(word string) bool {
	if word == '' || !word[0].is_letter() || word[0].is_capital() {
		return false
	}
	for c in word {
		if !(c.is_digit() || c == `-` || (c.is_letter() && !c.is_capital())) {
			return false
		}
	}
	return true
}

// is_placeholder is whether a usage line's word stands for a value: capitals,
// digits and underscores, before any `=`, as FILE, N and NAME=y are.
fn is_placeholder(word string) bool {
	name := word.all_before('=')
	if name == '' || !name[0].is_capital() {
		return false
	}
	for c in name {
		if !(c.is_capital() || c.is_digit() || c == `_`) {
			return false
		}
	}
	return true
}

fn is_flag_token(word string) bool {
	return word.len > 1 && word[0] == `-` && !is_negative_number(word)
}

// is_negative_number is whether a token that starts with a dash is a number,
// -1 or -.5, and so a value rather than an option.
fn is_negative_number(token string) bool {
	return token.len > 1 && token[0] == `-`
		&& (token[1].is_digit() || (token.len > 2 && token[1] == `.` && token[2].is_digit()))
}
